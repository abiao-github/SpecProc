
"""
Simple grating order tracing for spectral reduction pipeline.

Handles order detection and tracing for grating spectrographs,
where orders run horizontally (left to right) across the detector.
"""

import numpy as np
import logging
from typing import Optional
from scipy.signal import find_peaks
from scipy.ndimage import gaussian_filter1d, sobel
from scipy.interpolate import UnivariateSpline
from numpy.polynomial import Chebyshev, Polynomial
from src.core.data_structures import ApertureSet, ApertureLocation
from src.config.config_manager import ConfigManager

logger = logging.getLogger(__name__)


def find_grating_orders_simple(data: np.ndarray, mask: Optional[np.ndarray] = None,
                             config: Optional[ConfigManager] = None,
                             threshold: float = 0.3, trace_degree: int = 3,
                             separation: float = 30.0) -> ApertureSet:
    """
    Find the positions of grating orders on a CCD image.

    This is a simplified algorithm for grating spectrographs where orders
    run horizontally (left to right) across the detector.

    Args:
        data: Image data (2D array)
        mask: Optional bad pixel mask (same shape as data)
        config: Configuration manager
        threshold: Detection threshold (relative to max on collapsed profile)
        trace_degree: Polynomial degree used to fit order trace
        separation: Estimated order separation (in pixels) along Y direction

    Returns:
        ApertureSet with detected orders
    """
    if mask is None:
        mask = np.zeros_like(data, dtype=np.int32)

    h, w = data.shape

    logger.info(f"Finding grating orders in image of shape {data.shape}")

    # Build a trace image directly from combined master flat.
    # Use mild log compression so bright orders and faint order wings are both trackable.
    img = np.asarray(data, dtype=np.float64)
    p10 = np.nanpercentile(img, 10)
    img = np.clip(img - p10, 0.0, None)
    trace_img = np.log1p(img)

    x_center = w // 2
    half_band = max(5, w // 20)
    x1 = max(0, x_center - half_band)
    x2 = min(w, x_center + half_band + 1)
    center_profile = np.median(trace_img[:, x1:x2], axis=1)
    center_profile = gaussian_filter1d(center_profile, sigma=1.5)

    # Seed detection from central columns only (robust enough for this dataset).
    center_p50 = np.nanmedian(center_profile)
    center_p95 = np.nanpercentile(center_profile, 95)
    center_scale = max(1e-6, center_p95 - center_p50)
    seed_profile = (center_profile - center_p50) / center_scale

    # Configurable tracing controls.
    seed_threshold = threshold
    prominence_scale = 0.5
    search_half_scale = 0.45
    step_denominator = 220
    fill_missing_orders = True
    gap_fill_factor = 1.35
    fit_method = 'polynomial'
    edge_degree = 3
    bspline_smooth = 0.2
    if config is not None:
        seed_threshold = config.get_float('reduce.trace', 'seed_threshold', threshold)
        prominence_scale = config.get_float('reduce.trace', 'prominence_scale', 0.5)
        search_half_scale = config.get_float('reduce.trace', 'search_half_scale', 0.45)
        step_denominator = config.get_int('reduce.trace', 'step_denominator', 220)
        fill_missing_orders = config.get_bool('reduce.trace', 'fill_missing_orders', True)
        gap_fill_factor = config.get_float('reduce.trace', 'gap_fill_factor', 1.35)
        fit_method = config.get('reduce.trace', 'fit_method', 'polynomial').strip().lower()
        edge_degree = config.get_int('reduce.trace', 'edge_degree', 3)
        bspline_smooth = config.get_float('reduce.trace', 'bspline_smooth', 0.2)

    # Find seed peaks on center profile using two-pass detection (strong + weak peaks).
    # This keeps obvious peaks while recovering dim orders that still stand out from background.
    prof_max = float(np.max(seed_profile))
    baseline = float(np.median(seed_profile))
    noise = float(np.std(seed_profile - baseline))
    min_dist = max(6, int(separation * 0.72))

    strong_height = max(seed_threshold * prof_max, 0.16 * prof_max)
    strong_prom = max(0.015 * prof_max, prominence_scale * noise)
    peaks_strong, _ = find_peaks(
        seed_profile,
        height=strong_height,
        distance=min_dist,
        prominence=strong_prom,
    )

    weak_height = max(0.08 * prof_max, baseline + 0.8 * noise)
    weak_prom = max(0.008 * prof_max, 0.35 * prominence_scale * noise)
    peaks_weak, _ = find_peaks(
        seed_profile,
        height=weak_height,
        distance=max(4, int(0.6 * min_dist)),
        prominence=weak_prom,
    )

    # Merge strong/weak peaks while enforcing a minimum spacing.
    merged = np.asarray(sorted(np.unique(np.concatenate([peaks_strong, peaks_weak]))), dtype=int)
    peaks = []
    for p in merged:
        if not peaks:
            peaks.append(int(p))
            continue
        if p - peaks[-1] >= min_dist:
            peaks.append(int(p))
        else:
            # Keep the stronger peak when too close.
            if seed_profile[p] > seed_profile[peaks[-1]]:
                peaks[-1] = int(p)
    peaks = np.asarray(peaks, dtype=int)

    # Prune obvious false seeds (typically faint noise peaks near tails).
    if len(peaks) >= 4:
        snr = (seed_profile[peaks] - baseline) / (noise + 1e-6)
        diffs = np.diff(peaks)
        med_sep = np.median(diffs) if len(diffs) > 0 else separation
        keep = np.zeros_like(peaks, dtype=bool)

        for i, p in enumerate(peaks):
            if snr[i] >= 1.35:
                keep[i] = True
                continue

            # Keep marginal peaks only when spacing is consistent with neighbors.
            left_ok = (
                i > 0
                and 0.55 * med_sep <= (p - peaks[i - 1]) <= 1.45 * med_sep
            )
            right_ok = (
                i < len(peaks) - 1
                and 0.55 * med_sep <= (peaks[i + 1] - p) <= 1.45 * med_sep
            )
            if snr[i] >= 1.0 and (left_ok or right_ok):
                keep[i] = True

        dropped = int(np.sum(~keep))
        if dropped > 0:
            peaks = peaks[keep]
            logger.info(f"Pruned {dropped} low-confidence seed peak(s)")

    logger.info(
        f"Found {len(peaks)} seed peaks in center profile "
        f"(strong={len(peaks_strong)}, weak={len(peaks_weak)})"
    )

    # Fill missed weak orders by checking abnormal gaps in seed sequence.
    # Use a relatively sensitive default so weak orders are not skipped easily.
    peaks = np.asarray(sorted(peaks), dtype=int)
    if fill_missing_orders and len(peaks) >= 2:
        diffs = np.diff(peaks)
        med_sep = np.median(diffs)
        if np.isfinite(med_sep) and med_sep > 2:
            inserted = []
            for i in range(len(peaks) - 1):
                left, right = peaks[i], peaks[i + 1]
                gap = right - left
                if gap > gap_fill_factor * med_sep:
                    n_missing = int(round(gap / med_sep)) - 1
                    for k in range(1, n_missing + 1):
                        guess = int(round(left + k * gap / (n_missing + 1)))
                        win = max(5, int(0.55 * med_sep))
                        y1 = max(0, guess - win)
                        y2 = min(h, guess + win + 1)
                        if y2 - y1 < 3:
                            continue
                        seg = seed_profile[y1:y2]
                        noise = np.std(seg - np.median(seg))
                        # Accept weaker local maxima as long as they are above local noise.
                        loc = int(np.argmax(seg))
                        cand = y1 + loc
                        cand_val = seg[loc]
                        if noise > 0 and cand_val < (np.median(seg) + 1.2 * noise):
                            continue
                        if np.all(np.abs(peaks - cand) > max(3, int(0.30 * med_sep))):
                            inserted.append(cand)

            # Extrapolate possible missed dim orders near both ends.
            edge_win = max(5, int(0.60 * med_sep))
            min_sep = max(3, int(0.30 * med_sep))
            left_ref = int(peaks[0])
            while left_ref - med_sep > 2:
                guess = int(round(left_ref - med_sep))
                y1 = max(0, guess - edge_win)
                y2 = min(h, guess + edge_win + 1)
                if y2 - y1 < 3:
                    break
                seg = seed_profile[y1:y2]
                cand = y1 + int(np.argmax(seg))
                seg_noise = float(np.std(seg - np.median(seg)))
                if seg_noise > 0 and seed_profile[cand] >= (np.median(seg) + 1.0 * seg_noise):
                    if np.all(np.abs(peaks - cand) > min_sep):
                        inserted.append(cand)
                        left_ref = cand
                        continue
                break

            right_ref = int(peaks[-1])
            while right_ref + med_sep < h - 2:
                guess = int(round(right_ref + med_sep))
                y1 = max(0, guess - edge_win)
                y2 = min(h, guess + edge_win + 1)
                if y2 - y1 < 3:
                    break
                seg = seed_profile[y1:y2]
                cand = y1 + int(np.argmax(seg))
                seg_noise = float(np.std(seg - np.median(seg)))
                if seg_noise > 0 and seed_profile[cand] >= (np.median(seg) + 1.0 * seg_noise):
                    if np.all(np.abs(peaks - cand) > min_sep):
                        inserted.append(cand)
                        right_ref = cand
                        continue
                break

            if inserted:
                peaks = np.asarray(sorted(np.concatenate([peaks, np.array(inserted, dtype=int)])), dtype=int)
                logger.info(f"Inserted {len(inserted)} missing-order seed(s) by gap filling")

    if len(peaks) == 0:
        logger.warning("No order peaks detected; lower threshold or check flat exposures")
        return ApertureSet()

    def refine_peak_y(col: np.ndarray, y_guess: float, search_half: int) -> Optional[float]:
        y0 = int(round(y_guess))
        y1 = max(0, y0 - search_half)
        y2 = min(h, y0 + search_half + 1)
        if y2 - y1 < 5:
            return None

        segment = col[y1:y2]
        loc = int(np.argmax(segment))
        peak_val = segment[loc]
        if not np.isfinite(peak_val):
            return None

        # local centroid around peak for sub-pixel refinement
        c1 = max(0, loc - 2)
        c2 = min(segment.size, loc + 3)
        wv = segment[c1:c2]
        idx = np.arange(c1, c2, dtype=float)
        s = np.sum(wv)
        if s <= 0:
            return float(y1 + loc)
        return float(y1 + np.sum(idx * wv) / s)

    def trace_order(seed_y: float) -> tuple:
        # Trace both directions from center with continuity + short-gap tolerance.
        step = max(1, w // max(40, step_denominator))
        search_half = max(6, int(separation * search_half_scale))

        xs = [x_center]
        ys = [float(seed_y)]

        # right
        y_prev = float(seed_y)
        miss_count = 0
        right_x = [float(x_center)]
        right_y = [float(seed_y)]
        for x in range(x_center + step, w, step):
            if len(right_x) >= 4:
                recent_x = np.asarray(right_x[-4:], dtype=float)
                recent_y = np.asarray(right_y[-4:], dtype=float)
                coef1 = np.polyfit(recent_x, recent_y, deg=1)
                y_guess = float(np.polyval(coef1, x))
            else:
                y_guess = y_prev

            y_new = refine_peak_y(trace_img[:, x], y_guess, search_half)
            if y_new is None:
                miss_count += 1
                if miss_count > 6:
                    break
                continue

            if abs(y_new - y_guess) > search_half:
                miss_count += 1
                if miss_count > 6:
                    break
                continue

            xs.append(x)
            ys.append(y_new)
            right_x.append(float(x))
            right_y.append(float(y_new))
            y_prev = y_new
            miss_count = 0

        # left
        y_prev = float(seed_y)
        miss_count = 0
        left_x = [float(x_center)]
        left_y = [float(seed_y)]
        for x in range(x_center - step, -1, -step):
            if len(left_x) >= 4:
                recent_x = np.asarray(left_x[-4:], dtype=float)
                recent_y = np.asarray(left_y[-4:], dtype=float)
                coef1 = np.polyfit(recent_x, recent_y, deg=1)
                y_guess = float(np.polyval(coef1, x))
            else:
                y_guess = y_prev

            y_new = refine_peak_y(trace_img[:, x], y_guess, search_half)
            if y_new is None:
                miss_count += 1
                if miss_count > 6:
                    break
                continue

            if abs(y_new - y_guess) > search_half:
                miss_count += 1
                if miss_count > 6:
                    break
                continue

            xs.append(x)
            ys.append(y_new)
            left_x.append(float(x))
            left_y.append(float(y_new))
            y_prev = y_new
            miss_count = 0

        xs = np.asarray(xs, dtype=float)
        ys = np.asarray(ys, dtype=float)
        idx = np.argsort(xs)
        return xs[idx], ys[idx]

    def estimate_order_width(xs: np.ndarray, ys: np.ndarray) -> float:
        """Estimate aperture width from local FWHM samples along the traced order."""
        if xs.size == 0:
            return float(separation)

        sample_idx = np.linspace(0, xs.size - 1, num=min(24, xs.size), dtype=int)
        widths = []
        search_half = max(8, int(0.8 * separation))

        for i in sample_idx:
            x = int(round(xs[i]))
            yc = float(ys[i])
            y1 = max(0, int(np.floor(yc - search_half)))
            y2 = min(h, int(np.ceil(yc + search_half + 1)))
            if y2 - y1 < 7:
                continue

            prof = trace_img[y1:y2, x]
            pmax = np.max(prof)
            if not np.isfinite(pmax) or pmax <= 0:
                continue
            half = 0.5 * pmax

            iy = int(np.argmax(prof))
            l = iy
            while l > 0 and prof[l] > half:
                l -= 1
            r = iy
            while r < prof.size - 1 and prof[r] > half:
                r += 1

            fwhm = float(r - l)
            if 2.0 <= fwhm <= 2.5 * separation:
                widths.append(fwhm)

        if not widths:
            return float(separation)

        # Use aperture full width ~ 1.8*FWHM, clipped to sensible range.
        width = 1.8 * float(np.median(widths))
        width = float(np.clip(width, 6.0, 2.2 * separation))
        return width

    def _fit_values_with_method(x_obs: np.ndarray, y_obs: np.ndarray, x_eval: np.ndarray,
                                method: str, degree: int, bspline_s: float) -> np.ndarray:
        """Fit y(x) with selected method and evaluate on x_eval."""
        x = np.asarray(x_obs, dtype=float)
        y = np.asarray(y_obs, dtype=float)
        xe = np.asarray(x_eval, dtype=float)
        if x.size < 3:
            return np.interp(xe, x, y) if x.size > 1 else np.full_like(xe, y[0] if x.size == 1 else 0.0)

        degree = max(1, int(min(degree, max(1, len(x) - 2))))
        m = (method or 'polynomial').strip().lower()

        if m == 'chebyshev':
            cheb = Chebyshev.fit(x, y, deg=degree, domain=[x.min(), x.max()])
            return cheb(xe)

        if m == 'bspline':
            k = min(3, max(1, len(x) - 1))
            s_val = max(0.0, float(bspline_s)) * len(x)
            spl = UnivariateSpline(x, y, k=k, s=s_val)
            return spl(xe)

        coef = np.polyfit(x, y, deg=degree)
        return np.polyval(coef, xe)

    def _fit_center_coeff(x_vals: np.ndarray, y_vals: np.ndarray, method: str, degree: int,
                          bspline_s: float) -> np.ndarray:
        """Fit center and return poly coeffs (descending) used by downstream code."""
        x = np.asarray(x_vals, dtype=float)
        xe = np.asarray(x_vals, dtype=float)
        y_fit = _fit_values_with_method(x, y_vals, xe, method, degree, bspline_s)
        return np.array(np.polyfit(xe, y_fit, deg=max(1, int(min(degree, max(1, len(xe) - 2))))), dtype=float)

    def _fit_edge_coeff_with_center_prior(x_edge: np.ndarray, y_edge: np.ndarray,
                                          center_coef: np.ndarray,
                                          method: str,
                                          edge_deg: int,
                                          bspline_s: float,
                                          side: str) -> np.ndarray:
        """
        Fit one edge with full parameters, using center trajectory only as prior init.

        This performs complete curve fitting on edge points; center trajectory is
        used to build an initial model and robustly select inliers in early iterations.
        """
        x = np.asarray(x_edge, dtype=float)
        y = np.asarray(y_edge, dtype=float)
        if x.size < max(6, edge_deg + 2):
            return None

        deg = max(1, int(min(edge_deg, max(1, len(x) - 2))))

        # Prior initialization from center coefficients with vertical shift.
        init_coef = np.array(center_coef, dtype=float).copy()
        if init_coef.size != deg + 1:
            # Refit center prior to target degree for consistent parameter size.
            y_center = np.polyval(center_coef, x)
            init_coef = np.polyfit(x, y_center, deg=deg)

        y_center = np.polyval(init_coef, x)
        d_med = np.median(y - y_center)
        init_coef[-1] += d_med

        # Robust inlier selection seeded by center-prior model.
        keep = np.ones_like(x, dtype=bool)
        model = np.polyval(init_coef, x)
        resid = y - model
        s0 = np.std(resid) if resid.size > 1 else 0.0
        if s0 > 0:
            keep = np.abs(resid) < 3.5 * s0

        m = (method or 'polynomial').strip().lower()
        coef = init_coef.copy()
        for _ in range(3):
            xk = x[keep]
            yk = y[keep]
            if xk.size < max(6, deg + 2):
                xk = x
                yk = y

            if m == 'chebyshev':
                cheb = Chebyshev.fit(xk, yk, deg=deg, domain=[x.min(), x.max()])
                poly = cheb.convert(kind=Polynomial)
                coef = np.array(poly.coef[::-1], dtype=float)
            elif m == 'bspline':
                k = min(3, max(1, len(xk) - 1))
                s_val = max(0.0, float(bspline_s)) * len(xk)
                spl = UnivariateSpline(xk, yk, k=k, s=s_val)
                y_smooth = spl(x)
                coef = np.array(np.polyfit(x, y_smooth, deg=deg), dtype=float)
            else:
                coef = np.array(np.polyfit(xk, yk, deg=deg), dtype=float)

            model = np.polyval(coef, x)
            resid = y - model
            s = np.std(resid[keep]) if np.any(keep) else np.std(resid)
            if s <= 0:
                break
            new_keep = np.abs(resid) < 3.0 * s
            if np.array_equal(new_keep, keep):
                break
            keep = new_keep

        return coef

    def _detect_order_edges(center_coef: np.ndarray, width_guess: float,
                            step: int) -> tuple:
        """
        Detect lower/upper edge points using Sobel positive/negative peaks.

        Lower edge: strongest positive Sobel peak below center.
        Upper edge: strongest negative Sobel peak above center.
        """
        grad = sobel(trace_img, axis=0, mode='nearest')
        x_samples = np.arange(0, w, max(1, step))
        y_center = np.polyval(center_coef, x_samples)
        half = max(6, int(round(0.7 * width_guess)))

        xl, yl, xu, yu = [], [], [], []
        for x, yc in zip(x_samples, y_center):
            y0 = int(round(yc))
            low1 = max(0, y0 - half)
            low2 = min(h, y0 + 1)
            up1 = max(0, y0)
            up2 = min(h, y0 + half + 1)

            if low2 - low1 >= 4:
                g_low = grad[low1:low2, int(x)]
                i_low = int(np.argmax(g_low))
                xl.append(float(x))
                yl.append(float(low1 + i_low))

            if up2 - up1 >= 4:
                g_up = -grad[up1:up2, int(x)]
                i_up = int(np.argmax(g_up))
                xu.append(float(x))
                yu.append(float(up1 + i_up))

        return np.array(xl), np.array(yl), np.array(xu), np.array(yu)

    def _sanitize_edge_coeffs(center_coef: np.ndarray,
                              lower_coef: np.ndarray,
                              upper_coef: np.ndarray,
                              width_guess: float) -> tuple:
        """
        Enforce physical edge constraints so boundaries remain around center trajectory.

        If fitted edges drift/cross or become implausibly wide/narrow, fall back to a
        stable center +/- offset model inferred from valid portions.
        """
        x_eval = np.linspace(0, w - 1, num=min(200, w), dtype=float)
        yc = np.polyval(center_coef, x_eval)
        yl = np.polyval(lower_coef, x_eval)
        yu = np.polyval(upper_coef, x_eval)

        target_half = max(3.0, 0.5 * float(width_guess))
        min_half = max(2.0, 0.45 * target_half)
        max_half = max(min_half + 1.0, 2.0 * target_half)

        d_low = yc - yl
        d_up = yu - yc

        valid_low = np.isfinite(d_low) & (d_low > min_half) & (d_low < max_half)
        valid_up = np.isfinite(d_up) & (d_up > min_half) & (d_up < max_half)
        valid_both = valid_low & valid_up

        frac_valid = float(np.mean(valid_both)) if valid_both.size > 0 else 0.0

        # If edge fits are mostly reasonable, keep them.
        if frac_valid >= 0.75:
            return lower_coef, upper_coef

        # Otherwise build robust symmetric-like edges from valid medians.
        if np.any(valid_low):
            low_off = float(np.median(d_low[valid_low]))
        else:
            low_off = target_half

        if np.any(valid_up):
            up_off = float(np.median(d_up[valid_up]))
        else:
            up_off = target_half

        low_off = float(np.clip(low_off, min_half, max_half))
        up_off = float(np.clip(up_off, min_half, max_half))

        lower_fix = center_coef.copy()
        upper_fix = center_coef.copy()
        lower_fix[-1] -= low_off
        upper_fix[-1] += up_off

        logger.debug(
            "Edge coefficients sanitized: frac_valid=%.2f, low_off=%.2f, up_off=%.2f",
            frac_valid, low_off, up_off
        )
        return lower_fix, upper_fix

    def _fit_boundary_coeff(x_vals: np.ndarray, y_vals: np.ndarray,
                            method: str, degree: int, bspline_s: float) -> np.ndarray:
        """Fit a smooth shared inter-order boundary y(x)."""
        x = np.asarray(x_vals, dtype=float)
        y = np.asarray(y_vals, dtype=float)
        keep = np.isfinite(x) & np.isfinite(y)
        x = x[keep]
        y = y[keep]
        if x.size < max(6, degree + 2):
            return None

        deg = max(1, int(min(degree, max(1, x.size - 2))))
        coef = np.polyfit(x, y, deg=deg)
        for _ in range(3):
            model = np.polyval(coef, x)
            resid = y - model
            std = np.std(resid)
            if std <= 0:
                break
            good = np.abs(resid) < 3.0 * std
            if np.sum(good) < max(6, deg + 2):
                break
            new_coef = np.polyfit(x[good], y[good], deg=deg)
            if np.allclose(new_coef, coef, rtol=1e-5, atol=1e-5):
                coef = new_coef
                break
            coef = new_coef
        return np.array(coef, dtype=float)

    def _build_shared_boundary_coeffs(cands: list, method: str,
                                      edge_deg: int, bspline_s: float) -> list:
        """Build aperture edges from shared inter-order derivative/valley boundaries.

        For each adjacent order pair, determine a single shared boundary by
        locating the first-derivative zero crossing (negative -> positive) between
        their center traces in sampled columns. If no robust zero crossing exists,
        fall back to the local minimum and then to the center midline.
        This is more stable than fitting each order wing independently when
        adjacent orders are bright, steep, and separated by only 1-2 pixels,
        while naturally allowing the boundary spacing to widen in later orders.
        """
        if len(cands) == 0:
            return cands
        if len(cands) == 1:
            return cands

        x_samples = np.arange(0, w, max(1, w // max(40, step_denominator)), dtype=float)
        smooth_img = gaussian_filter1d(trace_img, sigma=1.0, axis=0)

        shared_bounds = []
        for i in range(len(cands) - 1):
            lower_c = cands[i]
            upper_c = cands[i + 1]
            xb = []
            yb = []

            yc1 = np.polyval(lower_c['center_coef'], x_samples)
            yc2 = np.polyval(upper_c['center_coef'], x_samples)
            for x, y1c, y2c in zip(x_samples, yc1, yc2):
                if (not np.isfinite(y1c)) or (not np.isfinite(y2c)):
                    continue
                ya = float(min(y1c, y2c))
                ybnd = float(max(y1c, y2c))
                if ybnd - ya < 2.0:
                    continue

                y1i = max(0, int(np.floor(ya)))
                y2i = min(h - 1, int(np.ceil(ybnd)))
                if y2i - y1i < 2:
                    continue

                col = smooth_img[:, int(round(x))]
                seg = col[y1i:y2i + 1]
                if seg.size < 3 or not np.any(np.isfinite(seg)):
                    boundary = 0.5 * (ya + ybnd)
                else:
                    # Ignore one pixel around each center so the boundary is not
                    # pulled back onto the order core when valleys are narrow.
                    trim_lo = 1 if seg.size >= 5 else 0
                    trim_hi = seg.size - 1 if seg.size >= 5 else seg.size
                    inner = seg[trim_lo:trim_hi]
                    if inner.size >= 2 and np.any(np.isfinite(inner)):
                        grad_inner = np.gradient(inner.astype(float))
                        zero_crossings = []
                        for j in range(len(grad_inner) - 1):
                            g0 = float(grad_inner[j])
                            g1 = float(grad_inner[j + 1])
                            if (not np.isfinite(g0)) or (not np.isfinite(g1)):
                                continue
                            # Valley condition: derivative changes from negative to positive.
                            if g0 <= 0.0 and g1 >= 0.0 and (g1 - g0) > 0.0:
                                denom = (g1 - g0)
                                frac = 0.0 if abs(denom) < 1e-12 else (-g0 / denom)
                                frac = float(np.clip(frac, 0.0, 1.0))
                                y_zero = float(y1i + trim_lo + j + frac)
                                # Use linear interpolation of intensity at zero crossing.
                                f0 = float(inner[j])
                                f1 = float(inner[j + 1])
                                f_zero = (1.0 - frac) * f0 + frac * f1
                                zero_crossings.append((f_zero, y_zero))

                        if zero_crossings:
                            # Pick the lowest-intensity zero crossing between adjacent peaks.
                            _, boundary = min(zero_crossings, key=lambda t: t[0])
                        else:
                            # Fallback to local minimum inside the inter-order interval.
                            boundary = float(y1i + trim_lo + int(np.nanargmin(inner)))
                    else:
                        boundary = 0.5 * (ya + ybnd)

                # Keep the shared boundary strictly between the two centers.
                boundary = float(np.clip(boundary, ya + 0.35, ybnd - 0.35))
                xb.append(float(x))
                yb.append(boundary)

            coef = _fit_boundary_coeff(np.asarray(xb), np.asarray(yb), method, edge_deg, bspline_s)
            if coef is None:
                # Midline fallback if valley points are too sparse.
                mid_y = 0.5 * (
                    np.polyval(lower_c['center_coef'], x_samples) +
                    np.polyval(upper_c['center_coef'], x_samples)
                )
                coef = _fit_boundary_coeff(x_samples, mid_y, 'polynomial', edge_deg, 0.0)
            shared_bounds.append(coef)

        refined = []
        for i, cand in enumerate(cands):
            center_coef = cand['center_coef']
            x_eval = np.linspace(0, w - 1, num=min(200, w), dtype=float)
            yc = np.polyval(center_coef, x_eval)

            if i == 0:
                upper_coef = shared_bounds[0]
                yu = np.polyval(upper_coef, x_eval)
                up_off = np.nanmedian(yu - yc)
                if not np.isfinite(up_off) or up_off <= 1:
                    up_off = max(3.0, 0.5 * cand['width'])
                lower_coef = center_coef.copy()
                lower_coef[-1] -= float(up_off)
            elif i == len(cands) - 1:
                lower_coef = shared_bounds[i - 1]
                yl = np.polyval(lower_coef, x_eval)
                low_off = np.nanmedian(yc - yl)
                if not np.isfinite(low_off) or low_off <= 1:
                    low_off = max(3.0, 0.5 * cand['width'])
                upper_coef = center_coef.copy()
                upper_coef[-1] += float(low_off)
            else:
                lower_coef = shared_bounds[i - 1]
                upper_coef = shared_bounds[i]

            lower_coef, upper_coef = _sanitize_edge_coeffs(center_coef, lower_coef, upper_coef, cand['width'])

            yl = np.polyval(lower_coef, x_eval)
            yu = np.polyval(upper_coef, x_eval)
            width_med = float(np.nanmedian(yu - yl)) if np.any(np.isfinite(yu - yl)) else float(cand['width'])
            width_med = float(np.clip(width_med, 4.0, 2.2 * separation))

            new_cand = dict(cand)
            new_cand['lower_coef'] = lower_coef
            new_cand['upper_coef'] = upper_coef
            new_cand['width'] = width_med
            refined.append(new_cand)

        return refined

    aperture_set = ApertureSet()
    candidates = []
    for order_idx, seed in enumerate(peaks):
        x_positions, y_positions = trace_order(float(seed))

        if len(x_positions) < max(10, trace_degree + 4):
            logger.warning(f"Order {order_idx + 1}: Not enough traced points ({len(x_positions)})")
            continue

        try:
            good = np.ones_like(x_positions, dtype=bool)
            for _ in range(4):
                coef = np.polyfit(x_positions[good], y_positions[good], deg=trace_degree)
                model = np.polyval(coef, x_positions)
                resid = y_positions - model
                std = np.std(resid[good]) if np.any(good) else 0.0
                if std <= 0:
                    break
                new_good = np.abs(resid) < 2.8 * std
                if np.sum(new_good) < max(10, trace_degree + 4):
                    break
                if np.array_equal(new_good, good):
                    break
                good = new_good

            coef = _fit_center_coeff(
                x_positions[good], y_positions[good], fit_method, trace_degree, bspline_smooth
            )

            # Reject likely false traces: too short across dispersion or too noisy around fit.
            x_good = x_positions[good]
            y_good = y_positions[good]
            if x_good.size < max(10, trace_degree + 4):
                logger.warning(f"Order {order_idx + 1}: Rejected (insufficient robust points)")
                continue

            x_coverage = float(np.max(x_good) - np.min(x_good))
            if x_coverage < 0.40 * w:
                logger.warning(
                    f"Order {order_idx + 1}: Rejected (short trace coverage {x_coverage:.0f}px)"
                )
                continue

            y_fit_good = np.polyval(coef, x_good)
            fit_rms = float(np.std(y_good - y_fit_good)) if y_good.size > 1 else 0.0
            if fit_rms > max(1.2, 0.22 * separation):
                logger.warning(
                    f"Order {order_idx + 1}: Rejected (fit rms {fit_rms:.2f}px too large)"
                )
                continue

            width = estimate_order_width(x_positions[good], y_positions[good])

            # Edge coefficients are initialized here and then refined globally
            # from shared inter-order boundaries after duplicate merging.
            half_width = width / 2.0
            lower_coef = coef.copy()
            lower_coef[-1] -= half_width
            upper_coef = coef.copy()
            upper_coef[-1] += half_width

            candidates.append({
                'center_coef': coef,
                'lower_coef': lower_coef,
                'upper_coef': upper_coef,
                'width': float(width),
                'n_good': int(np.sum(good)),
                'fit_rms': float(fit_rms),
                'x_coverage': float(x_coverage),
                'y_center': float(np.polyval(coef, x_center)),
            })

            logger.info(
                f"Candidate from seed {order_idx + 1}: points={np.sum(good)}, "
                f"fit={fit_method}, width={width:.1f}px"
            )

        except Exception as e:
            logger.warning(f"Order {order_idx + 1}: Failed polynomial tracing: {e}")
            continue

    # Merge duplicate candidates that represent the same physical order,
    # then assign continuous order IDs.
    if candidates:
        candidates_sorted = sorted(candidates, key=lambda c: c['y_center'])
        merge_sep = max(3.0, 0.45 * float(separation))
        unique = []

        for cand in candidates_sorted:
            if not unique:
                unique.append(cand)
                continue

            prev = unique[-1]
            if abs(cand['y_center'] - prev['y_center']) < merge_sep:
                prev_score = prev['n_good'] - 2.0 * prev['fit_rms'] + 0.002 * prev['x_coverage']
                cand_score = cand['n_good'] - 2.0 * cand['fit_rms'] + 0.002 * cand['x_coverage']
                if cand_score > prev_score:
                    unique[-1] = cand
            else:
                unique.append(cand)

        unique = _build_shared_boundary_coeffs(unique, fit_method, edge_degree, bspline_smooth)

        for new_id, cand in enumerate(unique, start=1):
            aperture_loc = ApertureLocation(
                aperture=new_id,
                order=new_id,
                center_coef=cand['center_coef'],
                lower_coef=cand['lower_coef'],
                upper_coef=cand['upper_coef'],
                width=cand['width'],
            )
            aperture_set.add_aperture(aperture_loc)

    logger.info(f"Detected {aperture_set.norders} orders (simple grating tracing)")

    return aperture_set
