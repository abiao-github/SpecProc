"""
Scattered light subtraction stage of spectral reduction pipeline.

Estimates and removes inter-order scattered light (stray light) from
echelle spectrograph CCD images.  Equivalent to IRAF's ``apscatter`` task.
"""

import numpy as np
import logging
from pathlib import Path
from typing import Tuple, Optional
from scipy.ndimage import median_filter
import matplotlib.pyplot as plt
import matplotlib.cm as mcm
from astropy.convolution import convolve_fft, Gaussian2DKernel
from src.utils.image_processing import estimate_background_2d, estimate_background_region
from src.config.config_manager import ConfigManager
from src.core.data_structures import ApertureSet
from src.utils.fits_io import write_fits_image

logger = logging.getLogger(__name__)


class BackgroundRemover:
    """Handles background/stray light estimation and correction."""

    def __init__(self):
        """
        Initialize background remover.
        """
        self.background_2d = None

    def _build_aperture_mask(self, image_shape: Tuple[int, int], apertures: ApertureSet,
                             margin_pixels: int = 3) -> np.ndarray:
        """Build boolean mask using traced lower/upper edges plus an outer margin."""
        rows, cols = image_shape
        mask = np.zeros((rows, cols), dtype=bool)
        x_coords = np.arange(cols)

        for aperture in apertures.apertures.values():
            center = aperture.get_position(x_coords)
            lower = aperture.get_lower(x_coords)
            upper = aperture.get_upper(x_coords)

            width_guess = float(aperture.width if aperture.width is not None else 8.0)
            half_w_max = max(4.0, min(0.12 * rows, 2.2 * width_guess))

            for x in x_coords:
                c = float(center[x])
                lo = float(lower[x])
                hi = float(upper[x])

                if (not np.isfinite(c)) or c < -margin_pixels or c > (rows - 1 + margin_pixels):
                    continue
                if (not np.isfinite(lo)) or (not np.isfinite(hi)):
                    continue

                y_lo = min(lo, hi) - margin_pixels
                y_hi = max(lo, hi) + margin_pixels

                # Guard against pathological edge fits that produce unrealistically
                # wide masks and wipe out inter-order regions.
                if (y_hi - y_lo) > 2.0 * half_w_max:
                    y_lo = c - half_w_max
                    y_hi = c + half_w_max

                y1 = max(0, int(np.floor(y_lo)))
                y2 = min(rows, int(np.ceil(y_hi + 1)))
                if y2 > y1:
                    mask[y1:y2, x] = True

        return mask

    @staticmethod
    def _estimate_background_convolution(image: np.ndarray, order_mask: np.ndarray,
                                         kernel_sigma_x: float = 15.0,
                                         kernel_sigma_y: float = 15.0,
                                         sigma_clip: float = 3.0,
                                         maxiters: int = 3,
                                         clip_mode: str = 'upper'
                                         ) -> np.ndarray:
        """Estimate scattered-light background via astropy NaN-aware 2D Gaussian convolution.

        Masked (order) pixels are set to NaN and astropy.convolution.convolve
        automatically ignores them, using only valid inter-order pixels for the
        weighted average under the kernel.  This naturally reproduces local
        structures (ripples, fringes) that global polynomial fits cannot capture.

        Parameters
        ----------
        image : 2D array — science or flat image.
        order_mask : bool mask, True = inside an order (to be ignored).
        kernel_sigma_x, kernel_sigma_y : Gaussian kernel std-dev in pixels
            (x = dispersion direction, y = cross-dispersion).  Use an
            anisotropic kernel: large x (dispersion, ~40) for broad smoothing,
            smaller y (cross-dispersion, ~15) to preserve ripple structure.
        sigma_clip : iterative outlier rejection threshold (MAD-based).
        maxiters : max sigma-clipping iterations.
        clip_mode : 'upper' / 'both' / 'lower'.
        """
        kernel = Gaussian2DKernel(x_stddev=kernel_sigma_x, y_stddev=kernel_sigma_y)
        # Pad size = kernel half-extent so edge-padding produces the same
        # result as boundary='extend' (not available in convolve_fft).
        pad_y = kernel.shape[0] // 2
        pad_x = kernel.shape[1] // 2

        current_mask = order_mask.copy()
        model = None
        for it in range(max(1, maxiters)):
            data = image.astype(np.float64).copy()
            data[current_mask] = np.nan
            # Edge-pad image and mask, then FFT-convolve, then crop back.
            data_pad = np.pad(data, ((pad_y, pad_y), (pad_x, pad_x)), mode='edge')
            model_pad = convolve_fft(data_pad, kernel, boundary='fill',
                                     fill_value=np.nan,
                                     nan_treatment='interpolate',
                                     allow_huge=True)
            model = model_pad[pad_y:-pad_y, pad_x:-pad_x]

            if sigma_clip <= 0:
                break

            # Sigma-clip on the valid (non-order) pixels only.
            valid = ~order_mask
            resid = image.astype(np.float64) - model
            valid_resid = resid[valid]
            med = float(np.nanmedian(valid_resid))
            mad = float(np.nanmedian(np.abs(valid_resid - med)))
            sig = 1.4826 * mad if mad > 0 else float(np.nanstd(valid_resid))
            if sig <= 0:
                break

            if clip_mode == 'upper':
                outlier = (resid - med) > sigma_clip * sig
            elif clip_mode == 'lower':
                outlier = (resid - med) < -sigma_clip * sig
            else:
                outlier = np.abs(resid - med) > sigma_clip * sig

            new_mask = order_mask | outlier
            if np.array_equal(new_mask, current_mask):
                break
            current_mask = new_mask

        return model

    @staticmethod
    def _estimate_background_column_spline(image: np.ndarray, order_mask: np.ndarray,
                                           spline_smooth_factor: float = 1.0,
                                           post_smooth_sigma_x: float = 5.0,
                                           sigma_clip: float = 3.0,
                                           maxiters: int = 3,
                                           clip_mode: str = 'upper'
                                           ) -> np.ndarray:
        """Estimate scattered-light background via column-wise 1D smoothing splines.

        For each column (or every 2nd column for speed on wide images), extract
        unmasked (inter-order) pixels, fit a cubic smoothing spline, then
        evaluate over the full column.  Skipped columns are filled by linear
        interpolation.  A final light horizontal Gaussian smooth removes any
        residual column-to-column discontinuity.

        Parameters
        ----------
        image : 2D array — science or flat image.
        order_mask : bool mask, True = inside an order (to be ignored).
        spline_smooth_factor : Multiplier for the spline smoothing parameter *s*.
            Computed as ``s = factor × N_bg × noise_variance``.  factor=1.0
            gives statistically optimal smoothing; values < 1 produce a tighter
            fit that follows ripple peaks/troughs.  Typical 0.3–2.0.
        post_smooth_sigma_x : Gaussian σ (pixels) for a light horizontal
            smoothing pass applied after assembling all columns.
            Set to 0 to disable.
        sigma_clip : iterative outlier rejection threshold (MAD-based).
        maxiters : max sigma-clipping iterations.
        clip_mode : 'upper' / 'both' / 'lower'.
        """
        from scipy.interpolate import splrep, splev
        from scipy.ndimage import gaussian_filter

        rows, cols = image.shape
        background = np.full((rows, cols), np.nan, dtype=np.float64)
        y_all = np.arange(rows, dtype=np.float64)

        # Process every step-th column; interpolate the rest.
        step = 2 if cols > 2000 else 1
        fit_cols = list(range(0, cols, step))
        if fit_cols[-1] != cols - 1:
            fit_cols.append(cols - 1)

        log_step = max(1, len(fit_cols) // 10)
        for ci, ix in enumerate(fit_cols):
            if ci % log_step == 0:
                logger.info(f"  column_spline: column {ix}/{cols}")
            col_data = image[:, ix].astype(np.float64)
            col_mask = order_mask[:, ix]
            valid = ~col_mask
            y_bg = y_all[valid]
            f_bg = col_data[valid]

            if len(y_bg) < 4:
                background[:, ix] = np.nanmedian(f_bg) if len(f_bg) > 0 else 0.0
                continue

            # Estimate per-column noise for proper s scaling.
            noise_est = float(np.median(np.abs(np.diff(f_bg)))) / 0.6745
            var_est = max(noise_est ** 2, 1.0)

            # Data range for clamping extrapolation artifacts
            f_min = float(np.min(f_bg))
            f_max = float(np.max(f_bg))
            f_range = max(f_max - f_min, 1.0)
            clamp_lo = f_min - 0.5 * f_range
            clamp_hi = f_max + 0.5 * f_range

            # Iterative sigma-clip on 1D residuals
            current_valid = np.ones(len(y_bg), dtype=bool)
            spl_tck = None
            for _it in range(max(1, maxiters)):
                yv = y_bg[current_valid]
                fv = f_bg[current_valid]
                if len(yv) < 4:
                    break
                s_val = spline_smooth_factor * len(yv) * var_est
                try:
                    spl_tck = splrep(yv, fv, s=s_val, k=3)
                except Exception:
                    spl_tck = None
                    break
                resid = f_bg - splev(y_bg, spl_tck)
                med = float(np.median(resid[current_valid]))
                mad = float(np.median(np.abs(resid[current_valid] - med)))
                sig = 1.4826 * mad if mad > 0 else float(np.std(resid[current_valid]))
                if sig <= 0 or sigma_clip <= 0:
                    break
                if clip_mode == 'upper':
                    outlier = (resid - med) > sigma_clip * sig
                elif clip_mode == 'lower':
                    outlier = (resid - med) < -sigma_clip * sig
                else:
                    outlier = np.abs(resid - med) > sigma_clip * sig
                new_valid = current_valid & ~outlier
                if np.array_equal(new_valid, current_valid):
                    break
                current_valid = new_valid

            # Final spline with cleaned data
            yv = y_bg[current_valid]
            fv = f_bg[current_valid]
            if len(yv) < 4:
                background[:, ix] = np.nanmedian(col_data[valid])
                continue
            s_val = spline_smooth_factor * len(yv) * var_est
            try:
                spl_tck = splrep(yv, fv, s=s_val, k=3)
                col_bg = splev(y_all, spl_tck)
                np.clip(col_bg, clamp_lo, clamp_hi, out=col_bg)
                background[:, ix] = col_bg
            except Exception:
                background[:, ix] = np.nanmedian(fv)

        # Interpolate skipped columns linearly
        if step > 1:
            fit_arr = np.array(fit_cols, dtype=np.float64)
            for iy in range(rows):
                row_vals = background[iy, fit_cols]
                good = np.isfinite(row_vals)
                if np.sum(good) >= 2:
                    background[iy, :] = np.interp(
                        np.arange(cols), fit_arr[good], row_vals[good])
                elif np.any(good):
                    background[iy, :] = row_vals[good][0]

        # Light horizontal (dispersion-direction) smoothing to remove
        # column-to-column discontinuities.
        if post_smooth_sigma_x > 0:
            background = gaussian_filter(background, sigma=[0, post_smooth_sigma_x])

        return background

    def _estimate_background_split_halves(self, image: np.ndarray, order_mask: np.ndarray,
                                          order: int, method: str,
                                          smooth_sigma: float, sigma_clip: float,
                                          maxiters: int, bspline_smooth: float,
                                          split_row: int = 0,
                                          clip_mode: str = 'both',
                                          kernel_sigma_x: float = 15.0,
                                          kernel_sigma_y: float = 15.0,
                                          spline_smooth_factor: float = 1.0,
                                          spline_post_smooth_x: float = 5.0,
                                          ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Estimate background separately for upper/lower detector halves.

        Args:
            split_row: 1-based detector row dividing the two amplifier halves.
                       0 (default) falls back to rows//2.
            clip_mode: 'both' (symmetric), 'upper' (reject bright outliers only),
                       or 'lower' (reject faint outliers only).
            kernel_sigma_x, kernel_sigma_y: Gaussian kernel σ for 'convolution' method.
            spline_smooth_factor: Smoothing factor multiplier for 'column_spline' method.
            spline_post_smooth_x: Horizontal Gaussian σ for post-smoothing in 'column_spline' method.
        """
        rows, cols = image.shape

        # --- Convolution method: runs on the full image (no split) ---
        if method == 'convolution':
            logger.info(f"Background convolution: Gaussian kernel σ_x={kernel_sigma_x:.1f}, "
                        f"σ_y={kernel_sigma_y:.1f}")
            background = self._estimate_background_convolution(
                image, order_mask,
                kernel_sigma_x=kernel_sigma_x, kernel_sigma_y=kernel_sigma_y,
                sigma_clip=sigma_clip, maxiters=maxiters, clip_mode=clip_mode,
            )
            # Use the configured detector split row (1-based → 0-based index).
            split = (split_row - 1) if 0 < split_row <= rows else rows // 2
            split = max(1, min(split, rows - 1))
            model_lower = background[:split, :]
            model_upper = background[split:, :]
            return background, model_lower, model_upper

        # --- Column-wise 1D spline method: runs on the full image (no split) ---
        if method == 'column_spline':
            logger.info(f"Background column_spline: smooth_factor={spline_smooth_factor:.3f}, "
                        f"post_smooth_σ_x={spline_post_smooth_x:.1f}")
            background = self._estimate_background_column_spline(
                image, order_mask,
                spline_smooth_factor=spline_smooth_factor,
                post_smooth_sigma_x=spline_post_smooth_x,
                sigma_clip=sigma_clip, maxiters=maxiters, clip_mode=clip_mode,
            )
            split = (split_row - 1) if 0 < split_row <= rows else rows // 2
            split = max(1, min(split, rows - 1))
            model_lower = background[:split, :]
            model_upper = background[split:, :]
            return background, model_lower, model_upper

        # --- Polynomial / B-spline / smooth: per-half fitting ---
        # Use the configured detector split row (1-based → 0-based index).
        split = (split_row - 1) if 0 < split_row <= rows else rows // 2
        split = max(1, min(split, rows - 1))
        background = np.zeros_like(image, dtype=np.float64)
        model_lower = np.zeros((split, cols), dtype=np.float64)
        model_upper = np.zeros((rows - split, cols), dtype=np.float64)

        halves = [
            (0, split, "lower"),
            (split, rows, "upper"),
        ]
        for y0, y1, label in halves:
            sub = image[y0:y1, :].astype(np.float64)
            sub_mask = order_mask[y0:y1, :]
            valid = ~sub_mask

            sampled = sub.copy()
            sampled[sub_mask] = np.nan

            if np.mean(valid) < 0.03:
                # Fallback in very dense order region: smooth heavily on unmasked
                # values while leaving masked zones as NaN constraints.
                logger.warning(
                    f"Inter-order fraction too small in {label} half ({np.mean(valid):.4f}); "
                    "using robust median-smoothed fallback background."
                )
                fill = np.nanmedian(sampled)
                if not np.isfinite(fill):
                    fill = float(np.nanmedian(sub))
                tmp = np.where(np.isfinite(sampled), sampled, fill)
                sub_bg = median_filter(tmp, size=31).astype(np.float64)
            else:
                sub_bg = estimate_background_region(
                    sampled,
                    order=order,
                    valid_mask=valid,
                    method=method,
                    smooth_sigma=smooth_sigma,
                    sigma_clip=sigma_clip,
                    maxiters=maxiters,
                    bspline_smooth=bspline_smooth,
                    clip_mode=clip_mode,
                ).astype(np.float64)

            background[y0:y1, :] = sub_bg
            if label == "lower":
                model_lower[:, :] = sub_bg
            else:
                model_upper[:, :] = sub_bg

        return background, model_lower, model_upper

    def estimate_background(self, image: np.ndarray, order: int = 2) -> np.ndarray:
        """
        Estimate 2D background using polynomial fitting.

        Args:
            image: 2D image array
            order: Polynomial order

        Returns:
            Background model
        """
        logger.info(f"Estimating 2D background (polyorder={order})...")

        background = estimate_background_2d(image, order=order)
        self.background_2d = background

        logger.info(f"Background estimated: min={np.min(background):.1f}, "
                   f"max={np.max(background):.1f}")

        return background

    def estimate_background_interorder(self, image: np.ndarray, apertures: ApertureSet,
                                       order: int = 2,
                                       mask_margin_pixels: int = 1,
                                       method: str = 'chebyshev',
                                       smooth_sigma: float = 20.0,
                                       sigma_clip: float = 3.0,
                                       maxiters: int = 4,
                                       bspline_smooth: float = 1.0) -> Tuple[np.ndarray, np.ndarray]:
        """
        Estimate 2D background only from inter-order pixels.

        Steps:
        1) Build aperture mask from traced orders.
        2) Keep only inter-order pixels as fit constraints.
        3) Fit smooth 2D background surface.
        """
        logger.info(
            f"Estimating inter-order background (method={method}, order={order}, margin_px={mask_margin_pixels})..."
        )

        order_mask = self._build_aperture_mask(image.shape, apertures, margin_pixels=mask_margin_pixels)
        inter_order = ~order_mask

        sampled = image.astype(np.float64).copy()
        sampled[order_mask] = np.nan

        # Use inter-order pixels directly as constraints with robust clipping.
        background = estimate_background_2d(
            sampled,
            order=order,
            valid_mask=inter_order,
            method=method,
            smooth_sigma=smooth_sigma,
            sigma_clip=sigma_clip,
            maxiters=maxiters,
            bspline_smooth=bspline_smooth,
        )
        self.background_2d = background

        logger.info(
            f"Inter-order background estimated ({method}): min={np.min(background):.1f}, max={np.max(background):.1f}, "
            f"inter-order fraction={np.mean(inter_order):.3f}"
        )
        return background, order_mask

    def apply_background_correction(self, image: np.ndarray,
                                   background: Optional[np.ndarray] = None) -> np.ndarray:
        """
        Apply background correction to image.

        Args:
            image: Original science image
            background: Background model (uses self.background_2d if None)

        Returns:
            Background-corrected image
        """
        bkg = background if background is not None else self.background_2d

        if bkg is None:
            logger.warning("No background model available, returning uncorrected image")
            return image

        if image.shape != bkg.shape:
            raise ValueError(f"Image shape {image.shape} does not match "
                           f"background shape {bkg.shape}")

        corrected = image.astype(float) - bkg.astype(float)
        logger.debug(f"Background corrected: min={np.min(corrected):.1f}, "
                    f"max={np.max(corrected):.1f}")

        return corrected

    def extract_interorder_background(self, image: np.ndarray,
                                     aperture_positions: np.ndarray) -> np.ndarray:
        """
        Extract inter-order background by masking order regions.

        Args:
            image: 2D science image
            aperture_positions: Aperture center positions

        Returns:
            Inter-order background map
        """
        logger.info("Extracting inter-order background...")

        # Create aperture mask
        rows, cols = image.shape
        aperture_mask = np.zeros((rows, cols), dtype=bool)

        # Mark pixels within orders
        for x in range(cols):
            for aperture in aperture_positions:
                center = aperture['center'] + aperture.get('center_slope', 0) * x
                width = aperture.get('width', 20)
                y_min = max(0, int(center - width/2))
                y_max = min(rows, int(center + width/2))
                aperture_mask[y_min:y_max, x] = True

        # Image outside orders is inter-order light
        interorder_pixels = image[~aperture_mask]

        if len(interorder_pixels) > 0:
            median_bkg = np.median(interorder_pixels)
            logger.info(f"Median inter-order background: {median_bkg:.1f}")

            return median_bkg
        else:
            return 0.0

    def save_background(self, output_path: str):
        """Save background model to FITS file."""
        if self.background_2d is None:
            raise RuntimeError("No background model to save")

        from astropy.io import fits
        from pathlib import Path

        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        hdu = fits.PrimaryHDU(data=self.background_2d.astype(np.float32))
        hdu.header['EXTNAME'] = '2D BACKGROUND'
        hdu.header['ORIGIN'] = 'SpecProc'
        hdul = fits.HDUList([hdu])
        hdul.writeto(str(output_path), overwrite=True)

        logger.info(f"Saved background model to {output_path}")


def _extract_boundaries_from_labels(order_labels):
    """Extract per-order per-column lower/upper boundaries from a labeled mask.

    Parameters
    ----------
    order_labels : np.ndarray[int16], shape (rows, cols)
        Pixel value = order ID, 0 = background.

    Returns
    -------
    order_ids : list[int]
        Sorted order IDs (ascending).
    y_lo : np.ndarray, shape (n_orders, cols)
        Lower boundary row per order per column (NaN where absent).
    y_hi : np.ndarray, shape (n_orders, cols)
        Upper boundary row per order per column (NaN where absent).
    """
    rows, cols = order_labels.shape
    unique_ids = sorted(int(v) for v in np.unique(order_labels) if v != 0)
    n = len(unique_ids)
    y_lo = np.full((n, cols), np.nan, dtype=np.float64)
    y_hi = np.full((n, cols), np.nan, dtype=np.float64)

    for i, oid in enumerate(unique_ids):
        order_mask_2d = (order_labels == oid)
        for xi in range(cols):
            col_rows = np.where(order_mask_2d[:, xi])[0]
            if col_rows.size > 0:
                y_lo[i, xi] = float(col_rows[0])
                y_hi[i, xi] = float(col_rows[-1])

    return unique_ids, y_lo, y_hi


def _extrapolate_virtual_boundaries(order_ids, y_lo, y_hi, n_below, n_above, rows,
                                    n_ref=8):
    """Add virtual order boundaries at detector edges by spacing extrapolation.

    Uses the nearest *n_ref* real orders to fit a low-degree polynomial trend
    of inter-order spacing and aperture width, then extrapolates outward.
    Degree is automatically chosen: quadratic when ≥4 data points, linear
    when 2–3 points, constant when only 1 point.  This captures the nonlinear
    spacing variation dictated by the grating equation (∝ 1/m²).

    Returns
    -------
    all_ids : list[int]
    y_lo_ext : np.ndarray, shape (n_orders + n_below + n_above, cols)
    y_hi_ext : np.ndarray, shape (n_orders + n_below + n_above, cols)
    """
    n_orders = len(order_ids)
    cols = y_lo.shape[1]

    # Per-column center and width
    centers = 0.5 * (y_lo + y_hi)   # (n_orders, cols)
    widths = y_hi - y_lo            # (n_orders, cols)

    def _polyfit_cols(j_indices, data_2d, max_degree=2):
        """Polynomial fit per column.  Returns coefficient array (deg+1, cols).

        Automatically reduces degree if there are too few data points.
        Coefficients are in np.polyval order (highest power first).
        """
        j = np.asarray(j_indices, dtype=float)
        N = len(j)
        deg = min(max_degree, N - 1)
        if deg < 0:
            return np.zeros((1, data_2d.shape[1]))
        # Fit each column independently
        coefs = np.empty((deg + 1, data_2d.shape[1]), dtype=np.float64)
        for ci in range(data_2d.shape[1]):
            coefs[:, ci] = np.polyfit(j, data_2d[:, ci], deg)
        return coefs

    def _polyeval_cols(coefs, j_val):
        """Evaluate polynomial coefs (deg+1, cols) at scalar j_val → (cols,)."""
        deg = coefs.shape[0] - 1
        result = coefs[0].copy()
        for p in range(1, deg + 1):
            result = result * j_val + coefs[p]
        return result

    # --- below (extrapolate below the first order) ---
    below_ids = []
    below_lo = []
    below_hi = []
    if n_below > 0 and n_orders >= 2:
        min_id = min(order_ids)
        nr = min(n_ref, n_orders)

        # Inter-order spacings of the nr bottom orders
        sp_2d = np.diff(centers[:nr], axis=0)      # (nr-1, cols)
        j_sp = np.arange(sp_2d.shape[0], dtype=float)
        sp_coefs = _polyfit_cols(j_sp, sp_2d)

        # Widths of the nr bottom orders
        w_2d = widths[:nr]                          # (nr, cols)
        j_w = np.arange(nr, dtype=float)
        w_coefs = _polyfit_cols(j_w, w_2d)

        c_prev = centers[0].copy()

        for k in range(1, n_below + 1):
            sp_k = _polyeval_cols(sp_coefs, -k)
            sp_k = np.maximum(sp_k, 3.0)
            vc = c_prev - sp_k
            wk = _polyeval_cols(w_coefs, -k)
            wk = np.maximum(wk, 3.0)
            vlo = vc - 0.5 * wk
            vhi = vc + 0.5 * wk

            vid = min_id - k
            below_ids.insert(0, vid)
            below_lo.insert(0, vlo)
            below_hi.insert(0, vhi)
            c_prev = vc

    # --- above (extrapolate above the last order) ---
    above_ids = []
    above_lo = []
    above_hi = []
    if n_above > 0 and n_orders >= 2:
        max_id = max(order_ids)
        nr = min(n_ref, n_orders)

        sp_2d = np.diff(centers[-nr:], axis=0)      # (nr-1, cols)
        j_sp = np.arange(sp_2d.shape[0], dtype=float)
        sp_coefs = _polyfit_cols(j_sp, sp_2d)

        w_2d = widths[-nr:]
        j_w = np.arange(nr, dtype=float)
        w_coefs = _polyfit_cols(j_w, w_2d)

        last_sp_j = float(sp_2d.shape[0] - 1)
        c_prev = centers[-1].copy()

        for k in range(1, n_above + 1):
            sp_k = _polyeval_cols(sp_coefs, last_sp_j + k)
            sp_k = np.maximum(sp_k, 3.0)
            vc = c_prev + sp_k
            wk = _polyeval_cols(w_coefs, float(nr - 1) + k)
            wk = np.maximum(wk, 3.0)
            vlo = vc - 0.5 * wk
            vhi = vc + 0.5 * wk

            vid = max_id + k
            above_ids.append(vid)
            above_lo.append(vlo)
            above_hi.append(vhi)
            c_prev = vc

    # Concatenate
    all_ids = below_ids + list(order_ids) + above_ids
    parts_lo = below_lo + [y_lo[i] for i in range(n_orders)] + above_lo
    parts_hi = below_hi + [y_hi[i] for i in range(n_orders)] + above_hi
    y_lo_ext = np.vstack(parts_lo) if parts_lo else y_lo
    y_hi_ext = np.vstack(parts_hi) if parts_hi else y_hi

    n_virt = len(below_ids) + len(above_ids)
    if n_virt > 0:
        logger.info(f"Virtual orders: {len(below_ids)} below ({below_ids}), "
                    f"{len(above_ids)} above ({above_ids})")

    return all_ids, y_lo_ext, y_hi_ext


def build_mask_from_labels(output_dir_base: str,
                           flat_data=None,
                           n_mask_below: int = 2,
                           n_mask_above: int = 1,
                           mask_margin_pixels: int = 1):
    """Build background boolean mask from the Step 2 labeled order mask.

    Reads ``Orders_mask.fits`` from step2_trace, optionally adds virtual
    orders at the edges, expands boundaries by *mask_margin_pixels*, and
    resolves overlaps using the valley minimum in the flat cross-section.

    Parameters
    ----------
    output_dir_base : str
        Base directory containing step2_trace/Orders_mask.fits
    flat_data : np.ndarray or None
        Master flat image used for valley detection.  If None, overlapping
        orders are split at the midpoint.
    n_mask_below : int
        Number of virtual orders to extend below hardware orders.
    n_mask_above : int
        Number of virtual orders to extend above hardware orders.
    mask_margin_pixels : int
        Pixels to expand the order boundaries by before background fitting.

    Returns
    -------
    mask : np.ndarray[bool]  — True = inside order (to be masked).
    order_ids : list[int]    — All order IDs in the mask (real + virtual).
    y_lo, y_hi : np.ndarray  — Per-order per-column boundaries after expansion.
    """
    from src.utils.fits_io import read_fits_image

    base = Path(output_dir_base) / 'step2_trace'
    labels_path = base / 'Orders_mask.fits'
    if not labels_path.exists():
        logger.error(f"Orders_mask.fits not found at {labels_path}; run Step 2 first")
        return None, [], None, None

    order_labels, _ = read_fits_image(str(labels_path))
    order_labels = order_labels.astype(np.int32)
    rows, cols = order_labels.shape
    logger.info(f"Loaded labeled order mask ({rows}×{cols}) from {labels_path.name}")

    # Extract per-order boundaries
    order_ids, y_lo, y_hi = _extract_boundaries_from_labels(order_labels)
    logger.info(f"Extracted boundaries for {len(order_ids)} orders "
                f"(IDs {min(order_ids)}..{max(order_ids)})")

    # Add virtual orders at edges
    order_ids, y_lo, y_hi = _extrapolate_virtual_boundaries(
        order_ids, y_lo, y_hi, n_mask_below, n_mask_above, rows)

    # Expand by margin
    margin = mask_margin_pixels
    n_orders = len(order_ids)

    y_lo_exp = y_lo - margin
    y_hi_exp = y_hi + margin

    # Resolve overlaps between adjacent orders (including virtual ones).
    # For each adjacent pair, if the expanded gap < 1 pixel, find the valley
    # (minimum intensity in the flat cross-section between the two orders)
    # and split there, keeping at least 1 unmasked background pixel at the valley.
    if n_orders >= 2:
        n_resolved = 0
        for i in range(n_orders - 1):
            for xi in range(cols):
                gap = y_lo_exp[i + 1, xi] - y_hi_exp[i, xi]
                if gap >= 1.0:
                    continue
                # The expanded masks overlap or nearly touch.
                # Raw (un-expanded) boundaries define the search region.
                raw_hi = y_hi[i, xi]
                raw_lo_next = y_lo[i + 1, xi]
                if not (np.isfinite(raw_hi) and np.isfinite(raw_lo_next)):
                    mid = int(np.round(0.5 * (y_hi_exp[i, xi] + y_lo_exp[i + 1, xi])))
                    y_hi_exp[i, xi] = mid - 1
                    y_lo_exp[i + 1, xi] = mid + 1
                    continue

                r0 = max(0, int(np.ceil(raw_hi)))
                r1 = min(rows, int(np.floor(raw_lo_next)) + 1)
                if flat_data is not None and r1 > r0 + 1:
                    seg = flat_data[r0:r1, xi].astype(np.float64)
                    valley = r0 + int(np.argmin(seg))
                else:
                    valley = int(np.round(0.5 * (raw_hi + raw_lo_next)))

                # Keep the valley pixel unmasked: upper mask stops 1 px below,
                # lower mask starts 1 px above.
                y_hi_exp[i, xi] = valley - 1
                y_lo_exp[i + 1, xi] = valley + 1
                n_resolved += 1

        if n_resolved > 0:
            logger.info(f"Resolved {n_resolved} overlapping margin positions "
                        f"using valley detection (margin={margin})")

    # Clip to detector and build boolean mask
    y_lo_int = np.clip(np.floor(y_lo_exp), 0, rows).astype(int)
    y_hi_int = np.clip(np.ceil(y_hi_exp) + 1, 0, rows).astype(int)

    mask = np.zeros((rows, cols), dtype=bool)
    for i in range(n_orders):
        for xi in range(cols):
            lo = y_lo_int[i, xi]
            hi = y_hi_int[i, xi]
            if hi > lo:
                mask[lo:hi, xi] = True

    bg_frac = np.mean(~mask)
    if bg_frac < 0.05:
        logger.warning(f"Order mask leaves only {bg_frac*100:.1f}% background pixels "
                       f"(margin={margin}). Consider reducing mask_margin_pixels.")
    logger.info(f"Background mask: {np.mean(mask)*100:.1f}% masked, "
                f"{bg_frac*100:.1f}% background")

    return mask, order_ids, y_lo_exp, y_hi_exp


def _build_envelope_aware_mask(image_shape, apertures, mask_margin_pixels=0,
                               fill_below_edge=False, fill_above_edge=False):
    """Build an aperture mask using traced profile-root boundaries plus an optional margin.

    The margin expands each order's mask beyond its profile-root edges.  To prevent
    adjacent masks from completely filling the inter-order gap (leaving no background
    pixels for fitting), the margin is clamped **per-column per-pair** so that at
    least 1 pixel of background always remains between every neighbouring pair.

    Parameters
    ----------
    fill_below_edge : bool
        If True, mask from row 0 up to the lowest order's lower boundary
        at each column.  Used when virtual orders extend below the detector.
    fill_above_edge : bool
        If True, mask from the highest order's upper boundary to the last
        row at each column.  Used when virtual orders extend above the detector.

    Returns
    -------
    mask : np.ndarray[bool]  – True for pixels *inside* orders (to be masked).
    """
    rows, cols = image_shape
    mask = np.zeros((rows, cols), dtype=bool)

    # Sort by physical position (center y at mid-column) so that pair clamping
    # operates on genuinely adjacent orders.  Sorting by aperture ID breaks
    # when virtual orders (negative IDs) are present: top-of-detector and
    # bottom-of-detector virtuals become interleaved in ID order.
    mid_col = cols // 2
    ordered_apertures = sorted(
        apertures.apertures.items(),
        key=lambda x: float(x[1].get_position(mid_col)))
    x_coords = np.arange(cols, dtype=float)

    n_ap = len(ordered_apertures)
    # Pre-evaluate per-aperture lower/upper arrays once (shape: n_ap × cols).
    all_lower = np.empty((n_ap, cols), dtype=np.float64)
    all_upper = np.empty((n_ap, cols), dtype=np.float64)
    all_center = np.empty((n_ap, cols), dtype=np.float64)
    for i, (_aid, ap) in enumerate(ordered_apertures):
        all_center[i] = ap.get_position(x_coords)
        all_lower[i] = ap.get_lower(x_coords)
        all_upper[i] = ap.get_upper(x_coords)

    margin = mask_margin_pixels

    # Compute effective y1 (lower edge - margin) and y2 (upper edge + margin + 1)
    # for each order at each column, then clamp adjacent pairs.
    # y1[i, xi] = floor(lower[i, xi] - margin)
    # y2[i, xi] = ceil(upper[i, xi] + margin + 1)
    y1_arr = np.floor(np.minimum(all_lower, all_upper) - margin).astype(np.float64)
    y2_arr = np.ceil(np.maximum(all_lower, all_upper) + margin + 1).astype(np.float64)

    # Per-column per-pair clamping: ensure at least 1 pixel of background between
    # adjacent orders.  If upper-edge+margin of order i >= lower-edge-margin of
    # order i+1, split the gap at the midpoint minus 0.5 / plus 0.5 so 1 pixel
    # remains.
    if margin > 0 and n_ap >= 2:
        n_clamped = 0
        for i in range(n_ap - 1):
            # Overlap columns: y2 of order i > y1 of order i+1 - 1
            overlap = y2_arr[i] > y1_arr[i + 1] - 1
            if not np.any(overlap):
                continue
            # Midpoint between the raw (un-margined) upper of i and lower of i+1
            mid = 0.5 * (all_upper[i] + all_lower[i + 1])
            # Clamp: order i's mask stops at floor(mid), order i+1 starts at ceil(mid)+1
            y2_arr[i, overlap] = np.floor(mid[overlap])
            y1_arr[i + 1, overlap] = np.ceil(mid[overlap]) + 1
            n_clamped += int(np.sum(overlap))
        if n_clamped > 0:
            logger.info(f"Mask margin clamped at {n_clamped} column-pair positions "
                        f"to preserve background pixels (margin={margin})")

    # Clip to detector bounds
    y1_arr = np.clip(y1_arr, 0, rows).astype(int)
    y2_arr = np.clip(y2_arr, 0, rows).astype(int)

    # Paint the mask — use clipped boundaries; no need to check center position.
    for i in range(n_ap):
        for xi in range(cols):
            lo = y1_arr[i, xi]
            hi = y2_arr[i, xi]
            if not np.isfinite(all_lower[i, xi]) or not np.isfinite(all_upper[i, xi]):
                cy = float(all_center[i, xi])
                if not np.isfinite(cy):
                    continue
                lo = max(0, int(cy - 4 - margin))
                hi = min(rows, int(cy + 4 + margin + 1))
            if hi > lo:
                mask[lo:hi, xi] = True

    # Edge fill: mask detector regions beyond the outermost virtual orders.
    # When virtual orders are present near the detector edge, we fill the
    # strip between the detector edge and the outermost virtual order's
    # boundary so the background fitter ignores edge artefacts.
    if fill_below_edge or fill_above_edge:
        if fill_below_edge:
            # Use the physically lowest aperture's lower boundary (this is
            # typically the outermost virtual order below the detector edge).
            lo_edge = y1_arr[0]  # lowest aperture in physical order
            n_filled = 0
            for xi in range(cols):
                if lo_edge[xi] > 0:
                    mask[0:lo_edge[xi], xi] = True
                    n_filled += 1
            logger.info(f"Edge fill (below): masked rows 0..{int(np.median(lo_edge))} "
                        f"(median) below outermost order ({n_filled}/{cols} columns)")
        if fill_above_edge:
            # Use the physically highest aperture's upper boundary.
            hi_edge = y2_arr[-1]  # highest aperture in physical order
            n_filled = 0
            for xi in range(cols):
                if hi_edge[xi] < rows:
                    mask[hi_edge[xi]:rows, xi] = True
                    n_filled += 1
            logger.info(f"Edge fill (above): masked rows {int(np.median(hi_edge))}..{rows-1} "
                        f"(median) above outermost order ({n_filled}/{cols} columns)")

    bg_frac = np.mean(~mask)
    if bg_frac < 0.05:
        logger.warning(f"Order mask leaves only {bg_frac*100:.1f}% background pixels "
                       f"(margin={margin}). Consider reducing mask_margin_pixels.")

    return mask


def build_and_save_order_mask(output_dir_base: str,
                              image_shape: tuple,
                              apertures: ApertureSet = None,
                              output_subdir: str = 'step3_scatterlight',
                              flat_data: np.ndarray = None,
                              n_mask_below: int = 2,
                              n_mask_above: int = 1,
                              mask_margin_pixels: int = 1,
                              save_plots: bool = True,
                              fig_format: str = 'png') -> Optional[np.ndarray]:
    """Build the MasterFlat order mask once and save it to disk.

    Uses the labeled order mask (Orders_mask.fits) produced by Step 2.
    Virtual orders are extrapolated at the edges, boundaries are expanded
    by *mask_margin_pixels*, and overlapping margins are resolved using the
    valley minimum of the flat cross-section.

    Parameters
    ----------
    flat_data : np.ndarray or None
        Master flat image for valley detection when margins overlap.

    Returns
    -------
    order_mask : np.ndarray[bool] or None
        True = inside an order (to be masked).  None on failure.
    """
    result = build_mask_from_labels(
        output_dir_base, flat_data=flat_data,
        n_mask_below=n_mask_below, n_mask_above=n_mask_above,
        mask_margin_pixels=mask_margin_pixels
    )
    if result[0] is None:
        # Fallback: if Orders_mask.fits not available, try legacy aperture-based mask.
        if apertures is not None and apertures.norders > 0:
            logger.warning("Orders_mask.fits not found, falling back to aperture-based mask")
            order_mask = _build_envelope_aware_mask(
                image_shape, apertures,
                mask_margin_pixels=mask_margin_pixels,
                fill_below_edge=(n_mask_below == 0),
                fill_above_edge=(n_mask_above == 0))
        else:
            logger.warning("build_and_save_order_mask: no labels or apertures available")
            return None
    else:
        order_mask = result[0]
        all_order_ids = result[1]
        y_lo_exp = result[2]
        y_hi_exp = result[3]

    out_dir = Path(output_dir_base) / output_subdir
    out_dir.mkdir(parents=True, exist_ok=True)

    # FITS
    mask_fits = out_dir / 'MasterFlat_orders_mask.fits'
    write_fits_image(str(mask_fits), order_mask.astype(np.int16), dtype='int16')
    logger.info(f"Saved MasterFlat order mask to {mask_fits}")

    # PNG diagnostic
    if save_plots:
        rows, cols = image_shape if order_mask.ndim == 2 else order_mask.shape
        rows, cols = order_mask.shape
        fig, ax = plt.subplots(figsize=(12, 10))
        ax.imshow(order_mask.astype(float), aspect='auto', origin='lower', cmap='gray')
        # Add order number labels on the right side
        x_label = cols - 1
        if result[0] is not None and y_lo_exp is not None:
            # Use label-based boundaries for label positions
            for idx, oid in enumerate(all_order_ids):
                yc = 0.5 * (y_lo_exp[idx, x_label] + y_hi_exp[idx, x_label])
                if np.isfinite(yc) and 0 <= yc < rows:
                    ax.text(x_label + cols * 0.005, yc, str(oid),
                            color='black', fontsize=3.5, ha='left', va='center',
                            clip_on=False)
        elif apertures is not None:
            for ap_id in sorted(apertures.get_orders()):
                ap = apertures.get_aperture(ap_id)
                yc = float(ap.get_position(x_label))
                if 0 <= yc < rows:
                    ax.text(x_label + cols * 0.005, yc, str(ap_id),
                            color='black', fontsize=3.5, ha='left', va='center',
                            clip_on=False)
        ax.set_xlim(0, cols - 1)
        ax.set_title(
            f'MasterFlat Order Mask  (label-based boundaries ± {mask_margin_pixels} px margin)\n'
            f'1=masked order region, 0=inter-order background'
        )
        ax.set_xlabel('Pixel (X)')
        ax.set_ylabel('Pixel (Y)')
        fig.colorbar(ax.images[0], ax=ax, label='Mask (1=order, 0=inter-order)')
        fig.tight_layout()
        fig.savefig(str(out_dir / f'MasterFlat_orders_mask.{fig_format}'), dpi=150, bbox_inches='tight')
        plt.close(fig)
        logger.info(f"Saved MasterFlat order mask plot to {out_dir / f'MasterFlat_orders_mask.{fig_format}'})")

    return order_mask


def process_background_stage(science_image: np.ndarray,
                             output_dir_base: str,
                             output_subdir: str = 'step3_scatterlight',
                             output_tag: str = '',
                             apertures: Optional[ApertureSet] = None,
                             mask_margin_scale: float = 1.2,
                             sci_header=None,
                             order_mask: Optional[np.ndarray] = None,
                             poly_order: int = 3,
                             bg_method: str = 'convolution',
                             smooth_sigma: float = 20.0,
                             sigma_clip_val: float = 3.0,
                             maxiters: int = 4,
                             bspline_smooth: float = 1.0,
                             mask_margin_pixels: int = 1,
                             n_mask_below: int = 2,
                             n_mask_above: int = 1,
                             clip_mode: str = 'upper',
                             split_row: int = 2068,
                             kernel_sigma_x: float = 40.0,
                             kernel_sigma_y: float = 15.0,
                             spline_smooth_factor: float = 1.0,
                             spline_post_smooth_x: float = 5.0,
                             save_plots: bool = True,
                             fig_format: str = 'png') -> np.ndarray:
    """
    Execute scattered light subtraction stage (Step 3).

    The order mask (MasterFlat_orders_mask) should be built once via
    build_and_save_order_mask() and passed in via the ``order_mask`` parameter.
    If omitted, the mask is built on-the-fly from ``apertures`` (legacy path,
    discouraged because it re-saves the mask file for every image).

    Returns:
        Background model (2D array, same shape as science_image)
    """
    remover = BackgroundRemover()

    if bg_method == '2d_poly':
        bg_method = 'chebyshev'
    elif bg_method == 'smooth':
        bg_method = 'gaussian_smooth'

    # Use the pre-built mask when provided; otherwise build on-the-fly (legacy).
    if order_mask is not None:
        # Mask was built once externally from the master flat — just use it.
        logger.info("Step 3: using pre-built MasterFlat order mask for scattered light fitting.")
        sampled_background = science_image.astype(np.float64).copy()
        sampled_background[order_mask] = np.nan

        background, model_lower, model_upper = remover._estimate_background_split_halves(
            science_image, order_mask,
            order=poly_order, method=bg_method,
            smooth_sigma=smooth_sigma, sigma_clip=sigma_clip_val,
            maxiters=maxiters, bspline_smooth=bspline_smooth,
            split_row=split_row, clip_mode=clip_mode,
            kernel_sigma_x=kernel_sigma_x, kernel_sigma_y=kernel_sigma_y,
            spline_smooth_factor=spline_smooth_factor,
            spline_post_smooth_x=spline_post_smooth_x,
        )
        remover.background_2d = background
        logger.info(
            f"Split-half scattered light estimated ({bg_method}): "
            f"min={np.nanmin(background):.1f}, max={np.nanmax(background):.1f}, "
            f"inter-order fraction={np.mean(~order_mask):.3f}"
        )
    elif apertures is not None and apertures.norders > 0:
        # Legacy path: build and save mask here (one image at a time).
        logger.info("Step 3: building envelope-aware aperture mask on-the-fly (prefer "
                    "calling build_and_save_order_mask() once before the image loop).")
        order_mask = _build_envelope_aware_mask(
            science_image.shape, apertures,
            mask_margin_pixels=mask_margin_pixels,
            fill_below_edge=(n_mask_below == 0),
            fill_above_edge=(n_mask_above == 0)
        )
        sampled_background = science_image.astype(np.float64).copy()
        sampled_background[order_mask] = np.nan

        background, model_lower, model_upper = remover._estimate_background_split_halves(
            science_image, order_mask,
            order=poly_order, method=bg_method,
            smooth_sigma=smooth_sigma, sigma_clip=sigma_clip_val,
            maxiters=maxiters, bspline_smooth=bspline_smooth,
            split_row=split_row, clip_mode=clip_mode,
            kernel_sigma_x=kernel_sigma_x, kernel_sigma_y=kernel_sigma_y,
            spline_smooth_factor=spline_smooth_factor,
            spline_post_smooth_x=spline_post_smooth_x,
        )
        remover.background_2d = background
        logger.info(
            f"Split-half scattered light estimated ({bg_method}): "
            f"min={np.nanmin(background):.1f}, max={np.nanmax(background):.1f}, "
            f"inter-order fraction={np.mean(~order_mask):.3f}"
        )
    else:
        logger.warning("No apertures available in Step 3; fitting without order mask.")
        order_mask = np.zeros_like(science_image, dtype=bool)
        sampled_background = science_image.astype(np.float64).copy()
        background, model_lower, model_upper = remover._estimate_background_split_halves(
            science_image, order_mask,
            order=poly_order, method=bg_method,
            smooth_sigma=smooth_sigma, sigma_clip=sigma_clip_val,
            maxiters=maxiters, bspline_smooth=bspline_smooth,
            split_row=split_row, clip_mode=clip_mode,
            kernel_sigma_x=kernel_sigma_x, kernel_sigma_y=kernel_sigma_y,
            spline_smooth_factor=spline_smooth_factor,
            spline_post_smooth_x=spline_post_smooth_x,
        )
        remover.background_2d = background

    # ------------------------------------------------------------------
    # Save outputs
    # ------------------------------------------------------------------
    tag = output_tag or ''
    out_dir = Path(output_dir_base) / output_subdir
    out_dir.mkdir(parents=True, exist_ok=True)

    # Derive the science-image stem from tag (tag = '{sci_name}_' by convention).
    sci_stem = tag.rstrip('_') if tag else 'output'

    # Compute residual now so all three arrays are ready for the combined plot.
    residual = sampled_background.astype(np.float64) - background.astype(np.float64)
    residual[order_mask] = np.nan

    # Combined 3-panel diagnostic: bkg_orders_masked | bkg_model | bkg_residual
    # All three panels share the same y-axis and a single shared colorbar whose
    # scale is determined by the 1-99th percentile of the data/model arrays
    # (panels 1 & 2), while the residual uses a symmetric scale centred on 0.
    if save_plots:
        data_arr  = np.array(sampled_background, dtype=float)
        model_arr = np.array(background,         dtype=float)
        resid_arr = np.array(residual,           dtype=float)

        # data_arr and resid_arr already have NaN in masked order regions.
        # model_arr is intentionally NOT masked so the full background surface
        # (including under orders) is visible, showing overall variation.

        # Shared vmin/vmax for panels 1 & 2 (data and model on the same scale)
        finite_dm = np.isfinite(data_arr) | np.isfinite(model_arr)
        if np.any(finite_dm):
            vmin_dm = float(np.nanpercentile(np.where(np.isfinite(data_arr),  data_arr,  np.nan), 1))
            vmax_dm = float(np.nanpercentile(np.where(np.isfinite(data_arr),  data_arr,  np.nan), 99))
        else:
            vmin_dm, vmax_dm = 0.0, 1.0

        # Symmetric scale for residual panel (centred on 0)
        finite_r = np.isfinite(resid_arr)
        if np.any(finite_r):
            rlim = float(np.nanpercentile(np.abs(resid_arr[finite_r]), 99))
        else:
            rlim = 1.0
        vmin_r, vmax_r = -rlim, rlim

        rows, cols = data_arr.shape
        # Width ratio: 3 image panels + 2 narrow colorbars (one shared, one residual)
        fig, axes = plt.subplots(
            1, 3,
            figsize=(18, 6),
            sharey=True,
        )
        fig.subplots_adjust(wspace=0.04)

        titles = [
            'Scattered Light Data',
            'Scattered Light Model',
            'Residual (data − model)',
        ]
        arrays  = [data_arr,  model_arr, resid_arr]
        vmins   = [vmin_dm,   vmin_dm,   vmin_r]
        vmaxs   = [vmax_dm,   vmax_dm,   vmax_r]
        cmap_names = ['viridis', 'viridis', 'seismic']

        # Build colormap copies: panels 1 & 3 have NaN (masked orders) shown as
        # white (blank); panel 2 (model) has no NaN, keep default.
        cmaps_used = []
        for i, cn in enumerate(cmap_names):
            cm = mcm.get_cmap(cn).copy()
            if i != 1:  # panels 0 (data) and 2 (residual)
                cm.set_bad(color='white')
            cmaps_used.append(cm)

        ims = []
        for ax, arr, title, vmin, vmax, cmap in zip(
                axes, arrays, titles, vmins, vmaxs, cmaps_used):
            im = ax.imshow(arr, origin='lower', aspect='auto',
                           cmap=cmap, vmin=vmin, vmax=vmax)
            ax.set_title(title, fontsize=10)
            ax.set_xlabel('Pixel (X)')
            ims.append(im)

        axes[0].set_ylabel('Pixel (Y)')
        for ax in axes[1:]:
            ax.tick_params(labelleft=False)

        # Shared colorbar for panels 1 & 2 (attached to panel 2)
        fig.colorbar(ims[1], ax=axes[1], fraction=0.035, pad=0.02,
                     label='Counts (shared scale)')
        # Separate colorbar for residual panel
        fig.colorbar(ims[2], ax=axes[2], fraction=0.035, pad=0.02,
                     label='Residual (counts)')

        fig.suptitle('Background Diagnostics (orders masked)', fontsize=11)
        plt.savefig(str(out_dir / f'{sci_stem}_bkg_diagnostic.{fig_format}'),
                    dpi=150, bbox_inches='tight')
        plt.close(fig)
        logger.info(f"Saved background diagnostic panel to "
                    f"{out_dir / f'{sci_stem}_bkg_diagnostic.{fig_format}'}")

    # 4) Corrected image FITS only (no PNG — full science image, not a diagnostic)
    corrected = (science_image.astype(np.float64) - background.astype(np.float64)).astype(np.float32)
    corrected_fits = out_dir / f'{sci_stem}.fits'
    write_fits_image(str(corrected_fits), corrected, header=sci_header, dtype='float32')

    return background
