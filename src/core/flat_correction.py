"""
2D flat-field correction stage — Step 4 of the spectral reduction pipeline.

Applies the pixel-to-pixel flat correction map (built in Step 2 / order_tracing)
to a science image, and persists all flat-model artefacts to the step4 output
directory.

Main entry point: process_flat_correction_stage()
"""

import logging
import pickle
from pathlib import Path
from typing import Optional

import numpy as np



from scipy import ndimage
from scipy.signal import savgol_filter
from scipy.interpolate import UnivariateSpline
from scipy.ndimage import gaussian_filter1d
from typing import Dict, Tuple, List
from src.core.data_structures import FlatField, ApertureSet
from src.utils.fits_io import write_fits_image
from src.plotting.spectra_plotter import plot_2d_image_to_file

logger = logging.getLogger(__name__)


def save_flat_correction_products(flat_field: FlatField,
                                  apertures: Optional[ApertureSet],
                                  out_dir: Path,
                                  processor=None,
                                  fig_format: str = 'png',
                                  save_plots: bool = True) -> None:
    """
    Persist Step-4 flat-model artefacts (blaze profiles, sensitivity map,
    pixel flat, 2D correction map, response map, diagnostic plots).

    Parameters
    ----------
    processor : FlatCorrectionModelBuilder, optional
        If provided, re-use its already-computed response_map and diagnostics
        instead of rebuilding from scratch (saves ~5 s).
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    if flat_field.blaze_profiles:
        with open(out_dir / 'blaze_profiles.pkl', 'wb') as f:
            pickle.dump(flat_field.blaze_profiles, f)

    if flat_field.cross_profiles:
        with open(out_dir / 'cross_profiles.pkl', 'wb') as f:
            pickle.dump(flat_field.cross_profiles, f)

    # If processor is provided, save its diagnostics directly.
    if processor is not None:
        # Compute clean_flat for diagnostics
        clean_flat = flat_field.flat_data
        if flat_field.scattered_light is not None:
            clean_flat = np.clip(
                flat_field.flat_data.astype(np.float32) - flat_field.scattered_light.astype(np.float32),
                1e-6,
                None,
            )

        clean_flat_corrected = None
        flat_corr = flat_field.pixel_flat if flat_field.pixel_flat is not None else flat_field.flat_corr_2d
        if flat_corr is not None:
            safe_flat = flat_corr.astype(np.float32)
            bad = (~np.isfinite(safe_flat)) | (safe_flat <= 0.05)
            if np.any(bad):
                safe_flat = safe_flat.copy()
                safe_flat[bad] = 1.0
            clean_flat_corrected = clean_flat.astype(np.float32) / safe_flat

        processor.save_step4_diagnostics(
            out_dir,
            clean_flat=clean_flat.astype(np.float32),
            clean_flat_corrected=clean_flat_corrected,
            fig_format=fig_format,
            save_plots=save_plots
        )



class FlatCorrectionModelBuilder:
    def __init__(self):
        self.flat_data = None
        self.flat_mask = None
        self.response_map = None
        self.smoothed_model = None
        self.pixel_flat = None
        self.illumination_flat = None
        self.flat_corr_2d = None
        self.flat_sens = None
        self.order_diagnostics = {}

    def _smooth_order_response(self, flux_1d: np.ndarray, method: str = 'savgol',
                               window: int = 21, bspline_smooth: float = 0.5) -> np.ndarray:
        """Smooth per-order 1D blaze extraction with selectable method."""
        x = np.arange(flux_1d.size)
        y = np.asarray(flux_1d, dtype=float)

        valid = np.isfinite(y) & (y > 0)
        if np.sum(valid) < max(20, int(0.2 * y.size)):
            med = np.nanmedian(y[valid]) if np.any(valid) else 1.0
            return np.full_like(y, med if np.isfinite(med) and med > 0 else 1.0)

        y_interp = np.interp(x, x[valid], y[valid])
        win = max(5, int(window))
        win = min(win, y.size - 1 if y.size % 2 == 0 else y.size)
        if win < 7:
            win = 7
        if win % 2 == 0:
            win += 1

        m = (method or 'savgol').strip().lower()
        if m == 'median':
            smooth = ndimage.median_filter(y_interp, size=win, mode='nearest')
            smooth = savgol_filter(smooth, window_length=win, polyorder=min(3, win - 2), mode='interp')
        elif m == 'bspline':
            k = min(3, max(1, len(x) - 1))
            s = max(0.0, float(bspline_smooth)) * len(x)
            smooth = UnivariateSpline(x, y_interp, k=k, s=s)(x)
        else:
            smooth = savgol_filter(y_interp, window_length=win, polyorder=min(3, win - 2), mode='interp')

        # Iterative clipping similar to apnormalize-style robust fitting.
        fit = smooth.copy()
        keep = np.ones_like(y_interp, dtype=bool)
        for _ in range(3):
            resid = y_interp - fit
            sigma = np.std(resid[keep]) if np.any(keep) else 0.0
            if sigma <= 0:
                break
            new_keep = np.abs(resid) < 3.0 * sigma
            if np.array_equal(new_keep, keep):
                break
            keep = new_keep
            if np.sum(keep) < 10:
                break
            fit = np.interp(x, x[keep], y_interp[keep])
            fit = savgol_filter(fit, window_length=win, polyorder=3, mode='interp')

        # Safety: no negative or zero values
        fit = np.maximum(fit, 1e-6)
        return fit

    def _smooth_blaze_auto_bspline(self, flux_1d: np.ndarray,
                                   smooth_factor: float = 1.0,
                                   prefilter_size: int = 15) -> np.ndarray:
        """Automated B-spline blaze extraction.
        
        Uses SciPy's automated knot selection driven by a dynamic smoothing condition
        derived from local noise variance. This eliminates the need for manual knot
        spacing and perfectly tracks complex order shapes natively.
        """
        from scipy.interpolate import splrep, splev

        x = np.arange(flux_1d.size)
        y = np.asarray(flux_1d, dtype=float)

        valid = np.isfinite(y) & (y > 0)
        n_valid = int(np.sum(valid))
        if n_valid < max(20, int(0.2 * y.size)):
            med = np.nanmedian(y[valid]) if np.any(valid) else 1.0
            return np.full_like(y, med if np.isfinite(med) and med > 0 else 1.0)

        actual_prefilter = min(prefilter_size, max(3, n_valid // 5))
        if actual_prefilter % 2 == 0:
            actual_prefilter += 1

        # Pre-clean: median filter kills local outliers
        y_clean = ndimage.median_filter(
            np.interp(x, x[valid], y[valid]), size=actual_prefilter, mode='nearest')

        # Fit region: > 1% of peak to allow knots to reach the very edges naturally
        peak = np.max(y_clean)
        sig_mask = y_clean > (peak * 0.01)
        x_fit = x[sig_mask].astype(float)
        y_fit = y_clean[sig_mask]
        
        if len(x_fit) < 10:
            return np.maximum(y_clean, 1.0)
            
        # BUG FIX: Estimate noise from the RAW unsmoothed data!
        # Using y_fit (which is median-filtered) yields a variance near 0,
        # causing 's' to be tiny and forcing B-spline to overfit all wiggles.
        y_raw = np.interp(x, x[valid], y[valid])
        y_raw_fit = y_raw[sig_mask]
        
        # Estimate true photon noise for 's' parameter (MAD of raw differences)
        noise_est = np.nanmedian(np.abs(np.diff(y_raw_fit))) / 0.6745
        var_est = max(noise_est ** 2, 1.0)
        
        # s = N * variance is the statistical sweet spot. 
        s_val = smooth_factor * len(x_fit) * var_est
        
        # Cubic B-spline fit with automatic knot selection!
        tck = splrep(x_fit, y_fit, k=3, s=s_val)
        
        # Natural extrapolation to the full array
        blaze_smooth = splev(x, tck)
        
        # --- NEW: Wing Deviation Check & Local Smoothing ---
        # Check outer 10% of the array for wild extrapolations (up/down curls)
        n_pts = len(x)
        wing_len = int(0.1 * n_pts)
        
        if wing_len > 10:
            trans_len = max(1, int(0.05 * n_pts))  # Transition zone is the inner 5% of the wing
            flat_len = wing_len - trans_len        # Absolute local smoothing zone is the outer 5%
            
            # Threshold: 10x noise or 5% of peak, whichever is larger
            dev_threshold = max(10.0 * noise_est, 0.05 * peak)
            
            macro_trend = None  # Compute lazily to save time if not needed
            
            # 1. Left wing check
            left_wing_spline = blaze_smooth[:wing_len]
            left_wing_raw = y_clean[:wing_len]
            left_dev = np.max(np.abs(left_wing_spline - left_wing_raw))
            
            # Trigger if deviation is too large OR if the spline dives below 1.0 (causing a flat cutoff)
            if left_dev > dev_threshold or np.min(left_wing_spline) < 1.0:
                if macro_trend is None:
                    macro_trend = ndimage.gaussian_filter1d(y_clean, sigma=max(20.0, wing_len / 4.0))
                
                # Cosine blend: 1.0 at x=0, smoothly transitions to 0.0 at x=wing_len
                w_l = np.ones(wing_len, dtype=float)
                x_trans = np.arange(trans_len, dtype=float)
                w_l[flat_len:] = 0.5 * (1.0 + np.cos(np.pi * x_trans / trans_len))
                blaze_smooth[:wing_len] = w_l * macro_trend[:wing_len] + (1.0 - w_l) * left_wing_spline
                
            # 2. Right wing check
            right_wing_spline = blaze_smooth[-wing_len:]
            right_wing_raw = y_clean[-wing_len:]
            right_dev = np.max(np.abs(right_wing_spline - right_wing_raw))
            
            if right_dev > dev_threshold or np.min(right_wing_spline) < 1.0:
                if macro_trend is None:
                    macro_trend = ndimage.gaussian_filter1d(y_clean, sigma=max(20.0, wing_len / 4.0))
                
                # Cosine blend: 0.0 at x=start_idx, smoothly transitions to 1.0 at the end
                w_r = np.ones(wing_len, dtype=float)
                x_trans = np.arange(trans_len, dtype=float)
                w_r[:trans_len] = 0.5 * (1.0 - np.cos(np.pi * x_trans / trans_len))
                
                start_idx = n_pts - wing_len
                blaze_smooth[start_idx:] = w_r * macro_trend[start_idx:] + (1.0 - w_r) * right_wing_spline

        # Safety: no negative or zero values
        return blaze_smooth

    def _smooth_blaze_chebyshev(self, flux_1d: np.ndarray,
                                degree: int = 5,
                                prefilter_size: int = 51) -> np.ndarray:
        """Smooth blaze using a low-order Chebyshev polynomial.
        
        This enforces a rigid global bow-shape curve, which is perfect for very bright 
        orders where B-splines might accidentally overfit high-frequency interference 
        fringes (ripples).
        """
        from numpy.polynomial.chebyshev import chebfit, chebval
        
        x = np.arange(flux_1d.size)
        y = np.asarray(flux_1d, dtype=float)

        valid = np.isfinite(y) & (y > 0)
        n_valid = int(np.sum(valid))
        if n_valid < max(20, int(0.2 * y.size)):
            med = np.nanmedian(y[valid]) if np.any(valid) else 1.0
            return np.full_like(y, med if np.isfinite(med) and med > 0 else 1.0)

        actual_prefilter = min(prefilter_size, max(3, n_valid // 5))
        if actual_prefilter % 2 == 0:
            actual_prefilter += 1

        # Pre-clean: median filter kills local outliers
        y_clean = ndimage.median_filter(
            np.interp(x, x[valid], y[valid]), size=actual_prefilter, mode='nearest')

        # Determine the valid (signal) region: > 5% of peak
        peak = np.max(y_clean)
        sig_mask = y_clean > (peak * 0.05)
        sig_idx = np.where(sig_mask)[0]
        if sig_idx.size < degree + 2:
            return y_clean.copy()
            
        x_fit = x[sig_mask].astype(float)
        y_fit = y_clean[sig_mask]
        
        # Map x to canonical [-1, 1] domain for Chebyshev
        x_min, x_max = x_fit[0], x_fit[-1]
        span = max(x_max - x_min, 1e-6)
        x_mapped = 2.0 * (x_fit - x_min) / span - 1.0
        
        # Robust iterative fitting with sigma-clipping
        keep = np.ones(len(x_fit), dtype=bool)
        coef = None
        for _ in range(5):
            if np.sum(keep) < degree + 2:
                break
            coef = chebfit(x_mapped[keep], y_fit[keep], degree)
            resid = y_fit - chebval(x_mapped, coef)
            std = np.std(resid[keep])
            if std < 1e-5:
                break
            # Symmetric clip is ideal for symmetric interference ripples
            new_keep = np.abs(resid) < 2.5 * std
            if np.array_equal(new_keep, keep):
                break
            keep = new_keep
            
        if coef is None:
            coef = chebfit(x_mapped, y_fit, degree)
            
        # Evaluate on all pixels, but clamp domain to [-1, 1] to avoid wild polynomial extrapolation
        x_all_mapped = 2.0 * (x.astype(float) - x_min) / span - 1.0
        x_all_mapped_clamped = np.clip(x_all_mapped, -1.0, 1.0)
        blaze_smooth = chebval(x_all_mapped_clamped, coef)
        
        # Safety clamp: no negative or zero values
        blaze_smooth = np.maximum(blaze_smooth, 1.0)
        
        return blaze_smooth

    def _smooth_blaze_upper_envelope_bspline(self, flux_1d: np.ndarray,
                                             knot_spacing: int = 200,
                                             edge_nknots: int = 6,
                                             prefilter_size: int = 15) -> np.ndarray:
        """Smooth blaze using an Iterative Asymmetric Upper-Envelope B-spline.
        
        This method forces the spline to rest on the 'peaks' of the interference fringes
        (the upper envelope) by assigning high weights to points above the fit and extremely 
        low weights to points in the fringe valleys.
        """
        from scipy.interpolate import splrep, splev
        
        x = np.arange(flux_1d.size)
        y = np.asarray(flux_1d, dtype=float)

        valid = np.isfinite(y) & (y > 0)
        n_valid = int(np.sum(valid))
        if n_valid < max(20, int(0.2 * y.size)):
            med = np.nanmedian(y[valid]) if np.any(valid) else 1.0
            return np.full_like(y, med if np.isfinite(med) and med > 0 else 1.0)

        actual_prefilter = min(prefilter_size, max(3, n_valid // 5))
        if actual_prefilter % 2 == 0:
            actual_prefilter += 1

        # Pre-clean: median filter kills sharp shot-noise so peaks are reliable
        y_clean = ndimage.median_filter(
            np.interp(x, x[valid], y[valid]), size=actual_prefilter, mode='nearest')

        # Determine the valid (signal) region: > 5% of peak
        peak = np.max(y_clean)
        sig_mask = y_clean > (peak * 0.05)
        sig_idx = np.where(sig_mask)[0]
        if sig_idx.size < 50:
            return y_clean.copy()
            
        x_fit = x[sig_mask].astype(float)
        y_fit = y_clean[sig_mask]
        
        # --- Build knot vector ---
        start, end = float(x_fit[0]), float(x_fit[-1])
        span = end - start
        edge_zone = max(50.0, span * 0.10)
        
        mid_knots = np.arange(start + edge_zone, end - edge_zone, knot_spacing, dtype=float)
        left_offsets = np.geomspace(10, edge_zone, edge_nknots)
        right_offsets = np.geomspace(10, edge_zone, edge_nknots)
        
        all_knots = np.sort(np.unique(np.concatenate([start + left_offsets, mid_knots, end - right_offsets])))
        all_knots = all_knots[(all_knots > x_fit[0]) & (all_knots < x_fit[-1])]
        
        # --- Iterative Asymmetric Fitting (Upper Envelope) ---
        weights = np.ones(len(x_fit), dtype=float)
        tck = None
        
        for iteration in range(5):
            try:
                if all_knots.size >= 2:
                    tck = splrep(x_fit, y_fit, w=weights, t=all_knots, k=3, task=-1)
                else:
                    tck = splrep(x_fit, y_fit, w=weights, k=3, s=len(x_fit))
            except Exception:
                tck = splrep(x_fit, y_fit, w=weights, k=3, s=len(x_fit)*2.0)
                
            current_fit = splev(x_fit, tck)
            resid = y_fit - current_fit
            
            # Asymmetric weighting: Points ABOVE the fit (peaks) get weight 1.0.
            # Points BELOW the fit (fringe valleys) get weight 0.01 to ignore them.
            weights = np.where(resid > 0, 1.0, 0.01)
            
            # Force high weights at the very edges to anchor the spline from flying away
            weights[:20] = 1.0
            weights[-20:] = 1.0
            
        blaze_smooth = splev(x, tck)
        
        return np.maximum(blaze_smooth, 1.0)

    def _smooth_blaze_upper_envelope_chebyshev(self, flux_1d: np.ndarray,
                                               degree: int = 4,
                                               prefilter_size: int = 21) -> np.ndarray:
        """Smooth blaze using an Iterative Asymmetric Upper-Envelope Chebyshev polynomial.
        
        Combines the extreme rigidity of a low-order Chebyshev polynomial (to bridge over 
        deep local pits and fringes) with the asymmetric weighting of an upper-envelope 
        fit (to perfectly rest on the wave peaks without being dragged down by the valleys).
        """
        from numpy.polynomial.chebyshev import chebfit, chebval
        
        x = np.arange(flux_1d.size)
        y = np.asarray(flux_1d, dtype=float)

        valid = np.isfinite(y) & (y > 0)
        n_valid = int(np.sum(valid))
        if n_valid < max(20, int(0.2 * y.size)):
            med = np.nanmedian(y[valid]) if np.any(valid) else 1.0
            return np.full_like(y, med if np.isfinite(med) and med > 0 else 1.0)

        actual_prefilter = min(prefilter_size, max(3, n_valid // 5))
        if actual_prefilter % 2 == 0:
            actual_prefilter += 1

        # Pre-clean: median filter kills sharp shot-noise so peaks are reliable
        y_clean = ndimage.median_filter(
            np.interp(x, x[valid], y[valid]), size=actual_prefilter, mode='nearest')

        # Determine the valid (signal) region: > 5% of peak
        peak = np.max(y_clean)
        sig_mask = y_clean > (peak * 0.05)
        sig_idx = np.where(sig_mask)[0]
        if sig_idx.size < degree + 2:
            return y_clean.copy()
            
        x_fit = x[sig_mask].astype(float)
        y_fit = y_clean[sig_mask]
        
        x_min, x_max = x_fit[0], x_fit[-1]
        span = max(x_max - x_min, 1e-6)
        x_mapped = 2.0 * (x_fit - x_min) / span - 1.0
        
        # --- Iterative Asymmetric Fitting (Upper Envelope) ---
        weights = np.ones(len(x_fit), dtype=float)
        coef = None
        
        edge_len = min(20, len(x_fit) // 10)
        
        for iteration in range(5):
            coef = chebfit(x_mapped, y_fit, degree, w=weights)
            current_fit = chebval(x_mapped, coef)
            resid = y_fit - current_fit
            
            # Asymmetric weighting: Points ABOVE the fit (peaks) get weight 1.0.
            # Points BELOW the fit (fringe valleys) get weight 0.01.
            weights = np.where(resid > 0, 1.0, 0.01)
            
            # Force high weights at the very edges to anchor the polynomial
            if edge_len > 0:
                weights[:edge_len] = 1.0
                weights[-edge_len:] = 1.0
                
        # Evaluate on all pixels, clamping the domain to [-1, 1] to avoid wild extrapolation
        x_all_mapped = 2.0 * (x.astype(float) - x_min) / span - 1.0
        x_all_mapped_clamped = np.clip(x_all_mapped, -1.0, 1.0)
        blaze_smooth = chebval(x_all_mapped_clamped, coef)
        
        return np.maximum(blaze_smooth, 1.0)

    def build_order_response_map(self, apertures: ApertureSet,
                                 source_image: Optional[np.ndarray] = None,
                                 blaze_smooth_factor: float = 1.0,
                                 width_smooth_window: int = 41,
                                 profile_bin_step: float = 0.01,
                                 n_profile_segments: int = 100,
                                 profile_smooth_sigma: float = 6.0,
                                 pixel_flat_min: float = 0.5,
                                 pixel_flat_max: float = 1.5,
                                 fringe_orders: int = 20
                                 ) -> Tuple[np.ndarray, Dict[int, np.ndarray], Dict[int, np.ndarray], np.ndarray, np.ndarray, Optional[np.ndarray]]:
        """
        Build per-order response model and convert it to a 2D sensitivity map.

        This follows the echelle flat-field philosophy used in IRAF/PypeIt/REDUCE:
        estimate order-by-order blaze/response first, then construct a 2D correction map.

        Improvements over a single global profile:
        1. **Segmented cross-dispersion profiles**: The profile is built in
           ``n_profile_segments`` bins along X, then linearly interpolated at
           each column, so the model tracks PSF shape changes (aberrations,
           defocus, vignetting).
        2. **Sparse-knots B-spline blaze fitting**: Interior knots are placed
           at wide spacing (``blaze_knot_spacing``) to crush fringe ripples,
           while dense knots near edges track the steep blaze roll-off
           — so the 2D model is fringe-free and the resulting pixel_flat
           absorbs fringe structure for correction.
        """
        if self.flat_data is None and source_image is None:
            raise RuntimeError("No flat data")

        image = source_image if source_image is not None else self.flat_data
        image = np.asarray(image, dtype=np.float64)
        height, width = image.shape
        x_coords = np.arange(width)

        # Clamp segments to a sane range
        n_profile_segments = max(1, min(n_profile_segments, 200))

        response_map = np.zeros_like(image, dtype=np.float32)
        coverage = np.zeros_like(image, dtype=bool)

        blaze_profiles: Dict[int, np.ndarray] = {}
        cross_profiles: Dict[int, np.ndarray] = {}
        order_diag: Dict[int, Dict[str, np.ndarray]] = {}

        # --- First Pass: Extract raw blaze for all orders and compute intensity percentiles ---
        extracted_raw_blaze = {}
        extracted_width_raw = {}
        extracted_bounds = {}

        for aperture_id, aperture in apertures.apertures.items():
            lower = aperture.get_lower(x_coords)
            upper = aperture.get_upper(x_coords)

            # 1D aperture sum extraction inside traced boundaries.
            raw_blaze = np.full(width, np.nan, dtype=float)
            width_raw = np.full(width, np.nan, dtype=float)
            per_col_bounds = []

            for x in x_coords:
                y_lo = min(lower[x], upper[x])
                y_hi = max(lower[x], upper[x])
                y1 = max(0, int(np.floor(y_lo)))
                y2 = min(height, int(np.ceil(y_hi + 1)))
                if y2 - y1 < 3:
                    per_col_bounds.append(None)
                    continue

                stripe = image[y1:y2, x]
                yy = np.arange(y1, y2)
                if self.flat_mask is not None:
                    local_bad = self.flat_mask[y1:y2, x] > 0
                    if np.all(local_bad):
                        per_col_bounds.append(None)
                        continue
                    stripe = stripe[~local_bad]
                    yy = yy[~local_bad]

                if stripe.size < 3:
                    per_col_bounds.append(None)
                    continue

                raw_blaze[x] = np.sum(stripe)
                width_raw[x] = max(1.0, y_hi - y_lo)
                per_col_bounds.append((yy, stripe))

            extracted_raw_blaze[aperture_id] = raw_blaze
            extracted_width_raw[aperture_id] = width_raw
            extracted_bounds[aperture_id] = per_col_bounds

        # Determine cutoff for fringe correction based on user parameter
        oids = sorted(extracted_raw_blaze.keys())
        if oids and fringe_orders > 0:
            n_fringed = min(len(oids), fringe_orders)
            fringed_oids = set(oids[:n_fringed])
            is_fringed = {oid: (oid in fringed_oids) for oid in oids}
            logger.info(f"Applying Upper Envelope Chebyshev Fit to the first {n_fringed} red-end orders.")
        else:
            is_fringed = {}
            n_fringed = 0

        # --- Second Pass: Smooth blaze and reconstruct 2D model with dynamic parameters ---
        for aperture_id, aperture in apertures.apertures.items():
            center = aperture.get_position(x_coords)
            lower = aperture.get_lower(x_coords)
            upper = aperture.get_upper(x_coords)
            
            raw_blaze = extracted_raw_blaze[aperture_id]
            width_raw = extracted_width_raw[aperture_id]
            per_col_bounds = extracted_bounds[aperture_id]

            # Dynamic rule based on auto-detected fringes
            if is_fringed.get(aperture_id, False):
                # Use Iterative Asymmetric Upper-Envelope Chebyshev to track peaks with rigid bow shape.
                blaze_1d = self._smooth_blaze_upper_envelope_chebyshev(
                    raw_blaze,
                    degree=6,
                    prefilter_size=21
                )
                fit_method_label = 'Upper Env Chebyshev'
            else:
                # Medium & Dark blue orders: use automated B-spline.
                # Create a two-step "transition zone" just past the fringe zone:
                # +1 to +5 orders: multiply smooth factor by 50.0
                # +6 to +15 orders: multiply smooth factor by 15.0
                # This irons out small-scale wiggles gradually before reverting to normal flexibility.
                current_factor = blaze_smooth_factor
                order_idx = oids.index(aperture_id) if aperture_id in oids else -1
                
                if n_fringed <= order_idx < n_fringed + 5:
                    current_factor *= 50.0
                    fit_method_label = 'Rigid Auto B-spline'
                elif n_fringed + 5 <= order_idx < n_fringed + 15:
                    current_factor *= 15.0
                    fit_method_label = 'Rigid Auto B-spline'
                else:
                    fit_method_label = 'Auto B-spline'

                blaze_1d = self._smooth_blaze_auto_bspline(
                    raw_blaze,
                    smooth_factor=current_factor,
                    prefilter_size=21
                )

            width_smooth = self._smooth_order_response(
                width_raw,
                method='savgol',
                window=width_smooth_window,
                bspline_smooth=0.0,
            )
            width_smooth = np.clip(width_smooth, 1.0, None)

            norm = np.nanmedian(blaze_1d[np.isfinite(blaze_1d) & (blaze_1d > 0)])
            if not np.isfinite(norm) or norm <= 0:
                norm = 1.0
            blaze = blaze_1d / norm
            blaze = np.where(np.isfinite(blaze) & (blaze > 1e-6), blaze, 1.0)
            blaze_profiles[aperture_id] = blaze.astype(np.float32)

            valid_cols = np.where(np.isfinite(raw_blaze))[0]
            if valid_cols.size == 0:
                logger.warning(f"Order {aperture_id}: no valid columns — skipped")
                continue

            sample_cols = set()
            if valid_cols.size > 0:
                sample_cols = set(np.unique(np.linspace(valid_cols.min(), valid_cols.max(), num=min(n_profile_segments, valid_cols.size), dtype=int)))

            # ==============================================================
            # Build segmented cross-dispersion profiles along X
            # ==============================================================
            step = max(0.002, float(profile_bin_step))
            n_seg = n_profile_segments

            # Determine column ranges for each segment
            x_lo_seg = valid_cols.min()
            x_hi_seg = valid_cols.max()
            seg_edges = np.linspace(x_lo_seg, x_hi_seg + 1, n_seg + 1).astype(int)
            seg_centers = 0.5 * (seg_edges[:-1] + seg_edges[1:])  # X center of each segment

            # Build one profile per segment
            seg_profiles_x = []   # list of (prof_x_i,) arrays — normalized y coordinates
            seg_profiles_y = []   # list of (prof_y_i,) arrays — profile values
            common_yn_grid = None  # will be unified after first pass

            for si in range(n_seg):
                xs0, xs1 = int(seg_edges[si]), int(seg_edges[si + 1])
                yn_seg = []
                fn_seg = []
                for x in range(xs0, min(xs1, width)):
                    if per_col_bounds[x] is None:
                        continue
                    yy, stripe = per_col_bounds[x]
                    bval = blaze_1d[x]
                    wval = width_smooth[x]
                    if (not np.isfinite(bval)) or bval <= 0 or (not np.isfinite(wval)) or wval <= 0:
                        continue
                    yn = (yy.astype(float) - center[x]) / wval
                    fn = stripe / bval
                    good = np.isfinite(yn) & np.isfinite(fn)
                    if np.any(good):
                        yn_seg.append(yn[good])
                        fn_seg.append(fn[good])

                if not yn_seg:
                    seg_profiles_x.append(None)
                    seg_profiles_y.append(None)
                    continue

                yn_cat = np.concatenate(yn_seg)
                fn_cat = np.concatenate(fn_seg)

                yn_min = max(-1.2, float(np.nanpercentile(yn_cat, 1)))
                yn_max_v = min(1.2, float(np.nanpercentile(yn_cat, 99)))
                if yn_max_v - yn_min < 0.2:
                    yn_min, yn_max_v = -0.8, 0.8
                bins = np.arange(yn_min, yn_max_v + step, step)
                if bins.size < 8:
                    bins = np.linspace(-0.8, 0.8, 81)

                idx = np.digitize(yn_cat, bins) - 1
                px, py = [], []
                for i in range(len(bins) - 1):
                    pick = idx == i
                    if np.sum(pick) < 6:
                        continue
                    px.append(0.5 * (bins[i] + bins[i + 1]))
                    py.append(np.nanmedian(fn_cat[pick]))

                if len(px) < 8:
                    seg_profiles_x.append(None)
                    seg_profiles_y.append(None)
                    continue

                px = np.asarray(px, dtype=float)
                py = np.asarray(py, dtype=float)
                py = np.where(np.isfinite(py), py, 1.0)
                py = np.clip(py, 1e-6, None)

                # Smooth the binned profile
                py = self._smooth_order_response(
                    py,
                    method='savgol',
                    window=max(11, int(11 * max(1.0, step / 0.01))),
                    bspline_smooth=0.0,
                )

                seg_profiles_x.append(px)
                seg_profiles_y.append(py)

            # Count valid segments
            valid_seg_idx = [i for i in range(n_seg) if seg_profiles_x[i] is not None]

            if not valid_seg_idx:
                logger.warning(f"Order {aperture_id}: no valid profile segments — skipped")
                continue

            # Build a common y_norm grid that covers all segments
            all_yn_min = min(seg_profiles_x[i].min() for i in valid_seg_idx)
            all_yn_max = max(seg_profiles_x[i].max() for i in valid_seg_idx)
            common_yn = np.arange(all_yn_min, all_yn_max + step, step)

            # Resample each segment profile onto the common grid
            seg_on_grid = np.full((n_seg, common_yn.size), np.nan, dtype=float)
            for i in valid_seg_idx:
                seg_on_grid[i] = np.interp(common_yn, seg_profiles_x[i], seg_profiles_y[i],
                                           left=1e-6, right=1e-6)

            # Fill missing segments by nearest-neighbour copy
            for i in range(n_seg):
                if i in valid_seg_idx:
                    continue
                # Find nearest valid segment
                dists = [abs(i - vi) for vi in valid_seg_idx]
                nearest = valid_seg_idx[int(np.argmin(dists))]
                seg_on_grid[i] = seg_on_grid[nearest]

            # ----------------------------------------------------------
            # Smooth profiles along dispersion direction (axis 0 = segment index)
            # so that CSP varies continuously from segment to segment.
            # ----------------------------------------------------------
            if n_seg >= 3 and profile_smooth_sigma > 0:
                for j in range(seg_on_grid.shape[1]):
                    col = seg_on_grid[:, j]
                    ok = np.isfinite(col)
                    if np.sum(ok) < 3:
                        continue
                    # Interpolate NaN gaps before smoothing
                    idx_ok = np.where(ok)[0]
                    col_filled = np.interp(np.arange(n_seg), idx_ok, col[idx_ok])
                    seg_on_grid[:, j] = gaussian_filter1d(
                        col_filled, sigma=profile_smooth_sigma, mode='nearest')

            # Store the composite profile for diagnostics (average of all segments)
            with np.errstate(all='ignore'):
                mean_prof_y = np.nanmean(seg_on_grid, axis=0)
            mean_prof_y = np.where(np.isfinite(mean_prof_y), mean_prof_y, 1e-6)
            csum = np.nansum(mean_prof_y)
            if np.isfinite(csum) and csum > 0:
                cross_profiles[aperture_id] = (mean_prof_y / csum).astype(np.float32)

            sampled_realspace = []

            # ==============================================================
            # Reconstruct 2D model: Model(X,Y) = P(y_norm; X) × Blaze(X)
            # with P interpolated between segment profiles
            # ==============================================================
            for x in x_coords:
                y_lo = min(lower[x], upper[x])
                y_hi = max(lower[x], upper[x])
                y1 = max(0, int(np.floor(y_lo)))
                y2 = min(height, int(np.ceil(y_hi + 1)))
                if y2 - y1 < 3:
                    continue

                bval = blaze_1d[x]
                wval = width_smooth[x]
                if (not np.isfinite(bval)) or bval <= 0 or (not np.isfinite(wval)) or wval <= 0:
                    continue

                yy = np.arange(y1, y2, dtype=float)
                yreq = (yy - center[x]) / wval

                # Interpolate profile between the two nearest segments
                if n_seg == 1:
                    pval = np.interp(yreq, common_yn, seg_on_grid[0], left=1e-6, right=1e-6)
                else:
                    # Find the fractional segment index for this column
                    frac = np.clip((float(x) - seg_centers[0]) / max(1.0, seg_centers[-1] - seg_centers[0]),
                                   0.0, 1.0)
                    seg_f = frac * (n_seg - 1)
                    si_lo = max(0, min(int(np.floor(seg_f)), n_seg - 2))
                    si_hi = si_lo + 1
                    w_hi = seg_f - si_lo
                    w_lo = 1.0 - w_hi

                    p_lo = np.interp(yreq, common_yn, seg_on_grid[si_lo], left=1e-6, right=1e-6)
                    p_hi = np.interp(yreq, common_yn, seg_on_grid[si_hi], left=1e-6, right=1e-6)
                    pval = w_lo * p_lo + w_hi * p_hi

                model_col = pval * bval

                if x in sample_cols:
                    sampled_realspace.append({
                        'x': int(x),
                        'y': yy.astype(np.float32),
                        'ynorm': yreq.astype(np.float32),
                        'profile_val': pval.astype(np.float32),
                        'model_flux': model_col.astype(np.float32),
                    })

                response_map[y1:y2, x] = model_col.astype(np.float32)
                coverage[y1:y2, x] = True

            order_diag[aperture_id] = {
                'raw_blaze': raw_blaze.astype(np.float32),
                'blaze_fit': blaze_1d.astype(np.float32),
                'width_smooth': width_smooth.astype(np.float32),
                'profile_x': common_yn.astype(np.float32),
                'profile_y': mean_prof_y.astype(np.float32),
                'seg_on_grid': seg_on_grid.astype(np.float32),
                'seg_centers': seg_centers.astype(np.float32),
                'sampled_realspace': sampled_realspace,
                'n_profile_segments': n_seg,
                'n_valid_segments': len(valid_seg_idx),
                'blaze_method': fit_method_label,
            }

        # Smoothed model keeps blaze + cross-order envelope, without pixel noise.
        smoothed_model = response_map.astype(np.float32)

        # Pixel-to-pixel flat: high-frequency component around unity.

        pixel_flat = np.ones_like(image, dtype=np.float32)
        pixel_flat[coverage] = image[coverage] / (smoothed_model[coverage] + 1e-10)
        bad_cov = coverage & (~np.isfinite(pixel_flat) | (pixel_flat < pixel_flat_min) | (pixel_flat > pixel_flat_max))
        if np.any(bad_cov):
            pixel_flat[bad_cov] = 1.0
        if np.any(coverage):
            pf_med = np.nanmedian(pixel_flat[coverage][np.isfinite(pixel_flat[coverage])])
            if np.isfinite(pf_med) and pf_med > 0:
                pixel_flat[coverage] = pixel_flat[coverage] / pf_med
        pixel_flat = np.where(np.isfinite(pixel_flat), pixel_flat, 1.0)
        pixel_flat[~coverage] = 1.0
        
        # 将级次上下边界 1 pixel 的区域强制设为 1.0，防止边缘极低信噪比带来的除法突变
        from scipy.ndimage import binary_erosion
        eroded_coverage = binary_erosion(coverage, structure=np.ones((3, 1), dtype=bool))
        pixel_flat[coverage & (~eroded_coverage)] = 1.0
        
        # 进一步消除级次内部的异常突变：将信噪比低不可靠的极端像元限制得更严格些
        edge_mask = coverage & ((pixel_flat < 0.8) | (pixel_flat > 1.2))
        if np.any(edge_mask):
            pixel_flat[edge_mask] = 1.0
        pixel_flat = np.clip(pixel_flat, 0.1, 10.0).astype(np.float32)

        # Illumination flat is intentionally not used/output in pixel-to-pixel-only workflow.
        illumination_flat = None

        # Final 2D correction map equals pixel flat in this workflow.
        flat_corr_2d = np.clip(pixel_flat, 0.1, 10.0).astype(np.float32)

        self.response_map = response_map
        self.smoothed_model = smoothed_model
        self.pixel_flat = pixel_flat
        self.illumination_flat = illumination_flat
        self.flat_corr_2d = flat_corr_2d
        self.flat_sens = None
        self.order_diagnostics = order_diag
        logger.info("Built Step4 flat model via 1D blaze + normalized-profile reconstruction")
        return None, blaze_profiles, cross_profiles, smoothed_model, pixel_flat, illumination_flat

    def save_step4_diagnostics(self, output_dir: Path,
                               clean_flat: Optional[np.ndarray] = None,
                               clean_flat_corrected: Optional[np.ndarray] = None,
                               fig_format: str = 'png',
                               save_plots: bool = True):
        """Save step4 diagnostic artifacts: blaze/profile/inverse-map/model/pixel-flat and corrected flat."""
        out_dir = Path(output_dir)
        out_dir.mkdir(parents=True, exist_ok=True)

        # fig_format and save_plots are now passed via method parameters
        # (kept as instance defaults for backward compat with save_step4_diagnostics signature)

        # ---- FITS outputs ----
        if self.response_map is not None:
            write_fits_image(str(out_dir / 'flat_smoothed_model.fits'), self.response_map.astype(np.float32), dtype='float32')
        if self.pixel_flat is not None:
            write_fits_image(str(out_dir / 'flat_pixel_2d.fits'), self.pixel_flat.astype(np.float32), dtype='float32')
        if clean_flat_corrected is not None:
            write_fits_image(str(out_dir / 'MasterFlat.fits'), clean_flat_corrected.astype(np.float32), dtype='float32')

        # ---- 2D overview PNG plots ----
        if save_plots:
            if self.response_map is not None:
                plot_2d_image_to_file(self.response_map, str(out_dir / f'flat_smoothed_model.{fig_format}'), 'Reconstructed 2D Flat Model')
            if self.pixel_flat is not None:
                plot_2d_image_to_file(self.pixel_flat, str(out_dir / f'flat_pixel_2d.{fig_format}'),
                                      'Pixel-to-Pixel Flat', vmin=0.85, vmax=1.15)
            if clean_flat_corrected is not None:
                plot_2d_image_to_file(clean_flat_corrected, str(out_dir / f'MasterFlat.{fig_format}'), 'Master Flat After 2D Flat Correction')

        if not self.order_diagnostics:
            return

        import matplotlib.pyplot as plt
        from matplotlib.cm import ScalarMappable
        from matplotlib.colors import Normalize

        blaze_dir = out_dir / 'blaze_per_order'
        profile_dir = out_dir / 'profile_per_order'
        realspace_dir = out_dir / 'realspace_profile_per_order'
        blaze_dir.mkdir(parents=True, exist_ok=True)
        profile_dir.mkdir(parents=True, exist_ok=True)
        realspace_dir.mkdir(parents=True, exist_ok=True)

        # Collect data for combined overview plots.
        all_order_ids = sorted(self.order_diagnostics.keys())

        # ---- Per-order diagnostic loop ----
        for aperture_id, diag in self.order_diagnostics.items():
            raw_blaze = diag.get('raw_blaze')
            blaze_fit = diag.get('blaze_fit')
            blaze_method = diag.get('blaze_method', 'B-spline Fit')
            profile_x = diag.get('profile_x')
            profile_y = diag.get('profile_y')
            seg_on_grid = diag.get('seg_on_grid')
            seg_centers = diag.get('seg_centers')
            sampled_realspace = diag.get('sampled_realspace', [])
            n_seg = diag.get('n_profile_segments', 1)
            n_valid_seg = diag.get('n_valid_segments', 0)

            # Save numerical data
            np.savez(
                blaze_dir / f'order_{aperture_id:03d}_blaze.npz',
                raw_blaze=raw_blaze,
                blaze_fit=blaze_fit,
                width_smooth=diag.get('width_smooth'),
            )
            np.savez(
                profile_dir / f'order_{aperture_id:03d}_profile.npz',
                y_norm=profile_x,
                profile_mean=profile_y,
                seg_on_grid=seg_on_grid,
                seg_centers=seg_centers,
            )

            if not save_plots:
                continue

            # --- Per-order blaze plot ---
            if raw_blaze is not None and blaze_fit is not None:
                plt.figure(figsize=(10, 4))
                x = np.arange(raw_blaze.size)
                plt.plot(x, raw_blaze, color='0.6', linewidth=0.8, label='Raw Blaze 1D')
                plt.plot(x, blaze_fit, color='tab:red', linewidth=1.4, label=blaze_method)
                plt.xlabel('X (column)')
                plt.ylabel('Flux')
                plt.title(f'Order {aperture_id} — Blaze (segments={n_seg}, valid={n_valid_seg})')
                plt.legend()
                plt.grid(alpha=0.2)
                plt.tight_layout()
                plt.savefig(blaze_dir / f'order_{aperture_id:03d}_blaze.{fig_format}', dpi=150)
                plt.close()

            # --- Per-order profile plot (2D heatmap + mean) ---
            if profile_x is not None and profile_y is not None:
                fig = plt.figure(figsize=(10, 4.5))
                
                if seg_on_grid is not None and seg_centers is not None and n_seg > 1:
                    ax1 = plt.subplot(121)
                    # seg_on_grid is shape (n_seg, n_y). We transpose to (n_y, n_seg)
                    img_data = seg_on_grid.T
                    # Origin='lower' so y goes up.
                    extent = [seg_centers[0], seg_centers[-1], profile_x[0], profile_x[-1]]
                    im = ax1.imshow(img_data, aspect='auto', origin='lower', extent=extent, cmap='viridis')
                    ax1.set_xlabel('Dispersion X (pixel)')
                    ax1.set_ylabel('Spatial y_norm')
                    ax1.set_title('Profile Variation along Dispersion')
                    fig.colorbar(im, ax=ax1, label='Normalized Flux')
                    
                    ax2 = plt.subplot(122)
                    
                    # Draw segment profiles as colored lines
                    cmap = plt.cm.plasma
                    norm = Normalize(vmin=seg_centers[0], vmax=seg_centers[-1])
                    for i in range(n_seg):
                        color = cmap(norm(seg_centers[i]))
                        ax2.plot(profile_x, seg_on_grid[i], color=color, linewidth=0.8, alpha=0.6)
                        
                    # Add colorbar for lines
                    sm = ScalarMappable(cmap=cmap, norm=norm)
                    sm.set_array([])
                    fig.colorbar(sm, ax=ax2, label='Dispersion X (pixel)')
                else:
                    ax2 = plt.subplot(111)
                    
                ax2.plot(profile_x, profile_y, color='black', linewidth=1.8, label='Mean profile', zorder=10)
                ax2.set_xlabel('Spatial y_norm')
                ax2.set_ylabel('Normalized Flux')
                ax2.set_title(f'Order {aperture_id} — Cross-Dispersion Profile')
                ax2.grid(alpha=0.2)
                ax2.legend()
                
                plt.suptitle(f'Order {aperture_id} — Cross-Dispersion Profile ({n_seg} segments)')
                plt.tight_layout()
                plt.savefig(profile_dir / f'order_{aperture_id:03d}_profile.{fig_format}', dpi=150)
                plt.close()

            # --- Per-order real-space profile plot ---
            if sampled_realspace:
                fig = plt.figure(figsize=(10, 4.5))
                sampled_realspace.sort(key=lambda d: d['x'])
                x_vals = [d['x'] for d in sampled_realspace]
                
                cmap = plt.cm.plasma
                norm = Normalize(vmin=min(x_vals), vmax=max(x_vals))
                
                if len(sampled_realspace) > 1:
                    ax1 = plt.subplot(121)
                    xs_arr, ys_arr, fs_arr = [], [], []
                    for prof in sampled_realspace:
                        xs_arr.extend([prof['x']] * len(prof['y']))
                        ys_arr.extend(prof['y'])
                        fs_arr.extend(prof['model_flux'])
                        
                    sc = ax1.scatter(xs_arr, ys_arr, c=fs_arr, cmap='viridis', s=5, marker='s', edgecolors='none')
                    ax1.set_xlabel('Dispersion X (pixel)')
                    ax1.set_ylabel('Physical Y (pixel)')
                    ax1.set_title('Real-Space Model Heatmap')
                    fig.colorbar(sc, ax=ax1, label='Model Flux')
                    
                    ax2 = plt.subplot(122)
                else:
                    ax2 = plt.subplot(111)
                    
                for prof in sampled_realspace:
                    color = cmap(norm(prof['x']))
                    ax2.plot(prof['y'], prof['model_flux'], color=color, linewidth=0.8, alpha=0.6)
                    
                if len(sampled_realspace) > 1:
                    sm = ScalarMappable(cmap=cmap, norm=norm)
                    sm.set_array([])
                    fig.colorbar(sm, ax=ax2, label='Dispersion X (pixel)')
                    
                ax2.set_xlabel('Physical Y (pixel)')
                ax2.set_ylabel('Model Flux')
                ax2.set_title(f'Order {aperture_id} — Real-Space Profiles')
                ax2.grid(alpha=0.2)
                
                plt.suptitle(f'Order {aperture_id} — Real-Space Profiles')
                plt.tight_layout()
                plt.savefig(realspace_dir / f'order_{aperture_id:03d}_realspace_profiles.{fig_format}', dpi=150)
                plt.close()

        # ================================================================
        # Combined overview: all orders' blaze in one figure
        # ================================================================
        if save_plots and all_order_ids:
            n_orders = len(all_order_ids)
            cmap = plt.cm.turbo
            # Reverse: low order IDs (bottom, long wavelength) → red,
            #          high order IDs (top, short wavelength) → blue.
            norm_c = Normalize(vmin=0, vmax=max(1, n_orders - 1))

            def _order_color(i):
                return cmap(norm_c(n_orders - 1 - i))

            # --- All-orders blaze overview ---
            fig, ax = plt.subplots(figsize=(14, 6), constrained_layout=True)
            for i, oid in enumerate(all_order_ids):
                diag = self.order_diagnostics[oid]
                bfit = diag.get('blaze_fit')
                if bfit is None:
                    continue
                x = np.arange(bfit.size)
                clr = _order_color(i)
                ax.plot(x, bfit, color=clr, linewidth=0.6, alpha=0.8)
            ax.set_xlabel('X (column)')
            ax.set_ylabel('Heavy-Smooth Blaze')
            ax.set_title(f'All Orders Blaze Overview ({n_orders} orders)')
            ax.grid(alpha=0.15)
            cmap_r = plt.cm.turbo_r
            sm = ScalarMappable(cmap=cmap_r, norm=norm_c)
            sm.set_array([])
            cbar = fig.colorbar(sm, ax=ax, shrink=0.8, pad=0.02)
            order_ids_arr = np.array(all_order_ids)
            cbar.set_ticks(np.linspace(0, n_orders - 1, min(10, n_orders)))
            cbar.set_ticklabels([str(order_ids_arr[int(t)]) for t in np.linspace(0, n_orders - 1, min(10, n_orders))])
            cbar.set_label('Order ID')
            fig.savefig(out_dir / f'all_orders_blaze_overview.{fig_format}', dpi=150)
            plt.close(fig)

            # --- All-orders profile overview ---
            fig, ax = plt.subplots(figsize=(10, 7))
            for i, oid in enumerate(all_order_ids):
                diag = self.order_diagnostics[oid]
                px = diag.get('profile_x')
                py = diag.get('profile_y')
                if px is None or py is None:
                    continue
                clr = _order_color(i)
                ax.plot(px, py, color=clr, linewidth=0.7, alpha=0.8)
            ax.set_xlabel('y_norm')
            ax.set_ylabel('Normalized Flux')
            ax.set_title(f'All Orders Cross-Dispersion Profile Overview ({n_orders} orders)')
            ax.grid(alpha=0.15)
            sm = ScalarMappable(cmap=cmap_r, norm=norm_c)
            sm.set_array([])
            cbar = fig.colorbar(sm, ax=ax, shrink=0.7, pad=0.02)
            cbar.set_ticks(np.linspace(0, n_orders - 1, min(10, n_orders)))
            cbar.set_ticklabels([str(order_ids_arr[int(t)]) for t in np.linspace(0, n_orders - 1, min(10, n_orders))])
            cbar.set_label('Order ID')
            fig.tight_layout()
            fig.savefig(out_dir / f'all_orders_profile_overview.{fig_format}', dpi=150)
            plt.close(fig)

            logger.info(f"Step 4 diagnostics saved to {out_dir} "
                        f"({n_orders} orders, format={fig_format})")



def process_flat_correction_stage(science_image: np.ndarray,
                                   flat_field: FlatField,
                                   output_dir_base: str,
                                   apertures: Optional[ApertureSet] = None,
                                   science_name: str = 'science',
                                   # Blaze / profile parameters for on-demand model build
                                   blaze_smooth_factor: float = 1.0,
                                   width_smooth_window: int = 41,
                                   profile_bin_step: float = 0.01,
                                   n_profile_segments: int = 100,
                                   profile_smooth_sigma: float = 6.0,
                                   pixel_flat_min: float = 0.5,
                                   pixel_flat_max: float = 1.5,
                                   fringe_orders: int = 20,
                                   save_plots: bool = True,
                                   fig_format: str = 'png') -> np.ndarray:
    """
    Step 4: apply the 2D pixel-flat correction map to a science image.
    """

    processor = None
    if apertures is None:
        raise RuntimeError("Apertures are required to build the flat model, but none were provided.")
    logger.info("  Building 2D flat model (blaze × profile → pixel flat)...")
    processor = FlatCorrectionModelBuilder()
    processor.flat_data = flat_field.flat_data
    processor.flat_mask = flat_field.flat_mask

    clean_flat = flat_field.flat_data
    if flat_field.scattered_light is not None:
        clean_flat = np.clip(
            flat_field.flat_data.astype(np.float32) - flat_field.scattered_light.astype(np.float32),
            1e-6,
            None,
        )

    flat_sens, blaze_profiles, cross_profiles, smoothed_model, pixel_flat, illum = (
        processor.build_order_response_map(
            apertures, source_image=clean_flat,
            blaze_smooth_factor=blaze_smooth_factor,
            width_smooth_window=width_smooth_window,
            profile_bin_step=profile_bin_step,
            n_profile_segments=n_profile_segments,
            profile_smooth_sigma=profile_smooth_sigma,
            pixel_flat_min=pixel_flat_min,
            pixel_flat_max=pixel_flat_max,
            fringe_orders=fringe_orders)
    )
    flat_field.flat_sens = flat_sens
    flat_field.blaze_profiles = blaze_profiles
    flat_field.cross_profiles = cross_profiles
    flat_field.smoothed_model = smoothed_model
    flat_field.pixel_flat = pixel_flat
    flat_field.illumination_flat = illum

    flat_corr = flat_field.pixel_flat
    if flat_corr is None:
        flat_corr = flat_field.flat_corr_2d
    if flat_corr is None:
        raise RuntimeError("No flat correction map available. Run Step 2 first.")

    if flat_corr.shape != science_image.shape:
        raise RuntimeError(
            f"Science image shape {science_image.shape} does not match "
            f"flat correction shape {flat_corr.shape}"
        )

    safe_flat = flat_corr.astype(np.float32)
    bad = (~np.isfinite(safe_flat)) | (safe_flat <= 0.05)
    if np.any(bad):
        safe_flat = safe_flat.copy()
        safe_flat[bad] = 1.0

    corrected = science_image.astype(np.float32) / safe_flat

    out_dir = Path(output_dir_base) / 'step4_flat_corrected'
    out_dir.mkdir(parents=True, exist_ok=True)

    # Save all flat-model artefacts (passing processor if it was freshly built here to save its diagnostics).
    save_flat_correction_products(flat_field, apertures, out_dir, processor=processor,
                                  fig_format=fig_format, save_plots=save_plots)

    write_fits_image(
        str(out_dir / f'{science_name}.fits'),
        corrected,
        dtype='float32',
    )

    if save_plots:
        plot_2d_image_to_file(
            corrected,
            str(out_dir / f'{science_name}_flat2d_corrected.{fig_format}'),
            'Science After 2D Flat Correction',
        )

    logger.info(f"✓ 2D flat correction complete: {science_name}")
    return corrected
