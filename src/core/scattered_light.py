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

def create_widened_mask(apertures: ApertureSet, flat_data: np.ndarray, margin: int,
                        n_mask_below: int = 0, n_mask_above: int = 0) -> Tuple[np.ndarray, np.ndarray, np.ndarray, list]:
    """
    Create a binary mask widened by margin, strictly preserving valleys between orders.
    Extrapolates virtual orders along columns if n_mask_below/above > 0.
    Returns: (binary_mask, lo_traces, hi_traces, order_ids)
    """
    from scipy.interpolate import UnivariateSpline

    height, width = flat_data.shape
    x_all = np.arange(width, dtype=float)
    
    all_ap_ids = sorted(list(apertures.apertures.keys()))
    n_ap = len(all_ap_ids)
    
    # 1. Evaluate real orders
    real_lo = np.zeros((n_ap, width))
    real_hi = np.zeros((n_ap, width))
    real_centers = np.zeros((n_ap, width))
    
    for i, ap_id in enumerate(all_ap_ids):
        ap = apertures.get_aperture(ap_id)
        real_centers[i] = ap.get_position(x_all)
        real_lo[i] = ap.get_lower(x_all) - margin
        real_hi[i] = ap.get_upper(x_all) + margin
        
    # 2. Extrapolate virtual orders along columns
    def extrapolate_2d(data_2d, n_extrap, direction='below'):
        if n_ap < 2 or n_extrap <= 0:
            return np.zeros((0, width))
            
        res = np.zeros((n_extrap, width))
        
        # 只使用靠近外推方向的一半级次数据，避免另一端的趋势带来远端牵制偏差
        n_half = max(4, n_ap // 2) if n_ap >= 4 else n_ap
        k_deg = min(3, n_half - 1)
        
        if direction == 'below':
            x_eval = np.arange(-n_extrap, 0, dtype=float)
            fit_indices = np.arange(n_half, dtype=float)
            data_subset = data_2d[:n_half, :]
        else:
            x_eval = np.arange(n_ap, n_ap + n_extrap, dtype=float)
            fit_indices = np.arange(n_ap - n_half, n_ap, dtype=float)
            data_subset = data_2d[-n_half:, :]
            
        for xi in range(width):
            col_data = data_subset[:, xi]
            # 沿列方向利用半数探测级次进行 Spline 平滑
            spl = UnivariateSpline(fit_indices, col_data, k=k_deg, s=n_half, ext=0)
            res[:, xi] = spl(x_eval)
            
        return res

    lo_below = extrapolate_2d(real_lo, n_mask_below, 'below')
    hi_below = extrapolate_2d(real_hi, n_mask_below, 'below')
    c_below  = extrapolate_2d(real_centers, n_mask_below, 'below')

    lo_above = extrapolate_2d(real_lo, n_mask_above, 'above')
    hi_above = extrapolate_2d(real_hi, n_mask_above, 'above')
    c_above  = extrapolate_2d(real_centers, n_mask_above, 'above')

    # Combine real and virtual
    full_lo = np.vstack((lo_below, real_lo, lo_above)) if n_mask_below or n_mask_above else real_lo
    full_hi = np.vstack((hi_below, real_hi, hi_above)) if n_mask_below or n_mask_above else real_hi
    full_centers = np.vstack((c_below, real_centers, c_above)) if n_mask_below or n_mask_above else real_centers
    
    base_id_min, base_id_max = (min(all_ap_ids), max(all_ap_ids)) if all_ap_ids else (0, 0)
    full_ids = [base_id_min - n_mask_below + i for i in range(n_mask_below)] + all_ap_ids + \
               [base_id_max + 1 + i for i in range(n_mask_above)]
    n_total = len(full_ids)
        
    # 3. Enforce valley gap strictly
    for i in range(n_total - 1):
        for xi in range(width):
            if full_hi[i, xi] >= full_lo[i+1, xi] - 1.0:
                c1 = int(round(full_centers[i, xi]))
                c2 = int(round(full_centers[i+1, xi]))
                c1 = max(0, min(c1, height - 1))
                c2 = max(0, min(c2, height - 1))
                valley = c1 + np.argmin(flat_data[c1:c2, xi]) if c2 > c1 + 2 else int(round(0.5 * (c1 + c2)))
                full_hi[i, xi] = valley - 0.5
                full_lo[i+1, xi] = valley + 0.5
                
    # 4. Draw Mask
    step3_mask = np.zeros((height, width), dtype=np.uint8)
    for i in range(n_total):
        for xi in range(width):
            if np.isfinite(full_lo[i, xi]) and np.isfinite(full_hi[i, xi]):
                r0 = max(0, int(np.ceil(full_lo[i, xi])))
                r1 = min(height, int(np.floor(full_hi[i, xi])) + 1)
                if r1 > r0:
                    step3_mask[r0:r1, xi] = 1
                    
    return step3_mask, full_lo, full_hi, full_ids


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
                             n_mask_below: int = 4,
                             n_mask_above: int = 4,
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

    if order_mask is not None:
        # Ensure order_mask is boolean (it may be loaded as int32 with order IDs)
        if order_mask.dtype != bool:
            order_mask = order_mask > 0
        logger.info("Step 3: using pre-built widened order mask for scattered light fitting.")
    else:
        logger.warning("No order_mask provided in Step 3; fitting without order mask.")
        order_mask = np.zeros_like(science_image, dtype=bool)
        
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
