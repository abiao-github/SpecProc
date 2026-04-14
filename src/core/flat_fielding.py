"""
Flat fielding stage of spectral reduction pipeline.

Handles flat field combination, order detection (tracing),
and sensitivity map extraction.
"""

import numpy as np
import logging
from pathlib import Path
from typing import List, Tuple, Dict, Optional
from scipy import ndimage
from scipy.signal import savgol_filter
from scipy.ndimage import gaussian_filter1d
from scipy.interpolate import UnivariateSpline
from src.utils.fits_io import read_fits_image, write_fits_image
from src.utils.image_processing import combine_images, normalize_flat, find_bad_pixels, estimate_background_2d
from src.core.data_structures import ApertureSet, ApertureLocation, FlatField
from src.config.config_manager import ConfigManager
from src.plotting.spectra_plotter import plot_2d_image_to_file
from src.core.grating_trace_simple import find_grating_orders_simple

logger = logging.getLogger(__name__)


class FlatFieldProcessor:
    """Handles flat field processing and order tracing."""

    def __init__(self, config: ConfigManager):
        """
        Initialize flat field processor.

        Args:
            config: Configuration manager
        """
        self.config = config
        self.flat_data = None
        self.flat_norm = None
        self.flat_sens = None
        self.flat_mask = None
        self.response_map = None
        self.scattered_light = None
        self.smoothed_model = None
        self.pixel_flat = None
        self.illumination_flat = None
        self.flat_corr_2d = None
        self.order_diagnostics = {}

    def combine_flat_frames(self, flat_filenames: List[str]) -> Tuple[np.ndarray, np.ndarray]:
        """
        Combine multiple flat frames to create master flat.

        Args:
            flat_filenames: List of flat FITS file paths

        Returns:
            Tuple of (combined_flat, flat_mask)
        """
        logger.info(f"Combining {len(flat_filenames)} flat frames...")

        if not flat_filenames:
            raise ValueError("No flat frames provided")

        # Read all flat images
        flat_images = []
        flat_exptimes = []
        for filename in flat_filenames:
            try:
                img, header = read_fits_image(filename)
                flat_images.append(img)

                exptime = np.nan
                if header is not None:
                    exptime = header.get('EXPTIME', np.nan)
                flat_exptimes.append(float(exptime) if exptime is not None else np.nan)
            except Exception as e:
                logger.warning(f"Failed to read flat frame {filename}: {e}")
                continue

        if not flat_images:
            raise RuntimeError("Could not read any flat frames")

        # Convert to array
        images_array = np.array(flat_images, dtype=np.float32)
        logger.info(f"Successfully read {len(images_array)} flat frames")

        # If flat exposures differ a lot, scale to a common exposure level first.
        # This avoids biasing the combined flat toward longest-exposure frames.
        exptime_array = np.array(flat_exptimes[:len(images_array)], dtype=np.float64)
        if len(exptime_array) == len(images_array):
            valid = np.isfinite(exptime_array) & (exptime_array > 0)
            if np.all(valid) and len(exptime_array) > 1:
                exptime_ratio = np.max(exptime_array) / np.min(exptime_array)
                if exptime_ratio > 1.1:
                    ref_exptime = float(np.median(exptime_array))
                    images_array = images_array / exptime_array[:, None, None] * ref_exptime
                    logger.info(
                        f"Flat EXPTIME range is large (min={np.min(exptime_array):.2f}, "
                        f"max={np.max(exptime_array):.2f}); scaled to reference {ref_exptime:.2f}s"
                    )

        # Combine (default workflow: median + sigma clipping for cosmic/bad pixels)
        method = self.config.get('reduce.flat', 'combine_method', 'median')
        if method == 'median':
            from astropy.stats import sigma_clip
            sigma = self.config.get_float('reduce.bias', 'combine_sigma', 3.0)
            clipped = sigma_clip(images_array, sigma=sigma, axis=0, masked=True)
            combined = np.ma.median(clipped, axis=0).filled(np.nan)
            # Fallback for any fully masked pixels
            nanpix = ~np.isfinite(combined)
            if np.any(nanpix):
                combined[nanpix] = np.median(images_array, axis=0)[nanpix]
            combined = combined.astype(np.float32)
        else:
            combined, _ = combine_images(images_array, method=method)

        # Detect bad pixels
        mosaic_maxcount = self.config.get_float('reduce.flat', 'mosaic_maxcount', 65535)
        bad_mask = combined > mosaic_maxcount

        self.flat_data = combined
        self.flat_mask = bad_mask.astype(np.int16)

        logger.info(f"Master flat created: shape {combined.shape}, "
                   f"saturated pixels: {np.sum(bad_mask)}")

        return combined, self.flat_mask

    def normalize_flat(self) -> np.ndarray:
        """
        Normalize flat field to unity.

        Returns:
            Normalized flat field
        """
        if self.flat_data is None:
            raise RuntimeError("No flat data to normalize")

        # Create normalized version for order tracing
        self.flat_norm = normalize_flat(self.flat_data, self.flat_mask)

        logger.info(f"Normalized flat field: min={np.min(self.flat_norm):.3f}, "
                   f"max={np.max(self.flat_norm):.3f}")

        return self.flat_norm

    def extract_sensitivity_map(self) -> np.ndarray:
        """
        Extract pixel-by-pixel sensitivity map from flat field.

        Returns:
            Sensitivity map (2D array)
        """
        if self.flat_data is None:
            raise RuntimeError("No flat data")

        # Simple method: normalize by smoothed version
        from scipy.ndimage import median_filter
        smoothed = median_filter(self.flat_data, size=31)
        sensitivity = self.flat_data / (smoothed + 1e-10)

        # Clip extreme values
        sensitivity = np.clip(sensitivity, 0.5, 2.0)

        self.flat_sens = sensitivity.astype(np.float32)

        logger.info(f"Sensitivity map extracted: shape {sensitivity.shape}")

        return self.flat_sens

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

        fit = np.where(fit > 0, fit, np.nanmedian(fit[fit > 0]) if np.any(fit > 0) else 1.0)
        return fit

    def _build_order_mask(self, apertures: ApertureSet, margin_pixels: int = 3) -> np.ndarray:
        """Build order-region mask from traced lower/upper edges plus a wing margin."""
        if self.flat_data is None:
            raise RuntimeError("No flat data")

        height, width = self.flat_data.shape
        mask = np.zeros((height, width), dtype=bool)
        x_coords = np.arange(width)

        for aperture in apertures.apertures.values():
            lower = aperture.get_lower(x_coords)
            upper = aperture.get_upper(x_coords)
            for x in x_coords:
                y_lo = min(lower[x], upper[x]) - margin_pixels
                y_hi = max(lower[x], upper[x]) + margin_pixels
                y1 = max(0, int(np.floor(y_lo)))
                y2 = min(height, int(np.ceil(y_hi + 1)))
                if y2 > y1:
                    mask[y1:y2, x] = True
        return mask

    def estimate_scattered_light(self, apertures: ApertureSet, poly_order: int = 3,
                                 method: str = 'chebyshev',
                                 margin_pixels: int = 3,
                                 smooth_sigma: float = 20.0,
                                 sigma_clip: float = 3.0,
                                 maxiters: int = 4,
                                 bspline_smooth: float = 1.0) -> np.ndarray:
        """
        Estimate scattered light from inter-order pixels and fit smooth 2D surface.

        This follows the common echelle strategy used by IRAF/PypeIt/REDUCE:
        use inter-order regions as constraints and fit a low-order 2D background.
        """
        if self.flat_data is None:
            raise RuntimeError("No flat data")

        order_mask = self._build_order_mask(apertures, margin_pixels=margin_pixels)
        interorder = ~order_mask

        # Mask order regions so only inter-order pixels are sampled.
        work = self.flat_data.astype(np.float64).copy()
        work[order_mask] = np.nan

        scattered = estimate_background_2d(
            work,
            order=poly_order,
            valid_mask=interorder,
            method=method,
            smooth_sigma=smooth_sigma,
            sigma_clip=sigma_clip,
            maxiters=maxiters,
            bspline_smooth=bspline_smooth,
        )
        scattered = np.clip(scattered, 0.0, np.percentile(self.flat_data, 99.9)).astype(np.float32)

        self.scattered_light = scattered
        logger.info("Estimated scattered light surface from inter-order pixels")
        return scattered

    def build_order_response_map(self, apertures: ApertureSet,
                                 source_image: Optional[np.ndarray] = None
                                 ) -> Tuple[np.ndarray, Dict[int, np.ndarray], Dict[int, np.ndarray], np.ndarray, np.ndarray, Optional[np.ndarray]]:
        """
        Build per-order response model and convert it to a 2D sensitivity map.

        This follows the echelle flat-field philosophy used in IRAF/PypeIt/REDUCE:
        estimate order-by-order blaze/response first, then construct a 2D correction map.
        """
        if self.flat_data is None and source_image is None:
            raise RuntimeError("No flat data")

        image = source_image if source_image is not None else self.flat_data
        image = np.asarray(image, dtype=np.float64)
        height, width = image.shape
        x_coords = np.arange(width)

        blaze_smooth_method = self.config.get('reduce.flat', 'blaze_smooth_method', 'savgol').strip().lower()
        blaze_smooth_window = self.config.get_int('reduce.flat', 'blaze_smooth_window', 21)
        blaze_bspline_smooth = self.config.get_float('reduce.flat', 'blaze_bspline_smooth', 0.5)
        width_smooth_window = self.config.get_int('reduce.flat', 'width_smooth_window', 41)
        profile_bin_step = self.config.get_float('reduce.flat', 'profile_bin_step', 0.01)
        pixel_flat_min = self.config.get_float('reduce.flat', 'pixel_flat_min', 0.5)
        pixel_flat_max = self.config.get_float('reduce.flat', 'pixel_flat_max', 1.5)

        response_map = np.ones_like(image, dtype=np.float32)
        coverage = np.zeros_like(image, dtype=bool)

        blaze_profiles: Dict[int, np.ndarray] = {}
        cross_profiles: Dict[int, np.ndarray] = {}
        order_diag: Dict[int, Dict[str, np.ndarray]] = {}

        for aperture_id, aperture in apertures.apertures.items():
            center = aperture.get_position(x_coords)
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

            blaze_1d = self._smooth_order_response(
                raw_blaze,
                method=blaze_smooth_method,
                window=blaze_smooth_window,
                bspline_smooth=blaze_bspline_smooth,
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
            sample_cols = set()
            if valid_cols.size > 0:
                sample_cols = set(np.unique(np.linspace(valid_cols.min(), valid_cols.max(), num=min(6, valid_cols.size), dtype=int)))

            # Build universal cross-order profile in normalized coordinates y_norm.
            yn_all = []
            fn_all = []
            for x in x_coords:
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
                    yn_all.append(yn[good])
                    fn_all.append(fn[good])

            if not yn_all:
                logger.warning(f"Order {aperture_id}: insufficient normalized points, fallback to simple model")
                continue

            yn_all = np.concatenate(yn_all)
            fn_all = np.concatenate(fn_all)

            step = max(0.002, float(profile_bin_step))
            yn_min = max(-1.2, float(np.nanpercentile(yn_all, 1)))
            yn_max = min(1.2, float(np.nanpercentile(yn_all, 99)))
            if yn_max - yn_min < 0.2:
                yn_min, yn_max = -0.8, 0.8
            bins = np.arange(yn_min, yn_max + step, step)
            if bins.size < 8:
                bins = np.linspace(-0.8, 0.8, 81)

            idx = np.digitize(yn_all, bins) - 1
            prof_x = []
            prof_y = []
            for i in range(len(bins) - 1):
                pick = idx == i
                if np.sum(pick) < 6:
                    continue
                prof_x.append(0.5 * (bins[i] + bins[i + 1]))
                prof_y.append(np.nanmedian(fn_all[pick]))

            if len(prof_x) < 8:
                logger.warning(f"Order {aperture_id}: insufficient binned profile points, fallback to simple model")
                continue

            prof_x = np.asarray(prof_x, dtype=float)
            prof_y = np.asarray(prof_y, dtype=float)
            prof_y = np.where(np.isfinite(prof_y), prof_y, 1.0)
            prof_y = np.clip(prof_y, 1e-6, None)

            # Smooth the empirical profile for high-SNR model reconstruction.
            prof_y = self._smooth_order_response(
                prof_y,
                method='savgol',
                window=max(11, int(11 * max(1.0, step / 0.01))),
                bspline_smooth=0.0,
            )
            csum = np.nansum(prof_y)
            if np.isfinite(csum) and csum > 0:
                cross_profiles[aperture_id] = (prof_y / csum).astype(np.float32)

            sampled_realspace = []

            # Reconstruct 2D model for this order: Model(X,Y)=Profile(y_norm)*Blaze(X).
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
                pval = np.interp(yreq, prof_x, prof_y, left=1.0, right=1.0)
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
                'smooth_blaze': blaze_1d.astype(np.float32),
                'width_smooth': width_smooth.astype(np.float32),
                'profile_x': prof_x.astype(np.float32),
                'profile_y': prof_y.astype(np.float32),
                'sampled_realspace': sampled_realspace,
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
        pixel_flat = np.clip(pixel_flat, 0.1, 10.0).astype(np.float32)

        # Illumination flat is intentionally not used/output in pixel-to-pixel-only workflow.
        illumination_flat = None

        # Final 2D correction map equals pixel flat in this workflow.
        flat_corr_2d = np.clip(pixel_flat, 0.1, 10.0).astype(np.float32)

        sens = flat_corr_2d.copy()

        # Normalize sensitivity map to unity on valid pixels.
        if np.any(coverage):
            med = np.median(sens[coverage][np.isfinite(sens[coverage])])
            if np.isfinite(med) and med > 0:
                sens[coverage] /= med

        sens = np.clip(sens, 0.25, 4.0).astype(np.float32)

        self.response_map = response_map
        self.smoothed_model = smoothed_model
        self.pixel_flat = pixel_flat
        self.illumination_flat = illumination_flat
        self.flat_corr_2d = flat_corr_2d
        self.flat_sens = sens
        self.order_diagnostics = order_diag
        logger.info("Built Step4 flat model via 1D blaze + normalized-profile reconstruction")
        return sens, blaze_profiles, cross_profiles, smoothed_model, pixel_flat, illumination_flat

    def save_step4_diagnostics(self, output_dir: Path,
                               clean_flat: Optional[np.ndarray] = None,
                               clean_flat_corrected: Optional[np.ndarray] = None):
        """Save step4 diagnostic artifacts: blaze/profile/inverse-map/model/pixel-flat and corrected flat."""
        out_dir = Path(output_dir)
        out_dir.mkdir(parents=True, exist_ok=True)

        fig_format = self.config.get('reduce', 'fig_format', 'png')
        save_plots = self.config.get_bool('reduce', 'save_plots', True)

        if self.response_map is not None:
            write_fits_image(str(out_dir / 'model_flat_2d.fits'), self.response_map.astype(np.float32), dtype='float32')
        if self.pixel_flat is not None:
            write_fits_image(str(out_dir / 'pixel_flat_2d.fits'), self.pixel_flat.astype(np.float32), dtype='float32')
        if self.flat_sens is not None:
            write_fits_image(str(out_dir / 'sensitivity_map_2d.fits'), self.flat_sens.astype(np.float32), dtype='float32')
        if clean_flat is not None:
            write_fits_image(str(out_dir / 'clean_master_flat.fits'), clean_flat.astype(np.float32), dtype='float32')
        if clean_flat_corrected is not None:
            write_fits_image(str(out_dir / 'clean_master_flat_flat2d_corrected.fits'), clean_flat_corrected.astype(np.float32), dtype='float32')

        if save_plots:
            if self.response_map is not None:
                plot_2d_image_to_file(self.response_map, str(out_dir / f'model_flat_2d.{fig_format}'), 'Reconstructed 2D Flat Model')
            if self.pixel_flat is not None:
                plot_2d_image_to_file(self.pixel_flat, str(out_dir / f'pixel_flat_2d.{fig_format}'), 'Pixel Flat (Sensitivity Map)')
            if clean_flat is not None:
                plot_2d_image_to_file(clean_flat, str(out_dir / f'clean_master_flat.{fig_format}'), 'Clean Master Flat')
            if clean_flat_corrected is not None:
                plot_2d_image_to_file(clean_flat_corrected, str(out_dir / f'clean_master_flat_flat2d_corrected.{fig_format}'), 'Clean Master Flat After 2D Flat Correction')

        if not self.order_diagnostics:
            return

        import matplotlib.pyplot as plt

        blaze_dir = out_dir / 'blaze_per_order'
        profile_dir = out_dir / 'universal_profile_per_order'
        realspace_dir = out_dir / 'realspace_profile_per_order'
        blaze_dir.mkdir(parents=True, exist_ok=True)
        profile_dir.mkdir(parents=True, exist_ok=True)
        realspace_dir.mkdir(parents=True, exist_ok=True)

        for aperture_id, diag in self.order_diagnostics.items():
            raw_blaze = diag.get('raw_blaze')
            smooth_blaze = diag.get('smooth_blaze')
            profile_x = diag.get('profile_x')
            profile_y = diag.get('profile_y')
            sampled_realspace = diag.get('sampled_realspace', [])

            np.savez(
                blaze_dir / f'order_{aperture_id:03d}_blaze.npz',
                raw_blaze=raw_blaze,
                smooth_blaze=smooth_blaze,
                width_smooth=diag.get('width_smooth'),
            )
            np.savez(
                profile_dir / f'order_{aperture_id:03d}_universal_profile.npz',
                y_norm=profile_x,
                profile=profile_y,
            )

            if save_plots and raw_blaze is not None and smooth_blaze is not None:
                plt.figure(figsize=(10, 4))
                x = np.arange(raw_blaze.size)
                plt.plot(x, raw_blaze, color='0.6', linewidth=0.8, label='Raw Blaze 1D')
                plt.plot(x, smooth_blaze, color='tab:red', linewidth=1.2, label='Smooth Blaze 1D')
                plt.xlabel('X (column)')
                plt.ylabel('Flux')
                plt.title(f'Order {aperture_id} Blaze Function')
                plt.legend()
                plt.grid(alpha=0.2)
                plt.tight_layout()
                plt.savefig(blaze_dir / f'order_{aperture_id:03d}_blaze.{fig_format}', dpi=150)
                plt.close()

            if save_plots and profile_x is not None and profile_y is not None:
                plt.figure(figsize=(7, 4))
                plt.plot(profile_x, profile_y, color='tab:blue', linewidth=1.4)
                plt.xlabel('y_norm')
                plt.ylabel('Normalized Flux')
                plt.title(f'Order {aperture_id} Universal 1D Profile')
                plt.grid(alpha=0.2)
                plt.tight_layout()
                plt.savefig(profile_dir / f'order_{aperture_id:03d}_universal_profile.{fig_format}', dpi=150)
                plt.close()

            if save_plots and sampled_realspace:
                plt.figure(figsize=(8, 5))
                for prof in sampled_realspace:
                    y = prof['y']
                    m = prof['model_flux']
                    xcol = prof['x']
                    plt.plot(y, m, linewidth=1.0, label=f'X={xcol}')
                plt.xlabel('Y (physical pixel)')
                plt.ylabel('Model Flux')
                plt.title(f'Order {aperture_id} Inverse-Mapped Real-Space Profiles')
                if len(sampled_realspace) <= 8:
                    plt.legend(fontsize=8)
                plt.grid(alpha=0.2)
                plt.tight_layout()
                plt.savefig(realspace_dir / f'order_{aperture_id:03d}_realspace_profiles.{fig_format}', dpi=150)
                plt.close()

    def detect_orders(self, threshold: float = 0.3, trace_degree: int = 3) -> ApertureSet:
        """
        Detect grating orders in flat field with simple tracing.

        For grating spectrographs, orders run horizontally (left to right) across
        the detector. The algorithm collapses the image along the dispersion direction
        (X) to find initial order centroids along Y, then traces each order across
        the image using polynomial fitting.

        Args:
            threshold: Detection threshold (relative to max on collapsed profile)
            trace_degree: Polynomial degree used to fit order trace.

        Returns:
            ApertureSet with detected orders
        """
        if self.flat_data is None:
            raise RuntimeError("No flat data")

        logger.info("Detecting grating orders (simple tracing)...")

        threshold = self.config.get_float('reduce.trace', 'seed_threshold', threshold)
        trace_degree = self.config.get_int('reduce.trace', 'degree', trace_degree)

        # Use the simplified grating order tracing method
        apertures = find_grating_orders_simple(
            data=self.flat_data,
            mask=self.flat_mask,
            config=self.config,
            threshold=threshold,
            trace_degree=trace_degree,
            separation=self.config.get_float('reduce.trace', 'separation', 30.0)
        )

        return apertures

    def extract_blaze_profiles(self, apertures: ApertureSet) -> Dict[int, np.ndarray]:
        """
        Extract blaze profiles for each order from the flat field.

        The blaze profile represents the wavelength-dependent efficiency
        of the grating along each order.

        Args:
            apertures: Detected aperture set

        Returns:
            Dictionary of blaze profiles per order
        """
        if self.flat_norm is None:
            self.normalize_flat()

        logger.info("Extracting blaze profiles for de-blazing...")

        blaze_profiles = {}
        height, width = self.flat_norm.shape

        for aperture_id, aperture in apertures.apertures.items():
            # Evaluate trace center directly as y(x).
            x_coords = np.arange(width)
            center_pos = aperture.get_position(x_coords)

            # Extract profile along the order
            profile = np.zeros(width)
            for x in range(width):
                y_center = center_pos[x]
                y_min = max(0, int(y_center - 5))
                y_max = min(height, int(y_center + 6))

                if y_max > y_min:
                    profile[x] = np.mean(self.flat_norm[y_min:y_max, x])

            # Smooth the profile to reduce noise
            from scipy.ndimage import median_filter
            smoothed_profile = median_filter(profile, size=21)

            # Normalize to avoid division by zero
            smoothed_profile = np.where(smoothed_profile > 0, smoothed_profile, 1.0)

            blaze_profiles[aperture_id] = smoothed_profile.astype(np.float32)

        logger.info(f"Extracted blaze profiles for {len(blaze_profiles)} orders")
        return blaze_profiles

    def _find_peaks(self, array: np.ndarray, min_height: float = 0.0,
                   distance: int = 10, prominence: float = 0.0) -> list:
        """Find local maxima in 1D array."""
        from scipy.signal import find_peaks
        peaks, _ = find_peaks(array, height=min_height, distance=distance, prominence=prominence)
        return peaks.tolist()

    def save_flat_field(self, output_path: str, config: ConfigManager):
        """Save flat field data to FITS file."""
        if self.flat_data is None:
            raise RuntimeError("No flat data to save")

        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        from astropy.io import fits

        # Create multi-HDU FITS with processing headers
        primary_hdu = fits.PrimaryHDU(data=self.flat_data.astype(np.float32))
        primary_hdu.header['FLATCR'] = (True, 'Flat fielding completed')

        hdul = fits.HDUList([
            primary_hdu,
            fits.ImageHDU(data=self.flat_mask, name='MASK'),
            fits.ImageHDU(data=self.flat_norm.astype(np.float32), name='NORM'),
        ])

        if self.flat_sens is not None:
            hdul.append(fits.ImageHDU(data=self.flat_sens, name='SENSITIVITY'))
        if self.response_map is not None:
            hdul.append(fits.ImageHDU(data=self.response_map.astype(np.float32), name='RESPONSE'))
        if self.scattered_light is not None:
            hdul.append(fits.ImageHDU(data=self.scattered_light.astype(np.float32), name='SCATTER'))
        if self.smoothed_model is not None:
            hdul.append(fits.ImageHDU(data=self.smoothed_model.astype(np.float32), name='SMOOTH_MODEL'))
        if self.pixel_flat is not None:
            hdul.append(fits.ImageHDU(data=self.pixel_flat.astype(np.float32), name='PIXEL_FLAT'))
        if self.flat_corr_2d is not None:
            hdul.append(fits.ImageHDU(data=self.flat_corr_2d.astype(np.float32), name='FLAT_CORR_2D'))

        hdul.writeto(str(output_path), overwrite=True)
        logger.info(f"Saved flat field to {output_path}")

        # Save diagnostic plots if enabled
        save_plots = config.get_bool('reduce', 'save_plots', True)
        if save_plots:
            out_dir = Path(output_path).parent
            fig_format = config.get('reduce', 'fig_format', 'png')
            # Plot master flat
            plot_2d_image_to_file(self.flat_data, str(out_dir / f'master_flat.{fig_format}'), "Master Flat Field")


def _fit_background_envelope_1d(cross_section, order_centers_y, half_win=4):
    """Fit a smooth background baseline envelope using inter-order valley samples.

    At each inter-order midpoint (and detector edges) takes the local minimum inside a
    ±half_win pixel window as a valley sample, then fits a smoothing cubic spline through
    all samples to produce a continuous background baseline across the whole column.

    Parameters
    ----------
    cross_section   : 1-D array – the CCD column cross-section
    order_centers_y : list[float] – center pixel of each order (sorted ascending)
    half_win        : int – half-width of local-minimum window at each sample point

    Returns
    -------
    envelope : np.ndarray, shape (n,) – background baseline at every pixel
    """
    from scipy.interpolate import UnivariateSpline

    n  = len(cross_section)
    cs = cross_section.astype(float)
    centers = sorted(float(c) for c in order_centers_y)

    def _local_min(pos):
        i0 = max(0, int(round(pos)) - half_win)
        i1 = min(n, int(round(pos)) + half_win + 1)
        seg = cs[i0:i1]
        return float(np.nanmin(seg)) if seg.size > 0 else float(cs[int(np.clip(round(pos), 0, n-1))])

    sample_y = []
    sample_v = []

    # Lower detector edge
    gap0 = (centers[1] - centers[0]) if len(centers) > 1 else 20.0
    edge_lo = max(0.0, centers[0] - gap0 * 0.6)
    sample_y.append(edge_lo);  sample_v.append(_local_min(edge_lo))

    # Midpoints between consecutive orders
    for k in range(len(centers) - 1):
        mid = (centers[k] + centers[k + 1]) * 0.5
        sample_y.append(mid);  sample_v.append(_local_min(mid))

    # Upper detector edge
    gap1 = (centers[-1] - centers[-2]) if len(centers) > 1 else 20.0
    edge_hi = min(float(n - 1), centers[-1] + gap1 * 0.6)
    sample_y.append(edge_hi);  sample_v.append(_local_min(edge_hi))

    sy = np.array(sample_y, dtype=float)
    sv = np.array(sample_v, dtype=float)

    # Deduplicate
    _, ui = np.unique(sy, return_index=True)
    sy, sv = sy[ui], sv[ui]

    x_all = np.arange(n, dtype=float)
    if len(sy) >= 4:
        spl = UnivariateSpline(sy, sv, k=3, s=0, ext=3)
        envelope = spl(x_all)
    elif len(sy) >= 2:
        envelope = np.interp(x_all, sy, sv)
    else:
        envelope = np.full(n, sv[0] if sv.size > 0 else float(np.nanmedian(cs)))

    # Clamp: envelope must not exceed the actual cross-section
    return np.clip(envelope, np.nanmin(cs), np.nanmax(cs))


def _find_profile_roots_1d(cross_section, center_idx, search_lo, search_hi,
                           threshold_frac=0.05, bg_envelope=None,
                           noise_floor=None):
    """Find where an order profile's wings blend into the background envelope.

    Walks outward from *center_idx*; stops when ``cs[j] - bg_envelope[j]`` drops to
    ``threshold_frac * (peak - bg_at_peak)``.  Position is subpixel-interpolated.

    Parameters
    ----------
    cross_section  : 1-D array
    center_idx     : int
    search_lo      : float – lower search limit (inter-order midpoint)
    search_hi      : float – upper search limit
    threshold_frac : float – fraction of peak amplitude above envelope for the root
    bg_envelope    : 1-D array or None – smooth background baseline; if None a flat
                     local-percentile estimate is used (legacy behaviour)
    noise_floor    : float or None – minimum excess above envelope in counts

    Returns
    -------
    root_lo, root_hi : float – subpixel crossing positions
    bg_at_peak       : float – envelope value at the order center
    threshold        : float – absolute intensity level of the crossing
    """
    n   = len(cross_section)
    cs  = cross_section.astype(float)
    ci  = int(np.clip(center_idx,  0, n - 1))
    slo = int(np.clip(float(search_lo), 0, n - 1))
    shi = int(np.clip(float(search_hi), 0, n - 1))

    peak = cs[ci]

    if bg_envelope is None:
        region = cs[slo:shi + 1]
        flat_bg = float(np.nanpercentile(region, 10)) if region.size > 3 else float(np.nanmin(region) if region.size else 0.0)
        bg_envelope = np.full(n, flat_bg)

    bg_at_peak    = float(bg_envelope[ci])
    amp_at_peak   = max(float(peak) - bg_at_peak, 1.0)
    excess_thresh = threshold_frac * amp_at_peak

    # Hybrid threshold: percentage of peak OR noise floor, whichever is larger.
    if noise_floor is not None and np.isfinite(noise_floor):
        excess_thresh = max(excess_thresh, float(noise_floor))

    # Keep threshold physically meaningful for this order.
    excess_thresh = float(np.clip(excess_thresh, 0.0, 0.9 * amp_at_peak))

    def _excess(j):
        return cs[j] - float(bg_envelope[int(np.clip(j, 0, n - 1))])

    # --- lower crossing: walk ci → slo ---
    root_lo = float(slo)
    for j in range(ci, slo - 1, -1):
        if _excess(j) <= excess_thresh:
            if j + 1 <= ci:
                e0, e1 = _excess(j), _excess(j + 1)   # e0 ≤ thresh < e1
                denom = e1 - e0
                t = float(np.clip((excess_thresh - e0) / denom, 0.0, 1.0)) if denom > 0 else 0.0
                root_lo = float(j) + t
            else:
                root_lo = float(j)
            break

    # --- upper crossing: walk ci → shi ---
    root_hi = float(shi)
    for j in range(ci, shi + 1):
        if _excess(j) <= excess_thresh:
            if j - 1 >= ci:
                e0, e1 = _excess(j - 1), _excess(j)   # e0 > thresh ≥ e1
                denom = e0 - e1
                t = float(np.clip((e0 - excess_thresh) / denom, 0.0, 1.0)) if denom > 0 else 0.0
                root_hi = float(j - 1) + t
            else:
                root_hi = float(j)
            break

    threshold = bg_at_peak + excess_thresh
    return root_lo, root_hi, bg_at_peak, threshold


def process_flat_stage(config: ConfigManager, flat_filenames: List[str],
                      midpath: str = './midpath') -> Tuple[FlatField, ApertureSet]:
    """
    Execute flat fielding stage.

    Args:
        config: Configuration manager
        flat_filenames: List of flat frame paths
        midpath: Intermediate output directory

    Returns:
        Tuple of (FlatField, ApertureSet)
    """
    processor = FlatFieldProcessor(config)

    # Combine flat frames
    flat_data, flat_mask = processor.combine_flat_frames(flat_filenames)

    # Normalize
    logger.info("Normalizing flat field...")
    flat_norm = processor.normalize_flat()
    logger.info("Flat field normalization completed")

    # Save master flat field as FITS and PNG
    logger.info("Saving master flat field...")
    base_output_path = config.get_output_path()
    flat_file = Path(base_output_path) / 'step2_trace' / 'master_flat.fits'
    flat_file.parent.mkdir(parents=True, exist_ok=True)
    processor.save_flat_field(str(flat_file), config)
    logger.info(f"Saved master flat field to {flat_file}")

    # Keep normalized flat in memory for internal modeling only.
    from src.plotting.spectra_plotter import plot_2d_image_to_file
    base_output_path = config.get_output_path()

    # Detect orders
    logger.info("Detecting grating orders...")
    apertures = processor.detect_orders(threshold=0.3)
    logger.info(f"Order detection completed: {apertures.norders} orders found")

    # Step 2 only traces orders and builds an initial flat model.
    # Step 3 will generate and subtract scattered-light models for master flat and science.
    logger.info("Building initial flat model (scattered-light subtraction deferred to Step 3)...")
    flat_sens, blaze_profiles, cross_profiles, smoothed_model, pixel_flat, illumination_flat = (
        processor.build_order_response_map(apertures, source_image=flat_data)
    )
    logger.info("Initial order-wise flat/blaze modeling completed")

    # Rewrite master flat to include newly derived SENSITIVITY/RESPONSE extensions.
    processor.save_flat_field(str(flat_file), config)
    logger.info("Updated master flat FITS with response/sensitivity extensions")

    # Save aperture information
    import json
    apertures_file = Path(base_output_path) / 'step2_trace' / 'apertures.json'
    with open(apertures_file, 'w') as f:
        # Convert aperture data to serializable format
        apertures_dict = {
            'norders': apertures.norders,
            'direct_axis': apertures.direct_axis,
            'apertures': {}
        }
        for aperture_id, aperture in apertures.apertures.items():
            apertures_dict['apertures'][aperture_id] = {
                'aperture': aperture.aperture,
                'order': aperture.order,
                'center_coef': aperture.center_coef.tolist(),
                'lower_coef': aperture.lower_coef.tolist(),
                'upper_coef': aperture.upper_coef.tolist(),
                'width': aperture.width
            }
        json.dump(apertures_dict, f, indent=2)
    logger.info(f"Saved aperture information to {apertures_file}")

    # Step2 should only emit tracing diagnostics.
    # Step4 writes flat-model and 2D flat-correction artifacts.
    save_plots = config.get_bool('reduce', 'save_plots', True)
    if save_plots:
        logger.info("Generating diagnostic plots...")
        out_dir = Path(base_output_path) / 'step2_trace'
        fig_format = config.get('reduce', 'fig_format', 'png')
        import matplotlib.pyplot as plt

        # Plot order traces on combined master flat (no normalized output dependency)
        logger.info("Plotting order traces on combined master flat...")
        plt.figure(figsize=(12, 10))
        display_flat = flat_data
        vmin = np.percentile(display_flat, 1)
        vmax = np.percentile(display_flat, 99)
        plt.imshow(display_flat, aspect='auto', origin='lower', cmap='viridis', vmin=vmin, vmax=vmax)
        plt.colorbar(label='Counts')

        # Plot order traces (for grating spectrographs, orders run left to right)
        height, width = flat_data.shape
        x_coords = np.arange(width)
        center_col = width // 2
        q25_col = int(np.clip(round(width * 0.25), 0, width - 1))
        q50_col = int(np.clip(round(width * 0.50), 0, width - 1))
        q75_col = int(np.clip(round(width * 0.75), 0, width - 1))
        diag_cols = [q25_col, q50_col, q75_col]
        root_frac = config.get_float('reduce.trace', 'aperture_root_fraction', 0.03)
        noise_floor_sigma = config.get_float('reduce.trace', 'aperture_noise_floor_sigma', 3.0)

        def _estimate_noise_floor(cross_section_1d, bg_envelope_1d, centers_y):
            valley_mask = np.zeros(height, dtype=bool)
            for k in range(len(centers_y) - 1):
                c0, c1 = centers_y[k], centers_y[k + 1]
                mid = 0.5 * (c0 + c1)
                half = max(1.0, 0.18 * abs(c1 - c0))
                i0 = int(np.clip(np.floor(mid - half), 0, height - 1))
                i1 = int(np.clip(np.ceil(mid + half), 0, height - 1))
                valley_mask[i0:i1 + 1] = True

            residual = cross_section_1d - bg_envelope_1d
            valley_residual = residual[valley_mask] if np.any(valley_mask) else residual
            if valley_residual.size >= 5:
                med = float(np.nanmedian(valley_residual))
                mad = float(np.nanmedian(np.abs(valley_residual - med)))
                sigma_bg = max(1.4826 * mad, 1e-6)
            else:
                sigma_bg = max(float(np.nanstd(residual)), 1e-6)
            return max(noise_floor_sigma * sigma_bg, 0.0)

        ordered_apertures = sorted(apertures.apertures.items(), key=lambda x: x[0])
        valid_orders = []
        for order_idx, (aperture_id, aperture) in enumerate(ordered_apertures, start=1):
            center_pos = aperture.get_position(x_coords)

            # Keep order_traces numbering/selection aligned with apertures diagnostics:
            # ignore edge-only partial orders whose center is out of frame at center column.
            center_mid = float(aperture.get_position(center_col))
            if center_mid < -5 or center_mid > height + 5:
                continue

            # Skip apertures where the center trajectory is mostly outside the image
            # (prevents phantom boundary lines with no visible center trace).
            n_in_frame = int(np.sum((center_pos >= 0) & (center_pos < height)))
            if n_in_frame < 0.15 * width:
                continue

            valid_center = np.isfinite(center_pos) & (center_pos >= 0) & (center_pos < height)
            if np.sum(valid_center) < max(10, int(0.10 * width)):
                continue

            plt.plot(x_coords[valid_center], center_pos[valid_center], 'r-', linewidth=1, label=f'Order {order_idx}')

            # Order index annotation, using the same numbering rule as apertures.png
            x_valid = x_coords[valid_center]
            y_valid = center_pos[valid_center]
            if x_valid.size > 0:
                x_anno = int(x_valid[-1])
                y_anno = float(y_valid[-1])
                plt.text(x_anno + 3, y_anno,
                         str(order_idx),
                         color='red', fontsize=5, ha='left', va='center')

            valid_orders.append({
                'order_idx': order_idx,
                'aperture': aperture,
                'valid_center': valid_center,
            })

        # Compute/order boundaries using the SAME root-crossing rule as apertures diagnostics,
        # and draw boundaries only for already-accepted valid orders.
        if len(valid_orders) >= 1:
            y_lo_all = np.full((len(valid_orders), width), np.nan, dtype=float)
            y_hi_all = np.full((len(valid_orders), width), np.nan, dtype=float)

            for xi, x in enumerate(x_coords):
                cross_section_col = flat_data[:, x].astype(float)

                centers = []
                active_rows = []
                for row_idx, item in enumerate(valid_orders):
                    cy = float(item['aperture'].get_position(x))
                    if not np.isfinite(cy) or cy < -5 or cy > height + 5:
                        continue
                    centers.append(float(np.clip(cy, 0, height - 1)))
                    active_rows.append((row_idx, cy))

                if len(centers) == 0:
                    continue

                bg_env_col = _fit_background_envelope_1d(cross_section_col, centers)
                noise_floor_col = _estimate_noise_floor(cross_section_col, bg_env_col, centers)

                for local_i, (row_idx, cy_raw) in enumerate(active_rows):
                    cy_y = float(np.clip(cy_raw, 0, height - 1))
                    cy_idx = int(round(cy_y))

                    if local_i > 0:
                        search_lo = 0.5 * (centers[local_i] + centers[local_i - 1])
                    else:
                        gap = (centers[1] - centers[0]) if len(centers) > 1 else 20.0
                        search_lo = max(0.0, centers[local_i] - 0.8 * gap)

                    if local_i < len(centers) - 1:
                        search_hi = 0.5 * (centers[local_i] + centers[local_i + 1])
                    else:
                        gap = (centers[-1] - centers[-2]) if len(centers) > 1 else 20.0
                        search_hi = min(float(height - 1), centers[local_i] + 0.8 * gap)

                    root_lo, root_hi, _bg_at_peak, _threshold = _find_profile_roots_1d(
                        cross_section_col, cy_idx, search_lo, search_hi,
                        root_frac, bg_env_col,
                        noise_floor=noise_floor_col)

                    y_lo_all[row_idx, xi] = float(np.clip(root_lo, 0, height - 1))
                    y_hi_all[row_idx, xi] = float(np.clip(root_hi, 0, height - 1))

            for row_idx, item in enumerate(valid_orders):
                valid_center = item['valid_center']
                y_lo = y_lo_all[row_idx, :]
                y_hi = y_hi_all[row_idx, :]
                valid_bounds = valid_center & np.isfinite(y_lo) & np.isfinite(y_hi)
                if np.sum(valid_bounds) < max(10, int(0.10 * width)):
                    continue
                plt.plot(x_coords[valid_bounds], y_lo[valid_bounds], color='white', linestyle='--', linewidth=1.1, alpha=0.55)
                plt.plot(x_coords[valid_bounds], y_hi[valid_bounds], color='white', linestyle='--', linewidth=1.1, alpha=0.55)
                plt.plot(x_coords[valid_bounds], y_lo[valid_bounds], color='#00E5FF', linestyle='--', linewidth=0.8, alpha=0.95)
                plt.plot(x_coords[valid_bounds], y_hi[valid_bounds], color='#00E5FF', linestyle='--', linewidth=0.8, alpha=0.95)

                # Mark q25/q50/q75 positions to visually cross-check against
                # apertures_q25/q50/q75 diagnostics.
                ap = item['aperture']
                for xc in diag_cols:
                    if not valid_center[xc] or not np.isfinite(y_lo[xc]) or not np.isfinite(y_hi[xc]):
                        continue
                    yc = float(np.clip(ap.get_position(xc), 0, height - 1))
                    plt.scatter([xc], [yc], s=9, c='red', marker='o', alpha=0.95, zorder=6)
                    plt.scatter([xc], [y_lo[xc]], s=8, c='#00E5FF', marker='o', alpha=0.95, zorder=6)
                    plt.scatter([xc], [y_hi[xc]], s=8, c='#00E5FF', marker='o', alpha=0.95, zorder=6)

        plt.xlabel('Pixel (X)')
        plt.ylabel('Pixel (Y)')
        plt.title('Order Traces on Combined Master Flat')

        if len(valid_orders) <= 10:  # Only show legend if not too many displayed orders
            plt.legend()
        plt.grid(False)
        traces_plot_file = out_dir / f'order_traces.{fig_format}'
        plt.savefig(str(traces_plot_file), dpi=150, bbox_inches='tight')
        plt.close()
        logger.info(f"Saved order traces plot to {traces_plot_file}")

        # Plot apertures diagnostics at multiple X columns: 25%, 50%, 75%
        logger.info("Plotting aperture cross-sections with envelope-based profile-root boundaries...")
        # Reuse the same threshold parameters in apertures diagnostics.

        def _plot_apertures_cross_section(col_idx: int, out_name: str):
            plt.figure(figsize=(10, 8))
            cross_section = flat_data[:, col_idx].astype(float)
            y_coords = np.arange(height)
            plt.plot(y_coords, cross_section, color='0.25', linewidth=0.75, label='Cross-section')
            cs_max = float(np.nanmax(cross_section))
            cs_span = max(cs_max - float(np.nanmin(cross_section)), 1e-6)

            ordered_local = sorted(apertures.apertures.items(), key=lambda x: x[0])

            # Pre-collect center positions (clipped to detector) for envelope + window computation
            order_centers_y = []
            for _ap_id, _ap in ordered_local:
                _cy = float(_ap.get_position(col_idx))
                order_centers_y.append(float(np.clip(_cy, 0, height - 1)))

            # Fit smooth background envelope through inter-order valley minima
            bg_envelope = _fit_background_envelope_1d(cross_section, order_centers_y)
            plt.plot(y_coords, bg_envelope, color='orange', linewidth=0.8,
                     linestyle='--', alpha=0.70, label='BG envelope', zorder=2)

            # Estimate background noise on inter-order regions using robust MAD.
            noise_floor = _estimate_noise_floor(cross_section, bg_envelope, order_centers_y)

            for order_idx, (aperture_id, aperture) in enumerate(ordered_local, start=1):
                center_y = float(aperture.get_position(col_idx))

                # Skip orders whose center is off the detector at this column
                if center_y < -5 or center_y > height + 5:
                    continue

                cy_y = float(np.clip(center_y, 0, height - 1))
                cy_idx = int(round(cy_y))
                i = order_idx - 1  # 0-based index in order_centers_y

                # Search window: midpoint to each neighbouring order (or detector edge)
                if i > 0:
                    search_lo = (order_centers_y[i] + order_centers_y[i - 1]) * 0.5
                else:
                    gap = (order_centers_y[1] - order_centers_y[0]) if len(order_centers_y) > 1 else 20.0
                    search_lo = max(0.0, order_centers_y[i] - gap * 0.8)

                if i < len(order_centers_y) - 1:
                    search_hi = (order_centers_y[i] + order_centers_y[i + 1]) * 0.5
                else:
                    gap = (order_centers_y[-1] - order_centers_y[-2]) if len(order_centers_y) > 1 else 20.0
                    search_hi = min(float(height - 1), order_centers_y[i] + gap * 0.8)

                root_lo, root_hi, bg_at_peak, threshold = _find_profile_roots_1d(
                    cross_section, cy_idx, search_lo, search_hi,
                    root_frac, bg_envelope,
                    noise_floor=noise_floor)

                lo_y = float(np.clip(root_lo, 0, height - 1))
                hi_y = float(np.clip(root_hi, 0, height - 1))
                peak_val = float(cross_section[cy_idx])

                # Horizontal error bar at the envelope crossing level
                plt.errorbar(
                    cy_y, threshold,
                    xerr=[[cy_y - lo_y], [hi_y - cy_y]],
                    fmt='none', color='#00E5FF',
                    capsize=2.5, capthick=0.9, elinewidth=0.9, alpha=0.90,
                    zorder=4,
                )

                # Short red vertical tick: fixed length for all orders in this panel
                tick_gap = 0.012 * cs_span
                tick_len = 0.040 * cs_span
                tick_bot = peak_val + tick_gap
                tick_top = tick_bot + tick_len
                plt.vlines(cy_y, tick_bot, tick_top,
                           colors='red', linewidth=1.1, alpha=0.95, zorder=5)

                # Order label just above the tick
                plt.text(cy_y, tick_top + 0.008 * cs_span,
                         str(order_idx),
                         color='red', fontsize=5, ha='center', va='bottom')

            plt.xlabel('Pixel (Y)')
            plt.ylabel('Counts')
            plt.title(f'Apertures Cross-section at X={col_idx}\n'
                      f'Orange = BG envelope, Cyan = max(frac*peak, noise floor), Red = center')
            if np.isfinite(noise_floor):
                plt.text(0.01, 0.985,
                         f'frac={root_frac:.3f}, noise floor={noise_floor:.2f} counts',
                         transform=plt.gca().transAxes,
                         fontsize=7, color='#00AACC', ha='left', va='top')
            plt.legend(fontsize=7)
            plt.grid(True, alpha=0.3)
            out_file = out_dir / out_name
            plt.savefig(str(out_file), dpi=150, bbox_inches='tight')
            plt.close()
            logger.info(f"Saved apertures plot to {out_file}")

        # Output apertures diagnostics at 25%, 50%, 75% columns.
        _plot_apertures_cross_section(q25_col, f'apertures_q25.{fig_format}')
        _plot_apertures_cross_section(q50_col, f'apertures_q50.{fig_format}')
        _plot_apertures_cross_section(q75_col, f'apertures_q75.{fig_format}')

    # Create FlatField object
    logger.info("Creating FlatField object...")
    flat_field = FlatField(
        flat_data=flat_data,
        flat_sens=flat_sens,
    flat_norm=flat_norm,
    flat_mask=flat_mask,
        scattered_light=None,
        smoothed_model=smoothed_model,
        pixel_flat=pixel_flat,
    illumination_flat=None,
        flat_corr_2d=processor.flat_corr_2d,
        aperture_set=apertures,
        cross_profiles=cross_profiles,
        blaze_profiles=blaze_profiles
    )

    logger.info("Step 2 order tracing completed successfully")
    return flat_field, apertures
