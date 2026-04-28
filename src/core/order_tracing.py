"""
Order tracing stage — Step 2 of the spectral reduction pipeline.

This module combines master flat frames, detects and traces echelle or grating
orders, defines their boundaries, and saves the results for subsequent
pipeline stages. The main entry point is `process_order_tracing_stage`.

The main entry point is process_order_tracing_stage().  FlatFieldProcessor is the
class that wraps all processing state.
"""

import numpy as np
import logging
from pathlib import Path
from typing import List, Tuple, Dict, Optional, Any
from scipy.signal import find_peaks
from scipy.ndimage import gaussian_filter1d
from scipy.interpolate import UnivariateSpline
from numpy.polynomial import chebyshev as cheb
from concurrent.futures import ThreadPoolExecutor
import json

from src.utils.fits_io import read_fits_image, write_fits_image
from src.utils.image_processing import combine_images, find_bad_pixels, estimate_background_2d
from src.core.data_structures import ApertureSet, ApertureLocation, FlatField

logger = logging.getLogger(__name__)


class FlatFieldProcessor:
    """Handles flat field processing and order tracing."""

    def __init__(self):
        """
        Initialize flat field processor.
        """
        self.flat_data: Optional[np.ndarray] = None
        self.flat_mask: Optional[np.ndarray] = None
        self.flat_norm: Optional[np.ndarray] = None
        self.scattered_light: Optional[np.ndarray] = None
        self.order_diagnostics = {}
        self.base_header = None

    def combine_flat_frames(self, flat_filenames: List[str],
                            combine_method: str = 'median',
                            combine_sigma: float = 3.0,
                            mosaic_maxcount: float = 65535) -> Tuple[np.ndarray, np.ndarray]:
        """
        Combine multiple flat frames to create master flat.

        Args:
            flat_filenames: List of flat FITS file paths
            combine_method: Combination method
            combine_sigma: Sigma for clipping
            mosaic_maxcount: Cutoff for bad pixels

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
                if self.base_header is None and header is not None:
                    self.base_header = header.copy()

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
        method = combine_method
        if method == 'median':
            from astropy.stats import sigma_clip
            sigma = combine_sigma
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
        bad_mask = combined > mosaic_maxcount

        self.flat_data = combined
        self.flat_mask = bad_mask.astype(np.int16)

        logger.info(f"Master flat created: shape {combined.shape}, "
                   f"saturated pixels: {np.sum(bad_mask)}")

        return combined, self.flat_mask

    def detect_orders(self, snr_threshold: float = 5.0,
                      min_trace_coverage: float = 0.2,
                      trace_degree: int = 4,
                      output_dir_base: str = '') -> ApertureSet:
        """
        Detect grating orders in flat field with simple tracing.

        For grating spectrographs, orders run horizontally (left to right) across
        the detector. The algorithm collapses the image along the dispersion direction
        (X) to find initial order centroids along Y, then traces each order across
        the image using polynomial fitting.

        Args:
            snr_threshold: Detection threshold (multiples of sigma above baseline)
            gap_fill_snr: SNR threshold for gap fill.
            min_trace_coverage: Minimum fraction of detector width a traced order must cover.
            trace_degree: Polynomial degree for fitting the trace.
            output_dir_base: Output directory for diagnostics.
        Returns:
            ApertureSet with detected orders
        """
        if self.flat_data is None:
            raise RuntimeError("No flat data")

        logger.info("Detecting grating orders (simple tracing)...")

        tracer = GratingOrderTracer(
            data=self.flat_data,
            mask=self.flat_mask,
            snr_threshold=snr_threshold,
            min_trace_coverage=min_trace_coverage,
            trace_degree=trace_degree,
        )
        apertures = tracer.trace_orders()
        self.order_diagnostics = tracer.get_diagnostics()

        return apertures

    def save_flat_field(self, output_path: str,
                        save_plots: bool = True,
                        fig_format: str = 'png'):
        """Save flat field data to FITS file.

        MasterFlat.fits contains only the combined flat image (PRIMARY) and
        the bad-pixel MASK extension so its size is comparable to a single
        input flat frame.
        """
        if self.flat_data is None:
            raise RuntimeError("No flat data to save")

        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        from astropy.io import fits

        # ---- MasterFlat.fits: combined image + mask only ----
        if self.base_header is not None:
            if isinstance(self.base_header, fits.Header):
                header = self.base_header.copy()
            else:
                header = fits.Header()
                for k, v in self.base_header.items():
                    try:
                        header[k] = v
                    except Exception:
                        pass
        else:
            header = fits.Header()
        header['EXTNAME'] = 'MASTER FLAT'
        primary_hdu = fits.PrimaryHDU(data=self.flat_data.astype(np.float32), header=header)
        hdul = fits.HDUList([
            primary_hdu,
            fits.ImageHDU(data=self.flat_mask, name='MASK'),
        ])
        hdul.writeto(str(output_path), overwrite=True)
        logger.info(f"Saved master flat to {output_path}")

        # Save diagnostic plots if enabled
        if save_plots:
            from src.plotting.spectra_plotter import plot_2d_image_to_file
            out_dir = Path(output_path).parent
            plot_2d_image_to_file(self.flat_data, str(out_dir / f'MasterFlat.{fig_format}'), "Master Flat Field")


class GratingOrderTracer:
    """
    Detects and traces orders for grating spectrographs.
    """
    def __init__(self, data: np.ndarray, mask: Optional[np.ndarray] = None, **kwargs):
        self.data = np.asarray(data, dtype=np.float64)
        self.mask = mask if mask is not None else np.zeros_like(data, dtype=bool)
        self.h, self.w = self.data.shape
        self.params = {
            'snr_threshold': 5.0,
            'min_trace_coverage': 0.20,
            'trace_degree': 4,
        }
        self.params.update(kwargs)
        self.diagnostics = {}

    def _get_seed_profile(self) -> Tuple[np.ndarray, dict]:
        """Extracts and analyzes a 1D cross-section to find order seeds."""
        x_center = self.w // 2
        half_band = max(5, self.w // 20)
        x1, x2 = max(0, x_center - half_band), min(self.w, x_center + half_band + 1)
        profile = np.median(self.data[:, x1:x2], axis=1)
        profile = gaussian_filter1d(profile, sigma=1.5)

        # Estimate dynamic order separation
        rough_peaks, _ = find_peaks(gaussian_filter1d(profile, 2.0), distance=5, prominence=np.percentile(profile, 99) * 0.02)
        if len(rough_peaks) >= 5:
            diffs = np.diff(rough_peaks)
            y_mids = 0.5 * (rough_peaks[:-1] + rough_peaks[1:])
            med_diff = float(np.median(diffs))
            # Filter out anomalously large diffs (e.g. missing orders) before fitting
            valid = (diffs > 0.5 * med_diff) & (diffs < 1.5 * med_diff)
            if np.sum(valid) >= 3:
                deg = min(2, np.sum(valid) - 1)
                coefs = np.polyfit(y_mids[valid], diffs[valid], deg)
                med_sep = med_diff
                local_sep_func = lambda y: float(np.clip(np.polyval(coefs, y), med_sep * 0.4, med_sep * 2.5))
            else:
                med_sep = med_diff
                local_sep_func = lambda y: float(med_sep)
        else:
            med_sep = 30.0
            local_sep_func = lambda y: float(med_sep)

        # Background estimation
        from scipy.ndimage import minimum_filter1d
        bg_size = int(2.5 * med_sep)
        envelope = gaussian_filter1d(minimum_filter1d(profile, size=bg_size), sigma=bg_size / 4.0)
        seed_pure = profile - envelope

        # Noise estimation
        valid_prof = seed_pure[np.isfinite(seed_pure)]
        s_med = np.median(np.sort(valid_prof)[:len(valid_prof)//4])
        s_sig = np.median(np.abs(np.diff(valid_prof))) / 0.9539
        noise, baseline = max(s_sig, 0.1), s_med

        diags = {'profile': profile, 'envelope': envelope, 'noise': noise, 'baseline': baseline,
                 'local_sep_func': local_sep_func, 'med_sep': med_sep, 'x1': x1, 'x2': x2}
        return seed_pure, diags

    def _find_seeds(self, seed_profile: np.ndarray, diags: dict) -> np.ndarray:
        """Finds initial peak locations (seeds) in the 1D profile."""
        det_threshold = diags['baseline'] + self.params['snr_threshold'] * diags['noise']
        prominence = max(0.5 * diags['noise'], 0.1 * self.params['snr_threshold'] * diags['noise'])
        
        # Initial peak finding with a small distance to not miss close orders
        peaks, props = find_peaks(seed_profile, height=det_threshold, distance=3, prominence=prominence)
        
        if len(peaks) < 5:
            logger.info(f"Found {len(peaks)} seed peaks (SNR > {self.params['snr_threshold']:.1f})")
            return np.asarray(sorted(peaks), dtype=int)

        # --- Gap detection and seed insertion ---
        sorted_peaks = np.asarray(sorted(peaks), dtype=int)
        diffs = np.diff(sorted_peaks)
        local_sep_func = diags.get('local_sep_func')
        if not local_sep_func:
            med_sep = diags.get('med_sep', np.median(diffs) if len(diffs) > 0 else 30.0)
            local_sep_func = lambda y: float(med_sep)

        new_seeds = list(sorted_peaks)
        for i in range(len(diffs)):
            gap = diffs[i]
            mid_y = (sorted_peaks[i] + sorted_peaks[i+1]) / 2.0
            expected_sep = local_sep_func(mid_y)
            gap_thresh = expected_sep * 1.5 # A gap is > 1.5x the local expected separation
            
            if gap > gap_thresh:
                # Found a gap, try to insert seeds
                n_missing = int(round(gap / expected_sep)) - 1
                if n_missing > 0 and n_missing <= 5: # Limit insertions
                    logger.info(f"Gap detected between {sorted_peaks[i]} and {sorted_peaks[i+1]} (gap={gap:.1f}, expected_sep={expected_sep:.1f}). Trying to insert {n_missing} seed(s).")
                    for j in range(1, n_missing + 1):
                        # Predict seed position
                        predicted_y = int(round(sorted_peaks[i] + j * gap / (n_missing + 1)))
                        
                        # Search for a peak near the predicted position, even if it's below the main SNR threshold
                        search_win = max(2, int(expected_sep * 0.2))
                        y1 = max(0, predicted_y - search_win)
                        y2 = min(len(seed_profile), predicted_y + search_win + 1)
                        
                        if y2 > y1:
                            local_profile = seed_profile[y1:y2]
                            # Use a lower threshold for gap filling
                            gap_fill_thresh = diags['baseline'] + 1.0 * diags['noise'] # Lower SNR for gaps
                            local_peaks, _ = find_peaks(local_profile, height=gap_fill_thresh, distance=3)
                            if len(local_peaks) > 0:
                                best_local_peak_idx = np.argmax(local_profile[local_peaks])
                                inserted_y = y1 + local_peaks[best_local_peak_idx]
                                new_seeds.append(inserted_y)
                                logger.info(f"  -> Inserted seed at y={inserted_y} (predicted at {predicted_y})")
            
        final_seeds = np.asarray(sorted(list(set(new_seeds))), dtype=int)
        logger.info(f"Found {len(peaks)} initial peaks, {len(final_seeds)} seeds after gap filling.")
        return final_seeds

    def _trace_one_order(self, seed_y: float, diags: dict) -> Optional[Tuple[np.ndarray, np.ndarray]]:
        """Traces a single order starting from a seed."""
        xs, ys = [], []

        def refine_peak(x, y_guess, search_half):
            y0 = int(round(y_guess))
            y1, y2 = max(0, y0 - search_half), min(self.h, y0 + search_half + 1)
            if y2 - y1 < 3: return None
            col = self.data[y1:y2, x]
            if not np.any(np.isfinite(col)): return None
            
            loc = int(np.nanargmax(col))
            
            # A true peak must be a local maximum, not just the edge of the window
            if loc == 0 or loc == len(col) - 1:
                return None, 0.0
            
            bg_local = np.nanmin(col)
            flux_bg_sub = col[loc] - bg_local
            # To prevent losing the track over absorption lines, we lower the hard cutoff
            # to 3.0 sigma. If it's between 3.0 and snr_threshold, we will use it for
            # steering but NOT for final trace fitting.
            if flux_bg_sub < 3.0 * diags['noise']:
                return None, 0.0
            
            # Sub-pixel centroiding
            hw = 2
            c1, c2 = max(0, loc - hw), min(len(col), loc + hw + 1)
            sub_col = col[c1:c2]
            valid = np.isfinite(sub_col)
            if np.sum(valid) >= 3:
                yy = np.arange(c1, c2)[valid]
                bg = np.min(sub_col[valid])
                weights = np.maximum(0, sub_col[valid] - bg)
                sum_w = np.sum(weights)
                if sum_w > 0:
                    peak_y_subpix = float(y1 + np.sum(yy * weights) / sum_w)
                    peak_snr = flux_bg_sub / diags['noise']
                    return peak_y_subpix, peak_snr
            peak_snr = flux_bg_sub / diags['noise']
            return float(y1 + loc), peak_snr # Fallback to integer pixel if sub-pixel fails

        # Start from the center, calculate dynamic search window
        initial_sep = diags['local_sep_func'](seed_y)
        initial_search_half = max(3, int(initial_sep * 0.25))
        
        res = refine_peak(self.w // 2, seed_y, initial_search_half)
        if res[0] is None: return None
        start_y, start_snr = res
        
        xs.append(float(self.w // 2))
        ys.append(start_y)

        # Trace right and left
        for direction in [1, -1]:
            hist_x = [float(self.w // 2)]
            hist_y = [start_y]
            miss_count = 0 # Counter for consecutive columns with signal below threshold or no peak
            
            x_range = range(self.w // 2 + direction, self.w if direction == 1 else -1, direction)
            for x in x_range:
                # Extrapolate y_guess using a polynomial fit to the recent trace history
                if len(hist_x) >= 15:
                    # Use quadratic fit for better curve handling
                    n_pts, deg = 15, 2
                elif len(hist_x) >= 5:
                    # Use linear fit for stability
                    n_pts, deg = 5, 1
                else:
                    n_pts, deg = 0, 0

                if n_pts > 0:
                    fit_x = np.array(hist_x[-n_pts:])
                    fit_y = np.array(hist_y[-n_pts:])
                    # Shift x to be around 0 for better numerical stability
                    x0 = fit_x[-1]
                    try:
                        coeffs = np.polyfit(fit_x - x0, fit_y, deg)
                        y_guess = np.polyval(coeffs, x - x0)
                    except (np.linalg.LinAlgError, ValueError):
                        y_guess = hist_y[-1] # Fallback
                elif len(hist_x) >= 2:
                    slope = (hist_y[-1] - hist_y[-2]) / (hist_x[-1] - hist_x[-2])
                    y_guess = hist_y[-1] + slope * (x - hist_x[-1])
                else:
                    y_guess = hist_y[-1]
                
                # Dynamically adjust search window based on the smooth separation curve
                local_sep = diags['local_sep_func'](y_guess)
                search_half = max(2, int(local_sep * 0.25))
                
                res = refine_peak(x, y_guess, search_half)
                if res[0] is not None and abs(res[0] - y_guess) <= search_half:
                    y_new, peak_snr = res
                    # Update steering history even if SNR is low (e.g. inside absorption line)
                    # so we don't drift off the track.
                    hist_x.append(float(x))
                    hist_y.append(y_new)
                    
                    if peak_snr >= self.params['snr_threshold']:
                        # Strong peak. Save to final trace points.
                        xs.append(float(x))
                        ys.append(y_new)
                        miss_count = 0
                    else:
                        # Weak peak (absorption line). Steer, but don't save to final trace.
                        miss_count += 1
                else:
                    # Completely lost the peak (e.g., edge of detector)
                    miss_count += 1 
                
                # Allow a coasting distance of 50 pixels to jump over wide absorption lines.
                max_miss = 50
                if miss_count > max_miss:
                    break

        if len(xs) < self.params['min_trace_coverage'] * self.w:
            return None

        xs, ys = np.array(xs), np.array(ys)
        order = np.argsort(xs)
        return xs[order], ys[order]

    def trace_orders(self) -> ApertureSet:
        """Main method to perform order tracing."""
        seed_profile, diags = self._get_seed_profile()
        self.diagnostics.update(diags)
        self.diagnostics['seed_profile'] = seed_profile
        seeds = self._find_seeds(seed_profile, diags)
        self.diagnostics['seeds'] = seeds

        aperture_set = ApertureSet()
        order_id_counter = 1
        for seed in seeds:
            trace_result = self._trace_one_order(float(seed), diags)
            if trace_result is None:
                continue

            x_pts, y_pts = trace_result
            if len(x_pts) < 10: continue

            # Fit with Chebyshev polynomial
            deg = min(self.params['trace_degree'], len(x_pts) - 2)
            if deg < 1: continue
            domain = (0.0, float(self.w - 1))
            x_mapped = 2.0 * (x_pts - domain[0]) / (domain[1] - domain[0]) - 1.0
            center_coef = cheb.chebfit(x_mapped, y_pts, deg)

            # Estimate width
            local_sep = diags['local_sep_func'](float(seed))
            width = 0.4 * local_sep # Simplified width estimation

            ap = ApertureLocation(
                aperture=order_id_counter, order=order_id_counter,
                center_coef=center_coef, width=width,
                is_chebyshev=True, domain=domain
            )
            aperture_set.add_aperture(ap)
            order_id_counter += 1

        logger.info(f"Successfully traced {aperture_set.norders} orders.")
        return aperture_set

    def get_diagnostics(self) -> dict:
        return self.diagnostics


def assign_order_indices(apertures: ApertureSet, gap_fill_factor: float = 1.35) -> dict:
    """
    Assign integer order indices *m* to detected apertures, detecting gaps.

    Walks ordered apertures from bottom (small Y) to top (large Y) and increments *m*
    by 1 for each normal step. When the gap between two consecutive centres
    exceeds *gap_fill_factor × local_expected_sep*, it skips the appropriate number
    of *m* values, reserving slots for missing orders.

    The expected separation is modelled as a smooth polynomial of sequential
    index to adapt to monotonically changing spacing (e.g., in echelle spectra).

    Args:
        apertures: The set of detected apertures.
        gap_fill_factor: Factor to identify a gap. A separation > factor * expected
                         is considered a gap.

    Returns:
        A dictionary mapping original aperture_id to the new gap-aware index `m`.
    """
    order_ids = apertures.get_orders()
    if len(order_ids) < 2:
        return {oid: i for i, oid in enumerate(order_ids)}

    x_ref = apertures.get_aperture(order_ids[0]).domain[1] / 2.0
    centres = {oid: float(apertures.get_aperture(oid).get_position(x_ref)) for oid in order_ids}

    sorted_ids = sorted(order_ids, key=lambda oid: centres[oid])
    diffs = np.array([centres[sorted_ids[i + 1]] - centres[sorted_ids[i]] for i in range(len(sorted_ids) - 1)])
    med_sep = float(np.median(diffs)) if len(diffs) > 0 else 30.0

    # Build a smooth, robust model for the expected local separation.
    n_diffs = len(diffs)
    if n_diffs >= 6:
        idx_arr = np.arange(n_diffs, dtype=float)
        fit_deg = min(3, n_diffs - 2)
        good = np.ones(n_diffs, dtype=bool)

        # Iteratively clip outliers to get a robust fit for separation trend
        for _ in range(4):
            if np.sum(good) < fit_deg + 2:
                break
            coef_tmp = np.polyfit(idx_arr[good], diffs[good], fit_deg)
            model_tmp = np.polyval(coef_tmp, idx_arr)
            resid = diffs - model_tmp
            std_r = np.std(resid[good])
            if std_r < 1e-6:
                break
            new_good = np.abs(resid) < 2.5 * std_r
            if np.array_equal(new_good, good):
                break
            good = new_good

        sep_coef = np.polyfit(idx_arr[good], diffs[good], fit_deg)
        expected_seps = np.polyval(sep_coef, idx_arr)
        inlier_min = np.min(diffs[good]) if np.any(good) else med_sep * 0.5
        inlier_max = np.max(diffs[good]) if np.any(good) else med_sep * 2.0
        expected_seps = np.clip(expected_seps, inlier_min * 0.5, inlier_max * 2.0)
    else:
        expected_seps = np.full(n_diffs, med_sep)

    mapping = {}
    m = 0
    mapping[sorted_ids[0]] = m
    for i in range(len(diffs)):
        gap = diffs[i]
        local_sep = float(expected_seps[i])
        ratio = gap / max(local_sep, 1.0)
        n_skip = int(round(ratio)) if ratio >= gap_fill_factor else 1
        m += n_skip
        mapping[sorted_ids[i + 1]] = m

    logger.info(f"Order index assignment: {len(mapping)} apertures, m_max={m}, "
                f"med_sep={med_sep:.1f}px, gap_factor={gap_fill_factor:.2f}")
    return mapping


def renumber_apertures_with_gaps(apertures: ApertureSet, gap_fill_factor: float = 1.35) -> ApertureSet:
    """
    Renumber aperture IDs to reflect physical order indices, including gaps.

    After initial detection, apertures might have sequential IDs (1, 2, 3, ...).
    This function re-assigns IDs based on the gap-aware index `m` from
    `assign_order_indices`, so that missing orders create gaps in the IDs
    (e.g., 1, 2, 4, 5 if order 3 is missing). The new ID is `m + 1`.

    Args:
        apertures: The original ApertureSet.
        gap_fill_factor: Factor for gap detection, passed to `assign_order_indices`.

    Returns:
        A new ApertureSet with renumbered apertures.
    """
    if apertures.norders == 0:
        return apertures

    order_map = assign_order_indices(apertures, gap_fill_factor)
    new_set = ApertureSet()
    for old_id, m in order_map.items():
        ap = apertures.get_aperture(old_id)
        new_id = m + 1
        # Create a new ApertureLocation with the updated ID
        new_ap = ApertureLocation(
            aperture=new_id,
            order=new_id,
            center_coef=ap.center_coef.copy() if ap.center_coef is not None else None,
            lower_coef=ap.lower_coef.copy() if ap.lower_coef is not None else None,
            upper_coef=ap.upper_coef.copy() if ap.upper_coef is not None else None,
            width=ap.width,
            is_chebyshev=ap.is_chebyshev,
            domain=ap.domain,
            is_interpolated=ap.is_interpolated,
        )
        # Preserve custom attributes if they exist
        if getattr(ap, 'center_arr', None) is not None:
            new_ap.center_arr = ap.center_arr.copy()
        if getattr(ap, 'w_up_cheb_coef', None) is not None:
            new_ap.w_up_cheb_coef = ap.w_up_cheb_coef.copy()
        if getattr(ap, 'w_low_cheb_coef', None) is not None:
            new_ap.w_low_cheb_coef = ap.w_low_cheb_coef.copy()

        new_set.add_aperture(new_ap)

    if set(apertures.get_orders()) != set(new_set.get_orders()):
        logger.info(f"Renumbered apertures with gap-aware IDs: {apertures.norders} -> {new_set.norders} orders.")
    return new_set


def _refine_predicted_trace(
    flat_data: np.ndarray,
    pred_center: np.ndarray,
    local_sep: float,
    trace_degree: int = 3,
    step: int = 20
) -> Optional[np.ndarray]:
    """
    Refine a predicted trace by finding peaks in the actual flat image.

    Args:
        flat_data: The 2D flat field image.
        pred_center: The predicted center trace array.
        local_sep: The local order separation in pixels.
        trace_degree: The polynomial degree for fitting the refined trace.
        step: The column step for sampling peaks.

    Returns:
        A refined center trace array, or None if refinement fails.
    """
    h, w = flat_data.shape
    search_half = max(6, int(0.4 * local_sep))
    xs, ys = [], []

    for x in range(0, w, step):
        y_guess = pred_center[x]
        y0 = int(round(y_guess))
        y1 = max(0, y0 - search_half)
        y2 = min(h, y0 + search_half + 1)
        if y2 - y1 < 5:
            continue

        col = flat_data[y1:y2, x].astype(float)
        if not np.any(np.isfinite(col)):
            continue

        loc = int(np.nanargmax(col))
        peak_val = col[loc]
        if not np.isfinite(peak_val) or peak_val <= 0:
            continue

        # Centroid refinement
        hw = 3
        g1, g2 = max(0, loc - hw), min(col.size, loc + hw + 1)
        sub = col[g1:g2]
        valid_sub = np.isfinite(sub)
        if np.sum(valid_sub) < 3:
            refined_y = float(y1 + loc)
        else:
            yy = np.arange(g1, g2)[valid_sub]
            weights = sub[valid_sub]
            bg = np.min(weights)
            weights = np.maximum(0, weights - bg)
            sum_w = np.sum(weights)
            if sum_w > 0:
                refined_y = float(y1 + np.sum(yy * weights) / sum_w)
            else:
                refined_y = float(y1 + loc)

        if abs(refined_y - y_guess) <= search_half:
            xs.append(float(x))
            ys.append(refined_y)

    if len(xs) < max(5, int(0.15 * (w / step))):
        return None

    xs_arr, ys_arr = np.asarray(xs), np.asarray(ys)
    deg = min(trace_degree, max(1, len(xs_arr) - 2))

    # Robust polynomial fit with sigma-clipping
    for _ in range(3):
        if len(xs_arr) < deg + 2:
            return None
        coef = np.polyfit(xs_arr, ys_arr, deg)
        resid = ys_arr - np.polyval(coef, xs_arr)
        sigma = max(0.5, float(np.std(resid)))
        keep = np.abs(resid) < 2.5 * sigma
        if np.all(keep):
            break
        xs_arr, ys_arr = xs_arr[keep], ys_arr[keep]

    final_coef = np.polyfit(xs_arr, ys_arr, deg)
    return np.polyval(final_coef, np.arange(w, dtype=float))


def fill_missing_orders_by_interpolation(
    apertures: ApertureSet,
    image_width: int,
    image_height: int,
    trace_degree: int = 3,
    n_extend_below: int = 0,
    n_extend_above: int = 0,
    flat_data: Optional[np.ndarray] = None,
    gap_fill_factor: float = 1.35
) -> ApertureSet:
    """
    Fill interior gaps and extend edges of an ApertureSet.

    Interior gaps are filled by interpolating traces between bounding real orders.
    Edge orders are extrapolated based on the spacing of the nearest real orders.
    If `flat_data` is provided, interpolated traces are refined against the image.

    Args:
        apertures: The ApertureSet, assumed to be renumbered with gaps.
        image_width: The width of the detector image.
        image_height: The height of the detector image.
        trace_degree: Polynomial degree for fitting new traces.
        n_extend_below: Number of orders to extrapolate below the first real order.
        n_extend_above: Number of orders to extrapolate above the last real order.
        flat_data: Optional 2D flat image for refining interpolated traces.
        gap_fill_factor: Factor used for gap detection.

    Returns:
        A new ApertureSet containing original, interpolated, and extrapolated orders.
    """
    if apertures.norders < 2:
        return apertures.copy()

    order_map = assign_order_indices(apertures, gap_fill_factor)
    m_to_oid = {m: oid for oid, m in order_map.items()}
    sorted_m = sorted(m_to_oid.keys())

    filled_set = ApertureSet()
    x_all = np.arange(image_width, dtype=float)
    x_mapped = 2.0 * x_all / (image_width - 1.0) - 1.0

    # --- 1. Fill interior gaps by interpolation ---
    for i in range(len(sorted_m) - 1):
        m_lo, m_hi = sorted_m[i], sorted_m[i + 1]
        oid_lo, oid_hi = m_to_oid[m_lo], m_to_oid[m_hi]

        # Add the lower-bound real order
        filled_set.add_aperture(apertures.get_aperture(oid_lo))

        n_gap = m_hi - m_lo - 1
        if n_gap <= 0:
            continue

        ap_lo = apertures.get_aperture(oid_lo)
        ap_hi = apertures.get_aperture(oid_hi)
        traces = {
            "cen": (ap_lo.get_position(x_all), ap_hi.get_position(x_all)),
            "low": (ap_lo.get_lower(x_all), ap_hi.get_lower(x_all)),
            "upp": (ap_lo.get_upper(x_all), ap_hi.get_upper(x_all)),
        }

        for k in range(1, n_gap + 1):
            frac = k / (n_gap + 1.0)
            pred_traces = {key: val[0] + frac * (val[1] - val[0]) for key, val in traces.items()}

            yc_mid = pred_traces["cen"][image_width // 2]
            if not (0 <= yc_mid < image_height):
                continue

            final_center = pred_traces["cen"]
            if flat_data is not None:
                local_sep = float(np.median(traces["cen"][1] - traces["cen"][0])) / (n_gap + 1)
                refined = _refine_predicted_trace(flat_data, pred_traces["cen"], local_sep, trace_degree)
                if refined is not None:
                    final_center = refined
                    logger.info(f"  Refined interpolated order m={m_lo + k}")

            # Reconstruct boundaries relative to the (potentially refined) center
            half_width_low = pred_traces["cen"] - pred_traces["low"]
            half_width_upp = pred_traces["upp"] - pred_traces["cen"]
            final_lower = final_center - half_width_low
            final_upper = final_center + half_width_upp

            new_id = m_lo + k + 1
            ap = ApertureLocation(
                aperture=new_id, order=new_id,
                center_coef=cheb.chebfit(x_mapped, final_center, trace_degree),
                lower_coef=cheb.chebfit(x_mapped, final_lower, trace_degree),
                upper_coef=cheb.chebfit(x_mapped, final_upper, trace_degree),
                width=float(np.median(final_upper - final_lower)),
                is_chebyshev=True,
                domain=(0.0, float(image_width - 1.0)),
                is_interpolated=True,
            )
            filled_set.add_aperture(ap)
            logger.info(f"  Filled gap order m={m_lo + k} (ID={new_id}) at y~{yc_mid:.1f}")

    # Add the last real order
    filled_set.add_aperture(apertures.get_aperture(m_to_oid[sorted_m[-1]]))

    # --- 2. Extrapolate orders at the edges ---
    if n_extend_below > 0:
        _extrapolate_edge_orders(filled_set, 'below', n_extend_below, image_width, image_height, trace_degree)
    if n_extend_above > 0:
        _extrapolate_edge_orders(filled_set, 'above', n_extend_above, image_width, image_height, trace_degree)

    logger.info(f"Interpolation/extrapolation complete: {apertures.norders} -> {filled_set.norders} orders.")
    return filled_set


def _extrapolate_edge_orders(
    apertures: ApertureSet,
    side: str,
    n_orders: int,
    image_width: int,
    image_height: int,
    trace_degree: int,
    n_ref: int = 8
):
    """Helper to extrapolate and add virtual orders to an ApertureSet in-place."""
    real_oids = sorted(apertures.get_orders())
    if len(real_oids) < 2:
        return

    if side == 'below':
        ref_oids = real_oids[:min(n_ref, len(real_oids))]
        direction = -1
        anchor_ap = apertures.get_aperture(ref_oids[0])
        start_idx = -1
    else: # 'above'
        ref_oids = real_oids[max(0, len(real_oids) - n_ref):]
        direction = 1
        anchor_ap = apertures.get_aperture(ref_oids[-1])
        start_idx = len(ref_oids) - 2

    x_all = np.arange(image_width, dtype=float)
    ref_centers = np.array([apertures.get_aperture(oid).get_position(x_all) for oid in ref_oids])
    spacings = np.diff(ref_centers, axis=0)

    # Model spacing trend with a simple linear fit
    if spacings.shape[0] >= 2:
        j_sp = np.arange(spacings.shape[0], dtype=float)
        b_sp = np.polyfit(j_sp, spacings, 1)[0] # slope only
    else:
        b_sp = np.zeros(image_width)

    last_center = anchor_ap.get_position(x_all)
    last_lower = anchor_ap.get_lower(x_all)
    last_upper = anchor_ap.get_upper(x_all)
    last_sep = spacings[start_idx]

    for k in range(1, n_orders + 1):
        # Extrapolate spacing, center, and width
        current_sep = last_sep + k * b_sp * direction
        new_center = last_center + current_sep * direction
        width = last_upper - last_lower
        new_lower, new_upper = new_center - width / 2, new_center + width / 2

        if (side == 'below' and np.max(new_upper) < 0) or \
           (side == 'above' and np.min(new_lower) >= image_height):
            break # Stop if order is completely off-detector

        x_mapped = 2.0 * x_all / (image_width - 1.0) - 1.0

        new_id = anchor_ap.aperture + k * direction
        ap = ApertureLocation(
            aperture=new_id, order=new_id,
            center_coef=cheb.chebfit(x_mapped, new_center, trace_degree),
            lower_coef=cheb.chebfit(x_mapped, new_lower, trace_degree),
            upper_coef=cheb.chebfit(x_mapped, new_upper, trace_degree),
            width=float(np.median(width)),
            is_chebyshev=True,
            domain=(0.0, float(image_width - 1.0)),
            is_interpolated=True,
        )
        apertures.add_aperture(ap)
        logger.info(f"  Extrapolated virtual order (ID={new_id}) at y~{new_center[image_width//2]:.1f}")


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


def process_order_tracing_stage(flat_filenames: List[str],
                      output_dir_base: str = './midpath',
                      # Flat combination parameters
                      combine_method: str = 'median',
                      combine_sigma: float = 3.0,
                      mosaic_maxcount: float = 65535,
                      # Order detection parameters
                      snr_threshold: float = 5.0,
                      gap_fill_factor: float = 1.6,
                      min_trace_coverage: float = 0.20,
                      trace_degree: int = 4,
                      # Profile boundary parameters
                      width_cheb_degree: int = 3,
                      aperture_boundary_snr: float = 3.0,
                      # Gap fill / extend parameters
                      n_extend_below: int = 0,
                      n_extend_above: int = 0,
                      gap_fill_factor_interp: float = 1.35,
                      # Diagnostic parameters
                      save_plots: bool = True,
                      fig_format: str = 'png') -> Tuple[FlatField, ApertureSet]:
    """_
    Execute flat fielding stage.

    Args:
        flat_filenames: List of flat frame paths
        output_dir_base: Base output directory

    Returns:
        Tuple of (FlatField, ApertureSet)
    """
    processor = FlatFieldProcessor()

    # Combine flat frames
    flat_data, flat_mask = processor.combine_flat_frames(
        flat_filenames,
        combine_method=combine_method,
        combine_sigma=combine_sigma,
        mosaic_maxcount=mosaic_maxcount)

    # Save master flat field as FITS and PNG (combined image + mask only)
    logger.info("Saving master flat field...")
    base_output_path = output_dir_base
    flat_file = Path(base_output_path) / 'step2_trace' / 'MasterFlat.fits'
    flat_file.parent.mkdir(parents=True, exist_ok=True)
    processor.save_flat_field(str(flat_file), save_plots=save_plots, fig_format=fig_format)
    logger.info(f"Saved master flat field to {flat_file}")

    # Detect orders — uses raw combined flat, no normalization needed
    logger.info("Detecting grating orders...")
    apertures = processor.detect_orders( # Pass relevant params from config
        snr_threshold=snr_threshold,
        min_trace_coverage=min_trace_coverage,
        trace_degree=trace_degree,
        output_dir_base=output_dir_base)
    logger.info(f"Order detection completed: {apertures.norders} orders found")

    # Renumber apertures with gap-aware physical order indices.
    # If spatial analysis detects missing orders, aperture IDs skip the
    # missing positions (e.g. 1,...,49,51,...,100) so that the downstream
    # 2D surface fitting uses correct order index m for every aperture.
    apertures = renumber_apertures_with_gaps(apertures, gap_fill_factor)

    # Pre-compute profile-root-based boundaries for order mask + diagnostic plots.
    # _find_profile_roots_1d finds where the signal drops to aperture_boundary_snr × σ
    # above background, which is always strictly inside the inter-order midpoint.  This
    # guarantees adjacent order masks never overlap even for densely packed echelle orders.
    logger.info("Pre-computing profile-root boundaries for order mask and diagnostics...")
    height, width = flat_data.shape
    x_coords_all = np.arange(width)
    center_col = width // 2

    def _estimate_noise_floor(cross_section_1d, bg_envelope_1d, centers_y):
        """Estimate background noise via iterative MAD sigma-clipping.""" 
        residual = cross_section_1d - bg_envelope_1d
        valid = np.isfinite(residual)
        res = residual[valid]

        if res.size < 5:
            return float(np.nanstd(residual))

        sorted_res = np.sort(res)
        bg_pixels = sorted_res[:max(5, len(sorted_res) // 4)]
        med = float(np.nanmedian(bg_pixels))
        mad = float(np.nanmedian(np.abs(bg_pixels - med)))
        sigma_est = max(1.4826 * mad, 1e-6)
        keep = res < (med + 3.0 * sigma_est)

        for _ in range(5):
            if np.sum(keep) < 5: break
            sub = res[keep]
            med = float(np.nanmedian(sub))
            mad = float(np.nanmedian(np.abs(sub - med)))
            sigma_est = max(1.4826 * mad, 1e-6)
            new_keep = res < (med + 3.0 * sigma_est)
            if np.array_equal(new_keep, keep): break
            keep = new_keep

        if np.sum(keep) >= 5:
            sub = res[keep]
            med = float(np.nanmedian(sub))
            mad = float(np.nanmedian(np.abs(sub - med)))
            sigma_bg = max(1.4826 * mad, 1e-6)
        else:
            sigma_bg = max(sigma_est, 1e-6)

        return sigma_bg * aperture_boundary_snr

    ordered_apertures = sorted(apertures.apertures.items(), key=lambda x: x[0])
    valid_orders = []
    for order_idx, (aperture_id, aperture) in enumerate(ordered_apertures, start=1):
        center_pos = aperture.get_position(x_coords_all)
        valid_center = np.isfinite(center_pos) & (center_pos >= 0) & (center_pos < height)
        if np.sum(valid_center) < 10:
            continue
        valid_orders.append({
            'order_idx': order_idx,
            'aperture_id': aperture_id,
            'aperture': aperture,
            'valid_center': valid_center,
        })

    y_lo_all = np.full((len(valid_orders), width), np.nan, dtype=float)
    y_hi_all = np.full((len(valid_orders), width), np.nan, dtype=float)
    if len(valid_orders) >= 1:
        def process_column(xi):
            x = x_coords_all[xi]
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
                return xi, None
            bg_env_col = _fit_background_envelope_1d(cross_section_col, centers)
            noise_floor_col = _estimate_noise_floor(cross_section_col, bg_env_col, centers)
            results = []
            for local_i, (row_idx, cy_raw) in enumerate(active_rows):
                cy_y = float(np.clip(cy_raw, 0, height - 1))
                cy_idx = int(round(cy_y))
                if local_i > 0:
                    search_lo = 0.5 * (centers[local_i] + centers[local_i - 1])
                else:
                    gap = (centers[1] - centers[0]) if len(centers) > 1 else 20.0
                    search_lo = max(0.0, centers[local_i] - 0.5 * gap)
                if local_i < len(centers) - 1:
                    search_hi = 0.5 * (centers[local_i] + centers[local_i + 1])
                else:
                    gap = (centers[-1] - centers[-2]) if len(centers) > 1 else 20.0
                    search_hi = min(float(height - 1), centers[local_i] + 0.5 * gap)
                root_lo, root_hi, _bg_at_peak, _threshold = _find_profile_roots_1d(
                    cross_section_col, cy_idx, search_lo, search_hi,
                    threshold_frac=0.0, bg_envelope=bg_env_col, noise_floor=noise_floor_col)
                results.append((row_idx, root_lo, root_hi))
            return xi, results

        with ThreadPoolExecutor() as executor:
            for xi, results in executor.map(process_column, range(width)):
                if results is None:
                    continue
                for row_idx, root_lo, root_hi in results:
                    y_lo_all[row_idx, xi] = float(np.clip(root_lo, 0, height - 1))
                    y_hi_all[row_idx, xi] = float(np.clip(root_hi, 0, height - 1))
                    
    logger.info("Filtering abnormal boundary widths via cross-order interpolation...")
    centers_all = np.full((len(valid_orders), width), np.nan)
    for row_idx, item in enumerate(valid_orders):
        centers_all[row_idx, :] = item['aperture'].get_position(x_coords_all)

    w_up_all = y_hi_all - centers_all
    w_low_all = centers_all - y_lo_all

    from scipy.ndimage import median_filter
    for xi in range(width):
        for w_arr in (w_up_all, w_low_all):
            col_w = w_arr[:, xi]
            valid = np.isfinite(col_w) & (col_w > 0)
            if np.sum(valid) >= 3:
                idx_valid = np.where(valid)[0]
                val_valid = col_w[valid]
                
                if len(val_valid) >= 5:
                    trend = median_filter(val_valid, size=5, mode='nearest')
                else:
                    trend = np.full_like(val_valid, np.median(val_valid))
                
                # 识别异常偏宽的边界：宽度 > 趋势宽度的 1.4 倍，或者比趋势宽出 4 个像素
                outlier_mask = (val_valid > 1.4 * trend) | (val_valid > trend + 4.0)
                
                if np.any(outlier_mask):
                    good_mask = ~outlier_mask
                    if np.sum(good_mask) >= 2:
                        val_valid[outlier_mask] = np.interp(idx_valid[outlier_mask], idx_valid[good_mask], val_valid[good_mask])
                    elif np.sum(good_mask) == 1:
                        val_valid[outlier_mask] = val_valid[good_mask][0]
                    col_w[idx_valid] = val_valid
                    
    y_hi_all = centers_all + w_up_all
    y_lo_all = centers_all - w_low_all

    # Resolve overlaps in raw boundaries before fitting to guide the polynomials
    for i in range(len(valid_orders) - 1):
        for xi in range(width):
            if y_hi_all[i, xi] >= y_lo_all[i+1, xi] - 1.0:
                c1 = int(round(centers_all[i, xi]))
                c2 = int(round(centers_all[i+1, xi]))
                c1 = max(0, min(c1, height - 1))
                c2 = max(0, min(c2, height - 1))
                if c2 > c1 + 2:
                    valley = c1 + np.argmin(flat_data[c1:c2, xi])
                else:
                    valley = int(round(0.5 * (c1 + c2)))
                y_hi_all[i, xi] = valley - 0.5
                y_lo_all[i+1, xi] = valley + 0.5

    x_all_f = x_coords_all.astype(float)

    logger.info(f"Hybrid smoothing: Chebyshev center (deg {trace_degree}) + Chebyshev width (deg {width_cheb_degree})")

    _trace_coefs = {}

    for row_idx, item in enumerate(valid_orders):
        aperture = item['aperture']
        valid_center = item['valid_center']
        center_raw = aperture.get_position(x_all_f)
        y_lo_raw = y_lo_all[row_idx, :]
        y_hi_raw = y_hi_all[row_idx, :]
        valid_bnd = valid_center & np.isfinite(y_lo_raw) & np.isfinite(y_hi_raw)
        n_valid = int(np.sum(valid_bnd))
        if n_valid < 10:
            continue

        center_smooth = center_raw
        w_up_raw  = y_hi_raw[valid_bnd] - center_smooth[valid_bnd]
        w_low_raw = center_smooth[valid_bnd] - y_lo_raw[valid_bnd]

        x_bnd = x_all_f[valid_bnd]
        span = x_bnd[-1] - x_bnd[0] if len(x_bnd) > 0 else 0
        deg = min(width_cheb_degree, n_valid - 1) if span >= 0.6 * width else (min(1, width_cheb_degree, n_valid - 1) if span >= 0.3 * width else 0)

        coef_w_up  = cheb.chebfit(x_bnd, w_up_raw,  deg)
        coef_w_low = cheb.chebfit(x_bnd, w_low_raw, deg)

        w_up_smooth  = np.maximum(cheb.chebval(x_all_f, coef_w_up), 0.0)
        w_low_smooth = np.maximum(cheb.chebval(x_all_f, coef_w_low), 0.0)

        hi_final = np.clip(center_smooth + w_up_smooth, 0, height - 1)
        lo_final = np.clip(center_smooth - w_low_smooth, 0, height - 1)

        y_lo_all[row_idx, :] = np.where(valid_center, lo_final, np.nan)
        y_hi_all[row_idx, :] = np.where(valid_center, hi_final, np.nan)

        _trace_coefs[item['aperture_id']] = {
            'center_arr': center_smooth.tolist(),
            'w_up_cheb':  coef_w_up.tolist(),
            'w_low_cheb': coef_w_low.tolist(),
        }

    logger.info("Hybrid smoothing complete (Chebyshev center + Chebyshev width → traces)")

    for item in valid_orders:
        aid = item['aperture_id']
        if aid in _trace_coefs:
            bc = _trace_coefs[aid]
            ap = item['aperture']
            ap.center_arr = np.asarray(bc['center_arr'], dtype=float)
            ap.w_up_cheb_coef = np.asarray(bc['w_up_cheb'], dtype=float)
            ap.w_low_cheb_coef = np.asarray(bc['w_low_cheb'], dtype=float)

    filled_apertures = fill_missing_orders_by_interpolation(
        apertures, width, height, trace_degree=trace_degree,
        n_extend_below=n_extend_below, n_extend_above=n_extend_above,
        flat_data=flat_data, gap_fill_factor=gap_fill_factor_interp)

    if filled_apertures and filled_apertures.norders > apertures.norders:
        apertures = filled_apertures
        logger.info(f"Gap-fill added {filled_apertures.norders - apertures.norders} predicted order(s)")

    coefs_path = Path(base_output_path) / 'step2_trace' / 'Orders_trace_coefs.json'
    with open(coefs_path, 'w') as f:
        json.dump({'orders': _trace_coefs}, f, indent=2)
    logger.info(f"Saved trace coefs ({len(_trace_coefs)} orders) to {coefs_path.name}")

    all_ap_ids = sorted(list(apertures.apertures.keys()))
    n_ap = len(all_ap_ids)
    lo_traces = np.zeros((n_ap, width))
    hi_traces = np.zeros((n_ap, width))
    centers = np.zeros((n_ap, width))

    for i, ap_id in enumerate(all_ap_ids):
        ap = apertures.get_aperture(ap_id)
        centers[i] = ap.get_position(x_all_f)
        lo_traces[i] = ap.get_lower(x_all_f)
        hi_traces[i] = ap.get_upper(x_all_f)

    # Enforce valley gap for the final mask and plot
    for i in range(n_ap - 1):
        for xi in range(width):
            if hi_traces[i, xi] >= lo_traces[i+1, xi] - 1.0:
                c1 = int(round(centers[i, xi]))
                c2 = int(round(centers[i+1, xi]))
                c1 = max(0, min(c1, height - 1))
                c2 = max(0, min(c2, height - 1))
                if c2 > c1 + 2:
                    valley = c1 + np.argmin(flat_data[c1:c2, xi])
                else:
                    valley = int(round(0.5 * (c1 + c2)))
                hi_traces[i, xi] = valley - 0.5
                lo_traces[i+1, xi] = valley + 0.5

    # Use uint8 for boolean mask (1 = order, 0 = background)
    order_mask = np.zeros((height, width), dtype=np.uint8)
    for i, ap_id in enumerate(all_ap_ids):
        for xi in range(width):
            if np.isfinite(lo_traces[i, xi]) and np.isfinite(hi_traces[i, xi]):
                r0 = max(0, int(np.ceil(lo_traces[i, xi])))
                r1 = min(height, int(np.floor(hi_traces[i, xi])) + 1)
                if r1 > r0:
                    order_mask[r0:r1, xi] = 1

    mask_path = Path(base_output_path) / 'step2_trace' / 'Orders_mask.fits'
    write_fits_image(str(mask_path), order_mask, dtype='uint8')
    logger.info(f"Saved order boolean mask to {mask_path}")

    if save_plots and apertures.norders > 0:
        logger.info("Generating diagnostic plots...")
        out_dir = Path(base_output_path) / 'step2_trace'
        from src.plotting.spectra_plotter import plot_2d_image_to_file
        import matplotlib.pyplot as plt
        plt.figure(figsize=(12, 10))
        vmin, vmax = np.percentile(flat_data, [1, 99])
        plt.imshow(flat_data, aspect='auto', origin='lower', cmap='viridis', vmin=vmin, vmax=vmax)
        plt.colorbar(label='Counts')

        for i, ap_id in enumerate(all_ap_ids):
            ap = apertures.get_aperture(ap_id)
            center_pos = centers[i]
            lower_pos = lo_traces[i]
            upper_pos = hi_traces[i]
            
            is_filled = getattr(ap, 'is_interpolated', False)
            alpha_val = 0.7 if is_filled else 1.0
            
            # 中心轨迹线用间隔大一点的虚线，上下边缘用点线
            # 颜色与背景光谱的 color map 对比明显
            plt.plot(x_coords_all, center_pos, color='red', linestyle='--', linewidth=0.5, alpha=alpha_val)
            plt.plot(x_coords_all, lower_pos, color='orange', linestyle=':', linewidth=0.4, alpha=alpha_val)
            plt.plot(x_coords_all, upper_pos, color='magenta', linestyle=':', linewidth=0.4, alpha=alpha_val)
            
            if np.any(np.isfinite(center_pos)):
                x_anno = x_coords_all[np.isfinite(center_pos)][-1]
                y_anno = center_pos[np.isfinite(center_pos)][-1]
                plt.text(x_anno + 5, y_anno, str(ap_id), color='red', fontsize=5, ha='left', va='center', alpha=alpha_val, clip_on=False)

        plt.xlim(0, width - 1)
        plt.ylim(bottom=0)
        plt.xlabel('Pixel (X)')
        plt.ylabel('Pixel (Y)')
        plt.title('Order Traces on Combined Master Flat')
        plt.grid(False)
        traces_plot_file = out_dir / f'order_traces.{fig_format}'
        plt.savefig(str(traces_plot_file), dpi=150, bbox_inches='tight')
        plt.close()
        logger.info(f"Saved order traces plot to {traces_plot_file}")

        # Plot seed peaks diagnostic
        diags = processor.order_diagnostics
        if 'seed_profile' in diags and 'seeds' in diags:
            logger.info("Plotting seed peaks diagnostic...")
            seed_profile = diags['seed_profile']
            seeds = diags['seeds']
            envelope = diags.get('envelope')

            fig, ax = plt.subplots(figsize=(12, 8))
            y_coords = np.arange(len(seed_profile))
            if envelope is not None and 'profile' in diags:
                ax.plot(y_coords, diags['profile'], 'k-', linewidth=0.7, alpha=1.0, label='Raw Profile')
                ax.plot(y_coords, envelope, 'g--', linewidth=0.7, label='Background Envelope')
                
                det_threshold_flux = snr_threshold * diags.get('noise', 1.0)
                ax.plot(y_coords, envelope + det_threshold_flux, color='magenta', linestyle=':', linewidth=0.8, alpha=0.8, label=f'{snr_threshold} $\\sigma$ Detection Threshold')

            # Add markings for each aperture
            max_prof_val = np.nanmax(diags['profile']) if 'profile' in diags else 1000
            x1_col = diags.get('x1', max(0, width // 2 - 5))
            x2_col = diags.get('x2', min(width, width // 2 + 6))
            x_cols_eval = np.arange(x1_col, x2_col, dtype=float)
            for ap_id, ap in apertures.apertures.items():
                center_y = np.nanmedian(ap.get_position(x_cols_eval))
                lower_y = np.nanmedian(ap.get_lower(x_cols_eval))
                upper_y = np.nanmedian(ap.get_upper(x_cols_eval))

                if np.isfinite(center_y):
                    if len(seeds) > 0:
                        seed_idx = np.argmin(np.abs(seeds - center_y))
                        if np.abs(seeds[seed_idx] - center_y) < diags.get('med_sep', 30) * 0.75:
                            peak_y_coord = seeds[seed_idx]
                            if 'profile' in diags and envelope is not None and 0 <= peak_y_coord < len(diags['profile']):
                                peak_top = diags['profile'][peak_y_coord]
                                
                                # 红色短线和级次序号标记在 peaks 顶部
                                v_offset = 0.015 * max_prof_val
                                ax.vlines(peak_y_coord, peak_top + v_offset, peak_top + v_offset + 0.025 * max_prof_val, color='red', linewidth=1.5, zorder=10)
                                ax.text(peak_y_coord, peak_top + v_offset + 0.03 * max_prof_val, str(ap_id), color='red', ha='center', va='bottom', fontsize=7)
                                
                                # 青色水平 error bar，两边的 cap 底部卡在实际的级次边缘流量对应的位置
                                err_lower = center_y - lower_y
                                err_upper = upper_y - center_y
                                if err_lower > 0 and err_upper > 0:
                                    idx_lower = int(np.clip(round(lower_y), 0, len(diags['profile']) - 1))
                                    idx_upper = int(np.clip(round(upper_y), 0, len(diags['profile']) - 1))
                                    edge_flux = (diags['profile'][idx_lower] + diags['profile'][idx_upper]) / 2.0
                                    ax.errorbar(peak_y_coord, edge_flux, xerr=[[err_lower], [err_upper]], fmt='none', ecolor='cyan', capsize=2.0, elinewidth=0.6, capthick=0.4, zorder=10)

            ax.set_xlabel('Pixel (Y)')
            ax.set_ylabel('Flux (Counts)')
            ax.set_title('Order Seeds Detection')
            ax.legend()
            ax.grid(True, alpha=0.3)
            seed_peaks_plot_file = out_dir / f'order_seed_peaks.{fig_format}'
            plt.savefig(str(seed_peaks_plot_file), dpi=150, bbox_inches='tight')
            plt.close(fig)
            logger.info(f"Saved seed peaks plot to {seed_peaks_plot_file}")

            # Plot peak SNR diagnostic
            logger.info("Plotting peak SNR diagnostic...")
            noise = diags.get('noise', 1.0)
            local_sep_func = diags.get('local_sep_func')
            fig, ax1 = plt.subplots(figsize=(12, 8))
            
            plot_data = []
            for ap_id, ap in apertures.apertures.items():
                center_y = np.nanmedian(ap.get_position(x_cols_eval))
                if np.isfinite(center_y) and len(seeds) > 0:
                    seed_idx = np.argmin(np.abs(seeds - center_y))
                    if np.abs(seeds[seed_idx] - center_y) < diags.get('med_sep', 30) * 0.75:
                        s_y = seeds[seed_idx]
                        if 0 <= s_y < len(seed_profile):
                            snr_val = seed_profile[s_y] / noise
                            plot_data.append((s_y, snr_val, ap_id))

            plot_data.sort(key=lambda x: x[0])

            if len(plot_data) > 0:
                plot_seeds, plot_snrs, plot_labels = zip(*plot_data)
                ax1.plot(plot_seeds, plot_snrs, 'b-o', linewidth=1.0, markersize=4, label='Peak SNR')
                
                # 根据图像截面长度（完整空间方向跨度）稍微左右宽展几列
                x_margin = len(seed_profile) * 0.02
                ax1.set_xlim(-x_margin, len(seed_profile) - 1 + x_margin)
                
                # 提高 Y 轴上限，为顶部居中的图例留出空间
                snr_max = max(plot_snrs)
                ax1.set_ylim(bottom=-snr_max * 0.05, top=max(snr_threshold * 1.5, snr_max * 1.20))
            
            ax1.axhline(y=snr_threshold, color='magenta', linestyle=':', linewidth=1.2, alpha=0.8, label=f'{snr_threshold} $\\sigma$ Threshold')
            
            ax1.set_xlabel('Pixel (Y)')
            ax1.set_ylabel('Signal-to-Noise Ratio (SNR)', color='b')
            ax1.tick_params(axis='y', labelcolor='b')
            ax1.set_title('Peak SNR and Inter-Order Separation')
            ax1.grid(True, alpha=0.3, axis='x')
            
            if local_sep_func and len(plot_data) > 0:
                ax2 = ax1.twinx()
                separations = [local_sep_func(y) for y in plot_seeds]
                
                # Inter-order Separation 只画点
                ax2.plot(plot_seeds, separations, 'rs', markersize=3, alpha=0.7, label='Inter-order Separation')
                
                # 画出前面多项式拟合得到的光滑曲线
                x_min, x_max = ax1.get_xlim()
                y_eval = np.linspace(max(0, x_min), min(len(seed_profile)-1, x_max), 300)
                sep_eval = [local_sep_func(y) for y in y_eval]
                ax2.plot(y_eval, sep_eval, 'r-', linewidth=1.0, alpha=0.4, label='Smooth Fit Curve')
                
                ax2.set_ylabel('Inter-order Separation (pixels)', color='r')
                ax2.tick_params(axis='y', labelcolor='r')
                
                # 提高右侧 Y 轴上限，防止级次间距曲线撞到顶部图例
                sep_max = max(sep_eval) if len(sep_eval) > 0 else 10.0
                ax2.set_ylim(bottom=-sep_max * 0.05, top=sep_max * 1.20)
                
                snr_range = max(plot_snrs) - min(plot_snrs) if len(plot_snrs) > 1 else 10.0
                y_range = max(sep_eval) - min(sep_eval) if len(sep_eval) > 1 else 10.0
                for y, sep, snr, label in zip(plot_seeds, separations, plot_snrs, plot_labels):
                    ax1.text(y, snr + snr_range * 0.015, str(label), va='bottom', ha='center', color='blue', fontsize=7, rotation=90)
                    ax2.text(y, sep + y_range * 0.015, str(label), va='bottom', ha='center', color='red', fontsize=7, rotation=90)
                
                # 合并两个图的图例到一起
                lines1, labels1 = ax1.get_legend_handles_labels()
                lines2, labels2 = ax2.get_legend_handles_labels()
                ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper center', ncol=3)
            elif len(plot_data) > 0:
                snr_range = max(plot_snrs) - min(plot_snrs) if len(plot_snrs) > 1 else 10.0
                for y, snr, label in zip(plot_seeds, plot_snrs, plot_labels):
                    ax1.text(y, snr + snr_range * 0.015, str(label), va='bottom', ha='center', color='blue', fontsize=7, rotation=90)
                ax1.legend(loc='upper center', ncol=2)
            
            fig.tight_layout()
            snr_plot_file = out_dir / f'order_seed_snr.{fig_format}'
            plt.savefig(str(snr_plot_file), dpi=150, bbox_inches='tight')
            plt.close(fig)
            logger.info(f"Saved peak SNR plot to {snr_plot_file}")

    logger.info("Creating FlatField object...")
    flat_field = FlatField(
        flat_data=processor.flat_data,
        flat_mask=processor.flat_mask,
        flat_norm=None,
        flat_sens=None,
        scattered_light=None,
        smoothed_model=None,
        pixel_flat=None,
        illumination_flat=None,
        aperture_set=apertures,
    )

    logger.info("Step 2 order tracing completed successfully")
    return flat_field, apertures
