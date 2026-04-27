"""
Order tracing stage — Step 2 of the spectral reduction pipeline.

Combines master flat frames, detects echelle orders on the detector,
fits polynomial traces, and builds the pixel-flat / blaze-profile models
that are consumed by later pipeline stages.

The main entry point is process_order_tracing_stage().  FlatFieldProcessor is the
class that wraps all processing state.
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

from src.plotting.spectra_plotter import plot_2d_image_to_file

logger = logging.getLogger(__name__)


class FlatFieldProcessor:
    """Handles flat field processing and order tracing."""

    def __init__(self):
        """
        Initialize flat field processor.
        """
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
        self.flat_sens = None
        self.order_diagnostics = {}

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

    def detect_orders(self, snr_threshold: float = 5.0,
                      step_denominator: int = 20,
                      gap_fill_factor: float = 1.35,
                      gap_fill_snr: float = 2.5,
                      min_trace_coverage: float = 0.20,
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
            trace_degree: Polynomial degree used to fit order trace.
            separation: Estimated order separation in pixels.
            gap_fill_factor: Factor for gap fill detection.
            gap_fill_snr: SNR threshold for gap fill.
            min_trace_coverage: Minimum fraction of detector width a traced
                order must cover to be accepted during Step 2 detection.
            trace_degree: Polynomial degree for center tracing.
            output_dir_base: Output directory for diagnostics.

        Returns:
            ApertureSet with detected orders
        """
        if self.flat_data is None:
            raise RuntimeError("No flat data")

        logger.info("Detecting grating orders (simple tracing)...")

        # Use the simplified grating order tracing method
        apertures = _find_grating_orders_impl(
            data=self.flat_data,
            mask=self.flat_mask,
            snr_threshold=snr_threshold,
            step_denominator=step_denominator,
            gap_fill_factor=gap_fill_factor,
            gap_fill_snr=gap_fill_snr, 
            min_trace_coverage=min_trace_coverage,
            trace_degree=trace_degree,
            output_dir_base=output_dir_base
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

    def save_flat_field(self, output_path: str,
                        save_plots: bool = True,
                        fig_format: str = 'png'):
        """Save flat field data to FITS file.

        MasterFlat.fits contains only the combined flat image (PRIMARY) and
        the bad-pixel MASK extension so its size is comparable to a single
        input flat frame.  All derived products are saved as separate FITS
        files alongside MasterFlat.fits.
        """
        if self.flat_data is None:
            raise RuntimeError("No flat data to save")

        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        out_dir = output_path.parent

        from astropy.io import fits

        # ---- MasterFlat.fits: combined image + mask only ----
        primary_hdu = fits.PrimaryHDU(data=self.flat_data.astype(np.float32))
        primary_hdu.header['FLATCR'] = (True, 'Flat fielding completed')
        hdul = fits.HDUList([
            primary_hdu,
            fits.ImageHDU(data=self.flat_mask, name='MASK'),
        ])
        hdul.writeto(str(output_path), overwrite=True)
        logger.info(f"Saved master flat to {output_path}")

        # ---- derived products as separate files ----
        def _save_ext(data, name: str):
            if data is None:
                return
            dest = out_dir / name
            hdu = fits.PrimaryHDU(data=data.astype(np.float32) if hasattr(data, 'astype') else data)
            hdu.writeto(str(dest), overwrite=True)
            logger.info(f"Saved {name}")

        if self.flat_norm is not None:
            _save_ext(self.flat_norm, 'master_flat_norm.fits')
        if self.flat_sens is not None:
            _save_ext(self.flat_sens, 'master_flat_sensitivity.fits')
        if self.response_map is not None:
            _save_ext(self.response_map, 'master_flat_response.fits')
        if self.scattered_light is not None:
            _save_ext(self.scattered_light, 'master_flat_scatter.fits')
        if self.smoothed_model is not None:
            _save_ext(self.smoothed_model, 'master_flat_smooth_model.fits')
        if self.pixel_flat is not None:
            _save_ext(self.pixel_flat, 'master_flat_pixel_flat.fits')
        if self.flat_corr_2d is not None:
            _save_ext(self.flat_corr_2d, 'master_flat_flat_corr_2d.fits')

        # Save diagnostic plots if enabled
        if save_plots:
            out_dir = Path(output_path).parent
            # Plot master flat
            plot_2d_image_to_file(self.flat_data, str(out_dir / f'MasterFlat.{fig_format}'), "Master Flat Field")


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


def process_order_tracing_stage(flat_filenames: List[str],
                      output_dir_base: str = './midpath',
                      # Flat combination parameters
                      combine_method: str = 'median',
                      combine_sigma: float = 3.0,
                      mosaic_maxcount: float = 65535,
                      # Order detection parameters
                      snr_threshold: float = 5.0,
                      step_denominator: int = 20,
                      gap_fill_factor: float = 1.6,
                      gap_fill_snr: float = 2.5,
                      min_trace_coverage: float = 0.20,
                      trace_degree: int = 4,
                      # Profile boundary parameters
                      boundary_frac: float = 0.02,
                      fwhm_scale: float = 1.5,
                      width_cheb_degree: int = 3,
                      # Gap fill / extend parameters
                      n_extend_below: int = 0,
                      n_extend_above: int = 0,
                      gap_fill_factor_interp: float = 1.35,
                      # Diagnostic parameters
                      save_plots: bool = True,
                      fig_format: str = 'png') -> Tuple[FlatField, ApertureSet]:
    """
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

    from src.plotting.spectra_plotter import plot_2d_image_to_file

    # Detect orders — uses raw combined flat, no normalization needed
    logger.info("Detecting grating orders...")
    apertures = processor.detect_orders(
        snr_threshold=snr_threshold,
        step_denominator=step_denominator,
        gap_fill_factor=gap_fill_factor,
        gap_fill_snr=gap_fill_snr,
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
        """Estimate background noise via iterative MAD sigma-clipping.

        After subtracting the smooth background envelope, the residual
        contains near-zero background pixels and positive order-peak
        pixels.  We iteratively clip bright outliers (order peaks) to
        converge on the pure background noise level.

        Steps:
          1. residual = data − envelope
          2. Initial aggressive cut: discard pixels > median + 2×MAD
             (removes obvious order peaks)
          3. Iterate 3σ upper clipping until the pixel set stabilises
          4. σ_bg = 1.4826 × MAD of the surviving background pixels
        """
        residual = cross_section_1d - bg_envelope_1d
        valid = np.isfinite(residual)
        res = residual[valid]

        if res.size < 5:
            # Legacy SNR-based boundary logic removed
            return float(np.nanstd(residual))

        # Use bottom 25% of pixels to seed the background estimate robustly
        sorted_res = np.sort(res)
        bg_pixels = sorted_res[:max(5, len(sorted_res) // 4)]
        med = float(np.nanmedian(bg_pixels))
        mad = float(np.nanmedian(np.abs(bg_pixels - med)))
        sigma_est = max(1.4826 * mad, 1e-6)
        keep = res < (med + 3.0 * sigma_est)

        # --- Iterative 3σ upper clipping (converge on pure background) ---
        for _ in range(5):
            if np.sum(keep) < 5:
                break
            sub = res[keep]
            med = float(np.nanmedian(sub))
            mad = float(np.nanmedian(np.abs(sub - med)))
            sigma_est = max(1.4826 * mad, 1e-6)
            new_keep = res < (med + 3.0 * sigma_est)
            if np.array_equal(new_keep, keep):
                break
            keep = new_keep

        # Final noise estimate from surviving background pixels
        if np.sum(keep) >= 5:
            sub = res[keep]
            med = float(np.nanmedian(sub))
            mad = float(np.nanmedian(np.abs(sub - med)))
            sigma_bg = max(1.4826 * mad, 1e-6)
        else:
            sigma_bg = max(sigma_est, 1e-6)

        # Legacy SNR-based boundary logic removed
        return sigma_bg

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
        from concurrent.futures import ThreadPoolExecutor

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

    # Float column coordinates (needed by hybrid smoothing below)
    x_all_f = x_coords_all.astype(float)

    # ---------------------------------------------------------------------------
    # Hybrid trace smoothing:
    #   1. Center trace  Y_center(X)  — B-spline with fixed knots (sub-pixel
    #      precision tracking of optical dispersion curve).
    #   2. Aperture widths  W_up(X) = Y_upper_raw − Y_center_smooth,
    #      W_low(X) = Y_center_smooth − Y_lower_raw  — low-order Chebyshev
    #      (rigid, ignores edge-detection noise; captures only slow optical
    #      aberration like astigmatism / FRD).
    #   3. Reconstruct final smooth boundaries:
    #        Y_upper_final = Y_center_bspline + W_up_chebyshev
    #        Y_lower_final = Y_center_bspline − W_low_chebyshev
    # ---------------------------------------------------------------------------
    from scipy.interpolate import splrep, splev
    from numpy.polynomial import chebyshev as cheb

    # Config knobs
    # Config knobs (now passed as function parameters)

    logger.info(f"Hybrid smoothing: Chebyshev center (deg {trace_degree}) + Chebyshev width (deg {width_cheb_degree})")

    # Collect Chebyshev width coefs + Chebyshev-smoothed centres per order for saving.
    _trace_coefs = {}   # aperture_id → {center_arr, w_up_cheb, w_low_cheb}

    for row_idx, item in enumerate(valid_orders):
        aperture = item['aperture']
        valid_center = item['valid_center']

        # --- (a) Raw center positions from the polynomial trace ---
        center_raw = aperture.get_position(x_all_f)   # shape (width,)

        # Columns where we have both a valid center *and* valid boundaries
        y_lo_raw = y_lo_all[row_idx, :]
        y_hi_raw = y_hi_all[row_idx, :]
        valid_bnd = valid_center & np.isfinite(y_lo_raw) & np.isfinite(y_hi_raw)
        n_valid = int(np.sum(valid_bnd))
        if n_valid < 10:
            continue

        # ---- Step 1: Center trace is already smooth and regularized ----
        # Do not re-fit it on a partial valid_bnd segment, which causes extrapolation explosion!
        center_smooth = center_raw

        # ---- Step 2: Compute widths relative to smooth centre ----
        w_up_raw  = y_hi_raw[valid_bnd] - center_smooth[valid_bnd]   # positive
        w_low_raw = center_smooth[valid_bnd] - y_lo_raw[valid_bnd]   # positive

        # ---- Step 3: Chebyshev fit the widths ----
        x_bnd = x_all_f[valid_bnd]
        span = x_bnd[-1] - x_bnd[0] if len(x_bnd) > 0 else 0
        # If valid boundaries cover less than half the image, restrict width fit degree
        # to prevent widths from exploding outwards in the faint wings.
        if span < 0.3 * width:
            deg = 0
        elif span < 0.6 * width:
            deg = min(1, width_cheb_degree, n_valid - 1)
        else:
            deg = min(width_cheb_degree, n_valid - 1)

        coef_w_up  = cheb.chebfit(x_bnd, w_up_raw,  deg)
        coef_w_low = cheb.chebfit(x_bnd, w_low_raw, deg)

        w_up_smooth  = cheb.chebval(x_all_f, coef_w_up)
        w_low_smooth = cheb.chebval(x_all_f, coef_w_low)

        # Ensure widths are non-negative
        w_up_smooth  = np.maximum(w_up_smooth, 0.0)
        w_low_smooth = np.maximum(w_low_smooth, 0.0)

        # ---- Step 4: Reconstruct final boundaries ----
        hi_final = center_smooth + w_up_smooth
        lo_final = center_smooth - w_low_smooth

        # Clamp to image extent
        hi_final = np.clip(hi_final, 0, height - 1)
        lo_final = np.clip(lo_final, 0, height - 1)

        # Write back (only where centre is valid)
        y_lo_all[row_idx, :] = np.where(valid_center, lo_final, np.nan)
        y_hi_all[row_idx, :] = np.where(valid_center, hi_final, np.nan)

        # Collect coefs for saving
        _trace_coefs[item['aperture_id']] = {
            'center_arr': center_smooth.tolist(),
            'w_up_cheb':  coef_w_up.tolist(),
            'w_low_cheb': coef_w_low.tolist(),
        }

    logger.info("Hybrid smoothing complete (Chebyshev center + Chebyshev width → traces)")

    # Write hybrid-smoothed results back to aperture objects so that
    # downstream gap-fill interpolation uses the precise narrow boundaries
    # (not the wide initial-detection polynomials).
    for item in valid_orders:
        aid = item['aperture_id']
        if aid not in _trace_coefs:
            continue
        bc = _trace_coefs[aid]
        ap = item['aperture']
        ap.center_arr = np.asarray(bc['center_arr'], dtype=float)
        ap.w_up_cheb_coef = np.asarray(bc['w_up_cheb'], dtype=float)
        ap.w_low_cheb_coef = np.asarray(bc['w_low_cheb'], dtype=float)

    # ---- Fill interior gaps & extend edges (spacing-based) ----
    # Gap-filled orders are treated as real orders and included in the mask
    # and boundary coefs output files.
    
    filled_apertures = fill_missing_orders_by_interpolation(
        apertures, width, height, trace_degree=trace_degree,
        n_extend_below=n_extend_below, n_extend_above=n_extend_above,
        flat_data=flat_data, gap_fill_factor=gap_fill_factor_interp,
        step_denominator=step_denominator)

    if filled_apertures.norders > apertures.norders:
        n_new = filled_apertures.norders - apertures.norders
        logger.info(f"Gap-fill added {n_new} predicted order(s): "
                    f"{apertures.norders} → {filled_apertures.norders}")

        # Identify the newly added order IDs.
        existing_ids = set(apertures.get_orders())
        new_ids = sorted(set(filled_apertures.get_orders()) - existing_ids)

        # For each new order compute profile-root boundaries and apply hybrid
        # smoothing, identical to what was done for originally detected orders.
        for new_oid in new_ids:
            new_ap = filled_apertures.get_aperture(new_oid)
            center_new = new_ap.get_position(x_all_f)
            valid_c = np.isfinite(center_new) & (center_new >= 0) & (center_new < height)

            # Compute profile-root boundaries for the filled order
            lo_new = np.full(width, np.nan, dtype=float)
            hi_new = np.full(width, np.nan, dtype=float)

            for xi in range(width):
                cy = float(center_new[xi])
                if not valid_c[xi] or cy < 0 or cy >= height:
                    continue
                # Use the polynomial lower/upper as initial boundaries
                lo_poly = float(new_ap.get_lower(float(xi)))
                hi_poly = float(new_ap.get_upper(float(xi)))
                if np.isfinite(lo_poly) and np.isfinite(hi_poly):
                    lo_new[xi] = float(np.clip(lo_poly, 0, height - 1))
                    hi_new[xi] = float(np.clip(hi_poly, 0, height - 1))

            valid_bnd = valid_c & np.isfinite(lo_new) & np.isfinite(hi_new)
            n_valid_new = int(np.sum(valid_bnd))
            if n_valid_new < 10:
                # Not enough valid columns — store polynomial-derived boundaries
                _trace_coefs[new_oid] = {
                    'center_arr': center_new.tolist(),
                    'w_up_cheb':  cheb.chebfit(x_all_f[valid_c], np.maximum(hi_new[valid_c] - center_new[valid_c], 0), min(width_cheb_degree, max(int(np.sum(valid_c))-1, 0))).tolist() if np.sum(valid_c) > 0 else [0.0],
                    'w_low_cheb': cheb.chebfit(x_all_f[valid_c], np.maximum(center_new[valid_c] - lo_new[valid_c], 0), min(width_cheb_degree, max(int(np.sum(valid_c))-1, 0))).tolist() if np.sum(valid_c) > 0 else [0.0],
                }
                # Add to mask arrays
                valid_orders.append({
                    'aperture_id': new_oid,
                    'aperture': new_ap,
                    'valid_center': valid_c,
                    'gap_filled': True,
                })
                new_row = np.where(valid_c, lo_new, np.nan)
                y_lo_all = np.vstack([y_lo_all, new_row[np.newaxis, :]])
                new_row_hi = np.where(valid_c, hi_new, np.nan)
                y_hi_all = np.vstack([y_hi_all, new_row_hi[np.newaxis, :]])
                continue

            # ---- Center trace is already smooth (interpolated from valid orders) ----
            c_smooth_new = center_new

            # ---- Chebyshev fit the widths ----
            w_up_r = hi_new[valid_bnd] - c_smooth_new[valid_bnd]
            w_lo_r = c_smooth_new[valid_bnd] - lo_new[valid_bnd]
            x_bnd_n = x_all_f[valid_bnd]
            span_n = x_bnd_n[-1] - x_bnd_n[0] if len(x_bnd_n) > 0 else 0
            if span_n < 0.3 * width:
                deg_n = 0
            elif span_n < 0.6 * width:
                deg_n = min(1, width_cheb_degree, n_valid_new - 1)
            else:
                deg_n = min(width_cheb_degree, n_valid_new - 1)

            coef_wu = cheb.chebfit(x_bnd_n, w_up_r, deg_n)
            coef_wl = cheb.chebfit(x_bnd_n, w_lo_r, deg_n)

            w_up_s = np.maximum(cheb.chebval(x_all_f, coef_wu), 0.0)
            w_lo_s = np.maximum(cheb.chebval(x_all_f, coef_wl), 0.0)

            hi_f = np.clip(c_smooth_new + w_up_s, 0, height - 1)
            lo_f = np.clip(c_smooth_new - w_lo_s, 0, height - 1)

            # Store coefs
            _trace_coefs[new_oid] = {
                'center_arr': c_smooth_new.tolist(),
                'w_up_cheb':  coef_wu.tolist(),
                'w_low_cheb': coef_wl.tolist(),
            }

            # Add to mask arrays
            valid_orders.append({
                'aperture_id': new_oid,
                'aperture': new_ap,
                'valid_center': valid_c,
                'gap_filled': True,
            })
            y_lo_all = np.vstack([y_lo_all, np.where(valid_c, lo_f, np.nan)[np.newaxis, :]])
            y_hi_all = np.vstack([y_hi_all, np.where(valid_c, hi_f, np.nan)[np.newaxis, :]])

        apertures = filled_apertures
        logger.info(f"Processed {len(new_ids)} gap-filled orders for mask & coefs")
    else:
        logger.info("No missing interior orders detected")

    # ---------------------------------------------------------------------------
    # Save Chebyshev boundary coefficients to companion JSON file.
    # ---------------------------------------------------------------------------
    import json as _json
    coefs_path = Path(base_output_path) / 'step2_trace' / 'Orders_trace_coefs.json'
    with open(coefs_path, 'w') as f:
        _json.dump({'orders': _trace_coefs}, f)
    logger.info(f"Saved trace coefs ({len(_trace_coefs)} orders) to {coefs_path.name}")

    # ---------------------------------------------------------------------------
    # Build and save the labeled order mask (Orders_mask.fits).
    # Pixel value = order ID (int32), 0 = background.
    # Includes both detected and gap-filled orders.
    # ---------------------------------------------------------------------------
    order_labels = np.zeros((height, width), dtype=np.int32)
    for row_idx, item in enumerate(valid_orders):
        oid = int(item['aperture_id'])
        for xi in range(width):
            lo = y_lo_all[row_idx, xi]
            hi = y_hi_all[row_idx, xi]
            if np.isfinite(lo) and np.isfinite(hi):
                r0 = max(0, int(np.floor(lo)))
                r1 = min(height, int(np.ceil(hi)) + 1)
                if r1 > r0:
                    order_labels[r0:r1, xi] = oid

    labels_path = Path(base_output_path) / 'step2_trace' / 'Orders_mask.fits'
    write_fits_image(str(labels_path), order_labels, dtype='int32')
    logger.info(f"Saved labeled order mask to {labels_path} "
                f"({len(valid_orders)} orders, "
                f"{np.count_nonzero(order_labels)}/{order_labels.size} px labeled)")

    # Step 3 reads Orders_mask.fits for background subtraction; Step 4/5
    # reconstruct ApertureSet from the mask + coefs via load_boundary_coefs().

    # Step 2 should only emit tracing diagnostics.
    rows, cols = height, width
    # Step4 writes flat-model and 2D flat-correction artifacts.
    if save_plots:
        logger.info("Generating diagnostic plots...")
        out_dir = Path(base_output_path) / 'step2_trace'
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
        # height, width, center_col, aperture_boundary_snr, _estimate_noise_floor,
        # valid_orders, y_lo_all, and y_hi_all are all pre-computed before mask building.
        q25_col = int(np.clip(round(width * 0.25), 0, width - 1))
        q50_col = int(np.clip(round(width * 0.50), 0, width - 1))
        q75_col = int(np.clip(round(width * 0.75), 0, width - 1))
        diag_cols = [q25_col, q50_col, q75_col]

        has_interpolated = False
        for idx, item in enumerate(valid_orders):
            aperture = item['aperture']
            valid_center = item['valid_center']
            center_pos = aperture.get_position(x_coords_all)
            is_filled = getattr(aperture, 'is_interpolated', False) or item.get('gap_filled', False)
            snr_val = None
            if 'peaks' in locals() and idx < len(peaks):
                snr_val = peaks[idx]
            # SNR大于阈值的order一律画实线且不降透明度
            if is_filled and snr_val is not None and snr_val >= snr_threshold:
                ls = '-'
                alpha_val = 1.0
            else:
                ls = '--' if is_filled else '-'
                alpha_val = 1.0
            plt.plot(x_coords_all[valid_center], center_pos[valid_center],
                     color='red', linestyle=ls, linewidth=0.6, alpha=alpha_val,
                     label=f'Order {item["aperture_id"]}')
            x_valid = x_coords_all[valid_center]
            y_valid = center_pos[valid_center]
            if x_valid.size > 0:
                x_anno = int(x_valid[-1])
                y_anno = float(y_valid[-1])
                plt.text(x_anno + 3, y_anno,
                         str(item['aperture_id']),
                         color='red', fontsize=3.5, ha='left', va='center', alpha=alpha_val)

        # Draw pre-computed smoothed profile-root boundaries for all valid orders.
        # y_lo_all / y_hi_all have already been Chebyshev-smoothed above.
        if len(valid_orders) >= 1:
            for row_idx, item in enumerate(valid_orders):
                aperture = item['aperture']
                valid_center = item['valid_center']
                y_lo = y_lo_all[row_idx, :]
                y_hi = y_hi_all[row_idx, :]
                valid_bounds = valid_center & np.isfinite(y_lo) & np.isfinite(y_hi)
                if np.sum(valid_bounds) < max(10, int(0.10 * width)):
                    continue

                y_lo_plot = np.where(valid_bounds, y_lo, np.nan)
                y_hi_plot = np.where(valid_bounds, y_hi, np.nan)

                is_filled = getattr(aperture, 'is_interpolated', False) or item.get('gap_filled', False)
                ls = '--' if is_filled else '-'
                alpha_white = 0.25 if is_filled else 0.45
                alpha_cyan = 0.50 if is_filled else 0.90
                
                plt.plot(x_coords_all[valid_bounds], y_lo_plot[valid_bounds], color='white', linestyle=ls, linewidth=0.7, alpha=alpha_white)
                plt.plot(x_coords_all[valid_bounds], y_hi_plot[valid_bounds], color='white', linestyle=ls, linewidth=0.7, alpha=alpha_white)
                plt.plot(x_coords_all[valid_bounds], y_lo_plot[valid_bounds], color='#00E5FF', linestyle=ls, linewidth=0.5, alpha=alpha_cyan)
                plt.plot(x_coords_all[valid_bounds], y_hi_plot[valid_bounds], color='#00E5FF', linestyle=ls, linewidth=0.5, alpha=alpha_cyan)

                # Mark q25/q50/q75 positions to visually cross-check against
                # apertures_q25/q50/q75 diagnostics.
                ap = item['aperture']
                for xc in diag_cols:
                    if not valid_center[xc] or not np.isfinite(y_lo[xc]) or not np.isfinite(y_hi[xc]):
                        continue
                    yc = float(np.clip(ap.get_position(xc), 0, height - 1))
                    y_lo_c = float(np.clip(y_lo[xc], 0, height - 1))
                    y_hi_c = float(np.clip(y_hi[xc], 0, height - 1))
                    # Short horizontal lines for lower and upper boundaries
                    plt.plot([xc - 6, xc + 6], [y_lo_c, y_lo_c],
                             color='#FF00FF', linewidth=0.9, alpha=0.90, zorder=6)
                    plt.plot([xc - 6, xc + 6], [y_hi_c, y_hi_c],
                             color='lime', linewidth=0.9, alpha=0.90, zorder=6)
                    # Short horizontal line marks the center position
                    plt.plot([xc - 6, xc + 6], [yc, yc],
                             color='#FFA500', linewidth=0.9, alpha=0.95, zorder=7)

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
        valid_ap_ids = {item['aperture_id'] for item in valid_orders}

        def _plot_apertures_cross_section(col_idx: int, out_name: str):
            plt.figure(figsize=(10, 8))
            half_step = max(1, step_denominator // 2)
            c1 = max(0, col_idx - half_step)
            c2 = min(width, col_idx + half_step + 1)
            cross_section = np.median(flat_data[:, c1:c2], axis=1).astype(float)
            y_coords = np.arange(height)
            plt.plot(y_coords, cross_section, color='0.25', linewidth=0.75, label='Cross-section')
            cs_max = float(np.nanmax(cross_section))
            cs_span = max(cs_max - float(np.nanmin(cross_section)), 1e-6)

            ordered_local = sorted(apertures.apertures.items(), key=lambda x: x[0])

            # Pre-collect center positions (clipped to detector) for envelope + window computation
            order_centers_y = []
            for _ap_id, _ap in ordered_local:
                _cy = float(_ap.get_position(col_idx))
                if np.isfinite(_cy) and -5 <= _cy <= height + 5:
                    order_centers_y.append(float(np.clip(_cy, 0, height - 1)))

            if len(order_centers_y) == 0:
                # No orders on detector at this column — just save blank cross-section
                plt.xlabel('Pixel (Y)')
                plt.ylabel('Counts')
                plt.title(f'Apertures Cross-section at X~{col_idx} (median of {c2-c1} cols, no orders)')
                plt.grid(True, alpha=0.3)
                out_file = out_dir / out_name
                plt.savefig(str(out_file), dpi=150, bbox_inches='tight')
                plt.close()
                return

            # Fit smooth background envelope through inter-order valley minima
            bg_envelope = _fit_background_envelope_1d(cross_section, order_centers_y)
            plt.plot(y_coords, bg_envelope, color='orange', linewidth=0.9,
                     linestyle='--', alpha=0.85, label='BG envelope', zorder=2)

            # Estimate background noise on inter-order regions using robust MAD.
            noise_floor = _estimate_noise_floor(cross_section, bg_envelope, order_centers_y)

            # Draw detection threshold curve (used for seeding, not for boundaries).
            if np.isfinite(noise_floor) and noise_floor > 0:
                noise_sigma_raw = noise_floor
                det_thresh = bg_envelope + snr_threshold * noise_sigma_raw
                plt.plot(y_coords, det_thresh, color='#FF00FF', linewidth=0.9,
                         linestyle=':', alpha=0.85,
                         label=f'Blind Detection ({snr_threshold:.1f}σ)', zorder=3)

            # Build a set of predicted (surface-filled) aperture IDs.
            predicted_ap_ids = {ap_id for ap_id in apertures.get_orders()
                                if ap_id not in valid_ap_ids and ap_id >= 0}

            # Build a list of (aperture_id, aperture, index_in_order_centers_y)
            # that are actually on-detector, so the search-window logic uses the
            # correct neighbour positions from order_centers_y.
            vis_items = []
            seq_counter = 0
            for aperture_id, aperture in ordered_local:
                center_y = float(aperture.get_position(col_idx))
                if not np.isfinite(center_y) or center_y < -5 or center_y > height + 5:
                    continue
                vis_items.append((seq_counter, aperture_id, aperture, center_y))
                seq_counter += 1

            for local_j, (seq_i, aperture_id, aperture, center_y) in enumerate(vis_items):
                cy_y = float(np.clip(center_y, 0, height - 1))
                cy_idx = int(round(cy_y))

                # Moffat boundaries stored in the aperture
                lo_y = float(np.clip(aperture.get_lower(col_idx), 0, height - 1))
                hi_y = float(np.clip(aperture.get_upper(col_idx), 0, height - 1))
                peak_val = float(cross_section[cy_idx])

                # Place error bar at the average cross-section value at the boundaries
                bnd_level = 0.5 * (float(np.interp(lo_y, y_coords, cross_section))
                                   + float(np.interp(hi_y, y_coords, cross_section)))

                # Horizontal error bar spanning center ± 1.5 FWHM (Moffat boundary)
                xerr_lo = max(0.0, cy_y - lo_y)
                xerr_hi = max(0.0, hi_y - cy_y)
                plt.errorbar(
                    cy_y, bnd_level,
                    xerr=[[xerr_lo], [xerr_hi]],
                    fmt='none', color='#00E5FF',
                    capsize=2.5, capthick=0.9, elinewidth=0.9, alpha=0.90,
                    zorder=4,
                )

                # Short vertical tick: solid for detected, dashed for predicted (same red color)
                tick_gap = 0.012 * cs_span
                tick_len = 0.040 * cs_span
                tick_bot = peak_val + tick_gap
                tick_top = tick_bot + tick_len
                is_predicted = aperture_id in predicted_ap_ids
                tick_color = 'red'
                # SNR大于阈值的order一律画实线
                tick_style = '-'
                if is_predicted:
                    snr_val = None
                    if 'ap_ids' in locals() and aperture_id in ap_ids:
                        idx = ap_ids.index(aperture_id)
                        if 'peaks' in locals() and idx < len(peaks):
                            snr_val = peaks[idx]
                    if snr_val is not None and snr_val >= snr_threshold:
                        tick_style = '-'
                    else:
                        tick_style = '--'
                plt.plot([cy_y, cy_y], [tick_bot, tick_top],
                         color=tick_color, linewidth=1.1, linestyle=tick_style,
                         alpha=0.95, zorder=5)

                # Order label just above the tick
                plt.text(cy_y, tick_top + 0.008 * cs_span,
                         str(aperture_id),
                         color=tick_color, fontsize=5, ha='center', va='bottom')

            plt.xlabel('Pixel (Y)')
            plt.ylabel('Counts')
            plt.title(f'Apertures Cross-section at X~{col_idx} (median of {c2-c1} cols)\n'
                      f'Orange = BG env, Magenta = Blind SNR, Cyan = Gap SNR\n'
                      f'Cyan = Moffat bnd (±1.5 FWHM), Red = center')
            if np.isfinite(noise_floor):
                plt.text(0.01, 0.985,
                         f'detection SNR={snr_threshold:.1f}, noise floor={noise_floor:.2f} counts',
                         transform=plt.gca().transAxes,
                         fontsize=7, color='#00AACC', ha='left', va='top')
            plt.legend(fontsize=7)
            plt.grid(True, alpha=0.3)
            out_file = out_dir / out_name
            plt.savefig(str(out_file), dpi=150, bbox_inches='tight')
            plt.close()
            logger.info(f"Saved apertures plot to {out_file}")

        # Output apertures diagnostics at 25%, 50%, 75%
        _plot_apertures_cross_section(q25_col, f'apertures_q25.{fig_format}')
        _plot_apertures_cross_section(q50_col, f'apertures_q50.{fig_format}')
        _plot_apertures_cross_section(q75_col, f'apertures_q75.{fig_format}')

        def _plot_order_metrics(col_idx: int, out_name: str):
            ordered_local = sorted(apertures.apertures.items(), key=lambda x: x[0])
            if not ordered_local:
                return

            ap_ids = []
            peaks = []
            centers = []
            
            half_step = max(1, step_denominator // 2)
            c1 = max(0, col_idx - half_step)
            c2 = min(width, col_idx + half_step + 1)
            cross_section = np.median(flat_data[:, c1:c2], axis=1).astype(float)
            
            for ap_id, ap in ordered_local:
                cy = float(ap.get_position(col_idx))
                if cy < 0 or cy >= height:
                    continue
                ap_ids.append(ap_id)
                centers.append(cy)

            if len(ap_ids) < 2:
                return

            bg_envelope = _fit_background_envelope_1d(cross_section, centers)
            noise_threshold = _estimate_noise_floor(cross_section, bg_envelope, centers)
            # Reverse the threshold multiplier to get raw noise sigma for standard SNR plot
            noise_sigma = noise_threshold
            
            if not np.isfinite(noise_sigma) or noise_sigma <= 0:
                noise_sigma = 1.0

            for cy in centers:
                cy_idx = int(round(cy))
                peak_val = cross_section[cy_idx]
                bg_val = bg_envelope[cy_idx]
                # Formula: (Peak - Background) / Noise_Sigma
                snr = (peak_val - bg_val) / noise_sigma
                peaks.append(snr)
            
            if len(ap_ids) < 2:
                return

            gaps = []
            for i in range(len(centers) - 1):
                gaps.append(abs(centers[i+1] - centers[i]))
            # Pad the last gap with the previous one to match array lengths
            gaps.append(gaps[-1])

            # Read detection SNR threshold from config for the threshold line
            det_snr = snr_threshold

            fig, ax1 = plt.subplots(figsize=(10, 6))

            color1 = 'tab:blue'
            ax1.set_xlabel('Order ID')
            ax1.set_ylabel('Peak SNR (above BG)', color=color1)
            ax1.plot(ap_ids, peaks, marker='o', color=color1, markersize=4,
                     linewidth=0.8, label='Peak SNR')
            for ap_id, snr in zip(ap_ids, peaks):
                ax1.annotate(str(ap_id), (ap_id, snr),
                             xytext=(0, 4), textcoords='offset points',
                             fontsize=5, color=color1,
                             ha='center', va='bottom')

            # Draw detection threshold line
            ax1.axhline(y=det_snr, color='orange', linestyle='--', linewidth=1.2,
                        alpha=0.85, label=f'Detection Threshold ({det_snr:.1f}σ)')
            ax1.text(ap_ids[0], det_snr * 1.08,
                     f'{det_snr:.1f}σ threshold', color='orange',
                     fontsize=7, va='bottom')

            ax1.tick_params(axis='y', labelcolor=color1)
            ax1.grid(True, alpha=0.3)

            ax2 = ax1.twinx()  
            color2 = 'tab:red'
            ax2.set_ylabel('Inter-order Gap (px)', color=color2)
            # Plot gap data points (严格只画N-1个)
            real_ap_ids = ap_ids[:-1]
            real_gaps = gaps[:-1]
            ax2.plot(real_ap_ids, real_gaps, marker='x', linestyle='none', color=color2, markersize=4)
            # Label each gap point with其order ID
            for ap_id_lbl, gap_val in zip(real_ap_ids, real_gaps):
                ax2.annotate(str(ap_id_lbl), (ap_id_lbl, gap_val),
                             xytext=(0, 5), textcoords='offset points',
                             fontsize=5, color=color2,
                             ha='center', va='bottom')
            # Fit and plot a smooth polynomial curve through the gap data
            if len(real_ap_ids) >= 4:
                try:
                    _x = np.array(real_ap_ids, dtype=float)
                    _y = np.array(real_gaps, dtype=float)
                    _deg = min(3, len(_x) - 1)
                    _coeffs = np.polyfit(_x, _y, _deg)
                    _x_smooth = np.linspace(_x[0], _x[-1], 400)
                    ax2.plot(_x_smooth, np.polyval(_coeffs, _x_smooth),
                             color='salmon', linewidth=1.2, linestyle='-', alpha=0.85)
                except Exception:
                    pass
            ax2.tick_params(axis='y', labelcolor=color2)

            ax1.set_title(
                f'Order Metrics at X~{col_idx} (median of {c2-c1} cols)  |  '
                f'noise σ={noise_sigma:.1f} counts  |  '
                f'threshold={det_snr:.1f}σ',
                fontsize=9,
            )
            fig.tight_layout()
            
            out_file = out_dir / out_name
            plt.savefig(str(out_file), dpi=150, bbox_inches='tight')
            plt.close()
            logger.info(f"Saved order metrics plot to {out_file}")

        _plot_order_metrics(q50_col, f'order_metrics_q50.{fig_format}')

    # Create FlatField object
    logger.info("Creating FlatField object...")
    flat_field = FlatField(
        flat_data=flat_data,
        flat_sens=None,
        flat_norm=None,
        flat_mask=flat_mask,
        scattered_light=None,
        smoothed_model=None,
        pixel_flat=None,
        illumination_flat=None,
        flat_corr_2d=None,
        aperture_set=apertures,
        cross_profiles=None,
        blaze_profiles=None,
    )

    logger.info("Step 2 order tracing completed successfully")
    return flat_field, apertures


# ======================================================================
# MERGED ALGORITHM KERNEL FROM grating_trace_simple.py
# ======================================================================


"""
Simple grating order tracing for spectral reduction pipeline.

Handles order detection and tracing for grating spectrographs,
where orders run horizontally (left to right) across the detector.
"""

import numpy as np
import logging
from typing import Optional
from scipy.signal import find_peaks, peak_widths
from scipy.ndimage import gaussian_filter1d, sobel
from scipy.interpolate import UnivariateSpline
from numpy.polynomial import Chebyshev, Polynomial
from src.core.data_structures import ApertureSet, ApertureLocation

logger = logging.getLogger(__name__)


def find_grating_orders_simple(data: np.ndarray, mask: Optional[np.ndarray] = None,
                             snr_threshold: float = 5.0,
                             step_denominator: int = 20,
                             gap_fill_factor: float = 1.35,
                             gap_fill_snr: float = 2.5,
                             min_trace_coverage: float = 0.20,
                             trace_degree: int = 4,
                             output_dir_base: str = '') -> ApertureSet:
    """
    Find the positions of grating orders on a CCD image.

    This is a simplified algorithm for grating spectrographs where orders
    run horizontally (left to right) across the detector.

    Args:
        data: Image data (2D array)
        mask: Optional bad pixel mask (same shape as data)
        snr_threshold: Detection threshold (multiples of sigma above profile baseline)
        gap_fill_factor: Factor for gap detection (gap > factor * predicted → missing).
        gap_fill_snr: SNR threshold for faint-order gap filling.
        min_trace_coverage: Minimum fraction of detector width a traced
            order must cover to be accepted.
        trace_degree: Polynomial degree for center tracing.
        output_dir_base: Output directory for diagnostic plots.

    Returns:
        ApertureSet with detected orders
    """
    if mask is None:
        mask = np.zeros_like(data, dtype=np.int32)

    h, w = data.shape

    logger.info(f"Finding grating orders in image of shape {data.shape}")

    img = np.asarray(data, dtype=np.float64)

    x_center = w // 2
    half_band = max(5, w // 20)
    x1 = max(0, x_center - half_band)
    x2 = min(w, x_center + half_band + 1)
    # Use LINEAR profile for seeding so the SNR threshold perfectly matches
    # the physical ADU SNR definition.
    center_profile = np.median(img[:, x1:x2], axis=1)
    center_profile = gaussian_filter1d(center_profile, sigma=1.5)

    # -------------------------------------------------------------------------
    # Auto-estimate dynamic order separation (Blind Peak Finding & Curve Fit)
    # -------------------------------------------------------------------------
    smooth_prof = gaussian_filter1d(center_profile, sigma=2.0)
    p99 = np.nanpercentile(smooth_prof, 99)
    rough_peaks, _ = find_peaks(smooth_prof, distance=5, prominence=max(5.0, p99 * 0.02))
    
    med_sep = 30.0
    def _fallback_sep(y): return med_sep
    local_sep_func = _fallback_sep

    if len(rough_peaks) >= 5:
        diffs = np.diff(rough_peaks)
        mids = 0.5 * (rough_peaks[:-1] + rough_peaks[1:])
        med_sep = float(np.median(diffs))

        # --- UnivariateSpline整体平滑 ---
        from scipy.interpolate import UnivariateSpline

        # 对间距做整体平滑
        s_factor = 0.5 * len(diffs)  # 平滑度参数，可根据实际调整
        try:
            spline = UnivariateSpline(mids, diffs, s=s_factor)
            def _dynamic_sep(y):
                return float(spline(y))
            local_sep_func = _dynamic_sep
            _inlier_min = float(np.min(diffs))
            _inlier_max = float(np.max(diffs))
            logger.info(f"Auto-estimated dynamic separation (UnivariateSpline): median={med_sep:.1f} px, s={s_factor}")
        except Exception as e:
            logger.warning(f"Spline fit failed: {e}, fallback to median separation.")
            local_sep_func = lambda y: med_sep
            _inlier_min = med_sep * 0.5
            _inlier_max = med_sep * 2.0
    else:
        logger.warning(f"Auto-estimation of separation failed, using fallback: {med_sep:.1f} px")

    # 1. Tracing controls (all explicit parameters now)
    fill_missing_orders = True

    # 2. Background Envelope Estimation
    # Use a large-scale rolling minimum to capture the background base.
    from scipy.ndimage import minimum_filter1d
    bg_size = int(2.5 * med_sep)
    envelope = minimum_filter1d(center_profile, size=bg_size)
    envelope = gaussian_filter1d(envelope, sigma=bg_size/4.0) # Smooth it
    
    # 3. Background-subtracted profile for seeding
    seed_pure = center_profile - envelope
    
    # Robust noise estimation on the residual
    def _estimate_prof_noise(prof):
        valid_prof = prof[np.isfinite(prof)]
        if len(valid_prof) < 10:
            return 1.0, 0.0
            
        # 1. Baseline estimation: bottom 25%
        sorted_prof = np.sort(valid_prof)
        bg_pixels = sorted_prof[:max(10, len(sorted_prof) // 4)]
        s_med = float(np.nanmedian(bg_pixels))
        
        # 2. Intrinsic noise estimation via MAD of differences
        diffs = np.diff(valid_prof)
        s_sig = float(np.nanmedian(np.abs(diffs)) / 0.9539) # sqrt(2)*0.6745
        
        return max(s_sig, 0.1), s_med

    noise, baseline = _estimate_prof_noise(seed_pure)
    seed_profile = seed_pure # Use the pure signal for peak finding
    
    prof_max = float(np.max(seed_profile))

    # Seeding threshold: single robust SNR-based gate
    det_threshold = baseline + snr_threshold * noise
    min_dist = max(3, int(med_sep * 0.20)) 

    peaks, _ = find_peaks(
        seed_profile,
        height=det_threshold,
        distance=min_dist,
        prominence=max(0.8 * noise, 0.15 * snr_threshold * noise),
    )

    logger.info(
        f"Found {len(peaks)} seed peaks (SNR > {snr_threshold:.1f}, "
        f"baseline={baseline:.4f}, noise={noise:.4f})"
    )
    
    if len(peaks) > 0:
        # Diagnostic: check the faintest detected peak and first peak region
        pk_vals = seed_profile[peaks]
        min_pk_idx = np.argmin(pk_vals)
        first_pk = int(peaks[0])
        
        logger.info(
            f"Faintest seed: y={peaks[min_pk_idx]}, val={pk_vals[min_pk_idx]:.4f}, "
            f"SNR={(pk_vals[min_pk_idx]-baseline)/(noise+1e-6):.2f}"
        )
        # Show what's below the first peak.
        check_y = max(0, first_pk - int(med_sep))
        for yy in range(check_y, first_pk):
            if seed_profile[yy] > baseline + 0.3 * noise:
                logger.info(
                    f"  Below-first: y={yy}, val={seed_profile[yy]:.4f}, "
                    f"SNR={(seed_profile[yy]-baseline)/(noise+1e-6):.2f}"
                )

    peaks = np.asarray(sorted(peaks), dtype=int)

    # -------------------------------------------------------------------------
    # 4. Fill missed weak orders in 1D cross-section
    # -------------------------------------------------------------------------
    if fill_missing_orders and len(peaks) >= 2:
        # User requirement: use 2.5 sigma for faint blind guessing.
        gap_fill_snr = max(2.5, gap_fill_snr)

        logger.info(
            f"Gap-fill analysis: {len(peaks)} seeds, global median_sep={med_sep:.1f}px, "
            f"gap_fill_factor=1.30, gap_fill_snr={gap_fill_snr:.2f}"
        )

        # Internal gap filling using curve
        if np.isfinite(med_sep) and med_sep > 2:
            inserted = []
            for i in range(len(peaks) - 1):
                left, right = peaks[i], peaks[i + 1]
                gap = right - left
                expected_sep = local_sep_func(0.5 * (left + right))
                ratio = gap / expected_sep
                
                # Gap > 1.3 curve prediction -> Missing orders
                if gap > 1.30 * expected_sep:
                    n_missing = int(round(gap / expected_sep)) - 1
                    if n_missing < 1: 
                        n_missing = 1
                    logger.info(
                        f"Gap-fill: gap between y={left} and y={right} "
                        f"is {gap}px ({ratio:.2f}x expected {expected_sep:.1f}), inserting {n_missing} candidate(s)"
                    )
                    for k in range(1, n_missing + 1):
                        guess = int(round(left + k * gap / (n_missing + 1)))
                        win = max(5, int(0.55 * expected_sep))
                        y1 = max(0, guess - win)
                        y2 = min(h, guess + win + 1)
                        if y2 - y1 < 3:
                            continue
                        seg = seed_profile[y1:y2]
                        seg_med = float(np.median(seg))
                        seg_mad = float(np.median(np.abs(seg - seg_med)))
                        seg_noise = 1.4826 * seg_mad if seg_mad > 0 else float(np.std(seg - seg_med))
                        
                        loc = int(np.argmax(seg))
                        cand = y1 + loc
                        cand_val = seg[loc]
                        
                        if seg_noise > 0 and cand_val < (seg_med + gap_fill_snr * seg_noise):
                            logger.debug(
                                f"  Gap-fill candidate y={cand}: rejected "
                                f"(SNR={(cand_val-seg_med)/(seg_noise+1e-6):.1f} < {gap_fill_snr})"
                            )
                            continue
                            
                        # Avoid duplicates
                        if np.all(np.abs(peaks - cand) > max(3, int(0.30 * expected_sep))):
                            inserted.append(cand)
                            logger.info(f"  Gap-fill candidate y={cand}: ACCEPTED (SNR={(cand_val-seg_med)/(seg_noise+1e-6):.1f})")

            # Edge extrapolation
            min_sep = max(3, int(0.30 * med_sep))

            left_ref = int(peaks[0])
            local_sep_bottom = local_sep_func(left_ref)
            edge_win_bottom = max(5, int(0.80 * local_sep_bottom))
            logger.info(
                f"Edge-fill (bottom) starting: left_ref=y{left_ref}, "
                f"local_sep={local_sep_bottom:.1f}, edge_win={edge_win_bottom}"
            )
            while left_ref - local_sep_func(left_ref) * 0.5 > 0:
                step_sep = local_sep_func(left_ref)
                edge_win = max(5, int(0.80 * step_sep))
                guess = int(round(left_ref - step_sep))
                if guess < 0:
                    guess = max(0, left_ref // 2)
                y1 = max(0, guess - edge_win)
                y2 = min(h, guess + edge_win + 1)
                if y2 - y1 < 3:
                    break
                seg = seed_profile[y1:y2]
                cand = y1 + int(np.argmax(seg))
                seg_med = float(np.median(seg))
                seg_mad = float(np.median(np.abs(seg - seg_med)))
                seg_noise = 1.4826 * seg_mad if seg_mad > 0 else float(np.std(seg - seg_med))
                cand_val = float(seed_profile[cand])
                threshold = seg_med + snr_threshold * seg_noise
                logger.info(
                    f"  Edge-fill (bottom): guess=y{guess}, window=[{y1},{y2}], "
                    f"argmax→y{cand}, val={cand_val:.4f}, "
                    f"seg_med={seg_med:.4f}, seg_noise={seg_noise:.4f}, "
                    f"threshold={threshold:.4f} ({snr_threshold:.1f}σ), left_ref=y{left_ref}"
                )
                if seg_noise > 0 and cand_val >= threshold:
                    if np.all(np.abs(peaks - cand) > min_sep):
                        inserted.append(cand)
                        logger.info(
                            f"Edge-fill (bottom): y={cand} ACCEPTED "
                            f"(val={cand_val:.3f} >= {threshold:.3f})"
                        )
                        left_ref = cand
                        continue
                logger.debug(
                    f"Edge-fill (bottom): y={cand} rejected "
                    f"(val={cand_val:.3f} < {threshold:.3f})"
                )
                break

            right_ref = int(peaks[-1])
            local_sep_top = local_sep_func(right_ref)
            edge_win_top = max(5, int(0.80 * local_sep_top))
            logger.info(
                f"Edge-fill (top) starting: right_ref=y{right_ref}, "
                f"local_sep={local_sep_top:.1f}, edge_win={edge_win_top}"
            )
            while right_ref + local_sep_func(right_ref) * 0.5 < h:
                step_sep = local_sep_func(right_ref)
                edge_win = max(5, int(0.80 * step_sep))
                guess = int(round(right_ref + step_sep))
                if guess >= h:
                    break
                y1 = max(0, guess - edge_win)
                y2 = min(h, guess + edge_win + 1)
                if y2 - y1 < 3:
                    break
                seg = seed_profile[y1:y2]
                cand = y1 + int(np.argmax(seg))
                seg_med = float(np.median(seg))
                seg_mad = float(np.median(np.abs(seg - seg_med)))
                seg_noise = 1.4826 * seg_mad if seg_mad > 0 else float(np.std(seg - seg_med))
                cand_val = float(seed_profile[cand])
                threshold = seg_med + snr_threshold * seg_noise
                if seg_noise > 0 and cand_val >= threshold:
                    if np.all(np.abs(peaks - cand) > min_sep):
                        inserted.append(cand)
                        logger.info(
                            f"Edge-fill (top): y={cand} ACCEPTED "
                            f"(val={cand_val:.3f} >= {threshold:.3f}, {snr_threshold:.1f}σ)"
                        )
                        right_ref = cand
                        continue
                logger.debug(
                    f"Edge-fill (top): y={cand} rejected "
                    f"(val={cand_val:.3f} < {threshold:.3f})"
                )
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
        peak_val = float(segment[loc])
        if not np.isfinite(peak_val) or peak_val <= 0:
            return None

        # Linear space intensity-weighted centroid
        bg = 0.5 * (float(segment[0]) + float(segment[-1]))
        sub_flux = segment.astype(float) - bg
        core_mask = sub_flux > 0.2 * (peak_val - bg)
        
        if np.sum(core_mask) >= 3:
            yy = np.arange(y1, y2, dtype=float)[core_mask]
            weights = sub_flux[core_mask]
            sum_w = np.sum(weights)
            if sum_w > 0:
                return float(np.sum(yy * weights) / sum_w)

        # Fallback
        return float(y1 + loc)

    def trace_order(seed_y: float) -> tuple:
        # Trace every column so the downstream sigma-clipping and Chebyshev fit
        # have dense coverage; faint-order wings use prediction + rejection to
        # avoid locking onto neighbouring orders.
        local_sep = local_sep_func(seed_y)
        search_half = max(4, int(local_sep * 0.45))

        # Refine the starting seed position to sub-pixel accuracy before tracing
        refined_y = refine_peak_y(img[:, x_center], seed_y, search_half)
        start_y = refined_y if refined_y is not None else float(seed_y)

        xs = [float(x_center)]
        ys = [start_y]

        # right
        y_prev = start_y
        miss_count = 0
        right_x = [float(x_center)]
        right_y = [start_y]
        for x in range(x_center + 1, w):
            if len(right_x) >= 4:
                recent_x = np.asarray(right_x[-4:], dtype=float)
                recent_y = np.asarray(right_y[-4:], dtype=float)
                coef1 = np.polyfit(recent_x, recent_y, deg=1)
                y_guess = float(np.polyval(coef1, x))
            else:
                y_guess = y_prev

            y_new = refine_peak_y(img[:, x], y_guess, search_half)
            if y_new is None:
                miss_count += 1
                if miss_count > 10:
                    break
                continue

            if abs(y_new - y_guess) > search_half:
                miss_count += 1
                if miss_count > 10:
                    break
                continue

            xs.append(x)
            ys.append(y_new)
            right_x.append(float(x))
            right_y.append(float(y_new))
            y_prev = y_new
            miss_count = 0

        # left
        y_prev = start_y
        miss_count = 0
        left_x = [float(x_center)]
        left_y = [start_y]
        for x in range(x_center - 1, -1, -1):
            if len(left_x) >= 4:
                recent_x = np.asarray(left_x[-4:], dtype=float)
                recent_y = np.asarray(left_y[-4:], dtype=float)
                coef1 = np.polyfit(recent_x, recent_y, deg=1)
                y_guess = float(np.polyval(coef1, x))
            else:
                y_guess = y_prev

            y_new = refine_peak_y(img[:, x], y_guess, search_half)
            if y_new is None:
                miss_count += 1
                if miss_count > 10:
                    break
                continue

            if abs(y_new - y_guess) > search_half:
                miss_count += 1
                if miss_count > 10:
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

    def estimate_order_width(xs: np.ndarray, ys: np.ndarray, local_sep: float) -> float:
        """Estimate aperture width from local FWHM samples along the traced order."""
        if xs.size == 0:
            return float(local_sep)

        sample_idx = np.linspace(0, xs.size - 1, num=min(24, xs.size), dtype=int)
        widths = []
        search_half = max(8, int(0.8 * local_sep))

        for i in sample_idx:
            x = int(round(xs[i]))
            yc = float(ys[i])
            y1 = max(0, int(np.floor(yc - search_half)))
            y2 = min(h, int(np.ceil(yc + search_half + 1)))
            if y2 - y1 < 7:
                continue

            prof = img[y1:y2, x]
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
            if 2.0 <= fwhm <= 2.5 * local_sep:
                widths.append(fwhm)

        if not widths:
            return float(local_sep)

        # Use aperture full width ~ 1.8*FWHM, clipped to sensible range.
        width = 1.8 * float(np.median(widths))
        width = float(np.clip(width, 6.0, 2.2 * local_sep))
        return width

    def _fit_center_to_coefs(x_vals: np.ndarray, y_vals: np.ndarray, 
                             degree: int, domain: tuple) -> np.ndarray:
        """Fit center using Chebyshev and return coefficients (ascending)."""
        x = np.asarray(x_vals, dtype=float)
        y = np.asarray(y_vals, dtype=float)
        
        if x.size < 4:
            return None

        deg = max(1, int(min(degree, max(1, x.size - 2))))
        cheb = Chebyshev.fit(x, y, deg=deg, domain=domain)
        return np.array(cheb.coef, dtype=float)
    def _fit_boundary_coeff(x_vals: np.ndarray, y_vals: np.ndarray,
                            degree: int, domain: tuple) -> np.ndarray:
        """Fit a smooth boundary using Chebyshev polynomials with robust outlier rejection."""
        x = np.asarray(x_vals, dtype=float)
        y = np.asarray(y_vals, dtype=float)
        if x.size < 4:
            return None
        
        deg = max(1, int(min(degree, max(1, x.size - 2))))
        try:
            # Iterative fit to handle outliers near detector edges
            cheb_obj = Chebyshev.fit(x, y, deg=deg, domain=domain)
            for _ in range(4):
                res = y - cheb_obj(x)
                s = np.std(res)
                if s <= 0.02: break
                good = np.abs(res) < 3.0 * s
                if np.sum(good) < deg + 2: break
                cheb_obj = Chebyshev.fit(x[good], y[good], deg=deg, domain=domain)
            return np.array(cheb_obj.coef, dtype=float)
        except:
            return None

    # 3. Tracing controls (all explicit parameters now)
    fill_missing_orders = True

    def _build_moffat_boundaries(cands: list,
                                   flux_fraction: float = 0.02,
                                   fwhm_scale: float = 1.5) -> None:
        """Determine order boundaries via two-step Moffat fitting.

        Step 1: For each sample column, use the traced center position as the
                peak.  Walk outward from that center in the log-compressed
                cross-section until the signal drops to *flux_fraction* × peak.
                This gives an initial lower/upper bracket.
        Step 2: Background = average of profile values at the two bracket
                positions.  Fit a Moffat profile to the background-subtracted
                cross-section within the bracket.
                Final boundary = fitted_center ± fwhm_scale × FWHM.
                Falls back to the Step-1 bracket when fitting fails.
        """
        from scipy.optimize import curve_fit

        def _moffat1d(y, y_c, I0, gamma, beta):
            return I0 * (1.0 + ((y - y_c) / gamma) ** 2) ** (-beta)

        if not cands:
            return

        x_samples = np.linspace(0, w - 1, num=min(128, w), dtype=int)

        for cand in cands:
            c_coef = cand['center_coef']
            c_model = Chebyshev(c_coef, domain=(0, w - 1))
            local_sep = local_sep_func(cand['y_center'])
            is_reg = cand.get('is_regularized', False)

            xl, yl, xu, yu = [], [], [], []
            xc, yc = [], []
            half_win = max(8, int(0.6 * local_sep))

            for x in x_samples:
                yc_float = float(c_model(x))          # traced center = peak
                y0 = int(round(yc_float))
                y1 = max(0, y0 - half_win)
                y2 = min(h, y0 + half_win + 1)
                if y2 - y1 < 5:
                    continue

                prof = img[y1:y2, x].astype(float)
                n = len(prof)
                mid_idx = max(0, min(n - 1, y0 - y1))  # center in local coords
                peak_val = float(prof[mid_idx])

                # ── Step 1: bracket at flux_fraction × peak ───────────────────
                threshold = flux_fraction * peak_val

                low_idx = 0
                for j in range(mid_idx, -1, -1):
                    if prof[j] <= threshold:
                        low_idx = j
                        break

                up_idx = n - 1
                for j in range(mid_idx, n):
                    if prof[j] <= threshold:
                        up_idx = j
                        break

                # ── Step 2: Moffat fit within bracket ─────────────────────────
                bg = 0.5 * (float(prof[low_idx]) + float(prof[up_idx]))
                seg_y = np.arange(low_idx, up_idx + 1, dtype=float)
                seg_f = np.maximum(prof[low_idx:up_idx + 1] - bg, 0.0)

                fwhm_final = None
                yc_fit = float(mid_idx)

                try:
                    if len(seg_y) >= 5 and (peak_val - bg) > 0:
                        I0_guess = max(float(peak_val - bg), 1e-6)
                        gamma_guess = max(1.0, (up_idx - low_idx) / 4.0)
                        p0 = [float(mid_idx), I0_guess, gamma_guess, 2.5]
                        b_lo = [float(low_idx), 0.0, 0.3, 0.5]
                        b_hi = [float(up_idx), I0_guess * 5.0,
                                float(up_idx - low_idx), 10.0]
                        popt, _ = curve_fit(
                            _moffat1d, seg_y, seg_f,
                            p0=p0, bounds=(b_lo, b_hi),
                            maxfev=3000
                        )
                        yc_fit, _I0f, gamma_fit, beta_fit = popt
                        if beta_fit > 0.5 and gamma_fit > 0.3:
                            fwhm_final = 2.0 * gamma_fit * np.sqrt(
                                2.0 ** (1.0 / beta_fit) - 1.0)
                except Exception:
                    pass

                if fwhm_final is not None and 1.0 <= fwhm_final <= local_sep:
                    # Moffat fit is good, trust its center even for faint orders
                    y_center_actual = float(y1 + yc_fit)
                    half_bnd = fwhm_scale * fwhm_final
                    y_lo = y_center_actual - half_bnd
                    y_hi = y_center_actual + half_bnd
                else:
                    # Moffat fit failed or is unreliable, stick to the polynomial trace
                    y_center_actual = yc_float
                    # Fallback to Step-1 bracket, forced symmetric for faint orders
                    if is_reg:
                        half_bnd = (up_idx - low_idx) / 2.0
                        y_lo = y_center_actual - half_bnd
                        y_hi = y_center_actual + half_bnd
                    else:
                        y_lo = float(y1 + low_idx)
                        y_hi = float(y1 + up_idx)

                xl.append(float(x))
                yl.append(y_lo)
                xu.append(float(x))
                yu.append(y_hi)
                # Always collect the best available center
                xc.append(float(x))
                yc.append(y_center_actual)

            # Upgrade the center trace using the new, more accurate centers
            if len(xc) >= 5:
                new_c_coef = _fit_boundary_coeff(np.array(xc), np.array(yc), trace_degree, (0, w - 1))
                if new_c_coef is not None:
                    cand['center_coef'] = new_c_coef

            cand['lower_coef'] = _fit_boundary_coeff(
                np.array(xl), np.array(yl), trace_degree, (0, w - 1))
            cand['upper_coef'] = _fit_boundary_coeff(
                np.array(xu), np.array(yu), trace_degree, (0, w - 1))

    def _plot_seeding_diagnostics(out_path):
        """Plot the seeding profile with envelope and thresholds restored to original state."""
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(12, 6))
        y_idx = np.arange(len(center_profile))
        
        ax.plot(y_idx, center_profile, color='gray', alpha=0.5, label='Raw Profile (log)')
        ax.plot(y_idx, envelope, color='blue', linestyle='--', alpha=0.8, label='BG Envelope')
        
        det_line = envelope + snr_threshold * noise
        root_line = envelope + 3.0 * noise
        
        ax.plot(y_idx, det_line, color='orange', label=f'Detection ({snr_threshold}σ)')
        ax.plot(y_idx, root_line, color='green', alpha=0.6, label='Aperture Boundary (3σ)')
        
        # Highlight detected peaks
        if len(peaks) > 0:
            ax.scatter(peaks, center_profile[peaks], color='red', s=15, zorder=5, label='Seeds')
            
        ax.set_title("Order Seeding & Threshold Diagnostics (Restored with Envelope)")
        ax.set_xlabel("Detector Y (px)")
        ax.set_ylabel("Intensity (log units)")
        ax.legend(loc='upper right', fontsize='small')
        ax.grid(True, alpha=0.2)
        
        plt.tight_layout()
        plt.savefig(out_path, dpi=150)
        plt.close()
        logger.info(f"Saved seeding diagnostics plot to {out_path}")

 # Helpers are defined, execution will happen at the end of the function.

    aperture_set = ApertureSet()
    candidates = []

    # 读取谱级间距容差参数
    import configparser
    spacing_tol = 0.3
    try:
        from src.config.config_manager import ConfigManager
        cfg = ConfigManager()
        spacing_tol = cfg.get_float('reduce.trace', 'spacing_tol', 0.3)
    except Exception:
        pass

    for order_idx, seed in enumerate(peaks):
        local_sep = local_sep_func(float(seed))
        x_positions, y_positions = trace_order(float(seed))

        if len(x_positions) < 15:
            logger.warning(f"Seed {order_idx + 1} (y~{seed:.1f}): Rejected - Not enough traced points ({len(x_positions)} < 15)")
            continue

        try:
            # 1. Iterative outlier rejection using一个简单多项式
            good = np.ones_like(x_positions, dtype=bool)
            deg_rough = 3
            for _ in range(3):
                if np.sum(good) < deg_rough + 2: break
                coef = np.polyfit(x_positions[good], y_positions[good], deg=deg_rough)
                resid = y_positions - np.polyval(coef, x_positions)
                std = np.std(resid[good])
                if std <= 0.05: break
                good = np.abs(resid) < 3.0 * std

            if np.sum(good) < 15:
                logger.warning(f"Seed {order_idx + 1} (y~{seed:.1f}): Rejected - Too few valid points after clipping ({np.sum(good)} < 15)")
                continue

            # 2. Final Center fitting (Always Chebyshev)
            center_deg = trace_degree
            center_coef = _fit_center_to_coefs(x_positions[good], y_positions[good], center_deg, (0, float(w - 1)))
            if center_coef is None:
                logger.warning(f"Seed {order_idx + 1} (y~{seed:.1f}): Rejected - Chebyshev center fit failed")
                continue

            center_model = Chebyshev(center_coef, domain=(0, w - 1))
            yc_mid = center_model(x_center)

            x_coverage = float(np.max(x_positions[good]) - np.min(x_positions[good]))
            min_cov_px = float(min_trace_coverage) * w
            if x_coverage < min_cov_px:
                logger.warning(f"Seed {order_idx + 1} (y~{seed:.1f}): Rejected - Trace coverage too short ({x_coverage:.0f}px < min {min_cov_px:.0f}px)")
                continue


            width = estimate_order_width(x_positions[good], y_positions[good], local_sep)
            half_width = width / 2.0
            lower_coef = np.array([yc_mid - half_width], dtype=float)
            upper_coef = np.array([yc_mid + half_width], dtype=float)

            candidates.append({
                'center_coef': center_coef,
                'lower_coef': lower_coef,
                'upper_coef': upper_coef,
                'width': float(width),
                'n_good': int(np.sum(good)),
                'fit_rms': float(np.std(y_positions[good] - center_model(x_positions[good]))),
                'x_coverage': float(x_coverage),
                'y_center': float(yc_mid),
                'seed_y': float(seed),
            })

            logger.info(
                f"Candidate from seed {order_idx + 1}: points={np.sum(good)}, "
                f"fit=chebyshev, width={width:.1f}px"
            )

        except Exception as e:
            logger.warning(f"Order {order_idx + 1}: Failed polynomial tracing: {e}")
            continue

    # Merge duplicate candidates that represent the same physical order,
    # then assign continuous order IDs.
    # Ensure boundary_frac and fwhm_scale are available in this scope (from function arguments)
    boundary_frac = boundary_frac if 'boundary_frac' in locals() else 0.02
    fwhm_scale = fwhm_scale if 'fwhm_scale' in locals() else 1.5
    if candidates:
        candidates_sorted = sorted(candidates, key=lambda c: c['y_center'])
        unique = []

        # 恢复并增强谱级筛选判据
        if candidates:
            candidates_sorted = sorted(candidates, key=lambda c: c['y_center'])
            unique = []

            for idx, cand in enumerate(candidates_sorted):
                if not unique:
                    unique.append(cand)
                    continue

                prev = unique[-1]
                pred_sep = local_sep_func(0.5 * (cand['y_center'] + prev['y_center']))
                actual_sep = abs(cand['y_center'] - prev['y_center'])
                # 输出调试信息
                logger.info(f"Seed spacing debug: idx={idx}, y1={prev['y_center']:.2f}, y2={cand['y_center']:.2f}, actual={actual_sep:.2f}, predicted={pred_sep:.2f}, ratio={actual_sep/pred_sep:.2f}")
                # 判据：实际间距小于0.3×预测间距则合并（保留分数高的），大于10.0×预测间距才排除
                if actual_sep < 0.3 * pred_sep:
                    prev_score = prev['n_good'] - 2.0 * prev['fit_rms'] + 0.002 * prev['x_coverage']
                    cand_score = cand['n_good'] - 2.0 * cand['fit_rms'] + 0.002 * cand['x_coverage']
                    if cand_score > prev_score:
                        unique[-1] = cand
                elif actual_sep > 10.0 * pred_sep:
                    logger.warning(f"Seed at y={cand['y_center']:.1f} rejected: spacing {actual_sep:.2f} > 10.0×predicted {pred_sep:.2f}")
                    continue
                else:
                    unique.append(cand)

        # ---------------------------------------------------------------------
        # Global Shape Regularization (Prevents faint orders from crossing)
        # ---------------------------------------------------------------------
        if len(unique) >= 3:
            def is_robust_shape(c):
                # Only use extremely well-traced orders across most of the CCD as global shape anchors.
                # Loose RMS allows wiggly noise-chasing traces to corrupt the global shape.
                return (c['x_coverage'] > 0.85 * w) and (c['fit_rms'] < max(1.5, 0.08 * local_sep_func(c['y_center'])))

            robust = [c for c in unique if is_robust_shape(c)]
            if len(robust) >= 3:
                max_len = max(len(c['center_coef']) for c in robust)
                # Ensure robust is strictly sorted by y_center so np.interp works correctly
                robust_sorted = sorted(robust, key=lambda c: c['y_center'])
                y_c_robust = np.array([c['y_center'] for c in robust_sorted])
                coefs_robust = np.array([np.pad(c['center_coef'], (0, max_len - len(c['center_coef']))) for c in robust_sorted])
                
                for c in unique:
                    # 如果该级次本身就很亮且轨迹稳健，则保留其自身真实的光学畸变轨迹
                    if is_robust_shape(c):
                        c['is_regularized'] = False
                        continue

                    c_model_old = Chebyshev(c['center_coef'], domain=(0, w - 1))
                    # Anchor directly to the robust photometric seed peak found in the median profile
                    y_anchor = c['y_center']
                    
                    new_coef = np.zeros(max_len)
                    for d in range(1, max_len):
                        # np.interp guarantees flat extrapolation at edges:
                        # Faint edge orders will exactly inherit the shape of the outermost bright order.
                        new_coef[d] = np.interp(y_anchor, y_c_robust, coefs_robust[:, d])
                        
                    c_model_new_temp = Chebyshev(new_coef, domain=(0, w - 1))
                    y_mid_new_temp = c_model_new_temp(x_center)
                    
                    # T_0(x) is 1 everywhere, so adjusting c[0] shifts the whole curve exactly
                    new_coef[0] = y_anchor - y_mid_new_temp
                    c['center_coef'] = new_coef
                    c['y_center'] = y_anchor
                    c['is_regularized'] = True
                    
                logger.info(f"Applied global shape regularization to {len(unique) - len(robust)} faint orders based on {len(robust)} robust traces.")

        # Apply Moffat boundary detection to the unique orders
        _build_moffat_boundaries(unique, flux_fraction=boundary_frac, fwhm_scale=fwhm_scale)

        # ---------------------------------------------------------------------
        # 2D Coefficient Interpolation for Completely Missing Orders
        # ---------------------------------------------------------------------
        if len(unique) >= 2:
            final_cands = []
            final_cands.append(unique[0])
            
            for i in range(len(unique) - 1):
                left_cand = unique[i]
                right_cand = unique[i + 1]
                gap = right_cand['y_center'] - left_cand['y_center']
                
                expected_sep = local_sep_func(0.5 * (left_cand['y_center'] + right_cand['y_center']))
                
                if gap > 1.30 * expected_sep:
                    n_missing = int(round(gap / expected_sep)) - 1
                    if n_missing > 0:
                        logger.info(
                            f"2D Interpolation: gap between y={left_cand['y_center']:.1f} "
                            f"and y={right_cand['y_center']:.1f} is {gap:.1f}px "
                            f"(expected {expected_sep:.1f}), synthesising {n_missing} completely blind order(s)."
                        )
                        for k in range(1, n_missing + 1):
                            alpha = float(k) / (n_missing + 1.0)
                            interp_cand = {}
                            # Pad to same length if degrees somehow differ (though they should be standard center_deg)
                            len_c = max(len(left_cand['center_coef']), len(right_cand['center_coef']))
                            lc_pad = np.pad(left_cand['center_coef'], (0, len_c - len(left_cand['center_coef'])))
                            rc_pad = np.pad(right_cand['center_coef'], (0, len_c - len(right_cand['center_coef'])))
                            interp_cand['center_coef'] = (1.0 - alpha) * lc_pad + alpha * rc_pad
                            
                            len_l = max(len(left_cand['lower_coef']), len(right_cand['lower_coef']))
                            ll_pad = np.pad(left_cand['lower_coef'], (0, len_l - len(left_cand['lower_coef'])))
                            rl_pad = np.pad(right_cand['lower_coef'], (0, len_l - len(right_cand['lower_coef'])))
                            interp_cand['lower_coef'] = (1.0 - alpha) * ll_pad + alpha * rl_pad
                            
                            len_u = max(len(left_cand['upper_coef']), len(right_cand['upper_coef']))
                            lu_pad = np.pad(left_cand['upper_coef'], (0, len_u - len(left_cand['upper_coef'])))
                            ru_pad = np.pad(right_cand['upper_coef'], (0, len_u - len(right_cand['upper_coef'])))
                            interp_cand['upper_coef'] = (1.0 - alpha) * lu_pad + alpha * ru_pad
                            
                            interp_cand['width'] = (1.0 - alpha) * left_cand['width'] + alpha * right_cand['width']
                            interp_cand['y_center'] = (1.0 - alpha) * left_cand['y_center'] + alpha * right_cand['y_center']
                            interp_cand['n_good'] = 0
                            interp_cand['fit_rms'] = 0.0
                            interp_cand['x_coverage'] = float(w)
                            interp_cand['is_interpolated'] = True
                            final_cands.append(interp_cand)
                            
                final_cands.append(right_cand)
            
            unique = final_cands

        # Save diagnostic plot (showing raw profile + envelope + thresholds)
        # NOTE: seeding diagnostics are now merged into apertures_q*.png via
        # _plot_apertures_cross_section(); do not regenerate a separate file here.

        for new_id, cand in enumerate(unique, start=1):
            aperture_loc = ApertureLocation(
                aperture=new_id,
                order=new_id,
                center_coef=cand.get('center_coef'),
                width=cand['width'],
                is_chebyshev=True, # native chebyshev for both center and boundaries
                domain=(0, float(w - 1)),
                is_interpolated=cand.get('is_interpolated', False)
            )
            aperture_set.add_aperture(aperture_loc)

    logger.info(f"Detected {aperture_set.norders} orders (3-sigma boundary mode)")
    return aperture_set


# ---------------------------------------------------------------------------
# Order index assignment & gap-aware interpolation
# ---------------------------------------------------------------------------


def assign_order_indices(apertures: ApertureSet, gap_fill_factor: float = 1.35) -> dict:
    """Assign integer order indices *m* to detected apertures.

    Walk ordered apertures from bottom (small Y) to top (large Y) and increment *m*
    by 1 for each normal step.  When the gap between two consecutive centres
    exceeds *gap_fill_factor × local_expected_sep*, skip the appropriate number
    of *m* values so that the missing orders have reserved slots.

    The expected separation is modelled as a smooth polynomial of sequential
    index so that it adapts to the monotonically increasing spacing of
    echelle orders (short-wavelength orders are wider apart).

    Returns
    -------
    mapping : dict[int, int]
        {aperture_id: m_index}
    """
    order_ids = apertures.get_orders()
    if len(order_ids) < 2:
        return {oid: i for i, oid in enumerate(order_ids)}

    # Evaluate centre Y at image midpoint (X = 2048 as typical reference).
    x_ref = 2048.0
    centres = {}
    for oid in order_ids:
        ap = apertures.get_aperture(oid)
        centres[oid] = float(ap.get_position(np.array([x_ref]))[0])

    sorted_ids = sorted(order_ids, key=lambda oid: centres[oid])
    diffs = np.array([centres[sorted_ids[i + 1]] - centres[sorted_ids[i]]
                       for i in range(len(sorted_ids) - 1)])
    med_sep = float(np.median(diffs)) if len(diffs) else 30.0

    # Build a smooth LOCAL expected-separation model.
    # For echelle spectra the separation increases monotonically; a low-order
    # polynomial captures this trend so that the gap detector doesn't treat
    # the naturally wider upper-order spacing as missing orders.
    #
    # IMPORTANT: the polynomial fit must be ROBUST to gap outliers.
    # Otherwise a large gap (e.g. 191 px among ~60 px diffs) skews the
    # polynomial and inflates the local expected separation, causing the
    # gap to be under-counted.  We use iterative sigma-clipping: fit the
    # polynomial, reject outliers, refit on inliers only, then evaluate
    # the clean polynomial at all positions.
    n = len(diffs)
    if n >= 6:
        idx_arr = np.arange(n, dtype=float)
        fit_deg = min(3, n - 1)
        good = np.ones(n, dtype=bool)
        for _clip_iter in range(4):
            coef_tmp = np.polyfit(idx_arr[good], diffs[good], fit_deg)
            model_tmp = np.polyval(coef_tmp, idx_arr)
            resid = diffs - model_tmp
            std_r = np.std(resid[good])
            if std_r <= 0:
                break
            new_good = np.abs(resid) < 2.5 * std_r
            if np.sum(new_good) < max(fit_deg + 2, 6):
                break  # too few inliers — keep previous good set
            if np.array_equal(new_good, good):
                break
            good = new_good
        # Final fit on clean data, evaluate everywhere.
        sep_coef = np.polyfit(idx_arr[good], diffs[good], fit_deg)
        expected_seps = np.polyval(sep_coef, idx_arr)
        # Clamp to reasonable range (based on INLIER diffs only).
        inlier_min = np.min(diffs[good])
        inlier_max = np.max(diffs[good])
        expected_seps = np.clip(expected_seps, inlier_min * 0.5,
                                inlier_max * 2.0)
    else:
        expected_seps = np.full(n, med_sep)

    mapping = {}
    m = 0
    mapping[sorted_ids[0]] = m
    for i in range(1, len(sorted_ids)):
        gap = centres[sorted_ids[i]] - centres[sorted_ids[i - 1]]
        local_sep = float(expected_seps[i - 1])
        ratio = gap / max(local_sep, 1.0)
        if ratio >= gap_fill_factor:
            n_skip = max(2, int(round(ratio)))
        else:
            n_skip = 1
        m += n_skip
        mapping[sorted_ids[i]] = m

    logger.info(f"Order index assignment: {len(mapping)} apertures, m_max={m}, "
                f"med_sep={med_sep:.1f}px, gap_fill_factor={gap_fill_factor:.2f}")
    # Log any gaps.
    prev_m = -1
    for oid in sorted_ids:
        cur_m = mapping[oid]
        if prev_m >= 0 and cur_m - prev_m > 1:
            n_miss = cur_m - prev_m - 1
            logger.info(f"  Gap detected: m={prev_m}→{cur_m} ({n_miss} missing order(s))")
        prev_m = cur_m

    return mapping


def renumber_apertures_with_gaps(apertures: ApertureSet, gap_fill_factor: float = 1.35) -> ApertureSet:
    """Renumber aperture IDs to reflect physical order indices with gaps.

    After initial detection, apertures have sequential IDs (1, 2, 3, ...).
    This function re-assigns IDs so that missing orders cause ID gaps
    (e.g. 1, 2, ..., 49, 51, ..., 100 if order 50 is missing).
    New ID = m + 1 where m is the gap-aware index from assign_order_indices().
    """
    order_map = assign_order_indices(apertures, gap_fill_factor)

    new_set = ApertureSet()
    for old_id, m in order_map.items():
        ap = apertures.get_aperture(old_id)
        new_id = m + 1
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
        new_set.add_aperture(new_ap)

    # Log renumbered IDs if any gaps exist.
    old_ids = sorted(apertures.get_orders())
    new_ids = sorted(new_set.get_orders())
    if old_ids != new_ids:
        logger.info(f"Renumbered apertures with gap-aware IDs: "
                    f"max_id {old_ids[-1]}→{new_ids[-1]} ({len(new_ids)} orders)")
    return new_set




def _refine_predicted_trace(flat_data: np.ndarray, pred_center: np.ndarray,
                            image_width: int, local_sep: float,
                            trace_degree: int = 3,
                            min_coverage: float = 0.15,
                            step_denominator: int = 20) -> Optional[np.ndarray]:
    """Refine a surface-predicted centre trace using actual flat image peaks.

    For a grid of x positions across the image, searches for the nearest peak
    to the predicted center and performs Gaussian sub-pixel refinement.  If
    enough valid peaks are found, fits a polynomial to produce a refined trace.

    Returns refined centre array (length = image_width) or None if refinement
    failed (too few valid peaks).
    """
    h, w = flat_data.shape
    search_half = max(6, int(0.4 * local_sep))
    step = max(1, step_denominator)  # use step_denominator as step size

    xs, ys = [], []
    for x in range(0, w, step):
        y_guess = pred_center[x]
        y0 = int(round(y_guess))
        y1 = max(0, y0 - search_half)
        y2 = min(h, y0 + search_half + 1)
        if y2 - y1 < 5:
            continue

        col = flat_data[y1:y2, x].astype(float)
        loc = int(np.argmax(col))
        peak_val = col[loc]
        if not np.isfinite(peak_val) or peak_val <= 0:
            continue

        # Linear space intensity-weighted centroid
        hw = 3
        g1 = max(0, loc - hw)
        g2 = min(col.size, loc + hw + 1)
        sub = col[g1:g2].astype(float)
        
        bg = 0.5 * (float(sub[0]) + float(sub[-1]))
        sub_flux = sub - bg
        pos = sub_flux > 0.2 * (peak_val - bg)

        refined_y = None
        if np.sum(pos) >= 3:
            yy = np.arange(g1, g2, dtype=float)[pos]
            weights = sub_flux[pos]
            sum_w = np.sum(weights)
            if sum_w > 0:
                refined_y = float(np.sum(yy * weights) / sum_w)
                
        if refined_y is None:
            refined_y = float(y1 + loc)

        # Reject if too far from prediction (likely jumped to a neighbour).
        if abs(refined_y - y_guess) > search_half:
            continue

        xs.append(float(x))
        ys.append(refined_y)

    if len(xs) < max(5, int(min_coverage * (w / step))):
        return None

    xs_arr = np.asarray(xs)
    ys_arr = np.asarray(ys)

    # Sigma-clip outliers (2.5σ, 3 iterations) before polynomial fit.
    deg = min(trace_degree, max(1, len(xs_arr) - 2))
    for _ in range(3):
        if len(xs_arr) < deg + 2:
            break
        coef = np.polyfit(xs_arr, ys_arr, deg)
        resid = ys_arr - np.polyval(coef, xs_arr)
        sigma = max(0.5, float(np.std(resid)))
        keep = np.abs(resid) < 2.5 * sigma
        if np.all(keep):
            break
        xs_arr = xs_arr[keep]
        ys_arr = ys_arr[keep]

    if len(xs_arr) < deg + 2:
        return None

    coef = np.polyfit(xs_arr, ys_arr, deg)
    x_all = np.arange(image_width, dtype=float)
    return np.polyval(coef, x_all)


def fill_missing_orders_by_interpolation(
        apertures: ApertureSet,
        image_width: int, image_height: int,
        trace_degree: int = 3,
        n_extend_below: int = 0,
        n_extend_above: int = 0,
        flat_data: Optional[np.ndarray] = None,
        gap_fill_factor: float = 1.35,
        step_denominator: int = 20) -> ApertureSet:
    """Fill interior gaps by interpolating between neighbouring real orders.

    Uses :func:`assign_order_indices` to detect gaps, then for each gap
    linearly interpolates centre/lower/upper traces between the bounding
    orders (per-column).  When *flat_data* is provided, interpolated centres
    are refined by tracing actual peaks in the flat image.

    Edge extension (``n_extend_below`` / ``n_extend_above``) is handled by
    :func:`_extrapolate_virtual_orders` (spacing-based, not surface-based).

    Parameters
    ----------
    apertures : original ApertureSet from tracing (with gap-aware IDs).
    image_width, image_height : detector dimensions.
    trace_degree : polynomial degree for per-order coefficients.
    n_extend_below, n_extend_above : extra orders to add beyond edges.
    flat_data : optional 2D flat image for image-based refinement.
    gap_fill_factor : forwarded to assign_order_indices for gap detection.

    Returns
    -------
    filled_apertures : new ApertureSet with originals + interpolated orders.
    """
    order_map = assign_order_indices(apertures, gap_fill_factor)
    if not order_map:
        return apertures

    m_to_oid = {m: oid for oid, m in order_map.items()}
    m_min = min(order_map.values())
    m_max = max(order_map.values())

    # Sort real orders by m-index for neighbour lookup.
    sorted_m = sorted(m_to_oid.keys())

    filled = ApertureSet()
    x_all = np.arange(image_width, dtype=float)
    n_filled = 0
    n_refined = 0

    # --- Fill interior gaps by interpolation ---
    for i in range(len(sorted_m) - 1):
        m_lo = sorted_m[i]
        m_hi = sorted_m[i + 1]
        # Copy the lower-bound real order.
        filled.add_aperture(apertures.get_aperture(m_to_oid[m_lo]))

        n_gap = m_hi - m_lo - 1
        if n_gap <= 0:
            continue

        # Get traces of bounding orders at every column.
        ap_lo = apertures.get_aperture(m_to_oid[m_lo])
        ap_hi = apertures.get_aperture(m_to_oid[m_hi])
        cen_lo = ap_lo.get_position(x_all)
        cen_hi = ap_hi.get_position(x_all)
        lo_lo = ap_lo.get_lower(x_all)
        lo_hi = ap_hi.get_lower(x_all)
        up_lo = ap_lo.get_upper(x_all)
        up_hi = ap_hi.get_upper(x_all)

        for k in range(1, n_gap + 1):
            frac = k / (n_gap + 1)
            pred_cen = cen_lo + frac * (cen_hi - cen_lo)
            pred_lower = lo_lo + frac * (lo_hi - lo_lo)
            pred_upper = up_lo + frac * (up_hi - up_lo)

            yc_mid = pred_cen[image_width // 2]
            if yc_mid < 0 or yc_mid >= image_height:
                continue

            # Try image-based refinement if flat data available.
            refined_center = None
            if flat_data is not None:
                local_sep = float(np.median(cen_hi - cen_lo)) / (n_gap + 1)
                refined_center = _refine_predicted_trace(
                    flat_data, pred_cen, image_width, local_sep,
                    trace_degree=trace_degree,
                    step_denominator=step_denominator)

            if refined_center is not None:
                half_w_up = pred_upper - pred_cen
                half_w_low = pred_cen - pred_lower
                cen_coef = np.polyfit(x_all, refined_center, trace_degree)
                low_coef = np.polyfit(x_all, refined_center - half_w_low, trace_degree)
                up_coef = np.polyfit(x_all, refined_center + half_w_up, trace_degree)
                width = float(np.median(half_w_up + half_w_low))
                shift = float(np.median(refined_center - pred_cen))
                logger.info(f"  Refined interpolated m={m_lo + k}: shifted {shift:+.1f}px")
                n_refined += 1
            else:
                cen_coef = np.polyfit(x_all, pred_cen, trace_degree)
                low_coef = np.polyfit(x_all, pred_lower, trace_degree)
                up_coef = np.polyfit(x_all, pred_upper, trace_degree)
                width = float(np.median(pred_upper - pred_lower))

            new_id = m_lo + k + 1  # consistent with renumber_apertures_with_gaps
            ap = ApertureLocation(
                aperture=new_id, order=new_id,
                center_coef=cen_coef, lower_coef=low_coef,
                upper_coef=up_coef, width=max(4.0, width),
            )
            filled.add_aperture(ap)
            logger.info(f"  Filled gap order m={m_lo + k}: id={new_id}, "
                        f"y_center≈{yc_mid:.1f}, width≈{width:.1f}px")
            n_filled += 1

    # Copy the last real order.
    filled.add_aperture(apertures.get_aperture(m_to_oid[sorted_m[-1]]))

    # --- Edge extension via spacing extrapolation ---
    if n_extend_below > 0:
        virts = _extrapolate_virtual_orders(
            filled, image_width, image_height,
            n_orders=n_extend_below, side='below',
            trace_degree=trace_degree)
        for ap in virts:
            # Use positive IDs below the current minimum.
            min_id = min(filled.get_orders())
            new_id = min_id - 1
            if new_id <= 0:
                new_id = min_id - 1
            ap.aperture = new_id
            ap.order = new_id
            filled.add_aperture(ap)
            n_filled += 1

    if n_extend_above > 0:
        virts = _extrapolate_virtual_orders(
            filled, image_width, image_height,
            n_orders=n_extend_above, side='above',
            trace_degree=trace_degree)
        for ap in virts:
            max_id = max(filled.get_orders())
            new_id = max_id + 1
            ap.aperture = new_id
            ap.order = new_id
            filled.add_aperture(ap)
            n_filled += 1

    logger.info(f"Interpolation gap-fill: {n_filled} order(s) added, total {filled.norders} orders")
    return filled


# ======================================================================
# MERGED ALGORITHM KERNEL FROM grating_trace_simple.py
# ======================================================================


"""
Simple grating order tracing for spectral reduction pipeline.

Handles order detection and tracing for grating spectrographs,
where orders run horizontally (left to right) across the detector.
"""

from scipy.signal import find_peaks, peak_widths


def _find_grating_orders_impl(data: np.ndarray, mask: Optional[np.ndarray] = None,
                             snr_threshold: float = 5.0,
                             step_denominator: int = 20,
                             gap_fill_factor: float = 1.35,
                             gap_fill_snr: float = 2.5,
                             min_trace_coverage: float = 0.20,
                             trace_degree: int = 4,
                             output_dir_base: str = '') -> ApertureSet:
    """
    Find the positions of grating orders on a CCD image.

    This is a simplified algorithm for grating spectrographs where orders
    run horizontally (left to right) across the detector.

    Args:
        data: Image data (2D array)
        mask: Optional bad pixel mask (same shape as data)
        snr_threshold: Detection threshold (multiples of sigma above profile baseline)
        gap_fill_factor: Factor for gap detection (gap > factor * predicted → missing).
        gap_fill_snr: SNR threshold for faint-order gap filling.
        min_trace_coverage: Minimum fraction of detector width a traced
            order must cover to be accepted.
        trace_degree: Polynomial degree for center tracing.
        output_dir_base: Output directory for diagnostic plots.

    Returns:
        ApertureSet with detected orders
    """
    if mask is None:
        mask = np.zeros_like(data, dtype=np.int32)

    h, w = data.shape

    logger.info(f"Finding grating orders in image of shape {data.shape}")

    img = np.asarray(data, dtype=np.float64)

    x_center = w // 2
    half_band = max(5, w // 20)
    x1 = max(0, x_center - half_band)
    x2 = min(w, x_center + half_band + 1)
    # Use LINEAR profile for seeding so the SNR threshold perfectly matches
    # the physical ADU SNR definition.
    center_profile = np.median(img[:, x1:x2], axis=1)
    center_profile = gaussian_filter1d(center_profile, sigma=1.5)

    # -------------------------------------------------------------------------
    # Auto-estimate dynamic order separation (Blind Peak Finding & Curve Fit)
    # -------------------------------------------------------------------------
    smooth_prof = gaussian_filter1d(center_profile, sigma=2.0)
    p99 = np.nanpercentile(smooth_prof, 99)
    rough_peaks, _ = find_peaks(smooth_prof, distance=5, prominence=max(5.0, p99 * 0.02))
    
    med_sep = 30.0
    def _fallback_sep(y): return med_sep
    local_sep_func = _fallback_sep

    if len(rough_peaks) >= 5:
        diffs = np.diff(rough_peaks)
        mids = 0.5 * (rough_peaks[:-1] + rough_peaks[1:])
        med_sep = float(np.median(diffs))

        # --- UnivariateSpline整体平滑 ---
        from scipy.interpolate import UnivariateSpline

        # 对间距做整体平滑
        s_factor = 0.5 * len(diffs)  # 平滑度参数，可根据实际调整
        try:
            spline = UnivariateSpline(mids, diffs, s=s_factor)
            def _dynamic_sep(y):
                return float(spline(y))
            local_sep_func = _dynamic_sep
            _inlier_min = float(np.min(diffs))
            _inlier_max = float(np.max(diffs))
            logger.info(f"Auto-estimated dynamic separation (UnivariateSpline): median={med_sep:.1f} px, s={s_factor}")
        except Exception as e:
            logger.warning(f"Spline fit failed: {e}, fallback to median separation.")
            local_sep_func = lambda y: med_sep
            _inlier_min = med_sep * 0.5
            _inlier_max = med_sep * 2.0
    else:
        logger.warning(f"Auto-estimation of separation failed, using fallback: {med_sep:.1f} px")

    # 1. Tracing controls (all explicit parameters now)
    fill_missing_orders = True

    # 2. Background Envelope Estimation
    # Use a large-scale rolling minimum to capture the background base.
    from scipy.ndimage import minimum_filter1d
    bg_size = int(2.5 * med_sep)
    envelope = minimum_filter1d(center_profile, size=bg_size)
    envelope = gaussian_filter1d(envelope, sigma=bg_size/4.0) # Smooth it
    
    # 3. Background-subtracted profile for seeding
    seed_pure = center_profile - envelope
    
    # Robust noise estimation on the residual
    def _estimate_prof_noise(prof):
        valid_prof = prof[np.isfinite(prof)]
        if len(valid_prof) < 10:
            return 1.0, 0.0
            
        # 1. Baseline estimation: bottom 25%
        sorted_prof = np.sort(valid_prof)
        bg_pixels = sorted_prof[:max(10, len(sorted_prof) // 4)]
        s_med = float(np.nanmedian(bg_pixels))
        
        # 2. Intrinsic noise estimation via MAD of differences
        diffs = np.diff(valid_prof)
        s_sig = float(np.nanmedian(np.abs(diffs)) / 0.9539) # sqrt(2)*0.6745
        
        return max(s_sig, 0.1), s_med

    noise, baseline = _estimate_prof_noise(seed_pure)
    seed_profile = seed_pure # Use the pure signal for peak finding
    
    prof_max = float(np.max(seed_profile))

    # Seeding threshold: single robust SNR-based gate
    det_threshold = baseline + snr_threshold * noise
    min_dist = max(3, int(med_sep * 0.20)) 

    peaks, _ = find_peaks(
        seed_profile,
        height=det_threshold,
        distance=min_dist,
        prominence=max(0.8 * noise, 0.15 * snr_threshold * noise),
    )

    logger.info(
        f"Found {len(peaks)} seed peaks (SNR > {snr_threshold:.1f}, "
        f"baseline={baseline:.4f}, noise={noise:.4f})"
    )
    
    if len(peaks) > 0:
        # Diagnostic: check the faintest detected peak and first peak region
        pk_vals = seed_profile[peaks]
        min_pk_idx = np.argmin(pk_vals)
        first_pk = int(peaks[0])
        
        logger.info(
            f"Faintest seed: y={peaks[min_pk_idx]}, val={pk_vals[min_pk_idx]:.4f}, "
            f"SNR={(pk_vals[min_pk_idx]-baseline)/(noise+1e-6):.2f}"
        )
        # Show what's below the first peak.
        check_y = max(0, first_pk - int(med_sep))
        for yy in range(check_y, first_pk):
            if seed_profile[yy] > baseline + 0.3 * noise:
                logger.info(
                    f"  Below-first: y={yy}, val={seed_profile[yy]:.4f}, "
                    f"SNR={(seed_profile[yy]-baseline)/(noise+1e-6):.2f}"
                )

    peaks = np.asarray(sorted(peaks), dtype=int)

    # -------------------------------------------------------------------------
    # 4. Fill missed weak orders in 1D cross-section
    # -------------------------------------------------------------------------
    if fill_missing_orders and len(peaks) >= 2:
        # User requirement: use 2.5 sigma for faint blind guessing.
        gap_fill_snr = max(2.5, gap_fill_snr)

        logger.info(
            f"Gap-fill analysis: {len(peaks)} seeds, global median_sep={med_sep:.1f}px, "
            f"gap_fill_factor=1.30, gap_fill_snr={gap_fill_snr:.2f}"
        )

        # Internal gap filling using curve
        if np.isfinite(med_sep) and med_sep > 2:
            inserted = []
            for i in range(len(peaks) - 1):
                left, right = peaks[i], peaks[i + 1]
                gap = right - left
                expected_sep = local_sep_func(0.5 * (left + right))
                ratio = gap / expected_sep
                
                # Gap > 1.3 curve prediction -> Missing orders
                if gap > 1.30 * expected_sep:
                    n_missing = int(round(gap / expected_sep)) - 1
                    if n_missing < 1: 
                        n_missing = 1
                    logger.info(
                        f"Gap-fill: gap between y={left} and y={right} "
                        f"is {gap}px ({ratio:.2f}x expected {expected_sep:.1f}), inserting {n_missing} candidate(s)"
                    )
                    for k in range(1, n_missing + 1):
                        guess = int(round(left + k * gap / (n_missing + 1)))
                        win = max(5, int(0.55 * expected_sep))
                        y1 = max(0, guess - win)
                        y2 = min(h, guess + win + 1)
                        if y2 - y1 < 3:
                            continue
                        seg = seed_profile[y1:y2]
                        seg_med = float(np.median(seg))
                        seg_mad = float(np.median(np.abs(seg - seg_med)))
                        seg_noise = 1.4826 * seg_mad if seg_mad > 0 else float(np.std(seg - seg_med))
                        
                        loc = int(np.argmax(seg))
                        cand = y1 + loc
                        cand_val = seg[loc]
                        
                        if seg_noise > 0 and cand_val < (seg_med + gap_fill_snr * seg_noise):
                            logger.debug(
                                f"  Gap-fill candidate y={cand}: rejected "
                                f"(SNR={(cand_val-seg_med)/(seg_noise+1e-6):.1f} < {gap_fill_snr})"
                            )
                            continue
                            
                        # Avoid duplicates
                        if np.all(np.abs(peaks - cand) > max(3, int(0.30 * expected_sep))):
                            inserted.append(cand)
                            logger.info(f"  Gap-fill candidate y={cand}: ACCEPTED (SNR={(cand_val-seg_med)/(seg_noise+1e-6):.1f})")

            # Edge extrapolation
            min_sep = max(3, int(0.30 * med_sep))

            left_ref = int(peaks[0])
            local_sep_bottom = local_sep_func(left_ref)
            edge_win_bottom = max(5, int(0.80 * local_sep_bottom))
            logger.info(
                f"Edge-fill (bottom) starting: left_ref=y{left_ref}, "
                f"local_sep={local_sep_bottom:.1f}, edge_win={edge_win_bottom}"
            )
            while left_ref - local_sep_func(left_ref) * 0.5 > 0:
                step_sep = local_sep_func(left_ref)
                edge_win = max(5, int(0.80 * step_sep))
                guess = int(round(left_ref - step_sep))
                if guess < 0:
                    guess = max(0, left_ref // 2)
                y1 = max(0, guess - edge_win)
                y2 = min(h, guess + edge_win + 1)
                if y2 - y1 < 3:
                    break
                seg = seed_profile[y1:y2]
                cand = y1 + int(np.argmax(seg))
                seg_med = float(np.median(seg))
                seg_mad = float(np.median(np.abs(seg - seg_med)))
                seg_noise = 1.4826 * seg_mad if seg_mad > 0 else float(np.std(seg - seg_med))
                cand_val = float(seed_profile[cand])
                threshold = seg_med + snr_threshold * seg_noise
                logger.info(
                    f"  Edge-fill (bottom): guess=y{guess}, window=[{y1},{y2}], "
                    f"argmax→y{cand}, val={cand_val:.4f}, "
                    f"seg_med={seg_med:.4f}, seg_noise={seg_noise:.4f}, "
                    f"threshold={threshold:.4f} ({snr_threshold:.1f}σ), left_ref=y{left_ref}"
                )
                if seg_noise > 0 and cand_val >= threshold:
                    if np.all(np.abs(peaks - cand) > min_sep):
                        inserted.append(cand)
                        logger.info(
                            f"Edge-fill (bottom): y={cand} ACCEPTED "
                            f"(val={cand_val:.3f} >= {threshold:.3f})"
                        )
                        left_ref = cand
                        continue
                logger.debug(
                    f"Edge-fill (bottom): y={cand} rejected "
                    f"(val={cand_val:.3f} < {threshold:.3f})"
                )
                break

            right_ref = int(peaks[-1])
            local_sep_top = local_sep_func(right_ref)
            edge_win_top = max(5, int(0.80 * local_sep_top))
            logger.info(
                f"Edge-fill (top) starting: right_ref=y{right_ref}, "
                f"local_sep={local_sep_top:.1f}, edge_win={edge_win_top}"
            )
            while right_ref + local_sep_func(right_ref) * 0.5 < h:
                step_sep = local_sep_func(right_ref)
                edge_win = max(5, int(0.80 * step_sep))
                guess = int(round(right_ref + step_sep))
                if guess >= h:
                    break
                y1 = max(0, guess - edge_win)
                y2 = min(h, guess + edge_win + 1)
                if y2 - y1 < 3:
                    break
                seg = seed_profile[y1:y2]
                cand = y1 + int(np.argmax(seg))
                seg_med = float(np.median(seg))
                seg_mad = float(np.median(np.abs(seg - seg_med)))
                seg_noise = 1.4826 * seg_mad if seg_mad > 0 else float(np.std(seg - seg_med))
                cand_val = float(seed_profile[cand])
                threshold = seg_med + snr_threshold * seg_noise
                if seg_noise > 0 and cand_val >= threshold:
                    if np.all(np.abs(peaks - cand) > min_sep):
                        inserted.append(cand)
                        logger.info(
                            f"Edge-fill (top): y={cand} ACCEPTED "
                            f"(val={cand_val:.3f} >= {threshold:.3f}, {snr_threshold:.1f}σ)"
                        )
                        right_ref = cand
                        continue
                logger.debug(
                    f"Edge-fill (top): y={cand} rejected "
                    f"(val={cand_val:.3f} < {threshold:.3f})"
                )
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
        peak_val = float(segment[loc])
        if not np.isfinite(peak_val) or peak_val <= 0:
            return None

        # Linear space intensity-weighted centroid
        bg = 0.5 * (float(segment[0]) + float(segment[-1]))
        sub_flux = segment.astype(float) - bg
        core_mask = sub_flux > 0.2 * (peak_val - bg)
        
        if np.sum(core_mask) >= 3:
            yy = np.arange(y1, y2, dtype=float)[core_mask]
            weights = sub_flux[core_mask]
            sum_w = np.sum(weights)
            if sum_w > 0:
                return float(np.sum(yy * weights) / sum_w)

        # Fallback
        return float(y1 + loc)

    def trace_order(seed_y: float) -> tuple:
        # Trace every column so the downstream sigma-clipping and Chebyshev fit
        # have dense coverage; faint-order wings use prediction + rejection to
        # avoid locking onto neighbouring orders.
        local_sep = local_sep_func(seed_y)
        search_half = max(4, int(local_sep * 0.45))

        # Refine the starting seed position to sub-pixel accuracy before tracing
        refined_y = refine_peak_y(img[:, x_center], seed_y, search_half)
        start_y = refined_y if refined_y is not None else float(seed_y)

        xs = [float(x_center)]
        ys = [start_y]

        # right
        y_prev = start_y
        miss_count = 0
        right_x = [float(x_center)]
        right_y = [start_y]
        for x in range(x_center + 1, w):
            if len(right_x) >= 4:
                recent_x = np.asarray(right_x[-4:], dtype=float)
                recent_y = np.asarray(right_y[-4:], dtype=float)
                coef1 = np.polyfit(recent_x, recent_y, deg=1)
                y_guess = float(np.polyval(coef1, x))
            else:
                y_guess = y_prev

            y_new = refine_peak_y(img[:, x], y_guess, search_half)
            if y_new is None:
                miss_count += 1
                if miss_count > 10:
                    break
                continue

            if abs(y_new - y_guess) > search_half:
                miss_count += 1
                if miss_count > 10:
                    break
                continue

            xs.append(x)
            ys.append(y_new)
            right_x.append(float(x))
            right_y.append(float(y_new))
            y_prev = y_new
            miss_count = 0

        # left
        y_prev = start_y
        miss_count = 0
        left_x = [float(x_center)]
        left_y = [start_y]
        for x in range(x_center - 1, -1, -1):
            if len(left_x) >= 4:
                recent_x = np.asarray(left_x[-4:], dtype=float)
                recent_y = np.asarray(left_y[-4:], dtype=float)
                coef1 = np.polyfit(recent_x, recent_y, deg=1)
                y_guess = float(np.polyval(coef1, x))
            else:
                y_guess = y_prev

            y_new = refine_peak_y(img[:, x], y_guess, search_half)
            if y_new is None:
                miss_count += 1
                if miss_count > 10:
                    break
                continue

            if abs(y_new - y_guess) > search_half:
                miss_count += 1
                if miss_count > 10:
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

    def estimate_order_width(xs: np.ndarray, ys: np.ndarray, local_sep: float) -> float:
        """Estimate aperture width from local FWHM samples along the traced order."""
        if xs.size == 0:
            return float(local_sep)

        sample_idx = np.linspace(0, xs.size - 1, num=min(24, xs.size), dtype=int)
        widths = []
        search_half = max(8, int(0.8 * local_sep))

        for i in sample_idx:
            x = int(round(xs[i]))
            yc = float(ys[i])
            y1 = max(0, int(np.floor(yc - search_half)))
            y2 = min(h, int(np.ceil(yc + search_half + 1)))
            if y2 - y1 < 7:
                continue

            prof = img[y1:y2, x]
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
            if 2.0 <= fwhm <= 2.5 * local_sep:
                widths.append(fwhm)

        if not widths:
            return float(local_sep)

        # Use aperture full width ~ 1.8*FWHM, clipped to sensible range.
        width = 1.8 * float(np.median(widths))
        width = float(np.clip(width, 6.0, 2.2 * local_sep))
        return width

    def _fit_center_to_coefs(x_vals: np.ndarray, y_vals: np.ndarray, 
                             degree: int, domain: tuple) -> np.ndarray:
        """Fit center using Chebyshev and return coefficients (ascending)."""
        x = np.asarray(x_vals, dtype=float)
        y = np.asarray(y_vals, dtype=float)
        
        if x.size < 4:
            return None

        deg = max(1, int(min(degree, max(1, x.size - 2))))
        cheb = Chebyshev.fit(x, y, deg=deg, domain=domain)
        return np.array(cheb.coef, dtype=float)
    def _fit_boundary_coeff(x_vals: np.ndarray, y_vals: np.ndarray,
                            degree: int, domain: tuple) -> np.ndarray:
        """Fit a smooth boundary using Chebyshev polynomials with robust outlier rejection."""
        x = np.asarray(x_vals, dtype=float)
        y = np.asarray(y_vals, dtype=float)
        if x.size < 4:
            return None
        
        deg = max(1, int(min(degree, max(1, x.size - 2))))
        try:
            # Iterative fit to handle outliers near detector edges
            cheb_obj = Chebyshev.fit(x, y, deg=deg, domain=domain)
            for _ in range(4):
                res = y - cheb_obj(x)
                s = np.std(res)
                if s <= 0.02: break
                good = np.abs(res) < 3.0 * s
                if np.sum(good) < deg + 2: break
                cheb_obj = Chebyshev.fit(x[good], y[good], deg=deg, domain=domain)
            return np.array(cheb_obj.coef, dtype=float)
        except:
            return None

    # 3. Tracing controls (all explicit parameters now)
    fill_missing_orders = True

    def _build_moffat_boundaries(cands: list,
                                   flux_fraction: float = 0.02,
                                   fwhm_scale: float = 1.5) -> None:
        """Determine order boundaries via two-step Moffat fitting.

        Step 1: For each sample column, use the traced center position as the
                peak.  Walk outward from that center in the log-compressed
                cross-section until the signal drops to *flux_fraction* × peak.
                This gives an initial lower/upper bracket.
        Step 2: Background = average of profile values at the two bracket
                positions.  Fit a Moffat profile to the background-subtracted
                cross-section within the bracket.
                Final boundary = fitted_center ± fwhm_scale × FWHM.
                Falls back to the Step-1 bracket when fitting fails.
        """
        from scipy.optimize import curve_fit

        def _moffat1d(y, y_c, I0, gamma, beta):
            return I0 * (1.0 + ((y - y_c) / gamma) ** 2) ** (-beta)

        if not cands:
            return

        x_samples = np.linspace(0, w - 1, num=min(128, w), dtype=int)

        for cand in cands:
            c_coef = cand['center_coef']
            c_model = Chebyshev(c_coef, domain=(0, w - 1))
            local_sep = local_sep_func(cand['y_center'])
            is_reg = cand.get('is_regularized', False)

            xl, yl, xu, yu = [], [], [], []
            xc, yc = [], []
            half_win = max(8, int(0.6 * local_sep))

            for x in x_samples:
                yc_float = float(c_model(x))          # traced center = peak
                y0 = int(round(yc_float))
                y1 = max(0, y0 - half_win)
                y2 = min(h, y0 + half_win + 1)
                if y2 - y1 < 5:
                    continue

                prof = img[y1:y2, x].astype(float)
                n = len(prof)
                mid_idx = max(0, min(n - 1, y0 - y1))  # center in local coords
                peak_val = float(prof[mid_idx])

                # ── Step 1: bracket at flux_fraction × peak ───────────────────
                threshold = flux_fraction * peak_val

                low_idx = 0
                for j in range(mid_idx, -1, -1):
                    if prof[j] <= threshold:
                        low_idx = j
                        break

                up_idx = n - 1
                for j in range(mid_idx, n):
                    if prof[j] <= threshold:
                        up_idx = j
                        break

                # ── Step 2: Moffat fit within bracket ─────────────────────────
                bg = 0.5 * (float(prof[low_idx]) + float(prof[up_idx]))
                seg_y = np.arange(low_idx, up_idx + 1, dtype=float)
                seg_f = np.maximum(prof[low_idx:up_idx + 1] - bg, 0.0)

                fwhm_final = None
                yc_fit = float(mid_idx)

                try:
                    if len(seg_y) >= 5 and (peak_val - bg) > 0:
                        I0_guess = max(float(peak_val - bg), 1e-6)
                        gamma_guess = max(1.0, (up_idx - low_idx) / 4.0)
                        p0 = [float(mid_idx), I0_guess, gamma_guess, 2.5]
                        b_lo = [float(low_idx), 0.0, 0.3, 0.5]
                        b_hi = [float(up_idx), I0_guess * 5.0,
                                float(up_idx - low_idx), 10.0]
                        popt, _ = curve_fit(
                            _moffat1d, seg_y, seg_f,
                            p0=p0, bounds=(b_lo, b_hi),
                            maxfev=3000
                        )
                        yc_fit, _I0f, gamma_fit, beta_fit = popt
                        if beta_fit > 0.5 and gamma_fit > 0.3:
                            fwhm_final = 2.0 * gamma_fit * np.sqrt(
                                2.0 ** (1.0 / beta_fit) - 1.0)
                except Exception:
                    pass

                if fwhm_final is not None and 1.0 <= fwhm_final <= local_sep:
                    # Moffat fit is good, trust its center even for faint orders
                    y_center_actual = float(y1 + yc_fit)
                    half_bnd = fwhm_scale * fwhm_final
                    y_lo = y_center_actual - half_bnd
                    y_hi = y_center_actual + half_bnd
                else:
                    # Moffat fit failed or is unreliable, stick to the polynomial trace
                    y_center_actual = yc_float
                    # Fallback to Step-1 bracket, forced symmetric for faint orders
                    if is_reg:
                        half_bnd = (up_idx - low_idx) / 2.0
                        y_lo = y_center_actual - half_bnd
                        y_hi = y_center_actual + half_bnd
                    else:
                        y_lo = float(y1 + low_idx)
                        y_hi = float(y1 + up_idx)

                xl.append(float(x))
                yl.append(y_lo)
                xu.append(float(x))
                yu.append(y_hi)
                # Always collect the best available center
                xc.append(float(x))
                yc.append(y_center_actual)

            # Upgrade the center trace using the new, more accurate centers
            if len(xc) >= 5:
                new_c_coef = _fit_boundary_coeff(np.array(xc), np.array(yc), trace_degree, (0, w - 1))
                if new_c_coef is not None:
                    cand['center_coef'] = new_c_coef

            cand['lower_coef'] = _fit_boundary_coeff(
                np.array(xl), np.array(yl), trace_degree, (0, w - 1))
            cand['upper_coef'] = _fit_boundary_coeff(
                np.array(xu), np.array(yu), trace_degree, (0, w - 1))

    def _plot_seeding_diagnostics(out_path):
        """Plot the seeding profile with envelope and thresholds restored to original state."""
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(12, 6))
        y_idx = np.arange(len(center_profile))
        
        ax.plot(y_idx, center_profile, color='gray', alpha=0.5, label='Raw Profile (log)')
        ax.plot(y_idx, envelope, color='blue', linestyle='--', alpha=0.8, label='BG Envelope')
        
        det_line = envelope + snr_threshold * noise
        root_line = envelope + 3.0 * noise
        
        ax.plot(y_idx, det_line, color='orange', label=f'Detection ({snr_threshold}σ)')
        ax.plot(y_idx, root_line, color='green', alpha=0.6, label='Aperture Boundary (3σ)')
        
        # Highlight detected peaks
        if len(peaks) > 0:
            ax.scatter(peaks, center_profile[peaks], color='red', s=15, zorder=5, label='Seeds')
            
        ax.set_title("Order Seeding & Threshold Diagnostics (Restored with Envelope)")
        ax.set_xlabel("Detector Y (px)")
        ax.set_ylabel("Intensity (log units)")
        ax.legend(loc='upper right', fontsize='small')
        ax.grid(True, alpha=0.2)
        
        plt.tight_layout()
        plt.savefig(out_path, dpi=150)
        plt.close(fig)
        logger.info(f"Saved seeding diagnostics plot to {out_path}")

 # Helpers are defined, execution will happen at the end of the function.

    aperture_set = ApertureSet()
    candidates = []

    # 读取谱级间距容差参数
    import configparser
    spacing_tol = 0.3
    try:
        from src.config.config_manager import ConfigManager
        cfg = ConfigManager()
        spacing_tol = cfg.get_float('reduce.trace', 'spacing_tol', 0.3)
    except Exception:
        pass

    for order_idx, seed in enumerate(peaks):
        local_sep = local_sep_func(float(seed))
        x_positions, y_positions = trace_order(float(seed))

        if len(x_positions) < 15:
            logger.warning(f"Seed {order_idx + 1} (y~{seed:.1f}): Rejected - Not enough traced points ({len(x_positions)} < 15)")
            continue

        try:
            # 1. Iterative outlier rejection using一个简单多项式
            good = np.ones_like(x_positions, dtype=bool)
            deg_rough = 3
            for _ in range(3):
                if np.sum(good) < deg_rough + 2: break
                coef = np.polyfit(x_positions[good], y_positions[good], deg=deg_rough)
                resid = y_positions - np.polyval(coef, x_positions)
                std = np.std(resid[good])
                if std <= 0.05: break
                good = np.abs(resid) < 3.0 * std

            if np.sum(good) < 15:
                logger.warning(f"Seed {order_idx + 1} (y~{seed:.1f}): Rejected - Too few valid points after clipping ({np.sum(good)} < 15)")
                continue

            # 2. Final Center fitting (Always Chebyshev)
            center_deg = trace_degree
            center_coef = _fit_center_to_coefs(x_positions[good], y_positions[good], center_deg, (0, float(w - 1)))
            if center_coef is None:
                logger.warning(f"Seed {order_idx + 1} (y~{seed:.1f}): Rejected - Chebyshev center fit failed")
                continue

            center_model = Chebyshev(center_coef, domain=(0, w - 1))
            yc_mid = center_model(x_center)

            x_coverage = float(np.max(x_positions[good]) - np.min(x_positions[good]))
            min_cov_px = float(min_trace_coverage) * w
            if x_coverage < min_cov_px:
                logger.warning(f"Seed {order_idx + 1} (y~{seed:.1f}): Rejected - Trace coverage too short ({x_coverage:.0f}px < min {min_cov_px:.0f}px)")
                continue


            width = estimate_order_width(x_positions[good], y_positions[good], local_sep)
            half_width = width / 2.0
            lower_coef = np.array([yc_mid - half_width], dtype=float)
            upper_coef = np.array([yc_mid + half_width], dtype=float)

            candidates.append({
                'center_coef': center_coef,
                'lower_coef': lower_coef,
                'upper_coef': upper_coef,
                'width': float(width),
                'n_good': int(np.sum(good)),
                'fit_rms': float(np.std(y_positions[good] - center_model(x_positions[good]))),
                'x_coverage': float(x_coverage),
                'y_center': float(yc_mid),
                'seed_y': float(seed),
            })

            logger.info(
                f"Candidate from seed {order_idx + 1}: points={np.sum(good)}, "
                f"fit=chebyshev, width={width:.1f}px"
            )

        except Exception as e:
            logger.warning(f"Order {order_idx + 1}: Failed polynomial tracing: {e}")
            continue

    # Merge duplicate candidates that represent the same physical order,
    # then assign continuous order IDs.
    # Ensure boundary_frac and fwhm_scale are available in this scope (from function arguments)
    boundary_frac = boundary_frac if 'boundary_frac' in locals() else 0.02
    fwhm_scale = fwhm_scale if 'fwhm_scale' in locals() else 1.5
    if candidates:
        candidates_sorted = sorted(candidates, key=lambda c: c['y_center'])
        unique = []

        # 恢复并增强谱级筛选判据
        if candidates:
            candidates_sorted = sorted(candidates, key=lambda c: c['y_center'])
            unique = []

            for idx, cand in enumerate(candidates_sorted):
                if not unique:
                    unique.append(cand)
                    continue

                prev = unique[-1]
                pred_sep = local_sep_func(0.5 * (cand['y_center'] + prev['y_center']))
                actual_sep = abs(cand['y_center'] - prev['y_center'])
                # 输出调试信息
                logger.info(f"Seed spacing debug: idx={idx}, y1={prev['y_center']:.2f}, y2={cand['y_center']:.2f}, actual={actual_sep:.2f}, predicted={pred_sep:.2f}, ratio={actual_sep/pred_sep:.2f}")
                # 判据：实际间距小于0.3×预测间距则合并（保留分数高的），大于10.0×预测间距才排除
                if actual_sep < 0.3 * pred_sep:
                    prev_score = prev['n_good'] - 2.0 * prev['fit_rms'] + 0.002 * prev['x_coverage']
                    cand_score = cand['n_good'] - 2.0 * cand['fit_rms'] + 0.002 * cand['x_coverage']
                    if cand_score > prev_score:
                        unique[-1] = cand
                elif actual_sep > 10.0 * pred_sep:
                    logger.warning(f"Seed at y={cand['y_center']:.1f} rejected: spacing {actual_sep:.2f} > 10.0×predicted {pred_sep:.2f}")
                    continue
                else:
                    unique.append(cand)

        # ---------------------------------------------------------------------
        # Global Shape Regularization (Prevents faint orders from crossing)
        # ---------------------------------------------------------------------
        if len(unique) >= 3:
            def is_robust_shape(c):
                # Only use extremely well-traced orders across most of the CCD as global shape anchors.
                # Loose RMS allows wiggly noise-chasing traces to corrupt the global shape.
                return (c['x_coverage'] > 0.85 * w) and (c['fit_rms'] < max(1.5, 0.08 * local_sep_func(c['y_center'])))

            robust = [c for c in unique if is_robust_shape(c)]
            if len(robust) >= 3:
                max_len = max(len(c['center_coef']) for c in robust)
                # Ensure robust is strictly sorted by y_center so np.interp works correctly
                robust_sorted = sorted(robust, key=lambda c: c['y_center'])
                y_c_robust = np.array([c['y_center'] for c in robust_sorted])
                coefs_robust = np.array([np.pad(c['center_coef'], (0, max_len - len(c['center_coef']))) for c in robust_sorted])
                
                for c in unique:
                    # 如果该级次本身就很亮且轨迹稳健，则保留其自身真实的光学畸变轨迹
                    if is_robust_shape(c):
                        c['is_regularized'] = False
                        continue

                    c_model_old = Chebyshev(c['center_coef'], domain=(0, w - 1))
                    # Anchor directly to the robust photometric seed peak found in the median profile
                    y_anchor = c['y_center']
                    
                    new_coef = np.zeros(max_len)
                    for d in range(1, max_len):
                        # np.interp guarantees flat extrapolation at edges:
                        # Faint edge orders will exactly inherit the shape of the outermost bright order.
                        new_coef[d] = np.interp(y_anchor, y_c_robust, coefs_robust[:, d])
                        
                    c_model_new_temp = Chebyshev(new_coef, domain=(0, w - 1))
                    y_mid_new_temp = c_model_new_temp(x_center)
                    
                    # T_0(x) is 1 everywhere, so adjusting c[0] shifts the whole curve exactly
                    new_coef[0] = y_anchor - y_mid_new_temp
                    c['center_coef'] = new_coef
                    c['y_center'] = y_anchor
                    c['is_regularized'] = True
                    
                logger.info(f"Applied global shape regularization to {len(unique) - len(robust)} faint orders based on {len(robust)} robust traces.")

        # Apply Moffat boundary detection to the unique orders
        _build_moffat_boundaries(unique, flux_fraction=boundary_frac, fwhm_scale=fwhm_scale)

        # ---------------------------------------------------------------------
        # 2D Coefficient Interpolation for Completely Missing Orders
        # ---------------------------------------------------------------------
        if len(unique) >= 2:
            final_cands = []
            final_cands.append(unique[0])
            
            for i in range(len(unique) - 1):
                left_cand = unique[i]
                right_cand = unique[i + 1]
                gap = right_cand['y_center'] - left_cand['y_center']
                
                expected_sep = local_sep_func(0.5 * (left_cand['y_center'] + right_cand['y_center']))
                
                if gap > 1.30 * expected_sep:
                    n_missing = int(round(gap / expected_sep)) - 1
                    if n_missing > 0:
                        logger.info(
                            f"2D Interpolation: gap between y={left_cand['y_center']:.1f} "
                            f"and y={right_cand['y_center']:.1f} is {gap:.1f}px "
                            f"(expected {expected_sep:.1f}), synthesising {n_missing} completely blind order(s)."
                        )
                        for k in range(1, n_missing + 1):
                            alpha = float(k) / (n_missing + 1.0)
                            interp_cand = {}
                            # Pad to same length if degrees somehow differ (though they should be standard center_deg)
                            len_c = max(len(left_cand['center_coef']), len(right_cand['center_coef']))
                            lc_pad = np.pad(left_cand['center_coef'], (0, len_c - len(left_cand['center_coef'])))
                            rc_pad = np.pad(right_cand['center_coef'], (0, len_c - len(right_cand['center_coef'])))
                            interp_cand['center_coef'] = (1.0 - alpha) * lc_pad + alpha * rc_pad
                            
                            len_l = max(len(left_cand['lower_coef']), len(right_cand['lower_coef']))
                            ll_pad = np.pad(left_cand['lower_coef'], (0, len_l - len(left_cand['lower_coef'])))
                            rl_pad = np.pad(right_cand['lower_coef'], (0, len_l - len(right_cand['lower_coef'])))
                            interp_cand['lower_coef'] = (1.0 - alpha) * ll_pad + alpha * rl_pad
                            
                            len_u = max(len(left_cand['upper_coef']), len(right_cand['upper_coef']))
                            lu_pad = np.pad(left_cand['upper_coef'], (0, len_u - len(left_cand['upper_coef'])))
                            ru_pad = np.pad(right_cand['upper_coef'], (0, len_u - len(right_cand['upper_coef'])))
                            interp_cand['upper_coef'] = (1.0 - alpha) * lu_pad + alpha * ru_pad
                            
                            interp_cand['width'] = (1.0 - alpha) * left_cand['width'] + alpha * right_cand['width']
                            interp_cand['y_center'] = (1.0 - alpha) * left_cand['y_center'] + alpha * right_cand['y_center']
                            interp_cand['n_good'] = 0
                            interp_cand['fit_rms'] = 0.0
                            interp_cand['x_coverage'] = float(w)
                            interp_cand['is_interpolated'] = True
                            final_cands.append(interp_cand)
                            
                final_cands.append(right_cand)
            
            unique = final_cands

        # Save diagnostic plot (showing raw profile + envelope + thresholds)
        # NOTE: seeding diagnostics are now merged into apertures_q*.png via
        # _plot_apertures_cross_section(); do not regenerate a separate file here.

        for new_id, cand in enumerate(unique, start=1):
            aperture_loc = ApertureLocation(
                aperture=new_id,
                order=new_id,
                center_coef=cand.get('center_coef'),
                width=cand['width'],
                is_chebyshev=True, # native chebyshev for both center and boundaries
                domain=(0, float(w - 1)),
                is_interpolated=cand.get('is_interpolated', False)
            )
            aperture_set.add_aperture(aperture_loc)

    logger.info(f"Detected {aperture_set.norders} orders (3-sigma boundary mode)")
    return aperture_set


# ---------------------------------------------------------------------------
#  Spacing-based extrapolation for virtual edge orders
# ---------------------------------------------------------------------------

def _extrapolate_virtual_orders(apertures: ApertureSet,
                                image_width: int, image_height: int,
                                n_orders: int,
                                side: str = 'below',
                                n_ref: int = 8,
                                trace_degree: int = 3) -> list:
    """Predict virtual orders by extrapolating inter-order spacing.

    Instead of evaluating a global 2D Chebyshev surface (which diverges at
    the edges), this function uses the **local spacing pattern** of the
    nearest ``n_ref`` real orders to predict where the next orders would be.

    For each column independently:
      1. Evaluate centre / lower / upper of the nearest ``n_ref`` orders.
      2. Compute per-column inter-order spacings Δy(m).
      3. Fit Δy as a low-degree polynomial of m (degree 1 or 2) to capture
         the smooth trend dictated by the grating equation.
      4. Extrapolate centre positions one at a time:
         y_{next} = y_last + Δy_extrapolated.
      5. Similarly propagate upper/lower boundary widths.

    Parameters
    ----------
    apertures : ApertureSet with real (positive-ID) orders only.
    image_width, image_height : detector dimensions.
    n_orders : how many virtual orders to generate.
    side : 'below' (extrapolate below first order) or 'above' (above last).
    n_ref : number of reference orders from the edge to use (default 8).
    trace_degree : polynomial degree for per-order trace coefficients.

    Returns
    -------
    list of ApertureLocation objects (with negative IDs) in order of
    increasing distance from the detector edge.  The list may be shorter
    than ``n_orders`` if all remaining virtual orders would be entirely
    off-detector.
    """
    # Collect real (positive-ID) orders sorted by centre Y at mid-column.
    mid_col = image_width // 2
    real_oids = [oid for oid in apertures.get_orders() if oid > 0]
    real_oids.sort(key=lambda oid: float(apertures.get_aperture(oid).get_position(mid_col)))

    if len(real_oids) < 3:
        return []

    # Select reference orders from the relevant edge.
    if side == 'below':
        ref_oids = real_oids[:min(n_ref, len(real_oids))]
    else:
        ref_oids = real_oids[max(0, len(real_oids) - n_ref):]

    n_r = len(ref_oids)
    x_all = np.arange(image_width, dtype=float)

    # Evaluate centre/lower/upper for each reference order at every column.
    ref_cen = np.empty((n_r, image_width), dtype=np.float64)  # shape (n_r, W)
    ref_lo  = np.empty((n_r, image_width), dtype=np.float64)
    ref_hi  = np.empty((n_r, image_width), dtype=np.float64)
    for i, oid in enumerate(ref_oids):
        ap = apertures.get_aperture(oid)
        ref_cen[i] = ap.get_position(x_all)
        ref_lo[i]  = ap.get_lower(x_all)
        ref_hi[i]  = ap.get_upper(x_all)

    # Per-column spacings between consecutive reference orders.
    # spacings[j, xi] = centre[j+1, xi] - centre[j, xi]  (j = 0..n_r-2)
    spacings = np.diff(ref_cen, axis=0)                # (n_r-1, W)
    wup_arr  = ref_hi - ref_cen                         # half-width above centre
    wlo_arr  = ref_cen - ref_lo                         # half-width below centre

    # ----- Vectorised linear fit across all columns at once -----
    # For spacing: fit Δy(j) ≈ a_sp + b_sp * j  (j = 0..n_r-2)
    # For width-up / width-lo: fit w(j) ≈ a_w + b_w * j  (j = 0..n_r-1)
    # Linear regression: b = (N*Σ(xy) - Σx*Σy) / (N*Σ(x²) - (Σx)²),
    #                    a = (Σy - b*Σx) / N
    def _vectorised_linfit(j, data):
        """Linear fit y = a + b*j for data shape (len(j), W). Returns (a, b) each (W,)."""
        N = len(j)
        Sx  = j.sum()
        Sx2 = (j * j).sum()
        Sy  = data.sum(axis=0)                          # (W,)
        Sxy = (j[:, None] * data).sum(axis=0)           # (W,)
        denom = N * Sx2 - Sx * Sx
        if abs(denom) < 1e-12:
            return Sy / max(N, 1), np.zeros_like(Sy)
        b = (N * Sxy - Sx * Sy) / denom
        a = (Sy - b * Sx) / N
        return a, b                                     # each shape (W,)

    j_sp = np.arange(n_r - 1, dtype=float)
    j_w  = np.arange(n_r, dtype=float)

    if n_r >= 4:
        a_sp, b_sp = _vectorised_linfit(j_sp, spacings)   # (W,) each
        a_wu, b_wu = _vectorised_linfit(j_w,  wup_arr)
        a_wl, b_wl = _vectorised_linfit(j_w,  wlo_arr)
    else:
        a_sp = np.median(spacings, axis=0)
        b_sp = np.zeros(image_width)
        a_wu = np.median(wup_arr, axis=0)
        b_wu = np.zeros(image_width)
        a_wl = np.median(wlo_arr, axis=0)
        b_wl = np.zeros(image_width)

    virtual_orders = []

    for k in range(1, n_orders + 1):
        # Accumulate predicted spacings from the anchor order.
        if side == 'below':
            # Extrapolate to j = -1, -2, ..., -k.
            cum_sp = np.zeros(image_width, dtype=np.float64)
            for kk in range(1, k + 1):
                cum_sp += a_sp + b_sp * (-float(kk))
            pred_cen = ref_cen[0] - np.abs(cum_sp)
            w_target = -float(k)
        else:
            cum_sp = np.zeros(image_width, dtype=np.float64)
            for kk in range(1, k + 1):
                cum_sp += a_sp + b_sp * float(n_r - 2 + kk)
            pred_cen = ref_cen[-1] + np.abs(cum_sp)
            w_target = float(n_r - 1 + k)

        pred_wup = np.maximum(2.0, a_wu + b_wu * w_target)
        pred_wlo = np.maximum(2.0, a_wl + b_wl * w_target)

        pred_upper = pred_cen + pred_wup
        pred_lower = pred_cen - pred_wlo

        # Check whether ANY part of this order overlaps the detector.
        if side == 'below':
            if np.max(pred_upper) < 0:
                break
        else:
            if np.min(pred_lower) >= image_height:
                break

        # Fit polynomial coefficients for the trace.
        cen_coef = np.polyfit(x_all, pred_cen, trace_degree)
        low_coef = np.polyfit(x_all, pred_lower, trace_degree)
        up_coef  = np.polyfit(x_all, pred_upper, trace_degree)
        width = float(np.median(pred_wup + pred_wlo))

        yc_mid = pred_cen[mid_col]
        ap = ApertureLocation(
            aperture=0,    # placeholder; caller sets the final ID
            order=0,
            center_coef=cen_coef,
            lower_coef=low_coef,
            upper_coef=up_coef,
            width=max(4.0, width),
        )
        virtual_orders.append(ap)
        logger.info(f"  Virtual order ({side}) k={k}: y_center≈{yc_mid:.1f}, "
                    f"width≈{width:.1f}")

    return virtual_orders


def extend_apertures_for_background(apertures: ApertureSet,
                                    image_width: int, image_height: int,
                                    n_below: int = 2,
                                    n_above: int = 1,
                                    trace_degree: int = 3) -> ApertureSet:
    """Create an extended ApertureSet with virtual orders beyond the edges.

    Used by Step 3 to build a better background mask: *n_below* virtual orders
    are appended below the first detected order and *n_above* above the last,
    predicted by **local spacing extrapolation** from the nearest real orders.

    Virtual orders that only partially overlap the detector are still added —
    the mask builder clips boundaries to [0, rows) so partial coverage works.

    Parameters
    ----------
    apertures : ApertureSet (possibly already gap-filled).
    image_width, image_height : detector dimensions.
    n_below : number of virtual orders to add below the first order.
    n_above : number of virtual orders to add above the last order.
    trace_degree : polynomial degree for per-order coefficients.

    Returns
    -------
    extended : new ApertureSet including virtual edge orders (marked with
               negative aperture IDs to distinguish from real ones).
    """
    # Start with a copy of all existing apertures.
    extended = ApertureSet()
    for oid in apertures.get_orders():
        extended.add_aperture(apertures.get_aperture(oid))

    real_ids = sorted(apertures.get_orders())
    min_real_id = real_ids[0]
    max_real_id = real_ids[-1]
    n_added = 0

    # Below the first order — IDs count down from min_real_id - 1.
    if n_below > 0:
        virts = _extrapolate_virtual_orders(
            apertures, image_width, image_height,
            n_orders=n_below, side='below',
            trace_degree=trace_degree)
        for k, ap in enumerate(virts, start=1):
            vid = min_real_id - k
            ap.aperture = vid
            ap.order = vid
            extended.add_aperture(ap)
            n_added += 1

    # Above the last order — IDs count up from max_real_id + 1.
    if n_above > 0:
        virts = _extrapolate_virtual_orders(
            apertures, image_width, image_height,
            n_orders=n_above, side='above',
            trace_degree=trace_degree)
        for k, ap in enumerate(virts, start=1):
            vid = max_real_id + k
            ap.aperture = vid
            ap.order = vid
            extended.add_aperture(ap)
            n_added += 1

    logger.info(f"Extended apertures for background: +{n_added} virtual orders, "
                f"total {extended.norders}")
    return extended
