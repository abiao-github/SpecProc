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
from src.utils.fits_io import read_fits_image, write_fits_image
from src.utils.image_processing import combine_images, normalize_flat, find_bad_pixels
from src.core.data_structures import ApertureSet, ApertureLocation, FlatField
from src.config.config_manager import ConfigManager

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
        for filename in flat_filenames:
            try:
                img, _ = read_fits_image(filename)
                flat_images.append(img)
            except Exception as e:
                logger.warning(f"Failed to read flat frame {filename}: {e}")
                continue

        if not flat_images:
            raise RuntimeError("Could not read any flat frames")

        # Convert to array
        images_array = np.array(flat_images, dtype=np.float32)
        logger.info(f"Successfully read {len(images_array)} flat frames")

        # Combine
        method = self.config.get('reduce.flat', 'combine_method', 'median')
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

    def detect_orders(self, threshold: float = 0.3, trace_degree: int = 3) -> ApertureSet:
        """
        Detect echelle orders in flat field with robust tracing.

        Uses collapsed profile peak finding followed by center tracing across image rows
        and polynomial fitting for curved orders.

        Args:
            threshold: Detection threshold (relative to max on collapsed profile)
            trace_degree: Polynomial degree used to fit order trace.

        Returns:
            ApertureSet with detected orders
        """
        if self.flat_norm is None:
            self.normalize_flat()

        logger.info("Detecting echelle orders (robust tracing)...")

        image = self.flat_norm
        height, width = image.shape

        # Collapse spatial direction (Y) to find initial order centroids along X
        profile = np.mean(image, axis=0)
        peak_threshold = threshold * np.max(profile)

        # Find candidate order peaks
        peaks = self._find_peaks(profile,
                                 min_height=peak_threshold,
                                 distance=30,
                                 prominence=0.05 * np.max(profile))

        if not peaks:
            logger.warning("No order peaks detected; lower threshold or check flat exposures")
            return ApertureSet()

        apertures = ApertureSet()
        order_num = 1

        for x_peak in peaks:
            # Trace each order in Y direction by finding local maximum near x_peak each row
            y_positions = []
            x_positions = []

            for y in range(0, height, 4):
                x_min = int(max(0, x_peak - 20))
                x_max = int(min(width, x_peak + 20))
                if x_max <= x_min:
                    continue
                strip = image[y, x_min:x_max]
                if strip.size == 0:
                    continue

                local_index = np.argmax(strip)
                x_center = x_min + local_index

                # require a minimum local flux to ignore noise
                if strip[local_index] < 0.2 * np.max(profile):
                    continue

                y_positions.append(y)
                x_positions.append(x_center)

            if len(y_positions) < 8:
                # Too few points to define stable trace
                continue

            # Fit polynomial to the traced centers
            coeff = np.polyfit(y_positions, x_positions, trace_degree)
            center_coef = np.array(coeff)

            # Determine width from local profile (approx median FWHM)
            profile_width = self.config.get_float('reduce.trace', 'separation', 30.0)
            half_width = profile_width / 2.0

            lower_coef = center_coef.copy()
            lower_coef[-1] -= half_width
            upper_coef = center_coef.copy()
            upper_coef[-1] += half_width

            loc = ApertureLocation(
                aperture=order_num,
                order=order_num,
                center_coef=center_coef,
                lower_coef=lower_coef,
                upper_coef=upper_coef,
                width=profile_width
            )
            apertures.add_aperture(loc)
            order_num += 1

        logger.info(f"Detected {apertures.norders} orders (optimized)")
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
            # Get order trace positions
            y_coords = np.arange(height)
            center_pos = aperture.get_position(y_coords)

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

    def save_flat_field(self, output_path: str):
        """Save flat field data to FITS file."""
        if self.flat_data is None:
            raise RuntimeError("No flat data to save")

        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        from astropy.io import fits

        # Create multi-HDU FITS
        hdul = fits.HDUList([
            fits.PrimaryHDU(data=self.flat_data.astype(np.float32)),
            fits.ImageHDU(data=self.flat_mask, name='MASK'),
            fits.ImageHDU(data=self.flat_norm.astype(np.float32), name='NORM'),
        ])

        if self.flat_sens is not None:
            hdul.append(fits.ImageHDU(data=self.flat_sens, name='SENSITIVITY'))

        hdul.writeto(str(output_path), overwrite=True)
        logger.info(f"Saved flat field to {output_path}")


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
    flat_norm = processor.normalize_flat()

    # Extract sensitivity
    flat_sens = processor.extract_sensitivity_map()

    # Detect orders
    apertures = processor.detect_orders(threshold=0.3)

    # Extract blaze profiles
    blaze_profiles = processor.extract_blaze_profiles(apertures)

    # Create output directory
    Path(midpath).mkdir(parents=True, exist_ok=True)

    # Save flat field
    flat_file = Path(midpath) / 'flat.fits'
    processor.save_flat_field(str(flat_file))

    # Create FlatField object
    flat_field = FlatField(
        flat_data=flat_data,
        flat_sens=flat_sens,
        flat_norm=flat_norm,
        flat_mask=flat_mask,
        aperture_set=apertures,
        blaze_profiles=blaze_profiles,
        nrows=flat_data.shape[0],
        ncols=flat_data.shape[1]
    )

    return flat_field, apertures
