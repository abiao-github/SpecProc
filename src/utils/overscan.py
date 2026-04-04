"""
Overscan correction for CCD images.

Handles removal of overscan regions and associated bias estimation.
"""

import numpy as np
import logging
from typing import Tuple, Optional
from scipy.ndimage import median_filter

logger = logging.getLogger(__name__)


class OverscanCorrector:
    """Handles overscan bias correction."""

    def __init__(self, overscan_regions: Optional[dict] = None):
        """
        Initialize overscan corrector.

        Args:
            overscan_regions: Dict with 'left', 'right', 'top', 'bottom' pixel ranges
                Example: {'left': (0, 50), 'right': (4110, 4160),
                         'top': (0, 10), 'bottom': (4120, 4136)}
        """
        # Default for generic 4160x4136 CCD
        self.overscan_regions = overscan_regions or {
            'left': None,      # No left overscan
            'right': (4110, 4160),   # Right overscan region
            'top': None,       # No top overscan
            'bottom': None,    # No bottom overscan
        }
        self.overscan_bias = None

    def extract_overscan_regions(self, image: np.ndarray) -> dict:
        """
        Extract overscan regions from image.

        Args:
            image: 2D CCD image array

        Returns:
            Dictionary with extracted overscan arrays
        """
        regions = {}
        height, width = image.shape

        # Extract right overscan
        if self.overscan_regions.get('right'):
            x1, x2 = self.overscan_regions['right']
            if x2 is None:
                x2 = width  # To image edge
            if x1 < x2 and x2 <= width:
                regions['right'] = image[:, x1:x2]

        # Extract left overscan
        if self.overscan_regions.get('left'):
            x1, x2 = self.overscan_regions['left']
            if x2 is None:
                x2 = width
            if x1 < x2 and x2 <= width:
                regions['left'] = image[:, x1:x2]

        # Extract top overscan
        if self.overscan_regions.get('top'):
            y1, y2 = self.overscan_regions['top']
            if y2 is None:
                y2 = height
            if y1 < y2 and y2 <= height:
                regions['top'] = image[y1:y2, :]

        # Extract bottom overscan
        if self.overscan_regions.get('bottom'):
            y1, y2 = self.overscan_regions['bottom']
            if y2 is None:
                y2 = height
            if y1 < y2 and y2 <= height:
                regions['bottom'] = image[y1:y2, :]

        return regions

    def estimate_overscan_bias(self, image: np.ndarray,
                              method: str = 'median') -> np.ndarray:
        """
        Estimate overscan bias from image.

        Args:
            image: 2D CCD image array
            method: 'median' or 'polynomial'

        Returns:
            2D overscan bias map
        """
        logger.info("Estimating overscan bias...")

        regions = self.extract_overscan_regions(image)

        if not regions:
            logger.warning("No overscan regions found, using median of image edges")
            # Fallback: use image edges
            edges = np.concatenate([
                image[0, :].flatten(),
                image[-1, :].flatten(),
                image[:, 0].flatten(),
                image[:, -1].flatten()
            ])
            overscan_level = np.median(edges)
            bias_map = np.full_like(image, overscan_level, dtype=float)
            return bias_map

        if method == 'median':
            # Use median of overscan regions
            overscan_values = []
            for region in regions.values():
                overscan_values.extend(region.flatten())

            overscan_level = np.median(np.array(overscan_values))
            logger.info(f"Median overscan level: {overscan_level:.1f} counts")

            # Create 2D bias map (flat)
            bias_map = np.full_like(image, overscan_level, dtype=float)

        elif method == 'polynomial':
            # Fit polynomial to overscan columns/rows
            bias_map = self._fit_polynomial_overscan(image, regions)

        else:
            raise ValueError(f"Unknown method: {method}")

        self.overscan_bias = bias_map
        return bias_map

    def _fit_polynomial_overscan(self, image: np.ndarray,
                                regions: dict) -> np.ndarray:
        """
        Fit polynomial to overscan regions along dispersion axis.

        Args:
            image: 2D image array
            regions: Extracted overscan regions

        Returns:
            2D polynomial overscan bias map
        """
        height, width = image.shape

        # Use right overscan if available
        if 'right' in regions:
            overscan_cols = regions['right']
            # Median across overscan columns, then fit along rows
            overscan_1d = np.median(overscan_cols, axis=1)
        elif 'left' in regions:
            overscan_cols = regions['left']
            overscan_1d = np.median(overscan_cols, axis=1)
        else:
            # Fallback
            overscan_1d = np.median(image, axis=1)

        # Fit polynomial to overscan variation
        y = np.arange(len(overscan_1d))
        coeffs = np.polyfit(y, overscan_1d, deg=3)
        overscan_fit = np.polyval(coeffs, y)

        # Create 2D map (constant across columns)
        bias_map = np.zeros_like(image, dtype=float)
        for i in range(height):
            bias_map[i, :] = overscan_fit[i]

        return bias_map

    def apply_overscan_correction(self, image: np.ndarray,
                                 bias_map: Optional[np.ndarray] = None) -> np.ndarray:
        """
        Apply overscan correction to image.

        Args:
            image: Raw CCD image
            bias_map: Optional precomputed overscan bias map

        Returns:
            Overscan-corrected image
        """
        bias = bias_map if bias_map is not None else self.overscan_bias

        if bias is None:
            logger.warning("No overscan bias available, estimating...")
            bias = self.estimate_overscan_bias(image)

        if image.shape != bias.shape:
            raise ValueError(f"Image shape {image.shape} != bias shape {bias.shape}")

        corrected = image.astype(float) - bias
        logger.debug(f"Overscan correction applied: "
                    f"image range [{np.min(corrected):.1f}, {np.max(corrected):.1f}]")

        return corrected

    def trim_image(self, image: np.ndarray) -> np.ndarray:
        """
        Trim overscan regions from image.

        Args:
            image: Image array (may include overscan)

        Returns:
            Trimmed image with overscan removed
        """
        height, width = image.shape

        # Define science region (exclude overscan)
        x1, x2 = 0, width
        y1, y2 = 0, height

        # Trim if overscan regions are defined
        if self.overscan_regions.get('right'):
            x_overscan_start = self.overscan_regions['right'][0]
            x2 = min(x2, x_overscan_start)

        if self.overscan_regions.get('left'):
            x_overscan_end = self.overscan_regions['left'][1]
            x1 = max(x1, x_overscan_end)

        if self.overscan_regions.get('bottom'):
            y_overscan_start = self.overscan_regions['bottom'][0]
            y2 = min(y2, y_overscan_start)

        if self.overscan_regions.get('top'):
            y_overscan_end = self.overscan_regions['top'][1]
            y1 = max(y1, y_overscan_end)

        trimmed = image[y1:y2, x1:x2]
        logger.info(f"Image trimmed: {image.shape} → {trimmed.shape}")

        return trimmed


def process_overscan_correction(image: np.ndarray,
                               overscan_config: Optional[dict] = None,
                               trim: bool = True) -> np.ndarray:
    """
    Convenience function for overscan correction.

    Args:
        image: Raw CCD image
        overscan_config: Overscan region configuration
        trim: Whether to trim overscan regions after correction

    Returns:
        Overscan-corrected (and optionally trimmed) image
    """
    corrector = OverscanCorrector(overscan_config)

    # Correct overscan
    corrected = corrector.apply_overscan_correction(image)

    # Optionally trim
    if trim:
        corrected = corrector.trim_image(corrected)

    return corrected
