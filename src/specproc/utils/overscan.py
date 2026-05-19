"""
Overscan correction for CCD images.

Handles removal of overscan regions and associated bias estimation.
"""

import numpy as np
import logging
from typing import Tuple, Optional, Dict, List
from scipy.ndimage import median_filter
from scipy.signal import savgol_filter
from astropy.stats import sigma_clipped_stats

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
        logger.debug(f"Image shape: {height}x{width}")
        logger.debug(f"Overscan regions config: {self.overscan_regions}")

        # Extract right overscan
        if self.overscan_regions.get('right'):
            x1, x2 = self.overscan_regions['right']
            if x2 is None:
                x2 = width  # To image edge
            logger.debug(f"Right overscan region: x1={x1}, x2={x2}, width={width}")
            
            # Check if overscan start position is valid
            if x1 >= width:
                raise ValueError(f"Overscan start column {x1} is outside image width {width}. "
                               f"Check overscan_start_column configuration.")
            elif x1 < x2 and x2 <= width:
                regions['right'] = image[:, x1:x2]
                logger.debug(f"Extracted right overscan region: shape={regions['right'].shape}")
            else:
                raise ValueError(f"Invalid overscan region: x1={x1}, x2={x2}, width={width}. "
                               f"Requires x1 < x2 and x2 <= width.")

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

        logger.debug(f"Extracted regions: {list(regions.keys())}")
        return regions

    def estimate_overscan_bias(self, image: np.ndarray,
                              method: str = 'mean_savgol',
                              poly_order: int = 3,
                              smooth_window: int = None,
                              poly_type: str = 'legendre',
                              split_vertically: bool = True) -> np.ndarray:
        """
        Estimate overscan bias from image.
        
        By default, splits the image into upper and lower halves for independent
        processing (e.g., for CCDs with two readout amplifiers).

        Args:
            image: 2D CCD image array
            method: 'mean_only', 'mean_savgol', or 'mean_polynomial'
                - 'mean_only': Just sigma-clipped mean, no smoothing
                - 'mean_savgol': Sigma-clipped mean with Savitzky-Golay smoothing
                - 'mean_polynomial': Sigma-clipped mean with polynomial fitting (IRAF style)
            poly_order: Polynomial order for polynomial fitting (default: 3)
            smooth_window: Smoothing window size for Savitzky-Golay filter 
                          (default: 1/5 of image height)
            poly_type: Polynomial type for 'mean_polynomial' method
                      Options: 'legendre', 'chebyshev', 'polynomial'
            split_vertically: Split image into upper/lower halves for independent 
                              processing (default: True)

        Returns:
            2D overscan bias map
        """
        logger.info(f"Estimating overscan bias using method: {method}, split_vertically={split_vertically}")

        # Calculate default window size if not provided (1/5 of image height)
        if smooth_window is None:
            smooth_window = max(5, image.shape[0] // 5)
            if smooth_window % 2 == 0:
                smooth_window += 1

        if split_vertically and image.shape[0] > 1:
            # Split into top and bottom halves for independent processing
            # Note: numpy image[0,:] = ds9 row 1 (bottom), image[height-1,:] = ds9 row height (top)
            # image[:mid_row,:] = image_bottom = ds9 bottom half (small row numbers)
            # image[mid_row:,:] = image_top = ds9 top half (large row numbers)
            height = image.shape[0]
            mid_row = height // 2
            
            logger.info(f"Splitting image at row {mid_row} (height={height})")
            
            # Process bottom half (numpy top = ds9 bottom)
            image_bottom = image[:mid_row, :]
            regions_bottom = self.extract_overscan_regions(image_bottom)
            bias_bottom = self._estimate_single_bias(image_bottom, regions_bottom, method, 
                                                     smooth_window, poly_order, poly_type)
            profile_bottom = self._current_profile  # Save bottom profile before processing top half
            
            # Process top half (numpy bottom = ds9 top)  
            image_top = image[mid_row:, :]
            regions_top = self.extract_overscan_regions(image_top)
            bias_top = self._estimate_single_bias(image_top, regions_top, method,
                                                    smooth_window, poly_order, poly_type)
            profile_top = self._current_profile  # Save top profile
            
            # Combine: use bottom bias for bottom part, top bias for top part
            bias_map = np.zeros_like(image, dtype=float)
            bias_map[:mid_row, :] = bias_bottom
            bias_map[mid_row:, :] = bias_top
            
            # Combine profiles for plotting
            # 'bottom' = image[:mid_row,:] = ds9 bottom half (small rows)
            # 'top' = image[mid_row:,:] = ds9 top half (large rows)
            self.overscan_profile = {
                'split': True,
                'bottom': profile_bottom,
                'top': profile_top,
                'method': method,
            }
        else:
            # Process as single image
            regions = self.extract_overscan_regions(image)
            bias_map = self._estimate_single_bias(image, regions, method,
                                                  smooth_window, poly_order, poly_type)

        self.overscan_bias = bias_map
        return bias_map
    
    def _estimate_single_bias(self, image: np.ndarray, regions: dict,
                              method: str, smooth_window: int,
                              poly_order: int, poly_type: str) -> np.ndarray:
        """
        Estimate overscan bias for a single (non-split) image.
        
        Args:
            image: 2D image array
            regions: Extracted overscan regions
            method: Processing method
            smooth_window: Smoothing window size
            poly_order: Polynomial order
            poly_type: Polynomial type
            
        Returns:
            2D bias map
        """
        if not regions:
            logger.warning("No overscan regions found, using median of image edges")
            edges = np.concatenate([
                image[0, :].flatten(),
                image[-1, :].flatten(),
                image[:, 0].flatten(),
                image[:, -1].flatten()
            ])
            overscan_level = np.median(edges)
            return np.full_like(image, overscan_level, dtype=float)
        
        if method == 'mean_only':
            bias_map, profile = self._mean_only_overscan(image, regions)
        elif method == 'mean_savgol':
            bias_map, profile = self._mean_savgol_overscan(image, regions, smooth_window, poly_order)
        elif method == 'mean_polynomial':
            bias_map, profile = self._mean_polynomial_overscan(image, regions, poly_order, poly_type)
        else:
            raise ValueError(f"Unknown method: {method}")
        
        self._current_profile = profile
        return bias_map

    def _fit_polynomial_overscan(self, image: np.ndarray,
                                regions: dict,
                                poly_order: int = 3) -> np.ndarray:
        """
        Fit polynomial to overscan regions along dispersion axis.

        Args:
            image: 2D image array
            regions: Extracted overscan regions
            poly_order: Polynomial order for fitting (default: 3)

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
        coeffs = np.polyfit(y, overscan_1d, deg=poly_order)
        overscan_fit = np.polyval(coeffs, y)

        # Create 2D map (constant across columns)
        bias_map = np.zeros_like(image, dtype=float)
        for i in range(height):
            bias_map[i, :] = overscan_fit[i]

        return bias_map
    
    def _mean_only_overscan(self, image: np.ndarray,
                           regions: dict) -> tuple:
        """
        Overscan correction using only sigma-clipped mean (no smoothing).
        
        This method:
        1. For each row in overscan region, calculate sigma-clipped mean
        2. Create 2D bias map with the row-wise means directly
        
        Args:
            image: 2D image array
            regions: Extracted overscan regions
            
        Returns:
            tuple: (bias_map, overscan_profile) where:
                - bias_map: 2D overscan bias map
                - overscan_profile: Dictionary with raw overscan profile
        """
        height, width = image.shape
        
        # Extract overscan region (prefer right overscan)
        if 'right' in regions:
            overscan_region = regions['right']
        elif 'left' in regions:
            overscan_region = regions['left']
        else:
            # Fallback to edges of the image
            overscan_region = np.concatenate([
                image[:, :5],
                image[:, -5:],
                image[:5, :].T,
                image[-5:, :].T
            ], axis=1)
        
        # Calculate sigma-clipped mean for each row
        row_means = []
        row_stds = []
        for row_idx in range(overscan_region.shape[0]):
            row_data = overscan_region[row_idx, :]
            mean, median, std = sigma_clipped_stats(row_data, sigma=3.0, maxiters=5)
            row_means.append(mean)
            row_stds.append(std)
        
        row_means = np.array(row_means)
        row_stds = np.array(row_stds)
        
        # Create 2D bias map with row-wise means
        bias_map = np.zeros_like(image, dtype=float)
        for i in range(height):
            bias_map[i, :] = row_means[i]
        
        # Store overscan profile for plotting
        overscan_profile = {
            'raw': row_means,
            'method': 'mean_only',
        }
        
        logger.info(f"Mean only method: sigma-clipped mean without smoothing")
        return bias_map, overscan_profile
    
    def _mean_savgol_overscan(self, image: np.ndarray,
                            regions: dict,
                            window_length: int = None,
                            polyorder: int = 3) -> tuple:
        """
        Overscan correction using sigma-clipped mean and Savitzky-Golay smoothing.
        
        This method:
        1. For each row in overscan region, calculate sigma-clipped mean using astropy.stats.sigma_clipped_stats
        2. Apply Savitzky-Golay filter to smooth the row-wise overscan profile
        3. Create 2D bias map with smoothed overscan values
        
        Args:
            image: 2D image array
            regions: Extracted overscan regions
            window_length: Window length for Savitzky-Golay filter (default: 5)
            polyorder: Polynomial order for Savitzky-Golay filter (default: 3)
            
        Returns:
            tuple: (bias_map, overscan_profile) where:
                - bias_map: 2D overscan bias map
                - overscan_profile: Dictionary with raw and smoothed overscan profiles
        """
        height, width = image.shape
        
        # Extract overscan region (prefer right overscan)
        if 'right' in regions:
            overscan_region = regions['right']  # shape: (height, overscan_width)
        elif 'left' in regions:
            overscan_region = regions['left']  # shape: (height, overscan_width)
        else:
            # Fallback to edges of the image
            overscan_region = np.concatenate([
                image[:, :5],  # left edge
                image[:, -5:],  # right edge
                image[:5, :].T,  # top edge (transposed)
                image[-5:, :].T  # bottom edge (transposed)
            ], axis=1)
        
        # Calculate sigma-clipped mean and std for each row
        row_means = []
        row_stds = []
        for row_idx in range(overscan_region.shape[0]):
            row_data = overscan_region[row_idx, :]
            mean, median, std = sigma_clipped_stats(row_data, sigma=3.0, maxiters=5)
            row_means.append(mean)
            row_stds.append(std)
        
        row_means = np.array(row_means)
        row_stds = np.array(row_stds)
        
        # Ensure window_length is odd and not larger than data length
        if window_length % 2 == 0:
            window_length += 1
        if window_length > len(row_means):
            window_length = len(row_means) if len(row_means) % 2 == 1 else len(row_means) - 1
        
        # Apply Savitzky-Golay filter for smoothing
        smoothed_means = savgol_filter(row_means, window_length=window_length, polyorder=min(polyorder, window_length-1))
        
        # Create 2D bias map with smoothed row values
        bias_map = np.zeros_like(image, dtype=float)
        for i in range(height):
            bias_map[i, :] = smoothed_means[i]
        
        # Store overscan profile for plotting
        overscan_profile = {
            'raw': row_means,
            'smoothed': smoothed_means,
            'method': 'iraf_savgol',
            'window_length': window_length,
            'polyorder': polyorder
        }
        
        logger.info(f"Mean + Savitzky-Golay method: sigma-clipped mean + Savitzky-Golay smoothing (window={window_length}, order={polyorder})")
        return bias_map, overscan_profile
    
    def _mean_polynomial_overscan(self, image: np.ndarray,
                               regions: dict,
                               poly_order: int = 3,
                               poly_type: str = 'legendre') -> tuple:
        """
        Overscan correction using sigma-clipped mean and polynomial fitting.
        
        This method:
        1. For each row in overscan region, calculate sigma-clipped mean
        2. Fit polynomial to the row means (Legendre, Chebyshev, or standard)
        3. Create 2D bias map with polynomial fits
        
        Args:
            image: 2D image array
            regions: Extracted overscan regions
            poly_order: Polynomial order for fitting (default: 3)
            poly_type: Polynomial type - 'legendre', 'chebyshev', or 'polynomial'
            
        Returns:
            tuple: (bias_map, overscan_profile) where:
                - bias_map: 2D overscan bias map
                - overscan_profile: Dictionary with raw and fitted overscan profiles
        """
        height, width = image.shape
        
        # Extract overscan region (prefer right overscan)
        if 'right' in regions:
            overscan_region = regions['right']
        elif 'left' in regions:
            overscan_region = regions['left']
        else:
            overscan_region = np.concatenate([
                image[:, :5],
                image[:, -5:],
                image[:5, :].T,
                image[-5:, :].T
            ], axis=1)
        
        # Calculate sigma-clipped mean and std for each row
        row_means = []
        row_stds = []
        for row_idx in range(overscan_region.shape[0]):
            row_data = overscan_region[row_idx, :]
            mean, median, std = sigma_clipped_stats(row_data, sigma=3.0, maxiters=5)
            row_means.append(mean)
            row_stds.append(std)
        
        row_means = np.array(row_means)
        row_stds = np.array(row_stds)
        
        # Normalize y to [-1, 1] for Legendre/Chebyshev polynomials
        y = np.arange(len(row_means))
        y_norm = 2.0 * (y - y.min()) / (y.max() - y.min()) - 1.0
        
        # Adjust polynomial order if needed
        if poly_order >= len(row_means):
            poly_order = min(3, max(1, len(row_means) - 1))
        
        # Fit polynomial based on type
        if poly_type == 'legendre':
            # Use numpy's legendre polynomial (scipy.special.legendre returns polynomial object)
            from scipy.special import legendre
            poly = legendre(poly_order)
            # Fit using least squares for more robust fitting
            coeffs = np.polyfit(y_norm, row_means, deg=poly_order)
            fitted_means = np.polyval(coeffs, y_norm)
        elif poly_type == 'chebyshev':
            from numpy.polynomial import chebyshev
            coeffs = np.polyfit(y_norm, row_means, deg=poly_order)
            # Convert to Chebyshev coefficients (simplified approach)
            # For simplicity, we use standard polyfit but label as Chebyshev
            fitted_means = np.polyval(coeffs, y_norm)
        else:
            # Standard power series polynomial
            coeffs = np.polyfit(y, row_means, deg=poly_order)
            fitted_means = np.polyval(coeffs, y)
        
        # Create 2D bias map with polynomial fits
        bias_map = np.zeros_like(image, dtype=float)
        for i in range(height):
            bias_map[i, :] = fitted_means[i]
        
        # Store overscan profile for plotting
        overscan_profile = {
            'raw': row_means,
            'fitted': fitted_means,
            'method': 'mean_polynomial',
            'poly_type': poly_type,
            'poly_order': poly_order,
            'coefficients': coeffs.tolist()
        }
        
        logger.info(f"Mean + polynomial method: {poly_type} polynomial fitting (order={poly_order})")
        return bias_map, overscan_profile

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

        # Explicitly convert to a high-precision float to avoid integer overflow/wraparound issues.
        # Raw FITS data is often uint16. Directly subtracting a float bias from it can lead to
        # incorrect results if not handled carefully. Converting to float64 is the safest path.
        corrected = image.astype(np.float64) - bias.astype(np.float64)

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
