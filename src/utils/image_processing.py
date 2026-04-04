"""
Image processing utilities for spectral data.

Provides functions for image operations like combining, filtering,
and detecting bad pixels.
"""

import numpy as np
from scipy import signal, ndimage
from typing import Tuple, Optional
import logging

logger = logging.getLogger(__name__)


def combine_images(images: np.ndarray, method: str = 'median',
                   sigma: float = 3.0) -> Tuple[np.ndarray, np.ndarray]:
    """
    Combine multiple 2D images.

    Args:
        images: Array of shape (n_images, height, width)
        method: 'mean', 'median', or 'sigma_clip'
        sigma: Sigma clipping threshold (for sigma_clip method)

    Returns:
        Tuple of (combined_image, uncertainty_map)
    """
    if method == 'mean':
        combined = np.mean(images, axis=0)
        uncertainty = np.std(images, axis=0) / np.sqrt(images.shape[0])
        return combined, uncertainty

    elif method == 'median':
        combined = np.median(images, axis=0)
        uncertainty = np.std(images, axis=0)
        return combined, uncertainty

    elif method == 'sigma_clip':
        combined = np.zeros_like(images[0], dtype=float)
        uncertainty = np.zeros_like(images[0], dtype=float)

        for i in range(images.shape[1]):
            for j in range(images.shape[2]):
                col = images[:, i, j]
                mean = np.mean(col)
                std = np.std(col)

                # Mask outliers
                mask = np.abs(col - mean) < sigma * std
                if np.sum(mask) > 0:
                    combined[i, j] = np.mean(col[mask])
                    uncertainty[i, j] = std
                else:
                    combined[i, j] = mean
                    uncertainty[i, j] = std

        return combined, uncertainty

    else:
        raise ValueError(f"Unknown combine method: {method}")


def savitzky_golay_2d(image: np.ndarray, window_size: int = 5,
                      order: int = 2) -> np.ndarray:
    """
    Apply 2D Savitzky-Golay filter.

    Args:
        image: 2D image array
        window_size: Filter window size (must be odd)
        order: Polynomial order

    Returns:
        Filtered image
    """
    if window_size % 2 == 0:
        window_size += 1

    # Create 1D filter
    coeffs = signal.savgol_coeffs(window_size, order)

    # Apply horizontal and vertical
    horizontal = ndimage.convolve1d(image, coeffs, axis=1)
    filtered = ndimage.convolve1d(horizontal, coeffs, axis=0)

    return filtered


def find_bad_pixels(image: np.ndarray, threshold: float = 3.0,
                    method: str = 'sigma') -> np.ndarray:
    """
    Find bad (outlier) pixels in an image.

    Args:
        image: 2D image array
        threshold: Detection threshold (sigma)
        method: 'sigma' for sigma-clipping, 'percentile' for percentile-based

    Returns:
        Boolean mask (True = bad pixel)
    """
    if method == 'sigma':
        mean = np.mean(image)
        std = np.std(image)
        bad_pixels = np.abs(image - mean) > threshold * std
        return bad_pixels

    elif method == 'percentile':
        lower = np.percentile(image, threshold)
        upper = np.percentile(image, 100 - threshold)
        bad_pixels = (image < lower) | (image > upper)
        return bad_pixels

    else:
        raise ValueError(f"Unknown method: {method}")


def estimate_background_2d(image: np.ndarray, order: int = 2) -> np.ndarray:
    """
    Estimate smooth 2D background using polynomial fitting.

    Args:
        image: 2D image array
        order: Polynomial order

    Returns:
        Background model
    """
    rows, cols = image.shape
    x = np.arange(cols)
    y = np.arange(rows)

    # Sample grid
    yy, xx = np.meshgrid(y, x, indexing='ij')

    # Flatten
    xx_flat = xx.flatten()
    yy_flat = yy.flatten()
    z_flat = image.flatten()

    # Build polynomial matrix
    ncoeff = (order + 1) * (order + 1)
    A = np.zeros((len(z_flat), ncoeff))
    idx = 0
    for i in range(order + 1):
        for j in range(order + 1):
            A[:, idx] = (xx_flat ** i) * (yy_flat ** j)
            idx += 1

    # Solve
    try:
        coeffs, _, _, _ = np.linalg.lstsq(A, z_flat, rcond=None)
        background = np.zeros_like(image)
        idx = 0
        for i in range(order + 1):
            for j in range(order + 1):
                background += coeffs[idx] * (xx ** i) * (yy ** j)
                idx += 1
        return background
    except Exception as e:
        logger.error(f"Error estimating background: {e}")
        return np.zeros_like(image)


def normalize_flat(flat_image: np.ndarray, mask: Optional[np.ndarray] = None) -> np.ndarray:
    """
    Normalize flat field to unity.

    Args:
        flat_image: Flat field image
        mask: Optional bad pixel mask

    Returns:
        Normalized flat (clipped to valid range)
    """
    if mask is not None:
        flat_array = flat_image.copy()
        flat_array[mask > 0] = np.nan
        mean_val = np.nanmean(flat_array)
    else:
        mean_val = np.mean(flat_image)

    normalized = flat_image / mean_val

    # Clip to avoid extreme values
    normalized = np.clip(normalized, 0.1, 10.0)

    return normalized


def extract_subimage(image: np.ndarray, x_start: int, x_end: int,
                    y_start: int, y_end: int) -> np.ndarray:
    """Extract rectangular sub-image."""
    return image[y_start:y_end, x_start:x_end]


def trim_image(image: np.ndarray, trim_left: int = 0, trim_right: int = 0,
              trim_top: int = 0, trim_bottom: int = 0) -> np.ndarray:
    """Trim edges of image."""
    rows, cols = image.shape
    return image[trim_top:rows-trim_bottom, trim_left:cols-trim_right]
