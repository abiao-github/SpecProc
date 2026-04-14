"""
Image processing utilities for spectral data.

Provides functions for image operations like combining, filtering,
and detecting bad pixels.
"""

import numpy as np
from typing import Optional, Tuple
from scipy.ndimage import gaussian_filter
from scipy import signal, ndimage
from scipy.interpolate import SmoothBivariateSpline
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


def _build_design_matrix(x: np.ndarray, y: np.ndarray, order: int,
                         method: str) -> np.ndarray:
    """Build 2D basis matrix for polynomial or Chebyshev fitting."""
    if method == 'chebyshev':
        vander = np.polynomial.chebyshev.chebvander2d(x, y, [order, order])
        return vander.reshape(x.size, -1)

    # Fallback: regular polynomial basis.
    ncoeff = (order + 1) * (order + 1)
    A = np.zeros((x.size, ncoeff), dtype=np.float64)
    idx = 0
    for i in range(order + 1):
        for j in range(order + 1):
            A[:, idx] = (x ** i) * (y ** j)
            idx += 1
    return A


def _robust_sigma(values: np.ndarray) -> float:
    """Robust sigma estimate via MAD; returns 0 for empty arrays."""
    if values.size == 0:
        return 0.0
    med = np.median(values)
    mad = np.median(np.abs(values - med))
    if mad <= 0:
        std = np.std(values)
        return float(std) if np.isfinite(std) else 0.0
    return float(1.4826 * mad)


def _estimate_single_background(image: np.ndarray, order: int = 2,
                                valid_mask: Optional[np.ndarray] = None,
                                method: str = 'chebyshev',
                                smooth_sigma: float = 20.0,
                                sigma_clip: float = 3.0,
                                maxiters: int = 4,
<<<<<<< HEAD
                                bspline_smooth: float = 1.0) -> np.ndarray:
    """Estimate smooth 2D background for a single (non-split) image."""
=======
                                bspline_smooth: float = 1.0,
                                clip_mode: str = 'both') -> np.ndarray:
    """Estimate smooth 2D background for a single (non-split) image.

    Parameters
    ----------
    clip_mode : {'both', 'upper', 'lower'}
        Sigma-clipping mode.  'both' rejects outliers on both sides
        (symmetric, the traditional behaviour).  'upper' rejects only
        positive outliers (pixels brighter than the model) which is
        recommended for flat-field background estimation where unmasked
        faint orders contaminate the inter-order region.  'lower' rejects
        only negative outliers.
    """
>>>>>>> cef6f04 (	modified:   README.md)
    rows, cols = image.shape
    yy, xx = np.meshgrid(np.arange(rows), np.arange(cols), indexing='ij')

    fit_mask = np.isfinite(image)
    if valid_mask is not None:
        fit_mask &= valid_mask.astype(bool)

    if np.sum(fit_mask) < max(100, (order + 1) ** 2):
        logger.warning("Too few valid samples for background fit; returning zeros")
        return np.zeros_like(image, dtype=np.float64)

    work_data = image.astype(np.float64)

<<<<<<< HEAD
    if method == 'smooth':
=======
    if method in ('smooth', 'gaussian_smooth'):
>>>>>>> cef6f04 (	modified:   README.md)
        current = fit_mask.copy()
        s = max(1.0, float(smooth_sigma))

        for _ in range(max(1, int(maxiters))):
            use = np.zeros_like(work_data, dtype=np.float64)
            wgt = np.zeros_like(work_data, dtype=np.float64)
            use[current] = work_data[current]
            wgt[current] = 1.0

            num = gaussian_filter(use, sigma=s, mode='nearest')
            den = gaussian_filter(wgt, sigma=s, mode='nearest')
            model = num / np.clip(den, 1e-8, None)

            if sigma_clip <= 0:
                return model

            resid = work_data[current] - model[current]
            sig = _robust_sigma(resid)
            if sig <= 0:
                break

<<<<<<< HEAD
            keep = np.abs(work_data - model) <= sigma_clip * sig
=======
            diff = work_data - model
            if clip_mode == 'upper':
                keep = diff <= sigma_clip * sig
            elif clip_mode == 'lower':
                keep = diff >= -sigma_clip * sig
            else:
                keep = np.abs(diff) <= sigma_clip * sig
>>>>>>> cef6f04 (	modified:   README.md)
            new_current = current & keep
            if np.array_equal(new_current, current):
                break
            if np.sum(new_current) < max(100, (order + 1) ** 2):
                break
            current = new_current

        return model

    x_norm = (xx.astype(np.float64) / max(cols - 1, 1)) * 2.0 - 1.0
    y_norm = (yy.astype(np.float64) / max(rows - 1, 1)) * 2.0 - 1.0

    xp = x_norm[fit_mask]
    yp = y_norm[fit_mask]
    zp = work_data[fit_mask]
    keep = np.ones_like(zp, dtype=bool)

    for _ in range(max(1, int(maxiters))):
        if np.sum(keep) < max(100, (order + 1) ** 2):
            break

        try:
            if method == 'bspline':
                kxy = max(1, min(3, int(order)))
                s = max(0.0, float(bspline_smooth)) * float(np.sum(keep))
                spline = SmoothBivariateSpline(
                    xp[keep], yp[keep], zp[keep], kx=kxy, ky=kxy, s=s
                )
                fit_all = spline.ev(xp, yp)
                model = spline.ev(x_norm.ravel(), y_norm.ravel()).reshape(rows, cols)
            else:
                fit_kind = 'chebyshev' if method not in ('poly', 'polynomial') else 'poly'
                A = _build_design_matrix(xp[keep], yp[keep], order, fit_kind)
                coeffs, _, _, _ = np.linalg.lstsq(A, zp[keep], rcond=None)
                A_all = _build_design_matrix(xp, yp, order, fit_kind)
                fit_all = A_all @ coeffs
                A_img = _build_design_matrix(x_norm.ravel(), y_norm.ravel(), order, fit_kind)
                model = (A_img @ coeffs).reshape(rows, cols)
        except Exception as e:
<<<<<<< HEAD
            logger.warning(f"Background {method} fit failed, fallback to smooth: {e}")
=======
            logger.warning(f"Background {method} fit failed, fallback to gaussian_smooth: {e}")
>>>>>>> cef6f04 (	modified:   README.md)
            return _estimate_single_background(
                image,
                order=order,
                valid_mask=valid_mask,
<<<<<<< HEAD
                method='smooth',
=======
                method='gaussian_smooth',
>>>>>>> cef6f04 (	modified:   README.md)
                smooth_sigma=smooth_sigma,
                sigma_clip=sigma_clip,
                maxiters=maxiters,
                bspline_smooth=bspline_smooth,
<<<<<<< HEAD
=======
                clip_mode=clip_mode,
>>>>>>> cef6f04 (	modified:   README.md)
            )

        if sigma_clip <= 0:
            return model

        resid = zp - fit_all
        sig = _robust_sigma(resid[keep])
        if sig <= 0:
            break

        center = np.median(resid[keep])
<<<<<<< HEAD
        new_keep = np.abs(resid - center) <= sigma_clip * sig
=======
        if clip_mode == 'upper':
            new_keep = (resid - center) <= sigma_clip * sig
        elif clip_mode == 'lower':
            new_keep = (resid - center) >= -sigma_clip * sig
        else:
            new_keep = np.abs(resid - center) <= sigma_clip * sig
>>>>>>> cef6f04 (	modified:   README.md)
        if np.array_equal(new_keep, keep):
            break
        keep = new_keep

    return model


<<<<<<< HEAD
=======
# Public alias: estimate a single region without any internal splitting.
# Use this when the caller handles detector-half splitting externally.
estimate_background_region = _estimate_single_background


>>>>>>> cef6f04 (	modified:   README.md)
def estimate_background_2d(image: np.ndarray, order: int = 2,
                           split_vertically: bool = True,
                           valid_mask: Optional[np.ndarray] = None,
                           method: str = 'chebyshev',
                           smooth_sigma: float = 20.0,
                           sigma_clip: float = 3.0,
                           maxiters: int = 4,
                           bspline_smooth: float = 1.0) -> np.ndarray:
    """
    Estimate smooth 2D background using polynomial fitting.
    
    By default, splits the image into upper and lower halves for independent
    processing (e.g., for CCDs with two readout amplifiers).

    Args:
        image: 2D image array
        order: Polynomial order
        split_vertically: Split image into upper/lower halves for independent 
                         background fitting (default: True)
        valid_mask: Optional boolean mask of valid sampling pixels (True = use)
<<<<<<< HEAD
        method: 'chebyshev', 'bspline', 'smooth', or 'poly'
=======
        method: 'chebyshev', 'bspline', 'gaussian_smooth', or 'poly'
>>>>>>> cef6f04 (	modified:   README.md)
        smooth_sigma: Gaussian sigma for smooth method
        sigma_clip: Sigma-clipping threshold for robust rejection
        maxiters: Maximum sigma-clipping iterations
        bspline_smooth: Smoothing factor scale for bspline method

    Returns:
        Background model
    """
    rows, cols = image.shape
    
    if split_vertically and rows > 1:
        # Split into top and bottom halves for independent processing
        # Note: numpy image[0,:] = ds9 row 1 (bottom), image[rows-1,:] = ds9 row rows (top)
        # image[:mid_row,:] = image_bottom = ds9 bottom half (small row numbers)
        # image[mid_row:,:] = image_top = ds9 top half (large row numbers)
        mid_row = rows // 2
        logger.info(f"Estimating background with vertical split at row {mid_row}")
        
        image_bottom = image[:mid_row, :]
        image_top = image[mid_row:, :]
        mask_bottom = valid_mask[:mid_row, :] if valid_mask is not None else None
        mask_top = valid_mask[mid_row:, :] if valid_mask is not None else None
        
        background_bottom = _estimate_single_background(
            image_bottom,
            order,
            valid_mask=mask_bottom,
            method=method,
            smooth_sigma=smooth_sigma,
            sigma_clip=sigma_clip,
            maxiters=maxiters,
            bspline_smooth=bspline_smooth,
        )
        background_top = _estimate_single_background(
            image_top,
            order,
            valid_mask=mask_top,
            method=method,
            smooth_sigma=smooth_sigma,
            sigma_clip=sigma_clip,
            maxiters=maxiters,
            bspline_smooth=bspline_smooth,
        )
        
        # Combine backgrounds
        background = np.zeros_like(image)
        background[:mid_row, :] = background_bottom
        background[mid_row:, :] = background_top
        
        return background
    else:
        return _estimate_single_background(
            image,
            order,
            valid_mask=valid_mask,
            method=method,
            smooth_sigma=smooth_sigma,
            sigma_clip=sigma_clip,
            maxiters=maxiters,
            bspline_smooth=bspline_smooth,
        )


def normalize_flat(flat_image: np.ndarray, mask: Optional[np.ndarray] = None) -> np.ndarray:
    """
    Normalize flat field to unity using sigma clipping.

    Args:
        flat_image: Flat field image
        mask: Optional bad pixel mask

    Returns:
        Normalized flat (clipped to valid range)
    """
    from astropy.stats import sigma_clipped_stats

    # Apply mask if provided
    if mask is not None:
        flat_array = flat_image.copy()
        flat_array[mask > 0] = np.nan
        # Use sigma clipping with masked values
        mean_val, _, _ = sigma_clipped_stats(flat_array, mask=np.isnan(flat_array), sigma=3.0, maxiters=5)
    else:
        # Use sigma clipping without mask
        mean_val, _, _ = sigma_clipped_stats(flat_image, sigma=3.0, maxiters=5)

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
