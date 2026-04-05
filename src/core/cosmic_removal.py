"""
Cosmic ray detection and correction stage.
"""

import numpy as np
import logging
from pathlib import Path
from typing import List, Optional
from src.utils.fits_io import read_fits_image, write_fits_image
from src.config.config_manager import ConfigManager

logger = logging.getLogger(__name__)


def _detect_cosmics(image: np.ndarray, sigma: float = 5.0, window: int = 7) -> np.ndarray:
    """Detect cosmic rays using local median + sigma threshold."""
    from scipy.ndimage import median_filter

    # local median and RMS estimate
    med = median_filter(image, size=window)
    diff = image - med

    # robust sigma estimate
    mad = np.median(np.abs(diff - np.median(diff)))
    std_est = 1.4826 * mad if mad > 0 else np.std(diff)

    threshold = sigma * max(std_est, 1e-6)
    cosmic_mask = diff > threshold

    return cosmic_mask


def _fix_cosmics(image: np.ndarray, cosmic_mask: np.ndarray, window: int = 7) -> np.ndarray:
    """Replace cosmic pixels with local median."""
    from scipy.ndimage import median_filter

    corrected = image.copy()
    local_med = median_filter(image, size=window)
    corrected[cosmic_mask] = local_med[cosmic_mask]

    return corrected


def process_cosmic_stage(config: ConfigManager, image_filenames: List[str], midpath: Optional[str] = None) -> List[str]:
    """
    Execute cosmic ray correction stage for a list of images.

    Args:
        config: Configuration manager
        image_filenames: List of FITS file paths
        midpath: Intermediate output directory

    Returns:
        List of corrected image file paths
    """
    logger.info("=" * 60)
    logger.info("STAGE 0.5: COSMIC RAY REMOVAL")
    logger.info("=" * 60)

    if midpath is None:
        out_path = config.get('reduce', 'out_path', './output')

    cosmic_path = Path(out_path) / 'step4_cosmic'
    cosmic_path.mkdir(parents=True, exist_ok=True)

    sigma = config.get_float('reduce', 'cosmic_sigma', 5.0)
    window = config.get_int('reduce', 'cosmic_window', 7)

    corrected_files = []

    for filename in image_filenames:
        try:
            input_path = Path(filename)
            output_file = cosmic_path / input_path.name

            image, header = read_fits_image(str(filename))

            mask = _detect_cosmics(image, sigma=sigma, window=window)
            corrected = _fix_cosmics(image, mask, window=window)

            write_fits_image(str(output_file), corrected, header=header, overwrite=True)

            corrected_files.append(str(output_file))
            logger.info(f"Cosmic correction: {filename} -> {output_file} ({np.sum(mask)} pixels fixed)")

        except Exception as e:
            logger.error(f"Error cosmic-correcting {filename}: {e}")
            corrected_files.append(str(filename))

    logger.info(f"✓ Cosmic correction complete: {len(corrected_files)} images processed")
    return corrected_files
