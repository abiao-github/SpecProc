"""
Overscan subtraction stage of spectral reduction pipeline.

This stage should be executed FIRST, before any other processing,
as it estimates and removes CCD readout bias from overscan regions.
"""

import numpy as np
import logging
from pathlib import Path
from typing import List, Tuple, Optional
from src.utils.overscan import OverscanCorrector
from src.utils.fits_io import read_fits_image, write_fits_image
from src.config.config_manager import ConfigManager

logger = logging.getLogger(__name__)


class OverscanCorrectionStage:
    """Handles overscan correction for raw images."""

    def __init__(self, config: ConfigManager):
        """
        Initialize overscan correction stage.

        Args:
            config: Configuration manager
        """
        self.config = config
        self.corrector = self._create_corrector()

    def _create_corrector(self) -> OverscanCorrector:
        """Create OverscanCorrector with configuration."""
        # Read overscan regions from config
        overscan_config = {}

        # Check if simplified overscan_start_column is provided (1-based)
        overscan_start_col = self.config.get_int('data', 'overscan_start_column', -1)
        if overscan_start_col > 0:
            # Simplified mode: overscan is a vertical strip from start_column to image right edge
            # Convert to 0-based indexing for internal use
            overscan_config['right'] = (overscan_start_col - 1, None)  # None means to image edge
            overscan_config['left'] = None
            overscan_config['top'] = None
            overscan_config['bottom'] = None
            logger.info(f"Simplified overscan configuration: start_column={overscan_start_col} (1-based)")
        else:
            # Legacy mode: read individual regions
            for region in ['left', 'right', 'top', 'bottom']:
                start = self.config.get_int('data', f'overscan_{region}_start', -1)
                end = self.config.get_int('data', f'overscan_{region}_end', -1)

                if start >= 0 and end > start:
                    overscan_config[region] = (start, end)
                else:
                    overscan_config[region] = None

        logger.info(f"Overscan configuration: {overscan_config}")
        return OverscanCorrector(overscan_config)

    def correct_image(self, image: np.ndarray,
                     trim: bool = True) -> np.ndarray:
        """
        Correct single image for overscan.

        Args:
            image: Raw CCD image
            trim: Whether to trim overscan regions

        Returns:
            Overscan-corrected image
        """
        # Check if overscan correction is enabled
        overscan_start_col = self.config.get_int('data', 'overscan_start_column', -1)
        has_overscan = overscan_start_col > 0

        if not has_overscan:
            logger.info("No overscan region configured, skipping overscan correction")
            return image

        # Estimate overscan bias
        method = self.config.get('data', 'overscan_method', 'median')
        overscan_bias = self.corrector.estimate_overscan_bias(image, method=method)

        # Apply correction
        corrected = self.corrector.apply_overscan_correction(image, overscan_bias)

        # Trim overscan regions if overscan is configured
        if trim and has_overscan:
            corrected = self.corrector.trim_image(corrected)

        return corrected

    def correct_file(self, input_path: str, output_path: Optional[str] = None) -> str:
        """
        Correct single FITS file for overscan.

        Args:
            input_path: Path to raw FITS file
            output_path: Optional output path (uses input_path if None)

        Returns:
            Path to corrected FITS file
        """
        logger.info(f"Reading raw image: {input_path}")
        raw_image, header = read_fits_image(input_path)

        # Correct overscan
        corrected = self.correct_image(raw_image, trim=True)

        # Save corrected image
        if output_path is None:
            output_path = input_path

        write_fits_image(output_path, corrected, header=header, overwrite=True)
        logger.info(f"Saved overscan-corrected image: {output_path}")

        return output_path


def process_overscan_stage(config: ConfigManager, image_filenames: List[str],
                          midpath: str = None) -> List[str]:
    """
    Execute overscan correction stage for multiple images.

    Args:
        config: Configuration manager
        image_filenames: List of raw FITS file paths
        midpath: Intermediate output directory

    Returns:
        List of corrected image file paths
    """
    logger.info("=" * 60)
    logger.info("STAGE 0: OVERSCAN CORRECTION")
    logger.info("=" * 60)

    # Use default midpath if not provided
    if midpath is None:
        out_path = config.get('reduce', 'out_path', './output')
        midpath = str(Path(out_path) / 'midpath')

    stage = OverscanCorrectionStage(config)

    # Create output directory
    Path(midpath).mkdir(parents=True, exist_ok=True)
    overscan_path = Path(midpath) / 'overscan_corrected'
    overscan_path.mkdir(parents=True, exist_ok=True)

    corrected_files = []

    for filename in image_filenames:
        try:
            # Create output filename
            input_path = Path(filename)
            output_file = overscan_path / input_path.name

            # Correct
            output_path = stage.correct_file(str(filename), str(output_file))
            corrected_files.append(output_path)

        except Exception as e:
            logger.error(f"Error correcting {filename}: {e}")
            # Return original if correction fails
            corrected_files.append(filename)

    logger.info(f"✓ Overscan correction complete: {len(corrected_files)} images processed")

    return corrected_files
