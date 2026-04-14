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
from src.plotting.spectra_plotter import plot_2d_image_to_file, plot_overscan_corrected_image

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
        logger.info(f"Overscan start column: {overscan_start_col}, has_overscan: {has_overscan}")

        if not has_overscan:
            logger.info("No overscan region configured, skipping overscan correction")
            return image

        # Estimate overscan bias
        method = self.config.get('data', 'overscan_method', 'mean_only')
        smooth_window = self.config.get_int('data', 'overscan_smooth_window', -1)
        poly_order = self.config.get_int('data', 'overscan_poly_order', 3)
        poly_type = self.config.get('data', 'overscan_poly_type', 'legendre')
        logger.info(f"Overscan correction method: {method}, smooth_window: {smooth_window}, poly_type: {poly_type}, poly_order: {poly_order}")
        overscan_bias = self.corrector.estimate_overscan_bias(
            image, method=method, 
            smooth_window=smooth_window if smooth_window > 0 else None,
            poly_order=poly_order,
            poly_type=poly_type
        )

        # Apply correction
        corrected = self.corrector.apply_overscan_correction(image, overscan_bias)

        # Trim overscan regions if overscan is configured
        if trim and has_overscan:
            corrected = self.corrector.trim_image(corrected)

        # Optional user crop (DS9 indexing, 1-based; inclusive bounds).
        corrected = self._apply_user_crop(corrected)

        return corrected

    def _apply_user_crop(self, image: np.ndarray) -> np.ndarray:
        """Apply optional crop bounds from config using DS9 1-based indexing."""
        rows, cols = image.shape

        x1 = self.config.get_int('data', 'trim_x_start', -1)
        x2 = self.config.get_int('data', 'trim_x_end', -1)
        y1 = self.config.get_int('data', 'trim_y_start', -1)
        y2 = self.config.get_int('data', 'trim_y_end', -1)

        # Default: disabled crop.
        if min(x1, x2, y1, y2) < 1:
            return image

        # Convert DS9 1-based inclusive indices to Python 0-based half-open slices.
        x1p = max(0, x1 - 1)
        x2p = min(cols, x2)
        y1p = max(0, y1 - 1)
        y2p = min(rows, y2)

        if x2p <= x1p or y2p <= y1p:
            logger.warning(
                f"Invalid crop range (DS9): x=[{x1},{x2}], y=[{y1},{y2}]. Skip cropping."
            )
            return image

        cropped = image[y1p:y2p, x1p:x2p]
        logger.info(
            f"Applied user crop (DS9): x=[{x1},{x2}], y=[{y1},{y2}] -> shape {cropped.shape}"
        )
        return cropped

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

        # Add processing marker to header (only in output file)
        header['OVRSCAN'] = (True, 'Overscan correction completed')

        # Save corrected image
        if output_path is None:
            output_path = input_path

        write_fits_image(output_path, corrected, header=header, overwrite=True, dtype='float32')
        logger.info(f"Saved overscan-corrected image: {output_path}")

        # Save diagnostic plots if enabled
        save_plots = self.config.get_bool('reduce', 'save_plots', True)
        if save_plots and output_path is not None and output_path != input_path:
            out_dir = Path(output_path).parent
            stem = Path(input_path).stem
            fig_format = self.config.get('reduce', 'fig_format', 'png')

            # Plot 0: Raw input image for before/after comparison.
            raw_plot_file = out_dir / f'{stem}_raw.{fig_format}'
            logger.debug(f"Saving raw image plot to: {raw_plot_file}")
            plot_2d_image_to_file(raw_image, str(raw_plot_file), "Raw Input Image")
            
            # Plot 1: Overscan-corrected image (no overscan region marking)
            plot_file = out_dir / f'{stem}_overscan_corrected.{fig_format}'
            logger.debug(f"Saving overscan-corrected image plot to: {plot_file}")
            # Get overscan region info for plotting
            overscan_start_col = self.config.get_int('data', 'overscan_start_column', -1)
            has_overscan = overscan_start_col > 0
            if has_overscan:
                # Get original image width from the raw image
                original_width = raw_image.shape[1]
                
                logger.debug(f"Overscan region: columns {overscan_start_col}-{original_width} (1-based)")
                # Plot corrected image (no overscan region marking)
                plot_overscan_corrected_image(
                    corrected, str(plot_file), "Overscan Corrected Image",
                    overscan_start_col=overscan_start_col,
                    original_width=original_width
                )
            else:
                plot_2d_image_to_file(corrected, str(plot_file), "Overscan Corrected Image")
            
            # Plot 2: Overscan profile (if using methods that generate profile)
            method = self.config.get('data', 'overscan_method', 'mean_savgol')
            if hasattr(self.corrector, 'overscan_profile') and self.corrector.overscan_profile:
                profile_file = out_dir / f'{stem}_overscan_profile.{fig_format}'
                logger.debug(f"Saving overscan profile plot to: {profile_file}")
                from src.plotting.spectra_plotter import plot_overscan_profile
                plot_overscan_profile(
                    self.corrector.overscan_profile, 
                    str(profile_file),
                    title=f"Overscan Profile - {stem}"
                )

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
        output_dir = config.get_output_path()

    stage = OverscanCorrectionStage(config)

    # Create output directory under unified Step1 basic preprocessing layout.
    overscan_path = Path(output_dir) / 'step1_basic' / 'overscan_corrected'
    overscan_path.mkdir(parents=True, exist_ok=True)

    corrected_files = []

    for filename in image_filenames:
        try:
            # Create output filename
            input_path = Path(filename)
            output_file = overscan_path / input_path.name

            # Correct
            corrected_file = stage.correct_file(str(filename), str(output_file))
            corrected_files.append(corrected_file)

        except Exception as e:
            logger.error(f"Error correcting {filename}: {e}")
            # Return original if correction fails
            corrected_files.append(filename)

    logger.info(f"✓ Overscan correction complete: {len(corrected_files)} images processed")

    return corrected_files
