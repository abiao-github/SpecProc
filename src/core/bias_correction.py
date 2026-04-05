"""
Bias subtraction stage of spectral reduction pipeline.

Handles bias frame combination and application to science images.
"""

import numpy as np
import logging
from pathlib import Path
from typing import List, Tuple, Optional
from src.utils.fits_io import read_fits_image, write_fits_image, combine_fits_images
from src.utils.image_processing import combine_images
from src.config.config_manager import ConfigManager
from src.plotting.spectra_plotter import plot_2d_image_to_file

logger = logging.getLogger(__name__)


class BiasCorrector:
    """Handles bias frame combination and correction."""

    def __init__(self, config: ConfigManager):
        """
        Initialize bias corrector.

        Args:
            config: Configuration manager
        """
        self.config = config
        self.master_bias = None
        self.bias_uncertainty = None

    def combine_bias_frames(self, bias_filenames: List[str]) -> Tuple[np.ndarray, np.ndarray]:
        """
        Combine multiple bias frames to create master bias.

        Args:
            bias_filenames: List of bias FITS file paths

        Returns:
            Tuple of (master_bias, uncertainty)
        """
        logger.info(f"Combining {len(bias_filenames)} bias frames...")

        if not bias_filenames:
            logger.warning("No bias frames provided")
            return None, None

        # Read all bias images
        bias_images = []
        for filename in bias_filenames:
            try:
                img, _ = read_fits_image(filename)
                bias_images.append(img)
            except Exception as e:
                logger.warning(f"Failed to read bias frame {filename}: {e}")
                continue

        if not bias_images:
            raise RuntimeError("Could not read any bias frames")

        # Convert to array
        images_array = np.array(bias_images)
        logger.info(f"Successfully read {len(images_array)} bias frames")

        # Combine using configured method
        method = self.config.get('reduce.bias', 'combine_method', 'median')
        sigma = self.config.get_float('reduce.bias', 'combine_sigma', 3.0)

        master_bias, uncertainty = combine_images(images_array, method=method, sigma=sigma)

        self.master_bias = master_bias
        self.bias_uncertainty = uncertainty

        logger.info(f"Master bias created: shape {master_bias.shape}, mean={np.mean(master_bias):.1f}")

        return master_bias, uncertainty

    def apply_bias_correction(self, science_image: np.ndarray,
                            master_bias: Optional[np.ndarray] = None) -> np.ndarray:
        """
        Apply bias correction to a science image.

        Args:
            science_image: Raw science image
            master_bias: Master bias frame (uses self.master_bias if None)

        Returns:
            Bias-corrected image
        """
        bias = master_bias if master_bias is not None else self.master_bias

        if bias is None:
            logger.warning("No master bias available, returning uncorrected image")
            return science_image

        if science_image.shape != bias.shape:
            raise ValueError(f"Image shape {science_image.shape} does not match "
                           f"bias shape {bias.shape}")

        corrected = science_image.astype(float) - bias.astype(float)
        logger.debug(f"Bias corrected image: mean={np.mean(corrected):.1f}")

        return corrected

    def save_master_bias(self, output_path: str):
        """Save master bias to FITS file."""
        if self.master_bias is None:
            raise RuntimeError("No master bias to save")

        header = {
            'EXTNAME': 'MASTER BIAS',
            'ORIGIN': 'SpecProc',
        }

        write_fits_image(output_path, self.master_bias, header=header)
        logger.info(f"Saved master bias to {output_path}")

        # Save diagnostic plot if enabled
        save_plots = self.config.get_bool('reduce', 'save_plots', True)
        if save_plots:
            out_dir = Path(output_path).parent
            fig_format = self.config.get('reduce', 'fig_format', 'png')
            plot_file = out_dir / f'master_bias.{fig_format}'
            plot_2d_image_to_file(self.master_bias, str(plot_file), "Master Bias Frame")

    def load_master_bias(self, filepath: str):
        """Load pre-computed master bias from file."""
        self.master_bias, header = read_fits_image(filepath)
        logger.info(f"Loaded master bias from {filepath}")


def process_bias_stage(config: ConfigManager, bias_filenames: List[str],
                      midpath: str = './midpath') -> Tuple[np.ndarray, str]:
    """
    Execute bias correction stage.

    Args:
        config: Configuration manager
        bias_filenames: List of bias frame paths
        midpath: Intermediate output directory

    Returns:
        Tuple of (master_bias, output_filepath)
    """
    corrector = BiasCorrector(config)

    # Combine bias frames
    master_bias, uncertainty = corrector.combine_bias_frames(bias_filenames)

    # Save master bias
    out_path = config.get('reduce', 'out_path', './output')
    output_path = Path(out_path) / 'step1_bias' / 'master_bias.fits'
    output_path.parent.mkdir(parents=True, exist_ok=True)
    corrector.save_master_bias(str(output_path))

    return master_bias, str(output_path)
