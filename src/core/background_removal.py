"""
Background subtraction stage of spectral reduction pipeline.

Handles inter-order background (stray light) estimation and removal.
"""

import numpy as np
import logging
from pathlib import Path
from typing import Tuple, Optional
from scipy.ndimage import median_filter
from src.utils.image_processing import estimate_background_2d
from src.config.config_manager import ConfigManager
from src.plotting.spectra_plotter import plot_background_residuals

logger = logging.getLogger(__name__)


class BackgroundRemover:
    """Handles background/stray light estimation and correction."""

    def __init__(self, config: ConfigManager):
        """
        Initialize background remover.

        Args:
            config: Configuration manager
        """
        self.config = config
        self.background_2d = None

    def estimate_background(self, image: np.ndarray, order: int = 2) -> np.ndarray:
        """
        Estimate 2D background using polynomial fitting.

        Args:
            image: 2D image array
            order: Polynomial order

        Returns:
            Background model
        """
        logger.info(f"Estimating 2D background (polyorder={order})...")

        background = estimate_background_2d(image, order=order)
        self.background_2d = background

        logger.info(f"Background estimated: min={np.min(background):.1f}, "
                   f"max={np.max(background):.1f}")

        return background

    def apply_background_correction(self, image: np.ndarray,
                                   background: Optional[np.ndarray] = None) -> np.ndarray:
        """
        Apply background correction to image.

        Args:
            image: Original science image
            background: Background model (uses self.background_2d if None)

        Returns:
            Background-corrected image
        """
        bkg = background if background is not None else self.background_2d

        if bkg is None:
            logger.warning("No background model available, returning uncorrected image")
            return image

        if image.shape != bkg.shape:
            raise ValueError(f"Image shape {image.shape} does not match "
                           f"background shape {bkg.shape}")

        corrected = image.astype(float) - bkg.astype(float)
        logger.debug(f"Background corrected: min={np.min(corrected):.1f}, "
                    f"max={np.max(corrected):.1f}")

        return corrected

    def extract_interorder_background(self, image: np.ndarray,
                                     aperture_positions: np.ndarray) -> np.ndarray:
        """
        Extract inter-order background by masking order regions.

        Args:
            image: 2D science image
            aperture_positions: Aperture center positions

        Returns:
            Inter-order background map
        """
        logger.info("Extracting inter-order background...")

        # Create aperture mask
        rows, cols = image.shape
        aperture_mask = np.zeros((rows, cols), dtype=bool)

        # Mark pixels within orders
        for x in range(cols):
            for aperture in aperture_positions:
                center = aperture['center'] + aperture.get('center_slope', 0) * x
                width = aperture.get('width', 20)
                y_min = max(0, int(center - width/2))
                y_max = min(rows, int(center + width/2))
                aperture_mask[y_min:y_max, x] = True

        # Image outside orders is inter-order light
        interorder_pixels = image[~aperture_mask]

        if len(interorder_pixels) > 0:
            median_bkg = np.median(interorder_pixels)
            logger.info(f"Median inter-order background: {median_bkg:.1f}")

            return median_bkg
        else:
            return 0.0

    def save_background(self, output_path: str):
        """Save background model to FITS file."""
        if self.background_2d is None:
            raise RuntimeError("No background model to save")

        from astropy.io import fits
        from pathlib import Path

        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        header = {
            'EXTNAME': '2D BACKGROUND',
            'ORIGIN': 'SpecProc',
        }

        hdu = fits.PrimaryHDU(data=self.background_2d.astype(np.float32), header=header)
        hdul = fits.HDUList([hdu])
        hdul.writeto(str(output_path), overwrite=True)

        logger.info(f"Saved background model to {output_path}")


def process_background_stage(config: ConfigManager, science_image: np.ndarray,
                            midpath: str = None) -> np.ndarray:
    """
    Execute background subtraction stage.

    Args:
        config: Configuration manager
        science_image: 2D science image
        midpath: Intermediate output directory

    Returns:
        Background model
    """
    # Use default midpath if not provided
    if midpath is None:
        base_output_path = config.get_output_path()
        midpath = str(Path(base_output_path) / 'midpath')

    remover = BackgroundRemover(config)

    # Get polynomial order from config
    poly_order = config.get_int('reduce.background', 'poly_order', 2)

    # Estimate background
    background = remover.estimate_background(science_image, order=poly_order)

    # Save background
    base_output_path = config.get_output_path()
    background_file = Path(base_output_path) / 'step3_background' / 'background_model.fits'
    background_file.parent.mkdir(parents=True, exist_ok=True)
    remover.save_background(str(background_file))

    # Save diagnostic plot if enabled
    save_plots = config.get_bool('reduce', 'save_plots', True)
    if save_plots:
        out_dir = background_file.parent
        fig_format = config.get('reduce', 'fig_format', 'png')
        plot_file = out_dir / f'background_residuals.{fig_format}'
        plot_background_residuals(science_image, background, str(plot_file))

    return background
