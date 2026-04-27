"""
Basic reduction stage — Step 1 of the spectral reduction pipeline.

Combines three sub-steps that are applied before order tracing:
  1. Overscan subtraction  (OverscanCorrectionStage / process_overscan_stage)
  2. Bias subtraction       (BiasCorrector / process_bias_stage)
  3. Cosmic-ray removal    (_detect_cosmics / _fix_cosmics / process_cosmic_stage)

The unified entry point is process_basic_reduction_stage(), which runs
whichever of the three sub-steps are enabled in the configuration.
"""


import numpy as np
import logging
import matplotlib.pyplot as plt
from pathlib import Path
from typing import List, Tuple, Optional

from src.utils.overscan import OverscanCorrector
from src.utils.fits_io import read_fits_image, write_fits_image, combine_fits_images
from src.utils.image_processing import combine_images
from src.config.config_manager import ConfigManager
from src.plotting.spectra_plotter import (
    plot_2d_image_to_file,
    plot_overscan_corrected_image,
)

# Use astroscrappy for L.A.Cosmic
try:
    import astroscrappy
except ImportError:
    astroscrappy = None

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Sub-step 1: Overscan correction
# ---------------------------------------------------------------------------

class OverscanCorrectionStage:
    """Handles overscan correction for raw images."""

    def __init__(self, 
                 overscan_start_column: int = -1,
                 overscan_left_start: int = -1, overscan_left_end: int = -1,
                 overscan_right_start: int = -1, overscan_right_end: int = -1,
                 overscan_top_start: int = -1, overscan_top_end: int = -1,
                 overscan_bottom_start: int = -1, overscan_bottom_end: int = -1,
                 overscan_method: str = 'mean_only',
                 overscan_smooth_window: int = -1,
                 overscan_poly_order: int = 3,
                 overscan_poly_type: str = 'legendre',
                 trim_x_start: int = -1, trim_x_end: int = -1,
                 trim_y_start: int = -1, trim_y_end: int = -1,
                 save_plots: bool = True,
                 fig_format: str = 'png'):
        self.overscan_start_column = overscan_start_column
        self.overscan_left_start = overscan_left_start
        self.overscan_left_end = overscan_left_end
        self.overscan_right_start = overscan_right_start
        self.overscan_right_end = overscan_right_end
        self.overscan_top_start = overscan_top_start
        self.overscan_top_end = overscan_top_end
        self.overscan_bottom_start = overscan_bottom_start
        self.overscan_bottom_end = overscan_bottom_end
        self.overscan_method = overscan_method
        self.overscan_smooth_window = overscan_smooth_window
        self.overscan_poly_order = overscan_poly_order
        self.overscan_poly_type = overscan_poly_type
        self.trim_x_start = trim_x_start
        self.trim_x_end = trim_x_end
        self.trim_y_start = trim_y_start
        self.trim_y_end = trim_y_end
        self.save_plots = save_plots
        self.fig_format = fig_format

        self.corrector = self._create_corrector()

    def _create_corrector(self) -> OverscanCorrector:
        """Create OverscanCorrector with configuration."""
        overscan_config = {}

        if self.overscan_start_column > 0:
            overscan_config['right'] = (self.overscan_start_column - 1, None)
            overscan_config['left'] = None
            overscan_config['top'] = None
            overscan_config['bottom'] = None
            logger.info(f"Simplified overscan configuration: start_column={self.overscan_start_column} (1-based)")
        else:
            regions = {
                'left': (self.overscan_left_start, self.overscan_left_end),
                'right': (self.overscan_right_start, self.overscan_right_end),
                'top': (self.overscan_top_start, self.overscan_top_end),
                'bottom': (self.overscan_bottom_start, self.overscan_bottom_end)
            }
            for region, (start, end) in regions.items():
                if start >= 0 and end > start:
                    overscan_config[region] = (start, end)
                else:
                    overscan_config[region] = None

        logger.info(f"Overscan configuration: {overscan_config}")
        return OverscanCorrector(overscan_config)

    def correct_image(self, image: np.ndarray, trim: bool = True) -> np.ndarray:
        """Correct single image for overscan."""
        has_overscan = self.overscan_start_column > 0
        logger.info(f"Overscan start column: {self.overscan_start_column}, has_overscan: {has_overscan}")

        if not has_overscan:
            logger.info("No overscan region configured, skipping overscan correction")
            return image

        logger.info(
            f"Overscan correction method: {self.overscan_method}, smooth_window: {self.overscan_smooth_window}, "
            f"poly_type: {self.overscan_poly_type}, poly_order: {self.overscan_poly_order}"
        )
        overscan_bias = self.corrector.estimate_overscan_bias(
            image, method=self.overscan_method,
            smooth_window=self.overscan_smooth_window if self.overscan_smooth_window > 0 else None,
            poly_order=self.overscan_poly_order,
            poly_type=self.overscan_poly_type,
        )

        corrected = self.corrector.apply_overscan_correction(image, overscan_bias)

        if trim and has_overscan:
            corrected = self.corrector.trim_image(corrected)

        corrected = self._apply_user_crop(corrected)
        return corrected

    def _apply_user_crop(self, image: np.ndarray) -> np.ndarray:
        """Apply optional crop bounds from config using DS9 1-based indexing."""
        rows, cols = image.shape

        x1 = self.trim_x_start
        x2 = self.trim_x_end
        y1 = self.trim_y_start
        y2 = self.trim_y_end

        if min(x1, x2, y1, y2) < 1:
            return image

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
        """Correct single FITS file for overscan."""
        logger.info(f"Reading raw image: {input_path}")
        raw_image, header = read_fits_image(input_path)

        corrected = self.correct_image(raw_image, trim=True)

        header['OVRSCAN'] = (True, 'Overscan correction completed')

        if output_path is None:
            output_path = input_path

        write_fits_image(output_path, corrected, header=header, overwrite=True, dtype='float32')
        logger.info(f"Saved overscan-corrected image: {output_path}")

        if self.save_plots and output_path != input_path:
            out_dir = Path(output_path).parent
            stem = Path(input_path).stem

            raw_plot_file = out_dir / f'{stem}_raw.{self.fig_format}'
            plot_2d_image_to_file(raw_image, str(raw_plot_file), "Raw Input Image")

            plot_file = out_dir / f'{stem}_overscan_corrected.{self.fig_format}'
            if self.overscan_start_column > 0:
                original_width = raw_image.shape[1]
                plot_overscan_corrected_image(
                    corrected, str(plot_file), "Overscan Corrected Image",
                    overscan_start_col=self.overscan_start_column,
                    original_width=original_width,
                )
            else:
                plot_2d_image_to_file(corrected, str(plot_file), "Overscan Corrected Image")

            if hasattr(self.corrector, 'overscan_profile') and self.corrector.overscan_profile:
                from src.plotting.spectra_plotter import plot_overscan_profile
                profile_file = out_dir / f'{stem}_overscan_profile.{self.fig_format}'
                plot_overscan_profile(
                    self.corrector.overscan_profile,
                    str(profile_file),
                    title=f"Overscan Profile - {stem}",
                )

        return output_path


def process_overscan_stage(image_filenames: List[str],
                           output_dir_base: str,
                           overscan_start_column: int = -1,
                           overscan_left_start: int = -1, overscan_left_end: int = -1,
                           overscan_right_start: int = -1, overscan_right_end: int = -1,
                           overscan_top_start: int = -1, overscan_top_end: int = -1,
                           overscan_bottom_start: int = -1, overscan_bottom_end: int = -1,
                           overscan_method: str = 'mean_only',
                           overscan_smooth_window: int = -1,
                           overscan_poly_order: int = 3,
                           overscan_poly_type: str = 'legendre',
                           trim_x_start: int = -1, trim_x_end: int = -1,
                           trim_y_start: int = -1, trim_y_end: int = -1,
                           save_plots: bool = True,
                           fig_format: str = 'png') -> List[str]:
    """Execute overscan correction for multiple images."""
    logger.info("=" * 60)
    logger.info("STAGE 0: OVERSCAN CORRECTION")
    logger.info("=" * 60)

    stage = OverscanCorrectionStage(
        overscan_start_column=overscan_start_column,
        overscan_left_start=overscan_left_start, overscan_left_end=overscan_left_end,
        overscan_right_start=overscan_right_start, overscan_right_end=overscan_right_end,
        overscan_top_start=overscan_top_start, overscan_top_end=overscan_top_end,
        overscan_bottom_start=overscan_bottom_start, overscan_bottom_end=overscan_bottom_end,
        overscan_method=overscan_method,
        overscan_smooth_window=overscan_smooth_window,
        overscan_poly_order=overscan_poly_order,
        overscan_poly_type=overscan_poly_type,
        trim_x_start=trim_x_start, trim_x_end=trim_x_end,
        trim_y_start=trim_y_start, trim_y_end=trim_y_end,
        save_plots=save_plots,
        fig_format=fig_format
    )

    overscan_path = Path(output_dir_base) / 'step1_basic' / 'overscan_corrected'
    overscan_path.mkdir(parents=True, exist_ok=True)

    corrected_files = []
    for filename in image_filenames:
        try:
            input_path = Path(filename)
            output_file = overscan_path / input_path.name
            corrected_file = stage.correct_file(str(filename), str(output_file))
            corrected_files.append(corrected_file)
        except Exception as e:
            logger.error(f"Error correcting {filename}: {e}")
            corrected_files.append(filename)

    logger.info(f"✓ Overscan correction complete: {len(corrected_files)} images processed")
    return corrected_files


# ---------------------------------------------------------------------------
# Sub-step 2: Bias correction
# ---------------------------------------------------------------------------

class BiasCorrector:
    """Handles bias frame combination and correction."""

    def __init__(self, 
                 combine_method: str = 'median',
                 combine_sigma: float = 3.0,
                 save_plots: bool = True,
                 fig_format: str = 'png'):
        self.combine_method = combine_method
        self.combine_sigma = combine_sigma
        self.save_plots = save_plots
        self.fig_format = fig_format
        self.master_bias = None
        self.bias_uncertainty = None

    def combine_bias_frames(self, bias_filenames: List[str]) -> Tuple[np.ndarray, np.ndarray]:
        """Combine multiple bias frames to create master bias."""
        logger.info(f"Combining {len(bias_filenames)} bias frames...")

        if not bias_filenames:
            logger.warning("No bias frames provided")
            return None, None

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

        images_array = np.array(bias_images)
        logger.info(f"Successfully read {len(images_array)} bias frames")

        master_bias, uncertainty = combine_images(images_array, method=self.combine_method, sigma=self.combine_sigma)

        self.master_bias = master_bias
        self.bias_uncertainty = uncertainty

        logger.info(f"Master bias created: shape {master_bias.shape}, mean={np.mean(master_bias):.1f}")
        return master_bias, uncertainty

    def apply_bias_correction(self, science_image: np.ndarray,
                              master_bias: Optional[np.ndarray] = None) -> np.ndarray:
        """Apply bias correction to a science image."""
        bias = master_bias if master_bias is not None else self.master_bias

        if bias is None:
            logger.warning("No master bias available, returning uncorrected image")
            return science_image

        if science_image.shape != bias.shape:
            raise ValueError(
                f"Image shape {science_image.shape} does not match bias shape {bias.shape}"
            )

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
            'BIASCR': (True, 'Master bias created'),
        }

        write_fits_image(output_path, self.master_bias, header=header, dtype='float32')
        logger.info(f"Saved master bias to {output_path}")

        if self.save_plots:
            out_dir = Path(output_path).parent
            plot_file = out_dir / f'master_bias.{self.fig_format}'
            plot_2d_image_to_file(self.master_bias, str(plot_file), "Master Bias Frame")

    def load_master_bias(self, filepath: str):
        """Load pre-computed master bias from file."""
        self.master_bias, header = read_fits_image(filepath)
        logger.info(f"Loaded master bias from {filepath}")


def process_bias_stage(bias_filenames: List[str],
                       output_dir_base: str,
                       science_file: str = None,
                       combine_method: str = 'median',
                       combine_sigma: float = 3.0,
                       save_plots: bool = True,
                       fig_format: str = 'png') -> Tuple[np.ndarray, str]:
    """Execute bias correction stage."""
    corrector = BiasCorrector(
        combine_method=combine_method,
        combine_sigma=combine_sigma,
        save_plots=save_plots,
        fig_format=fig_format
    )

    master_bias, uncertainty = corrector.combine_bias_frames(bias_filenames)

    output_path = Path(output_dir_base) / 'step1_basic' / 'bias_subtracted'
    output_path.mkdir(parents=True, exist_ok=True)
    master_bias_path = output_path / 'master_bias.fits'
    corrector.save_master_bias(str(master_bias_path))

    if science_file is not None:
        science_image, science_header = read_fits_image(science_file)
        logger.info(f"Applying bias correction to science image: {science_file}")

        bias_corrected = corrector.apply_bias_correction(science_image, master_bias)

        output_science_path = output_path / Path(science_file).name
        science_header['BIASCOR'] = (True, 'Bias correction applied')
        write_fits_image(str(output_science_path), bias_corrected,
                         header=science_header, dtype='float32')
        logger.info(f"Saved bias-corrected science image to {output_science_path}")

        if save_plots:
            plot_file = output_path / f'{Path(science_file).stem}_bias_corrected.{fig_format}'
            plot_2d_image_to_file(bias_corrected, str(plot_file), "Bias Corrected Science Image")

    return master_bias, str(master_bias_path)


# ---------------------------------------------------------------------------
# Sub-step 3: Cosmic-ray removal
# ---------------------------------------------------------------------------

def process_cosmic_stage(image_filenames: List[str],
                         output_dir_base: str,
                         output_subdir: str = 'step1_basic/cosmic_corrected',
                         cosmic_sigclip: float = 5.0,
                         cosmic_objlim: float = 5.0,
                         cosmic_gain: float = 1.0,
                         cosmic_readnoise: float = 5.0,
                         save_plots: bool = True,
                         fig_format: str = 'png') -> List[str]:
    """Execute cosmic ray correction stage for a list of images."""
    if astroscrappy is None:
        logger.error("`astroscrappy` package not found. Skipping cosmic ray removal.")
        logger.error("Please install it via: pip install astroscrappy")
        return image_filenames

    logger.info("=" * 60)
    logger.info("STEP 1: COSMIC RAY REMOVAL")
    logger.info("=" * 60)

    cosmic_path = Path(output_dir_base) / output_subdir
    cosmic_path.mkdir(parents=True, exist_ok=True)

    corrected_files = []

    for filename in image_filenames:
        try:
            input_path = Path(filename)
            output_file = cosmic_path / input_path.name

            image, header = read_fits_image(str(filename))

            logger.info(f"Using L.A.Cosmic with sigclip={cosmic_sigclip}, objlim={cosmic_objlim}")

            # Prepare arguments for astroscrappy. If gain/readnoise are 0,
            # do not pass them, allowing astroscrappy to auto-detect from header.
            lacosmic_kwargs = {
                'sigclip': cosmic_sigclip,
                'objlim': cosmic_objlim,
                'cleantype': 'medmask'
            }
            if cosmic_gain > 0:
                lacosmic_kwargs['gain'] = cosmic_gain
            if cosmic_readnoise > 0:
                lacosmic_kwargs['readnoise'] = cosmic_readnoise

            mask, corrected = astroscrappy.detect_cosmics(image, **lacosmic_kwargs)

            header['COSMICR'] = (True, 'Cosmic ray removal completed (L.A.Cosmic)')
            write_fits_image(str(output_file), corrected, header=header,
                             overwrite=True, dtype='float32')

            if save_plots:
                plot_file = cosmic_path / f"{input_path.stem}_cosmic_corrected.{fig_format}"
                plot_2d_image_to_file(corrected, str(plot_file), "Cosmic Corrected Image")

                marked_file = cosmic_path / f"{input_path.stem}_cosmic_marked.{fig_format}"
                fig, ax = plt.subplots(figsize=(10, 8))
                vmin = np.percentile(image, 1)
                vmax = np.percentile(image, 99)
                ax.imshow(image, cmap='viridis', origin='lower', vmin=vmin, vmax=vmax)
                ys, xs = np.where(mask)
                if xs.size > 0:
                    ax.scatter(
                        xs, ys, s=10,
                        facecolors='none', edgecolors='red',
                        linewidths=0.6, alpha=0.9,
                        label=f'Cosmic pixels: {xs.size}',
                    )
                    ax.legend(loc='upper right')
                ax.set_title('Detected Cosmic Rays (marked in red)')
                ax.set_xlabel('Column')
                ax.set_ylabel('Row')
                fig.tight_layout()
                fig.savefig(str(marked_file), dpi=150)
                plt.close(fig)

            corrected_files.append(str(output_file))
            logger.info(
                f"Cosmic correction: {filename} -> {output_file} ({np.sum(mask)} pixels fixed)"
            )

        except Exception as e:
            logger.error(f"Error cosmic-correcting {filename}: {e}")
            corrected_files.append(str(filename))

    logger.info(f"✓ Cosmic correction complete: {len(corrected_files)} images processed")
    return corrected_files


# ---------------------------------------------------------------------------
# Unified Step-1 entry point
# ---------------------------------------------------------------------------

def process_basic_reduction_stage(science_filenames: List[str],
                                  output_dir_base: str,
                                  bias_filenames: Optional[List[str]] = None,
                                  do_overscan: bool = True,
                                  do_bias: bool = True,
                                  do_cosmic: bool = True,
                                  overscan_kwargs: dict = None,
                                  bias_kwargs: dict = None, 
                                  cosmic_kwargs: dict = None) -> Tuple[List[str], Optional[str]]:
    """
    Run all enabled sub-steps of Step 1 (basic pre-processing) in order.

    Parameters
    ----------
    science_filenames : Raw (or partially-corrected) science FITS file list.
    output_dir_base   : Base directory for outputs.
    bias_filenames    : Bias frames for master-bias creation (needed if do_bias=True).
    do_overscan       : Run overscan correction.
    do_bias           : Run bias subtraction.
    do_cosmic         : Run cosmic-ray removal.

    Returns
    -------
    (processed_science_files, master_bias_path)
        master_bias_path is None when do_bias is False.
    """
    logger.info("=" * 60)
    logger.info("STEP 1: BASIC PRE-PROCESSING")
    logger.info("=" * 60)

    overscan_kwargs = overscan_kwargs or {}
    bias_kwargs = bias_kwargs or {}
    cosmic_kwargs = cosmic_kwargs or {}

    current_files = list(science_filenames)
    master_bias_path: Optional[str] = None

    if do_overscan:
        current_files = process_overscan_stage(current_files, output_dir_base, **overscan_kwargs)

    if do_bias and bias_filenames:
        # Apply bias correction to each science file individually.
        processed = []
        first = True
        for sci_file in current_files:
            if first:
                _, master_bias_path = process_bias_stage(bias_filenames, output_dir_base, sci_file, **bias_kwargs)
                first = False
            else:
                # Reuse the master bias already created; just apply it.
                corrector = BiasCorrector(**bias_kwargs)
                corrector.load_master_bias(master_bias_path)
                sci_img, sci_hdr = read_fits_image(sci_file)
                corrected = corrector.apply_bias_correction(sci_img)
                out_dir = Path(output_dir_base) / 'step1_basic' / 'bias_subtracted'
                out_dir.mkdir(parents=True, exist_ok=True)
                out_file = str(out_dir / Path(sci_file).name)
                sci_hdr['BIASCOR'] = (True, 'Bias correction applied')
                write_fits_image(out_file, corrected, header=sci_hdr, dtype='float32')
                processed.append(out_file)
                continue
            # First file: path was written inside process_bias_stage; reconstruct it.
            out_dir = Path(output_dir_base) / 'step1_basic' / 'bias_subtracted'
            processed.append(str(out_dir / Path(sci_file).name))
        current_files = processed

    if do_cosmic:
        current_files = process_cosmic_stage(current_files, output_dir_base, **cosmic_kwargs)

    logger.info(f"✓ Step 1 basic reduction complete: {len(current_files)} science files")
    return current_files, master_bias_path
