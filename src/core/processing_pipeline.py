"""
Main processing pipeline orchestrator.

Coordinates all processing stages and manages workflow execution.
"""

import logging
from pathlib import Path
from typing import List, Callable
import numpy as np
from src.config.config_manager import ConfigManager
from src.core.data_structures import ProcessingState, WaveCalib, SpectraSet
from src.core.overscan_correction import process_overscan_stage
from src.core.cosmic_removal import process_cosmic_stage
from src.core.bias_correction import process_bias_stage
from src.core.flat_fielding import process_flat_stage
from src.core.wave_calibration import WavelengthCalibrator, process_wavelength_stage
from src.core.background_removal import process_background_stage
from src.core.extraction import process_extraction_stage
from src.core.de_blazing import process_de_blazing_stage
from src.utils.fits_io import read_fits_image, write_fits_image

logger = logging.getLogger(__name__)


class ProcessingPipeline:
    """Main spectral reduction pipeline."""

    def __init__(self, config: ConfigManager):
        """
        Initialize processing pipeline.

        Args:
            config: Configuration manager
        """
        self.config = config
        self.state = ProcessingState()
        self.progress_callback = None

    def set_progress_callback(self, callback: Callable[[float, str], None]):
        """
        Set callback for progress updates.

        Args:
            callback: Function taking (progress_fraction, stage_name)
        """
        self.progress_callback = callback

    def _report_progress(self, progress: float, stage: str):
        """Report progress to callback."""
        self.state.progress = progress
        self.state.current_stage = stage
        if self.progress_callback:
            self.progress_callback(progress, stage)
            logger.info(f"{stage}: {progress*100:.1f}%")

    def stage_overscan_correction(self, raw_filenames: List[str]) -> List[str]:
        """
        Execute overscan correction stage (STAGE 0 - must be first!).

        Args:
            raw_filenames: List of raw image paths

        Returns:
            List of overscan-corrected image paths
        """
        logger.info("=" * 50)
        logger.info("STAGE 0: OVERSCAN CORRECTION")
        logger.info("=" * 50)

        self._report_progress(0.01, "Overscan Correction")

        try:
            corrected_files = process_overscan_stage(self.config, raw_filenames, None)

            self._report_progress(0.05, "Overscan Correction Complete")
            logger.info("✓ Overscan correction stage complete")

            return corrected_files

        except Exception as e:
            logger.error(f"Error in overscan correction: {e}")
            # Return original files if overscan correction fails
            return raw_filenames

    def stage_cosmic_correction(self, filenames: List[str]) -> List[str]:
        """
        Execute cosmic ray correction stage.

        Args:
            filenames: List of image paths to correct

        Returns:
            List of cosmic-corrected image paths
        """
        if not self.config.get_bool('reduce', 'cosmic_enabled', True):
            logger.info("Cosmic correction disabled, skipping stage.")
            return filenames

        logger.info("=" * 50)
        logger.info("STAGE 0.5: COSMIC RAY CORRECTION")
        logger.info("=" * 50)

        self._report_progress(0.06, "Cosmic Ray Correction")

        try:
            corrected_files = process_cosmic_stage(self.config, filenames, None)

            self._report_progress(0.08, "Cosmic Ray Correction Complete")
            logger.info("✓ Cosmic ray correction stage complete")

            return corrected_files

        except Exception as e:
            logger.error(f"Error in cosmic ray correction: {e}")
            return filenames

    def stage_bias_correction(self, bias_filenames: List[str]) -> np.ndarray:
        """
        Execute bias correction stage.

        Note: Bias frames are short exposure (0s), so no cosmic ray correction needed.
        Use mean/median combination directly.

        Args:
            bias_filenames: List of bias frame paths

        Returns:
            Master bias image
        """
        logger.info("=" * 50)
        logger.info("STAGE 1: BIAS CORRECTION")
        logger.info("=" * 50)

        self._report_progress(0.05, "Bias Correction")

        try:
            master_bias, bias_file = process_bias_stage(self.config, bias_filenames, None)

            self.state.bias_frame = master_bias
            self.state.bias_done = True

            self._report_progress(0.15, "Bias Correction Complete")
            logger.info("✓ Bias correction stage complete")

            return master_bias

        except Exception as e:
            logger.error(f"Error in bias correction: {e}")
            raise

    def stage_flat_fielding(self, flat_filenames: List[str]) -> tuple:
        """
        Execute flat fielding stage.

        Note: Flat frames should have cosmic ray correction applied before this stage.

        Args:
            flat_filenames: List of flat frame paths

        Returns:
            Tuple of (FlatField, ApertureSet)
        """
        logger.info("=" * 50)
        logger.info("STAGE 3: FLAT FIELDING & ORDER TRACING")
        logger.info("=" * 50)

        self._report_progress(0.15, "Flat Fielding")

        try:
            flat_field, apertures = process_flat_stage(self.config, flat_filenames, None)

            self.state.flat_field = flat_field
            self.state.apertures = apertures
            self.state.flat_done = True
            self.state.trace_done = True

            self._report_progress(0.30, "Flat Fielding Complete")
            logger.info("✓ Flat fielding stage complete")

            return flat_field, apertures

        except Exception as e:
            logger.error(f"Error in flat fielding: {e}")
            raise

    def stage_wavelength_calibration(self, calib_filename: str) -> WaveCalib:
        """
        Execute wavelength calibration stage on calibration frame.

        Args:
            calib_filename: Path to calibration (ThAr) frame

        Returns:
            WaveCalib object
        """
        logger.info("=" * 50)
        logger.info("STAGE 4: WAVELENGTH CALIBRATION (Calibration Frame)")
        logger.info("=" * 50)

        self._report_progress(0.30, "Wavelength Calibration")

        try:
            wave_calib = process_wavelength_stage(self.config, calib_filename, None)

            self.state.wavelength_calib = wave_calib
            self.state.wavelength_done = True

            self._report_progress(0.30, "Wavelength Calibration Complete")
            logger.info("✓ Wavelength calibration stage complete")

            return wave_calib

        except Exception as e:
            logger.error(f"Error in wavelength calibration: {e}")
            raise

    def stage_wavelength_calibration_on_spectra(self, spectra_set: SpectraSet,
                                                calib_filename: str) -> SpectraSet:
        """
        Apply wavelength calibration to extracted science spectra.

        This is the correct order: extract 1D spectra first, then apply wavelength calibration.
        Following the gamse standard pipeline.

        Args:
            spectra_set: Extracted spectra without wavelength information
            calib_filename: Path to calibration (ThAr) frame for reference

        Returns:
            SpectraSet with wavelength-calibrated spectra
        """
        logger.info("=" * 50)
        logger.info("STAGE 7: WAVELENGTH CALIBRATION (Science Spectra)")
        logger.info("=" * 50)

        self._report_progress(0.75, "Applying Wavelength Calibration")

        try:
            # First, calibrate the calibration frame if not already done
            if self.state.wavelength_calib is None:
                wave_calib = process_wavelength_stage(self.config, calib_filename, None)
                self.state.wavelength_calib = wave_calib
            else:
                wave_calib = self.state.wavelength_calib

            # Apply wavelength calibration to each extracted spectrum
            calibrator = WavelengthCalibrator(self.config)
            calibrator.wave_calib = wave_calib

            calibrated_spectra = SpectraSet()

            for spectrum in spectra_set.spectra:
                # Create pixel array
                n_pixels = len(spectrum.flux)
                pixel_array = np.arange(n_pixels)

                # Apply wavelength calibration
                wavelengths = calibrator.apply_wavelength_calibration(
                    spectrum.flux, pixel_array, aperture_y=float(spectrum.aperture)
                )

                # Update spectrum with wavelengths
                calibrated_spectrum = spectrum.copy()
                calibrated_spectrum.wavelength = wavelengths

                calibrated_spectra.add_spectrum(calibrated_spectrum)

            self._report_progress(0.85, "Wavelength Calibration Complete")
            logger.info(f"✓ Wavelength calibration applied to {len(calibrated_spectra.spectra)} spectra")

            return calibrated_spectra

        except Exception as e:
            logger.error(f"Error in wavelength calibration of spectra: {e}")
            raise

    def stage_background_removal(self, science_image: np.ndarray) -> np.ndarray:
        """
        Execute background subtraction stage.

        Args:
            science_image: 2D science image

        Returns:
            Background model
        """
        logger.info("=" * 50)
        logger.info("STAGE 5: BACKGROUND SUBTRACTION")
        logger.info("=" * 50)

        self._report_progress(0.45, "Background Removal")

        try:
            background = process_background_stage(self.config, science_image, None)

            self.state.background_done = True

            self._report_progress(0.60, "Background Removal Complete")
            logger.info("✓ Background removal stage complete")

            return background

        except Exception as e:
            logger.error(f"Error in background removal: {e}")
            raise

    def stage_cosmic_correction_science(self, science_image: np.ndarray) -> np.ndarray:
        """
        Execute cosmic ray correction stage for science images only.

        This stage is placed after background subtraction and before spectrum extraction.

        Why only science images need cosmic ray correction:
        - Bias frames: 0s exposure, no cosmic rays accumulate
        - Flat frames: Multiple frames combined, cosmic rays removed during combination
        - ThAr frames: Lamp spectra, typically no cosmic ray correction needed
        - Science frames: Long exposure (minutes to hours), single frame, cosmic rays must be removed

        Args:
            science_image: Background-corrected 2D science image

        Returns:
            Cosmic-ray corrected science image
        """
        logger.info("=" * 50)
        logger.info("STAGE 6: COSMIC RAY CORRECTION (Science Frame)")
        logger.info("=" * 50)

        self._report_progress(0.60, "Cosmic Ray Correction")

        try:
            # Apply cosmic ray correction directly to the image
            # Using sigma-threshold detection and interpolation
            cosmic_sigma = self.config.get_float('reduce', 'cosmic_sigma', 5.0)
            cosmic_window = self.config.get_int('reduce', 'cosmic_window', 7)

            # Simple cosmic ray detection: identify pixels > mean + N * std
            mean_val = np.mean(science_image)
            std_val = np.std(science_image)
            threshold = mean_val + cosmic_sigma * std_val

            # Create mask for cosmic rays
            cosmic_mask = science_image > threshold

            logger.info(f"Detected {np.sum(cosmic_mask)} cosmic ray affected pixels "
                       f"(threshold: {threshold:.2f}, sigma: {cosmic_sigma})")

            # Interpolate over cosmic ray pixels
            from scipy.ndimage import median_filter
            corrected_image = science_image.copy()

            if np.sum(cosmic_mask) > 0:
                # Use median filter to replace cosmic ray pixels
                # Get local median for each cosmic ray pixel
                filtered = median_filter(science_image, size=cosmic_window)
                corrected_image[cosmic_mask] = filtered[cosmic_mask]

            self._report_progress(0.65, "Cosmic Ray Correction Complete")
            logger.info(f"✓ Cosmic ray correction complete: "
                       f"{np.sum(cosmic_mask)} pixels corrected")

            return corrected_image

        except Exception as e:
            logger.error(f"Error in cosmic ray correction: {e}")
            raise

    def stage_extraction(self, science_image: np.ndarray) -> SpectraSet:
        """
        Execute spectrum extraction stage.

        Args:
            science_image: Background-corrected 2D science image

        Returns:
            SpectraSet with extracted 1D spectra (without wavelength calibration)
        """
        logger.info("=" * 50)
        logger.info("STAGE 7: SPECTRUM EXTRACTION")
        logger.info("=" * 50)

        self._report_progress(0.80, "Spectrum Extraction")

        if self.state.apertures is None:
            raise RuntimeError("Apertures not available. Run flat fielding stage first.")

        try:
            # Extract without wavelength calibration (will be applied in next stage)
            extracted_spectra = process_extraction_stage(
                self.config, science_image, self.state.apertures,
                None,  # No wavelength calibration yet
                self.state.flat_field, None
            )

            self.state.extracted_spectra = extracted_spectra
            self.state.extraction_done = True

            self._report_progress(0.75, "Spectrum Extraction Complete")
            logger.info("✓ Spectrum extraction stage complete")

            return extracted_spectra

        except Exception as e:
            logger.error(f"Error in spectrum extraction: {e}")
            raise

    def stage_de_blazing(self, spectra_set: SpectraSet) -> SpectraSet:
        """
        Execute de-blazing correction stage.

        Args:
            spectra_set: Extracted spectra

        Returns:
            De-blazed spectra
        """
        logger.info("=" * 50)
        logger.info("STAGE 9: DE-BLAZING CORRECTION")
        logger.info("=" * 50)

        self._report_progress(0.95, "De-blazing Correction")

        try:
            de_blazed_spectra = process_de_blazing_stage(
                self.config, spectra_set, self.state.flat_field
            )

            self.state.de_blazed_spectra = de_blazed_spectra
            self.state.de_blazing_done = True

            self._report_progress(1.0, "Pipeline Complete")
            logger.info("✓ De-blazing stage complete")

            return de_blazed_spectra

        except Exception as e:
            logger.error(f"Error in de-blazing: {e}")
            raise

    def run_full_pipeline(self, raw_image_path: str,
                         bias_filenames: List[str],
                         flat_filenames: List[str],
                         calib_filename: str) -> SpectraSet:
        """
        Execute complete spectral reduction pipeline.

        Args:
            raw_image_path: Path to science FITS file
            bias_filenames: List of bias frame paths
            flat_filenames: List of flat frame paths
            calib_filename: Path to wavelength calibration frame

        Returns:
            SpectraSet with final extracted spectra
        """
        logger.info("=" * 60)
        logger.info("STARTING FULL SPECTRAL REDUCTION PIPELINE")
        logger.info("=" * 60)

        # Reset state
        self.state.reset()

        # Create output directories
        self.config.create_directories()

        try:
            # Stage 0: Overscan correction (MUST be first! - apply to ALL images)
            all_raw_files = [raw_image_path] + bias_filenames + flat_filenames + [calib_filename]
            corrected_files = self.stage_overscan_correction(all_raw_files)

            # Extract corrected filenames
            raw_image_path = corrected_files[0]
            bias_filenames = corrected_files[1:1+len(bias_filenames)]
            flat_filenames = corrected_files[1+len(bias_filenames):1+len(bias_filenames)+len(flat_filenames)]
            calib_filename = corrected_files[-1]

            # Stage 1: Bias correction (bias is short exposure, use mean/median combine)
            master_bias = self.stage_bias_correction(bias_filenames)

            # Stage 2: Flat fielding
            flat_field, apertures = self.stage_flat_fielding(flat_filenames)

            # Read science image (already overscan-corrected)
            logger.info(f"Reading science image: {raw_image_path}")
            science_image, header = read_fits_image(raw_image_path)

            # Apply bias correction
            science_corrected = science_image.astype(float) - master_bias

            # Stage 3: Background subtraction
            background = self.stage_background_removal(science_corrected)
            science_final = science_corrected - background

            # Stage 4: Cosmic ray correction for science frame only
            # Only science images need cosmic ray correction:
            # - Bias: 0s exposure, no cosmic rays
            # - Flat: multiple frames combined, cosmic rays removed during combination
            # - ThAr: lamp spectra, typically no cosmic ray correction needed
            # - Science: long exposure, single frame, cosmic rays must be removed
            if self.config.get_bool('reduce', 'cosmic_enabled', True):
                logger.info("Applying cosmic ray correction to science frame...")
                science_final = self.stage_cosmic_correction_science(science_final)

            # Stage 5: Spectrum extraction (extract 1D spectra first)
            extracted_spectra = self.stage_extraction(science_final)

            # Stage 8: Wavelength calibration (applied to extracted 1D spectra)
            calibrated_spectra = self.stage_wavelength_calibration_on_spectra(extracted_spectra, calib_filename)

            # Stage 9: De-blazing correction
            final_spectra = self.stage_de_blazing(calibrated_spectra)

            # Save final spectra to output/step8_final_spectra/
            out_path = self.config.get('reduce', 'out_path', './output')
            spectra_dir = Path(out_path) / 'step8_final_spectra'
            spectra_dir.mkdir(parents=True, exist_ok=True)

            # Generate output filename from science image
            science_filename = Path(raw_image_path).stem
            final_spectra_file = spectra_dir / f'{science_filename}_final.fits'

            from src.core.de_blazing import save_deblazed_spectra
            save_deblazed_spectra(str(final_spectra_file), final_spectra)

            # Save diagnostic plots if enabled
            save_plots = self.config.get_bool('reduce', 'save_plots', True)
            if save_plots:
                from src.plotting.spectra_plotter import plot_spectrum_to_file
                fig_format = self.config.get('reduce', 'fig_format', 'png')
                for spectrum in final_spectra.spectra.values():
                    order = spectrum.aperture
                    plot_file = spectra_dir / f'final_order_{order:02d}.{fig_format}'
                    plot_spectrum_to_file(
                        spectrum.wavelength,
                        spectrum.flux,
                        str(plot_file),
                        spectrum.error if spectrum.error is not None else None,
                        f"Final Spectrum - Order {order}"
                    )

            self._report_progress(1.0, "Pipeline Complete")

            logger.info("=" * 60)
            logger.info("PIPELINE EXECUTION COMPLETED SUCCESSFULLY")
            logger.info("=" * 60)

            return final_spectra

        except Exception as e:
            logger.error(f"Pipeline execution failed: {e}")
            logger.error("=" * 60)
            raise

    def get_state(self) -> ProcessingState:
        """Get current processing state."""
        return self.state

    def reset(self):
        """Reset pipeline state."""
        self.state.reset()
