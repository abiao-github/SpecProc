"""
Main processing pipeline orchestrator.

Coordinates all processing stages and manages workflow execution.
"""

import logging
import pickle
from pathlib import Path
from typing import List, Callable, Optional
import numpy as np
from src.config.config_manager import ConfigManager
from src.core.data_structures import ProcessingState, WaveCalib, SpectraSet
from src.core.overscan_correction import process_overscan_stage
from src.core.cosmic_removal import process_cosmic_stage
from src.core.bias_correction import process_bias_stage
from src.core.flat_fielding import process_flat_stage, FlatFieldProcessor
from src.core.wave_calibration import WavelengthCalibrator, process_wavelength_stage
from src.core.background_removal import process_background_stage
from src.core.extraction import process_extraction_stage
from src.core.de_blazing import process_de_blazing_stage
from src.core.order_stitching import process_order_stitching_stage
from src.utils.fits_io import read_fits_image, write_fits_image
from src.plotting.spectra_plotter import plot_2d_image_to_file

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
        logger.info("STEP 1: OVERSCAN CORRECTION")
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
        logger.info("STEP 1: BASIC PRE-PROCESSING (BIAS)")
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

    def _apply_master_bias_to_files(self, input_files: List[str], master_bias: np.ndarray,
                                    output_subdir: str) -> List[str]:
        """Apply master bias to a list of files and save corrected copies."""
        if not input_files:
            return []

        output_dir = Path(self.config.get_output_path()) / 'step1_basic' / 'bias_subtracted'
        output_dir.mkdir(parents=True, exist_ok=True)
        save_plots = self.config.get_bool('reduce', 'save_plots', True)
        fig_format = self.config.get('reduce', 'fig_format', 'png')

        corrected_files = []
        for filename in input_files:
            image, header = read_fits_image(filename)

            if image.shape != master_bias.shape:
                raise ValueError(
                    f"Cannot apply bias to {filename}: image shape {image.shape} "
                    f"!= master bias shape {master_bias.shape}"
                )

            corrected = image.astype(np.float32) - master_bias.astype(np.float32)
            header['BIASCOR'] = (True, 'Bias correction applied')

            output_file = output_dir / Path(filename).name
            write_fits_image(str(output_file), corrected, header=header, dtype='float32')

            if save_plots:
                plot_file = output_dir / f"{Path(filename).stem}_bias_corrected.{fig_format}"
                plot_2d_image_to_file(
                    corrected,
                    str(plot_file),
                    f"Bias Corrected {output_subdir.capitalize()} Image"
                )

            corrected_files.append(str(output_file))

        logger.info(
            f"Applied master bias to {len(corrected_files)} frame(s) for {output_subdir}; "
            f"saved in {output_dir}"
        )
        return corrected_files

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
        logger.info("STEP 2: ORDER TRACING")
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
        logger.info("STEP 6: WAVELENGTH CALIBRATION (Science Spectra)")
        logger.info("=" * 50)

        self._report_progress(0.82, "Applying Wavelength Calibration")

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

            for spectrum in spectra_set.spectra.values():
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

            self._report_progress(0.88, "Wavelength Calibration Complete")
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
            background = process_background_stage(
                self.config,
                science_image,
                None,
                output_subdir='step2_scatterlight',
                output_tag='science_',
                apertures=self.state.apertures,
                mask_margin_scale=1.2,
            )

            self.state.background_done = True

            self._report_progress(0.60, "Background Removal Complete")
            logger.info("✓ Background removal stage complete")

            return background

        except Exception as e:
            logger.error(f"Error in background removal: {e}")
            raise

    def stage_cosmic_correction_science(self, science_image: np.ndarray,
                                       science_name: str = 'science') -> np.ndarray:
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
        logger.info("STEP 1: COSMIC RAY CORRECTION (Science Frame)")
        logger.info("=" * 50)

        self._report_progress(0.45, "Cosmic Ray Correction")

        try:
            # Reuse the shared Step 1 cosmic detector so selected-step execution
            # and full-pipeline execution behave consistently.
            from src.core.cosmic_removal import _detect_cosmics, _fix_cosmics

            cosmic_sigma = self.config.get_float('reduce', 'cosmic_sigma', 5.0)
            cosmic_window = self.config.get_int('reduce', 'cosmic_window', 5)
            cosmic_fine_sigma = self.config.get_float('reduce', 'cosmic_fine_sigma', 2.0)
            cosmic_line_sigma = self.config.get_float('reduce', 'cosmic_line_sigma', 1.5)
            cosmic_grow_sigma = self.config.get_float('reduce', 'cosmic_grow_sigma', 2.5)
            cosmic_maxsize = self.config.get_int('reduce', 'cosmic_maxsize', 8)

            cosmic_mask = _detect_cosmics(
                science_image,
                sigma=cosmic_sigma,
                window=cosmic_window,
                fine_sigma=cosmic_fine_sigma,
                line_sigma=cosmic_line_sigma,
                grow_sigma=cosmic_grow_sigma,
                maxsize=cosmic_maxsize,
            )

            logger.info(
                f"Detected {np.sum(cosmic_mask)} cosmic ray affected pixels "
                f"(sigma={cosmic_sigma}, window={cosmic_window}, fine_sigma={cosmic_fine_sigma}, "
                f"line_sigma={cosmic_line_sigma})"
            )

            corrected_image = _fix_cosmics(science_image, cosmic_mask, window=cosmic_window)

            out_dir = Path(self.config.get_output_path()) / 'step1_basic' / 'cosmic_corrected'
            out_dir.mkdir(parents=True, exist_ok=True)
            write_fits_image(
                str(out_dir / f'{science_name}_science_cosmic_corrected.fits'),
                corrected_image.astype(np.float32),
                dtype='float32'
            )

            if self.config.get_bool('reduce', 'save_plots', True):
                fig_format = self.config.get('reduce', 'fig_format', 'png')
                plot_2d_image_to_file(
                    corrected_image,
                    str(out_dir / f'{science_name}_science_cosmic_corrected.{fig_format}'),
                    'Science After Cosmic Ray Correction'
                )

            self._report_progress(0.48, "Cosmic Ray Correction Complete")
            logger.info(f"✓ Cosmic ray correction complete: "
                       f"{np.sum(cosmic_mask)} pixels corrected")

            return corrected_image

        except Exception as e:
            logger.error(f"Error in cosmic ray correction: {e}")
            raise

    def stage_extraction(self, science_image: np.ndarray) -> SpectraSet:
        """
        Execute spectrum extraction stage and output both sum/optimal products.

        Args:
            science_image: Background-corrected 2D science image

        Returns:
            SpectraSet with extracted 1D spectra (without wavelength calibration)
        """
        logger.info("=" * 50)
        logger.info("STEP 5: SPECTRUM EXTRACTION")
        logger.info("=" * 50)

        self._report_progress(0.74, "Spectrum Extraction")

        if self.state.apertures is None:
            raise RuntimeError("Apertures not available. Run flat fielding stage first.")

        try:
            # Step 5 requirement: output both aperture-sum and optimal extraction.
            extracted_sum = process_extraction_stage(
                self.config, science_image, self.state.apertures,
                None,  # No wavelength calibration yet
                self.state.flat_field, None,
                method_override='sum',
                output_filename='extracted_spectra_sum.fits',
                plot_prefix='extracted_sum',
            )

            extracted_optimal = process_extraction_stage(
                self.config, science_image, self.state.apertures,
                None,
                self.state.flat_field, None,
                method_override='optimal',
                output_filename='extracted_spectra_optimal.fits',
                plot_prefix='extracted_optimal',
            )

            # Keep optimal extraction as default downstream product.
            extracted_spectra = extracted_optimal
            self.state.extracted_spectra = extracted_spectra
            self.state.extraction_done = True

            self._report_progress(0.80, "Spectrum Extraction Complete")
            logger.info(
                f"✓ Spectrum extraction stage complete (sum={extracted_sum.norders} orders, "
                f"optimal={extracted_optimal.norders} orders)"
            )

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
        logger.info("STEP 7: DE-BLAZING CORRECTION")
        logger.info("=" * 50)

        self._report_progress(0.90, "De-blazing Correction")

        try:
            de_blazed_spectra = process_de_blazing_stage(
                self.config, spectra_set, self.state.flat_field
            )

            self.state.de_blazed_spectra = de_blazed_spectra
            self.state.de_blazing_done = True

            self._report_progress(0.94, "De-blazing Complete")
            logger.info("✓ De-blazing stage complete")

            return de_blazed_spectra

        except Exception as e:
            logger.error(f"Error in de-blazing: {e}")
            raise

    def _refresh_flat_model_after_scattered_subtraction(self, flat_clean: np.ndarray):
        """Rebuild flat response/pixel-flat products using scattered-subtracted master flat."""
        if self.state.flat_field is None or self.state.apertures is None:
            raise RuntimeError("Flat/aperture state missing before Step 3 flat refresh")

        processor = FlatFieldProcessor(self.config)
        processor.flat_data = self.state.flat_field.flat_data
        processor.flat_mask = self.state.flat_field.flat_mask

        flat_sens, blaze_profiles, cross_profiles, smoothed_model, pixel_flat, _ = (
            processor.build_order_response_map(self.state.apertures, source_image=flat_clean)
        )

        self.state.flat_field.scattered_light = np.clip(
            self.state.flat_field.flat_data.astype(np.float32) - flat_clean.astype(np.float32),
            0.0,
            None,
        )
        self.state.flat_field.flat_sens = flat_sens
        self.state.flat_field.smoothed_model = smoothed_model
        self.state.flat_field.pixel_flat = pixel_flat
        self.state.flat_field.flat_corr_2d = processor.flat_corr_2d
        self.state.flat_field.blaze_profiles = blaze_profiles
        self.state.flat_field.cross_profiles = cross_profiles

        out_dir = Path(self.config.get_output_path()) / 'step2_scatterlight'
        out_dir.mkdir(parents=True, exist_ok=True)
        write_fits_image(str(out_dir / 'master_flat_clean.fits'), flat_clean.astype(np.float32), dtype='float32')

    def _save_step4_flat_products(self, out_dir: Path):
        """Persist Step4 flat-model artifacts (generated from current flat-field state)."""
        if self.state.flat_field is None:
            return

        out_dir.mkdir(parents=True, exist_ok=True)

        flat_state = self.state.flat_field
        if flat_state.blaze_profiles:
            with open(out_dir / 'blaze_profiles.pkl', 'wb') as f:
                pickle.dump(flat_state.blaze_profiles, f)

        if flat_state.flat_sens is not None:
            write_fits_image(
                str(out_dir / 'flat_sensitivity.fits'),
                flat_state.flat_sens.astype(np.float32),
                dtype='float32',
            )

        if flat_state.smoothed_model is not None:
            write_fits_image(
                str(out_dir / 'flat_smoothed_model.fits'),
                flat_state.smoothed_model.astype(np.float32),
                dtype='float32',
            )

        if flat_state.pixel_flat is not None:
            write_fits_image(
                str(out_dir / 'pixel_to_pixel_flat.fits'),
                flat_state.pixel_flat.astype(np.float32),
                dtype='float32',
            )

        if flat_state.flat_corr_2d is not None:
            write_fits_image(
                str(out_dir / 'flat_correction_2d.fits'),
                flat_state.flat_corr_2d.astype(np.float32),
                dtype='float32',
            )

        # Build response map on demand from the current clean flat + apertures.
        if self.state.apertures is not None:
            clean_flat = flat_state.flat_data
            if flat_state.scattered_light is not None:
                clean_flat = np.clip(
                    flat_state.flat_data.astype(np.float32) - flat_state.scattered_light.astype(np.float32),
                    1e-6,
                    None,
                )

            processor = FlatFieldProcessor(self.config)
            processor.flat_data = flat_state.flat_data
            processor.flat_mask = flat_state.flat_mask
            processor.build_order_response_map(self.state.apertures, source_image=clean_flat)
            if processor.response_map is not None:
                write_fits_image(
                    str(out_dir / 'flat_response_map.fits'),
                    processor.response_map.astype(np.float32),
                    dtype='float32',
                )

            clean_flat_corrected = None
            if flat_state.flat_corr_2d is not None:
                safe_flat = flat_state.flat_corr_2d.astype(np.float32)
                bad = (~np.isfinite(safe_flat)) | (safe_flat <= 0.05)
                if np.any(bad):
                    safe_flat = safe_flat.copy()
                    safe_flat[bad] = 1.0
                clean_flat_corrected = clean_flat.astype(np.float32) / safe_flat

            processor.save_step4_diagnostics(
                out_dir,
                clean_flat=clean_flat.astype(np.float32),
                clean_flat_corrected=clean_flat_corrected,
            )

    def stage_scattered_light_subtraction_science(self, science_image: np.ndarray,
                                                  science_name: str = 'science') -> np.ndarray:
        """Step 3: build/subtract scattered-light models for master flat and current science frame."""
        logger.info("=" * 50)
        logger.info("STEP 3: BACKGROUND SCATTERED-LIGHT MODELING AND SUBTRACTION")
        logger.info("=" * 50)

        self._report_progress(0.50, "Step 3: Master Flat Background")

        if self.state.flat_field is None or self.state.apertures is None:
            raise RuntimeError("Flat/apertures unavailable for Step 3")

        flat_background = process_background_stage(
            self.config,
            self.state.flat_field.flat_data,
            None,
            output_subdir='step2_scatterlight',
            output_tag='master_flat_',
            apertures=self.state.apertures,
            mask_margin_scale=1.2,
        )
        flat_clean = np.clip(
            self.state.flat_field.flat_data.astype(np.float32) - flat_background.astype(np.float32),
            1e-6,
            None,
        )
        self._refresh_flat_model_after_scattered_subtraction(flat_clean)

        self._report_progress(0.56, "Step 3: Science Background")
        science_background = process_background_stage(
            self.config,
            science_image,
            None,
            output_subdir='step2_scatterlight',
            output_tag=f'{science_name}_',
            apertures=self.state.apertures,
            mask_margin_scale=1.2,
        )

        corrected = science_image.astype(np.float32) - science_background.astype(np.float32)

        sci_out_dir = Path(self.config.get_output_path()) / 'step2_scatterlight'
        sci_out_dir.mkdir(parents=True, exist_ok=True)
        write_fits_image(
            str(sci_out_dir / f'{science_name}_science_background_subtracted.fits'),
            corrected.astype(np.float32),
            dtype='float32'
        )

        self._report_progress(0.60, "Science Scattered Light Subtraction Complete")
        return corrected

    def stage_flat_fielding_science_2d(self, science_image: np.ndarray,
                                       science_name: str = 'science') -> np.ndarray:
        """Step 4: apply 2D pixel-flat correction map to science image."""
        logger.info("=" * 50)
        logger.info("STEP 4: 2D FLAT FIELD CORRECTION")
        logger.info("=" * 50)

        self._report_progress(0.62, "2D Flat Fielding")

        if self.state.flat_field is None:
            raise RuntimeError("Flat model is unavailable. Run Step 2 first.")

        # Follow Step 4 definition: prefer pixel-to-pixel flat correction map.
        flat_corr = self.state.flat_field.pixel_flat
        if flat_corr is None:
            flat_corr = self.state.flat_field.flat_corr_2d
        if flat_corr is None:
            raise RuntimeError("No flat correction map available. Run Step 2 first.")

        if flat_corr.shape != science_image.shape:
            raise RuntimeError(
                f"Science image shape {science_image.shape} does not match flat correction shape {flat_corr.shape}"
            )

        safe_flat = flat_corr.astype(np.float32)
        bad = (~np.isfinite(safe_flat)) | (safe_flat <= 0.05)
        if np.any(bad):
            safe_flat = safe_flat.copy()
            safe_flat[bad] = 1.0

        corrected = science_image.astype(np.float32) / safe_flat

        out_dir = Path(self.config.get_output_path()) / 'step3_flat'
        out_dir.mkdir(parents=True, exist_ok=True)

        # Step4 owns all flat-model persistence artifacts.
        self._save_step4_flat_products(out_dir)

        write_fits_image(str(out_dir / f'{science_name}_science_flat2d_corrected.fits'), corrected, dtype='float32')

        if self.config.get_bool('reduce', 'save_plots', True):
            fig_format = self.config.get('reduce', 'fig_format', 'png')
            plot_2d_image_to_file(corrected, str(out_dir / f'{science_name}_science_flat2d_corrected.{fig_format}'), 'Science After 2D Flat Correction')
            # Shared correction map for current run, not duplicated per science image.
            plot_2d_image_to_file(safe_flat, str(out_dir / f'flat_correction_used.{fig_format}'), '2D Flat Correction Used')
            write_fits_image(str(out_dir / 'flat_correction_used.fits'), safe_flat.astype(np.float32), dtype='float32')

        self._report_progress(0.70, "2D Flat Fielding Complete")
        return corrected

    def stage_order_stitching(self, spectra_set: SpectraSet) -> tuple:
        """Step 8: stitch overlapping neighboring orders into one continuous 1D spectrum."""
        logger.info("=" * 50)
        logger.info("STEP 8: ORDER STITCHING")
        logger.info("=" * 50)

        self._report_progress(0.96, "Order Stitching")
        stitched = process_order_stitching_stage(self.config, spectra_set, output_subdir='step7_stitching')
        self._report_progress(1.00, "Pipeline Complete")
        return stitched

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
            # Step 1: Overscan subtraction & trimming for all inputs
            all_raw_files = [raw_image_path] + bias_filenames + flat_filenames + [calib_filename]
            corrected_files = self.stage_overscan_correction(all_raw_files)

            # Extract corrected filenames
            raw_image_path = corrected_files[0]
            bias_filenames = corrected_files[1:1+len(bias_filenames)]
            flat_start = 1 + len(bias_filenames)
            flat_end = flat_start + len(flat_filenames)
            flat_filenames = corrected_files[flat_start:flat_end]
            calib_filename = corrected_files[-1]

            # Step 2: Bias subtraction
            master_bias = self.stage_bias_correction(bias_filenames)

            # Apply bias to flat/calib first.
            flat_filenames = self._apply_master_bias_to_files(flat_filenames, master_bias, 'flat')
            calib_corrected = self._apply_master_bias_to_files([calib_filename], master_bias, 'calib')
            if calib_corrected:
                calib_filename = calib_corrected[0]

            # Step 2: Order tracing and flat model products
            flat_field, apertures = self.stage_flat_fielding(flat_filenames)

            # Read science image (already overscan-corrected)
            logger.info(f"Reading science image: {raw_image_path}")
            science_image, header = read_fits_image(raw_image_path)

            # Step 1 (continued): apply bias to science
            science_corrected = science_image.astype(float) - master_bias

            # Step 1 (continued): cosmic-ray removal in basic preprocessing.
            if self.config.get_bool('reduce', 'cosmic_enabled', True):
                logger.info("Step 1: Applying cosmic ray correction to science frame...")
                science_name = Path(raw_image_path).stem
                science_corrected = self.stage_cosmic_correction_science(science_corrected, science_name=science_name)

            pre_dir = Path(self.config.get_output_path()) / 'step1_basic' / 'preprocessed'
            pre_dir.mkdir(parents=True, exist_ok=True)
            write_fits_image(str(pre_dir / 'science_preprocessed.fits'), science_corrected, dtype='float32')

            # Step 3: master-flat and science scattered-light subtraction
            science_name = Path(raw_image_path).stem
            science_scattered_sub = self.stage_scattered_light_subtraction_science(
                science_corrected,
                science_name=science_name,
            )

            # Step 4: 2D flat-field correction
            science_flat2d = self.stage_flat_fielding_science_2d(
                science_scattered_sub,
                science_name=science_name,
            )

            # Step 5: 1D extraction (sum + optimal outputs; optimal used downstream)
            extracted_spectra = self.stage_extraction(science_flat2d)

            # Step 6: apply 2D wavelength solution λ=f(pixel, order) to extracted orders.
            calibrated_spectra = self.stage_wavelength_calibration_on_spectra(extracted_spectra, calib_filename)

            # Step 7: de-blaze using 1D flat/blaze model
            final_spectra = self.stage_de_blazing(calibrated_spectra)

            # Step 8: order stitching to continuous 1D spectrum
            self.stage_order_stitching(final_spectra)

            # Save final per-order spectra for archival
            base_output_path = self.config.get_output_path()
            spectra_dir = Path(base_output_path) / 'step8_final_spectra'
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
