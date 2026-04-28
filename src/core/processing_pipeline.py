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
from src.core.basic_reduction import (
    process_overscan_stage,
    process_bias_stage,
    process_cosmic_stage,
)
from src.core.order_tracing import process_order_tracing_stage
from src.core.flat_correction import process_flat_correction_stage
from src.core.wave_calibration import WavelengthCalibrator, process_wavelength_stage
from src.core.scattered_light import process_background_stage
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
            # --- Extract config for overscan stage ---
            overscan_kwargs = {
                'overscan_start_column': self.config.get_int('data', 'overscan_start_column', -1),
                'overscan_method': self.config.get('data', 'overscan_method', 'mean_only'),
                'overscan_smooth_window': self.config.get_int('data', 'overscan_smooth_window', -1),
                'overscan_poly_order': self.config.get_int('data', 'overscan_poly_order', 3),
                'overscan_poly_type': self.config.get('data', 'overscan_poly_type', 'legendre'),
                'trim_x_start': self.config.get_int('data', 'trim_x_start', -1),
                'trim_x_end': self.config.get_int('data', 'trim_x_end', -1),
                'trim_y_start': self.config.get_int('data', 'trim_y_start', -1),
                'trim_y_end': self.config.get_int('data', 'trim_y_end', -1),
                'save_plots': self.config.get_bool('reduce', 'save_plots', True),
                'fig_format': self.config.get('reduce', 'fig_format', 'png')
            }
            corrected_files = process_overscan_stage(
                raw_filenames, self.config.get_output_path(), **overscan_kwargs
            )

            self._report_progress(0.05, "Overscan Correction Complete")
            logger.info("✓ Overscan correction stage complete")

            return corrected_files

        except Exception as e:
            logger.error(f"Error in overscan correction: {e}")
            # Return original files if overscan correction fails
            return raw_filenames

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
            # --- Extract config for bias stage ---
            bias_kwargs = {
                'combine_method': self.config.get('reduce.bias', 'combine_method', 'median'),
                'combine_sigma': self.config.get_float('reduce.bias', 'combine_sigma', 3.0),
                'save_plots': self.config.get_bool('reduce', 'save_plots', True),
                'fig_format': self.config.get('reduce', 'fig_format', 'png'),
            }
            master_bias, bias_file = process_bias_stage(
                bias_filenames, self.config.get_output_path(), **bias_kwargs
            )

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
            # --- Extract config for flat / order tracing stage ---
            output_dir_base = self.config.get_output_path()
            flat_kwargs = {
                'output_dir_base': output_dir_base,
                'combine_method': self.config.get('reduce.flat', 'combine_method', 'median'),
                'combine_sigma': self.config.get_float('reduce.bias', 'combine_sigma', 3.0),
                'mosaic_maxcount': self.config.get_float('reduce.flat', 'mosaic_maxcount', 65535),
                'snr_threshold': self.config.get_float('reduce.trace', 'snr_threshold', 5.0),
                'gap_fill_factor': self.config.get_float('reduce.trace', 'gap_fill_factor', 1.6),
                'min_trace_coverage': self.config.get_float('reduce.trace', 'min_trace_coverage', 0.20),
                'trace_degree': self.config.get_int('reduce.trace', 'degree', 4),
                'width_cheb_degree': self.config.get_int('reduce.trace', 'width_cheb_degree', 3),
                'n_extend_below': self.config.get_int('reduce.trace', 'n_extend_below', 0),
                'n_extend_above': self.config.get_int('reduce.trace', 'n_extend_above', 0),
                'aperture_boundary_snr': self.config.get_float('reduce.trace', 'aperture_boundary_snr', 3.0),
                'gap_fill_factor_interp': self.config.get_float('reduce.trace', 'gap_fill_factor', 1.35),
                'save_plots': self.config.get_bool('reduce', 'save_plots', True),
                'fig_format': self.config.get('reduce', 'fig_format', 'png'),
            }
            flat_field, apertures = process_order_tracing_stage(
                flat_filenames, **flat_kwargs
            )

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
            logger.info(f"Extracting ThAr lamp spectrum from {calib_filename} ...")
            lamp_img, _ = read_fits_image(calib_filename)
            
            extract_kwargs = {
                'optimal_sigma': self.config.get_float('reduce.extract', 'optimal_sigma', 3.0),
                'extraction_method': 'sum',
                'save_plots': self.config.get_bool('reduce', 'save_plots', True),
                'fig_format': self.config.get('reduce', 'fig_format', 'png'),
            }
            lamp_spectra = process_extraction_stage(
                lamp_img, self.state.apertures,
                output_dir_base=self.config.get_output_path(),
                wavelength_calib=None,
                flat_field=self.state.flat_field,
                method_override='sum',
                output_filename='thar_1D_sum.fits',
                plot_prefix='thar_1D_sum',
                **extract_kwargs
            )
            
            wave_calib = process_wavelength_stage(
                lamp_spectra=lamp_spectra,
                config=self.config,
                output_dir_base=self.config.get_output_path(),
                lamp_type=self.config.get('telescope.linelist', 'linelist_type', 'ThAr'),
                save_plots=self.config.get_bool('reduce', 'save_plots', True),
                fig_format=self.config.get('reduce', 'fig_format', 'png'),
            )

            self.state.wavelength_calib = wave_calib
            self.state.wavelength_done = True

            self._report_progress(0.30, "Wavelength Calibration Complete")
            logger.info("✓ Wavelength calibration stage complete")

            return wave_calib

        except Exception as e:
            logger.error(f"Error in wavelength calibration: {e}")
            raise

    def stage_wavelength_calibration_on_spectra(self, spectra_set: SpectraSet,
                                                calib_filename: str,
                                                science_name: str = 'science') -> SpectraSet:
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
        logger.info("STEP 7: WAVELENGTH CALIBRATION (Science Spectra)")
        logger.info("=" * 50)

        self._report_progress(0.82, "Applying Wavelength Calibration")

        try:
            # First, calibrate the calibration frame if not already done
            if self.state.wavelength_calib is None:
                logger.info(f"Extracting ThAr lamp spectrum from {calib_filename} ...")
                lamp_img, _ = read_fits_image(calib_filename)
                
                extract_kwargs = {
                    'optimal_sigma': self.config.get_float('reduce.extract', 'optimal_sigma', 3.0),
                    'extraction_method': 'sum',
                    'save_plots': self.config.get_bool('reduce', 'save_plots', True),
                    'fig_format': self.config.get('reduce', 'fig_format', 'png'),
                }
                lamp_spectra = process_extraction_stage(
                    lamp_img, self.state.apertures,
                    output_dir_base=self.config.get_output_path(),
                    wavelength_calib=None,
                    flat_field=self.state.flat_field,
                    method_override='sum',
                    output_filename='thar_1D_sum.fits',
                    plot_prefix='thar_1D_sum',
                    **extract_kwargs
                )
                
                wave_calib = process_wavelength_stage(
                    lamp_spectra=lamp_spectra,
                    config=self.config,
                    output_dir_base=self.config.get_output_path(),
                    lamp_type=self.config.get('telescope.linelist', 'linelist_type', 'ThAr'),
                    save_plots=self.config.get_bool('reduce', 'save_plots', True),
                    fig_format=self.config.get('reduce', 'fig_format', 'png'),
                )
                self.state.wavelength_calib = wave_calib
            else:
                wave_calib = self.state.wavelength_calib

            # --- Apply Order Offset (Renumbering) ---
            delta_m = getattr(wave_calib, 'delta_m', 0)
            if delta_m != 0:
                logger.info(f"Applying order offset delta_m = {delta_m} to renumber all extracted data...")
                if self.state.apertures:
                    self.state.apertures.shift_orders(delta_m)
                if self.state.flat_field:
                    self.state.flat_field.shift_orders(delta_m)
                if self.state.extracted_spectra:
                    self.state.extracted_spectra.shift_orders(delta_m)
                if self.state.de_blazed_spectra:
                    self.state.de_blazed_spectra.shift_orders(delta_m)
                if spectra_set is not self.state.de_blazed_spectra:
                    spectra_set.shift_orders(delta_m)

            # Apply wavelength calibration to each extracted spectrum
            calibrator = WavelengthCalibrator()
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

            # === 将波长定标后的科学光谱独立保存到 Step 7 对应目录 ===
            if self.config.get_bool('reduce.save_intermediate', 'save_wlcalib', True):
                out_dir = Path(self.config.get_output_path()) / 'step7_wavelength'
                out_dir.mkdir(parents=True, exist_ok=True)
                
                from src.core.de_blazing import save_deblazed_spectra
                save_deblazed_spectra(str(out_dir / f'{science_name}_1D_calibrated.fits'), calibrated_spectra)
                
                if self.config.get_bool('reduce', 'save_plots', True):
                    from src.plotting.spectra_plotter import plot_spectra_to_pdf
                    pdf_path = out_dir / f"{science_name}_1D_calibrated.pdf"
                    plot_spectra_to_pdf(calibrated_spectra, str(pdf_path), 
                                        title_prefix="Wavelength Calibrated Spectrum", xlabel=r"Wavelength ($\AA$)")

            self._report_progress(0.88, "Wavelength Calibration Complete")
            logger.info(f"✓ Wavelength calibration applied to {len(calibrated_spectra.spectra)} spectra")

            return calibrated_spectra

        except Exception as e:
            logger.error(f"Error in wavelength calibration of spectra: {e}")
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
        logger.info("STEP 1 (sub): COSMIC RAY CORRECTION (Science Frame)")
        logger.info("=" * 50)

        self._report_progress(0.45, "Cosmic Ray Correction")

        try:
            # Reuse the shared Step 1 cosmic detector so selected-step execution
            # and full-pipeline execution behave consistently.
            cosmic_kwargs = {
                'cosmic_sigclip': self.config.get_float('reduce', 'cosmic_sigclip', 5.0),
                'cosmic_objlim': self.config.get_float('reduce', 'cosmic_objlim', 5.0),
                'cosmic_gain': self.config.get_float('reduce', 'cosmic_gain', 1.0),
                'cosmic_readnoise': self.config.get_float('reduce', 'cosmic_readnoise', 5.0),
                'save_plots': self.config.get_bool('reduce', 'save_plots', True),
                'fig_format': self.config.get('reduce', 'fig_format', 'png'),
            }

            out_dir = Path(self.config.get_output_path()) / 'step1_basic' / 'cosmic_corrected'
            out_dir.mkdir(parents=True, exist_ok=True)
            
            # Create a temporary file to pass to process_cosmic_stage
            temp_input_path = out_dir / f"{science_name}_precosmic.fits"
            write_fits_image(str(temp_input_path), science_image, dtype='float32')

            corrected_files = process_cosmic_stage(
                [str(temp_input_path)],
                output_dir_base=self.config.get_output_path(),
                **cosmic_kwargs
            )
            
            corrected_image, _ = read_fits_image(corrected_files[0])
            temp_input_path.unlink() # Clean up temporary file

            if self.config.get_bool('reduce', 'save_plots', True):
                fig_format = self.config.get('reduce', 'fig_format', 'png')
                plot_2d_image_to_file(
                    corrected_image,
                    str(out_dir / f'{science_name}_science_cosmic_corrected.{fig_format}'),
                    'Science After Cosmic Ray Correction'
                )

            self._report_progress(0.48, "Cosmic Ray Correction Complete")
            logger.info("✓ Cosmic ray correction complete")

            return corrected_image

        except Exception as e:
            logger.error(f"Error in cosmic ray correction: {e}")
            raise

    def stage_extraction(self, science_image: np.ndarray, science_name: str = 'science') -> SpectraSet:
        """
        Execute spectrum extraction stage and output both sum/optimal products.

        Args:
            science_image: Background-corrected 2D science image
            science_name: Base name of the science image

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
            # --- Extract config for extraction stage ---
            output_dir_base = self.config.get_output_path()
            extract_kwargs = {
                'output_dir_base': output_dir_base,
                'optimal_sigma': self.config.get_float('reduce.extract', 'optimal_sigma', 3.0),
                'extraction_method': self.config.get('reduce.extract', 'method', 'optimal'),
                'save_plots': self.config.get_bool('reduce', 'save_plots', True),
                'fig_format': self.config.get('reduce', 'fig_format', 'png'),
            }
            # Step 5 requirement: output both aperture-sum and optimal extraction.
            extracted_sum = process_extraction_stage(
                science_image, self.state.apertures,
                wavelength_calib=None,
                flat_field=self.state.flat_field,
                method_override='sum',
                output_filename=f'{science_name}_1D_sum.fits',
                plot_prefix=f'{science_name}_1D_sum',
                **extract_kwargs,
            )

            extracted_optimal = process_extraction_stage(
                science_image, self.state.apertures,
                wavelength_calib=None,
                flat_field=self.state.flat_field,
                method_override='optimal',
                output_filename=f'{science_name}_1D_optimal.fits',
                plot_prefix=f'{science_name}_1D_optimal',
                **extract_kwargs,
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

    def stage_de_blazing(self, spectra_set: SpectraSet, science_name: str = 'science') -> SpectraSet:
        """
        Execute de-blazing correction stage.

        Args:
            spectra_set: Extracted spectra
            science_name: Base name of the science image

        Returns:
            De-blazed spectra
        """
        logger.info("=" * 50)
        logger.info("STEP 6: DE-BLAZING (Blaze Function Correction)")
        logger.info("=" * 50)

        self._report_progress(0.90, "De-blazing Correction")

        try:
            # --- Extract config for de-blazing stage ---
            de_blazed_spectra = process_de_blazing_stage(
                spectra_set,
                output_dir_base=self.config.get_output_path(),
                flat_field=self.state.flat_field,
                save_deblaze=self.config.get_bool('reduce.save_intermediate', 'save_deblaze', True),
                output_filename=f'{science_name}_1D_Deblaze.fits',
                plot_prefix=f'{science_name}_1D_Deblaze',
                save_plots=self.config.get_bool('reduce', 'save_plots', True),
                fig_format=self.config.get('reduce', 'fig_format', 'png'),
            )

            self.state.de_blazed_spectra = de_blazed_spectra
            self.state.de_blazing_done = True

            self._report_progress(0.94, "De-blazing Complete")
            logger.info("✓ De-blazing stage complete")

            return de_blazed_spectra

        except Exception as e:
            logger.error(f"Error in de-blazing: {e}")
            raise

    def _save_step4_flat_products(self, out_dir: Path):
        """Delegate Step-4 flat-model persistence to flat_correction module."""
        if self.state.flat_field is None:
            return
        from src.core.flat_correction import save_flat_correction_products
        save_flat_correction_products(
            self.state.flat_field, self.state.apertures, out_dir,
            fig_format=self.config.get('reduce', 'fig_format', 'png'),
            save_plots=self.config.get_bool('reduce', 'save_plots', True),
        )

    def stage_scattered_light_subtraction_science(self, science_image: np.ndarray,
                                                  science_name: str = 'science') -> np.ndarray:
        """Step 3: build/subtract scattered-light models for master flat and current science frame."""
        logger.info("=" * 50)
        logger.info("STEP 3: SCATTERED LIGHT MODELING AND SUBTRACTION")
        logger.info("=" * 50)

        self._report_progress(0.50, "Step 3: Master Flat Scattered Light")

        if self.state.flat_field is None or self.state.apertures is None:
            raise RuntimeError("Flat/apertures unavailable for Step 3")

        # --- Extract config for background stage ---
        output_dir_base = self.config.get_output_path()
        # 读取用户设置的mask参数
        mask_margin_pixels = self.config.get_int('reduce.background', 'mask_margin_pixels', 1)
        n_mask_below = self.config.get_int('reduce.trace', 'n_mask_below', 4)
        n_mask_above = self.config.get_int('reduce.trace', 'n_mask_above', 4)
        
        out_dir = Path(output_dir_base) / 'step3_scatterlight'
        out_dir.mkdir(parents=True, exist_ok=True)
        
        from src.core.scattered_light import create_widened_mask
        h, w = self.state.flat_field.flat_data.shape
        
        step3_mask, lo_traces, hi_traces, full_ids = create_widened_mask(
            self.state.apertures, self.state.flat_field.flat_data, mask_margin_pixels,
            n_mask_below, n_mask_above
        )
        
        write_fits_image(str(out_dir / 'Scattered_Light_Mask.fits'), step3_mask, dtype='uint8')
        
        if self.config.get_bool('reduce', 'save_plots', True):
            import matplotlib.pyplot as plt
            fig_format = self.config.get('reduce', 'fig_format', 'png')
            plt.figure(figsize=(12, 10))
            plt.imshow(step3_mask, aspect='auto', origin='lower', cmap='gray', vmin=0, vmax=1)
            
            x_anno = w - 1
            for i, ap_id in enumerate(full_ids):
                yc = 0.5 * (lo_traces[i, x_anno] + hi_traces[i, x_anno])
                if np.isfinite(yc) and 0 <= yc < h:
                    plt.text(x_anno + 5, yc, str(ap_id), color='red', fontsize=3.5, ha='left', va='center', clip_on=False)
                    
            plt.title(f'Step 3 Widened Mask (margin={mask_margin_pixels})\nWhite=Masked (Orders + Virtual), Black=Background')
            plt.xlabel('Pixel (X)')
            plt.ylabel('Pixel (Y)')
            plt.xlim(0, w - 1)
            plt.ylim(0, h - 1)
            plot_file = out_dir / f'scattered_light_mask.{fig_format}'
            plt.savefig(str(plot_file), dpi=150, bbox_inches='tight')
            plt.close()
        
        bg_kwargs = {
            'output_subdir': 'step3_scatterlight',
            'mask_margin_scale': 1.0,
            'mask_margin_pixels': 0, # Already widened by us
            'order_mask': step3_mask,
            'poly_order': self.config.get_int('reduce.background', 'poly_order', 3),
            'bg_method': self.config.get('reduce.background', 'method', 'convolution'),
            'sigma_clip_val': self.config.get_float('reduce.background', 'sigma_clip', 3.0),
            'maxiters': self.config.get_int('reduce.background', 'sigma_clip_maxiters', 4),
            'bspline_smooth': self.config.get_float('reduce.background', 'bspline_smooth', 1.0),
            'n_mask_below': n_mask_below,
            'n_mask_above': n_mask_above,
            'clip_mode': self.config.get('reduce.background', 'sigma_clip_mode', 'upper'),
            'split_row': self.config.get_int('data', 'detector_split_row', 2068),
            'kernel_sigma_x': self.config.get_float('reduce.background', 'kernel_sigma_x', 13.0),
            'kernel_sigma_y': self.config.get_float('reduce.background', 'kernel_sigma_y', 13.0),
            'spline_smooth_factor': self.config.get_float('reduce.background', 'spline_smooth_factor', 1.0),
            'spline_post_smooth_x': self.config.get_float('reduce.background', 'spline_post_smooth_x', 5.0),
            'save_plots': self.config.get_bool('reduce', 'save_plots', True),
            'fig_format': self.config.get('reduce', 'fig_format', 'png'),
        }
        flat_background = process_background_stage(
            self.state.flat_field.flat_data,
            output_dir_base,
            apertures=self.state.apertures,
            output_tag='MasterFlat_',
            **bg_kwargs,
        )
        flat_clean = np.clip(
            self.state.flat_field.flat_data.astype(np.float32) - flat_background.astype(np.float32),
            1e-6,
            None,
        )
        self.state.flat_field.scattered_light = flat_background
        
        master_flat_path = Path(output_dir_base) / 'step2_trace' / 'MasterFlat.fits'
        _, flat_header = read_fits_image(str(master_flat_path)) if master_flat_path.exists() else (None, None)
        if flat_header is not None:
            flat_header['BKGSCAT'] = (True, 'Scattered light background subtracted')
            
        write_fits_image(str(out_dir / 'MasterFlat.fits'), flat_clean.astype(np.float32), header=flat_header, dtype='float32')

        self._report_progress(0.56, "Step 3: Science Scattered Light")
        science_background = process_background_stage(
            science_image,
            output_dir_base,
            output_tag=f'{science_name}_',
            apertures=self.state.apertures,
            **bg_kwargs,
        )

        corrected = science_image.astype(np.float32) - science_background.astype(np.float32)

        # The corrected image with original filename is already saved by process_background_stage.

        self._report_progress(0.60, "Science Scattered Light Subtraction Complete")
        return corrected

    def stage_flat_fielding_science_2d(self, science_image: np.ndarray,
                                       science_name: str = 'science') -> np.ndarray:
        """Step 4: apply 2D pixel-flat correction map to science image."""
        self._report_progress(0.62, "2D Flat Fielding")
        # --- Extract config for flat correction stage ---
        flat_corr_kwargs = {
            'output_dir_base': self.config.get_output_path(),
            'blaze_knot_spacing': self.config.get_int('reduce.flat', 'blaze_knot_spacing', 500),
            'blaze_edge_nknots': self.config.get_int('reduce.flat', 'blaze_edge_nknots', 6),
            'width_smooth_window': self.config.get_int('reduce.flat', 'width_smooth_window', 41),
            'profile_bin_step': self.config.get_float('reduce.flat', 'profile_bin_step', 0.01),
            'n_profile_segments': self.config.get_int('reduce.flat', 'n_profile_segments', 100),
            'profile_smooth_sigma': self.config.get_float('reduce.flat', 'profile_smooth_sigma', 6.0),
            'pixel_flat_min': self.config.get_float('reduce.flat', 'pixel_flat_min', 0.5),
            'pixel_flat_max': self.config.get_float('reduce.flat', 'pixel_flat_max', 1.5),
            'save_plots': self.config.get_bool('reduce', 'save_plots', True),
            'fig_format': self.config.get('reduce', 'fig_format', 'png'),
        }
        corrected = process_flat_correction_stage(
            science_image, self.state.flat_field,
            apertures=self.state.apertures, science_name=science_name,
            **flat_corr_kwargs,
        )
        self._report_progress(0.70, "2D Flat Fielding Complete")
        return corrected

    def stage_order_stitching(self, spectra_set: SpectraSet) -> tuple:
        """Step 8: stitch overlapping neighboring orders into one continuous 1D spectrum."""
        logger.info("=" * 50)
        logger.info("STEP 8: ORDER STITCHING")
        logger.info("=" * 50)

        self._report_progress(0.96, "Order Stitching")
        stitched = process_order_stitching_stage(
            spectra_set,
            output_dir_base=self.config.get_output_path(),
            output_subdir='step8_stitching',
            save_plots=self.config.get_bool('reduce', 'save_plots', True),
            fig_format=self.config.get('reduce', 'fig_format', 'png'),
        )
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
            raw_image_path_corr = corrected_files[0]
            bias_filenames_corr = corrected_files[1:1+len(bias_filenames)]
            flat_start = 1 + len(bias_filenames)
            flat_end = flat_start + len(flat_filenames)
            flat_filenames_corr = corrected_files[flat_start:flat_end]
            calib_filename_corr = corrected_files[-1]

            # Step 2: Bias subtraction
            master_bias = self.stage_bias_correction(bias_filenames_corr)

            # Apply bias to flat/calib first.
            flat_filenames_bias = self._apply_master_bias_to_files(flat_filenames_corr, master_bias, 'flat')
            calib_corrected = self._apply_master_bias_to_files([calib_filename_corr], master_bias, 'calib')
            if calib_corrected:
                calib_filename_corr = calib_corrected[0]

            # Determine which flats to use for master flat:
            # If overscan/bias correction was run (i.e. flat_filenames_bias exists and is not empty), use those;
            # otherwise, fall back to raw flat_filenames.
            if flat_filenames_bias and all(Path(f).exists() for f in flat_filenames_bias):
                master_flat_files = flat_filenames_bias
            else:
                master_flat_files = flat_filenames

            # Step 2: Order tracing and flat model products
            flat_field, apertures = self.stage_flat_fielding(master_flat_files)

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
            extracted_spectra = self.stage_extraction(science_flat2d, science_name=science_name)

            # Step 6: de-blaze using 1D flat/blaze model (Blaze Function Correction)
            deblazed_spectra = self.stage_de_blazing(extracted_spectra, science_name=science_name)

            # Step 7: apply 2D wavelength solution λ=f(pixel, order) to de-blazed orders.
            final_spectra = self.stage_wavelength_calibration_on_spectra(
                deblazed_spectra, calib_filename_corr if 'calib_filename_corr' in locals() else calib_filename, science_name=science_name)

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
