"""
Main processing pipeline orchestrator.

Coordinates all processing stages and manages workflow execution.
"""

import logging
import pickle
from pathlib import Path
from typing import List, Callable, Optional, Union, Dict
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
            calib_name = Path(calib_filename).stem
            
            extract_kwargs = {
                'optimal_sigma': self.config.get_float('reduce.extract', 'optimal_sigma', 3.0),
                'save_plots': self.config.get_bool('reduce', 'save_plots', True),
                'fig_format': self.config.get('reduce', 'fig_format', 'png'),
            }
            lamp_spectra = process_extraction_stage(
                lamp_img, self.state.apertures,
                output_dir_base=self.config.get_output_path(),
                wavelength_calib=None,
                flat_field=self.state.flat_field,
                    method_override='optimal',
                    output_filename=f'{calib_name}_1D_optimal.fits',
                    plot_prefix=f'{calib_name}_1D_optimal',
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
                                                wave_calib: WaveCalib,
                                                science_name: str = 'science',
                                                method: str = 'optimal') -> SpectraSet:
        """
        Apply wavelength calibration to extracted science spectra.

        This is the correct order: extract 1D spectra first, then apply wavelength calibration.
        Following the gamse standard pipeline.

        Args:
            spectra_set: Extracted spectra without wavelength information
            wave_calib: Wavelength calibration model to apply

        Returns:
            SpectraSet with wavelength-calibrated spectra
        """
        logger.info("=" * 50)
        logger.info(f"STEP 7: WAVELENGTH CALIBRATION ({science_name} - {method})")
        logger.info("=" * 50)

        self._report_progress(0.82, "Applying Wavelength Calibration")

        try:
            # --- Apply Order Offset (Renumbering) ---
            delta_m = getattr(wave_calib, 'delta_m', 0)
            if delta_m != 0:
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
                save_deblazed_spectra(str(out_dir / f'{science_name}_1D_{method}_calibrated.fits'), calibrated_spectra)
                
                if self.config.get_bool('reduce', 'save_plots', True):
                    from src.plotting.spectra_plotter import plot_spectra_to_pdf
                    pdf_path = out_dir / f"{science_name}_1D_{method}_calibrated.pdf"
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

    def stage_extraction(self, science_image: np.ndarray, science_name: str = 'science',
                         calib_image_paths: Optional[List[str]] = None,
                         master_flat_path: Optional[str] = None) -> SpectraSet:
        """
        Execute spectrum extraction stage and output both sum/optimal products.

        Args:
            science_image: Background-corrected 2D science image
            science_name: Base name of the science image
            calib_image_paths: List of paths to calibration images
            master_flat_path: Path to the MasterFlat image

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
                'save_plots': self.config.get_bool('reduce', 'save_plots', True),
                'fig_format': self.config.get('reduce', 'fig_format', 'png'),
            }
            
            # 1. Science Extraction
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

            # 2. MasterFlat Extraction
            if master_flat_path and Path(master_flat_path).exists():
                try:
                    mf_img, _ = read_fits_image(master_flat_path)
                    mf_name = "MasterFlat"
                    process_extraction_stage(
                        mf_img, self.state.apertures,
                        wavelength_calib=None, flat_field=self.state.flat_field,
                        method_override='sum', output_filename=f'{mf_name}_1D_sum.fits',
                        plot_prefix=f'{mf_name}_1D_sum', **extract_kwargs
                    )
                    process_extraction_stage(
                        mf_img, self.state.apertures,
                        wavelength_calib=None, flat_field=self.state.flat_field,
                        method_override='optimal', output_filename=f'{mf_name}_1D_optimal.fits',
                        plot_prefix=f'{mf_name}_1D_optimal', **extract_kwargs
                    )
                    logger.info(f"✓ Extracted MasterFlat 1D spectra")
                except Exception as e:
                    logger.warning(f"Failed to extract MasterFlat: {e}")

            # 3. Calibration Extraction
            if calib_image_paths:
                for calib_image_path in calib_image_paths:
                    if Path(calib_image_path).exists():
                        try:
                            c_img, _ = read_fits_image(calib_image_path)
                            c_name = Path(calib_image_path).stem
                            process_extraction_stage(
                                c_img, self.state.apertures,
                                wavelength_calib=None, flat_field=self.state.flat_field,
                                method_override='sum', output_filename=f'{c_name}_1D_sum.fits',
                                plot_prefix=f'{c_name}_1D_sum', **extract_kwargs
                            )
                            process_extraction_stage(
                                c_img, self.state.apertures,
                                wavelength_calib=None, flat_field=self.state.flat_field,
                                method_override='optimal', output_filename=f'{c_name}_1D_optimal.fits',
                                plot_prefix=f'{c_name}_1D_optimal', **extract_kwargs
                            )
                            logger.info(f"✓ Extracted Calibration 1D spectra for {c_name}")
                        except Exception as e:
                            logger.warning(f"Failed to extract Calibration {Path(calib_image_path).name}: {e}")

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
                output_filename=f'{science_name}_1D_optimal_Deblaze.fits',
                plot_prefix=f'{science_name}_1D_optimal_Deblaze',
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
            'blaze_smooth_factor': self.config.get_float('reduce.flat', 'blaze_smooth_factor', 1.0),
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

    def stage_flat_fielding_calib_2d(self, calib_file: str) -> str:
        """Apply 2D pixel-flat correction map to calibration image."""
        self._report_progress(0.72, "2D Flat Fielding (Calibration)")
        
        calib_name = Path(calib_file).stem
        orig_name = Path(calib_file).name
        out_dir = Path(self.config.get_output_path()) / 'step4_flat_corrected'
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / orig_name
        
        # Priority: bias_subtracted -> overscan_corrected -> raw
        input_calib = calib_file
        bias_path = Path(self.config.get_output_path()) / 'step1_basic' / 'bias_subtracted' / orig_name
        overscan_path = Path(self.config.get_output_path()) / 'step1_basic' / 'overscan_corrected' / orig_name
        
        if bias_path.exists():
            input_calib = str(bias_path)
        elif overscan_path.exists():
            input_calib = str(overscan_path)
            
        logger.info(f"Applying 2D flat correction to calibration frame: {input_calib}")
        
        if self.state.flat_field is None:
            logger.warning("No flat field available, skipping calib 2D flat correction")
            return input_calib
            
        flat_corr = self.state.flat_field.pixel_flat
        if flat_corr is None:
            flat_corr = self.state.flat_field.flat_corr_2d
            
        if flat_corr is None:
            logger.warning("No flat correction map available in flat_field. Skipping calib 2D flat correction.")
            return input_calib
            
        calib_img, header = read_fits_image(input_calib)
            
        safe_flat = flat_corr.astype(np.float32)
        bad = (~np.isfinite(safe_flat)) | (safe_flat <= 0.05)
        if np.any(bad):
            safe_flat = safe_flat.copy()
            safe_flat[bad] = 1.0
            
        corrected = calib_img.astype(np.float32) / safe_flat
        if header is not None:
            header['FLATCOR'] = (True, '2D pixel flat correction applied')
            
        write_fits_image(str(out_path), corrected, header=header, dtype='float32')
        
        if self.config.get_bool('reduce', 'save_plots', True):
            from src.plotting.spectra_plotter import plot_2d_image_to_file
            fig_format = self.config.get('reduce', 'fig_format', 'png')
            plot_2d_image_to_file(corrected, str(out_dir / f"{calib_name}_flat2d_corrected.{fig_format}"), "Calibration Frame After 2D Flat Correction")
            
        return str(out_path)

    def get_best_calib_file(self, raw_calib_path: str) -> str:
        """Find the most processed version of the calibration file."""
        if not raw_calib_path:
            return raw_calib_path
        calib_name = Path(raw_calib_path).stem
        orig_name = Path(raw_calib_path).name
        out_dir = Path(self.config.get_output_path())
        
        flat_path = out_dir / 'step4_flat_corrected' / orig_name
        if flat_path.exists():
            return str(flat_path)
            
        bias_path = out_dir / 'step1_basic' / 'bias_subtracted' / orig_name
        if bias_path.exists():
            return str(bias_path)
            
        overscan_path = out_dir / 'step1_basic' / 'overscan_corrected' / orig_name
        if overscan_path.exists():
            return str(overscan_path)
            
        return raw_calib_path

    def run_full_pipeline(self, raw_image_paths: Union[str, List[str]],
                         bias_filenames: List[str],
                         flat_filenames: List[str],
                         calib_filenames: List[str]) -> Dict[str, SpectraSet]:
        """
        Execute complete spectral reduction pipeline.

        Args:
            raw_image_paths: List of paths to science FITS files
            bias_filenames: List of bias frame paths
            flat_filenames: List of flat frame paths
            calib_filenames: List of wavelength calibration frame paths

        Returns:
            Dict mapping science filenames to SpectraSet with final extracted spectra
        """
        if isinstance(raw_image_paths, str):
            raw_image_paths = [raw_image_paths]

        logger.info("=" * 60)
        logger.info("STARTING FULL SPECTRAL REDUCTION PIPELINE")
        logger.info("=" * 60)

        # Reset state
        self.state.reset()

        # Create output directories
        self.config.create_directories()

        try:
            # Step 1: Overscan subtraction & trimming for all inputs
            all_raw_files = [raw_image_path] + bias_filenames + flat_filenames + calib_filenames
            corrected_files = self.stage_overscan_correction(all_raw_files)

            # Extract corrected filenames
            raw_image_path_corr = corrected_files[0]
            bias_filenames_corr = corrected_files[1:1+len(bias_filenames)]
            flat_start = 1 + len(bias_filenames)
            flat_end = flat_start + len(flat_filenames)
            flat_filenames_corr = corrected_files[flat_start:flat_end]
            calib_filenames_corr = corrected_files[flat_end:]

            # Step 2: Bias subtraction
            master_bias = self.stage_bias_correction(bias_filenames_corr)

            # Apply bias to flat/calib first.
            flat_filenames_bias = self._apply_master_bias_to_files(flat_filenames_corr, master_bias, 'flat')
            calib_corrected = self._apply_master_bias_to_files(calib_filenames_corr, master_bias, 'calib')
            if calib_corrected:
                calib_filenames_corr = calib_corrected

            # Determine which flats to use for master flat:
            # If overscan/bias correction was run (i.e. flat_filenames_bias exists and is not empty), use those;
            # otherwise, fall back to raw flat_filenames.
            if flat_filenames_bias and all(Path(f).exists() for f in flat_filenames_bias):
                master_flat_files = flat_filenames_bias
            else:
                master_flat_files = flat_filenames

            # Step 2: Order tracing and flat model products
            flat_field, apertures = self.stage_flat_fielding(master_flat_files)

            # Step 3: Scattered light subtraction
            self._report_progress(0.40, "Scattered Light Subtraction")
            logger.info("=" * 50)
            logger.info("STEP 3: SCATTERED LIGHT MODELING AND SUBTRACTION")
            logger.info("=" * 50)
            
            out_dir = Path(self.config.get_output_path()) / 'step3_scatterlight'
            out_dir.mkdir(parents=True, exist_ok=True)
            from src.core.scattered_light import create_widened_mask, process_background_stage
            
            mask_margin_pixels = self.config.get_int('reduce.background', 'mask_margin_pixels', 1)
            n_mask_below = self.config.get_int('reduce.trace', 'n_mask_below', 4)
            n_mask_above = self.config.get_int('reduce.trace', 'n_mask_above', 4)
            
            step3_mask, _, _, _ = create_widened_mask(
                apertures, flat_field.flat_data, mask_margin_pixels, n_mask_below, n_mask_above
            )
            write_fits_image(str(out_dir / 'Scattered_Light_Mask.fits'), step3_mask, dtype='uint8')
            
            bg_kwargs = {
                'output_subdir': 'step3_scatterlight',
                'mask_margin_scale': 1.0,
                'mask_margin_pixels': 0,
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
            
            self._report_progress(0.45, "Step 3: Master Flat Scattered Light")
            flat_background = process_background_stage(
                flat_field.flat_data, self.config.get_output_path(),
                apertures=apertures, output_tag='MasterFlat_', **bg_kwargs
            )
            flat_clean = np.clip(
                flat_field.flat_data.astype(np.float32) - flat_background.astype(np.float32), 1e-6, None
            )
            flat_field.scattered_light = flat_background
            
            master_flat_path = Path(self.config.get_output_path()) / 'step2_trace' / 'MasterFlat.fits'
            _, flat_header = read_fits_image(str(master_flat_path)) if master_flat_path.exists() else (None, None)
            if flat_header is not None:
                flat_header['BKGSCAT'] = (True, 'Scattered light background subtracted')
            write_fits_image(str(out_dir / 'MasterFlat.fits'), flat_clean.astype(np.float32), header=flat_header, dtype='float32')

            self._report_progress(0.50, "Step 3: Science Scattered Light")
            science_scatter_files = []
            for sci_file in science_corrected_files:
                sci_img, sci_header = read_fits_image(sci_file)
                sci_name = Path(sci_file).stem
                sci_bg = process_background_stage(
                    sci_img, self.config.get_output_path(),
                    output_tag=f'{sci_name}_', apertures=apertures, **bg_kwargs
                )
                sci_clean = sci_img.astype(np.float32) - sci_bg.astype(np.float32)
                if sci_header is not None:
                    sci_header['BKGSCAT'] = (True, 'Scattered light background subtracted')
                out_path = out_dir / Path(sci_file).name
                write_fits_image(str(out_path), sci_clean, header=sci_header, dtype='float32')
                science_scatter_files.append(str(out_path))

            # Step 4: 2D flat-field correction
            self._report_progress(0.60, "2D Flat Fielding")
            
            # Step 4: 2D flat correction for calib files
            calib_filenames_final = []
            for c_file in calib_filenames_corr:
                calib_filenames_final.append(self.stage_flat_fielding_calib_2d(c_file))

            science_flat2d_files = []
            for sci_file in science_scatter_files:
                sci_img, sci_header = read_fits_image(sci_file)
                sci_name = Path(sci_file).stem
                corr = self.stage_flat_fielding_science_2d(sci_img, science_name=sci_name)
                
                out_dir_f = Path(self.config.get_output_path()) / 'step4_flat_corrected'
                out_path = out_dir_f / Path(sci_file).name
                if sci_header is not None:
                    sci_header['FLATCOR'] = (True, '2D pixel flat correction applied')
                write_fits_image(str(out_path), corr, header=sci_header, dtype='float32')
                science_flat2d_files.append(str(out_path))

            # Step 5: 1D Extraction
            self._report_progress(0.70, "Spectrum Extraction")
            logger.info("=" * 50)
            logger.info("STEP 5: SPECTRUM EXTRACTION")
            logger.info("=" * 50)
            extract_kwargs = {
                'output_dir_base': self.config.get_output_path(),
                'optimal_sigma': self.config.get_float('reduce.extract', 'optimal_sigma', 3.0),
                'save_plots': self.config.get_bool('reduce', 'save_plots', True),
                'fig_format': self.config.get('reduce', 'fig_format', 'png'),
            }
            
            mf_path = Path(self.config.get_output_path()) / 'step4_flat_corrected' / 'MasterFlat.fits'
            if mf_path.exists():
                mf_img, _ = read_fits_image(str(mf_path))
                process_extraction_stage(mf_img, apertures, wavelength_calib=None, flat_field=flat_field, method_override='sum', output_filename='MasterFlat_1D_sum.fits', plot_prefix='MasterFlat_1D_sum', **extract_kwargs)
                process_extraction_stage(mf_img, apertures, wavelength_calib=None, flat_field=flat_field, method_override='optimal', output_filename='MasterFlat_1D_optimal.fits', plot_prefix='MasterFlat_1D_optimal', **extract_kwargs)

            for c_file in calib_filenames_final:
                c_img, _ = read_fits_image(c_file)
                c_name = Path(c_file).stem
                process_extraction_stage(c_img, apertures, wavelength_calib=None, flat_field=flat_field, method_override='sum', output_filename=f'{c_name}_1D_sum.fits', plot_prefix=f'{c_name}_1D_sum', **extract_kwargs)
                process_extraction_stage(c_img, apertures, wavelength_calib=None, flat_field=flat_field, method_override='optimal', output_filename=f'{c_name}_1D_optimal.fits', plot_prefix=f'{c_name}_1D_optimal', **extract_kwargs)

            extracted_science_dict = {}
            for sci_file in science_flat2d_files:
                sci_img, _ = read_fits_image(sci_file)
                sci_name = Path(sci_file).stem
                process_extraction_stage(sci_img, apertures, wavelength_calib=None, flat_field=flat_field, method_override='sum', output_filename=f'{sci_name}_1D_sum.fits', plot_prefix=f'{sci_name}_1D_sum', **extract_kwargs)
                opt_spec = process_extraction_stage(sci_img, apertures, wavelength_calib=None, flat_field=flat_field, method_override='optimal', output_filename=f'{sci_name}_1D_optimal.fits', plot_prefix=f'{sci_name}_1D_optimal', **extract_kwargs)
                extracted_science_dict[sci_name] = opt_spec
            
            # Step 6: De-blazing
            self._report_progress(0.80, "De-blazing")
            logger.info("=" * 50)
            logger.info("STEP 6: DE-BLAZING")
            logger.info("=" * 50)
            def do_deblaze(target_name, spectra_obj=None):
                from src.core.extraction import load_extracted_spectra
                db_opt = None
                for method in ['sum', 'optimal']:
                    spectra = spectra_obj if (method == 'optimal' and spectra_obj is not None) else None
                    if spectra is None:
                        try: spectra = load_extracted_spectra(str(Path(self.config.get_output_path()) / 'step5_extraction' / f'{target_name}_1D_{method}.fits'))
                        except Exception: pass
                    if spectra is not None:
                        db_spec = process_de_blazing_stage(
                            spectra, self.config.get_output_path(), flat_field=flat_field,
                            save_deblaze=self.config.get_bool('reduce.save_intermediate', 'save_deblaze', True),
                            output_filename=f'{target_name}_1D_{method}_Deblaze.fits',
                            plot_prefix=f'{target_name}_1D_{method}_Deblaze',
                            save_plots=self.config.get_bool('reduce', 'save_plots', True),
                            fig_format=self.config.get('reduce', 'fig_format', 'png')
                        )
                        if method == 'optimal': db_opt = db_spec
                return db_opt

            do_deblaze("MasterFlat")
            for c_file in calib_filenames_final:
                do_deblaze(Path(c_file).stem)
            
            deblazed_science_dict = {}
            for sci_name, opt_spec in extracted_science_dict.items():
                db_opt = do_deblaze(sci_name, opt_spec)
                if db_opt is not None:
                    deblazed_science_dict[sci_name] = db_opt

            # Step 7: Wavelength Calibration
            self._report_progress(0.90, "Wavelength Calibration")
            
            import os
            from astropy.time import Time
            time_key = self.config.get('data', 'statime_key', 'DATE-OBS')
            def get_file_time(filepath):
                if not filepath or not Path(filepath).exists(): return 0
                try:
                    _, hdr = read_fits_image(filepath)
                    if hdr and time_key in hdr:
                        return Time(hdr[time_key]).unix
                except Exception: pass
                return os.path.getmtime(filepath)

            calibrations = []
            for c_file in (calib_filenames_final if calib_filenames_final else calib_filenames):
                if not c_file: continue
                wave_calib = self.stage_wavelength_calibration(c_file)
                c_time = get_file_time(c_file)
                calibrations.append((c_time, wave_calib, Path(c_file).stem))
            
            if not calibrations:
                raise RuntimeError("No wavelength calibrations could be generated.")

            calibrated_science_dict = {}
            
            targets = [("MasterFlat", None)]
            for c_file in (calib_filenames_final if calib_filenames_final else calib_filenames):
                if c_file: targets.append((Path(c_file).stem, c_file))
            for sci_name in extracted_science_dict.keys():
                sci_file_match = next((f for f in science_scatter_files if Path(f).stem == sci_name), None)
                if not sci_file_match:
                    sci_file_match = next((f for f in raw_image_paths if Path(f).stem == sci_name), None)
                targets.append((sci_name, sci_file_match))

            for target_name, target_file in targets:
                target_time = get_file_time(target_file)
                closest_calib = min(calibrations, key=lambda x: abs(x[0] - target_time))
                wave_calib = closest_calib[1]
                logger.info(f"Using calibration {closest_calib[2]} for {target_name} (closest in time)")

                calib_opt = None
                for method in ['sum', 'optimal']:
                    from src.core.extraction import load_extracted_spectra
                    spectra = None
                    db_path = Path(self.config.get_output_path()) / 'step6_deblazing' / f'{target_name}_1D_{method}_Deblaze.fits'
                    if db_path.exists():
                        try: spectra = load_extracted_spectra(str(db_path))
                        except Exception: pass
                    if spectra is None:
                        ex_path = Path(self.config.get_output_path()) / 'step5_extraction' / f'{target_name}_1D_{method}.fits'
                        if ex_path.exists():
                            try: spectra = load_extracted_spectra(str(ex_path))
                            except Exception: pass
                    
                    if spectra is not None:
                        calib_spec = self.stage_wavelength_calibration_on_spectra(
                            spectra, wave_calib, science_name=target_name, method=method)
                        if method == 'optimal':
                            calib_opt = calib_spec
                
                if calib_opt is not None and target_name in extracted_science_dict:
                    calibrated_science_dict[target_name] = calib_opt

            self._report_progress(1.0, "Pipeline Complete")

            logger.info("=" * 60)
            logger.info("PIPELINE EXECUTION COMPLETED SUCCESSFULLY")
            logger.info("=" * 60)

            return calibrated_science_dict

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
