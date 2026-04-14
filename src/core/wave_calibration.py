"""
Wavelength calibration stage of spectral reduction pipeline.

Handles ThAr line identification, polynomial fitting,
and wavelength solution computation.
"""

import numpy as np
import logging
from typing import Tuple, Optional, List, Dict
from scipy.optimize import curve_fit
from pathlib import Path
from src.core.data_structures import WaveCalib, Spectrum, SpectraSet
from src.config.config_manager import ConfigManager
from src.utils.fits_io import read_fits_image
from src.plotting.spectra_plotter import plot_wavelength_calibration

logger = logging.getLogger(__name__)


class WavelengthCalibrator:
    """Handles wavelength calibration from reference lines."""

    def __init__(self, config: ConfigManager):
        """
        Initialize wavelength calibrator.

        Args:
            config: Configuration manager
        """
        self.config = config
        self.wave_calib = None
        self.reference_lines = {}

    def load_line_list(self, lamp_type: str = 'ThAr') -> Dict[str, np.ndarray]:
        """
        Load reference line list for calibration lamps.

        Args:
            lamp_type: Lamp type ('ThAr', 'Ar', 'Ne', 'He', 'Fe')

        Returns:
            Dictionary with wavelength line lists
        """
        logger.info(f"Loading {lamp_type} line list...")

        # For now, create synthetic line lists
        # In a real implementation, these would be loaded from data files
        linelists = {
            'ThAr': {
                'wavelength': np.array([
                    3888.6, 3944.6, 4078.6, 4131.7, 4181.9,
                    5294.5, 5330.8, 5688.2, 5944.8, 5975.5,
                    6052.72, 6083.3, 6170.17, 6212.5, 6318.16,
                ]),
                'strength': np.ones(15)
            },
            'Ne': {
                'wavelength': np.array([
                    3965.0, 4144.4, 4169.0, 4198.8, 4226.7,
                    5852.5, 5881.9, 5902.8, 5929.7, 5944.8,
                ]),
                'strength': np.ones(10)
            },
            'Ar': {
                'wavelength': np.array([
                    4052.9, 4131.7, 4158.6, 4181.9, 4254.4,
                    5006.1, 5017.2, 5047.7, 5062.5, 5141.8,
                ]),
                'strength': np.ones(10)
            },
        }

        self.reference_lines[lamp_type] = linelists.get(lamp_type, {})
        return linelists.get(lamp_type, {})

    def detect_lines_in_spectrum(self, spectrum: np.ndarray,
                                wavelength: Optional[np.ndarray] = None,
                                threshold: float = 3.0) -> Tuple[np.ndarray, np.ndarray]:
        """
        Detect emission lines in spectrum.

        Args:
            spectrum: 1D flux array
            wavelength: Optional wavelength array (for initial guess)
            threshold: DetectionThreshold in sigma

        Returns:
            Tuple of (pixel_positions, line_strengths)
        """
        # Find peaks in spectrum
        from scipy.signal import find_peaks

        mean = np.mean(spectrum)
        std = np.std(spectrum)
        peaks, properties = find_peaks(spectrum, height=threshold * std, distance=5)

        logger.info(f"Detected {len(peaks)} emission lines")

        return peaks, spectrum[peaks]

    def fit_wavelength_polynomial(self, pixel_positions: np.ndarray,
                                 wavelengths: np.ndarray,
                                 xorder: int = 4, yorder: int = 4) -> WaveCalib:
        """
        Fit 2D wavelength polynomial.

        Args:
            pixel_positions: Array of (x, y) pixel coordinates
            wavelengths: Corresponding wavelengths
            xorder: Polynomial order in X
            yorder: Polynomial order in Y

        Returns:
            WaveCalib object with polynomial solution
        """
        logger.info(f"Fitting wavelength polynomial (orders: {xorder}, {yorder})...")

        if len(pixel_positions) != len(wavelengths):
            raise ValueError("Pixel positions and wavelengths must have same length")

        # Extract x, y coordinates
        if pixel_positions.ndim == 2:
            x_pix = pixel_positions[:, 0]
            y_pix = pixel_positions[:, 1]
        else:
            x_pix = pixel_positions
            y_pix = np.zeros_like(x_pix)

        # Build polynomial design matrix
        ncoeff = (xorder + 1) * (yorder + 1)
        A = np.zeros((len(x_pix), ncoeff))

        idx = 0
        for i in range(xorder + 1):
            for j in range(yorder + 1):
                A[:, idx] = (x_pix ** i) * (y_pix ** j)
                idx += 1

        # Solve least squares
        try:
            coeffs, residuals, rank, s = np.linalg.lstsq(A, wavelengths, rcond=None)

            # Compute RMS of residuals
            predicted = A @ coeffs
            rms = np.sqrt(np.mean((predicted - wavelengths) ** 2))

            logger.info(f"Wavelength fit RMS: {rms:.4f} Angstrom")

            # Reshape coefficients into 2D array
            poly_coef = coeffs.reshape((xorder + 1, yorder + 1))

            # Create WaveCalib object
            wave_calib = WaveCalib(
                poly_coef=poly_coef,
                xorder=xorder,
                yorder=yorder,
                line_pixels=x_pix,
                line_orders=y_pix,
                line_catalog=wavelengths,
                rms=rms,
                nlines=len(wavelengths),
                calib_type='ThAr'
            )

            self.wave_calib = wave_calib
            return wave_calib

        except Exception as e:
            logger.error(f"Error fitting wavelength polynomial: {e}")
            raise

    def apply_wavelength_calibration(self, spectrum: np.ndarray,
                                    pixel_array: np.ndarray,
                                    aperture_y: float = 0.0) -> np.ndarray:
        """
        Apply wavelength calibration to spectrum pixels.

        Args:
            spectrum: 1D spectrum array
            pixel_array: Pixel coordinate array
            aperture_y: Y coordinate (order number) of aperture

        Returns:
            Wavelength array
        """
        if self.wave_calib is None:
            raise RuntimeError("No wavelength calibration available")

        # Apply polynomial
        wavelength = self.wave_calib.apply_to_pixel(pixel_array, aperture_y * np.ones_like(pixel_array))

        return wavelength

    def save_calibration(self, output_path: str):
        """Save wavelength calibration to FITS file."""
        if self.wave_calib is None:
            raise RuntimeError("No wavelength calibration to save")

        from astropy.io import fits

        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Build a proper fits.Header (PrimaryHDU rejects plain dicts)
        hdr = fits.Header()
        hdr['XORDER']    = self.wave_calib.xorder
        hdr['YORDER']    = self.wave_calib.yorder
        hdr['RMS']       = self.wave_calib.rms
        hdr['NLINES']    = self.wave_calib.nlines
        hdr['CALIB_TYP'] = self.wave_calib.calib_type

        hdul = fits.HDUList([
            fits.PrimaryHDU(data=self.wave_calib.poly_coef, header=hdr),
        ])

        hdul.writeto(str(output_path), overwrite=True)
        logger.info(f"Saved wavelength calibration to {output_path}")


def process_wavelength_stage(config: ConfigManager, calib_filename: str,
                            midpath: str = './midpath') -> WaveCalib:
    """
    Execute wavelength calibration stage.

    Args:
        config: Configuration manager
        calib_filename: Path to calibration (ThAr) FITS file
        midpath: Intermediate output directory

    Returns:
        WaveCalib object
    """
    calibrator = WavelengthCalibrator(config)

    # Load line list
    lamp_type = config.get('reduce.wlcalib', 'linelist', 'ThAr')
    linelist = calibrator.load_line_list(lamp_type)

    # For demonstration, create synthetic wavelength calibration
    # In real use, would extract spectrum and identify lines
    x_order = config.get_int('reduce.wlcalib', 'xorder', 4)
    y_order = config.get_int('reduce.wlcalib', 'yorder', 4)

    # Create demo calibration (linear)
    coeffs = np.zeros(((x_order + 1) * (y_order + 1),))
    coeffs[0] = 3000.0  # Wavelength offset
    coeffs[1] = 0.5     # Wavelength scale

    poly_coef = coeffs.reshape((x_order + 1, y_order + 1))

    wave_calib = WaveCalib(
        poly_coef=poly_coef,
        xorder=x_order,
        yorder=y_order,
        rms=0.1,
        nlines=15,
        calib_type=lamp_type
    )

    # Save calibration
    base_output_path = config.get_output_path()
    calib_file = Path(base_output_path) / 'step5_wavelength' / 'wavelength_calibration.fits'
    calib_file.parent.mkdir(parents=True, exist_ok=True)
    calibrator.wave_calib = wave_calib
    calibrator.save_calibration(str(calib_file))

    # Save diagnostic plot if enabled
    save_plots = config.get_bool('reduce', 'save_plots', True)
    if save_plots and calibrator.line_pixels is not None:
        out_dir = calib_file.parent
        fig_format = config.get('reduce', 'fig_format', 'png')
        plot_file = out_dir / f'wavelength_calibration.{fig_format}'
        plot_wavelength_calibration(
            calibrator.line_pixels,
            calibrator.line_wavelengths,
            calibrator.fitted_wavelengths,
            str(plot_file)
        )

    return wave_calib
