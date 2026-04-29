"""
Spectrum extraction stage of spectral reduction pipeline.

Handles 1D spectrum extraction using sum or optimal (Gaussian) methods.
"""

import numpy as np
import logging
from pathlib import Path
from typing import Tuple, Optional
from scipy.optimize import curve_fit
from src.core.data_structures import Spectrum, SpectraSet, ApertureSet, WaveCalib, FlatField
from src.plotting.spectra_plotter import plot_spectrum_to_file, plot_spectra_to_pdf

logger = logging.getLogger(__name__)


class SpectrumExtractor:
    """Handles 1D spectrum extraction from 2D images."""

    def __init__(self, optimal_sigma: float = 3.0, method: str = 'optimal'):
        """
        Initialize spectrum extractor.

        Args:
            optimal_sigma: Gaussian sigma for optimal extraction profile.
            method: Default extraction method ('optimal' or 'sum').
        """
        self.optimal_sigma = optimal_sigma
        self.method = method
        self.extracted_spectra = None

    def sum_extraction(self, image: np.ndarray, aperture_positions: np.ndarray,
                      lower_bounds: np.ndarray = None, upper_bounds: np.ndarray = None) -> np.ndarray:
        """
        Extract spectrum using simple aperture sum.

        Args:
            image: 2D image array
            aperture_positions: Aperture center positions (per column)
            lower_bounds: Lower boundary positions (per column, absolute y)
            upper_bounds: Upper boundary positions (per column, absolute y)

        Returns:
            1D spectrum array
        """
        logger.debug("Performing sum extraction...")

        rows, cols = image.shape

        if len(aperture_positions) != cols:
            raise ValueError("Aperture positions length must match image width")

        spectrum = np.zeros(cols)

        for x in range(cols):
            center = aperture_positions[x]
            if lower_bounds is not None and upper_bounds is not None:
                y_min = max(0, int(np.floor(lower_bounds[x])))
                y_max = min(rows, int(np.ceil(upper_bounds[x])) + 1)
            else:
                y_min = max(0, int(center - 5))
                y_max = min(rows, int(center + 5) + 1)

            if y_max > y_min:
                spectrum[x] = np.sum(image[y_min:y_max, x])
            else:
                spectrum[x] = 0.0

        logger.debug(f"Sum extraction complete: {len(spectrum)} pixels")

        return spectrum

    def optimal_extraction(self, image: np.ndarray, aperture_positions: np.ndarray,
                          spatial_profile_2d: Optional[np.ndarray] = None,
                          lower_bounds: np.ndarray = None,
                          upper_bounds: np.ndarray = None) -> Tuple[np.ndarray, np.ndarray]:
        """
        Extract spectrum using optimal extraction.

        Uses a profile-weighted approach. If a flat-field-derived profile is
        provided, use it; otherwise fallback to Gaussian profile.

        Args:
            image: 2D image array
            aperture_positions: Aperture center positions (per column)
            profile: Optional spatial profile shape (relative weights)
            lower_bounds: Lower boundary positions (per column, absolute y)
            upper_bounds: Upper boundary positions (per column, absolute y)

        Returns:
            Tuple of (spectrum, error)
        """
        logger.debug("Performing optimal extraction...")

        rows, cols = image.shape

        if len(aperture_positions) != cols:
            raise ValueError("Aperture positions length must match image width")

        spectrum = np.zeros(cols)
        error = np.zeros(cols)

        for x in range(cols):
            center = aperture_positions[x]
            if lower_bounds is not None and upper_bounds is not None:
                y_min = max(0, int(np.floor(lower_bounds[x])))
                y_max = min(rows, int(np.ceil(upper_bounds[x])) + 1)
            else:
                y_min = max(0, int(center - 5))
                y_max = min(rows, int(center + 5) + 1)

            if y_max <= y_min:
                spectrum[x] = 0.0
                error[x] = 0.0
                continue

            sub_image = image[y_min:y_max, x]
            local_y = np.arange(y_min, y_max) - center

            weights = None
            if spatial_profile_2d is not None:
                p_col = spatial_profile_2d[y_min:y_max, x].astype(np.float64)
                sum_p = np.sum(p_col)
                if sum_p > 0:
                    weights = p_col / sum_p

            if weights is None:
                sigma = self.optimal_sigma
                weights = np.exp(-0.5 * (local_y / sigma) ** 2)
                weights /= np.sum(weights) if np.sum(weights) > 0 else 1.0

            weighted = sub_image * weights
            spectrum[x] = np.sum(weighted)
            error[x] = np.sqrt(np.sum((sub_image * weights) ** 2))

        logger.debug(f"Optimal extraction complete: {len(spectrum)} pixels")

        return spectrum, error

    def extract_aperture_set(self, image: np.ndarray, apertures: ApertureSet,
                           wavelength_calib: Optional[WaveCalib] = None,
                           flat_field: Optional[FlatField] = None,
                           method_override: Optional[str] = None) -> SpectraSet:
        """
        Extract all apertures from 2D image.

        Args:
            image: 2D science image
            apertures: ApertureSet with aperture definitions
            wavelength_calib: Optional wavelength calibration

        Returns:
            SpectraSet with extracted spectra
        """
        logger.info(f"Extracting {apertures.norders} apertures...")

        spectra_set = SpectraSet()

        method = method_override if method_override is not None else self.method

        _, cols = image.shape

        for aperture_id in apertures.get_orders():
            aperture = apertures.get_aperture(aperture_id)
            if aperture is None:
                continue

            x_pix = np.arange(cols)
            center_pos = aperture.get_position(x_pix)
            lower_bounds = aperture.get_lower(x_pix)
            upper_bounds = aperture.get_upper(x_pix)

            spatial_profile_2d = None
            if flat_field is not None:
                if flat_field.smoothed_model is not None:
                    spatial_profile_2d = flat_field.smoothed_model

            if method == 'optimal':
                flux, error = self.optimal_extraction(image, center_pos,
                                                      spatial_profile_2d=spatial_profile_2d,
                                                      lower_bounds=lower_bounds,
                                                      upper_bounds=upper_bounds)
            else:
                flux = self.sum_extraction(image, center_pos,
                                          lower_bounds=lower_bounds,
                                          upper_bounds=upper_bounds)
                error = np.sqrt(np.abs(flux))

            # Apply wavelength calibration if available
            if wavelength_calib is not None:
                wavelength = wavelength_calib.apply_to_pixel(x_pix,
                                                            aperture_id * np.ones_like(x_pix))
            else:
                wavelength = x_pix

            # Create Spectrum object
            spec = Spectrum(
                aperture=aperture_id,
                order=aperture.order,
                wavelength=wavelength,
                flux=flux,
                error=error,
                background=None,
                mask=None
            )

            spectra_set.add_spectrum(spec)

        self.extracted_spectra = spectra_set
        logger.info(f"Extraction complete: {spectra_set.norders} spectra")

        return spectra_set

    def save_spectra(self, output_path: str, spectra_set: SpectraSet):
        """Save extracted spectra to FITS file."""
        if spectra_set is None:
            raise RuntimeError("No spectra to save")

        from astropy.io import fits
        from pathlib import Path

        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Create structured array for table
        norders = spectra_set.norders
        max_pixels = max([s.npixel() for s in spectra_set.spectra.values()])


        # Create binary table
        cols = []
        for aperture_id in spectra_set.get_orders():
            spec = spectra_set.get_spectrum(aperture_id)

            # Pad arrays to same length
            wl = np.pad(spec.wavelength, (0, max_pixels - len(spec.wavelength)))
            fl = np.pad(spec.flux, (0, max_pixels - len(spec.flux)))

            cols.append(fits.Column(name=f'WAV_{aperture_id:02d}', format='D', array=wl))
            cols.append(fits.Column(name=f'FLUX_{aperture_id:02d}', format='D', array=fl))

            if spec.error is not None:
                err = np.pad(spec.error, (0, max_pixels - len(spec.error)))
                cols.append(fits.Column(name=f'ERR_{aperture_id:02d}', format='D', array=err))

        # Create HDUs
        primary_hdu = fits.PrimaryHDU()
        table_hdu = fits.BinTableHDU.from_columns(cols)
        table_hdu.header['NORDERS'] = norders

        hdul = fits.HDUList([primary_hdu, table_hdu])
        hdul.writeto(str(output_path), overwrite=True)

        logger.info(f"Saved {norders} spectra to {output_path}")


def process_extraction_stage(science_image: np.ndarray,
                            apertures: ApertureSet,
                            output_dir_base: str,
                            wavelength_calib: Optional[WaveCalib] = None,
                            flat_field: Optional[FlatField] = None,
                            method_override: Optional[str] = None,
                            output_filename: str = 'extracted_spectra.fits',
                            plot_prefix: str = 'extracted',
                            optimal_sigma: float = 3.0,
                            extraction_method: str = 'optimal',
                            save_plots: bool = True,
                            fig_format: str = 'png') -> SpectraSet:
    """
    Execute spectrum extraction stage.

    Args:
        science_image: 2D science image
        apertures: Detected apertures
        output_dir_base: Base output directory
        wavelength_calib: Wavelength calibration
        flat_field: Flat field data
        method_override: Override default extraction method
        output_filename: Name for output FITS file
        plot_prefix: Prefix for diagnostic plots
        optimal_sigma: Sigma for Gaussian profile in optimal extraction
        extraction_method: Default extraction method
        save_plots: Whether to save diagnostic plots
        fig_format: Format for diagnostic plots

    Returns:
        SpectraSet with extracted 1D spectra
    """
    extractor = SpectrumExtractor(optimal_sigma=optimal_sigma, method=extraction_method)

    # Extract spectra
    spectra_set = extractor.extract_aperture_set(science_image, apertures,
                                                wavelength_calib, flat_field,
                                                method_override=method_override)

    # Save spectra
    spectra_file = Path(output_dir_base) / 'step5_extraction' / output_filename
    spectra_file.parent.mkdir(parents=True, exist_ok=True)
    extractor.save_spectra(str(spectra_file), spectra_set)

    # Save diagnostic plots if enabled
    if save_plots:
        out_dir = spectra_file.parent
        if wavelength_calib is not None:
            pdf_path = out_dir / f"{plot_prefix}_all_orders.pdf"
            plot_spectra_to_pdf(spectra_set, str(pdf_path), title_prefix="Extracted Spectrum")
        else:
            # If no wavelength calibration, plot raw pixels
            pdf_name = f"{plot_prefix}.pdf"
            pdf_path = out_dir / pdf_name
            plot_spectra_to_pdf(spectra_set, str(pdf_path), title_prefix="Extracted Spectrum (pixels)", xlabel="Pixel")

    return spectra_set

def load_extracted_spectra(filepath: str) -> SpectraSet:
    """Load extracted spectra from FITS file."""
    from astropy.io import fits
    from src.core.data_structures import Spectrum, SpectraSet
    
    spectra_set = SpectraSet()
    
    with fits.open(filepath) as hdul:
        if len(hdul) < 2:
            raise ValueError("FITS file does not contain a binary table")
            
        table_hdu = hdul[1]
        
        # Extract unique aperture IDs from column names
        aperture_ids = set()
        for col in table_hdu.columns:
            if col.name.startswith('FLUX_'):
                try:
                    aperture_ids.add(int(col.name.split('_')[1]))
                except ValueError:
                    pass
                
        for ap_id in sorted(list(aperture_ids)):
            wav_col = f'WAV_{ap_id:02d}'
            flux_col = f'FLUX_{ap_id:02d}'
            err_col = f'ERR_{ap_id:02d}'
            
            if wav_col in table_hdu.data.columns.names and flux_col in table_hdu.data.columns.names:
                wav = table_hdu.data[wav_col]
                flux = table_hdu.data[flux_col]
                err = table_hdu.data[err_col] if err_col in table_hdu.data.columns.names else None
                
                if len(wav.shape) > 1 and wav.shape[0] == 1:
                    wav = wav[0]
                if len(flux.shape) > 1 and flux.shape[0] == 1:
                    flux = flux[0]
                if err is not None and len(err.shape) > 1 and err.shape[0] == 1:
                    err = err[0]
                
                spec = Spectrum(
                    aperture=ap_id,
                    order=ap_id, # order is same as aperture for now
                    wavelength=wav,
                    flux=flux,
                    error=err,
                    background=None,
                    mask=None
                )
                spectra_set.add_spectrum(spec)
                
    return spectra_set
