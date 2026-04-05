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
from src.config.config_manager import ConfigManager
from src.plotting.spectra_plotter import plot_spectrum_to_file

logger = logging.getLogger(__name__)


class SpectrumExtractor:
    """Handles 1D spectrum extraction from 2D images."""

    def __init__(self, config: ConfigManager):
        """
        Initialize spectrum extractor.

        Args:
            config: Configuration manager
        """
        self.config = config
        self.extracted_spectra = None

    def sum_extraction(self, image: np.ndarray, aperture_positions: np.ndarray,
                      lower_limit: float = -5.0, upper_limit: float = 5.0) -> np.ndarray:
        """
        Extract spectrum using simple aperture sum.

        Args:
            image: 2D image array
            aperture_positions: Aperture center positions (per column)
            lower_limit: Lower extraction limit relative to center (pixels)
            upper_limit: Upper extraction limit relative to center (pixels)

        Returns:
            1D spectrum array
        """
        logger.info("Performing sum extraction...")

        rows, cols = image.shape

        if len(aperture_positions) != cols:
            raise ValueError("Aperture positions length must match image width")

        spectrum = np.zeros(cols)

        for x in range(cols):
            center = aperture_positions[x]
            y_min = max(0, int(center + lower_limit))
            y_max = min(rows, int(center + upper_limit) + 1)

            if y_max > y_min:
                spectrum[x] = np.sum(image[y_min:y_max, x])
            else:
                spectrum[x] = 0.0

        logger.info(f"Sum extraction complete: {len(spectrum)} pixels")

        return spectrum

    def optimal_extraction(self, image: np.ndarray, aperture_positions: np.ndarray,
                          profile: Optional[np.ndarray] = None,
                          lower_limit: float = -5.0, upper_limit: float = 5.0) -> Tuple[np.ndarray, np.ndarray]:
        """
        Extract spectrum using optimal extraction.

        Uses a profile-weighted approach. If a flat-field-derived profile is
        provided, use it; otherwise fallback to Gaussian profile.

        Args:
            image: 2D image array
            aperture_positions: Aperture center positions (per column)
            profile: Optional spatial profile shape (relative weights)
            lower_limit: lower extraction limit relative to center (pixels)
            upper_limit: upper extraction limit relative to center (pixels)

        Returns:
            Tuple of (spectrum, error)
        """
        logger.info("Performing optimal extraction...")

        rows, cols = image.shape

        if len(aperture_positions) != cols:
            raise ValueError("Aperture positions length must match image width")

        spectrum = np.zeros(cols)
        error = np.zeros(cols)

        for x in range(cols):
            center = aperture_positions[x]
            y_min = max(0, int(center + lower_limit))
            y_max = min(rows, int(center + upper_limit) + 1)

            if y_max <= y_min:
                spectrum[x] = 0.0
                error[x] = 0.0
                continue

            sub_image = image[y_min:y_max, x]
            local_y = np.arange(y_min, y_max) - center

            if profile is not None and len(profile) >= (y_max - y_min):
                p0 = len(profile) // 2
                pstart = max(0, p0 - (y_max - y_min) // 2)
                pend = pstart + (y_max - y_min)
                weights = profile[pstart:pend].astype(np.float64)
                weights /= np.sum(weights) if np.sum(weights) > 0 else 1.0
            else:
                sigma = self.config.get_float('reduce.extract', 'optimal_sigma', 3.0)
                weights = np.exp(-0.5 * (local_y / sigma) ** 2)
                weights /= np.sum(weights) if np.sum(weights) > 0 else 1.0

            weighted = sub_image * weights
            spectrum[x] = np.sum(weighted)
            error[x] = np.sqrt(np.sum((sub_image * weights) ** 2))

        logger.info(f"Optimal extraction complete: {len(spectrum)} pixels")

        return spectrum, error

    def extract_aperture_set(self, image: np.ndarray, apertures: ApertureSet,
                           wavelength_calib: Optional[WaveCalib] = None) -> SpectraSet:
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

        method = self.config.get('reduce.extract', 'method', 'sum')
        lower_limit = self.config.get_float('reduce.extract', 'lower_limit', -5.0)
        upper_limit = self.config.get_float('reduce.extract', 'upper_limit', 5.0)

        _, cols = image.shape

        for aperture_id in apertures.get_orders():
            aperture = apertures.get_aperture(aperture_id)
            if aperture is None:
                continue

            x_pix = np.arange(cols)
            center_pos = aperture.get_position(x_pix)

            profile = None
            if flat_field is not None and flat_field.flat_norm is not None:
                y_mid = int(np.median(center_pos))
                profile_window = flat_field.flat_norm[max(y_mid-30,0):min(y_mid+30, flat_field.flat_norm.shape[0]), :]
                profile = np.mean(profile_window, axis=1)
                profile = profile / (np.sum(profile) if np.sum(profile) > 0 else 1.0)

            if method == 'optimal':
                flux, error = self.optimal_extraction(image, center_pos,
                                                      profile=profile,
                                                      lower_limit=lower_limit,
                                                      upper_limit=upper_limit)
            else:
                flux = self.sum_extraction(image, center_pos,
                                          lower_limit=lower_limit,
                                          upper_limit=upper_limit)
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

            cols.append(fits.Column(name=f'WAV_{aperture_id:02d}', format='D', array=[wl]))
            cols.append(fits.Column(name=f'FLUX_{aperture_id:02d}', format='D', array=[fl]))

            if spec.error is not None:
                err = np.pad(spec.error, (0, max_pixels - len(spec.error)))
                cols.append(fits.Column(name=f'ERR_{aperture_id:02d}', format='D', array=[err]))

        # Create HDUs
        primary_hdu = fits.PrimaryHDU()
        table_hdu = fits.BinTableHDU.from_columns(cols)
        table_hdu.header['NORDERS'] = norders

        hdul = fits.HDUList([primary_hdu, table_hdu])
        hdul.writeto(str(output_path), overwrite=True)

        logger.info(f"Saved {norders} spectra to {output_path}")


def process_extraction_stage(config: ConfigManager, science_image: np.ndarray,
                            apertures: ApertureSet,
                            wavelength_calib: Optional[WaveCalib] = None,
                            flat_field: Optional[FlatField] = None,
                            midpath: str = None) -> SpectraSet:
    """
    Execute spectrum extraction stage.

    Args:
        config: Configuration manager
        science_image: 2D science image
        apertures: Detected apertures
        wavelength_calib: Wavelength calibration
        midpath: Intermediate output directory

    Returns:
        SpectraSet with extracted 1D spectra
    """
    # Use default midpath if not provided
    if midpath is None:
        out_path = config.get('reduce', 'out_path', './output')
        midpath = str(Path(out_path) / 'spectra')

    extractor = SpectrumExtractor(config)

    # Extract spectra
    spectra_set = extractor.extract_aperture_set(science_image, apertures,
                                                wavelength_calib)

    # Save spectra
    out_path = config.get('reduce', 'out_path', './output')
    spectra_file = Path(out_path) / 'step5_extraction' / 'extracted_spectra.fits'
    spectra_file.parent.mkdir(parents=True, exist_ok=True)
    extractor.save_spectra(str(spectra_file), spectra_set)

    # Save diagnostic plots if enabled
    save_plots = config.get_bool('reduce', 'save_plots', True)
    if save_plots:
        out_dir = spectra_file.parent
        fig_format = config.get('reduce', 'fig_format', 'png')
        for spectrum in spectra_set.spectra.values():
            order = spectrum.aperture
            if wavelength_calib is not None and len(spectrum.wavelength) > 0:
                plot_file = out_dir / f'extracted_order_{order:02d}.{fig_format}'
                plot_spectrum_to_file(
                    spectrum.wavelength,
                    spectrum.flux,
                    str(plot_file),
                    spectrum.error if spectrum.error is not None else None,
                    f"Extracted Spectrum - Order {order}"
                )

    return spectra_set
