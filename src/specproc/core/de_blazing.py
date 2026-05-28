"""
De-blazing correction stage.

Corrects for the blaze function of the spectrograph grating
to flatten the spectral response.
"""

import numpy as np
import logging
from pathlib import Path
from typing import Optional
from specproc.core.data_structures import SpectraSet, FlatField
from specproc.plotting.spectra_plotter import plot_spectrum_to_file, plot_spectra_to_pdf

logger = logging.getLogger(__name__)


def apply_de_blazing(spectra_set: SpectraSet, flat_field: FlatField) -> SpectraSet:
    """
    Apply de-blazing correction to extracted spectra.

    Divides each spectrum by its corresponding blaze profile
    to remove wavelength-dependent efficiency variations.

    Args:
        spectra_set: Extracted spectra to correct
        flat_field: Flat field data containing blaze profiles

    Returns:
        De-blazed spectra set
    """
    if flat_field.blaze_profiles is None or len(flat_field.blaze_profiles) == 0:
        logger.warning("No blaze profiles available, skipping de-blazing")
        return spectra_set

    corrected_spectra = SpectraSet()

    for spectrum in spectra_set.spectra.values():
        aperture_id = spectrum.aperture

        if aperture_id in flat_field.blaze_profiles:
            blaze_profile = flat_field.blaze_profiles[aperture_id]

            # Ensure blaze profile matches spectrum length
            if len(blaze_profile) == len(spectrum.flux):
                blaze_peak = np.nanmax(blaze_profile)
                if not np.isfinite(blaze_peak) or blaze_peak <= 0:
                    blaze_peak = 1.0
                
                blaze_threshold = 0.02 * blaze_peak

                # --- FIX: Contiguous valid region ---
                # Only keep the contiguous region around the peak that is above the threshold.
                # This strictly ignores any artificial "up-swings" or polynomial ringing at the extreme edges.
                bad_blaze_mask = np.ones_like(blaze_profile, dtype=bool)
                peak_idx = int(np.nanargmax(blaze_profile))
                
                # Walk left from the peak
                left_bound = peak_idx
                while left_bound > 0 and blaze_profile[left_bound - 1] >= blaze_threshold:
                    left_bound -= 1
                    
                # Walk right from the peak
                right_bound = peak_idx
                while right_bound < len(blaze_profile) - 1 and blaze_profile[right_bound + 1] >= blaze_threshold:
                    right_bound += 1
                    
                # The contiguous region from left_bound to right_bound is good.
                bad_blaze_mask[left_bound:right_bound + 1] = False

                # Avoid division by zero or very small numbers by setting them to NaN
                blaze_for_division = np.where(bad_blaze_mask, np.nan, blaze_profile)

                # Apply de-blazing, ignoring division warnings for NaNs
                with np.errstate(divide='ignore', invalid='ignore'):
                    corrected_flux = spectrum.flux / blaze_for_division
                    corrected_error = spectrum.error / blaze_for_division if spectrum.error is not None else None

                # In regions of low blaze, the result is unreliable. Set flux and error to NaN.
                corrected_flux[bad_blaze_mask] = np.nan
                if corrected_error is not None:
                    corrected_error[bad_blaze_mask] = np.nan
                
                # Clean up any other non-finite values that might have resulted from the division
                corrected_flux[~np.isfinite(corrected_flux)] = np.nan
                if corrected_error is not None:
                    corrected_error[~np.isfinite(corrected_error)] = np.nan

                # Create corrected spectrum
                corrected_spectrum = spectrum.copy()
                corrected_spectrum.flux = corrected_flux
                corrected_spectrum.error = corrected_error
                corrected_spectra.add_spectrum(corrected_spectrum)
            else:
                logger.warning(f"Blaze profile length mismatch for aperture {aperture_id}, skipping")
                corrected_spectra.add_spectrum(spectrum)
        else:
            logger.warning(f"No blaze profile for aperture {aperture_id}, skipping de-blazing")
            corrected_spectra.add_spectrum(spectrum)

    logger.info(f"✓ De-blazing complete: {len(corrected_spectra.spectra)} spectra processed")
    return corrected_spectra


def save_deblazed_spectra(output_path: str, spectra_set: SpectraSet, 
                          content_desc: str = 'De-blazed calibrated spectra',
                          log_desc: str = 'de-blazed'):
    """Save de-blazed spectra to FITS file."""
    if spectra_set is None:
        raise RuntimeError("No de-blazed spectra to save")

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
        wl = np.pad(spec.wavelength, (0, max_pixels - len(spec.wavelength)), constant_values=np.nan)
        fl = np.pad(spec.flux, (0, max_pixels - len(spec.flux)), constant_values=np.nan)

        cols.append(fits.Column(name=f'WAV_{aperture_id:02d}', format='D', array=wl))
        cols.append(fits.Column(name=f'FLUX_{aperture_id:02d}', format='D', array=fl))

        if spec.error is not None:
            err = np.pad(spec.error, (0, max_pixels - len(spec.error)), constant_values=np.nan)
            cols.append(fits.Column(name=f'ERR_{aperture_id:02d}', format='D', array=err))

    coldefs = fits.ColDefs(cols)
    table_hdu = fits.BinTableHDU.from_columns(coldefs)

    # Create header
    header = fits.Header()
    header['NORDERS'] = norders
    header['NPX'] = max_pixels
    header['CONTENTS'] = content_desc
    header['COMMENT'] = 'Wavelength in Angstroms, Flux in ADU'

    primary_hdu = fits.PrimaryHDU(header=header)

    hdul = fits.HDUList([primary_hdu, table_hdu])
    hdul.writeto(str(output_path), overwrite=True)

    logger.info(f"Saved {norders} {log_desc} spectra to {output_path}")


def process_de_blazing_stage(spectra_set: SpectraSet,
                           output_dir_base: str,
                           flat_field: Optional[FlatField] = None,
                           save_output: bool = True,
                           save_deblaze: bool = True,
                           output_filename: str = 'deblazed_spectra.fits',
                           plot_prefix: str = 'deblazed',
                           save_plots: bool = True,
                           fig_format: str = 'png') -> SpectraSet:
    """
    Execute de-blazing correction stage.

    Args:
        spectra_set: Extracted spectra
        output_dir_base: Base output directory
        flat_field: Flat field data with blaze profiles
        save_output: Whether to save de-blazed spectra to file
        save_deblaze: Whether to save de-blazed intermediate products
        output_filename: Filename for the saved FITS file
        plot_prefix: Prefix for saved plot files
        save_plots: Whether to save diagnostic plots
        fig_format: Format for diagnostic plots

    Returns:
        De-blazed spectra set
    """
    if flat_field is None:
        logger.warning("No flat field provided, skipping de-blazing")
        return spectra_set

    deblazed_spectra = apply_de_blazing(spectra_set, flat_field)

    # Save de-blazed spectra if enabled
    if save_output and save_deblaze:
        deblazed_file = Path(output_dir_base) / 'step6_deblazing' / output_filename
        deblazed_file.parent.mkdir(parents=True, exist_ok=True)
        save_deblazed_spectra(str(deblazed_file), deblazed_spectra)

        # Save diagnostic plots if enabled
        if save_plots:
            out_dir = deblazed_file.parent
            pdf_name = f"{plot_prefix}.pdf"
            pdf_path = out_dir / pdf_name
            plot_spectra_to_pdf(deblazed_spectra, str(pdf_path), title_prefix="De-blazed Spectrum (pixels)", xlabel="Pixel")

    return deblazed_spectra