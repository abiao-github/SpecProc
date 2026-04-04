"""
De-blazing correction stage.

Corrects for the blaze function of the spectrograph grating
to flatten the spectral response.
"""

import numpy as np
import logging
from pathlib import Path
from typing import Optional
from src.core.data_structures import SpectraSet, FlatField
from src.config.config_manager import ConfigManager

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
    logger.info("=" * 60)
    logger.info("STAGE 6: DE-BLAZING CORRECTION")
    logger.info("=" * 60)

    if flat_field.blaze_profiles is None or len(flat_field.blaze_profiles) == 0:
        logger.warning("No blaze profiles available, skipping de-blazing")
        return spectra_set

    corrected_spectra = SpectraSet()

    for spectrum in spectra_set.spectra:
        aperture_id = spectrum.aperture

        if aperture_id in flat_field.blaze_profiles:
            blaze_profile = flat_field.blaze_profiles[aperture_id]

            # Ensure blaze profile matches spectrum length
            if len(blaze_profile) == len(spectrum.flux):
                # Apply de-blazing: divide by blaze function
                corrected_flux = spectrum.flux / (blaze_profile + 1e-10)  # Avoid division by zero

                # Create corrected spectrum
                corrected_spectrum = spectrum.copy()
                corrected_spectrum.flux = corrected_flux

                corrected_spectra.add_spectrum(corrected_spectrum)

                logger.info(f"De-blazed spectrum for aperture {aperture_id}")
            else:
                logger.warning(f"Blaze profile length mismatch for aperture {aperture_id}, skipping")
                corrected_spectra.add_spectrum(spectrum)
        else:
            logger.warning(f"No blaze profile for aperture {aperture_id}, skipping de-blazing")
            corrected_spectra.add_spectrum(spectrum)

    logger.info(f"✓ De-blazing complete: {len(corrected_spectra.spectra)} spectra processed")
    return corrected_spectra


def process_de_blazing_stage(config: ConfigManager, spectra_set: SpectraSet,
                           flat_field: Optional[FlatField] = None) -> SpectraSet:
    """
    Execute de-blazing correction stage.

    Args:
        config: Configuration manager
        spectra_set: Extracted spectra
        flat_field: Flat field data with blaze profiles

    Returns:
        De-blazed spectra set
    """
    if flat_field is None:
        logger.warning("No flat field provided, skipping de-blazing")
        return spectra_set

    return apply_de_blazing(spectra_set, flat_field)