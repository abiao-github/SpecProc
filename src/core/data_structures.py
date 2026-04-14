"""
Core data structures for spectral data representation.

Defines classes for storing and managing spectral data throughout the reduction pipeline.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Dict, Optional, Tuple, Any
import json


@dataclass
class ApertureLocation:
    """Single aperture (order) location and properties."""

    aperture: int              # Aperture number
    order: int                 # Physical diffraction order

    # Polynomial coefficients for aperture position
    # position(y) = c0 + c1*y + c2*y^2 + ... (along dispersion axis)
    center_coef: np.ndarray    # Aperture center position coefficients
    lower_coef: np.ndarray     # Lower boundary coefficients
    upper_coef: np.ndarray     # Upper boundary coefficients

    width: float = 30.0        # Aperture width (pixels)

    def get_position(self, y: np.ndarray) -> np.ndarray:
        """Get aperture center position at spatial coordinate y.

        Coefficients are in descending order: c0*y^n + c1*y^(n-1) + ... + cn
        """
        return np.polyval(self.center_coef, y)

    def get_lower(self, y: np.ndarray) -> np.ndarray:
        """Get lower boundary at spatial coordinate y.

        Coefficients are in descending order: c0*y^n + c1*y^(n-1) + ... + cn
        """
        return np.polyval(self.lower_coef, y)

    def get_upper(self, y: np.ndarray) -> np.ndarray:
        """Get upper boundary at spatial coordinate y.

        Coefficients are in descending order: c0*y^n + c1*y^(n-1) + ... + cn
        """
        return np.polyval(self.upper_coef, y)


@dataclass
class ApertureSet:
    """Collection of apertures (orders) for entire detector."""

    apertures: Dict[int, ApertureLocation] = field(default_factory=dict)
    norders: int = 0
    direct_axis: int = 0  # 0 for Y-direction, 1 for X-direction

    def add_aperture(self, aperture: ApertureLocation):
        """Add aperture to set."""
        self.apertures[aperture.aperture] = aperture
        self.norders = len(self.apertures)

    def get_aperture(self, aperture_id: int) -> Optional[ApertureLocation]:
        """Retrieve aperture by ID."""
        return self.apertures.get(aperture_id)

    def get_orders(self) -> list:
        """Return list of order numbers."""
        return sorted(self.apertures.keys())

    def __len__(self) -> int:
        """Return number of apertures (orders) in the set."""
        return len(self.apertures)


@dataclass
class Spectrum:
    """1D spectrum data for a single order."""

    aperture: int
    order: int

    # Core spectrum data
    wavelength: np.ndarray    # Wavelength array [Angstrom]
    flux: np.ndarray          # Flux array [counts or normalized]
    error: Optional[np.ndarray] = None      # Error/uncertainty
    background: Optional[np.ndarray] = None # Background/stray light
    mask: Optional[np.ndarray] = None       # Pixel mask (0=good, 1=bad)

    # Metadata
    exptime: float = 1.0
    bjd: float = 0.0          # Barycentric Julian Date
    ra: float = 0.0           # Right Ascension
    dec: float = 0.0          # Declination

    def npixel(self) -> int:
        """Number of pixels in spectrum."""
        return len(self.flux)

    def get_good_pixels(self) -> np.ndarray:
        """Get mask of good pixels."""
        if self.mask is None:
            return np.ones(len(self.flux), dtype=bool)
        return self.mask == 0

    def to_dict(self) -> dict:
        """Convert to dictionary for storage."""
        return {
            'aperture': self.aperture,
            'order': self.order,
            'wavelength': self.wavelength.tolist(),
            'flux': self.flux.tolist(),
            'error': self.error.tolist() if self.error is not None else None,
            'background': self.background.tolist() if self.background is not None else None,
            'mask': self.mask.tolist() if self.mask is not None else None,
            'exptime': self.exptime,
            'bjd': self.bjd,
        }

    def copy(self):
        """Return a deep copy of this Spectrum instance.

        Uses numpy's copy for array fields to avoid shared references.
        """
        return type(self)(
            aperture=self.aperture,
            order=self.order,
            wavelength=np.copy(self.wavelength) if self.wavelength is not None else None,
            flux=np.copy(self.flux) if self.flux is not None else None,
            error=np.copy(self.error) if self.error is not None else None,
            background=np.copy(self.background) if self.background is not None else None,
            mask=np.copy(self.mask) if self.mask is not None else None,
            exptime=self.exptime,
            bjd=self.bjd,
            ra=self.ra,
            dec=self.dec,
        )


@dataclass
class SpectraSet:
    """Collection of 1D spectra for all orders."""

    spectra: Dict[int, Spectrum] = field(default_factory=dict)
    norders: int = 0

    def add_spectrum(self, spec: Spectrum):
        """Add spectrum to set."""
        self.spectra[spec.aperture] = spec
        self.norders = len(self.spectra)

    def get_spectrum(self, aperture_id: int) -> Optional[Spectrum]:
        """Retrieve spectrum by aperture ID."""
        return self.spectra.get(aperture_id)

    def get_orders(self) -> list:
        """Return list of aperture numbers."""
        return sorted(self.spectra.keys())


@dataclass
class FlatField:
    """Flat field calibration data."""

    # Raw combined flat field image
    flat_data: np.ndarray

    # Sensitivity map (relative response per pixel)
    flat_sens: np.ndarray

    # Normalized flat (for order tracing)
    flat_norm: np.ndarray

    # Bad/saturated pixel mask
    flat_mask: np.ndarray

    # Scattered light model subtracted from master flat
    scattered_light: Optional[np.ndarray] = None

    # Smoothed 2D model preserving blaze/order envelope
    smoothed_model: Optional[np.ndarray] = None

    # Pixel-to-pixel flat (high-frequency component near 1.0)
    pixel_flat: Optional[np.ndarray] = None

    # Illumination correction map (large-scale spatial profile)
    illumination_flat: Optional[np.ndarray] = None

    # Final 2D flat correction map used on science images
    flat_corr_2d: Optional[np.ndarray] = None

    # Per-order 1D flat spectra (optional)
    flat_spec: Optional[SpectraSet] = None

    # Aperture traces from flat field
    aperture_set: Optional[ApertureSet] = None

    # Cross-order spatial profiles per order (optional)
    cross_profiles: Dict[int, np.ndarray] = field(default_factory=dict)

    # Blaze profiles per order (for de-blazing)
    blaze_profiles: Dict[int, np.ndarray] = field(default_factory=dict)

    # Metadata
    exptime: float = 1.0


@dataclass
class WaveCalib:
    """Wavelength calibration solution."""

    # 2D polynomial coefficients for wavelength
    # wavelength(x, y) = sum(p_ij * x^i * y^j)
    poly_coef: np.ndarray     # Polynomial coefficients, shape (xorder+1, yorder+1)

    xorder: int = 3           # Wavelength poly order in X
    yorder: int = 3           # Wavelength poly order in Y

    # Identified calibration lines
    line_catalog: Optional[np.ndarray] = None  # Reference wavelengths
    line_pixels: Optional[np.ndarray] = None   # Detected pixel positions
    line_orders: Optional[np.ndarray] = None   # Order numbers of lines

    # Quality metrics
    rms: float = 0.0          # RMS residual [Angstrom]
    nlines: int = 0           # Number of identified lines

    # Metadata
    calib_type: str = "ThAr"  # Calibration lamp type

    def apply_to_pixel(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        """Apply wavelength calibration to pixel coordinates.

        Args:
            x: X pixel coordinates
            y: Y pixel coordinates (order number)

        Returns:
            Wavelength array
        """
        wavelength = np.zeros_like(x, dtype=float)
        for i in range(self.poly_coef.shape[0]):
            for j in range(self.poly_coef.shape[1]):
                wavelength += self.poly_coef[i, j] * (x ** i) * (y ** j)
        return wavelength


@dataclass
class ProcessingConfig:
    """Configuration for a single processing stage."""

    name: str
    params: Dict[str, Any] = field(default_factory=dict)
    enabled: bool = True

    def get_param(self, key: str, default: Any = None) -> Any:
        """Get parameter value."""
        return self.params.get(key, default)

    def set_param(self, key: str, value: Any):
        """Set parameter value."""
        self.params[key] = value


@dataclass
class ProcessingState:
    """State tracking for reduction pipeline."""

    # Input data
    raw_images: list = field(default_factory=list)  # List of raw FITS filenames

    # Calibration data
    bias_frame: Optional[np.ndarray] = None
    flat_field: Optional[FlatField] = None
    apertures: Optional[ApertureSet] = None
    wavelength_calib: Optional[WaveCalib] = None

    # Pipeline stage flags
    bias_done: bool = False
    flat_done: bool = False
    trace_done: bool = False
    wavelength_done: bool = False
    background_done: bool = False
    extraction_done: bool = False
    de_blazing_done: bool = False

    # Output data
    extracted_spectra: Optional[SpectraSet] = None
    de_blazed_spectra: Optional[SpectraSet] = None

    # Progress tracking
    current_stage: str = ""
    progress: float = 0.0  # 0.0 to 1.0

    def reset(self):
        """Reset pipeline state."""
        self.bias_done = False
        self.flat_done = False
        self.trace_done = False
        self.wavelength_done = False
        self.background_done = False
        self.extraction_done = False
        self.de_blazing_done = False
        self.current_stage = ""
        self.progress = 0.0
