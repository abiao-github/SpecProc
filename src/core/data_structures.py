"""
Core data structures for spectral data representation.

Defines classes for storing and managing spectral data throughout the reduction pipeline.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Dict, Optional, Tuple, Any
import json
from numpy.polynomial.chebyshev import chebval


@dataclass
class ApertureLocation:
    """Single aperture (order) location and properties.

    Supports two storage modes:
    - **Polynomial**: ``center_coef`` / ``lower_coef`` / ``upper_coef`` are
      polynomial coefficient arrays evaluated via ``np.polyval``.
    - **Dense array**: ``center_arr`` / ``lower_arr`` / ``upper_arr`` are
      per-column value arrays (one value per detector column) evaluated via
      ``np.interp``.  This mode is used when boundaries come from a labeled
      mask image (``Orders_mask.fits``) rather than from a parametric fit.

    When Chebyshev width coefficients (``w_up_cheb_coef`` / ``w_low_cheb_coef``)
    are present, ``get_lower()`` / ``get_upper()`` reconstruct boundaries as::

        Y_upper(x) = center(x) + chebval(x, w_up_cheb_coef)
        Y_lower(x) = center(x) - chebval(x, w_low_cheb_coef)

    Otherwise dense arrays or polynomial coefficients are used as fallback.
    """

    aperture: int              # Aperture number
    order: int                 # Physical diffraction order

    # Polynomial coefficients (descending order: c0*x^n + … + cn)
    center_coef: np.ndarray = None
    lower_coef: np.ndarray = None
    upper_coef: np.ndarray = None

    # Dense per-column arrays (index = column number)
    center_arr: np.ndarray = None
    lower_arr: np.ndarray = None
    upper_arr: np.ndarray = None

    # Chebyshev width coefficients (from hybrid smoothing).
    # W_up(x)  = chebval(x, w_up_cheb_coef)   → upper = center + W_up
    # W_low(x) = chebval(x, w_low_cheb_coef)  → lower = center − W_low
    w_up_cheb_coef: np.ndarray = None
    w_low_cheb_coef: np.ndarray = None

    width: float = 30.0        # Aperture width (pixels)
    is_chebyshev: bool = False # Whether coefs use Chebyshev basis
    domain: Tuple[float, float] = (0, 4096) # Domain for Chebyshev evaluation
    
    # Advanced Missing Order Recovery Flags
    is_interpolated: bool = False # True if the entire order was mathematically interpolated
    
    # ------------------------------------------------------------------
    def _eval(self, coef, arr, x):
        """Return values at column positions *x*.

        Uses dense array when available (``np.interp``), otherwise uses
        Chebyshev polynomial evaluation (``chebval``). Standard power-series
        polynomials (np.polyval) are no longer supported for trace storage.
        """
        x = np.asarray(x, dtype=float)
        if arr is not None:
            # Dense array mode (evaluated via interpolation)
            idx = np.arange(len(arr), dtype=float)
            return np.interp(x, idx, arr)
        if coef is not None:
            # Chebyshev basis mode.
            # Coefficients produced by Chebyshev.fit(..., domain=(d0, d1)) live in
            # the canonical [-1, 1] window.  chebval() always evaluates in the
            # canonical window, so we must map raw column coordinates first.
            if self.is_chebyshev and self.domain is not None:
                d0, d1 = float(self.domain[0]), float(self.domain[1])
                if d1 > d0:
                    x_eval = 2.0 * (x - d0) / (d1 - d0) - 1.0
                else:
                    x_eval = x
                return chebval(x_eval, coef)
            return chebval(x, coef)
        return np.full_like(x, np.nan)

    def get_position(self, y: np.ndarray) -> np.ndarray:
<<<<<<< HEAD
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
=======
        """Get aperture center position at column coordinates *y*."""
        return self._eval(self.center_coef, self.center_arr, y)

    def get_lower(self, y: np.ndarray) -> np.ndarray:
        """Get lower boundary at column coordinates *y*."""
        if self.w_low_cheb_coef is not None:
            from numpy.polynomial.chebyshev import chebval
            x = np.asarray(y, dtype=float)
            center = self.get_position(x)
            w_low = np.maximum(chebval(x, self.w_low_cheb_coef), 0.0)
            return center - w_low
        return self._eval(self.lower_coef, self.lower_arr, y)

    def get_upper(self, y: np.ndarray) -> np.ndarray:
        """Get upper boundary at column coordinates *y*."""
        if self.w_up_cheb_coef is not None:
            from numpy.polynomial.chebyshev import chebval
            x = np.asarray(y, dtype=float)
            center = self.get_position(x)
            w_up = np.maximum(chebval(x, self.w_up_cheb_coef), 0.0)
            return center + w_up
        return self._eval(self.upper_coef, self.upper_arr, y)
>>>>>>> cef6f04 (	modified:   README.md)


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

<<<<<<< HEAD
=======
    def shift_orders(self, delta_m: int):
        """Shift all aperture IDs and order numbers by a constant offset."""
        if delta_m == 0:
            return
        new_apertures = {}
        for old_id, ap in self.apertures.items():
            new_id = old_id + delta_m
            ap.aperture = new_id
            ap.order = new_id
            new_apertures[new_id] = ap
        self.apertures = new_apertures

    @staticmethod
    def from_labeled_mask(order_labels: np.ndarray) -> 'ApertureSet':
        """Build an ApertureSet from a labeled order mask image.

        Parameters
        ----------
        order_labels : np.ndarray[int16], shape (rows, cols)
            Pixel value = order ID, 0 = background.

        Returns
        -------
        ApertureSet with dense-array ApertureLocation objects.
        """
        rows, cols = order_labels.shape
        unique_ids = sorted(int(v) for v in np.unique(order_labels) if v != 0)
        aps = ApertureSet()

        for oid in unique_ids:
            order_mask = (order_labels == oid)
            lo_arr = np.full(cols, np.nan, dtype=np.float64)
            hi_arr = np.full(cols, np.nan, dtype=np.float64)

            for xi in range(cols):
                col_rows = np.where(order_mask[:, xi])[0]
                if col_rows.size > 0:
                    lo_arr[xi] = float(col_rows[0])
                    hi_arr[xi] = float(col_rows[-1])

            center_arr = 0.5 * (lo_arr + hi_arr)
            valid = np.isfinite(lo_arr)
            w = float(np.nanmedian(hi_arr[valid] - lo_arr[valid])) if np.any(valid) else 10.0

            ap = ApertureLocation(
                aperture=oid, order=oid,
                center_arr=center_arr,
                lower_arr=lo_arr,
                upper_arr=hi_arr,
                width=w,
            )
            aps.add_aperture(ap)

        return aps

    def load_trace_coefs(self, coefs_path: str):
        """Load trace centers and Chebyshev width coefficients from a companion JSON file.

        For each order present in both the file and this ApertureSet, sets
        ``w_up_cheb_coef`` and ``w_low_cheb_coef`` on the ApertureLocation,
        and replaces ``center_arr`` with the B-spline-smoothed centre array
        stored in the file (sub-pixel precision).
        """
        import json as _json
        from pathlib import Path as _Path

        p = _Path(coefs_path)
        if not p.exists():
            return
        with open(p, 'r') as f:
            data = _json.load(f)
        for oid_str, entry in data.get('orders', {}).items():
            oid = int(oid_str)
            ap = self.apertures.get(oid)
            if ap is None:
                continue
            if 'center_arr' in entry:
                ap.center_arr = np.array(entry['center_arr'], dtype=np.float64)
            if 'w_up_cheb' in entry:
                ap.w_up_cheb_coef = np.array(entry['w_up_cheb'], dtype=np.float64)
            if 'w_low_cheb' in entry:
                ap.w_low_cheb_coef = np.array(entry['w_low_cheb'], dtype=np.float64)

>>>>>>> cef6f04 (	modified:   README.md)

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

    def shift_orders(self, delta_m: int):
        """Shift all aperture IDs and order numbers by a constant offset."""
        if delta_m == 0:
            return
        new_spectra = {}
        for old_id, spec in self.spectra.items():
            new_id = old_id + delta_m
            spec.aperture = new_id
            spec.order = new_id
            new_spectra[new_id] = spec
        self.spectra = new_spectra

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

    def shift_orders(self, delta_m: int):
        """Shift all order-keyed dictionaries by a constant offset."""
        if delta_m == 0:
            return
        if self.cross_profiles:
            self.cross_profiles = {k + delta_m: v for k, v in self.cross_profiles.items()}
        if self.blaze_profiles:
            self.blaze_profiles = {k + delta_m: v for k, v in self.blaze_profiles.items()}
        if self.flat_spec:
            self.flat_spec.shift_orders(delta_m)
        # aperture_set 统一在外部进行 shift，这里不重复处理

@dataclass
class WaveCalib:
    """Wavelength calibration solution."""

    # 2D polynomial coefficients for wavelength
    # wavelength(x, y) = sum(p_ij * x^i * y^j)
    poly_coef: np.ndarray     # Polynomial coefficients, shape (xorder+1, yorder+1)

    xorder: int = 3           # Wavelength poly order in X
    yorder: int = 3           # Wavelength poly order in Y
    poly_type: str = 'chebyshev'                  # Polynomial type: 'chebyshev', 'legendre', or 'polynomial'
    domain_x: Tuple[float, float] = (0.0, 4096.0) # Normalization domain for X (pixels)
    domain_y: Tuple[float, float] = (10.0, 100.0) # Normalization domain for Y (orders)

    # Identified calibration lines
    line_catalog: Optional[np.ndarray] = None  # Reference wavelengths
    line_pixels: Optional[np.ndarray] = None   # Detected pixel positions
    line_orders: Optional[np.ndarray] = None   # Order numbers of lines

    # Quality metrics
    rms: float = 0.0          # RMS residual [Angstrom]
    nlines: int = 0           # Number of identified lines

    # Metadata
    calib_type: str = "ThAr"  # Calibration lamp type
    delta_m: int = 0          # Order offset applied to match references

    def apply_to_pixel(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        """Apply wavelength calibration to pixel coordinates.

        Args:
            x: X pixel coordinates
            y: Y pixel coordinates (order number)

        Returns:
            Wavelength array
        """
        xmin, xmax = self.domain_x
        ymin, ymax = self.domain_y
        
        # 将坐标归一化到 [-1, 1]
        x_norm = 2.0 * (x - xmin) / max(xmax - xmin, 1e-6) - 1.0
        y_norm = 2.0 * (y - ymin) / max(ymax - ymin, 1e-6) - 1.0
        
        # 使用对应的 2D 多项式求 m * lambda
        if self.poly_type == 'chebyshev':
            from numpy.polynomial.chebyshev import chebval2d
            m_lambda = chebval2d(x_norm, y_norm, self.poly_coef)
        elif self.poly_type == 'legendre':
            from numpy.polynomial.legendre import legval2d
            m_lambda = legval2d(x_norm, y_norm, self.poly_coef)
        else:
            from numpy.polynomial.polynomial import polyval2d
            m_lambda = polyval2d(x_norm, y_norm, self.poly_coef)
            
        # 最终波长 lambda = m_lambda / m
        return m_lambda / y


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
