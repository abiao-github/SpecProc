"""
Helpers for identifying lamp lines in long-slit spectra.

Shared by the interactive identification window (``specproc.utils.identify``,
run as ``specproc identify``) and the slit pipeline. An *anchor* is a manually
identified (pixel, wavelength) pair; anchors seed the dispersion solution, which
is then extended automatically by matching all detected lamp peaks against the
full line catalogue.
"""

import itertools
import logging
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
from numpy.polynomial import Chebyshev, Legendre, Polynomial
from astropy.io import fits
from scipy import ndimage

from specproc.core.peak_fitting import fit_peak_center, fit_peak_shape

logger = logging.getLogger(__name__)

POLY_TYPES = ('polynomial', 'legendre', 'chebyshev')
_POLY_ALIASES = {
    'poly': 'polynomial',
    'ordinary': 'polynomial',
    'leg': 'legendre',
    'cheb': 'chebyshev',
}


def global_clipped_sigma(s: np.ndarray, clip: float = 3.0, n_iter: int = 20) -> float:
    """Robust global noise: MAD-based sigma with iterative sigma clipping."""
    s = np.asarray(s, dtype=float)
    use = np.isfinite(s)
    sigma = 0.0
    for _ in range(n_iter):
        med = np.median(s[use])
        sigma = 1.4826 * np.median(np.abs(s[use] - med))
        new = np.isfinite(s) & (np.abs(s - med) < clip * sigma)
        if (new == use).all() or new.sum() < 10:
            break
        use = new
    return float(sigma)


def find_arc_peaks(spec: np.ndarray, snr_threshold: float = 8.0, with_snr: bool = False,
                   noise_mode: str = 'local'):
    """Detect emission peaks above snr_threshold (peak height / background noise).

    Returns (centroid pixel (0-based), amplitude) tuples, or
    (pixel, amplitude, snr) tuples when with_snr is True.  The centroid is a
    local Gaussian fit (weighted centroid fallback).  noise_mode='global' uses one
    sigma for the whole spectrum: the MAD of the continuum-subtracted spectrum,
    iteratively 3-sigma clipped so that the emission lines do not inflate it.
    """
    spec = np.nan_to_num(np.asarray(spec, dtype=float))
    s = spec - ndimage.minimum_filter1d(ndimage.gaussian_filter1d(spec, 3), 61)
    s = ndimage.gaussian_filter1d(s, 1.0)
    glob = 1.4826 * np.median(np.abs(s - np.median(s)))
    # Local, one-sided noise: emission lines only add positive outliers, so the gap
    # between the running median and 16th percentile (~1 sigma) ignores them.
    size = min(101, len(s) // 2 * 2 - 1) if len(s) > 4 else 1
    local = ndimage.percentile_filter(s, 50, size=size) - ndimage.percentile_filter(s, 16, size=size)
    noise_arr = np.maximum(local, 0.5 * glob) + 1e-9
    if noise_mode == 'global':
        noise_arr = np.full(len(s), global_clipped_sigma(s) + 1e-9)
    pk = [i for i in range(2, len(s) - 2)
          if s[i] > s[i - 1] and s[i] >= s[i + 1] and s[i] > snr_threshold * noise_arr[i] and s[i] > s[i - 2]]
    out = []
    grid = np.arange(len(s), dtype=float)
    pk = np.asarray(pk, dtype=int)
    for k, i in enumerate(pk):
        left_width = right_width = 4
        if k > 0:
            left_width = min(4, max(2, int((i - pk[k - 1]) // 2)))
        if k + 1 < len(pk):
            right_width = min(4, max(2, int((pk[k + 1] - i) // 2)))
        centre = fit_peak_center(grid, s, int(i), half_width=4,
                                 left_width=left_width, right_width=right_width)
        entry = (centre, float(s[i]))
        out.append(entry + (float(s[i] / noise_arr[i]),) if with_snr else entry)
    return out


def measure_resolution(pixel, flux, peaks, fit,
                       min_fwhm: float = 0.5, max_fwhm: float = 15.0):
    """Estimate resolving power R = lambda / FWHM from strong arc lines.

    Each peak is re-fitted locally with a Gaussian; its FWHM in pixels is
    converted to Angstrom with the local derivative of the wavelength solution.
    The median over isolated, well-fitted lines is returned with the median
    FWHM.  Returns ``None`` when no reliable line is available.
    """
    pixel = np.asarray(pixel, dtype=float).reshape(-1)
    flux = np.asarray(flux, dtype=float).reshape(-1)
    if fit is None or len(pixel) != len(flux) or len(pixel) < 4:
        return None
    centres = []
    for p in peaks:
        value = getattr(p, 'pixel', p[0] if isinstance(p, (tuple, list)) else p)
        try:
            value = float(value)
        except (TypeError, ValueError):
            continue
        if np.isfinite(value):
            centres.append(value)
    centres = np.asarray(sorted(centres), dtype=float)
    if not len(centres):
        return None
    coef = dispersion_polyval_coef(fit)
    if coef is None or not np.all(np.isfinite(coef)):
        return None
    deriv = np.polyder(np.asarray(coef, dtype=float))
    results = []
    for k, centre in enumerate(centres):
        i = int(np.clip(round(centre), 0, len(pixel) - 1))
        left_width = right_width = 4
        if k > 0:
            left_width = min(4, max(2, int((centre - centres[k - 1]) // 2)))
        if k + 1 < len(centres):
            right_width = min(4, max(2, int((centres[k + 1] - centre) // 2)))
        shape = fit_peak_shape(pixel, flux, i, half_width=4,
                               left_width=left_width, right_width=right_width)
        if shape is None:
            continue
        c, sigma = shape
        fwhm_px = 2.354820045 * sigma
        if not (min_fwhm <= fwhm_px <= max_fwhm):
            continue
        slope = abs(float(np.polyval(deriv, c)))
        wave = float(fit(c))
        if not (np.isfinite(slope) and slope > 0 and np.isfinite(wave) and wave > 0):
            continue
        fwhm_ang = fwhm_px * slope
        if not np.isfinite(fwhm_ang) or fwhm_ang <= 0:
            continue
        results.append((wave, fwhm_px, fwhm_ang, wave / fwhm_ang))
    if not results:
        return None
    arr = np.asarray(results, dtype=float)
    # R normally varies smoothly with wavelength.  Reject outliers relative
    # to a robust linear trend in log(R) vs wavelength, not relative to a
    # single global median (which would reject the genuine trend ends).
    keep = np.ones(len(arr), dtype=bool)
    trend_waves = np.empty(0, dtype=float)
    trend_res = np.empty(0, dtype=float)
    if len(arr) >= 5:
        x = arr[:, 0]
        log_r = np.log10(np.clip(arr[:, 3], 1e-30, None))
        coef = None
        for _ in range(4):
            if keep.sum() < 4:
                break
            xk, yk = x[keep], log_r[keep]
            # Theil-Sen style median slope: an endpoint outlier must not tilt
            # the trend enough to protect itself.
            dx = xk[:, None] - xk[None, :]
            dy = yk[:, None] - yk[None, :]
            upper = np.triu_indices(len(xk), 1)
            valid = np.abs(dx[upper]) > 1e-9
            if not valid.any():
                break
            slopes = dy[upper][valid] / dx[upper][valid]
            slope = float(np.median(slopes))
            intercept = float(np.median(yk - slope * xk))
            resid = log_r - (slope * x + intercept)
            centre = float(np.median(resid[keep]))
            mad = 1.4826 * float(np.median(np.abs(resid[keep] - centre)))
            if not np.isfinite(mad) or mad <= 0:
                break
            new_keep = np.abs(resid - centre) <= 3.0 * mad
            if new_keep.sum() < 3:
                break
            if np.array_equal(new_keep, keep):
                break
            keep = new_keep
            coef = (slope, intercept)
        if coef is None:
            # A straight line through a small set of points.
            if len(arr) >= 2:
                c = np.polyfit(x, log_r, 1)
                coef = (float(c[0]), float(c[1]))
        if coef is not None:
            gx = np.linspace(float(x[keep].min()), float(x[keep].max()), 100)
            gy = coef[0] * gx + coef[1]
            trend_waves = gx
            trend_res = np.power(10.0, gy)
    kept = arr[keep]
    rejected = arr[~keep]
    resolution = float(np.median(kept[:, 3]))
    fwhm_ang = float(np.median(kept[:, 2]))
    fwhm_px = float(np.median(kept[:, 1]))
    return {'resolution': resolution,
            'fwhm_ang': fwhm_ang,
            'fwhm_px': fwhm_px,
            'n': int(len(kept)),
            'n_all': int(len(arr)),
            'n_rejected': int(len(rejected)),
            'waves': kept[:, 0].copy(),
            'fwhm_ang_lines': kept[:, 2].copy(),
            'resolutions': kept[:, 3].copy(),
            'rejected_waves': rejected[:, 0].copy(),
            'rejected_resolutions': rejected[:, 3].copy(),
            'trend_waves': trend_waves,
            'trend_resolutions': trend_res}


def load_linelist(path: str) -> Tuple[np.ndarray, List[str], np.ndarray]:
    """Read a catalogue: wavelength, label and relative strength (1 when absent).

    The first column is the wavelength; if the second is numeric it is the line
    strength, otherwise the remaining text is the label.
    """
    waves, labels, strengths = [], [], []
    for line in Path(path).read_text(errors='ignore').splitlines():
        text = line.strip()
        if not text or text.startswith('#'):
            continue
        parts = text.replace(',', ' ').split(None, 1)
        try:
            waves.append(float(parts[0]))
        except ValueError:
            continue
        rest = parts[1].strip() if len(parts) > 1 else ''
        try:
            strengths.append(float(rest.split()[0]) if rest else 1.0)
            labels.append('')
        except ValueError:
            strengths.append(1.0)
            labels.append(rest)
    order = np.argsort(waves)
    return (np.asarray(waves, dtype=float)[order], [labels[i] for i in order],
            np.asarray(strengths, dtype=float)[order])


def _collapse_image(img: np.ndarray, dispersion_axis: str) -> np.ndarray:
    """Sum a 2-D lamp image over a band around the brightest spatial position."""
    img = np.nan_to_num(np.asarray(img, dtype=float))
    if dispersion_axis == 'x':
        img = img.T
    profile = ndimage.gaussian_filter1d(img.sum(axis=0), 3)
    c = int(np.argmax(profile))
    lo, hi = max(c - 12, 0), min(c + 13, img.shape[1])
    return img[:, lo:hi].sum(axis=1)


def load_lamp_spectrum(path: str, dispersion_axis: str = 'y') -> np.ndarray:
    """Load a 1-D lamp spectrum (counts versus 0-based pixel).

    Accepts a SpecProc ``*_pixel.fits`` / ``*_wave.fits`` table (FLUX_SUM/FLUX_OPT), a 1-D or 2-D
    FITS image, or a text/CSV file (last column is the flux).
    """
    p = Path(path)
    if p.suffix.lower() in ('.fits', '.fit', '.fts'):
        with fits.open(p) as hdul:
            for hdu in hdul:
                names = [n.upper() for n in (getattr(hdu, 'columns', None).names
                                              if getattr(hdu, 'columns', None) is not None else [])]
                for col in ('FLUX_SUM', 'FLUX_OPT', 'FLUX'):
                    if col in names:
                        return np.nan_to_num(np.asarray(hdu.data[col], dtype=float))
            for hdu in hdul:
                data = getattr(hdu, 'data', None)
                if data is not None and getattr(data, 'ndim', 0) in (1, 2) and data.dtype.names is None:
                    data = np.asarray(data, dtype=float)
                    return np.nan_to_num(data) if data.ndim == 1 else _collapse_image(
                        data, dispersion_axis.lower())
        raise ValueError(f"No spectrum or image data found in {path}")
    rows = []
    for line in p.read_text(errors='ignore').splitlines():
        text = line.split('#')[0].strip()
        if not text:
            continue
        try:
            rows.append([float(v) for v in text.replace(',', ' ').split()])
        except ValueError:
            continue
    if not rows:
        raise ValueError(f"No numeric data found in {path}")
    if len(rows) == 1 and len(rows[0]) > 2:
        # Common CSV/TXT export: all flux values on a single comma/space line.
        return np.asarray(rows[0], dtype=float)
    return np.array([r[-1] for r in rows], dtype=float)


def load_anchors(path: str) -> np.ndarray:
    """Read anchors as an (N, 2) array of (0-based pixel, wavelength)."""
    rows = []
    for line in Path(path).read_text().splitlines():
        text = line.split('#')[0].strip()
        if not text:
            continue
        parts = text.replace(',', ' ').split()
        try:
            rows.append((float(parts[0]), float(parts[1])))
        except (ValueError, IndexError):
            continue
    return np.array(rows, dtype=float).reshape(-1, 2)


def save_anchors(path: str, pairs: np.ndarray, note: str = ''):
    """Write anchors as CSV, readable by the slit pipeline ``identification_file`` option."""
    pairs = np.asarray(pairs, dtype=float).reshape(-1, 2)
    pairs = pairs[np.argsort(pairs[:, 0])]
    lines = ['# SpecProc anchor lines for long-slit wavelength calibration',
             '# pixel (0-based along the dispersion axis), wavelength (Angstrom)']
    if note:
        lines.append('# ' + note)
    lines += ['%.3f, %.4f' % (x, w) for x, w in pairs]
    Path(path).write_text('\n'.join(lines) + '\n')


def fit_polynomial(pairs: np.ndarray, order: int) -> Optional[np.ndarray]:
    """Fit wavelength = poly(pixel); the order is limited by the number of pairs."""
    pairs = np.asarray(pairs, dtype=float).reshape(-1, 2)
    if len(pairs) < 2:
        return None
    return np.polyfit(pairs[:, 0], pairs[:, 1], min(order, len(pairs) - 1))


def normalize_poly_type(poly_type: Optional[str]) -> str:
    """Return one of ``polynomial``, ``legendre`` or ``chebyshev``."""
    name = str(poly_type or 'legendre').strip().lower()
    name = _POLY_ALIASES.get(name, name)
    if name not in POLY_TYPES:
        raise ValueError(f'Unknown polynomial type {poly_type!r}; '
                         f'choose one of {", ".join(POLY_TYPES)}')
    return name


def fit_dispersion(x, y, order: int, poly_type: str = 'legendre',
                   domain: Optional[Tuple[float, float]] = None):
    """Fit a 1-D pixel -> wavelength polynomial in the requested basis.

    Returns a :mod:`numpy.polynomial` object (``Legendre``, ``Chebyshev`` or
    ``Polynomial``) that can be called directly, or ``None`` when fewer than
    two finite points are available.  The order is limited by the number of
    points.  Fitting in an orthogonal basis is numerically better conditioned
    for degree 5 and above than ``numpy.polyfit`` on raw pixel values.
    """
    x = np.asarray(x, dtype=float).reshape(-1)
    y = np.asarray(y, dtype=float).reshape(-1)
    good = np.isfinite(x) & np.isfinite(y)
    x, y = x[good], y[good]
    if len(x) < 2:
        return None
    order = max(1, min(int(order), len(x) - 1))
    if domain is None:
        domain = (float(np.min(x)), float(np.max(x)))
    else:
        domain = (float(domain[0]), float(domain[1]))
    if not (domain[1] > domain[0]):
        domain = (domain[0] - 1.0, domain[1] + 1.0)

    kind = normalize_poly_type(poly_type)
    cls = {'polynomial': Polynomial, 'legendre': Legendre, 'chebyshev': Chebyshev}[kind]
    return cls.fit(x, y, order, domain=domain)


def dispersion_type(fit) -> str:
    """Return the basis name of a :func:`fit_dispersion` result."""
    if isinstance(fit, Legendre):
        return 'legendre'
    if isinstance(fit, Chebyshev):
        return 'chebyshev'
    return 'polynomial'


def dispersion_polyval_coef(fit) -> Optional[np.ndarray]:
    """Convert a fitted series to ``numpy.polyval`` (descending) coefficients."""
    if fit is None:
        return None
    return np.asarray(fit.convert(kind=Polynomial).coef[::-1], dtype=float)


def dispersion_formula(fit, precision: int = 6) -> str:
    """Human-readable formula of a fitted dispersion series."""
    if fit is None:
        return ''
    kind = dispersion_type(fit)
    degree = fit.degree()
    terms = ' '.join(f'{c:+.{precision}g}*{name}'
                     for c, name in zip(fit.coef, _basis_names(kind, degree)))
    return f'{kind} degree {degree}: {terms}'


def _basis_names(kind: str, degree: int) -> List[str]:
    if kind == 'legendre':
        return [f'P{k}(t)' for k in range(degree + 1)]
    if kind == 'chebyshev':
        return [f'T{k}(t)' for k in range(degree + 1)]
    return ['1'] + [f't^{k}' for k in range(1, degree + 1)]


def _match_nearest(pix: np.ndarray, catalog: np.ndarray, model, tol: float) -> np.ndarray:
    """Pair each peak with its nearest catalogue line within tol (each line used once)."""
    wl = model(pix)
    idx = np.clip(np.searchsorted(catalog, wl), 1, len(catalog) - 1)
    near = np.where(np.abs(catalog[idx - 1] - wl) < np.abs(catalog[idx] - wl), idx - 1, idx)
    dist = np.abs(catalog[near] - wl)
    best = {}
    for i in np.where(dist < tol)[0]:
        j = near[i]
        if j not in best or dist[i] < dist[best[j]]:
            best[j] = i
    sel = np.array(sorted(best.values()), dtype=int)
    if len(sel) == 0:
        return np.empty((0, 2))
    return np.column_stack([pix[sel], catalog[near[sel]]])


def pattern_match_linear(peak_pix, catalog, peak_amp=None, wave_range=None, pixel_range=None,
                         n_peaks: int = 15, max_catalog: int = 150, ratio_eps: float = 0.01):
    """Find a linear pixel->wavelength map from the pattern of line spacings.

    For three lines the ratio (middle-first)/(last-first) is unchanged by a linear
    dispersion, in either direction. Triples of the strongest lamp peaks are paired
    with catalogue triples of equal ratio; each pairing gives a linear hypothesis
    that is scored by how many lamp peaks land on a catalogue line. When
    ``wave_range`` and ``pixel_range`` are given, hypotheses whose wavelength
    coverage falls outside the range are discarded. Returns the coefficients
    ``[scale, offset]`` of the best hypothesis, or None.
    """
    pix = np.asarray(peak_pix, dtype=float)
    cat = np.sort(np.asarray(catalog, dtype=float))
    if len(pix) < 3 or len(cat) < 3:
        return None
    amp = np.ones(len(pix)) if peak_amp is None else np.asarray(peak_amp, dtype=float)
    key = pix[np.argsort(-amp)[:n_peaks]]
    key.sort()
    span = pix.max() - pix.min()
    if pixel_range is None:
        pixel_range = (float(pix.min()), float(pix.max()))
    # Catalogue lines used for hashing: evenly thinned if too many.
    sub = cat if len(cat) <= max_catalog else cat[np.linspace(0, len(cat) - 1, max_catalog).astype(int)]
    ia, ib, ic = np.array(list(itertools.combinations(range(len(sub)), 3))).T
    ca, cb, cc = sub[ia], sub[ib], sub[ic]
    cr = (cb - ca) / (cc - ca)
    order_ = np.argsort(cr)
    ca, cc, cr = ca[order_], cc[order_], cr[order_]

    # Collect the hypotheses first and score them in blocks.  Scoring one
    # triple at a time needs one searchsorted call per hypothesis; with a dense
    # catalogue that is tens of thousands of tiny calls and dominates run time.
    scales, offsets = [], []
    min_span = 0.15 * span
    for i, j, k in itertools.combinations(range(len(key)), 3):
        p1, p2, p3 = key[i], key[j], key[k]
        if p3 - p1 < min_span:
            continue
        r = (p2 - p1) / (p3 - p1)
        for reverse in (False, True):
            target = 1 - r if reverse else r
            lo, hi = np.searchsorted(cr, [target - ratio_eps, target + ratio_eps])
            if hi <= lo:
                continue
            w1, w3 = (cc[lo:hi], ca[lo:hi]) if reverse else (ca[lo:hi], cc[lo:hi])
            scale = (w3 - w1) / (p3 - p1)
            offset = w1 - scale * p1
            if wave_range is not None:
                wlo, whi = sorted(map(float, wave_range))
                margin = 0.05 * (whi - wlo)
                e0 = scale * pixel_range[0] + offset
                e1 = scale * pixel_range[1] + offset
                ok = (np.minimum(e0, e1) >= wlo - margin) & (np.maximum(e0, e1) <= whi + margin)
                scale, offset = scale[ok], offset[ok]
                if len(scale) == 0:
                    continue
            scales.append(scale)
            offsets.append(offset)

    if not scales:
        return None
    scale_all = np.concatenate(scales)
    offset_all = np.concatenate(offsets)
    best, best_key = None, (0, 0.0)
    block = 4096
    for start in range(0, len(scale_all), block):
        scale = scale_all[start:start + block]
        offset = offset_all[start:start + block]
        wl = scale[:, None] * pix[None, :] + offset[:, None]
        idx = np.searchsorted(cat, wl.ravel())
        idx = np.clip(idx, 1, len(cat) - 1).reshape(wl.shape)
        dist = np.minimum(np.abs(cat[idx] - wl), np.abs(cat[idx - 1] - wl))
        tol_w = np.maximum(1.5 * np.abs(scale), 1e-3)[:, None]
        hit = dist < tol_w
        score = hit.sum(axis=1)
        rms = np.sqrt((np.where(hit, dist, 0) ** 2).sum(axis=1) / np.maximum(score, 1))
        # np.argmax keeps the first maximum, so ties favour the earlier triple
        # exactly as the old per-hypothesis loop did.
        n = int(np.argmax(score - 1e-6 * rms))
        cand = (int(score[n]), -float(rms[n]))
        if cand > best_key:
            best_key, best = cand, np.array([scale[n], offset[n]])
    if best is None or best_key[0] < 4:
        return None
    return best


def refine_with_catalog(peak_pix, anchors: np.ndarray, catalog: np.ndarray, order: int,
                        tol: float = 3.0, min_tol: float = 0.4, max_iter: int = 6,
                        initial_range: Optional[Tuple[float, float]] = None,
                        initial_pixel_range: Optional[Tuple[float, float]] = None,
                        peak_amp=None, initial_coef=None,
                        poly_type: str = 'polynomial'):
    """Match detected lamp peaks to a catalogue, optionally without manual anchors.

    With anchors, they give the rough solution and are always retained. With no
    anchors, the starting solution comes from line-spacing pattern matching
    (``pattern_match_linear``). Matching then uses a shrinking tolerance while the
    polynomial order grows.  ``poly_type`` selects the basis used for the
    successive fits (``polynomial``, ``legendre`` or ``chebyshev``);
    ``initial_coef`` (ordinary polynomial coefficients, pixel -> wavelength)
    replaces the starting solution, e.g. a hand-aligned overlay.
    Returns (pairs, ordinary polynomial coefficients) so existing callers can
    keep using ``numpy.polyval``.
    """
    anchors = np.asarray(anchors, dtype=float).reshape(-1, 2)
    pix = np.asarray(peak_pix, dtype=float)
    catalog = np.asarray(catalog, dtype=float)
    catalog = np.sort(catalog[np.isfinite(catalog)])
    domain = ((float(np.min(pix)), float(np.max(pix))) if len(pix) > 1 else
              ((float(np.min(anchors[:, 0])), float(np.max(anchors[:, 0])))
               if len(anchors) > 1 else (0.0, 1.0)))
    if domain[1] <= domain[0]:
        domain = (domain[0] - 1.0, domain[1] + 1.0)

    def _fit(x, y, deg):
        return fit_dispersion(x, y, deg, poly_type, domain=domain)

    def _coef(fit):
        return dispersion_polyval_coef(fit) if fit is not None else None

    def _fallback():
        fit = _fit(anchors[:, 0], anchors[:, 1], order) if len(anchors) >= 2 else None
        coef = _coef(fit)
        return anchors, coef if coef is not None else fit_polynomial(anchors, order)

    if len(catalog) < 2 or len(pix) == 0:
        return _fallback()

    if len(anchors) == 0:
        seed = (np.asarray(initial_coef, dtype=float) if initial_coef is not None else
                pattern_match_linear(pix, catalog, peak_amp, initial_range,
                                     initial_pixel_range))
        if seed is None:
            return np.empty((0, 2)), None
        model = Polynomial(np.asarray(seed, dtype=float)[::-1])
        pairs = np.empty((0, 2))
        for it in range(max_iter):
            cur_tol = max(tol / (it + 1), min_tol)
            matched = _match_nearest(pix, catalog, model, cur_tol)
            if len(matched) < 2:
                break
            pairs = matched
            model = _fit(pairs[:, 0], pairs[:, 1], min(order, it + 1))
            if model is None:
                break
        if len(pairs) < 2 or model is None:
            return np.empty((0, 2)), None
        final_fit = _fit(pairs[:, 0], pairs[:, 1], order)
        if final_fit is None:
            return np.empty((0, 2)), None
        return pairs[np.argsort(pairs[:, 0])], _coef(final_fit)

    model = (Polynomial(np.asarray(initial_coef, dtype=float)[::-1])
             if initial_coef is not None
             else _fit(anchors[:, 0], anchors[:, 1], 1))
    if model is None:
        return _fallback()

    pairs = anchors
    for it in range(max_iter):
        cur_tol = max(tol / (it + 1), min_tol)
        matched = _match_nearest(pix, catalog, model, cur_tol)
        if len(matched):
            near_anchor = np.min(np.abs(matched[:, :1] - anchors[:, 0][None, :]), axis=1) < 1.5
            same_line = np.isin(np.round(matched[:, 1], 4), np.round(anchors[:, 1], 4))
            matched = matched[~(near_anchor | same_line)]
        candidate = np.vstack([anchors, matched])
        if len(candidate) < 2:
            break
        pairs = candidate
        model = _fit(pairs[:, 0], pairs[:, 1], min(order, it + 1))
        if model is None:
            break
    final_fit = _fit(pairs[:, 0], pairs[:, 1], order)
    if final_fit is None:
        return _fallback()
    return pairs[np.argsort(pairs[:, 0])], _coef(final_fit)
