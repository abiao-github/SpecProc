"""
Automatic anchor-line matching for long-slit lamp spectra.

The module matches the emission peaks of an extracted lamp spectrum (pixel
space) against a high-resolution reference lamp spectrum or line catalogue,
fits a low-order pixel-to-wavelength solution and writes a small, clean anchor
file that can seed the full wavelength calibration in the slit pipeline.

The recommended workflow for the Muztagh-Ata 1.93 m slit spectrograph is:

1. provide the approximate wavelength coverage from the instrument
   configuration (``[slit] wave_start`` / ``wave_end``);
2. detect the strong lamp peaks in the extracted 1-D spectrum;
3. match them against the strong peaks of the NOIRLAB/KPNO FeAr atlas
   (``calib_data/linelists/fear.fits``), which is a 0.2 A/pixel, absolutely
   wavelength-calibrated spectrum;
4. keep only the clean, isolated, well-spread matches as anchors;
5. write ``pixel, wavelength`` rows plus a diagnostic figure for inspection.

The functions are independent of the GUI so they can also be called from
notebooks or tests.  The command line entry point is
:mod:`specproc.match_anchors` (``specproc-match-anchors``).
"""

import logging
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

import numpy as np
from astropy.io import fits
from scipy import ndimage
from scipy.signal import find_peaks

from specproc.core.line_identification import load_linelist, pattern_match_linear

logger = logging.getLogger(__name__)


class Peak:
    """A detected emission peak (0-based pixel, prominence, S/N and width)."""

    __slots__ = ('pixel', 'amplitude', 'snr', 'width')

    def __init__(self, pixel: float, amplitude: float, snr: float, width: float):
        self.pixel = float(pixel)
        self.amplitude = float(amplitude)
        self.snr = float(snr)
        self.width = float(width)

    def __repr__(self) -> str:  # pragma: no cover - debugging aid
        return f'Peak(pixel={self.pixel:.2f}, snr={self.snr:.1f})'


class Reference:
    """Reference lamp spectrum or line catalogue in wavelength space."""

    def __init__(self, wave: np.ndarray, flux: Optional[np.ndarray] = None,
                 labels: Optional[Sequence[str]] = None, source: str = ''):
        self.wave = np.asarray(wave, dtype=float)
        self.flux = (np.ones_like(self.wave) if flux is None
                     else np.asarray(flux, dtype=float))
        self.labels = list(labels) if labels is not None else [''] * len(self.wave)
        self.source = source

    def __len__(self) -> int:
        return len(self.wave)

    def subset(self, wave_range: Optional[Sequence[float]],
               margin: float = 0.0) -> 'Reference':
        """Return the entries inside ``wave_range`` (optionally with a margin)."""
        if wave_range is None:
            return Reference(self.wave.copy(), self.flux.copy(),
                             list(self.labels), self.source)
        w0, w1 = sorted(float(v) for v in wave_range)
        m = (self.wave >= w0 - margin) & (self.wave <= w1 + margin)
        return Reference(self.wave[m], self.flux[m],
                         [self.labels[i] for i in np.where(m)[0]], self.source)


class Solution:
    """A pixel-to-wavelength solution and the pairs that produced it."""

    def __init__(self, coef: np.ndarray, pairs: np.ndarray, residuals: np.ndarray,
                 peak_indices: np.ndarray, catalog_indices: np.ndarray,
                 n_pixels: int):
        self.coef = np.asarray(coef, dtype=float)
        self.pairs = np.asarray(pairs, dtype=float).reshape(-1, 2)
        self.residuals = np.asarray(residuals, dtype=float)
        self.peak_indices = np.asarray(peak_indices, dtype=int)
        self.catalog_indices = np.asarray(catalog_indices, dtype=int)
        self.n_pixels = int(n_pixels)

    def wavelength(self, pixel) -> np.ndarray:
        return np.polyval(self.coef, pixel)

    @property
    def rms(self) -> float:
        return float(np.sqrt(np.mean(self.residuals ** 2))) if len(self.residuals) else np.nan

    @property
    def wave_range(self) -> Tuple[float, float]:
        return (float(self.wavelength(0.0)),
                float(self.wavelength(self.n_pixels - 1)))

    @property
    def n_lines(self) -> int:
        return len(self.pairs)


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------
def _as_float_array(data) -> np.ndarray:
    return np.nan_to_num(np.asarray(data, dtype=float), nan=0.0, posinf=0.0, neginf=0.0)


def _fits_wcs_wavelength(header) -> Optional[np.ndarray]:
    """Build a wavelength array from a simple 1-D FITS WCS, if present."""
    n = int(header.get('NAXIS1', 0) or 0)
    if n <= 0:
        return None
    crval = header.get('CRVAL1')
    cdelt = header.get('CDELT1', header.get('CD1_1'))
    if crval is None or cdelt is None:
        return None
    crpix = float(header.get('CRPIX1', 1.0))
    return float(crval) + (np.arange(n) + 1.0 - crpix) * float(cdelt)


def _wavelength_unit(header) -> str:
    return str(header.get('WAVEUNIT', header.get('CUNIT1', ''))).strip().lower()


def _is_angstrom(unit: str) -> bool:
    return unit.startswith(('ang', 'a')) or unit in ('', 'none')


def load_observed(path: str) -> Tuple[np.ndarray, np.ndarray, dict]:
    """Read an extracted 1-D lamp spectrum.

    Returns ``(flux, pixel, info)`` where ``pixel`` is 0-based.  Both SpecProc
    ``*_spec.fits`` tables and plain 1-D FITS arrays are accepted.  When the
    table carries an absolute wavelength column it is returned in
    ``info['wavelength']`` (and ``info['wave_range']``); the pixel axis is
    always returned as well because that is what the slit pipeline expects.
    """
    p = Path(path)
    if p.suffix.lower() not in ('.fits', '.fit', '.fts'):
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
            raise ValueError(f'No numeric spectrum found in {path}')
        arr = np.asarray(rows, dtype=float)
        flux = _as_float_array(arr[:, -1])
        return flux, np.arange(len(flux), dtype=float), {'path': str(p)}

    with fits.open(p) as hdul:
        for hdu in hdul:
            data = getattr(hdu, 'data', None)
            if data is None or getattr(data, 'dtype', None) is None or data.dtype.names is None:
                continue
            names = [n.upper() for n in data.dtype.names]
            flux_col = next((n for n in ('FLUX_SUM', 'FLUX_OPT', 'FLUX') if n in names), None)
            if flux_col is None:
                continue
            flux = _as_float_array(data[flux_col])
            pixel = np.arange(len(flux), dtype=float)
            info = {'path': str(p), 'wave_unit': _wavelength_unit(hdu.header),
                    'header': hdu.header}
            wave_col = next((n for n in ('WAVELENGTH', 'WAVE', 'LAMBDA') if n in names), None)
            if wave_col is not None:
                wave = _as_float_array(data[wave_col])
                unit = info['wave_unit']
                # SpecProc writes 1..N in WAVELENGTH when the spectrum is still
                # on a pixel axis; only trust it when the unit says Angstrom.
                if _is_angstrom(unit) and len(wave) > 2 and np.all(np.isfinite(wave)):
                    info['wavelength'] = wave
                    if wave[-1] > wave[0]:
                        info['wave_range'] = (float(wave[0]), float(wave[-1]))
            return flux, pixel, info
        for hdu in hdul:
            data = getattr(hdu, 'data', None)
            if data is None or getattr(data, 'ndim', 0) != 1 or data.dtype.names is not None:
                continue
            flux = _as_float_array(data)
            info = {'path': str(p), 'header': hdu.header}
            wave = _fits_wcs_wavelength(hdu.header)
            if wave is not None:
                info['wavelength'] = wave
                info['wave_range'] = (float(wave[0]), float(wave[-1]))
            return flux, np.arange(len(flux), dtype=float), info
    raise ValueError(f'No 1-D spectrum found in {path}')


def _read_text_reference(path: str) -> Reference:
    waves, fluxes, labels = [], [], []
    has_flux_column = None
    for line in Path(path).read_text(errors='ignore').splitlines():
        text = line.split('#')[0].strip()
        if not text:
            continue
        parts = text.replace(',', ' ').split()
        try:
            w = float(parts[0])
        except ValueError:
            continue
        rest = parts[1:]
        flux = None
        label = ''
        if rest:
            try:
                flux = float(rest[0])
                has_flux_column = True
            except ValueError:
                label = ' '.join(rest)
                has_flux_column = False if has_flux_column is None else has_flux_column
        waves.append(w)
        fluxes.append(flux)
        labels.append(label)
    if not waves:
        raise ValueError(f'No wavelength data found in {path}')
    if has_flux_column and all(f is not None for f in fluxes):
        flux_arr = np.asarray(fluxes, dtype=float)
    else:
        flux_arr = None
    return Reference(np.asarray(waves, dtype=float), flux_arr, labels, source=str(path))


def load_reference(path: str, linelist: Optional[str] = None) -> Reference:
    """Read a reference lamp spectrum or line catalogue.

    ``path`` may be a FITS 1-D spectrum (with a simple WCS), a FITS table with
    wavelength/flux columns, or a text file.  When ``linelist`` is given its
    wavelengths and ion labels are attached to the nearest reference entries.
    """
    p = Path(path)
    if p.suffix.lower() in ('.fits', '.fit', '.fts'):
        with fits.open(p) as hdul:
            ref = None
            for hdu in hdul:
                data = getattr(hdu, 'data', None)
                if data is None or getattr(data, 'dtype', None) is None or data.dtype.names is None:
                    continue
                names = [n.upper() for n in data.dtype.names]
                wave_col = next((n for n in ('WAVELENGTH', 'WAVE', 'LAMBDA') if n in names), None)
                flux_col = next((n for n in ('FLUX', 'FLUX_SUM', 'FLUX_OPT', 'COUNTS') if n in names), None)
                if wave_col is not None:
                    ref = Reference(_as_float_array(data[wave_col]),
                                    _as_float_array(data[flux_col]) if flux_col else None,
                                    source=str(p))
                    break
            if ref is None:
                for hdu in hdul:
                    data = getattr(hdu, 'data', None)
                    if data is None or getattr(data, 'ndim', 0) != 1 or data.dtype.names is not None:
                        continue
                    wave = _fits_wcs_wavelength(hdu.header)
                    if wave is None:
                        raise ValueError(
                            f'{path} has no wavelength WCS; supply --linelist or a '
                            'wavelength column')
                    ref = Reference(wave, _as_float_array(data), source=str(p))
                    break
            if ref is None:
                raise ValueError(f'No reference spectrum found in {path}')
    else:
        ref = _read_text_reference(str(p))

    if linelist:
        cat, labels, _ = load_linelist(linelist)
        if len(cat):
            # Attach the catalogue label to the closest reference entry.
            order = np.argsort(ref.wave)
            ref.wave = ref.wave[order]
            ref.flux = ref.flux[order]
            ref.labels = [ref.labels[i] for i in order]
            if len(ref.wave) > 1:
                idx = np.clip(np.searchsorted(ref.wave, cat), 1, len(ref.wave) - 1)
                near = np.where(np.abs(ref.wave[idx - 1] - cat) < np.abs(ref.wave[idx] - cat),
                                idx - 1, idx)
                for j, lab in zip(near, labels):
                    if not ref.labels[j]:
                        ref.labels[j] = lab
    order = np.argsort(ref.wave)
    ref.wave = ref.wave[order]
    ref.flux = ref.flux[order]
    ref.labels = [ref.labels[i] for i in order]
    return ref


# ---------------------------------------------------------------------------
# Peak detection and catalogue construction
# ---------------------------------------------------------------------------
def detect_peaks(flux: np.ndarray, snr: float = 8.0, min_width: float = 0.5,
                 max_width: float = 10.0, distance: int = 4,
                 baseline_size: int = 61) -> List[Peak]:
    """Detect narrow emission peaks on top of a smooth continuum.

    A median-filtered version of the smoothed spectrum is subtracted before
    peak finding, so broad blaze/continuum variations do not produce spurious
    peaks.  The returned pixel is a three-pixel weighted centroid.
    """
    flux = _as_float_array(flux)
    if len(flux) < 5:
        return []
    smooth = ndimage.gaussian_filter1d(flux, 1.0)
    size = min(baseline_size, len(smooth) // 2 * 2 + 1)
    baseline = ndimage.median_filter(smooth, size=max(size, 3))
    residual = smooth - baseline
    noise = 1.4826 * np.median(np.abs(residual - np.median(residual)))
    if not np.isfinite(noise) or noise <= 0:
        noise = np.std(residual) or 1.0
    idx, props = find_peaks(residual, prominence=max(snr * noise, 1e-9),
                            width=(min_width, max_width), distance=distance)
    peaks = []
    for k, i in enumerate(idx):
        lo, hi = max(i - 3, 0), min(i + 4, len(residual))
        w = np.clip(residual[lo:hi], 0, None)
        centroid = float((np.arange(lo, hi) * w).sum() / w.sum()) if w.sum() > 0 else float(i)
        peaks.append(Peak(centroid, props['prominences'][k],
                          props['prominences'][k] / noise, props['widths'][k]))
    return peaks


def _labels_for(waves: np.ndarray, label_wave: Optional[np.ndarray],
                label_text: Optional[Sequence[str]], tol: float = 0.8):
    """Attach the nearest catalogue species label to each reference peak."""
    labels = [''] * len(waves)
    if label_wave is None or label_text is None or len(label_wave) == 0:
        return labels
    label_wave = np.asarray(label_wave, dtype=float)
    order = np.argsort(label_wave)
    label_wave = label_wave[order]
    label_text = [label_text[i] for i in order]
    idx = np.clip(np.searchsorted(label_wave, waves), 1, len(label_wave) - 1)
    near = np.where(np.abs(label_wave[idx - 1] - waves) < np.abs(label_wave[idx] - waves),
                    idx - 1, idx)
    for k, j in enumerate(near):
        if abs(label_wave[j] - waves[k]) <= tol:
            labels[k] = label_text[j]
    return labels


def build_catalog(ref: Reference, wave_range: Optional[Sequence[float]] = None,
                  snr: float = 20.0, strength_frac: float = 0.01,
                  min_strength: float = 50.0, max_lines: int = 400,
                  linelist_wave: Optional[np.ndarray] = None,
                  linelist_labels: Optional[Sequence[str]] = None) -> Reference:
    """Build the strong-line catalogue used for the pattern match.

    When ``ref`` is a real spectrum its strongest peaks are used, with a
    threshold relative to the strongest line *inside* the requested range.
    Using a local threshold is important for the KPNO atlas: its throughput
    changes by more than an order of magnitude between echelle orders.
    """
    if ref.flux is None or np.allclose(ref.flux, 1.0):
        margin = 50.0
        if wave_range is not None:
            w0, w1 = sorted(float(v) for v in wave_range)
            margin = max(margin, 0.02 * (w1 - w0))
        cat = ref.subset(wave_range, margin=margin)
        labels = list(cat.labels)
        if linelist_wave is not None:
            labels = _labels_for(cat.wave, linelist_wave, linelist_labels)
        return Reference(cat.wave.copy(), np.ones_like(cat.wave), labels, ref.source)

    ref_peaks = detect_peaks(ref.flux, snr=snr, min_width=0.2, max_width=8.0,
                             distance=3, baseline_size=101)
    if not ref_peaks:
        raise ValueError('No peaks detected in the reference spectrum')
    cat_wave = np.interp([p.pixel for p in ref_peaks], np.arange(len(ref)), ref.wave)
    strength = np.array([p.snr for p in ref_peaks], dtype=float)
    labels = _labels_for(cat_wave, linelist_wave, linelist_labels)

    if wave_range is not None:
        w0, w1 = sorted(float(v) for v in wave_range)
        # A small margin absorbs an approximate wave_start/wave_end hint without
        # letting the far-red atlas lines (which are much brighter) back in.
        margin = max(50.0, 0.02 * (w1 - w0))
        m = (cat_wave >= w0 - margin) & (cat_wave <= w1 + margin)
        cat_wave, strength = cat_wave[m], strength[m]
        labels = [labels[i] for i in np.where(m)[0]]
    if len(cat_wave) == 0:
        raise ValueError('The reference catalogue is empty inside the requested range')

    threshold = max(min_strength, strength_frac * float(strength.max()))
    keep = strength >= threshold
    if keep.sum() > max_lines:
        order = np.argsort(-strength[keep])[:max_lines]
        sel = np.where(keep)[0][order]
    else:
        sel = np.where(keep)[0]
    sel = sel[np.argsort(cat_wave[sel])]
    cat = Reference(cat_wave[sel], strength[sel],
                    [labels[i] for i in sel], ref.source)
    logger.info('Reference catalogue: %d lines from %s (S/N >= %.0f)',
                len(cat), ref.source, threshold)
    return cat


# ---------------------------------------------------------------------------
# Dispersion solution
# ---------------------------------------------------------------------------
def initial_coefficients(n_pixels: int, wave_range: Sequence[float],
                         order: int = 1) -> np.ndarray:
    """Return polynomial coefficients mapping 0-based pixel to wavelength."""
    w0, w1 = sorted(float(v) for v in wave_range)
    if not (w1 > w0):
        raise ValueError('wave_range must be increasing')
    if order <= 1:
        return np.array([(w1 - w0) / max(n_pixels - 1, 1), w0], dtype=float)
    return np.polyfit([0, n_pixels - 1], [w0, w1], order)


def _match_peaks(pixel: np.ndarray, catalog: np.ndarray, coef: np.ndarray,
                 tol: float):
    """Pair peaks with the nearest catalogue line, at most one peak per line."""
    pixel = np.asarray(pixel, dtype=float)
    catalog = np.asarray(catalog, dtype=float)
    empty = (np.empty(0, dtype=float), np.empty(0, dtype=float), np.empty(0, dtype=int))
    if len(pixel) == 0 or len(catalog) == 0:
        return empty
    wave = np.polyval(coef, pixel)
    idx = np.clip(np.searchsorted(catalog, wave), 1, len(catalog) - 1)
    near = np.where(np.abs(catalog[idx - 1] - wave) < np.abs(catalog[idx] - wave),
                    idx - 1, idx)
    dist = np.abs(catalog[near] - wave)
    best = {}
    for i in np.where(dist < tol)[0]:
        j = int(near[i])
        if j not in best or dist[i] < dist[best[j]]:
            best[j] = int(i)
    if not best:
        return empty
    ci = np.array(sorted(best.keys()), dtype=int)
    pi = np.array([best[j] for j in ci], dtype=int)
    return pixel[pi], catalog[ci], pi


def _refine(pixel: np.ndarray, catalog: np.ndarray, coef: np.ndarray,
            order: int, tol0: float = 2.0, min_tol: float = 0.3,
            iterations: int = 12):
    """Iteratively match and fit, starting from a linear/quadratic guess."""
    coef = np.asarray(coef, dtype=float)
    px = np.empty(0)
    py = np.empty(0)
    pi = np.empty(0, dtype=int)
    for it in range(iterations):
        tol = max(tol0 * (0.75 ** it), min_tol)
        px, py, pi = _match_peaks(pixel, catalog, coef, tol)
        if len(px) < order + 2:
            break
        fit_order = min(order, 1) if it < 2 else order
        coef = np.polyfit(px, py, fit_order)
        if fit_order < order:
            coef = np.concatenate([np.zeros(order - fit_order), coef])
        # Robust clipping: blends and cosmic rays should not drag the fit.
        for _ in range(4):
            res = py - np.polyval(coef, px)
            mad = 1.4826 * np.median(np.abs(res - np.median(res)))
            keep = np.abs(res - np.median(res)) < max(3.0 * mad, 0.2)
            if keep.sum() < order + 2:
                break
            new = np.polyfit(px[keep], py[keep], order)
            if np.allclose(new, coef, rtol=0, atol=1e-12):
                coef = new
                break
            coef = new
    px, py, pi = _match_peaks(pixel, catalog, coef, max(min_tol * 1.5, 0.45))
    if len(px) < order + 2:
        return None
    res = py - np.polyval(coef, px)
    ci = np.array([int(np.argmin(np.abs(catalog - w))) for w in py], dtype=int)
    return coef, np.column_stack([px, py]), res, pi, ci


def _solution_score(coef, pixel, catalog, amp) -> Tuple[float, int, float]:
    """Score a solution by the number and quality of sub-arcsecond matches."""
    px, py, pi = _match_peaks(pixel, catalog, coef, 0.6)
    if len(px) == 0:
        return 0.0, 0, np.inf
    res = py - np.polyval(coef, px)
    weight = np.sqrt(np.clip(amp[pi] / max(amp.max(), 1e-30), 1e-6, 1.0))
    score = float(np.sum(weight)) + 0.25 * len(px)
    rms = float(np.sqrt(np.mean(res ** 2)))
    return score, len(px), rms


def solve_dispersion(peaks: Sequence[Peak], catalog: Reference, n_pixels: int,
                     wave_range: Optional[Sequence[float]] = None, order: int = 3,
                     match_tol: float = 2.0, min_tol: float = 0.3,
                     curvature_step: float = 2e-5, max_curvature: float = 2e-4,
                     n_peaks: int = 12, max_catalog: int = 80) -> Optional[Solution]:
    """Find the pixel-to-wavelength solution without manual anchors.

    A grating's dispersion is never exactly linear, so the line-spacing
    pattern match is run for a family of curvatures ``q = p + a*p^2``.  In
    ``q`` space the dispersion is nearly linear and :func:`pattern_match_linear`
    can find the correct alias; each seed is then refined with a polynomial in
    the original pixel coordinate.

    ``wave_range`` is only a *soft* hint here: it gives a small score bonus
    to solutions whose observed lines fall inside it.  The caller is expected
    to restrict the catalogue with :func:`build_catalog` when the approximate
    coverage is known; doing both keeps a wrong alias from winning without
    making the seed search fragile.
    """
    if len(peaks) < 4 or len(catalog) < 4:
        return None
    pixel = np.array([p.pixel for p in peaks], dtype=float)
    amp = np.array([p.amplitude for p in peaks], dtype=float)
    cat_wave = np.asarray(catalog.wave, dtype=float)

    alphas = np.arange(-max_curvature, max_curvature + 0.5 * curvature_step,
                       curvature_step)
    if not np.any(np.abs(alphas) < 1e-12):
        alphas = np.append(alphas, 0.0)

    candidates = []
    for alpha in alphas:
        q = pixel + alpha * pixel ** 2
        if np.any(np.diff(q) <= 0):
            continue
        try:
            seed = pattern_match_linear(
                q, cat_wave, amp, wave_range=(2500.0, 12500.0),
                pixel_range=(0.0, float(q.max())),
                n_peaks=min(n_peaks, len(pixel)),
                max_catalog=min(len(cat_wave), max_catalog),
                ratio_eps=0.02)
        except Exception:  # pragma: no cover - defensive
            seed = None
        if seed is None:
            continue
        scale, offset = float(seed[0]), float(seed[1])
        if not np.isfinite(scale) or scale <= 0:
            continue
        # lambda = scale*q + offset = (scale*a)*p^2 + scale*p + offset
        coef0 = np.array([scale * alpha, scale, offset], dtype=float)
        out = _refine(pixel, cat_wave, coef0, order, match_tol, min_tol)
        if out is None:
            continue
        # The shrinking tolerance can drop faint lines at the detector ends,
        # which then makes the polynomial extrapolate badly.  A second pass
        # starting from the first solution recovers them without loosening
        # the early iterations (where a loose tolerance would accept blends).
        out2 = _refine(pixel, cat_wave, out[0], order,
                       max(match_tol, 1.5), max(min_tol, 0.7))
        if out2 is not None and len(out2[1]) > len(out[1]):
            _, _, rms2 = _solution_score(out2[0], pixel, cat_wave, amp)
            if rms2 < 1.0:
                out = out2
        coef, pairs, res, pi, ci = out
        # A physical dispersion must be monotonic over the detector.
        probe = np.linspace(0, n_pixels - 1, 32)
        if np.any(np.diff(np.polyval(coef, probe)) <= 0):
            continue
        score, ngood, rms = _solution_score(coef, pixel, cat_wave, amp)
        if ngood < max(order + 2, 8):
            continue
        # Reward coverage: aliases typically match only 10-20 peaks, while the
        # true solution matches most of the strong lamp lines.
        score += 0.5 * ngood
        if wave_range is not None:
            w0, w1 = sorted(float(v) for v in wave_range)
            obs = np.polyval(coef, pixel)
            lo, hi = float(obs.min()), float(obs.max())
            overlap = max(0.0, min(w1, hi) - max(w0, lo))
            score += 2.0 * overlap / max(hi - lo, 1.0)
        candidates.append((score, ngood, -rms, coef, pairs, res, pi, ci))

    if not candidates:
        return None
    candidates.sort(key=lambda t: (t[0], t[1], t[2]), reverse=True)
    for k, cand in enumerate(candidates[:5]):
        coef = cand[3]
        logger.info('  candidate %d: %d matches, rms=%.3f A, range %.1f-%.1f A',
                    k + 1, cand[1], -cand[2], np.polyval(coef, 0),
                    np.polyval(coef, n_pixels - 1))
    score, ngood, neg_rms, coef, pairs, res, pi, ci = candidates[0]
    # Final recovery pass: the winning fit may have dropped faint lines at the
    # detector ends.  Add them back if they are consistent with the solution.
    out = _refine(pixel, cat_wave, coef, order,
                  max(match_tol, 1.5), max(min_tol, 0.7))
    if out is not None and len(out[1]) > len(pairs):
        coef2, pairs2, res2, pi2, ci2 = out
        if float(np.sqrt(np.mean(res2 ** 2))) < 0.6:
            coef, pairs, res, pi, ci = coef2, pairs2, res2, pi2, ci2
    logger.info('Dispersion solution: %d matches, rms=%.3f A, range %.1f-%.1f A',
                len(pairs), float(np.sqrt(np.mean(res ** 2))),
                np.polyval(coef, 0), np.polyval(coef, n_pixels - 1))
    return Solution(coef, pairs, res, pi, ci, n_pixels)


# ---------------------------------------------------------------------------
# Anchor selection and output
# ---------------------------------------------------------------------------
def select_anchors(peaks: Sequence[Peak], solution: Solution, ref: Reference,
                   min_snr: float = 8.0, max_residual: float = 0.6,
                   min_separation: float = 4.0, min_catalog_separation: float = 1.5,
                   max_anchors: int = 30, min_wavelength_gap: float = 8.0,
                   edge_min_snr: float = 5.0, resolution: float = 0.0):
    """Choose a clean, well-spread subset of matched lines as anchors.

    Strong, isolated lines are selected first.  The bluest and reddest
    candidates are added when they extend the wavelength coverage, because
    they constrain the slope much better than a purely S/N-ranked selection.
    """
    if solution is None or len(solution.pairs) == 0:
        return []
    pixel_all = np.array([p.pixel for p in peaks], dtype=float)
    cand = []
    for k, (x, w) in enumerate(solution.pairs):
        pi = int(solution.peak_indices[k])
        peak = peaks[pi]
        res = float(solution.residuals[k])
        if peak.snr < edge_min_snr or abs(res) > max_residual:
            continue
        others = np.abs(pixel_all - peak.pixel)
        others[pi] = np.inf
        if len(others) and others.min() < min_separation:
            continue
        j = int(solution.catalog_indices[k])
        label = ref.labels[j] if 0 <= j < len(ref.labels) else ''
        if 0 <= j < len(ref.wave):
            sep_min = float(min_catalog_separation)
            if resolution and resolution > 0:
                # Two reference lines inside one resolution element are a
                # blend in the observed spectrum and make poor anchors.
                sep_min = max(sep_min, 1.2 * float(w) / float(resolution))
            sep = np.abs(ref.wave - ref.wave[j])
            sep[j] = np.inf
            if len(sep) and sep.min() < sep_min:
                continue
        cand.append((peak.snr, float(x), float(w), res, label))
    if not cand:
        return []

    # Seed the selection with the bluest and reddest useful lines: they are
    # what constrains the slope, and a purely S/N-ranked list often stops
    # before the detector edges.  Then fill up with the strongest lines.
    by_wave = sorted(cand, key=lambda t: t[2])
    selected = []
    for edge in (by_wave[0], by_wave[-1]):
        if all(abs(edge[2] - s[2]) >= min_wavelength_gap for s in selected):
            selected.append(edge)

    strong = sorted((c for c in cand if c[0] >= min_snr), key=lambda t: -t[0])
    for item in strong:
        if len(selected) >= max_anchors:
            break
        if all(abs(item[2] - s[2]) >= min_wavelength_gap for s in selected):
            selected.append(item)
    selected.sort(key=lambda t: t[1])
    return selected


def _format_solution(coef: np.ndarray, precision: int = 8) -> Tuple[str, str]:
    degree = len(coef) - 1
    terms = []
    for i, c in enumerate(coef):
        power = degree - i
        if power == 0:
            terms.append(f'{c:.{precision}g}')
        elif power == 1:
            terms.append(f'{c:.{precision}g}*x')
        else:
            terms.append(f'{c:.{precision}g}*x^{power}')
    return ' + '.join(terms), '[' + ', '.join(f'{c:.{precision}g}' for c in coef) + ']'


def write_anchor_file(path: str, anchors, solution: Solution, lamp: str,
                      reference: str, note: str = '', created: str = ''):
    """Write anchors as ``pixel, wavelength`` CSV readable by the slit pipeline."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    formula, coef_text = _format_solution(solution.coef)
    lines = [
        '# SpecProc anchor lines for long-slit wavelength calibration',
        '# pixel (0-based along the dispersion axis), wavelength (Angstrom, air)',
        '# These are a small set of clean, isolated strong lines.  They seed the',
        '# full-catalogue matching; they are not the complete lamp line list.',
        f'# lamp       : {lamp}',
        f'# reference  : {reference}',
    ]
    if created:
        lines.append(f'# created    : {created}')
    lines += [
        f'# solution   : lambda(x) = {formula}',
        f'# coef       : {coef_text}  (np.polyval order)',
        f'# solution_rms: {solution.rms:.3f} A',
        f'# anchors    : {len(anchors)}',
    ]
    if note:
        lines.append(f'# note       : {note}')
    lines.append('# pixel, wavelength, ion, snr, residual')
    for snr, x, w, res, label in anchors:
        lines.append(f'{x:.3f}, {w:.4f}, {label or "?"}, {snr:.1f}, {res:+.2f}')
    p.write_text('\n'.join(lines) + '\n')
    return p


# ---------------------------------------------------------------------------
# Diagnostic figure
# ---------------------------------------------------------------------------
def _convolve_reference(ref: Reference, wave_grid: np.ndarray, fwhm: float) -> np.ndarray:
    """Resample the reference spectrum and smooth it to the observed resolution."""
    flux = np.interp(wave_grid, ref.wave, np.clip(ref.flux, 0, None), left=0, right=0)
    step = np.median(np.diff(wave_grid)) if len(wave_grid) > 1 else 0.0
    if fwhm > 0 and step > 0:
        flux = ndimage.gaussian_filter1d(flux, max(fwhm / 2.355 / step, 0.3))
    return flux


def _pick_zoom(anchors, wavelength: np.ndarray, width: float = 420.0) -> Tuple[float, float]:
    """Pick the window that contains the largest number of anchors."""
    if not anchors:
        return float(wavelength[0]), float(wavelength[-1])
    waves = np.array([a[2] for a in anchors])
    best = max(waves, key=lambda w: np.sum(np.abs(waves - w) < width / 2))
    lo, hi = best - width / 2, best + width / 2
    lo = max(lo, float(wavelength[0]))
    hi = min(hi, float(wavelength[-1]))
    if hi - lo < width:
        if lo <= wavelength[0]:
            hi = min(lo + width, float(wavelength[-1]))
        else:
            lo = max(hi - width, float(wavelength[0]))
    return float(lo), float(hi)


def plot_diagnostic(path: str, flux: np.ndarray, peaks: Sequence[Peak],
                    solution: Solution, anchors, ref: Optional[Reference] = None,
                    title: str = '', zoom: Optional[Tuple[float, float]] = None,
                    pixel_wavelength: Optional[np.ndarray] = None):
    """Draw the visual-check figure: spectrum, zoomed overlay and residuals."""
    import os
    import tempfile
    os.environ.setdefault('MPLCONFIGDIR', tempfile.mkdtemp(prefix='specproc-mpl-'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    pixel = np.arange(len(flux), dtype=float)
    if pixel_wavelength is None:
        pixel_wavelength = solution.wavelength(pixel)
    pixel_wavelength = np.asarray(pixel_wavelength, dtype=float)

    fig = plt.figure(figsize=(14, 10.5))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.35, 1.05, 0.85], hspace=0.30)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    ax3 = fig.add_subplot(gs[2])

    # --- full spectrum -----------------------------------------------------
    ax1.plot(pixel, np.clip(flux, 1, None), color='0.25', lw=0.55, zorder=2)
    ax1.set_yscale('log')
    ax1.set_xlim(0, len(flux) - 1)
    ax1.set_ylabel('Counts (log)')
    ax1.set_xlabel('Pixel (0-based)')
    ax1.set_title(title or 'Lamp spectrum and anchor lines', pad=10)
    ax1.grid(alpha=0.18)

    matched = solution.pairs[:, 0] if solution is not None and len(solution.pairs) else []
    for x in matched:
        ax1.axvline(x, color='0.60', lw=0.45, alpha=0.65, zorder=1)
    for snr, x, w, res, label in anchors:
        ax1.axvline(x, color='crimson', lw=0.9, ls='--', alpha=0.9, zorder=3)
    # Alternate two label rows so close pairs stay readable.
    for k, (snr, x, w, res, label) in enumerate(anchors):
        ypos = 0.985 if k % 2 == 0 else 0.905
        ax1.text(x, ypos, f'{w:.0f} {label or ""}', rotation=90, ha='center', va='top',
                 transform=ax1.get_xaxis_transform(), fontsize=6.2, color='crimson',
                 bbox=dict(fc='white', ec='none', alpha=0.80, pad=0.30))

    formula, _ = _format_solution(solution.coef, precision=6)
    w_first, w_last = solution.wave_range
    if anchors:
        a0 = min(a[2] for a in anchors)
        a1 = max(a[2] for a in anchors)
        coverage = f'anchors {a0:.0f}-{a1:.0f} Å'
    else:
        coverage = 'no anchors'
    ax1.text(0.995, 0.015,
             f'λ(x) = {formula}\nrange {w_first:.0f}-{w_last:.0f} Å'
             f'  ({coverage})\n{len(anchors)} anchors  |  rms {solution.rms:.2f} Å',
             transform=ax1.transAxes, ha='right', va='bottom', fontsize=7.3,
             color='0.20', bbox=dict(fc='white', ec='0.7', alpha=0.88, pad=2.5))

    # --- zoom on the densest anchor region ---------------------------------
    if zoom is None:
        zoom = _pick_zoom(anchors, pixel_wavelength)
    if zoom is not None:
        m = ((pixel_wavelength >= zoom[0]) & (pixel_wavelength <= zoom[1]))
        if m.sum() > 10:
            ax2b = ax2.twinx()
            ax2.plot(pixel[m], flux[m], color='black', lw=0.9, label='observed', zorder=3)
            if ref is not None and len(ref) and not np.allclose(ref.flux, 1.0):
                wgrid = pixel_wavelength[m]
                fwhm_px = np.median([p.width for p in peaks]) if peaks else 3.0
                mid = 0.5 * (float(pixel[m].min()) + float(pixel[m].max()))
                disp = abs(float(np.polyval(np.polyder(solution.coef), mid)))
                fwhm = max(fwhm_px * disp, 1.0)
                rflux = _convolve_reference(ref, wgrid, fwhm)
                # Match the reference to the observed window scale; the atlas
                # has no absolute flux calibration between echelle orders.
                scale = (np.nanpercentile(flux[m], 90) /
                         max(np.nanpercentile(rflux, 90), 1e-30))
                ax2b.plot(pixel[m], rflux * scale, color='tab:orange', lw=0.9,
                          alpha=0.75, label='reference', zorder=2)
                ax2b.set_ylabel('Reference (scaled)', color='tab:orange')
                ax2b.tick_params(axis='y', colors='tab:orange')
            for snr, x, w, res, label in anchors:
                if zoom[0] <= w <= zoom[1]:
                    ax2.axvline(x, color='crimson', lw=0.8, ls='--', alpha=0.75, zorder=1)
                    ax2.text(x, 0.97, f'{w:.1f}\n{label or ""}', rotation=90, ha='right',
                             va='top', transform=ax2.get_xaxis_transform(), fontsize=6.3,
                             color='crimson',
                             bbox=dict(fc='white', ec='none', alpha=0.72, pad=0.30))
            ax2.set_xlim(pixel[m].min(), pixel[m].max())
            ax2.set_ylabel('Observed counts')
            ax2.set_xlabel('Pixel (0-based)')
            ax2.grid(alpha=0.18)
            ax2.set_title(
                f'Zoom: {zoom[0]:.0f}-{zoom[1]:.0f} Å  '
                f'(pixel {pixel[m].min():.0f}-{pixel[m].max():.0f}; '
                'black: observed, orange: atlas)')
            h1, l1 = ax2.get_legend_handles_labels()
            h2, l2 = ax2b.get_legend_handles_labels()
            ax2.legend(h1 + h2, l1 + l2, loc='upper left', fontsize=8)

    # --- residuals ---------------------------------------------------------
    if solution is not None and len(solution.pairs):
        res = solution.residuals
        ax3.axhline(0, color='0.3', lw=0.8)
        ax3.axhline(solution.rms, color='tab:blue', lw=0.7, ls=':')
        ax3.axhline(-solution.rms, color='tab:blue', lw=0.7, ls=':')
        ax3.plot(solution.pairs[:, 1], res, 'o', ms=3.5, color='0.55',
                 label='all matches', zorder=2)
        if anchors:
            ax3.plot([a[2] for a in anchors], [a[3] for a in anchors], 'o', ms=5.5,
                     color='crimson', label='anchors', zorder=3)
        ax3.set_xlabel('Wavelength (Å)')
        ax3.set_ylabel('Residual (Å)')
        ax3.grid(alpha=0.18)
        ax3.legend(loc='best', fontsize=8)
        ax3.set_title(f'Residuals: rms = {solution.rms:.2f} Å, '
                      f'{len(anchors)} anchors')

    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=160, bbox_inches='tight')
    plt.close(fig)
    return out
