"""
Helpers for identifying lamp lines in long-slit spectra.

Shared by the interactive identification window (``specproc.gui.line_identifier``)
and the slit pipeline. An *anchor* is a manually identified (pixel, wavelength)
pair; anchors seed the dispersion solution, which is then extended automatically
by matching all detected lamp peaks against the full line catalogue.
"""

import itertools
import logging
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
from astropy.io import fits
from scipy import ndimage

logger = logging.getLogger(__name__)


def find_arc_peaks(spec: np.ndarray, snr_threshold: float = 8.0, with_snr: bool = False):
    """Detect emission peaks above snr_threshold (peak height / background noise).

    Returns (centroid pixel (0-based), amplitude) tuples, or
    (pixel, amplitude, snr) tuples when with_snr is True.
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
    pk = [i for i in range(2, len(s) - 2)
          if s[i] > s[i - 1] and s[i] >= s[i + 1] and s[i] > snr_threshold * noise_arr[i] and s[i] > s[i - 2]]
    out = []
    for i in pk:
        lo, hi = max(i - 3, 0), min(i + 4, len(s))
        w = np.clip(s[lo:hi], 0, None)
        entry = (float((np.arange(lo, hi) * w).sum() / w.sum()), float(s[i]))
        out.append(entry + (float(s[i] / noise_arr[i]),) if with_snr else entry)
    return out


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

    Accepts a SpecProc ``*_spec.fits`` table (FLUX_SUM/FLUX_OPT), a 1-D or 2-D
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


def _match_nearest(pix: np.ndarray, catalog: np.ndarray, coef: np.ndarray, tol: float) -> np.ndarray:
    """Pair each peak with its nearest catalogue line within tol (each line used once)."""
    wl = np.polyval(coef, pix)
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

    best, best_key = None, (0, 0.0)
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
            wl = scale[:, None] * pix[None, :] + offset[:, None]
            idx = np.clip(np.searchsorted(cat, wl), 1, len(cat) - 1)
            dist = np.minimum(np.abs(cat[idx] - wl), np.abs(cat[idx - 1] - wl))
            tol_w = np.maximum(1.5 * np.abs(scale), 1e-3)[:, None]
            hit = dist < tol_w
            score = hit.sum(axis=1)
            rms = np.sqrt((np.where(hit, dist, 0) ** 2).sum(axis=1) / np.maximum(score, 1))
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
                        peak_amp=None, initial_coef=None):
    """Match detected lamp peaks to a catalogue, optionally without manual anchors.

    With anchors, they give the rough solution and are always retained. With no
    anchors, the starting solution comes from line-spacing pattern matching
    (``pattern_match_linear``). Matching then uses a shrinking tolerance while the
    polynomial order grows. ``initial_coef`` (polynomial coefficients, pixel ->
    wavelength) replaces the starting solution, e.g. a hand-aligned overlay.
    Returns (pairs, coef).
    """
    anchors = np.asarray(anchors, dtype=float).reshape(-1, 2)
    pix = np.asarray(peak_pix, dtype=float)
    catalog = np.asarray(catalog, dtype=float)
    catalog = np.sort(catalog[np.isfinite(catalog)])
    if len(catalog) < 2 or len(pix) == 0:
        return anchors, fit_polynomial(anchors, order)

    if len(anchors) == 0:
        coef = (np.asarray(initial_coef, dtype=float) if initial_coef is not None else
                pattern_match_linear(pix, catalog, peak_amp, initial_range, initial_pixel_range))
        if coef is None:
            return np.empty((0, 2)), None
        pairs = np.empty((0, 2))
        for it in range(max_iter):
            cur_tol = max(tol / (it + 1), min_tol)
            matched = _match_nearest(pix, catalog, coef, cur_tol)
            if len(matched) < 2:
                break
            pairs = matched
            coef = fit_polynomial(pairs, min(order, it + 1))
        if len(pairs) < 2:
            return np.empty((0, 2)), None
        return pairs[np.argsort(pairs[:, 0])], fit_polynomial(pairs, order)

    coef = (np.asarray(initial_coef, dtype=float) if initial_coef is not None
            else fit_polynomial(anchors, 1))
    if coef is None:
        return anchors, fit_polynomial(anchors, order)

    pairs = anchors
    for it in range(max_iter):
        cur_tol = max(tol / (it + 1), min_tol)
        matched = _match_nearest(pix, catalog, coef, cur_tol)
        if len(matched):
            near_anchor = np.min(np.abs(matched[:, :1] - anchors[:, 0][None, :]), axis=1) < 1.5
            same_line = np.isin(np.round(matched[:, 1], 4), np.round(anchors[:, 1], 4))
            matched = matched[~(near_anchor | same_line)]
        candidate = np.vstack([anchors, matched])
        if len(candidate) < 2:
            break
        pairs = candidate
        coef = fit_polynomial(pairs, min(order, it + 1))
    coef = fit_polynomial(pairs, order)
    return pairs[np.argsort(pairs[:, 0])], coef
