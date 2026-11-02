"""Sub-pixel centres for emission-line peaks.

The detection routines return integer peak samples from
``scipy.signal.find_peaks``.  This module refines each sample with a local
three-parameter Gaussian (plus a constant local baseline), the standard
arc-lamp centring method used by IRAF ``center1d`` and PypeIt's
``fit_arcspec``.  If the fit fails, a weighted centroid over the same window
is used.
"""

import numpy as np
from scipy.optimize import curve_fit


def _gaussian(x, amplitude, centre, sigma, baseline):
    return baseline + amplitude * np.exp(-0.5 * ((x - centre) / sigma) ** 2)


def _centroid(x, y, baseline):
    weights = np.clip(np.asarray(y, dtype=float) - baseline, 0.0, None)
    total = float(weights.sum())
    if total <= 0 or not np.isfinite(total):
        return None
    return float(np.sum(np.asarray(x, dtype=float) * weights) / total)


def fit_peak_center(x, y, index: int, half_width: int = 3,
                    min_sigma: float = 0.05,
                    left_width: int = None, right_width: int = None) -> float:
    """Return a sub-pixel emission-line centre near ``index``.

    ``x`` is the coordinate of the samples (pixels or wavelengths), ``y`` the
    continuum-subtracted spectrum and ``index`` the sample nearest the peak.
    A four-parameter Gaussian (amplitude, centre, sigma, constant baseline) is
    fitted in a symmetric window; failure falls back to a weighted centroid.
    """
    x = np.asarray(x, dtype=float).reshape(-1)
    y = np.asarray(y, dtype=float).reshape(-1)
    if len(x) != len(y) or len(x) == 0:
        return float(index)
    index = int(np.clip(index, 0, len(x) - 1))
    half_width = max(1, int(half_width))
    left_width = half_width if left_width is None else max(1, int(left_width))
    right_width = half_width if right_width is None else max(1, int(right_width))
    lo = max(0, index - left_width)
    hi = min(len(x), index + right_width + 1)
    if hi - lo < 3:
        return float(x[index])

    xf = x[lo:hi]
    yf = y[lo:hi]
    finite = np.isfinite(xf) & np.isfinite(yf)
    if finite.sum() < 3:
        return float(x[index])
    xf, yf = xf[finite], yf[finite]

    baseline = float(np.min(yf))
    amplitude = float(np.max(yf) - baseline)
    if not np.isfinite(amplitude) or amplitude <= 0:
        centre = _centroid(xf, yf, baseline)
        return float(x[index]) if centre is None else centre

    spacing = np.median(np.diff(xf)) if len(xf) > 1 else 1.0
    sigma0 = max(abs(float(spacing)) * 1.5, 1e-6)
    p0 = [amplitude, float(x[index]), sigma0, baseline]
    lower = [0.0, float(np.min(xf)), min_sigma, -np.inf]
    upper = [np.inf, float(np.max(xf)), np.inf, np.inf]
    try:
        popt, _ = curve_fit(_gaussian, xf, yf, p0=p0, bounds=(lower, upper),
                            maxfev=3000)
        centre = float(popt[1])
        if np.isfinite(centre) and float(np.min(xf)) <= centre <= float(np.max(xf)):
            return centre
    except (RuntimeError, ValueError, TypeError):
        pass

    centre = _centroid(xf, yf, baseline)
    return float(x[index]) if centre is None else centre


def fit_peak_shape(x, y, index: int, half_width: int = 3,
                   min_sigma: float = 0.05,
                   left_width: int = None, right_width: int = None):
    """Return ``(centre, sigma)`` from a local Gaussian fit, or ``None``.

    The fit window has the same neighbour-aware bounds as
    :func:`fit_peak_center`; ``sigma`` is in the units of ``x``.
    """
    x = np.asarray(x, dtype=float).reshape(-1)
    y = np.asarray(y, dtype=float).reshape(-1)
    if len(x) != len(y) or len(x) == 0:
        return None
    index = int(np.clip(index, 0, len(x) - 1))
    half_width = max(1, int(half_width))
    left_width = half_width if left_width is None else max(1, int(left_width))
    right_width = half_width if right_width is None else max(1, int(right_width))
    lo = max(0, index - left_width)
    hi = min(len(x), index + right_width + 1)
    if hi - lo < 4:
        return None
    xf = x[lo:hi]
    yf = y[lo:hi]
    finite = np.isfinite(xf) & np.isfinite(yf)
    if finite.sum() < 4:
        return None
    xf, yf = xf[finite], yf[finite]
    baseline = float(np.min(yf))
    amplitude = float(np.max(yf) - baseline)
    if not np.isfinite(amplitude) or amplitude <= 0:
        return None
    spacing = np.median(np.diff(xf)) if len(xf) > 1 else 1.0
    sigma0 = max(abs(float(spacing)) * 1.5, 1e-6)
    p0 = [amplitude, float(x[index]), sigma0, baseline]
    lower = [0.0, float(np.min(xf)), min_sigma, -np.inf]
    upper = [np.inf, float(np.max(xf)), np.inf, np.inf]
    try:
        popt, _ = curve_fit(_gaussian, xf, yf, p0=p0, bounds=(lower, upper),
                            maxfev=3000)
    except (RuntimeError, ValueError, TypeError):
        return None
    centre, sigma = float(popt[1]), float(popt[2])
    if not (np.isfinite(centre) and np.isfinite(sigma) and sigma > 0):
        return None
    if not (float(np.min(xf)) <= centre <= float(np.max(xf))):
        return None
    return centre, sigma
