"""Tests for the automatic lamp-line anchor matcher."""

from pathlib import Path

import numpy as np
from astropy.io import fits

import specproc
from specproc.core.line_identification import load_anchors, pattern_match_linear
from specproc.core.anchor_matching import (
    Reference,
    _match_peaks,
    detect_peaks,
    load_observed,
    load_reference,
    select_anchors,
    solve_dispersion,
    write_anchor_file,
)


def _synthetic_lamp():
    """A lamp with a known cubic dispersion and 20 lines."""
    n = 2048
    true_coef = np.array([-2.0e-8, 8.0e-5, 2.1, 5000.0])
    x = np.arange(n, dtype=float)
    wave = np.polyval(true_coef, x)
    catalog = np.array([
        5200.0, 5380.0, 5610.0, 5840.0, 6030.0, 6250.0, 6480.0, 6710.0,
        6930.0, 7140.0, 7360.0, 7580.0, 7790.0, 8010.0, 8240.0, 8460.0,
        8680.0, 8890.0, 9110.0, 9320.0,
    ])
    peaks_px = np.interp(catalog, wave, x)  # wave is monotonic
    rng = np.random.default_rng(7)
    flux = np.full(n, 50.0)
    for p in peaks_px:
        flux += 1200.0 * np.exp(-0.5 * ((x - p) / 1.35) ** 2)
    flux += rng.normal(0.0, 2.0, n)
    return flux, catalog, true_coef


def test_solve_dispersion_recovers_known_solution():
    flux, catalog, true_coef = _synthetic_lamp()
    peaks = detect_peaks(flux, snr=8.0, min_width=0.5, max_width=8.0,
                         distance=5, baseline_size=101)
    assert len(peaks) >= 18

    solution = solve_dispersion(
        peaks, Reference(catalog), len(flux), wave_range=(4800.0, 9600.0),
        order=3, match_tol=3.0, min_tol=0.5)
    assert solution is not None
    assert solution.n_lines >= 15
    assert solution.rms < 0.5
    true_w0 = float(np.polyval(true_coef, 0.0))
    true_w1 = float(np.polyval(true_coef, len(flux) - 1))
    assert abs(solution.wave_range[0] - true_w0) < 2.0
    assert abs(solution.wave_range[1] - true_w1) < 2.0

    anchors = select_anchors(
        peaks, solution, Reference(catalog), min_snr=8.0, max_residual=0.4,
        min_separation=5.0, max_anchors=30, resolution=2000.0)
    assert len(anchors) >= 15
    assert anchors == sorted(anchors, key=lambda a: a[1])


def test_match_peaks_assigns_each_peak_once():
    # Two catalogue lines fall inside the tolerance of the first peak; without
    # the one-peak-per-line rule that peak would appear twice in the fit.
    pixel = np.array([10.0, 20.0])
    catalog = np.array([100.0, 100.5, 200.0])
    coef = np.array([10.0, 0.0])
    px, wl, pi = _match_peaks(pixel, catalog, coef, tol=1.0)
    assert len(px) == 2
    assert len(np.unique(px)) == 2
    assert list(wl) == [100.0, 200.0]


def test_load_observed_specproc_table(tmp_path):
    n = 64
    wave = np.arange(1, n + 1, dtype=float)
    flux_sum = np.linspace(100.0, 200.0, n)
    flux_opt = np.linspace(90.0, 190.0, n)
    cols = [
        fits.Column(name='WAVELENGTH', format='D', array=wave),
        fits.Column(name='FLUX_OPT', format='D', array=flux_opt),
        fits.Column(name='FLUX_SUM', format='D', array=flux_sum),
    ]
    hdu = fits.BinTableHDU.from_columns(cols)
    hdu.header['WAVEUNIT'] = 'pixel'
    path = tmp_path / 'lamp_spec.fits'
    fits.HDUList([fits.PrimaryHDU(), hdu]).writeto(path)

    flux, pixel, info = load_observed(str(path))
    assert len(flux) == n
    assert pixel[0] == 0.0 and pixel[-1] == n - 1
    assert np.allclose(flux, flux_sum)
    assert 'wavelength' not in info


def test_write_anchor_file(tmp_path):
    flux, catalog, _ = _synthetic_lamp()
    peaks = detect_peaks(flux, snr=8.0, min_width=0.5, max_width=8.0,
                         distance=5, baseline_size=101)
    solution = solve_dispersion(
        peaks, Reference(catalog), len(flux), wave_range=(4800.0, 9600.0),
        order=3, match_tol=3.0, min_tol=0.5)
    anchors = select_anchors(
        peaks, solution, Reference(catalog), min_snr=8.0, max_residual=0.4,
        min_separation=5.0, max_anchors=30, resolution=2000.0)

    out = tmp_path / 'anchors.csv'
    write_anchor_file(out, anchors, solution, 'lamp.fits', 'reference.fits',
                      note='test', created='2026-01-01')
    text = out.read_text()
    assert '# anchors' in text
    assert '# pixel, wavelength, ion, snr, residual' in text
    rows = [ln for ln in text.splitlines() if ln and not ln.startswith('#')]
    assert len(rows) == len(anchors)
    assert len(rows[0].split(',')) == 5


def test_load_reference_reads_fear_atlas_wcs():
    """fear.fits has a linear WCS: 3003.46 A + 0.2 A/pixel, 39709 pixels."""
    path = Path(specproc.__file__).resolve().parent / 'calib_data' / 'linelists' / 'fear.fits'
    assert path.exists()
    ref = load_reference(str(path))
    assert len(ref) == 39709
    assert abs(ref.wave[0] - 3003.4624) < 1e-3
    assert abs(ref.wave[-1] - (3003.4624 + 0.2 * 39708)) < 1e-2


def test_pattern_match_recovers_linear_dispersion():
    """The batched scoring in pattern_match_linear keeps the old behaviour."""
    catalog = np.linspace(5000.0, 7000.0, 40)
    scale, offset = 0.5, 4800.0
    rng = np.random.default_rng(11)
    pixel = (catalog - offset) / scale + rng.normal(0.0, 0.15, len(catalog))
    coef = pattern_match_linear(pixel, catalog, np.ones_like(pixel))
    assert coef is not None
    assert abs(coef[0] - scale) < 0.01
    assert abs(coef[1] - offset) < 5.0


def test_anchor_file_round_trip(tmp_path):
    flux, catalog, _ = _synthetic_lamp()
    peaks = detect_peaks(flux, snr=8.0, min_width=0.5, max_width=8.0,
                         distance=5, baseline_size=101)
    solution = solve_dispersion(
        peaks, Reference(catalog), len(flux), wave_range=(4800.0, 9600.0),
        order=3, match_tol=3.0, min_tol=0.5)
    anchors = select_anchors(
        peaks, solution, Reference(catalog), min_snr=8.0, max_residual=0.4,
        min_separation=5.0, max_anchors=30, resolution=2000.0)
    out = tmp_path / 'anchors.csv'
    write_anchor_file(out, anchors, solution, 'lamp.fits', 'reference.fits')

    loaded = load_anchors(str(out))
    assert loaded.shape == (len(anchors), 2)
    assert np.allclose(loaded[:, 0], sorted(a[1] for a in anchors))
    assert np.allclose(loaded[:, 1], [a[2] for a in sorted(anchors, key=lambda a: a[1])])
