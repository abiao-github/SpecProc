"""Tests for interactive/automatic identification and the anchor exporter."""

import os
from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits

from specproc.core.anchor_matching import (
    Reference,
    build_catalog,
    detect_peaks,
    load_observed,
    load_reference,
    solution_from_pairs,
    solve_from_anchors,
)
from specproc.core.line_identification import (
    load_lamp_spectrum,
    dispersion_polyval_coef,
    dispersion_type,
    fit_dispersion,
    measure_resolution,
    normalize_poly_type,
)
from specproc.utils.identify import _fits_header_meta, build_parser, main, solve_auto

TRUE_COEF = np.array([-2.0e-8, 8.0e-5, 2.1, 5000.0])
CATALOG = np.array([
    5200.0, 5380.0, 5610.0, 5840.0, 6030.0, 6250.0, 6480.0, 6710.0,
    6930.0, 7140.0, 7360.0, 7580.0, 7790.0, 8010.0, 8240.0, 8460.0,
    8680.0, 8890.0, 9110.0, 9320.0,
])


def _synthetic_lamp(seed: int = 7, reverse: bool = False):
    n = 2048
    x = np.arange(n, dtype=float)
    wave = np.polyval(TRUE_COEF, n - 1 - x) if reverse else np.polyval(TRUE_COEF, x)
    order = np.argsort(wave)
    peaks_px = np.interp(CATALOG, wave[order], x[order])
    rng = np.random.default_rng(seed)
    flux = np.full(n, 50.0)
    for p in peaks_px:
        flux += 1200.0 * np.exp(-0.5 * ((x - p) / 1.35) ** 2)
    flux += rng.normal(0.0, 2.0, n)
    return flux


@pytest.fixture()
def lamp_files(tmp_path):
    flux = _synthetic_lamp()
    lamp = tmp_path / 'lamp.txt'
    ref = tmp_path / 'ref.txt'
    np.savetxt(lamp, flux)
    np.savetxt(ref, CATALOG)
    return lamp, ref


@pytest.fixture()
def reverse_lamp_files(tmp_path):
    flux = _synthetic_lamp(reverse=True)
    lamp = tmp_path / 'lamp_reverse.txt'
    ref = tmp_path / 'ref.txt'
    np.savetxt(lamp, flux)
    np.savetxt(ref, CATALOG)
    return lamp, ref


def test_csv_and_txt_lamp_inputs(tmp_path):
    flux = np.linspace(1.0, 20.0, 12)
    csv = tmp_path / 'lamp.csv'
    csv.write_text('pixel,flux\n' + '\n'.join(
        f'{i},{value:.6f}' for i, value in enumerate(flux)))
    txt = tmp_path / 'lamp.txt'
    np.savetxt(txt, flux)
    single = tmp_path / 'single_row.csv'
    single.write_text(','.join(f'{value:.6f}' for value in flux))

    for path in (csv, txt, single):
        loaded = load_lamp_spectrum(str(path))
        assert np.allclose(loaded, flux)


def test_mstg_telescope_and_arm_from_fits_header(tmp_path):
    for arm, instrument in (('R', 'Camera Spectrum-R'),
                            ('B', 'Camera Spectrum-B')):
        path = tmp_path / f'lamp_{arm}.fits'
        hdu = fits.PrimaryHDU(np.ones(16))
        hdu.header['TELESCOP'] = 'MOST    '
        hdu.header['INSTRUME'] = instrument
        hdu.writeto(path)
        telescope, detected = _fits_header_meta(str(path))
        assert telescope == 'mstg193slit'
        assert detected == arm


def test_gui_uses_header_arm_for_default_name(tmp_path):
    QtWidgets = pytest.importorskip('PyQt5.QtWidgets')
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    from specproc.utils.identify import LineIdentifier

    lamp = tmp_path / 'lamp_B.fits'
    hdu = fits.PrimaryHDU(np.full(64, 10.0))
    hdu.header['TELESCOP'] = 'MOST    '
    hdu.header['INSTRUME'] = 'Camera Spectrum-B'
    hdu.writeto(lamp)
    ref = tmp_path / 'ref.txt'
    ref.write_text('5200\n5380\n')
    window = LineIdentifier(str(lamp), reference_path=str(ref))
    try:
        assert window.telescope_name == 'mstg193slit'
        assert window.arm_combo.currentText() == 'B'
        assert window.default_anchor_path().endswith('mstg193slit_anchors_B.csv')
        assert Path(window.default_anchor_path()).parent == Path.cwd()
    finally:
        window.close()
        app.processEvents()


def test_gui_auto_identify_shows_match_connectors_and_unobtrusive_legend(tmp_path):
    QtWidgets = pytest.importorskip('PyQt5.QtWidgets')
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    from specproc.utils.identify import LineIdentifier

    lamp = tmp_path / 'lamp.txt'
    reference = tmp_path / 'reference.txt'
    np.savetxt(lamp, _synthetic_lamp())
    np.savetxt(reference, CATALOG)
    window = LineIdentifier(str(lamp), reference_path=str(reference))
    try:
        for index, peak in enumerate(window.peaks):
            peak.snr = 20.0 if index < 3 else 40.0 if index < 6 else 100.0
        detected_peaks = list(window.peaks)

        def assert_matches_respect_snr(threshold):
            for pixel, _wave, _label in window.auto:
                nearest = min(detected_peaks, key=lambda peak: abs(peak.pixel - pixel))
                assert abs(nearest.pixel - pixel) < 1.0
                assert nearest.snr >= threshold

        window.auto_identify()
        app.processEvents()
        assert len(window.auto) >= 7
        assert_matches_respect_snr(50.0)
        high_threshold_count = len(window.auto)
        window.snr_spin.setValue(30.0)
        window.auto_identify()
        app.processEvents()
        assert len(window.auto) > high_threshold_count
        assert_matches_respect_snr(30.0)
        window.snr_spin.setValue(50.0)
        window.auto_identify()
        app.processEvents()
        assert len(window.auto) == high_threshold_count
        assert_matches_respect_snr(50.0)
        legend = window.ax.get_legend()
        assert legend is not None
        assert legend._loc == 1
        window.canvas.draw()
        legend_box = legend.get_window_extent()
        axes_box = window.ax.get_window_extent()
        assert axes_box.x0 <= legend_box.x0 < legend_box.x1 <= axes_box.x1
        assert axes_box.y0 <= legend_box.y0 < legend_box.y1 <= axes_box.y1
        black_connectors = [
            artist for artist in window.ax.collections
            if len(artist.get_colors())
            and np.allclose(artist.get_colors()[0, :3], [0.0, 0.0, 0.0])
        ]
        assert len(black_connectors) >= len(window.auto)
    finally:
        window.close()
        app.processEvents()


def test_parser_exposes_reference_and_auto():
    parser = build_parser()
    args = parser.parse_args(['--lamp', 'x.fits', '--ref', 'fear.fits', '--auto'])
    assert args.auto is True
    assert args.reference == 'fear.fits'
    assert args.lamp == 'x.fits'
    assert args.min_snr == 50.0
    assert args.flip_x is False
    assert args.poly_type == 'legendre'
    assert args.arm is None
    assert parser.parse_args(['--arm', 'B']).arm == 'B'
    assert parser.parse_args(['--flip-x']).flip_x is True
    assert parser.parse_args(['--poly-type', 'chebyshev']).poly_type == 'chebyshev'
    assert parser.parse_args(['--poly-order', '5']).poly_order == 5


def test_headless_identify_writes_anchors_and_plot(lamp_files, tmp_path):
    lamp, ref = lamp_files
    out = tmp_path / 'anchors.csv'
    plot = tmp_path / 'diagnostic.png'
    rc = main(['--auto', '--lamp', str(lamp), '--ref', str(ref),
               '--wave-range', '5000', '9500', '--min-snr', '5',
               '--poly-type', 'chebyshev', '--poly-order', '5',
               '--output', str(out), '--plot', str(plot)])
    assert rc == 0
    assert out.exists() and plot.exists() and plot.stat().st_size > 1000

    header = [ln for ln in out.read_text().splitlines() if ln.startswith('#')]
    data = [ln for ln in out.read_text().splitlines() if ln and not ln.startswith('#')]
    assert any('solution' in ln for ln in header)
    assert any('chebyshev degree 5' in ln for ln in header)
    assert len(data) >= 10
    first = [part.strip() for part in data[0].split(',')]
    assert len(first) == 3
    assert float(first[0]) > 0 and float(first[1]) > 0


def test_missing_inputs_fail_cleanly(tmp_path, capsys):
    rc = main(['--auto', '--output', str(tmp_path / 'x.csv')])
    assert rc == 2


def test_solve_from_anchors_refines(lamp_files):
    lamp, ref = lamp_files
    flux, _pixel, _info = load_observed(str(lamp))
    reference = load_reference(str(ref))
    catalog = build_catalog(reference)
    peaks = detect_peaks(flux, snr=5.0)

    # Three anchors taken from the known solution, then match the rest.
    pix = np.arange(len(flux), dtype=float)
    truth_px = np.interp(CATALOG, np.polyval(TRUE_COEF, pix), pix)
    anchors = [(float(truth_px[i]), float(CATALOG[i])) for i in (2, 8, 14)]
    solution = solve_from_anchors(peaks, catalog, anchors, len(flux))

    assert solution is not None
    assert solution.n_lines >= 10
    assert solution.rms < 0.5
    # the refinement must stay on the true alias
    assert abs(solution.wavelength(0) - np.polyval(TRUE_COEF, 0)) < 30


def test_solution_from_pairs():
    pairs = np.array([[10.0, 5000.0], [1010.0, 6000.0], [2010.0, 7000.0]])
    solution = solution_from_pairs(pairs, 2048, order=2)
    assert solution.n_lines == 3
    assert abs(solution.rms) < 1e-6
    assert abs(solution.wavelength(1010.0) - 6000.0) < 1e-6


def test_solution_rms_uses_inliers():
    x = np.arange(20, dtype=float) * 50.0
    y = 5000.0 + 1.1 * x
    y[-1] += 100.0
    pairs = np.column_stack([x, y])
    solution = solution_from_pairs(pairs, 1024, order=1)
    # The last point is a clear outlier and must not enter the reported RMS.
    assert solution.n_rms < solution.n_lines
    assert solution.rms < np.sqrt(np.mean(solution.residuals ** 2))
    explicit = solution_from_pairs(pairs, 1024, order=1,
                                   rms_mask=np.arange(20) < 19)
    assert explicit.n_rms == 19
    assert abs(explicit.rms - np.sqrt(np.mean(explicit.residuals[:19] ** 2))) < 1e-12


def test_reference_catalog_uses_wavelength_column():
    wave = np.linspace(5000.0, 6000.0, 1001)
    truth = np.array([5200.0, 5500.0, 5800.0])
    flux = np.full(len(wave), 100.0)
    for centre in truth:
        flux += 800.0 * np.exp(-0.5 * ((wave - centre) / 1.5) ** 2)
    catalog = build_catalog(Reference(wave, flux), snr=5.0,
                            strength_frac=0.0, min_strength=0.0)
    assert len(catalog) == len(truth)
    assert np.allclose(np.sort(catalog.wave), truth, atol=0.5)


def test_resolution_from_gaussian_line_widths():
    n = 1000
    x = np.arange(n, dtype=float)
    centres = np.array([200.0, 400.0, 600.0, 800.0])
    flux = np.full(n, 10.0)
    for centre in centres:
        flux += 1000.0 * np.exp(-0.5 * ((x - centre) / 2.0) ** 2)
    fit = fit_dispersion(centres, 5000.0 + 2.0 * centres, 3,
                         'legendre', domain=(0.0, float(n - 1)))
    info = measure_resolution(x, flux, centres, fit)
    assert info is not None
    # sigma=2 px, FWHM=4.709 px, dispersion=2 A/px -> FWHM ~9.42 A
    assert abs(info['fwhm_ang'] - 9.42) < 0.4
    assert 5800.0 / 9.42 * 0.9 < info['resolution'] < 5800.0 / 9.42 * 1.1
    assert len(info['waves']) == len(info['resolutions']) == info['n']


def test_resolution_rejects_obvious_outlier():
    n = 1000
    x = np.arange(n, dtype=float)
    centres = np.array([200.0, 400.0, 600.0, 800.0, 900.0])
    flux = np.full(n, 10.0)
    for centre in centres:
        sigma = 6.0 if centre == 900.0 else 2.0
        flux += 1000.0 * np.exp(-0.5 * ((x - centre) / sigma) ** 2)
    fit = fit_dispersion(centres, 5000.0 + 2.0 * centres, 3,
                         'legendre', domain=(0.0, float(n - 1)))
    info = measure_resolution(x, flux, centres, fit)
    assert info is not None
    assert info['n_all'] == 5
    assert info['n'] == 4
    assert info['n_rejected'] == 1
    assert np.any(np.isclose(info['rejected_waves'], 6800.0))


def test_reference_catalog_keeps_n_brightest():
    wave = np.linspace(5000.0, 6200.0, 1201)
    centres = np.array([5100.0, 5200.0, 5300.0, 5400.0, 5500.0, 5600.0])
    amplitudes = np.array([100.0, 800.0, 300.0, 900.0, 200.0, 700.0])
    flux = np.full(len(wave), 10.0)
    for centre, amp in zip(centres, amplitudes):
        flux += amp * np.exp(-0.5 * ((wave - centre) / 1.5) ** 2)
    catalog = build_catalog(Reference(wave, flux), snr=5.0,
                            strength_frac=0.0, min_strength=0.0, n_peaks=3)
    assert len(catalog) == 3
    # The three strongest centres are 5400, 5200 and 5600 Angstrom.
    assert np.allclose(np.sort(catalog.wave), [5200.0, 5400.0, 5600.0], atol=0.5)


def test_dispersion_bases_degree5():
    x = np.linspace(100.0, 1900.0, 20)
    y = np.polyval(TRUE_COEF, x)
    for kind in ('legendre', 'chebyshev'):
        fit = fit_dispersion(x, y, 5, kind, domain=(0.0, 2047.0))
        assert dispersion_type(fit) == kind
        assert fit.degree() == 5
        assert np.max(np.abs(fit(x) - y)) < 1e-6
        coef = dispersion_polyval_coef(fit)
        assert np.max(np.abs(np.polyval(coef, x) - y)) < 1e-6
    assert normalize_poly_type('cheb') == 'chebyshev'
    assert normalize_poly_type('leg') == 'legendre'


def test_solve_from_anchors_legendre_degree5(lamp_files):
    lamp, ref = lamp_files
    flux, _pixel, _info = load_observed(str(lamp))
    catalog = build_catalog(load_reference(str(ref)))
    peaks = detect_peaks(flux, snr=5.0)
    pix = np.arange(len(flux), dtype=float)
    truth_px = np.interp(CATALOG, np.polyval(TRUE_COEF, pix), pix)
    anchors = [(float(truth_px[i]), float(CATALOG[i])) for i in (2, 8, 14)]

    solution = solve_from_anchors(peaks, catalog, anchors, len(flux),
                                  order=5, poly_type='legendre')

    assert solution is not None
    assert solution.n_lines >= 15
    assert solution.rms < 0.5
    assert abs(solution.wavelength(0) - np.polyval(TRUE_COEF, 0)) < 30


def test_solve_auto_handles_reverse_dispersion(reverse_lamp_files):
    lamp, ref = reverse_lamp_files
    flux, _pixel, _info = load_observed(str(lamp))
    peaks = detect_peaks(flux, snr=8.0, min_width=0.5, max_width=10.0,
                         distance=4, baseline_size=61)

    solution = solve_auto(peaks, Reference(CATALOG), len(flux),
                          wave_range=(5000.0, 9500.0), order=3,
                          match_tol=2.0, min_tol=0.5)

    assert solution is not None
    assert solution.n_lines >= 15
    assert solution.rms < 0.5
    # wavelength decreases with pixel: the mirrored fallback found the alias
    assert solution.wave_range[0] > solution.wave_range[1]


def test_gui_matches_offscreen(lamp_files):
    pytest.importorskip('PyQt5')
    os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
    QtWidgets = pytest.importorskip('PyQt5.QtWidgets')
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    from specproc.utils.identify import LineIdentifier

    lamp, ref = lamp_files
    window = LineIdentifier(str(lamp), reference_path=str(ref),
                            wave_range=(5000.0, 9500.0))
    try:
        window.auto_identify()
        assert len(window.auto) >= 10
        assert window.display_flux().max() <= 1.0001
        assert window.ax.xaxis.get_label_position() == 'top'
        assert window.disp_fit is not None
        assert window.disp_fit.degree() == window.fit_order_spin.value()
        assert window.resolution_info is not None
        assert window.resolution_info['resolution'] > 0
        assert len(window.resolution_info['waves']) == window.resolution_info['n']
        assert len(window.res_ax.lines) >= 2
        residual_text = [text.get_text() for text in window.resid_ax.texts]
        resolution_text = [text.get_text() for text in window.res_ax.texts]
        assert not any('R ≈' in text for text in residual_text)
        assert any('median R' in text and 'FWHM' in text
                   for text in resolution_text)
        assert window.fit_type() == 'legendre'
        assert 'P' in window.fit_formula.text()
        assert not window.flip_x_box.isChecked()   # wavelength rises with pixel
        window.ax.set_xlim(100.0, 500.0)
        assert np.allclose(
            window.resid_ax.get_xlim(),
            window.disp_fit(np.asarray(window.ax.get_xlim())))
        wave_limits = window.disp_fit(np.array([200.0, 400.0]))
        window.resid_ax.set_xlim(*wave_limits)
        assert np.allclose(window.ax.get_xlim(), [200.0, 400.0], atol=0.1)
        window.ax.set_xlim(-0.5, len(window.flux) - 0.5)
        solution = window.build_solution()
        rows = window._anchor_rows(solution)
        assert len(rows) == len(window.auto)
        assert solution.rms < 0.5
        # the saved polynomial is the selected fit converted for the pipeline
        assert abs(float(window.disp_fit(500.0)) - solution.wavelength(500.0)) < 1e-6
        # the table carries one row per identified line
        assert window.table.rowCount() == len(rows)
        # the basis and degree are adjustable from the window
        window.fit_type_combo.setCurrentIndex(
            window.fit_type_combo.findData('chebyshev'))
        assert window.fit_type() == 'chebyshev'
        assert 'T' in window.fit_formula.text()
        window.fit_order_spin.setValue(2)
        assert window.disp_fit.degree() == 2
        # a decreasing model flips the display axis, and so does its opposite
        window._set_flip_x_from_endpoints(9000.0, 5000.0)
        assert window.flip_x_box.isChecked()
        assert np.allclose(
            window.resid_ax.get_xlim(),
            window.disp_fit(np.asarray(window.ax.get_xlim())))
        window._set_flip_x_from_endpoints(5000.0, 9000.0)
        assert not window.flip_x_box.isChecked()
        assert np.allclose(
            window.resid_ax.get_xlim(),
            window.disp_fit(np.asarray(window.ax.get_xlim())))
        assert window.snr_spin.value() == 50.0
        assert window.clip_mode() == 'sigma'
        assert window.fit_clip_spin.value() == 3.0
        assert [window.clip_mode_combo.itemData(i)
                for i in range(window.clip_mode_combo.count())] == ['sigma', 'absolute']
        window.clip_mode_combo.setCurrentIndex(1)
        assert window.clip_mode() == 'absolute'
        assert window.fit_clip_spin.value() == 0.4
        window.clip_mode_combo.setCurrentIndex(0)
        assert window.clip_mode() == 'sigma'
        assert window.fit_clip_spin.value() == 3.0
        # selection is a two-step state and Esc undoes one step at a time
        window.sel_peak = 500.0
        window.sel_ref = 6000.0
        window.cancel_last_selection()
        assert window.sel_peak == 500.0 and window.sel_ref is None
        window.cancel_last_selection()
        assert window.sel_peak is None
        # +/- buttons use the zoom-step slider
        window.ax.set_xlim(-0.5, len(window.flux) - 0.5)
        window.zoom_slider.setValue(10)
        window._zoom_button(1)
        assert window.ax.get_xlim()[1] - window.ax.get_xlim()[0] < len(window.flux)
        # A panned view keeps its centre, even when it extends past the detector.
        window.ax.set_xlim(5000.0, 6000.0)
        before_centre = sum(window.ax.get_xlim()) / 2.0
        window._zoom_button(1)
        after_centre = sum(window.ax.get_xlim()) / 2.0
        assert abs(after_centre - before_centre) < 1e-6
        # toolbar Home returns to the full detector view, not the last zoom
        window.ax.set_xlim(300.0, 600.0)
        window.toolbar.home()
        assert window.ax.get_xlim() == (-0.5, len(window.flux) - 0.5)
    finally:
        window.close()
        window.deleteLater()
        app.processEvents()


def test_gui_reference_axis_before_and_after_fit(lamp_files):
    pytest.importorskip('PyQt5')
    os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
    QtWidgets = pytest.importorskip('PyQt5.QtWidgets')
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    from specproc.utils.identify import LineIdentifier

    lamp, ref = lamp_files
    window = LineIdentifier(str(lamp), reference_path=str(ref),
                            wave_range=(5000.0, 9500.0))
    try:
        app.processEvents()
        # Before a fit the reference lives on its own wavelength axis.
        assert window.disp_fit is None
        assert window.ref_axis is not None
        first_wave = float(window.overlay_lines()[0][0])
        assert window._pick_reference(first_wave, window.ref_axis)
        # Two anchors create a model and the reference is mapped to pixels.
        pix = np.arange(len(window.flux), dtype=float)
        truth_px = np.interp(CATALOG[:2], np.polyval(TRUE_COEF, pix), pix)
        window.manual = [(float(truth_px[i]), float(CATALOG[i])) for i in (0, 1)]
        window.update_model()
        assert window.disp_fit is not None
        assert window.ref_axis is None
    finally:
        window.close()
        window.deleteLater()
        app.processEvents()


def test_gui_auto_flips_reverse_dispersion(reverse_lamp_files):
    pytest.importorskip('PyQt5')
    os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
    QtWidgets = pytest.importorskip('PyQt5.QtWidgets')
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    from specproc.utils.identify import LineIdentifier

    lamp, ref = reverse_lamp_files
    window = LineIdentifier(str(lamp), reference_path=str(ref),
                            wave_range=(5000.0, 9500.0))
    try:
        window.auto_identify()
        assert len(window.auto) >= 10
        assert window.flip_x_box.isChecked()
        x0, x1 = window.ax.get_xlim()
        assert x0 > x1
        assert np.allclose(
            window.resid_ax.get_xlim(),
            window.disp_fit(np.asarray(window.ax.get_xlim())))
        window.flip_x_box.setChecked(False)
        assert window.ax.get_xlim()[0] < window.ax.get_xlim()[1]
        assert np.allclose(
            window.resid_ax.get_xlim(),
            window.disp_fit(np.asarray(window.ax.get_xlim())))
        window.flip_x_box.setChecked(True)
        assert window.ax.get_xlim()[0] > window.ax.get_xlim()[1]
        solution = window.build_solution()
        assert solution.wave_range[0] > solution.wave_range[1]
    finally:
        window.close()
        window.deleteLater()
        app.processEvents()
