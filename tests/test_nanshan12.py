"""Tests for Nanshan 1.2 m resolution classification and per-order wavelength fits."""

import numpy as np

from specproc.core import nanshan12
from specproc.core.wave_calibration import WavelengthCalibrator


def _banded_image(spacing, ny=1200, nx=800, width=2.0):
    y = np.arange(ny)[:, None]
    img = np.zeros((ny, nx))
    for yc in np.arange(40, ny - 40, spacing):
        img += 1000.0 * np.exp(-0.5 * ((y - yc) / width) ** 2)
    img += np.random.default_rng(0).normal(0, 5, img.shape)
    return img


def test_resolution_from_order_spacing():
    assert nanshan12.classify_by_spacing(_banded_image(11)) == 'high'
    assert nanshan12.classify_by_spacing(_banded_image(17, width=3.0)) == 'mid'
    assert nanshan12.classify_by_spacing(_banded_image(90, width=8.0)) == 'low'


def test_per_order_fit_beats_global_when_orders_differ():
    rng = np.random.default_rng(1)
    calibrator = WavelengthCalibrator()
    pix, waves = [], []
    for m in range(40, 60):
        x = np.sort(rng.uniform(0, 2000, 40))
        # Order-specific distortion that no smooth f(order, x) of low degree can follow.
        shift = rng.normal(0, 0.3)
        lam = 6000.0 - 2.0 * (m - 40) + 0.2 * x + 1e-5 * x ** 2 + shift * (x / 2000.0)
        pix += [[xi, m] for xi in x]
        waves += list(lam)
    pix, waves = np.array(pix), np.array(waves)

    wc = calibrator.fit_wavelength_polynomial(pix, waves, xorder=3, yorder=3, max_residual=2.0)
    rms_global = wc.rms
    calibrator.fit_per_order(wc, degree=3, max_residual=2.0)

    assert wc.order_coefs and len(wc.order_coefs) == 20
    assert wc.rms < 0.1 * rms_global
    pred = wc.apply_to_pixel(pix[:, 0], pix[:, 1])
    assert np.sqrt(np.mean((pred - waves) ** 2)) < 1e-6


def test_per_order_solution_round_trip(tmp_path):
    rng = np.random.default_rng(2)
    calibrator = WavelengthCalibrator()
    x = np.tile(np.linspace(0, 1000, 30), 3)
    m = np.repeat([10, 11, 12], 30)
    waves = 5000.0 + 0.5 * x + 20.0 * m + rng.normal(0, 0.001, x.size)
    wc = calibrator.fit_wavelength_polynomial(np.column_stack([x, m]), waves, xorder=2, yorder=2)
    calibrator.fit_per_order(wc, degree=2, min_lines=8)
    calibrator.save_calibration(str(tmp_path / 'wcal.fits'))

    loaded = WavelengthCalibrator().load_calibration(str(tmp_path / 'wcal.fits'))
    assert set(loaded.order_coefs) == {10, 11, 12}
    np.testing.assert_allclose(loaded.apply_to_pixel(x, m), wc.apply_to_pixel(x, m))


def test_reject_false_traces_and_crop(tmp_path):
    from numpy.polynomial import chebyshev as cheb
    from specproc.core.data_structures import ApertureLocation, ApertureSet
    from specproc.core.order_tracing import reject_false_traces
    from specproc.core.basic_reduction import OverscanCorrectionStage

    w, h = 1000, 700
    xm = 2.0 * np.arange(w) / (w - 1) - 1.0
    aps = ApertureSet()
    for k in range(8):
        y = 300 + 12 * k - 20 * (xm + 1)          # all real orders bend alike
        aps.add_aperture(ApertureLocation(
            aperture=k + 1, order=k + 1, center_coef=cheb.chebfit(xm, y, 2),
            is_chebyshev=True, domain=(0.0, w - 1.0)))
    steep = 100 + 3000 * (xm + 1) ** 2             # leaves the detector
    aps.add_aperture(ApertureLocation(aperture=9, order=9, center_coef=cheb.chebfit(xm, steep, 2),
                                      is_chebyshev=True, domain=(0.0, w - 1.0)))
    odd = 380 + 80 * xm ** 2                       # in the image but bends differently
    aps.add_aperture(ApertureLocation(aperture=10, order=10, center_coef=cheb.chebfit(xm, odd, 2),
                                      is_chebyshev=True, domain=(0.0, w - 1.0)))
    kept = reject_false_traces(aps, w, h).get_orders()
    assert 9 not in kept and 10 not in kept
    assert kept == list(range(1, 9))

    stage = OverscanCorrectionStage(trim_y_start=600, trim_y_end=1300)
    assert stage._apply_user_crop(np.zeros((2048, 2048))).shape == (701, 2048)
    assert OverscanCorrectionStage(trim_y_end=1000)._apply_user_crop(np.zeros((2048, 100))).shape == (1000, 100)


def test_mark_peak_wavelengths_labels_lamp_peaks():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from specproc.plotting.spectra_plotter import mark_peak_wavelengths

    x = np.arange(1000.0)
    flux = np.random.default_rng(2).normal(0, 1, x.size)
    for c, a in ((200.3, 500), (500.0, 900), (800.6, 300)):
        flux += a * np.exp(-0.5 * ((x - c) / 1.5) ** 2)
    wave = 4000.0 + 2.0 * x
    fig, ax = plt.subplots()
    ax.plot(wave, flux)
    assert mark_peak_wavelengths(ax, wave, flux) == 3
    labels = sorted(t.get_text() for t in ax.texts)
    assert labels[0] == '4400.6' + labels[0][6:] and all(len(t.split('.')[1]) == 4 for t in labels)
    plt.close(fig)


def test_exclude_edge_apertures():
    from numpy.polynomial import chebyshev as cheb
    from specproc.core.data_structures import ApertureLocation, ApertureSet
    from specproc.core.order_tracing import _exclude_edge_apertures

    w = 500
    xm = 2.0 * np.arange(w) / (w - 1) - 1.0
    aps = ApertureSet()
    for k in range(5):
        aps.add_aperture(ApertureLocation(aperture=k + 1, order=k + 1,
                                          center_coef=cheb.chebfit(xm, 100 + 50 * k + 0 * xm, 1),
                                          is_chebyshev=True, domain=(0.0, w - 1.0)))
    assert _exclude_edge_apertures(aps, w, 1, 0).get_orders() == [2, 3, 4, 5]
    assert _exclude_edge_apertures(aps, w, 1, 2).get_orders() == [2, 3]
    assert _exclude_edge_apertures(aps, w, 0, 0).get_orders() == [1, 2, 3, 4, 5]


def test_estimate_snr_at_center():
    from specproc.plotting.spectra_plotter import estimate_snr_at_center
    rng = np.random.default_rng(3)
    flux = 1000 + rng.normal(0, 10, 2000)
    idx, snr = estimate_snr_at_center(flux)
    assert idx == 999 and 80 < snr < 125


def test_load_profile_replaces_settings_without_touching_packaged_file(tmp_path, monkeypatch):
    from specproc.config.config_manager import ConfigManager, resolve_telescope_profile

    monkeypatch.chdir(tmp_path)
    cfg = ConfigManager(str(resolve_telescope_profile('XL216_hrs')))
    cfg.load_profile(resolve_telescope_profile('nanshan12'))
    assert cfg.get_bool('data', 'flip_vertical', False)
    assert cfg.config_path == tmp_path / 'default_config_NS12_hrs.cfg'
    assert not cfg.config.has_option('data', 'direction')


def test_match_anchor_positions_uses_table_x_positions():
    import numpy as np
    from specproc.core.wave_calibration import WavelengthCalibrator

    table = np.array([
        [3, 20.0, 7000.0], [3, 100.0, 7100.0], [3, 180.0, 7200.0],
        [3, 260.0, 7300.0], [4, 40.0, 6500.0], [4, 130.0, 6600.0],
        [4, 220.0, 6700.0],
    ])
    peaks = {1: np.array([24.0, 104.0, 184.0, 264.0]),
             2: np.array([44.0, 134.0, 224.0])}
    delta_m, matched = WavelengthCalibrator().match_anchor_positions(peaks, peaks, table)
    assert delta_m == 2
    assert len(matched) == 7
    assert any(abs(x - 104.0) < 1e-6 and m == 3 and w == 7100.0 for x, m, w in matched)
