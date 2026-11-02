import numpy as np
from astropy.io import fits

from specproc.config.config_manager import ConfigManager
from specproc.core.slit_pipeline import SlitPipeline


def _pipeline_with_shift(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    pipeline = SlitPipeline(ConfigManager('MSTG193_slit'))
    ny, nx = 160, 80
    xref, xhalf = (nx - 1) / 2, nx / 2

    def design(x, y):
        return (((np.asarray(x) - xref) / xhalf)[..., None])

    pipeline._shift_model = {
        'design': design,
        'coef': np.array([0.8]),
        'band': (8, 72),
    }
    pipeline._flat_trace = (
        np.full(ny, xref),
        np.full(ny, 8.0),
        np.full(ny, 72.0),
    )
    pipeline.sky_gap = 8
    pipeline.sky_edge_margin = 2
    pipeline.sky_order = 1
    pipeline.sky_max_order = 1
    pipeline._band = (8, 72)
    return pipeline, ny, nx


def test_curvature_sky_subtraction_uses_common_wavelength_rows(tmp_path, monkeypatch):
    pipeline, ny, nx = _pipeline_with_shift(tmp_path, monkeypatch)
    rows = np.arange(ny, dtype=float)[:, None]
    columns = np.arange(nx, dtype=float)[None, :]
    shift = pipeline._shift(pipeline._shift_model, columns, rows)
    common_rows = rows - shift
    image = (20.0 + 0.1 * columns
             + 100.0 * np.exp(-0.5 * ((common_rows - 70.0) / 1.8) ** 2)
             + 65.0 * np.exp(-0.5 * ((common_rows - 112.0) / 2.2) ** 2))

    corrected_curved, _ = pipeline._subtract_sky(image, pipeline._flat_trace[0],
                                                  pipeline._band)
    curved_rms = np.std(corrected_curved[10:-10, 8:72])

    pipeline._shift_model = None
    corrected_rowwise, _ = pipeline._subtract_sky(image, pipeline._flat_trace[0],
                                                   pipeline._band)
    rowwise_rms = np.std(corrected_rowwise[10:-10, 8:72])
    assert curved_rms < 0.5
    assert curved_rms < rowwise_rms / 5


def test_debias_keeps_native_pixels_with_curvature_model(tmp_path, monkeypatch):
    pipeline, ny, nx = _pipeline_with_shift(tmp_path, monkeypatch)
    image = np.arange(ny * nx, dtype=float).reshape(ny, nx)
    bias_path = tmp_path / 'frame.fits'
    fits.PrimaryHDU(image).writeto(bias_path)

    corrected, header = pipeline._debias(str(bias_path), np.zeros_like(image), write=False)

    np.testing.assert_array_equal(corrected, image)
    assert 'RECTIFY' not in header


def test_source_trace_follows_slit_edge_trend_when_star_fades(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    pipeline = SlitPipeline(ConfigManager('MSTG193_slit'))
    pipeline.trace_nbins = 40
    ny, nx = 2000, 200
    rows = np.arange(ny, dtype=float)
    flat_centre = 100.0 - 0.005 * rows
    truth = flat_centre + 3.0
    rng = np.random.default_rng(1)
    amp = 50.0 * np.exp(-rows / 400.0)
    cols = np.arange(nx, dtype=float)[None, :]
    img = rng.normal(0, 1.0, (ny, nx))
    img += amp[:, None] * np.exp(-0.5 * ((cols - truth[:, None]) / 3.0) ** 2)
    trace = pipeline._trace_source(img, (5, nx - 5), flat_centre + 3.0)
    assert np.max(np.abs(trace - truth)) < 3.0
    assert np.ptp((trace - flat_centre)) < 3.0


def test_mstg_full_linelist_exists_and_loads(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    pipeline = SlitPipeline(ConfigManager('MSTG193_slit'))
    catalog = pipeline._load_catalog()
    assert len(catalog) > 100


def test_sky_fit_follows_curved_cross_slit_profile(tmp_path, monkeypatch):
    pipeline, ny, nx = _pipeline_with_shift(tmp_path, monkeypatch)
    pipeline._shift_model = None
    columns = np.arange(nx, dtype=float)[None, :]
    rows = np.arange(ny, dtype=float)[:, None]
    profile = 1.0 - 0.00045 * (columns - 40.0) ** 2
    image = 20.0 + 600.0 * profile * ((rows > 60) & (rows < 100))
    pipeline.sky_max_order = 1
    linear, _ = pipeline._subtract_sky(image, pipeline._flat_trace[0], pipeline._band)
    pipeline.sky_max_order = 4
    adaptive, _ = pipeline._subtract_sky(image, pipeline._flat_trace[0], pipeline._band)
    band = (slice(65, 95), slice(8, 72))
    assert np.abs(adaptive[band]).max() < 0.05 * np.abs(linear[band]).max()
