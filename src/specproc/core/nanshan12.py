"""Automatic resolution handling for the Nanshan 1.2 m fibre echelle spectrograph.

The instrument is operated in three resolution modes (high / mid / low) that cannot be
told apart from the FITS headers. When ``[data] auto_resolution = yes`` (the packaged
``nanshan12`` profile) the input frames are split into one group per resolution:

* flat / object / standard frames: from the order spacing of the bright bands, measured
  as the first peak of the autocorrelation of the cross-dispersion profile;
* lamp frames: by cross-correlating every extracted order with the calibrated reference
  lamp spectra of the instrument (``nanshan-thar-{high,mid}res.txt``); a lamp matching
  neither is assigned to the low-resolution mode, whose reference has no flux;
* bias frames: from the ``high`` token in the file name (bias frames carry no other cue).

Each group is then reduced by the standard echelle pipeline with resolution-specific
settings (reference spectrum for the initial wavelength solution, tracing parameters).
"""

import copy
import logging
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

logger = logging.getLogger(__name__)

RESOLUTIONS = ('high', 'mid', 'low')

_REFERENCE_DIR = 'calib_data/telescopes/nanshan12/'
_REFERENCE_FILES = {
    'high': _REFERENCE_DIR + 'nanshan-thar-highres.txt',
    'mid': _REFERENCE_DIR + 'nanshan-thar-medres.txt',
    'low': _REFERENCE_DIR + 'nanshan-thar-lowres.txt',
}

# (section, key) -> value applied on top of the profile for each resolution group.
_RESOLUTION_OVERRIDES = {
    'high': {},
    'mid': {},
    'low': {
        ('reduce.trace', 'snr_threshold'): 8.0,
        # Wide curved bands: follow each one with a narrow window so the tracker cannot
        # jump to a neighbouring order.
        ('reduce.trace', 'max_search_half'): 5,
        ('reduce.trace', 'min_trace_coverage'): 0.8,
        ('reduce.trace', 'seed_columns'): '0.15,0.5,0.85',
        ('reduce.trace', 'side_min_coverage'): 0.3,
        ('reduce.trace', 'degree'): 2,
        # The faintest, lowest band carries only noise.
        ('reduce.trace', 'exclude_lowest_apertures'): 1,
        # Few orders with ~10-20 identifiable lines each: fit every order independently with a
        # low-degree polynomial; orders that fit poorly fall back to the rigid 2D model.
        ('reduce.wlcalib', 'xorder'): 3,
        ('reduce.wlcalib', 'yorder'): 1,
        ('reduce.wlcalib', 'fit_mode'): 'per_order',
        ('reduce.wlcalib', 'order_degree'): 3,
        ('reduce.wlcalib', 'order_min_lines'): 6,
        ('reduce.wlcalib', 'order_min_span'): 0.4,
        ('reduce.wlcalib', 'order_max_rms'): 0.3,
    },
}

# Order spacing (pixels) boundaries between the modes; measured: high 10-12, mid 16-17, low 85-100.
_HIGH_MID_SPACING = 13.5
_MID_LOW_SPACING = 40.0

_MIN_LAMP_ORDERS = 5


def order_spacing(image: np.ndarray) -> float:
    """Return the dominant order spacing (rows) from the cross-dispersion profile, 0 if none."""
    from scipy.ndimage import median_filter
    from scipy.signal import find_peaks

    ny, nx = image.shape
    half = min(300, nx // 2)
    # The low percentile suppresses emission lines in lamp frames and cosmic rays everywhere.
    profile = np.percentile(image[:, nx // 2 - half:nx // 2 + half], 25, axis=1)
    profile = profile - median_filter(profile, size=151)
    profile = profile - profile.mean()
    if not np.any(profile):
        return 0.0
    ac = np.correlate(profile, profile, 'full')[ny - 1:]
    ac = ac / ac[0]
    peaks, _ = find_peaks(ac[6:400], prominence=0.05)
    return float(peaks[0] + 6) if len(peaks) else 0.0


def classify_by_spacing(image: np.ndarray) -> Optional[str]:
    """Classify a flat/object/standard frame as 'high', 'mid' or 'low' (None if undecided)."""
    spacing = order_spacing(image)
    if spacing <= 0:
        return None
    if spacing < _HIGH_MID_SPACING:
        return 'high'
    return 'mid' if spacing < _MID_LOW_SPACING else 'low'


def _extract_for_match(image: np.ndarray, apertures) -> Dict[int, np.ndarray]:
    """Fast sum extraction of ``image`` along traced apertures (used only for classification)."""
    from numpy.polynomial.chebyshev import chebval

    ny, nx = image.shape
    x = np.arange(nx)
    spectra = {}
    for ap in apertures.apertures.values():
        lo, hi = ap.domain
        yc = chebval(2.0 * (x - lo) / max(hi - lo, 1e-6) - 1.0, ap.center_coef)
        half = max(1, int(round(0.5 * ap.width)))
        offsets = np.arange(-half, half + 1)
        rows = np.clip(np.rint(yc[None, :] + offsets[:, None]).astype(int), 0, ny - 1)
        spectra[ap.aperture] = image[rows, x[None, :]].sum(axis=0)
    return spectra


def _group_apertures(flat_files: List[str], config, flip: bool):
    """Trace the orders of the median flat of a resolution group."""
    from ..utils.fits_io import read_fits_image
    from .order_tracing import GratingOrderTracer, renumber_apertures_with_gaps

    flats = [read_fits_image(f)[0].astype(float) for f in flat_files]
    master = np.median(np.stack(flats), axis=0)
    if flip:
        master = np.flipud(master)
    tracer = GratingOrderTracer(
        data=master,
        snr_threshold=config.get_float('reduce.trace', 'snr_threshold', 5.0),
        min_trace_coverage=config.get_float('reduce.trace', 'min_trace_coverage', 0.2),
        trace_degree=config.get_int('reduce.trace', 'degree', 3),
    )
    apertures = tracer.trace_orders()
    return renumber_apertures_with_gaps(apertures, config.get_float('reduce.trace', 'gap_fill_factor', 1.6))


def classify_lamp(lamp_file: str, apertures_by_res: Dict[str, object], flip: bool) -> Optional[str]:
    """Pick the resolution whose reference lamp spectrum best matches the lamp frame.

    Only 'high' and 'mid' have a flux-bearing reference. A lamp that correlates with
    neither (fewer than ``_MIN_LAMP_ORDERS`` orders) is returned as 'low'.
    """
    from ..utils.fits_io import read_fits_image
    from .wave_calibration import WavelengthCalibrator, _resolve_data_path

    image = read_fits_image(lamp_file)[0].astype(float)
    if flip:
        image = np.flipud(image)

    best, best_score = None, 0.0
    scores = {}
    for res, apertures in apertures_by_res.items():
        if res == 'low' or apertures is None or apertures.norders == 0:
            continue
        calibrator = WavelengthCalibrator()
        spectra = _extract_for_match(image, apertures)
        calibrator.match_reference_spectrum(spectra, _resolve_data_path(_REFERENCE_FILES[res]))
        n_corr, score = calibrator.last_match_orders, calibrator.last_match_score
        scores[res] = (n_corr, score)
        if n_corr >= _MIN_LAMP_ORDERS and score > best_score:
            best, best_score = res, score
    logger.info(f"Lamp {Path(lamp_file).name}: reference match (orders, score) {scores}")
    if best is None:
        return 'low'
    return best


def plan_groups(grouped: Dict[str, List[str]], config) -> Dict[str, Dict[str, List[str]]]:
    """Split classified input files into per-resolution reduction groups."""
    from ..utils.fits_io import read_fits_image

    flip = config.get_bool('data', 'flip_vertical', False)
    groups: Dict[str, Dict[str, List[str]]] = {
        res: {'bias': [], 'flat': [], 'calib': [], 'science': []} for res in RESOLUTIONS}

    for cat in ('flat', 'science'):
        for path in grouped[cat]:
            try:
                res = classify_by_spacing(read_fits_image(path)[0].astype(float))
            except Exception as exc:
                logger.warning(f"Cannot classify {Path(path).name}: {exc}")
                continue
            if res is None:
                logger.warning(f"Resolution of {Path(path).name} could not be determined; skipped.")
                continue
            groups[res][cat].append(path)
            logger.info(f"{cat:8s} {Path(path).name}: {res} resolution (order spacing)")

    apertures_by_res = {}
    for res in ('high', 'mid'):
        if groups[res]['flat']:
            try:
                apertures_by_res[res] = _group_apertures(groups[res]['flat'], config, flip)
            except Exception as exc:
                logger.warning(f"Tracing the {res}-resolution flats for lamp matching failed: {exc}")

    for path in grouped['calib']:
        res = classify_lamp(path, apertures_by_res, flip)
        if res is None or not groups[res]['flat']:
            logger.warning(f"Lamp {Path(path).name}: no flat for resolution '{res}'; skipped.")
            continue
        groups[res]['calib'].append(path)
        logger.info(f"calib    {Path(path).name}: {res} resolution (reference spectrum)")

    for path in grouped['bias']:
        high = 'high' in Path(path).name.lower()
        for res in RESOLUTIONS:
            if (res == 'high') == high:
                groups[res]['bias'].append(path)
    for res in RESOLUTIONS:
        if not groups[res]['bias']:
            groups[res]['bias'] = list(grouped['bias'])

    return {res: g for res, g in groups.items() if g['flat'] or g['calib'] or g['science']}


def group_config(config, resolution: str):
    """Return a config copy configured for one resolution group."""
    cfg = copy.deepcopy(config)
    cfg.set('telescope.linelist', 'reference_spectrum', _REFERENCE_FILES[resolution])
    cfg.set('telescope.linelist', 'anchor_file', '')
    for (section, key), value in _RESOLUTION_OVERRIDES[resolution].items():
        cfg.set(section, key, value)
    base = Path(config.get_output_path())
    cfg.set('reduce', 'output_path', (base / resolution).as_posix())
    return cfg


def run_by_resolution(config, grouped: Dict[str, List[str]], pipeline_cls) -> Dict[str, dict]:
    """Reduce every resolution group with its own configuration and output directory."""
    plan = plan_groups(grouped, config)
    results = {}
    failures = {}
    for res, files in plan.items():
        counts = ', '.join(f"{k}={len(v)}" for k, v in files.items())
        missing = [k for k in ('bias', 'flat', 'calib', 'science') if not files[k]]
        if missing:
            logger.warning(f"Resolution '{res}' skipped ({counts}); missing {', '.join(missing)}.")
            continue
        logger.info(f"=== Nanshan 1.2 m {res}-resolution group ({counts}) ===")
        try:
            pipeline = pipeline_cls(group_config(config, res))
            results[res] = pipeline.run_full_pipeline(
                raw_image_paths=files['science'], bias_filenames=files['bias'],
                flat_filenames=files['flat'], calib_filenames=files['calib'])
        except Exception as exc:
            # One failing mode must not discard the results of the others.
            failures[res] = exc
            logger.error(f"Nanshan {res}-resolution group failed: {exc}")
    if not results:
        raise RuntimeError('No Nanshan resolution group was reduced successfully'
                           + (f": {failures}" if failures else '.'))
    if failures:
        logger.warning('Groups that failed: ' + ', '.join(failures))
    return results
