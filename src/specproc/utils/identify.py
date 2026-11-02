"""
Interactive and automatic lamp-line identification for long-slit spectra.

The top panel shows the normalized observed lamp spectrum (upper curve, pixel
axis on top) and the reference lamp spectrum or line list (lower curve, blue,
wavelength axis at the bottom), with the reference mapped through the current
dispersion model.  Only observed peaks above the SNR threshold are marked, with
red inverted triangles.  The bottom panel shows the residuals of a Legendre or
Chebyshev series fitted to the identified anchor lines; its basis, degree and
residual-clip parameters can be changed in the window.  ``Flip x`` reverses the
displayed pixel direction so spectra whose wavelength decreases with pixel can
still be compared with the reference in the usual blue-to-red sense; automatic
identification sets it from the fitted dispersion.

Load a reference atlas (for example ``calib_data/linelists/fear.fits``) or a
text line list, press "Auto identify" to match all strong lamp peaks, adjust
the result by hand (click a lamp peak, then the matching reference line;
right-click to remove a line), and save the anchor CSV together with a
diagnostic figure.  The first two columns of the CSV (0-based pixel,
wavelength in Angstrom) are exactly the ``[slit] identification_file`` format.

The same matching engine runs without a window with ``specproc identify
--auto``.
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import List, Optional

import numpy as np
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import Qt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure

from specproc.core.anchor_matching import (
    Peak,
    Reference,
    build_catalog,
    detect_peaks,
    load_observed,
    load_reference,
    plot_diagnostic,
    select_anchors,
    solution_from_pairs,
    solve_dispersion,
    solve_from_anchors,
    write_anchor_file,
)
from specproc.core.line_identification import (
    dispersion_polyval_coef,
    dispersion_type,
    fit_dispersion,
    load_anchors,
    load_lamp_spectrum,
    load_linelist,
    measure_resolution,
    normalize_poly_type,
)

logger = logging.getLogger('specproc.identify')

MIN_DETECT_SNR = 3.0

# Colours chosen to stay readable on top of the (usually blue) lamp spectrum.
LAMP_COLOR = '#455a64'
PEAK_COLOR = '#d50000'
MANUAL_COLOR = '#d50000'
AUTO_COLOR = '#ef6c00'
REF_COLOR = '#1565c0'
REF_MARK_COLOR = '#8e24aa'
SELECT_MARK_COLOR = '#00b050'
# Vertical offsets in the normalized-flux panel: measured spectrum on top,
# reference spectrum just below it.
MEASURED_OFFSET = 1.12
REFERENCE_OFFSET = 0.0


# ---------------------------------------------------------------------------
# Command line: shared matching options and automatic mode.
# ---------------------------------------------------------------------------
def add_match_arguments(parser: argparse.ArgumentParser) -> None:
    """Add the input/matching options shared by the window and --auto mode."""
    parser.add_argument('--lamp',
                        help='extracted 1-D lamp spectrum (FITS table or text)')
    parser.add_argument('--ref', dest='reference',
                        help='reference lamp spectrum (e.g. fear.fits) or line list')
    parser.add_argument('--linelist',
                        help='optional line catalogue used for ion labels (e.g. fear.dat)')
    parser.add_argument('--anchors',
                        help='optional existing anchor file (pixel, wavelength) used to '
                             'seed the solution')
    parser.add_argument('--config',
                        help='optional SpecProc config/profile providing [slit] '
                             'wave_start/wave_end (e.g. MSTG193_slit or a .cfg path)')
    parser.add_argument('--wave-range', nargs=2, type=float, metavar=('W0', 'W1'),
                        help='approximate wavelength coverage in Angstrom (first/last '
                             'pixel).  It selects the reference lines used for the match '
                             'and seeds the curvature scan, so it should be within a few '
                             'percent of the true range.  Read from --config or the FITS '
                             'header when omitted; without any hint a blind pattern scan '
                             'is tried.')
    parser.add_argument('--output',
                        help='output anchor CSV (default: <lamp>_anchors.csv)')
    parser.add_argument('--plot',
                        help='output diagnostic figure (default: <lamp>_anchors.pdf; PNG is also allowed)')
    parser.add_argument('--poly-order', type=int, default=5,
                        choices=tuple(range(1, 10)),
                        help='degree of the dispersion polynomial (default 5)')
    parser.add_argument('--poly-type', default='legendre',
                        choices=('legendre', 'chebyshev', 'polynomial'),
                        help='polynomial basis for the dispersion fit '
                             '(default: legendre)')
    parser.add_argument('--resolution', type=float, default=2000.0,
                        help='approximate resolving power, used to reject blended anchors '
                             '(default 2000; 0 disables the check)')
    parser.add_argument('--min-snr', type=float, default=50.0,
                        help='minimum S/N of a lamp peak used as an anchor (default 50)')
    parser.add_argument('--match-tol', type=float, default=2.0,
                        help='initial matching tolerance in Angstrom (default 2.0)')
    parser.add_argument('--max-residual', type=float, default=0.6,
                        help='maximum residual of an anchor against the preliminary fit '
                             '(default 0.6 A; the wavelength itself comes from the catalogue)')
    parser.add_argument('--min-separation', type=float, default=4.0,
                        help='minimum separation from another lamp peak in pixels (default 4)')
    parser.add_argument('--max-anchors', type=int, default=30,
                        help='maximum number of anchors (default 30)')
    parser.add_argument('--catalog-snr', type=float, default=20.0,
                        help='peak-detection S/N in the reference spectrum (default 20)')
    parser.add_argument('--catalog-frac', type=float, default=0.01,
                        help='keep reference lines above this fraction of the strongest '
                             'line inside the requested range (default 0.01)')
    parser.add_argument('--catalog-min-snr', type=float, default=50.0,
                        help='absolute minimum reference-line S/N (default 50)')
    parser.add_argument('--zoom', nargs=2, type=float, metavar=('W0', 'W1'),
                        help='wavelength range for the zoom panel')
    parser.add_argument('--log-level', default='INFO',
                        choices=('DEBUG', 'INFO', 'WARNING', 'ERROR'))


def build_parser(prog: str = 'specproc.identify') -> argparse.ArgumentParser:
    """Parser for the identification tool (window plus --auto mode)."""
    parser = argparse.ArgumentParser(
        prog=prog,
        description='Interactively identify lamp lines against a reference atlas or line '
                    'list, or run the same matching automatically with --auto.')
    add_match_arguments(parser)
    parser.add_argument('--dispersion-axis', choices=['x', 'y'], default='y',
                        help='dispersion axis when --lamp is a 2-D image (default: y)')
    parser.add_argument('--auto', action='store_true',
                        help='run the automatic matcher without opening the window')
    parser.add_argument('--arm', choices=('R', 'B'), default=None,
                        help='MSTG 1.93m camera arm (R or B); by default it is read '
                             'from INSTRUME in the lamp FITS header')
    parser.add_argument('--flip-x', action='store_true',
                        help='start with the displayed pixel axis reversed (GUI; auto '
                             'identify still sets the direction from the fitted dispersion)')
    return parser


def _range_from_config(config_path: Optional[str]) -> Optional[tuple]:
    if not config_path:
        return None
    try:
        from specproc.config.config_manager import ConfigManager
        cfg = ConfigManager(config_path)
        w0 = cfg.get_float('slit', 'wave_start', 0.0)
        w1 = cfg.get_float('slit', 'wave_end', 0.0)
        if w0 > 0 and w1 > w0:
            return w0, w1
    except Exception as exc:  # pragma: no cover - configuration is optional
        logger.warning('Could not read wave_start/wave_end from %s: %s', config_path, exc)
    return None


def _mirror_peaks(peaks, n_pixels: int):
    """Mirror the pixel axis, keeping the peak list sorted by pixel."""
    mirrored = [Peak(n_pixels - 1 - p.pixel, p.amplitude, p.snr, p.width)
                for p in peaks]
    return sorted(mirrored, key=lambda p: p.pixel)


def solve_auto(peaks, catalog, n_pixels: int, wave_range=None, order: int = 3,
               match_tol: float = 2.0, min_tol: float = 0.3):
    """Blind pattern match, supporting either wavelength direction.

    The matcher evaluates both signs of the linear seed in one search. The
    mirrored fallback is retained for cases where the original-axis seed
    search cannot produce a valid solution.
    """
    solution = solve_dispersion(
        peaks, catalog, n_pixels, wave_range=wave_range, order=order,
        match_tol=match_tol, min_tol=min_tol)
    if solution is not None:
        return solution

    reversed_solution = solve_dispersion(
        _mirror_peaks(peaks, n_pixels), catalog, n_pixels,
        wave_range=wave_range, order=order, match_tol=match_tol, min_tol=min_tol)
    if reversed_solution is None:
        return None
    pairs = np.column_stack([n_pixels - 1 - reversed_solution.pairs[:, 0],
                             reversed_solution.pairs[:, 1]])
    pairs = pairs[np.argsort(pairs[:, 0])]
    logger.info('Dispersion increases towards smaller pixels; using the mirrored match')
    return solution_from_pairs(pairs, n_pixels, order=order, peaks=peaks,
                               catalog=catalog)


def _telescope_name(config_path: Optional[str]) -> str:
    if not config_path:
        return ''
    try:
        from specproc.config.config_manager import ConfigManager
        return ConfigManager(str(config_path)).get('telescope', 'name', '') or ''
    except Exception:
        return 'mstg193slit' if 'mstg193' in str(config_path).lower() else ''


def _fits_header_meta(path: Optional[str]):
    """Return (telescope_name, arm) from TELESCOP/INSTRUME when available."""
    if not path:
        return '', ''
    suffix = Path(path).suffix.lower()
    if suffix not in ('.fits', '.fit', '.fts', '.fz'):
        return '', ''
    try:
        from astropy.io import fits
        with fits.open(path, memmap=False) as hdul:
            header = hdul[0].header
            for hdu in hdul[1:]:
                if 'TELESCOP' in hdu.header or 'INSTRUME' in hdu.header:
                    header = hdu.header
                    break
        telescope = str(header.get('TELESCOP', '')).strip().upper()
        instrument = str(header.get('INSTRUME', '')).strip().upper()
    except Exception:
        return '', ''
    name = 'mstg193slit' if 'MOST' in telescope else ''
    arm = ''
    if 'SPECTRUM-R' in instrument or instrument.endswith('-R'):
        arm = 'R'
    elif 'SPECTRUM-B' in instrument or instrument.endswith('-B'):
        arm = 'B'
    return name, arm


def run_auto(args) -> int:
    """Run the automatic matching described by parsed ``args``."""
    logging.basicConfig(level=getattr(logging, args.log_level.upper()),
                        format='%(asctime)s - %(levelname)s - %(message)s')

    if not args.lamp or not args.reference:
        logger.error('Automatic matching needs --lamp and --ref. '
                     'Run `specproc identify` to load files interactively.')
        return 2

    try:
        flux, pixel, info = load_observed(args.lamp)
    except Exception as exc:
        logger.error('Cannot read lamp spectrum: %s', exc)
        return 2
    if len(flux) < 50:
        logger.error('Lamp spectrum is too short (%d pixels)', len(flux))
        return 2

    try:
        ref = load_reference(args.reference, args.linelist)
    except Exception as exc:
        logger.error('Cannot read reference: %s', exc)
        return 2

    label_wave = label_text = None
    if args.linelist:
        label_wave, label_text, _ = load_linelist(args.linelist)

    wave_range = tuple(args.wave_range) if args.wave_range else None
    if wave_range is None:
        wave_range = info.get('wave_range')
    if wave_range is None:
        wave_range = _range_from_config(args.config)
    if wave_range is not None:
        logger.info('Approximate wavelength coverage: %.1f-%.1f A',
                    wave_range[0], wave_range[1])
    else:
        logger.info('No wavelength hint; running the blind line-pattern scan.')

    try:
        if ref.flux is not None and not np.allclose(ref.flux, 1.0):
            catalog_ref = ref
            if wave_range is not None:
                w0, w1 = sorted(map(float, wave_range))
                margin = max(50.0, 0.05 * (w1 - w0))
                catalog_ref = ref.subset(wave_range, margin=margin)
            catalog = build_catalog(
                catalog_ref, snr=args.catalog_snr,
                strength_frac=args.catalog_frac,
                min_strength=args.catalog_min_snr, max_lines=400,
                linelist_wave=label_wave, linelist_labels=label_text)
        else:
            catalog = build_catalog(
                ref, wave_range=wave_range, snr=args.catalog_snr,
                strength_frac=args.catalog_frac,
                min_strength=args.catalog_min_snr,
                linelist_wave=label_wave, linelist_labels=label_text)
    except Exception as exc:
        logger.error('Cannot build the reference catalogue: %s', exc)
        return 2
    if len(catalog) < 5:
        logger.error('Too few strong reference lines (%d); lower --catalog-min-snr',
                     len(catalog))
        return 2

    detect_snr = min(args.min_snr, 8.0)
    peaks = detect_peaks(flux, snr=detect_snr, min_width=0.5, max_width=10.0,
                         distance=4, baseline_size=61)
    logger.info('Detected %d lamp peaks above S/N=%.1f', len(peaks), detect_snr)

    anchor_seed = None
    if args.anchors:
        try:
            anchor_seed = load_anchors(args.anchors)
        except Exception as exc:
            logger.error('Cannot read anchor seed file: %s', exc)
            return 2
        if len(anchor_seed) < 2:
            logger.warning('Anchor file %s has fewer than two lines; ignoring it',
                           args.anchors)
            anchor_seed = None

    peaks = [peak for peak in peaks if peak.snr >= args.min_snr]
    logger.info('%d lamp peaks meet --min-snr %.1f', len(peaks), args.min_snr)
    min_match_peaks = 4 if anchor_seed is not None else 7
    if len(peaks) < min_match_peaks:
        logger.error('Too few lamp peaks meet --min-snr for matching (%d; need %d)',
                     len(peaks), min_match_peaks)
        return 2

    if anchor_seed is not None:
        logger.info('Seeding the solution with %d manual anchors from %s',
                    len(anchor_seed), args.anchors)
        solution = solve_from_anchors(peaks, catalog, anchor_seed, len(flux),
                                      order=args.poly_order, match_tol=args.match_tol,
                                      min_tol=max(0.5, args.match_tol * 0.25),
                                      poly_type=args.poly_type)
    else:
        solution = solve_auto(
            peaks, catalog, len(flux), wave_range=wave_range,
            order=args.poly_order, match_tol=args.match_tol,
            min_tol=max(0.5, args.match_tol * 0.25))

    if solution is None:
        logger.error('No dispersion solution found.  Check --wave-range, lower '
                     '--catalog-min-snr/--catalog-frac, or try --poly-order 2.')
        return 2
    logger.info('Selected solution: %d matches, rms=%.3f A, range %.1f-%.1f A',
                solution.n_lines, solution.rms, *solution.wave_range)

    anchors = select_anchors(
        peaks, solution, catalog,
        min_snr=args.min_snr, max_residual=args.max_residual,
        min_separation=args.min_separation, max_anchors=args.max_anchors,
        resolution=args.resolution)
    if not anchors:
        logger.error('No anchors passed the selection cuts; relax --max-residual, '
                     '--min-separation or --min-snr')
        return 2
    logger.info('Selected %d anchors', len(anchors))
    for snr, x, w, res, label in anchors:
        logger.info('  %8.2f px  %9.3f A  %-5s  S/N=%8.0f  res=%+.2f A',
                    x, w, label or '?', snr, res)

    # Refit the selected anchors in the requested basis.  The coefficients are
    # converted back to ordinary polynomial form for the pipeline, while the
    # final residuals and diagnostic use the Legendre/Chebyshev fit.
    anchor_pairs = np.array([(a[1], a[2]) for a in anchors], dtype=float)
    final_fit = fit_dispersion(anchor_pairs[:, 0], anchor_pairs[:, 1],
                               args.poly_order, args.poly_type,
                               domain=(0.0, float(len(flux) - 1)))
    if final_fit is not None:
        coef = dispersion_polyval_coef(final_fit)
        rms_mask = np.array([
            bool(np.any((np.abs(anchor_pairs[:, 0] - x) < 1e-6)
                        & (np.abs(anchor_pairs[:, 1] - w) < 1e-6)))
            for x, w in solution.pairs], dtype=bool)
        solution = solution_from_pairs(
            solution.pairs, len(flux), order=args.poly_order, peaks=peaks,
            catalog=catalog, coef=coef,
            rms_mask=rms_mask if rms_mask.sum() >= 2 else None)
        anchors = [(snr, x, w, float(w - final_fit(x)), label)
                   for snr, x, w, _res, label in anchors]
        logger.info('Final fit: %s degree %d, %d anchors',
                    dispersion_type(final_fit), final_fit.degree(), len(anchors))
    else:
        logger.warning('Could not refit the anchors with %s; keeping the '
                       'preliminary polynomial', args.poly_type)

    lamp_path = Path(args.lamp)
    header_telescope, header_arm = _fits_header_meta(args.lamp)
    telescope = header_telescope or _telescope_name(args.config)
    arm = args.arm or header_arm or 'R'
    is_mstg = (telescope.lower() == 'mstg193slit'
               or 'mstg193slit' in str(args.reference).lower()
               or 'mstg193slit' in str(args.config).lower())
    if args.output:
        out = Path(args.output)
    elif is_mstg:
        out = Path.cwd() / f'mstg193slit_anchors_{arm}.csv'
    else:
        out = Path.cwd() / (lamp_path.stem + '_anchors.csv')
    fig = (Path(args.plot) if args.plot
           else Path.cwd() / (lamp_path.stem + '_anchors.pdf'))
    note = (f'seeded with {wave_range[0]:.1f}-{wave_range[1]:.1f} A'
            if wave_range else 'blind pattern match')
    fitted_kind = dispersion_type(final_fit) if final_fit is not None else args.poly_type
    fitted_degree = final_fit.degree() if final_fit is not None else args.poly_order
    note += f'; {fitted_kind} degree {fitted_degree}'
    written = write_anchor_file(out, anchors, solution, str(lamp_path),
                                str(args.reference), note=note,
                                created=np.datetime64('today').astype(str))
    title = f'{lamp_path.name}  |  {len(anchors)} anchors  |  rms={solution.rms:.2f} A'
    zoom = tuple(args.zoom) if args.zoom else None
    pixel_wavelength = (info.get('wavelength') if info.get('wavelength') is not None
                        else solution.wavelength(pixel))
    plot_diagnostic(fig, flux, peaks, solution, anchors, ref=ref, title=title,
                    zoom=zoom, pixel_wavelength=pixel_wavelength)
    logger.info('Wrote %s', written)
    logger.info('Wrote %s', fig)
    return 0


class SpectrumToolbar(NavigationToolbar2QT):
    """Toolbar with a pointer/select action plus continuous pan and zoom."""

    @staticmethod
    def _select_icon():
        """A plain pointer arrow in the same black line style as the toolbar."""
        pix = QtGui.QPixmap(24, 24)
        pix.fill(Qt.transparent)
        painter = QtGui.QPainter(pix)
        painter.setRenderHint(QtGui.QPainter.Antialiasing)
        colour = QtGui.QColor('#212121')
        painter.setPen(QtGui.QPen(colour, 1.4, Qt.SolidLine,
                                  Qt.RoundCap, Qt.RoundJoin))
        painter.setBrush(QtGui.QBrush(colour))
        points = [
            QtCore.QPointF(6.0, 3.0),
            QtCore.QPointF(6.0, 20.0),
            QtCore.QPointF(10.0, 16.0),
            QtCore.QPointF(13.0, 22.0),
            QtCore.QPointF(16.0, 20.5),
            QtCore.QPointF(13.0, 15.0),
            QtCore.QPointF(18.5, 15.0),
        ]
        painter.drawPolygon(QtGui.QPolygonF(points))
        painter.end()
        return QtGui.QIcon(pix)

    def __init__(self, canvas, parent, window):
        self.window = window
        super().__init__(canvas, parent)
        self.select_action = QtWidgets.QAction('Select', self)
        self.select_action.setCheckable(True)
        self.select_action.setChecked(True)
        self.select_action.setToolTip('Select peaks and reference lines for pairing')
        self.select_action.setIcon(self._select_icon())
        pan_action = next((action for action in self.actions()
                           if action.text().strip().lower() == 'pan'), None)
        if pan_action is not None:
            self.insertAction(pan_action, self.select_action)
        else:
            self.addAction(self.select_action)
        self.select_action.triggered.connect(self.window.set_select_mode)

    def pan(self, *args):
        super().pan(*args)
        if self.mode and self.window.select_mode and not self.window._changing_tool:
            self.window.set_select_mode(False, update_toolbar=False)

    def zoom(self, *args):
        super().zoom(*args)
        if self.mode and self.window.select_mode and not self.window._changing_tool:
            self.window.set_select_mode(False, update_toolbar=False)

    def home(self, *args):
        self.window.reset_view()
        try:
            self._nav_stack.clear()
            self.push_current()
            self.set_history_buttons()
        except Exception:
            pass


class LineIdentifier(QtWidgets.QMainWindow):
    """Measured and reference spectra on top, dispersion-fit residuals below."""

    def __init__(self, lamp_path: Optional[str] = None,
                 linelist_path: Optional[str] = None,
                 anchor_path: Optional[str] = None, dispersion_axis: str = 'y',
                 reference_path: Optional[str] = None,
                 wave_range: Optional[tuple] = None,
                 config_path: Optional[str] = None,
                 output_path: Optional[str] = None,
                 plot_path: Optional[str] = None,
                 flip_x: bool = False, poly_type: str = 'legendre',
                 arm: Optional[str] = None):
        super().__init__()
        self.setWindowTitle('SpecProc - Lamp Line Identification')
        self.resize(1500, 880)
        self.dispersion_axis = dispersion_axis
        self.flux = np.array([])
        self.peaks: List[Peak] = []
        self.ref: Optional[Reference] = None
        self.ref_is_spectrum = False
        self.ref_strengths = np.array([])
        self.ref_path = ''
        self.manual = []                # (pixel, wavelength) picked by hand
        self.auto = []                  # (pixel, wavelength, label) from auto matching
        self.disp_fit = None             # pixel -> wavelength dispersion model
        self.disp_residuals = np.empty(0)
        self.disp_mask = np.empty(0, dtype=bool)
        self.sel_peak: Optional[float] = None
        self.sel_ref: Optional[float] = None
        self.lamp_path = ''
        self.anchor_path = ''
        self.output_path = output_path or ''
        self.plot_path = plot_path or ''
        self.solution = None
        self.resolution_info = None
        self._plot_in_wave = False
        self._view_in_wave = False
        self._syncing = False
        self._ax_xlim_cid = None
        self._resid_xlim_cid = None
        self._ref_top: Optional[float] = None
        self._ref_texts = []
        self._wave_texts = []
        self._catalog_cache = None
        self.ref_axis = None
        self.telescope_name = ''
        self.select_mode = True
        self._changing_tool = False
        self._hidden_pair = None
        self._last_zoom_axis = None
        self._build_ui()
        index = self.fit_type_combo.findData(normalize_poly_type(poly_type))
        if index >= 0:
            self.fit_type_combo.blockSignals(True)
            self.fit_type_combo.setCurrentIndex(index)
            self.fit_type_combo.blockSignals(False)
        self._arm_locked = arm is not None
        if arm:
            arm_index = self.arm_combo.findText(str(arm).upper())
            if arm_index >= 0:
                self.arm_combo.blockSignals(True)
                self.arm_combo.setCurrentIndex(arm_index)
                self.arm_combo.blockSignals(False)
        if flip_x:
            self.flip_x_box.blockSignals(True)
            self.flip_x_box.setChecked(True)
            self.flip_x_box.blockSignals(False)
        if reference_path:
            self.load_reference(reference_path, linelist_path)
        elif linelist_path:
            self.load_linelist(linelist_path)
        if lamp_path:
            self.load_lamp(lamp_path)
        if anchor_path:
            self.load_anchor_file(anchor_path)
        if config_path:
            self._apply_config_range(config_path)
        if wave_range:
            self.apply_wave_range(*wave_range)

    # ------------------------------------------------------------------ UI
    def _build_ui(self):
        QW = QtWidgets
        tb = QW.QToolBar('Files')
        tb.setMovable(False)
        self.addToolBar(tb)
        for text, slot in (('Load reference atlas / line list…', self.choose_reference),
                           ('Load lamp spectrum…', self.choose_lamp),
                           ('Load identified lines…', self.choose_anchors),
                           ('Save identified lines…', self.save_dialog),
                           ('Save diagnostic figure…', self.save_plot_dialog)):
            tb.addAction(text).triggered.connect(slot)

        central = QW.QWidget()
        self.setCentralWidget(central)
        layout = QW.QHBoxLayout(central)

        left = QW.QVBoxLayout()
        self.fig = Figure(figsize=(10, 7.6))
        gs = self.fig.add_gridspec(2, 1, height_ratios=[3.2, 1.35], hspace=0.0)
        self.ax = self.fig.add_subplot(gs[0])
        self.resid_ax = self.fig.add_subplot(gs[1])
        self.canvas = FigureCanvasQTAgg(self.fig)
        self.fig_res = Figure(figsize=(10, 7.6))
        self.res_ax = self.fig_res.add_subplot(111)
        self.res_canvas = FigureCanvasQTAgg(self.fig_res)
        self.tabs = QW.QTabWidget()
        self.tabs.addTab(self.canvas, 'Spectrum and residuals')
        res_tab = QW.QWidget()
        res_layout = QW.QVBoxLayout(res_tab)
        res_layout.setContentsMargins(0, 0, 0, 0)
        self.res_toolbar = NavigationToolbar2QT(self.res_canvas, res_tab)
        res_layout.addWidget(self.res_toolbar)
        res_layout.addWidget(self.res_canvas, 1)
        self.tabs.addTab(res_tab, 'Resolution R(λ)')
        self.toolbar = SpectrumToolbar(self.canvas, self, self)
        left.addWidget(self.toolbar)
        left.addWidget(self.tabs, 1)
        self.canvas.mpl_connect('button_press_event', self.on_click)
        self.canvas.mpl_connect('motion_notify_event', self.on_move)
        self.canvas.mpl_connect('scroll_event', self.on_scroll)
        self.canvas.mpl_connect('resize_event', self._on_resize)
        self.canvas.setFocusPolicy(Qt.StrongFocus)

        opts = QW.QHBoxLayout()
        self.flip_x_box = QW.QCheckBox('Flip x')
        self.flip_x_box.setToolTip('Reverse the displayed pixel axis.  Use it when wavelength '
                                   'decreases with pixel so the measured spectrum and the '
                                   'reference read in the same direction; auto identify sets '
                                   'it from the fitted dispersion.')
        self.flip_x_box.stateChanged.connect(lambda *_: self._flip_x_changed())
        self.lo_spin = QW.QDoubleSpinBox()
        self.hi_spin = QW.QDoubleSpinBox()
        for sp in (self.lo_spin, self.hi_spin):
            sp.setRange(0.0, 100000.0)
            sp.setDecimals(1)
            sp.setSingleStep(50.0)
            sp.valueChanged.connect(lambda *_: self.redraw())
        opts.addWidget(self.flip_x_box)
        opts.addWidget(QW.QLabel('Camera arm'))
        self.arm_combo = QW.QComboBox()
        self.arm_combo.addItems(['R', 'B'])
        self.arm_combo.setToolTip('MSTG 1.93m camera arm, used for the default anchor '
                                  'filename')
        opts.addWidget(self.arm_combo)
        opts.addStretch(1)
        left.addLayout(opts)
        zoom = QW.QHBoxLayout()
        self.zoom_out_btn = QW.QPushButton('−')
        self.zoom_out_btn.setFixedWidth(32)
        self.zoom_out_btn.setToolTip('Zoom out the active panel')
        self.zoom_out_btn.clicked.connect(lambda: self._zoom_button(-1))
        self.zoom_in_btn = QW.QPushButton('+')
        self.zoom_in_btn.setFixedWidth(32)
        self.zoom_in_btn.setToolTip('Zoom in the active panel')
        self.zoom_in_btn.clicked.connect(lambda: self._zoom_button(1))
        self.zoom_slider = QW.QSlider(Qt.Horizontal)
        self.zoom_slider.setRange(1, 50)
        self.zoom_slider.setValue(10)
        self.zoom_slider.setFixedWidth(150)
        self.zoom_slider.setToolTip('Zoom step used by the +/− buttons and the scroll wheel')
        self.zoom_label = QW.QLabel('10%')
        self.zoom_label.setFixedWidth(38)
        self.zoom_slider.valueChanged.connect(self._zoom_slider_changed)
        for w in (self.zoom_out_btn, self.zoom_in_btn, QW.QLabel('Zoom step'),
                  self.zoom_slider, self.zoom_label):
            zoom.addWidget(w)
        zoom.addStretch(1)
        left.addLayout(zoom)
        ov = QW.QHBoxLayout()
        self.ov_box = QW.QCheckBox('Show reference')
        self.ov_box.setChecked(True)
        self.ov_box.setToolTip('Show the reference spectrum or line list below the measured '
                               'spectrum.')
        self.ov_box.stateChanged.connect(lambda *_: self.redraw())
        ov.addWidget(self.ov_box)
        ov.addStretch(1)
        left.addLayout(ov)
        rng = QW.QHBoxLayout()
        rng.addWidget(QW.QLabel('Reference wavelength range: low'))
        rng.addWidget(self.lo_spin)
        rng.addWidget(QW.QLabel('upper'))
        rng.addWidget(self.hi_spin)
        rng.addStretch(1)
        left.addLayout(rng)
        layout.addLayout(left, 1)

        right = QW.QVBoxLayout()
        self.pick_label = QW.QLabel()
        self.pick_label.setWordWrap(True)
        right.addWidget(self.pick_label)
        pick = QW.QHBoxLayout()
        cancel = QW.QPushButton('Cancel selection (Esc)')
        cancel.setToolTip('Clear the pending peak/reference selection')
        cancel.clicked.connect(self.cancel_selection)
        pick.addWidget(cancel)
        pick.addWidget(QW.QLabel('peak → reference adds pair'))
        pick.addStretch(1)
        right.addLayout(pick)

        right.addWidget(QW.QLabel('Identified lines (pixel is 0-based)'))
        self.table = QW.QTableWidget(0, 5)
        self.table.setHorizontalHeaderLabels(
            ['Type', 'Pixel', 'λ (Å)', 'Resid. (Å)', 'Ion'])
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.setSelectionBehavior(QW.QAbstractItemView.SelectRows)
        self.table.setEditTriggers(QW.QAbstractItemView.NoEditTriggers)
        self.table.verticalHeader().setDefaultSectionSize(22)
        self.table.setMinimumHeight(380)
        right.addWidget(self.table, 3)
        btns = QW.QHBoxLayout()
        for text, slot in (('Remove selected', self.remove_selected),
                           ('Clear manual', self.clear_manual),
                           ('Clear auto', self.clear_auto)):
            b = QW.QPushButton(text)
            b.clicked.connect(slot)
            btns.addWidget(b)
        right.addLayout(btns)

        box = QW.QGroupBox('Auto identification')
        form = QW.QFormLayout(box)
        self.snr_spin = QW.QDoubleSpinBox()
        self.snr_spin.setRange(MIN_DETECT_SNR, 1000.0)
        self.snr_spin.setValue(50.0)
        self.snr_spin.setToolTip(
            'Only measured lamp peaks at or above this S/N are used for automatic '
            'matching and anchor selection.')
        self.snr_spin.valueChanged.connect(lambda *_: self.redraw(keep_view=True))
        self.order_spin = QW.QSpinBox()
        self.order_spin.setRange(1, 7)
        self.order_spin.setValue(5)
        self.order_spin.valueChanged.connect(lambda *_: self.update_model())
        self.tol_spin = QW.QDoubleSpinBox()
        self.tol_spin.setRange(0.1, 100.0)
        self.tol_spin.setValue(3.0)
        self.tol_spin.setSuffix(' Å')
        form.addRow('Peak SNR threshold', self.snr_spin)
        form.addRow('Matching polynomial order', self.order_spin)
        form.addRow('Initial match tolerance', self.tol_spin)
        run = QW.QPushButton('Auto identify')
        run.setToolTip('Match the strong lamp peaks against the reference (manual pairs, '
                       'if any, seed the solution).')
        run.clicked.connect(self.auto_identify)
        form.addRow(run)
        right.addWidget(box)

        fitbox = QW.QGroupBox('Residual fit')
        fitform = QW.QFormLayout(fitbox)
        self.fit_type_combo = QW.QComboBox()
        self.fit_type_combo.addItem('Legendre', 'legendre')
        self.fit_type_combo.addItem('Chebyshev', 'chebyshev')
        self.fit_type_combo.addItem('Ordinary polynomial', 'polynomial')
        self.fit_type_combo.currentIndexChanged.connect(lambda *_: self.update_model())
        self.fit_order_spin = QW.QSpinBox()
        self.fit_order_spin.setRange(1, 9)
        self.fit_order_spin.setValue(5)
        self.fit_order_spin.setToolTip('Degree of the polynomial fitted to the identified '
                                       '(pixel, wavelength) anchors')
        self.fit_order_spin.valueChanged.connect(lambda *_: self.update_model())
        self.clip_mode_combo = QW.QComboBox()
        self.clip_mode_combo.addItem('Sigma clipping (nσ)', 'sigma')
        self.clip_mode_combo.addItem('Absolute (Å)', 'absolute')
        self.clip_mode_combo.currentIndexChanged.connect(self._clip_mode_changed)
        self.fit_clip_spin = QW.QDoubleSpinBox()
        self.fit_clip_spin.valueChanged.connect(lambda *_: self.update_model())
        self._configure_clip_control()
        refit = QW.QPushButton('Refit residuals')
        refit.clicked.connect(lambda *_: self.update_model())
        fitform.addRow('Basis', self.fit_type_combo)
        fitform.addRow('Polynomial degree', self.fit_order_spin)
        fitform.addRow('Clip mode', self.clip_mode_combo)
        fitform.addRow('Residual clip', self.fit_clip_spin)
        fitform.addRow(refit)
        self.resolution_label = QW.QLabel('Resolution R: identify at least two lines')
        self.resolution_label.setWordWrap(True)
        self.resolution_label.setToolTip('Estimated from the Gaussian FWHM of clean arc lines: '
                                         'R = lambda / FWHM.')
        fitform.addRow(self.resolution_label)
        self.fit_formula = QW.QLabel('')
        self.fit_formula.setWordWrap(True)
        self.fit_formula.setTextInteractionFlags(Qt.TextSelectableByMouse)
        self.fit_formula.setStyleSheet('font-family: monospace; font-size: 8pt;')
        fitform.addRow(self.fit_formula)
        right.addWidget(fitbox)

        self.fit_label = QW.QLabel('')
        self.fit_label.setWordWrap(True)
        right.addWidget(self.fit_label)
        help_txt = QW.QLabel(
            'Click a measured peak, then a reference peak/line; the pair is added '
            'automatically. Esc cancels the last selection. Pan/zoom continuously with the '
            'toolbar, then press Select to pick lines.')
        help_txt.setWordWrap(True)
        right.addWidget(help_txt)
        panel = QW.QWidget()
        panel.setLayout(right)
        scroll = QW.QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(panel)
        scroll.setFixedWidth(460)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        layout.addWidget(scroll)
        self._esc_shortcut = QW.QShortcut(Qt.Key_Escape, self)
        self._esc_shortcut.activated.connect(self.cancel_last_selection)
        self.update_pick_label()
        self.statusBar().showMessage(
            'Load a reference atlas (or line list) and a lamp spectrum to begin.')

    # ------------------------------------------------------------ file I/O
    def choose_reference(self):
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, 'Reference atlas spectrum or line list', '',
            'References (*.fits *.fit *.dat *.txt *.csv);;All files (*)')
        if path:
            self.load_reference(path)

    def choose_lamp(self):
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, 'Extracted 1-D lamp spectrum', '',
            'Spectra (*.fits *.fit *.txt *.csv *.dat);;All files (*)')
        if path:
            self.load_lamp(path)

    def choose_anchors(self):
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, 'Identified lines', '', 'CSV/text (*.csv *.txt);;All files (*)')
        if path:
            self.load_anchor_file(path)

    def load_reference(self, path: str, linelist_path: Optional[str] = None):
        try:
            ref = load_reference(path, linelist_path)
        except Exception as exc:
            QtWidgets.QMessageBox.warning(self, 'Reference', str(exc))
            return
        self._set_reference(ref, path)
        kind = 'atlas spectrum' if self.ref_is_spectrum else 'line list'
        self.statusBar().showMessage(
            f'{len(ref)} reference entries ({kind}) loaded from {path}')

    def load_linelist(self, path: str):
        try:
            wave, labels, strength = load_linelist(path)
        except Exception as exc:
            QtWidgets.QMessageBox.warning(self, 'Reference lines', str(exc))
            return
        ref = Reference(wave, None, labels, source=path)
        self.ref_strengths = np.asarray(strength, dtype=float)
        self._set_reference(ref, path)
        self.statusBar().showMessage(f'{len(wave)} reference lines loaded from {path}')

    def _set_reference(self, ref: Reference, path: str):
        self.ref = ref
        self.ref_path = path
        self.ref_is_spectrum = len(ref) > 0 and not np.allclose(ref.flux, 1.0)
        if not self.ref_is_spectrum and len(self.ref_strengths) != len(ref):
            self.ref_strengths = np.ones(len(ref))
        self._catalog_cache = None
        if len(ref):
            lo, hi = float(np.nanmin(ref.wave)), float(np.nanmax(ref.wave))
            for sp, value in ((self.lo_spin, lo), (self.hi_spin, hi)):
                if sp.value() <= 0:
                    sp.blockSignals(True)
                    sp.setValue(value)
                    sp.blockSignals(False)
        self.redraw()

    def _apply_config_range(self, config_path: str):
        try:
            from specproc.config.config_manager import (
                ConfigManager, resolve_telescope_profile)
            resolved = resolve_telescope_profile(config_path) or config_path
            cfg = ConfigManager(str(resolved))
            if not self.telescope_name:
                self.telescope_name = cfg.get('telescope', 'name', '') or ''
            w0 = cfg.get_float('slit', 'wave_start', 0.0)
            w1 = cfg.get_float('slit', 'wave_end', 0.0)
            if w0 > 0 and w1 > w0:
                self.apply_wave_range(w0, w1)
        except Exception as exc:  # pragma: no cover - configuration is optional
            self.statusBar().showMessage(f'Could not read wave range from config: {exc}')

    def apply_wave_range(self, w0: float, w1: float):
        """Set the reference wavelength range used for matching."""
        w0, w1 = sorted(float(v) for v in (w0, w1))
        for sp, value in ((self.lo_spin, w0), (self.hi_spin, w1)):
            sp.blockSignals(True)
            sp.setValue(value)
            sp.blockSignals(False)
        self.update_model()

    def load_lamp(self, path: str):
        try:
            self.flux = load_lamp_spectrum(path, self.dispersion_axis)
        except Exception as exc:
            QtWidgets.QMessageBox.warning(self, 'Lamp spectrum', str(exc))
            return
        self.lamp_path = path
        tel_name, arm = _fits_header_meta(path)
        if tel_name:
            self.telescope_name = tel_name
        if arm and not self._arm_locked:
            idx = self.arm_combo.findText(arm)
            if idx >= 0:
                self.arm_combo.blockSignals(True)
                self.arm_combo.setCurrentIndex(idx)
                self.arm_combo.blockSignals(False)
        self.peaks = detect_peaks(self.flux, MIN_DETECT_SNR, min_width=0.5,
                                  max_width=10.0, distance=4, baseline_size=61)
        self.sel_peak = self.sel_ref = None
        self.update_pick_label()
        self._update_fit()
        self.redraw()
        self.statusBar().showMessage(
            f'{len(self.flux)} pixels, {len(self._visible_peaks())} peaks '
            f'(SNR ≥ {self.snr_spin.value():g}): {path}')

    def load_anchor_file(self, path: str):
        try:
            pairs = load_anchors(path)
        except Exception as exc:
            QtWidgets.QMessageBox.warning(self, 'Identified lines', str(exc))
            return
        self.manual = [(float(x), float(w)) for x, w in pairs]
        self.auto = []
        self.anchor_path = path
        self.update_model(keep_view=False, auto_flip=True)
        self.statusBar().showMessage(f'{len(self.manual)} lines loaded as manual pairs from {path}')

    def default_anchor_path(self):
        """Suggested anchor filename, using TELESCOP/INSTRUME when available."""
        path_text = ' '.join(str(p) for p in
                             (self.lamp_path, self.ref_path, self.anchor_path) if p).lower()
        is_mstg = (self.telescope_name.lower() == 'mstg193slit'
                   or 'mstg193slit' in path_text)
        if self.output_path:
            return self.output_path
        if is_mstg:
            return str(Path.cwd() / f'mstg193slit_anchors_{self.arm_combo.currentText()}.csv')
        if self.lamp_path:
            return str(Path.cwd() / (Path(self.lamp_path).stem + '_anchors.csv'))
        return str(Path.cwd() / 'anchors.csv')

    def save_dialog(self):
        pairs = self.all_pairs()
        if len(pairs) < 2:
            QtWidgets.QMessageBox.information(
                self, 'Save', 'Identify at least two lines before saving.')
            return
        default = self.default_anchor_path()
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, 'Save identified lines', default, 'CSV (*.csv);;All files (*)')
        if not path:
            return
        try:
            solution = self.build_solution()
        except ValueError as exc:
            QtWidgets.QMessageBox.warning(self, 'Save', str(exc))
            return
        self.redraw(keep_view=True)
        anchors = self._anchor_rows(solution)
        note = f'{len(self.manual)} manual + {len(self.auto)} automatic lines'
        if self.disp_fit is not None:
            note += (f'; {dispersion_type(self.disp_fit)} '
                     f'degree {self.disp_fit.degree()}')
        write_anchor_file(path, anchors, solution,
                          lamp=self.lamp_path or '(none)',
                          reference=self.ref_path or '(none)',
                          note=note,
                          created=str(np.datetime64('today')))
        plot_file = (str(Path(self.plot_path).with_suffix('.pdf'))
                     if self.plot_path else str(Path(path).with_suffix('.pdf')))
        self.solution = solution
        self.anchor_path = path
        self.plot_path = plot_file
        self._render_diagnostic(plot_file, solution, anchors)
        self.statusBar().showMessage(
            f'{len(anchors)} lines saved to {path}; diagnostic figure: {plot_file}')

    def save_plot_dialog(self):
        pairs = self.all_pairs()
        if len(pairs) < 2:
            QtWidgets.QMessageBox.information(
                self, 'Diagnostic figure', 'Identify at least two lines first.')
            return
        if self.plot_path:
            default = str(Path(self.plot_path).with_suffix('.pdf'))
        elif self.anchor_path:
            default = str(Path.cwd() / (Path(self.anchor_path).stem + '.pdf'))
        elif self.lamp_path:
            default = str(Path.cwd() / (Path(self.lamp_path).stem + '_anchors.pdf'))
        else:
            default = str(Path.cwd() / 'anchors_diagnostic.pdf')
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, 'Save diagnostic figure', default,
            'PDF (*.pdf);;All files (*)')
        if not path:
            return
        try:
            solution = self.build_solution()
        except ValueError as exc:
            QtWidgets.QMessageBox.warning(self, 'Diagnostic figure', str(exc))
            return
        self.redraw(keep_view=True)
        anchors = self._anchor_rows(solution)
        self.plot_path = path
        self._render_diagnostic(path, solution, anchors)
        self.statusBar().showMessage(f'Diagnostic figure saved to {path}')

    def _render_diagnostic(self, path: str, solution, anchors):
        # Save the exact canvases the user is looking at.  PDF gets one page
        # for the spectrum/residual tab and a second for R(lambda).
        path = Path(path)
        if path.suffix.lower() != '.pdf':
            path = path.with_suffix('.pdf')
        from matplotlib.backends.backend_pdf import PdfPages
        with PdfPages(path) as pdf:
            pdf.savefig(self.fig, facecolor='white')
            pdf.savefig(self.fig_res, facecolor='white')
        return path

    # ------------------------------------------------------------- state
    def min_tol(self):
        """Final match tolerance for the automatic anchor extension."""
        return 0.4

    def display_flux(self):
        """Lamp flux normalized to 0-1 with the baseline removed."""
        if len(self.flux) == 0:
            return self.flux
        base = np.percentile(self.flux, 1)
        top = max(np.nanmax(self.flux) - base, 1e-12)
        return np.clip((self.flux - base) / top, 0.0, None)

    def ref_mask(self):
        if self.ref is None:
            return np.zeros(0, dtype=bool)
        return (self.ref.wave >= self.lo_spin.value()) & (self.ref.wave <= self.hi_spin.value())

    def all_pairs(self):
        pairs = [(x, w) for x, w in self.manual]
        pairs += [(x, w) for x, w, _ in self.auto]
        return sorted(pairs)

    def base_wavelength(self, pixel):
        """Dispersion model used for display, picking and the residual fit."""
        if self.disp_fit is None:
            return None
        return self.disp_fit(np.asarray(pixel, dtype=float))

    def wavelength_to_pixel(self, wave):
        """Map reference wavelengths to the displayed pixel axis (NaN if outside)."""
        model = self.overlay_model()
        if model is None or len(model) < 2:
            return np.full(np.shape(wave), np.nan, dtype=float)
        model = np.asarray(model, dtype=float)
        if not (np.all(np.diff(model) > 0) or np.all(np.diff(model) < 0)):
            return np.full(np.shape(wave), np.nan, dtype=float)
        pixel = np.arange(len(model), dtype=float)
        order = np.argsort(model)
        return np.interp(np.asarray(wave, dtype=float), model[order], pixel[order],
                         left=np.nan, right=np.nan)

    def fit_type(self) -> str:
        """Currently selected fitting basis."""
        return normalize_poly_type(self.fit_type_combo.currentData() or 'legendre')

    def clip_mode(self) -> str:
        return str(self.clip_mode_combo.currentData() or 'sigma')

    def _configure_clip_control(self):
        mode = self.clip_mode()
        self.fit_clip_spin.blockSignals(True)
        if mode == 'absolute':
            self.fit_clip_spin.setRange(0.0, 5.0)
            self.fit_clip_spin.setDecimals(2)
            self.fit_clip_spin.setSingleStep(0.05)
            self.fit_clip_spin.setSuffix(' Å')
            self.fit_clip_spin.setValue(0.40)
            self.fit_clip_spin.setToolTip(
                'Absolute residual cut in Angstrom (default 0.40 Å); 0 disables clipping.')
        else:
            self.fit_clip_spin.setRange(0.0, 10.0)
            self.fit_clip_spin.setDecimals(1)
            self.fit_clip_spin.setSingleStep(0.5)
            self.fit_clip_spin.setSuffix(' σ')
            self.fit_clip_spin.setValue(3.0)
            self.fit_clip_spin.setToolTip(
                'Iterative RMS clipping: keep residuals within nσ of the current fit. '
                'nσ is measured against the RMS of the current inliers; 1σ is exactly '
                'the ±1σ RMS line in the residual plot.')
        self.fit_clip_spin.setSpecialValueText('off')
        self.fit_clip_spin.blockSignals(False)

    def _clip_mode_changed(self, *_):
        self._configure_clip_control()
        self.update_model()

    def _clip_threshold(self, residuals, keep):
        value = float(self.fit_clip_spin.value())
        if value <= 0:
            return None
        if self.clip_mode() == 'absolute':
            return value
        res = np.asarray(residuals, dtype=float)[keep]
        if not len(res):
            return None
        rms = float(np.sqrt(np.mean(res ** 2)))
        return value * rms

    def _build_dispersion_fit(self):
        """Fit ``pixel -> wavelength`` with the selected basis and clip outliers."""
        pairs = np.asarray(self.all_pairs(), dtype=float).reshape(-1, 2)
        if len(pairs) < 2 or len(self.flux) < 2:
            return None, np.empty(0), np.zeros(0, dtype=bool)
        x, y = pairs[:, 0], pairs[:, 1]
        keep = np.isfinite(x) & np.isfinite(y)
        if keep.sum() < 2:
            return None, np.empty(0), np.zeros(0, dtype=bool)

        degree = int(self.fit_order_spin.value())
        domain = (0.0, float(len(self.flux) - 1))
        fit = None
        for _ in range(10):
            n_keep = int(keep.sum())
            if n_keep < 2:
                break
            degree = max(1, min(degree, n_keep - 1))
            fit = fit_dispersion(x[keep], y[keep], degree, self.fit_type(),
                                 domain=domain)
            if fit is None:
                break
            residuals = y - fit(x)
            threshold = self._clip_threshold(residuals, keep)
            if threshold is None:
                break
            new_keep = np.abs(residuals) <= threshold
            # Always keep enough points for the requested degree.
            if new_keep.sum() < max(2, degree + 1):
                idx = np.argsort(np.abs(residuals))[:max(2, degree + 1)]
                new_keep = np.zeros(len(x), dtype=bool)
                new_keep[idx] = True
            if np.array_equal(new_keep, keep):
                break
            keep = new_keep
        if fit is None:
            return None, np.empty(0), np.zeros(0, dtype=bool)
        return fit, y - fit(x), keep

    def _update_fit(self):
        self.disp_fit, self.disp_residuals, self.disp_mask = self._build_dispersion_fit()
        self.resolution_info = None
        if self.disp_fit is not None and len(self.flux) >= 4:
            try:
                self.resolution_info = measure_resolution(
                    np.arange(len(self.flux), dtype=float), self.flux,
                    self._visible_peaks(), self.disp_fit)
            except Exception:
                self.resolution_info = None
        if hasattr(self, 'resolution_label'):
            if self.resolution_info:
                info = self.resolution_info
                self.resolution_label.setText(
                    f'R ≈ {info["resolution"]:,.0f}   '
                    f'(FWHM {info["fwhm_ang"]:.2f} Å, n={info["n"]}'
                    f', rejected {info.get("n_rejected", 0)})')
            else:
                self.resolution_label.setText(
                    'Resolution R: identify at least two clean lines')

    def _poly_coefficients(self):
        """Selected fit expressed as ordinary polynomial coefficients (descending)."""
        return dispersion_polyval_coef(self.disp_fit)

    def _fit_formula_text(self):
        if self.disp_fit is None:
            return 'Dispersion fit: identify at least two lines'
        kind = dispersion_type(self.disp_fit)
        n = max(len(self.flux) - 1, 1)
        if kind == 'polynomial':
            terms = ' '.join(f'{c:+.6g}·t^{k}'
                             for k, c in enumerate(self.disp_fit.coef))
        else:
            symbol = 'P' if kind == 'legendre' else 'T'
            terms = ' '.join(f'{c:+.6g}·{symbol}{k}(t)'
                             for k, c in enumerate(self.disp_fit.coef))
        return f'{kind} degree {self.disp_fit.degree()}\nt = 2x/{n} - 1\nλ(t) = {terms}'

    def _set_flip_x_from_endpoints(self, w0: float, w1: float):
        if not np.isfinite(w0) or not np.isfinite(w1) or w0 == w1:
            return
        flip = bool(w1 < w0)
        if self.flip_x_box.isChecked() == flip:
            return
        self.flip_x_box.blockSignals(True)
        self.flip_x_box.setChecked(flip)
        self.flip_x_box.blockSignals(False)

    def _auto_flip_from_solution(self, solution):
        if solution is not None and len(self.flux) >= 2:
            w0, w1 = solution.wavelength([0.0, float(len(self.flux) - 1)])
            self._set_flip_x_from_endpoints(float(w0), float(w1))

    def _auto_flip_from_fit(self):
        if self.disp_fit is not None and len(self.flux) >= 2:
            w0, w1 = self.disp_fit([0.0, float(len(self.flux) - 1)])
            self._set_flip_x_from_endpoints(float(w0), float(w1))

    def _flip_x_changed(self):
        if self._syncing:
            return
        self.redraw()

    def _zoom_slider_changed(self, value):
        self.zoom_label.setText(f'{int(value)}%')

    def _zoom_button(self, direction: int):
        step = max(1, int(self.zoom_slider.value())) / 100.0
        factor = 1.0 + step
        if direction > 0:
            factor = 1.0 / factor
        # The +/- buttons always act on the measured spectrum and use its
        # current view, so panning first does not send the zoom back to an
        # earlier position.
        self._last_zoom_axis = self.ax
        self._zoom_active(factor)

    def _zoom_active(self, factor: float):
        axis = self._last_zoom_axis
        if axis is None or (axis is self.ref_axis and self.ref_axis is None):
            axis = self.ax
        x0, x1 = axis.get_xlim()
        if not (np.isfinite(x0) and np.isfinite(x1)) or x0 == x1:
            return
        centre = 0.5 * (x0 + x1)
        half = 0.5 * abs(x1 - x0) * factor
        a, b = centre - half, centre + half
        if b - a <= 1e-9:
            return
        # Keep the current centre exactly.  Clamping to the detector bounds
        # here makes a panned view jump back when +/− is pressed at an edge.
        if x1 >= x0:
            axis.set_xlim(a, b)
        else:
            axis.set_xlim(b, a)
        if axis is self.ref_axis:
            self.label_reference_lines()
        self.canvas.draw_idle()

    def _on_xlim_changed(self, _ax):
        if self._syncing:
            return
        self._sync_residual_xlim()
        if not self._plot_in_wave:
            self._draw_wave_labels()
        self.label_reference_lines()
        self.canvas.draw_idle()

    def _sync_residual_xlim(self):
        if self.disp_fit is None or len(self.flux) < 2:
            return
        x_limits = self.ax.get_xlim()
        wave_limits = self._model_at(np.asarray(x_limits, dtype=float))
        if np.all(np.isfinite(wave_limits)):
            self.resid_ax.set_xlim(float(wave_limits[0]), float(wave_limits[1]))

    def _on_residual_xlim_changed(self, _ax):
        if self._syncing or self.disp_fit is None or len(self.flux) < 2:
            return
        pixels = np.arange(len(self.flux), dtype=float)
        waves = np.asarray(self._model_at(pixels), dtype=float)
        if not (np.all(np.diff(waves) > 0) or np.all(np.diff(waves) < 0)):
            return
        wave_limits = self.resid_ax.get_xlim()
        if (min(wave_limits) < min(waves[0], waves[-1])
                or max(wave_limits) > max(waves[0], waves[-1])):
            return
        if waves[0] > waves[-1]:
            waves = waves[::-1]
            pixels = pixels[::-1]
        x_limits = np.interp(wave_limits, waves, pixels)
        self._syncing = True
        try:
            self.ax.set_xlim(float(x_limits[0]), float(x_limits[1]))
        finally:
            self._syncing = False
        if not self._plot_in_wave:
            self._draw_wave_labels()
        self.label_reference_lines()
        self.canvas.draw_idle()

    def _on_resize(self, _event):
        if self.ref_axis is None:
            return
        pos = self.ax.get_position()
        self.ref_axis.set_position(
            [pos.x0, pos.y0, pos.width, pos.height * 0.47])
        self.label_reference_lines()
        self.canvas.draw_idle()

    def set_select_mode(self, active, update_toolbar: bool = True):
        """Switch between line selection and continuous pan/zoom."""
        self.select_mode = bool(active)
        if update_toolbar or not self.select_mode:
            action = self.toolbar.select_action
            action.blockSignals(True)
            action.setChecked(self.select_mode)
            action.blockSignals(False)
        if self.select_mode:
            mode = str(self.toolbar.mode).upper()
            self._changing_tool = True
            try:
                if 'PAN/ZOOM' in mode:
                    self.toolbar.pan()
                elif 'ZOOM' in mode:
                    self.toolbar.zoom()
                elif 'PAN' in mode:
                    self.toolbar.pan()
            finally:
                self._changing_tool = False
            self.statusBar().showMessage(
                'Select mode: click a lamp peak, then a reference peak/line.')
        elif update_toolbar:
            self.statusBar().showMessage(
                'Pan/zoom mode: drag or scroll to align; click Select to pick lines.')

    def catalog(self, wave_range=None, for_matching: bool = False):
        """Strong-line catalogue for matching: atlas peaks or the line list."""
        ref = self.ref
        if ref is None or not len(ref):
            return Reference(np.empty(0), None, [])
        if not self.ref_is_spectrum:
            cat = ref.subset(wave_range)
            return Reference(cat.wave, np.ones_like(cat.wave), list(cat.labels), ref.source)
        # Matching uses a fixed atlas catalogue so changing the display/anchor
        # S/N threshold cannot change the wavelength solution. Overlay plots
        # remain uncluttered by showing only as many strong lines as visible
        # observed peaks.
        n_ref = None if for_matching else max(3, len(self._visible_peaks()))
        key = (bool(for_matching), n_ref)
        if wave_range is not None and wave_range[1] > wave_range[0]:
            key = (round(float(wave_range[0]), 1),
                   round(float(wave_range[1]), 1), bool(for_matching), n_ref)
        if self._catalog_cache is not None and self._catalog_cache[0] == key:
            return self._catalog_cache[1]
        try:
            cat = build_catalog(ref, wave_range=wave_range, snr=5.0,
                                strength_frac=0.0, min_strength=0.0,
                                max_lines=400, n_peaks=n_ref)
        except Exception as exc:
            self.statusBar().showMessage(
                f'Reference catalogue: {exc}; falling back to all reference lines')
            cat = ref.subset(wave_range)
        self._catalog_cache = (key, cat)
        return cat

    def overlay_lines(self):
        """Reference lines drawn on the lamp spectrum (strong lines for atlas)."""
        if self.ref is None:
            return np.empty(0), np.empty(0)
        lo, hi = self.lo_spin.value(), self.hi_spin.value()
        wave_range = (lo, hi) if hi > lo else None
        if self.ref_is_spectrum:
            cat = self.catalog(wave_range)
            return cat.wave, cat.flux
        m = self.ref_mask()
        return self.ref.wave[m], self.ref_strengths[m]

    def _reference_marker_data(self):
        """Wavelengths and display heights of the strong reference peaks."""
        waves, strengths = self.overlay_lines()
        waves = np.asarray(waves, dtype=float)
        if not len(waves):
            return np.empty(0), np.empty(0)
        if self.ref_is_spectrum:
            flux = np.clip(np.nan_to_num(np.asarray(self.ref.flux, dtype=float)), 0, None)
            top = self._ref_top or max(float(np.nanmax(flux)) if len(flux) else 1.0, 1e-12)
            y = np.interp(waves, self.ref.wave, flux, left=0.0, right=0.0) / top
            return waves, y + REFERENCE_OFFSET
        strengths = np.asarray(strengths, dtype=float)
        top = max(float(np.nanmax(strengths)) if len(strengths) else 1.0, 1e-12)
        y = REFERENCE_OFFSET + 0.05 + 0.22 * np.clip(strengths / top, 0, 1)
        return waves, y

    def _reference_y_at(self, wave_value: float) -> float:
        if self.ref_is_spectrum:
            flux = np.clip(np.nan_to_num(np.asarray(self.ref.flux, dtype=float)), 0, None)
            top = self._ref_top or max(float(np.nanmax(flux)) if len(flux) else 1.0, 1e-12)
            return float(np.interp(wave_value, self.ref.wave, flux,
                                   left=0.0, right=0.0)) / top + REFERENCE_OFFSET
        waves, strengths = self.overlay_lines()
        if len(waves):
            i = int(np.argmin(np.abs(np.asarray(waves, dtype=float) - wave_value)))
            top = max(float(np.nanmax(strengths)) if len(strengths) else 1.0, 1e-12)
            return (REFERENCE_OFFSET + 0.05
                    + 0.22 * float(np.clip(strengths[i] / top, 0, 1)))
        return REFERENCE_OFFSET + 0.05

    def _model_at(self, pixel):
        """Evaluate the dispersion model at arbitrary pixels."""
        return self.base_wavelength(pixel)

    def overlay_model(self):
        """Pixel -> wavelength model on the detector grid."""
        n = len(self.flux)
        if n < 2:
            return None
        return self._model_at(np.arange(n, dtype=float))

    def update_model(self, keep_view: bool = True, auto_flip: bool = False):
        self.manual.sort()
        self._update_fit()
        if auto_flip:
            self._auto_flip_from_fit()
        self.refresh_table()
        self.redraw(keep_view=keep_view)

    def refresh_table(self):
        rows = [('manual', x, w, self._label_for(x, w)) for x, w in self.manual]
        rows += [('auto', x, w, label) for x, w, label in self.auto]
        rows.sort(key=lambda r: r[1])
        self.table.setRowCount(len(rows))
        for i, (kind, x, w, label) in enumerate(rows):
            res = ('%.3f' % (w - float(self.disp_fit(x)))
                   if self.disp_fit is not None else '')
            for j, text in enumerate((kind, '%.2f' % x, '%.4f' % w, res, label)):
                self.table.setItem(i, j, QtWidgets.QTableWidgetItem(str(text)))
        if self.disp_fit is not None:
            mask = self.disp_mask
            if mask.any():
                rms = float(np.sqrt(np.mean(self.disp_residuals[mask] ** 2)))
            else:
                rms = float(np.sqrt(np.mean(self.disp_residuals ** 2)))
            self.fit_label.setText(
                f'{len(self.manual)} manual + {len(self.auto)} auto lines | '
                f'{dispersion_type(self.disp_fit)} degree {self.disp_fit.degree()} '
                f'| rms = {rms:.3f} Å'
                f'  ({int(mask.sum())}/{len(mask)} points in fit)')
            self.fit_formula.setText(self._fit_formula_text())
            self.fit_formula.setToolTip(self._fit_formula_text())
        elif self.auto:
            self.fit_label.setText(f'{len(self.auto)} auto-identified lines')
            self.fit_formula.setText('Dispersion fit: identify at least two lines')
        else:
            self.fit_label.setText(
                'Manual pairs are optional. Auto identify can start from the low/upper '
                'wavelength range and tries both dispersion directions.')
            self.fit_formula.setText('Dispersion fit: identify at least two lines')

    def _nearest_label(self, w: float) -> str:
        ref = self.ref
        if ref is None or not len(ref):
            return ''
        j = int(np.argmin(np.abs(ref.wave - w)))
        if abs(float(ref.wave[j]) - w) <= 2.0 and j < len(ref.labels):
            return str(ref.labels[j])
        return ''

    def _label_for(self, x: float, w: float) -> str:
        for ax, _aw, label in self.auto:
            if abs(ax - x) <= 0.5 and label:
                return label
        return self._nearest_label(w)

    def build_solution(self):
        """Build a Solution (with residuals) from the current manual + auto lines.

        The displayed fit uses the selected Legendre/Chebyshev basis; it is
        converted to ordinary
        polynomial coefficients here because :class:`Solution` and the slit
        pipeline evaluate the final model with ``numpy.polyval``.
        """
        pairs = np.array(self.all_pairs(), dtype=float).reshape(-1, 2)
        if len(pairs) < 2:
            raise ValueError('At least two identified lines are required')
        self._hidden_pair = None
        if self.disp_fit is None:
            self._update_fit()
        coef = self._poly_coefficients()
        mask = (self.disp_mask if len(self.disp_mask) == len(pairs)
                and self.disp_mask.any() else None)
        if coef is not None and np.all(np.isfinite(coef)):
            return solution_from_pairs(pairs, len(self.flux), peaks=self.peaks,
                                       catalog=self.ref, coef=coef,
                                       rms_mask=mask)
        order = max(1, min(self.order_spin.value(), len(pairs) - 1))
        return solution_from_pairs(pairs, len(self.flux), order=order,
                                   peaks=self.peaks, catalog=self.ref,
                                   rms_mask=mask)

    def _anchor_rows(self, solution):
        peak_px = np.array([p.pixel for p in self.peaks], dtype=float) if self.peaks else np.empty(0)
        peak_snr = np.array([p.snr for p in self.peaks], dtype=float) if self.peaks else np.empty(0)
        rows = []
        for (x, w), res in zip(solution.pairs, solution.residuals):
            snr = float(peak_snr[int(np.argmin(np.abs(peak_px - x)))]) if len(peak_px) else 0.0
            rows.append((snr, float(x), float(w), float(res), self._label_for(float(x), float(w))))
        return rows

    # -------------------------------------------------- picking and editing
    def update_pick_label(self):
        peak = 'none' if self.sel_peak is None else 'pixel %.2f' % self.sel_peak
        ref = 'none' if self.sel_ref is None else '%.4f Å' % self.sel_ref
        self.pick_label.setText(f'Selected lamp peak: {peak}\nSelected reference line: {ref}')

    def cancel_selection(self):
        self.sel_peak = self.sel_ref = None
        self.update_pick_label()
        self.redraw(keep_view=True)

    def cancel_last_selection(self):
        """Esc: undo the most recent pending peak/reference selection."""
        if self.sel_ref is not None:
            self.sel_ref = None
        elif self.sel_peak is not None:
            self.sel_peak = None
        else:
            return
        if self.sel_peak is None and self.sel_ref is None:
            self.update_pick_label()
        self.redraw(keep_view=True)

    def add_pair(self):
        if self.sel_peak is None or self.sel_ref is None:
            self.statusBar().showMessage('Select both a lamp peak and a reference line first.')
            return
        x, w = self.sel_peak, self.sel_ref
        self.manual = [m for m in self.manual if abs(m[0] - x) > 1.0 and abs(m[1] - w) > 1e-3]
        self.auto = [a for a in self.auto
                     if abs(a[0] - x) > 1.0 and abs(a[1] - w) > 1e-3]
        self.manual.append((x, w))
        self._hidden_pair = (x, w)
        # Keep the green selection markers visible until the next click; this
        # newest pair's anchor line stays hidden until the final fit.
        self.update_pick_label()
        self.update_model()

    def remove_selected(self):
        for r in sorted({i.row() for i in self.table.selectedIndexes()}):
            kind = self.table.item(r, 0).text()
            key = (float(self.table.item(r, 1).text()), float(self.table.item(r, 2).text()))
            lst = self.manual if kind == 'manual' else self.auto
            lst[:] = [p for p in lst if abs(p[0] - key[0]) > 0.01 or abs(p[1] - key[1]) > 1e-4]
        if self._hidden_pair not in self.manual:
            self._hidden_pair = None
        self.update_model()

    def clear_manual(self):
        self.manual = []
        self._hidden_pair = None
        self.update_model()

    def clear_auto(self):
        self.auto = []
        self.update_model()

    def auto_identify(self):
        if self.ref is None or not len(self.ref):
            self.statusBar().showMessage('Load a reference atlas or line list first.')
            return
        if not len(self.flux):
            self.statusBar().showMessage('Load a lamp spectrum first.')
            return
        lo, hi = self.lo_spin.value(), self.hi_spin.value()
        wave_range = (lo, hi) if hi > lo > 0 else None
        catalog = self.catalog(wave_range, for_matching=True)
        if len(catalog) < 4:
            self.statusBar().showMessage(
                'The reference catalogue has too few lines in this range.')
            return
        visible_peaks = self._visible_peaks()
        peaks = list(visible_peaks)
        if len(peaks) < 7 and len(self.manual) < 2:
            self.auto = []
            self.solution = None
            self.update_model()
            self.statusBar().showMessage(
                f'Only {len(peaks)} peaks meet S/N ≥ {self.snr_spin.value():g}; '
                'at least 7 are required for the fifth-order matching model.')
            return
        QtWidgets.QApplication.setOverrideCursor(Qt.WaitCursor)
        QtWidgets.QApplication.processEvents()
        try:
            if len(self.manual) >= 2:
                # Anchors first, then extend them to the full catalogue; this
                # is the same matching logic the slit pipeline applies to
                # [slit] identification_file.
                solution = solve_from_anchors(
                    peaks, catalog, self.manual, len(self.flux),
                    order=self.order_spin.value(), match_tol=self.tol_spin.value(),
                    min_tol=self.min_tol(), poly_type=self.fit_type())
                if solution is None:
                    solution = solution_from_pairs(
                        self.manual, len(self.flux),
                        order=self.order_spin.value(), peaks=peaks, catalog=catalog)
            else:
                solution = solve_auto(
                    peaks, catalog, len(self.flux), wave_range=wave_range,
                    order=self.order_spin.value(), match_tol=self.tol_spin.value(),
                    min_tol=self.min_tol())
        finally:
            QtWidgets.QApplication.restoreOverrideCursor()
        if solution is None:
            self.auto = []
            self.solution = None
            self.update_model()
            self.statusBar().showMessage(
                'No dispersion solution found. Check the wavelength range, SNR threshold, '
                'reference lines and initial match tolerance.')
            return
        labels = [str(catalog.labels[ci]) if 0 <= ci < len(catalog.labels) else ''
                  for ci in solution.catalog_indices]
        auto = [(float(x), float(w), lab) for (x, w), lab in zip(solution.pairs, labels)]
        if self.manual:
            manual_px = np.array([m[0] for m in self.manual], dtype=float)
            auto = [a for a in auto if float(np.min(np.abs(manual_px - a[0]))) > 1.0]
        self.auto = auto
        self.solution = solution
        self._hidden_pair = None
        self._auto_flip_from_solution(solution)
        self.update_model(keep_view=False)
        self.statusBar().showMessage(
            f'Auto identify: {len(self.auto)} lines matched with {len(self.manual)} manual '
            f'({len(peaks)} peaks at S/N ≥ {self.snr_spin.value():g}); '
            f'rms {solution.rms:.3f} Å'
            + ('; pixel axis flipped' if self.flip_x_box.isChecked() else ''))

    # ---------------------------------------------------------- plotting
    def reset_view(self):
        """Return both panels to the full detector view (toolbar Home)."""
        self.redraw()

    def redraw(self, keep_view: bool = False):
        view_top = (self.ax.get_xlim(), self.ax.get_ylim()) if keep_view and len(self.flux) else None
        has_ref = self.ref is not None and len(self.ref) > 0
        draw_ref = has_ref and self.ov_box.isChecked()
        self._syncing = True
        if self._ax_xlim_cid is not None:
            self.ax.callbacks.disconnect(self._ax_xlim_cid)
        if self._resid_xlim_cid is not None:
            self.resid_ax.callbacks.disconnect(self._resid_xlim_cid)
        # Remove every auxiliary reference axis, including any that may have
        # been left behind by an earlier draw.  Only the main spectrum and
        # residual axes belong to this figure permanently.
        if self.ref_axis is not None and self._last_zoom_axis is self.ref_axis:
            self._last_zoom_axis = None
        for extra in list(self.fig.axes):
            if extra is self.ax or extra is self.resid_ax:
                continue
            try:
                extra.remove()
            except Exception:
                pass
        self.ref_axis = None
        self.ax.clear()
        self.resid_ax.clear()
        self._ax_xlim_cid = self.ax.callbacks.connect(
            'xlim_changed', self._on_xlim_changed)
        self._resid_xlim_cid = self.resid_ax.callbacks.connect(
            'xlim_changed', self._on_residual_xlim_changed)
        self.res_ax.clear()
        self._ref_texts = []
        self._wave_texts = []
        self._ref_top = None
        thr = self.snr_spin.value()
        x = np.arange(len(self.flux), dtype=float)
        self._plot_in_wave = False
        plot_x = x
        measured_offset = MEASURED_OFFSET if draw_ref else 0.0
        if len(self.flux):
            disp = self.display_flux()
            self.ax.plot(plot_x, disp + measured_offset, color=LAMP_COLOR, lw=0.65,
                         label='measured lamp', zorder=2)
            self.ax.set_yscale('linear')
            self.ax.xaxis.set_ticks_position('top')
            self.ax.tick_params(axis='x', top=True, labeltop=True,
                                bottom=False, labelbottom=False)
            self.ax.xaxis.set_label_position('top')
            self.ax.set_xlabel('')
            visible = self._visible_peaks()
            if visible:
                px = np.array([p.pixel for p in visible])
                ypeak = np.array([disp[max(int(round(v)) - 2, 0):int(round(v)) + 3].max()
                                  for v in px]) + 0.015 + measured_offset
                peak_x = px
                self.ax.plot(peak_x, ypeak, 'v', color=PEAK_COLOR, ms=5.0,
                             ls='none', alpha=0.95, zorder=9,
                             label=f'peaks (SNR ≥ {thr:g})')
                if self.sel_peak is not None:
                    j = int(np.argmin(np.abs(px - self.sel_peak)))
                    self.ax.plot([peak_x[j]], [ypeak[j]], 'v', color=SELECT_MARK_COLOR,
                                 ms=8.0, ls='none', zorder=10)
                    self.ax.text(peak_x[j], ypeak[j] + 0.035, f'{px[j]:.1f}',
                                 ha='center', va='bottom', fontsize=8,
                                 color=SELECT_MARK_COLOR, zorder=11,
                                 bbox=dict(fc='white', ec='none', alpha=0.8, pad=0.3))
            if self._hidden_pair is None:
                manual_draw = list(self.manual)
            else:
                manual_draw = [pair for pair in self.manual
                               if pair != self._hidden_pair]
            anchor_lines = ([(x, w, 'manual') for x, w in manual_draw]
                            + [(x, w, 'auto') for x, w, _ in self.auto])
            for k, (rx, rw, kind) in enumerate(sorted(anchor_lines, key=lambda p: p[0])):
                colour = MANUAL_COLOR if kind == 'manual' else AUTO_COLOR
                ypos = float(np.interp(rx, x, disp)) + measured_offset
                # A short vertical tick just above the measured peak; the
                # wavelength label sits above the tick.
                tick_x = rx
                lift = 0.065 + (0.035 if k % 2 else 0.0)
                tick_top = ypos + lift + 0.105
                self.ax.plot([tick_x, tick_x], [ypos + lift, tick_top],
                             color=colour, lw=1.15, solid_capstyle='butt',
                             zorder=9)
                self.ax.text(tick_x, tick_top + 0.020, f'{rw:.4f}',
                             ha='center', va='bottom', rotation=90,
                             color=colour, fontsize=6.2, zorder=10,
                             bbox=dict(fc='white', ec='none', alpha=0.75, pad=0.20))
            if draw_ref:
                self._draw_reference(x)
        self.ax.set_ylabel('Normalized flux' + (' + offset' if draw_ref else ''))
        title = (Path(self.lamp_path).name if self.lamp_path else 'Lamp spectrum')
        self.ax.set_title(title)
        self.ax.set_xlabel('')

        if len(self.flux):
            if view_top:
                self.ax.set_xlim(*view_top[0])
                self.ax.set_ylim(*view_top[1])
            else:
                lo, hi = -0.5, len(self.flux) - 0.5
                if self.flip_x_box.isChecked():
                    lo, hi = hi, lo
                self.ax.set_xlim(lo, hi)
                top = MEASURED_OFFSET + 1.50 if draw_ref else 1.50
                self.ax.set_ylim(-0.05, top)
            if draw_ref:
                self._draw_reference(plot_x)
                if self.disp_fit is not None and anchor_lines:
                    for rx, rw, _kind in anchor_lines:
                        if rx < 0 or rx >= len(disp):
                            continue
                        ref_y = self._reference_y_at(rw)
                        lamp_y = float(np.interp(rx, x, disp)) + measured_offset
                        self.ax.vlines(
                            rx, ref_y, lamp_y, color='black', linewidth=0.55,
                            linestyles=(0, (2, 2)), alpha=0.75, zorder=5)
        self._draw_residuals()
        self._draw_resolution()
        self._sync_residual_xlim()
        self._syncing = False
        self._view_in_wave = self._plot_in_wave
        self.label_reference_lines()
        self._draw_top_legend()
        self.ax.xaxis.set_ticks_position('top')
        self.ax.tick_params(axis='x', top=True, labeltop=True,
                            bottom=False, labelbottom=False)
        self.resid_ax.xaxis.set_ticks_position('bottom')
        self.resid_ax.tick_params(axis='x', top=False, labeltop=False,
                                  bottom=True, labelbottom=True)
        self.canvas.draw_idle()
        self.res_canvas.draw_idle()

    def _draw_top_legend(self):
        handles, labels = self.ax.get_legend_handles_labels()
        if self.ref_axis is not None:
            h2, l2 = self.ref_axis.get_legend_handles_labels()
            handles += h2
            labels += l2
        if not handles:
            return
        pairs = list(dict(zip(labels, handles)).items())
        self.ax.legend([handle for _, handle in pairs],
                       [label for label, _ in pairs],
                       loc='upper right', bbox_to_anchor=(0.99, 0.99),
                       fontsize=5.5, framealpha=0.88, borderpad=0.25,
                       handlelength=1.0, handletextpad=0.35,
                       borderaxespad=0.0)

    def _draw_wave_labels(self):
        """Manual bottom wavelength tick labels for the fitted case."""
        for t in self._wave_texts:
            t.remove()
        self._wave_texts = []
        if self.disp_fit is None:
            return
        if not self.ov_box.isChecked() or len(self.flux) < 2:
            return
        model = self.overlay_model()
        if model is None or len(model) < 2:
            return
        if not (np.all(np.diff(model) > 0) or np.all(np.diff(model) < 0)):
            return
        x0, x1 = sorted(self.ax.get_xlim())
        x0 = max(x0, 0.0)
        x1 = min(x1, float(len(self.flux) - 1))
        if x1 <= x0:
            return
        pix_ticks = np.linspace(x0, x1, 7)
        waves = self._model_at(pix_ticks)
        if waves is None:
            return
        for px, wv in zip(pix_ticks, waves):
            if not np.isfinite(wv):
                continue
            self._wave_texts.append(
                self.ax.text(px, -0.028, f'{wv:.0f}', transform=self.ax.get_xaxis_transform(),
                             ha='center', va='top', color=REF_COLOR, fontsize=7,
                             clip_on=False, zorder=8))
        self._wave_texts.append(
            self.ax.text(0.5, -0.058, 'Wavelength (Å)',
                         transform=self.ax.transAxes, ha='center', va='top',
                         color=REF_COLOR, fontsize=8, clip_on=False, zorder=8))

    def _draw_reference(self, pixel):
        """Reference below the measured spectrum, before/after a fit."""
        if self.disp_fit is None:
            self._draw_reference_independent()
            return
        model = self.overlay_model()
        if model is None or len(model) < 2:
            return
        if not (np.all(np.diff(model) > 0) or np.all(np.diff(model) < 0)):
            return
        if self.ref_is_spectrum:
            wave = np.asarray(self.ref.wave, dtype=float)
            flux = np.clip(np.nan_to_num(np.asarray(self.ref.flux, dtype=float)), 0, None)
            lo, hi = float(np.min(model)), float(np.max(model))
            sel = (wave >= lo) & (wave <= hi) & np.isfinite(flux)
            if sel.sum() < 2:
                return
            order = np.argsort(model)
            px = np.interp(wave[sel], model[order], pixel[order])
            fl = flux[sel]
            ipx = np.full(len(px), -1, dtype=int)
            finite = np.isfinite(px)
            ipx[finite] = np.round(px[finite]).astype(int)
            good = finite & (ipx >= 0) & (ipx < len(pixel))
            prof = np.zeros(len(pixel), dtype=float)
            np.maximum.at(prof, ipx[good], fl[good])
            top = max(float(np.nanmax(prof)), 1e-12)
            self._ref_top = top
            self.ax.plot(pixel, prof / top + REFERENCE_OFFSET,
                         color=REF_COLOR, lw=0.65, alpha=0.9, zorder=1,
                         label='reference atlas')
            self._draw_reference_markers_aligned()
            return
        lines_wave, lines_strength = self.overlay_lines()
        if not len(lines_wave):
            return
        m = (lines_wave >= model.min()) & (lines_wave <= model.max())
        if not m.any():
            return
        order = np.argsort(model)
        px = np.interp(lines_wave[m], model[order], pixel[order])
        strength = np.asarray(lines_strength[m], dtype=float)
        top = max(float(np.nanmax(strength)) if len(strength) else 1.0, 1e-12)
        heights = 0.05 + 0.22 * np.clip(strength / top, 0, 1)
        self.ax.vlines(px, REFERENCE_OFFSET, REFERENCE_OFFSET + heights,
                       color=REF_COLOR, lw=0.65, alpha=0.8, zorder=1,
                       label='reference lines')
        self._draw_reference_markers_aligned()

    def _draw_reference_markers_aligned(self):
        if not self.ref_is_spectrum:
            if self.sel_ref is not None:
                px_sel = self.wavelength_to_pixel(self.sel_ref)
                if np.isfinite(px_sel):
                    y_sel = self._reference_y_at(self.sel_ref)
                    self.ax.vlines([px_sel], REFERENCE_OFFSET, y_sel + 0.02,
                                   color=SELECT_MARK_COLOR, lw=2.2, zorder=6)
                    self.ax.text(px_sel, y_sel + 0.045, f'{self.sel_ref:.2f}',
                                 ha='center', va='bottom', fontsize=8,
                                 color=SELECT_MARK_COLOR, zorder=7,
                                 bbox=dict(fc='white', ec='none', alpha=0.8, pad=0.3))
            return
        waves, ys = self._reference_marker_data()
        if not len(waves):
            return
        px = self.wavelength_to_pixel(waves)
        valid = np.isfinite(px)
        if valid.any():
            self.ax.plot(px[valid], (ys + 0.015)[valid], 'v',
                         color=REF_MARK_COLOR, ms=3.8, ls='none', alpha=0.9,
                         zorder=4, label='reference peaks')
        if self.sel_ref is not None:
            px_sel = self.wavelength_to_pixel(self.sel_ref)
            if np.isfinite(px_sel):
                y_sel = self._reference_y_at(self.sel_ref) + 0.015
                self.ax.plot([px_sel], [y_sel], 'v', color=SELECT_MARK_COLOR,
                             ms=8.0, ls='none', zorder=6)
                self.ax.text(px_sel, y_sel + 0.035, f'{self.sel_ref:.2f}',
                             ha='center', va='bottom', fontsize=8,
                             color=SELECT_MARK_COLOR, zorder=7,
                             bbox=dict(fc='white', ec='none', alpha=0.8, pad=0.3))

    def _draw_reference_independent(self):
        """Reference on its own wavelength axis (before a dispersion fit)."""
        pos = self.ax.get_position()
        self.ref_axis = self.fig.add_axes(
            [pos.x0, pos.y0, pos.width, pos.height * 0.47],
            zorder=self.ax.get_zorder() + 1)
        # Separate figure axes: an inset axes is not part of ``fig.axes`` and
        # therefore never receives mouse events.
        self.ref_axis.set_in_layout(False)
        self.ref_axis.patch.set_visible(False)
        self.ref_axis.set_yticks([])
        self.ref_axis.tick_params(axis='x', colors=REF_COLOR, labelsize=7)
        self.ref_axis.xaxis.set_ticks_position('bottom')
        self.ref_axis.xaxis.set_label_position('bottom')
        self.ref_axis.set_xlabel('Wavelength (Å)', color=REF_COLOR, fontsize=8)
        self.ref_axis.spines['bottom'].set_color(REF_COLOR)
        self.ref_axis.spines['top'].set_visible(False)
        self.ref_axis.spines['left'].set_visible(False)
        self.ref_axis.spines['right'].set_visible(False)

        lo, hi = sorted((self.lo_spin.value(), self.hi_spin.value()))
        if not (hi > lo > 0):
            if self.ref is None or not len(self.ref):
                return
            lo, hi = float(np.nanmin(self.ref.wave)), float(np.nanmax(self.ref.wave))
        self.ref_axis.set_xlim(lo, hi)
        self.ref_axis.set_ylim(-0.05, 1.05)

        if self.ref_is_spectrum:
            wave = np.asarray(self.ref.wave, dtype=float)
            flux = np.clip(np.nan_to_num(np.asarray(self.ref.flux, dtype=float)), 0, None)
            sel = (wave >= lo) & (wave <= hi) & np.isfinite(flux)
            if sel.sum() < 2:
                return
            top = max(float(np.nanmax(flux[sel])), 1e-12)
            self._ref_top = top
            self.ref_axis.plot(wave[sel], flux[sel] / top + REFERENCE_OFFSET,
                               color=REF_COLOR, lw=0.65, alpha=0.9, zorder=1,
                               label='reference atlas')
        else:
            lines_wave, lines_strength = self.overlay_lines()
            if not len(lines_wave):
                return
            strength = np.asarray(lines_strength, dtype=float)
            top = max(float(np.nanmax(strength)) if len(strength) else 1.0, 1e-12)
            heights = 0.05 + 0.22 * np.clip(strength / top, 0, 1)
            self.ref_axis.vlines(lines_wave, REFERENCE_OFFSET,
                                 REFERENCE_OFFSET + heights,
                                 color=REF_COLOR, lw=0.65, alpha=0.8, zorder=1,
                                 label='reference lines')
        if self.ref_is_spectrum:
            waves, ys = self._reference_marker_data()
            if len(waves):
                self.ref_axis.plot(waves, ys + 0.015, 'v', color=REF_MARK_COLOR,
                                   ms=3.8, ls='none', alpha=0.9, zorder=4,
                                   label='reference peaks')
            if self.sel_ref is not None:
                ysel = self._reference_y_at(self.sel_ref) + 0.015
                self.ref_axis.plot([self.sel_ref], [ysel], 'v',
                                   color=SELECT_MARK_COLOR, ms=8.0, ls='none',
                                   zorder=6)
                self.ref_axis.text(self.sel_ref, ysel + 0.035,
                                   f'{self.sel_ref:.2f}', ha='center', va='bottom',
                                   fontsize=8, color=SELECT_MARK_COLOR, zorder=7,
                                   bbox=dict(fc='white', ec='none', alpha=0.8, pad=0.3))
        elif self.sel_ref is not None:
            ysel = self._reference_y_at(self.sel_ref)
            self.ref_axis.vlines([self.sel_ref], REFERENCE_OFFSET, ysel + 0.02,
                                 color=SELECT_MARK_COLOR, lw=2.2, zorder=6)
            self.ref_axis.text(self.sel_ref, ysel + 0.045,
                               f'{self.sel_ref:.2f}', ha='center', va='bottom',
                               fontsize=8, color=SELECT_MARK_COLOR, zorder=7,
                               bbox=dict(fc='white', ec='none', alpha=0.8, pad=0.3))
        self.ref_axis.callbacks.connect('xlim_changed',
                                        lambda _ax: self.label_reference_lines())

    def _draw_resolution(self):
        ax = self.res_ax
        ax.set_ylabel('R = λ/Δλ')
        ax.set_xlabel('Wavelength (Å)')
        ax.grid(alpha=0.2)
        info = self.resolution_info
        if not info or not len(info.get('waves', [])):
            ax.text(0.5, 0.55, 'Identify clean lines to estimate R(λ)',
                    transform=ax.transAxes, ha='center', va='center',
                    fontsize=8, color='0.45')
            return
        waves = np.asarray(info['waves'], dtype=float)
        values = np.asarray(info['resolutions'], dtype=float)
        good = np.isfinite(waves) & np.isfinite(values) & (values > 0)
        if not good.any():
            return
        rejected_w = np.asarray(info.get('rejected_waves', []), dtype=float)
        rejected_r = np.asarray(info.get('rejected_resolutions', []), dtype=float)
        if len(rejected_w):
            ax.plot(rejected_w, rejected_r, 'x', ms=4.0, color='0.55',
                    mew=0.8, ls='none', zorder=2)
        ax.plot(waves[good], values[good], 'o', ms=4.0, color='#1565c0',
                mec='white', mew=0.4, ls='none', zorder=3)
        median = float(np.median(values[good]))
        # Scale the y axis to the bulk of the points so far-away outliers
        # (kept in the data, still reachable by zooming) do not flatten the plot.
        v = values[good]
        q1, q3 = np.percentile(v, [25, 75])
        span = max(q3 - q1, 0.1 * median)
        ylo, yhi = max(q1 - 2.5 * span, 0.0), q3 + 2.5 * span
        ax.set_ylim(ylo, yhi)
        n_out = int(np.sum((v < ylo) | (v > yhi)))
        if n_out:
            ax.text(0.005, 0.95, f'{n_out} outlier(s) outside the view',
                    transform=ax.transAxes, ha='left', va='top', fontsize=7, color='0.4')
        tw = np.asarray(info.get('trend_waves', []), dtype=float)
        tr = np.asarray(info.get('trend_resolutions', []), dtype=float)
        if len(tw) and len(tr):
            ax.plot(tw, tr, color='#1565c0', lw=1.0, ls='--',
                    label='trend')
        else:
            ax.axhline(median, color='#1565c0', lw=0.8, ls='--',
                       alpha=0.8)
        ax.text(0.995, 0.95,
                f'median R ≈ {median:,.0f}   (FWHM {info["fwhm_ang"]:.2f} Å, '
                f'n={good.sum()}, '
                f'rejected {info.get("n_rejected", 0)})',
                transform=ax.transAxes, ha='right', va='top', fontsize=7,
                color='#1565c0',
                bbox=dict(fc='white', ec='none', alpha=0.8, pad=1.2))

    def _draw_residuals(self):
        ax = self.resid_ax
        ax.axhline(0, color='0.35', lw=0.8)
        ax.set_ylabel('Residual (Å)')
        ax.set_xlabel('Wavelength (Å)' if self.disp_fit is not None else '')
        ax.xaxis.set_ticks_position('bottom')
        ax.tick_params(axis='x', top=False, labeltop=False,
                       bottom=True, labelbottom=True)
        ax.xaxis.set_label_position('bottom')
        ax.grid(alpha=0.2)
        if self.disp_fit is None:
            ax.text(0.5, 0.55,
                    'Identify at least two lines to show the fit residuals',
                    transform=ax.transAxes, ha='center', va='center',
                    fontsize=8, color='0.45')
            return
        fit = self.disp_fit
        mask = self.disp_mask
        rms = (float(np.sqrt(np.mean(self.disp_residuals[mask] ** 2)))
               if mask.any() else float(np.sqrt(np.mean(self.disp_residuals ** 2))))
        all_x = np.array([p[0] for p in self.all_pairs()], dtype=float)
        all_res = np.array([p[1] for p in self.all_pairs()], dtype=float) - fit(all_x)
        clip_keep = mask if len(mask) == len(all_x) and mask.any() else np.ones(len(all_x), bool)
        threshold = self._clip_threshold(all_res, clip_keep)
        ax.axhspan(-rms, rms, color='0.75', alpha=0.18, zorder=0)
        ax.axhline(rms, color='0.45', lw=0.7, ls=':', zorder=1)
        ax.axhline(-rms, color='0.45', lw=0.7, ls=':', zorder=1)
        for pairs, color, label in ((self.manual, MANUAL_COLOR, 'manual'),
                                    (self.auto, AUTO_COLOR, 'auto')):
            if not pairs:
                continue
            xs = np.array([p[0] for p in pairs], dtype=float)
            plot_xs = fit(xs)
            ys = np.array([p[1] for p in pairs], dtype=float)
            res = ys - fit(xs)
            clipped = (np.abs(res) > threshold if threshold is not None
                       else np.zeros(len(res), bool))
            ax.plot(plot_xs[~clipped], res[~clipped], 'o', ms=4.5, color=color,
                    mec='white', mew=0.4, ls='none', label=label, zorder=3)
            if clipped.any():
                ax.plot(plot_xs[clipped], res[clipped], 'o', ms=5.2, mfc='none',
                        mec=color, ls='none', label='_nolegend_', zorder=3)
        ax.legend(loc='best', fontsize=6.5, ncol=2, framealpha=0.85)
        value = float(self.fit_clip_spin.value())
        mode = self.clip_mode()
        if threshold is None:
            clip_text = 'clipping off'
        elif mode == 'absolute':
            clip_text = f'clip ±{value:g} Å'
        else:
            clip_text = f'clip {value:g}σ = ±{threshold:.4f} Å'
        ax.text(0.995, 0.95,
                f'{dispersion_type(fit)} degree {fit.degree()}   rms = {rms:.3f} Å'
                f'   (dotted ±1σ RMS; {clip_text})',
                transform=ax.transAxes, ha='right', va='top', fontsize=7,
                bbox=dict(fc='white', ec='none', alpha=0.75, pad=1.5))
        wlo, whi = sorted((float(fit(0.0)), float(fit(float(len(self.flux) - 1)))))
        if self.flip_x_box.isChecked():
            ax.set_xlim(whi, wlo)
        else:
            ax.set_xlim(wlo, whi)

    def label_reference_lines(self):
        """Label the reference lines currently visible in the top panel."""
        for t in self._ref_texts:
            t.remove()
        self._ref_texts = []
        if (self.ref is None or self.ref_is_spectrum or len(self.ref) == 0
                or not self.ov_box.isChecked()):
            return
        independent = self.ref_axis is not None and self.disp_fit is None
        if independent:
            x0, x1 = sorted(self.ref_axis.get_xlim())
            target_ax = self.ref_axis
        else:
            model = self.overlay_model()
            if model is None or len(model) < 2:
                return
            if not (np.all(np.diff(model) > 0) or np.all(np.diff(model) < 0)):
                return
            x0, x1 = sorted(self.ax.get_xlim())
            target_ax = self.ax
        m = self.ref_mask()
        wave = np.asarray(self.ref.wave[m], dtype=float)
        if not len(wave):
            return
        strength = (self.ref_strengths[m] / max(np.nanmax(self.ref_strengths), 1e-12)
                    if len(self.ref_strengths) == len(self.ref.wave) else np.ones(len(wave)))
        if independent or self._plot_in_wave:
            px = wave
            visible = (px >= x0) & (px <= x1)
        else:
            px = self.wavelength_to_pixel(wave)
            visible = np.isfinite(px) & (px >= x0) & (px <= x1)
        if visible.sum() > 50:
            return
        for i in np.where(visible)[0]:
            if self.sel_ref is not None and abs(float(wave[i]) - self.sel_ref) < 1e-6:
                continue
            text = '%.1f' % wave[i]
            label = self.ref.labels[i] if i < len(self.ref.labels) else ''
            if label:
                text += ' ' + str(label)
            ypos = REFERENCE_OFFSET + 0.30 + 0.10 * float(np.clip(strength[i], 0, 1))
            self._ref_texts.append(
                target_ax.text(px[i], ypos, text, rotation=90, ha='center', va='bottom',
                               fontsize=5.5, color=REF_COLOR, zorder=2,
                               bbox=dict(fc='white', ec='none', alpha=0.65, pad=0.2)))

    # ------------------------------------------------------------ events
    @staticmethod
    def _tol(ax):
        return 0.015 * abs(np.diff(ax.get_xlim())[0])

    def _visible_peaks(self):
        threshold = self.snr_spin.value()
        return [p for p in self.peaks if p.snr >= threshold]

    def _pick_peak(self, xdata):
        peaks = self._visible_peaks()
        if not peaks:
            return False
        px = np.array([p.pixel for p in peaks])
        i = int(np.argmin(np.abs(px - xdata)))
        tol = max(1.5, min(self._tol(self.ax), 12.0))
        if abs(px[i] - xdata) > tol:
            return False
        self.sel_peak = float(px[i])
        return True

    def _pick_reference(self, xdata, axis=None):
        if self.ref is None or not len(self.ref) or not self.ov_box.isChecked():
            return False
        wave, _strength = self.overlay_lines()
        if axis is self.ref_axis and self.ref_axis is not None:
            if len(wave):
                d = np.abs(np.asarray(wave, dtype=float) - xdata)
                i = int(np.argmin(d))
                tol = max(2.0, min(self._tol(self.ref_axis), 15.0))
                if d[i] <= tol:
                    self.sel_ref = float(wave[i])
                    return True
            return False
        if self.disp_fit is None:
            return False
        if len(wave):
            px = self.wavelength_to_pixel(wave)
            valid = np.isfinite(px)
            if valid.any():
                d = np.abs(np.where(valid, px, np.inf) - xdata)
                i = int(np.argmin(d))
                tol = max(2.0, min(self._tol(self.ax), 12.0))
                if d[i] <= tol:
                    self.sel_ref = float(wave[i])
                    return True
        return False

    def on_click(self, event):
        if event.xdata is not None and event.inaxes in (self.ax, self.ref_axis):
            # Remember which panel was last touched, including while panning.
            self._last_zoom_axis = event.inaxes
        if event.xdata is None or not self.select_mode:
            return
        if self.sel_peak is not None and self.sel_ref is not None:
            # The previous pair was already added; start a fresh selection.
            self.sel_peak = self.sel_ref = None
            self.update_pick_label()
        if self.sel_peak is None and self.sel_ref is None:
            # Reveal the just-completed pair while the next one is selected.
            self._hidden_pair = None
        if self.toolbar.mode:
            self.statusBar().showMessage(
                'Zoom/pan tool is active; click it again to return to line selection.')
            return
        xdata = float(event.xdata)
        if event.inaxes is self.ref_axis and self.ref_axis is not None:
            self._last_zoom_axis = self.ref_axis
            if event.button == 1 and self._pick_reference(xdata, self.ref_axis):
                self.update_pick_label()
                if self.sel_peak is not None:
                    self.add_pair()
                else:
                    self.redraw(keep_view=True)
            return
        if event.inaxes is not self.ax:
            return
        self._last_zoom_axis = self.ax
        if event.button == 3:
            self.remove_nearest(xdata)
            return
        if event.button != 1:
            return
        # Two-step selection on the same axis: after a lamp peak, the next
        # click is meant for the reference line (and vice versa).
        want_reference = self.sel_peak is not None and self.sel_ref is None
        if want_reference:
            picked = self._pick_reference(xdata, self.ax)
            if not picked:
                picked = self._pick_peak(xdata)
        else:
            picked = self._pick_peak(xdata)
            if not picked:
                picked = self._pick_reference(xdata, self.ax)
        if not picked:
            self.statusBar().showMessage(
                'No detected lamp peak or mapped reference line near the click; '
                'zoom in or align the reference first.')
            return
        self.update_pick_label()
        if self.sel_peak is not None and self.sel_ref is not None:
            self.add_pair()
        else:
            self.redraw(keep_view=True)

    def remove_nearest(self, xdata):
        tol = self._tol(self.ax)
        best = None
        for lst in (self.manual, self.auto):
            for p in lst:
                px = (float(self.disp_fit(p[0])) if self._plot_in_wave
                      else p[0])
                d = abs(px - xdata)
                if d <= tol and (best is None or d < best[0]):
                    best = (d, lst, p)
        if best:
            best[1].remove(best[2])
            self.update_model()

    def on_move(self, event):
        if event.xdata is None:
            return
        if event.inaxes is self.ref_axis and self.ref_axis is not None:
            self.statusBar().showMessage('reference λ = %.2f Å' % event.xdata)
            return
        if event.inaxes is self.ax:
            if self._plot_in_wave:
                msg = 'λ = %.2f Å' % event.xdata
                px = self.wavelength_to_pixel(event.xdata)
                if np.isfinite(px):
                    msg += '   pixel≈%.1f' % float(px)
            else:
                msg = 'pixel %.1f' % event.xdata
                wave = self.base_wavelength(event.xdata)
                if wave is not None and np.isfinite(wave):
                    msg += '   λ≈%.2f Å' % float(wave)
                    if self.flip_x_box.isChecked():
                        msg += '   (x flipped)'
            self.statusBar().showMessage(msg)
        elif event.inaxes is self.resid_ax and self.disp_fit is not None:
            pairs = np.asarray(self.all_pairs(), dtype=float).reshape(-1, 2)
            if len(pairs):
                i = int(np.argmin(np.abs(pairs[:, 0] - event.xdata)))
                res = float(pairs[i, 1] - self.disp_fit(pairs[i, 0]))
                self.statusBar().showMessage(
                    f'nearest line: pixel {pairs[i, 0]:.2f}, '
                    f'λ = {pairs[i, 1]:.3f} Å, residual = {res:+.3f} Å')

    def on_scroll(self, event):
        if event.xdata is None:
            return
        step = max(1, int(self.zoom_slider.value())) / 100.0
        scale = 1.0 / (1.0 + step) if event.button == 'up' else 1.0 + step
        if event.inaxes is self.ref_axis and self.ref_axis is not None:
            self._last_zoom_axis = self.ref_axis
            x0, x1 = self.ref_axis.get_xlim()
            a = event.xdata - (event.xdata - x0) * scale
            b = event.xdata + (x1 - event.xdata) * scale
            lo, hi = sorted((a, b))
            if x1 >= x0:
                self.ref_axis.set_xlim(lo, hi)
            else:
                self.ref_axis.set_xlim(hi, lo)
            self.label_reference_lines()
            self.canvas.draw_idle()
            return
        if event.inaxes not in (self.ax, self.resid_ax):
            return
        self._last_zoom_axis = self.ax
        x0, x1 = self.ax.get_xlim()
        a = event.xdata - (event.xdata - x0) * scale
        b = event.xdata + (x1 - event.xdata) * scale
        lo, hi = sorted((a, b))
        if len(self.flux):
            lo = max(lo, -0.5)
            hi = min(hi, len(self.flux) - 0.5)
        if x1 >= x0:
            self.ax.set_xlim(lo, hi)
        else:
            self.ax.set_xlim(hi, lo)
        self.canvas.draw_idle()


def main(argv=None, prog: str = 'specproc.identify'):
    args = build_parser(prog=prog).parse_args(argv)
    if args.auto:
        return run_auto(args)
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv[:1])
    win = LineIdentifier(args.lamp, args.linelist, args.anchors, args.dispersion_axis,
                         reference_path=args.reference, wave_range=args.wave_range,
                         config_path=args.config, output_path=args.output,
                         plot_path=args.plot, flip_x=args.flip_x,
                         poly_type=args.poly_type, arm=args.arm)
    win.show()
    return app.exec_()


if __name__ == '__main__':
    sys.exit(main())
