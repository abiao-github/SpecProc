"""
Interactive lamp-line identification for long-slit spectra.

The window shows the extracted 1-D lamp spectrum (top, versus pixel) above the
reference line list drawn as delta functions (bottom, versus wavelength, zero
baseline). Compare the two patterns by eye, pair a few obvious lines (click a
lamp peak, then the matching reference line), and press "Auto identify": the
manual pairs seed a dispersion solution and all lamp peaks above the SNR
threshold are matched to the catalogue. Save the result and use it as
``[slit] identification_file``.

Run with ``specproc-identify`` or ``python -m specproc.gui.line_identifier``.
"""

import argparse
import sys
from pathlib import Path
from typing import Optional

import numpy as np
from PyQt5 import QtWidgets
from PyQt5.QtCore import Qt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure

from specproc.core.line_identification import (
    find_arc_peaks, fit_polynomial, load_anchors, load_lamp_spectrum, load_linelist,
    refine_with_catalog, save_anchors)

MIN_DETECT_SNR = 3.0


class LineIdentifier(QtWidgets.QMainWindow):
    """Lamp spectrum above, reference delta-function spectrum below."""

    def __init__(self, lamp_path: Optional[str] = None, linelist_path: Optional[str] = None,
                 anchor_path: Optional[str] = None, dispersion_axis: str = 'y'):
        super().__init__()
        self.setWindowTitle('SpecProc - Lamp Line Identification')
        self.resize(1450, 850)
        self.dispersion_axis = dispersion_axis
        self.flux = np.array([])
        self.peaks = np.empty((0, 3))          # pixel, amplitude, snr
        self.ref_wave = np.array([])
        self.ref_label = []
        self.ref_strength = np.array([])
        self.manual = []                       # (pixel, wavelength) picked by hand
        self.auto = []                         # (pixel, wavelength) from auto identification
        self.coef = None
        self.sel_peak: Optional[float] = None
        self.sel_ref: Optional[float] = None
        self.lamp_path = ''
        self.anchor_path = ''
        self._syncing = False
        self._ref_texts = []
        self._build_ui()
        if linelist_path:
            self.load_linelist(linelist_path)
        if lamp_path:
            self.load_lamp(lamp_path)
        if anchor_path:
            self.load_anchor_file(anchor_path)

    # ------------------------------------------------------------------ UI
    def _build_ui(self):
        QW = QtWidgets
        tb = QW.QToolBar('Files')
        tb.setMovable(False)
        self.addToolBar(tb)
        for text, slot in (('Load reference lines…', self.choose_linelist),
                           ('Load lamp spectrum…', self.choose_lamp),
                           ('Load identified lines…', self.choose_anchors),
                           ('Save identified lines…', self.save_dialog)):
            tb.addAction(text).triggered.connect(slot)

        central = QW.QWidget()
        self.setCentralWidget(central)
        layout = QW.QHBoxLayout(central)

        left = QW.QVBoxLayout()
        self.fig = Figure(figsize=(10, 7))
        gs = self.fig.add_gridspec(2, 1, height_ratios=[3, 2], hspace=0.28)
        self.ax = self.fig.add_subplot(gs[0])
        self.ref_ax = self.fig.add_subplot(gs[1])
        self.canvas = FigureCanvasQTAgg(self.fig)
        self.toolbar = NavigationToolbar2QT(self.canvas, self)
        left.addWidget(self.toolbar)
        left.addWidget(self.canvas, 1)
        self.canvas.mpl_connect('button_press_event', self.on_click)
        self.canvas.mpl_connect('motion_notify_event', self.on_move)
        self.canvas.mpl_connect('scroll_event', self.on_scroll)
        self.canvas.mpl_connect('key_press_event', self.on_key)
        self.canvas.setFocusPolicy(Qt.StrongFocus)
        self.ax.callbacks.connect('xlim_changed', lambda ax: self.sync_reference_view())

        opts = QW.QHBoxLayout()
        self.log_box = QW.QCheckBox('Log y')
        self.log_box.setChecked(True)
        self.log_box.stateChanged.connect(lambda *_: self.redraw())
        self.norm_box = QW.QCheckBox('Normalize')
        self.norm_box.setToolTip('Show the lamp spectrum scaled to 0-1 (baseline removed, '
                                 'brightest line = 1)')
        self.norm_box.stateChanged.connect(lambda *_: self.redraw())
        self.sync_box = QW.QCheckBox('Sync reference view')
        self.sync_box.setChecked(True)
        self.sync_box.setToolTip('Once a dispersion model exists (≥2 pairs or a λ hint), '
                                 'zooming the lamp spectrum moves the reference panel too.')
        self.sync_box.stateChanged.connect(lambda *_: self.sync_reference_view())
        self.lo_spin = QW.QDoubleSpinBox()
        self.hi_spin = QW.QDoubleSpinBox()
        for sp in (self.lo_spin, self.hi_spin):
            sp.setRange(0.0, 100000.0)
            sp.setDecimals(1)
            sp.setSingleStep(50.0)
            sp.valueChanged.connect(lambda *_: self.redraw())
        self.w0_spin = QW.QDoubleSpinBox()
        self.w1_spin = QW.QDoubleSpinBox()
        for sp in (self.w0_spin, self.w1_spin):
            sp.setRange(0.0, 30000.0)
            sp.setDecimals(1)
            sp.setSpecialValueText('off')
            sp.valueChanged.connect(lambda *_: self.update_model())
        for w in (self.log_box, self.norm_box, self.sync_box,
                  QW.QLabel('λ hint at first/last pixel'), self.w0_spin, self.w1_spin):
            opts.addWidget(w)
        opts.addStretch(1)
        left.addLayout(opts)
        ov = QW.QHBoxLayout()
        self.ov_box = QW.QCheckBox('Overlay reference lines on lamp spectrum')
        self.ov_box.setToolTip('Draws the catalogue as vertical lines on the lamp spectrum using '
                               'the current model. Align them with shift/stretch (or the arrow '
                               'keys: ←/→ shift, ↑/↓ stretch), then press Fit overlay.')
        self.ov_box.stateChanged.connect(lambda *_: self.redraw(keep_view=True))
        self.flip_box = QW.QCheckBox('Flip')
        self.flip_box.setToolTip('Wavelength decreases with pixel (used before any model exists)')
        self.flip_box.stateChanged.connect(lambda *_: self.redraw(keep_view=True))
        self.shift_spin = QW.QDoubleSpinBox()
        self.shift_spin.setRange(-5000.0, 5000.0)
        self.shift_spin.setDecimals(2)
        self.shift_spin.setSingleStep(1.0)
        self.shift_spin.setSuffix(' Å')
        self.stretch_spin = QW.QDoubleSpinBox()
        self.stretch_spin.setRange(-50.0, 50.0)
        self.stretch_spin.setDecimals(2)
        self.stretch_spin.setSingleStep(0.1)
        self.stretch_spin.setSuffix(' %')
        for sp in (self.shift_spin, self.stretch_spin):
            sp.valueChanged.connect(lambda *_: self.redraw(keep_view=True))
        fit_btn = QW.QPushButton('Fit overlay')
        fit_btn.setToolTip('Match lamp peaks to the overlaid lines and fit the dispersion')
        fit_btn.clicked.connect(self.fit_overlay)
        reset_btn = QW.QPushButton('Reset shift/stretch')
        reset_btn.clicked.connect(self.reset_overlay)
        for w in (self.ov_box, self.flip_box, QW.QLabel('Shift'), self.shift_spin,
                  QW.QLabel('Stretch'), self.stretch_spin, reset_btn, fit_btn):
            ov.addWidget(w)
        ov.addStretch(1)
        left.addLayout(ov)
        rng = QW.QHBoxLayout()
        rng.addWidget(QW.QLabel('Reference wavelength range: low'))
        rng.addWidget(self.lo_spin)
        rng.addWidget(QW.QLabel('upper'))
        rng.addWidget(self.hi_spin)
        rng.addWidget(QW.QLabel('Å   R=λ/Δλ'))
        self.res_spin = QW.QDoubleSpinBox()
        self.res_spin.setRange(0.0, 1000000.0)
        self.res_spin.setDecimals(0)
        self.res_spin.setSingleStep(500.0)
        self.res_spin.setSpecialValueText('off')
        self.res_spin.setToolTip('Spectral resolution. Draws the reference lines broadened to '
                                 'this resolution and sets the final match tolerance to one FWHM.')
        self.res_spin.valueChanged.connect(lambda *_: self.redraw(keep_view=True))
        rng.addWidget(self.res_spin)
        rng.addStretch(1)
        left.addLayout(rng)
        layout.addLayout(left, 1)

        right = QW.QVBoxLayout()
        self.pick_label = QW.QLabel()
        self.pick_label.setWordWrap(True)
        right.addWidget(self.pick_label)
        pick = QW.QHBoxLayout()
        add_btn = QW.QPushButton('Add pair')
        add_btn.setToolTip('Pair the selected lamp peak with the selected reference line')
        add_btn.clicked.connect(self.add_pair)
        cancel = QW.QPushButton('Cancel selection')
        cancel.clicked.connect(self.cancel_selection)
        pick.addWidget(add_btn)
        pick.addWidget(cancel)
        right.addLayout(pick)

        right.addWidget(QW.QLabel('Identified lines (pixel is 0-based)'))
        self.table = QW.QTableWidget(0, 4)
        self.table.setHorizontalHeaderLabels(['Type', 'Pixel', 'λ (Å)', 'Resid. (Å)'])
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.setSelectionBehavior(QW.QAbstractItemView.SelectRows)
        self.table.setEditTriggers(QW.QAbstractItemView.NoEditTriggers)
        right.addWidget(self.table, 1)
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
        self.snr_spin.setValue(8.0)
        self.snr_spin.setToolTip('Only lamp peaks with peak/noise above this value are matched '
                                 '(also set [slit] peak_snr in the config).')
        self.snr_spin.valueChanged.connect(lambda *_: self.redraw(keep_view=True))
        self.order_spin = QW.QSpinBox()
        self.order_spin.setRange(1, 7)
        self.order_spin.setValue(3)
        self.order_spin.valueChanged.connect(lambda *_: self.update_model())
        self.tol_spin = QW.QDoubleSpinBox()
        self.tol_spin.setRange(0.1, 100.0)
        self.tol_spin.setValue(3.0)
        self.tol_spin.setSuffix(' Å')
        form.addRow('Peak SNR threshold', self.snr_spin)
        form.addRow('Polynomial order', self.order_spin)
        form.addRow('Initial match tolerance', self.tol_spin)
        run = QW.QPushButton('Auto identify')
        run.clicked.connect(self.auto_identify)
        form.addRow(run)
        right.addWidget(box)

        self.fit_label = QW.QLabel('')
        self.fit_label.setWordWrap(True)
        right.addWidget(self.fit_label)
        help_txt = QW.QLabel(
            'Click a lamp peak (top), then the matching reference line (bottom): the pair is '
            'added. Right-click near an identified line to remove it. Zoom/pan with the toolbar '
            '(clicks are ignored while a tool is active) or the mouse wheel.')
        help_txt.setWordWrap(True)
        right.addWidget(help_txt)
        wrap = QW.QWidget()
        wrap.setLayout(right)
        wrap.setFixedWidth(400)
        layout.addWidget(wrap)
        self.update_pick_label()
        self.statusBar().showMessage('Load a reference line list and a lamp spectrum to begin.')

    # ------------------------------------------------------------ file I/O
    def choose_linelist(self):
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, 'Reference line list', '', 'Line lists (*.dat *.txt *.csv);;All files (*)')
        if path:
            self.load_linelist(path)

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

    def load_linelist(self, path: str):
        try:
            self.ref_wave, self.ref_label, self.ref_strength = load_linelist(path)
        except Exception as exc:
            QtWidgets.QMessageBox.warning(self, 'Reference lines', str(exc))
            return
        for sp, v in ((self.lo_spin, self.ref_wave[0] - 1), (self.hi_spin, self.ref_wave[-1] + 1)):
            sp.blockSignals(True)
            sp.setValue(max(v, 0.0))
            sp.blockSignals(False)
        self.redraw()
        self.statusBar().showMessage(f'{len(self.ref_wave)} reference lines loaded from {path}')

    def load_lamp(self, path: str):
        try:
            self.flux = load_lamp_spectrum(path, self.dispersion_axis)
        except Exception as exc:
            QtWidgets.QMessageBox.warning(self, 'Lamp spectrum', str(exc))
            return
        self.lamp_path = path
        self.peaks = np.array(find_arc_peaks(self.flux, MIN_DETECT_SNR, with_snr=True),
                              dtype=float).reshape(-1, 3)
        self.sel_peak = self.sel_ref = None
        self.update_pick_label()
        self.coef = self.current_coef()
        self.redraw()
        self.statusBar().showMessage(f'{len(self.flux)} pixels, {len(self.peaks)} peaks '
                                     f'(SNR ≥ {MIN_DETECT_SNR:g}): {path}')

    def load_anchor_file(self, path: str):
        try:
            pairs = load_anchors(path)
        except Exception as exc:
            QtWidgets.QMessageBox.warning(self, 'Identified lines', str(exc))
            return
        self.manual = [(float(x), float(w)) for x, w in pairs]
        self.auto = []
        self.anchor_path = path
        self.update_model()
        self.statusBar().showMessage(f'{len(self.manual)} lines loaded as manual pairs from {path}')

    def save_dialog(self):
        pairs = self.all_pairs()
        if not pairs:
            QtWidgets.QMessageBox.information(self, 'Save', 'No identified lines to save.')
            return
        if self.anchor_path:
            default = self.anchor_path
        elif self.lamp_path:
            default = str(Path(self.lamp_path).with_suffix('')) + '_ids.csv'
        else:
            default = 'ids.csv'
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, 'Save identified lines', default, 'CSV (*.csv);;All files (*)')
        if not path:
            return
        note = (f'lamp: {Path(self.lamp_path).name}; {len(self.manual)} manual + '
                f'{len(self.auto)} auto (SNR>={self.snr_spin.value():g})')
        save_anchors(path, np.array(pairs), note=note)
        self.anchor_path = path
        self.statusBar().showMessage(f'{len(pairs)} lines saved to {path}')

    # ------------------------------------------------------------- state
    def min_tol(self):
        """Final match tolerance: one resolution FWHM at mid-range, else the default 0.4 Å."""
        R = self.res_spin.value()
        if R <= 0:
            return 0.4
        return max((self.lo_spin.value() + self.hi_spin.value()) / 2 / R, 0.05)

    def display_flux(self):
        """Lamp flux as plotted: raw counts, or 0-1 with the baseline removed."""
        if not self.norm_box.isChecked() or len(self.flux) == 0:
            return self.flux
        base = np.percentile(self.flux, 1)
        top = max(np.nanmax(self.flux) - base, 1e-12)
        out = (self.flux - base) / top
        return np.clip(out, 1e-3, None) if self.log_box.isChecked() else out

    def ref_mask(self):
        return (self.ref_wave >= self.lo_spin.value()) & (self.ref_wave <= self.hi_spin.value())

    def all_pairs(self):
        return sorted(self.manual + self.auto)

    def current_coef(self):
        """Dispersion model from identified pairs or the wavelength hints."""
        if len(self.manual) >= 2:
            return fit_polynomial(np.array(self.manual), self.order_spin.value())
        pairs = self.all_pairs()
        if len(pairs) >= 2:
            return fit_polynomial(np.array(pairs), self.order_spin.value())
        w0, w1 = self.w0_spin.value(), self.w1_spin.value()
        if w0 > 0 and w1 > 0 and len(self.flux) > 1:
            return np.polyfit([0, len(self.flux) - 1], [w0, w1], 1)
        return None

    def overlay_model(self):
        """Pixel -> wavelength model for the overlay: base model plus shift and stretch."""
        n = len(self.flux)
        if n < 2:
            return None
        pix = np.arange(n, dtype=float)
        if self.coef is not None:
            base = np.polyval(self.coef, pix)
        else:
            lo, hi = self.lo_spin.value(), self.hi_spin.value()
            if hi <= lo:
                return None
            base = np.linspace(hi, lo, n) if self.flip_box.isChecked() else np.linspace(lo, hi, n)
        mid = base[n // 2]
        return mid + (base - mid) * (1 + self.stretch_spin.value() / 100) + self.shift_spin.value()

    def reset_overlay(self):
        for sp in (self.shift_spin, self.stretch_spin):
            sp.blockSignals(True)
            sp.setValue(0)
            sp.blockSignals(False)
        self.redraw(keep_view=True)

    def fit_overlay(self):
        model = self.overlay_model()
        if model is None or len(self.peaks) == 0 or len(self.ref_wave) == 0:
            self.statusBar().showMessage('Load lamp and line list and set the wavelength range.')
            return
        pix = np.arange(len(model), dtype=float)
        start = np.polyfit(pix, model, min(self.order_spin.value(), 3))
        use = self.peaks[:, 2] >= self.snr_spin.value()
        pairs, _ = refine_with_catalog(self.peaks[use, 0], np.array(self.manual).reshape(-1, 2),
                                       self.ref_wave[self.ref_mask()],
                                       self.order_spin.value(), self.tol_spin.value(),
                                       min_tol=self.min_tol(), initial_coef=start)
        if len(pairs) < 2:
            self.statusBar().showMessage('Too few lines matched; align the overlay better or '
                                         'increase the match tolerance.')
            return
        manual_px = np.array([m[0] for m in self.manual])
        if len(manual_px):
            pairs = pairs[np.min(np.abs(pairs[:, :1] - manual_px[None, :]), axis=1) > 1.0]
        self.auto = [(float(x), float(w)) for x, w in pairs]
        for sp in (self.shift_spin, self.stretch_spin):
            sp.blockSignals(True)
            sp.setValue(0)
            sp.blockSignals(False)
        self.update_model()
        self.statusBar().showMessage(f'Fit overlay: {len(self.auto)} auto + {len(self.manual)} '
                                     f'manual lines; {self.fit_label.text()}')

    def on_key(self, event):
        key = event.key
        if key in ('left', 'right') and self.ov_box.isChecked():
            model = self.overlay_model()
            step = abs(np.median(np.diff(model))) if model is not None else 1.0
            self.shift_spin.setValue(self.shift_spin.value() + (step if key == 'right' else -step))
        elif key in ('up', 'down') and self.ov_box.isChecked():
            self.stretch_spin.setValue(self.stretch_spin.value() + (0.05 if key == 'up' else -0.05))

    def update_model(self):
        self.manual.sort()
        self.coef = self.current_coef()
        self.refresh_table()
        self.redraw(keep_view=True)
        self.sync_reference_view()

    def refresh_table(self):
        rows = [('manual', x, w) for x, w in self.manual] + [('auto', x, w) for x, w in self.auto]
        rows.sort(key=lambda r: r[1])
        pairs = np.array([(r[1], r[2]) for r in rows]).reshape(-1, 2)
        coef = fit_polynomial(pairs, self.order_spin.value()) if len(pairs) > 2 else None
        self.table.setRowCount(len(rows))
        for i, (kind, x, w) in enumerate(rows):
            res = '%.3f' % (w - np.polyval(coef, x)) if coef is not None else ''
            for j, text in enumerate((kind, '%.2f' % x, '%.4f' % w, res)):
                self.table.setItem(i, j, QtWidgets.QTableWidgetItem(text))
        if coef is not None:
            rms = np.std(pairs[:, 1] - np.polyval(coef, pairs[:, 0]))
            self.fit_label.setText(f'{len(self.manual)} manual + {len(self.auto)} auto lines, '
                                   f'order {len(coef) - 1}, rms = {rms:.3f} Å')
        else:
            if self.auto:
                self.fit_label.setText(f'{len(self.auto)} auto-identified lines')
            else:
                self.fit_label.setText(
                    'Manual pairs are optional. Auto identify can start from the low/upper '
                    'wavelength range and tries both dispersion directions.')

    # -------------------------------------------------- picking and editing
    def update_pick_label(self):
        peak = 'none' if self.sel_peak is None else 'pixel %.2f' % self.sel_peak
        ref = 'none' if self.sel_ref is None else '%.4f Å' % self.sel_ref
        self.pick_label.setText(f'Selected lamp peak: {peak}\nSelected reference line: {ref}')

    def cancel_selection(self):
        self.sel_peak = self.sel_ref = None
        self.update_pick_label()
        self.redraw(keep_view=True)

    def add_pair(self):
        if self.sel_peak is None or self.sel_ref is None:
            self.statusBar().showMessage('Select both a lamp peak and a reference line first.')
            return
        x, w = self.sel_peak, self.sel_ref
        self.manual = [m for m in self.manual if abs(m[0] - x) > 1.0 and abs(m[1] - w) > 1e-3]
        self.auto = [a for a in self.auto if abs(a[0] - x) > 1.0 and abs(a[1] - w) > 1e-3]
        self.manual.append((x, w))
        self.sel_peak = self.sel_ref = None
        self.update_pick_label()
        self.update_model()

    def remove_selected(self):
        for r in sorted({i.row() for i in self.table.selectedIndexes()}):
            kind = self.table.item(r, 0).text()
            key = (float(self.table.item(r, 1).text()), float(self.table.item(r, 2).text()))
            lst = self.manual if kind == 'manual' else self.auto
            lst[:] = [p for p in lst if abs(p[0] - key[0]) > 0.01 or abs(p[1] - key[1]) > 1e-4]
        self.update_model()

    def clear_manual(self):
        self.manual = []
        self.update_model()

    def clear_auto(self):
        self.auto = []
        self.update_model()

    def auto_identify(self):
        if len(self.ref_wave) == 0 or len(self.peaks) == 0:
            self.statusBar().showMessage('Load a reference line list and a lamp spectrum first.')
            return
        low, high = self.lo_spin.value(), self.hi_spin.value()
        if high <= low:
            self.statusBar().showMessage('Set a valid wavelength range: upper must exceed low.')
            return
        use = self.peaks[:, 2] >= self.snr_spin.value()
        pairs, _ = refine_with_catalog(self.peaks[use, 0], np.array(self.manual).reshape(-1, 2),
                                       self.ref_wave[self.ref_mask()],
                                       self.order_spin.value(), self.tol_spin.value(),
                                       min_tol=self.min_tol(),
                                       initial_range=(low, high) if not self.manual else None,
                                       initial_pixel_range=(0, len(self.flux) - 1),
                                       peak_amp=self.peaks[use, 1])
        if len(pairs) == 0:
            self.statusBar().showMessage(
                'No lines matched. Check the wavelength range, SNR threshold, line list, '
                'and initial match tolerance.')
            return
        manual_px = np.array([m[0] for m in self.manual])
        if len(manual_px):
            pairs = pairs[np.min(np.abs(pairs[:, :1] - manual_px[None, :]), axis=1) > 1.0]
        self.auto = [(float(x), float(w)) for x, w in pairs]
        self.update_model()
        self.statusBar().showMessage(
            f'Auto identify: {len(self.auto)} lines matched with {len(self.manual)} manual '
            f'({int(use.sum())} peaks with SNR ≥ {self.snr_spin.value():g})')

    # ---------------------------------------------------------- plotting
    def redraw(self, keep_view: bool = False):
        view_top = (self.ax.get_xlim(), self.ax.get_ylim()) if keep_view and len(self.flux) else None
        view_ref = self.ref_ax.get_xlim() if keep_view and len(self.ref_wave) else None
        self._syncing = True
        self.ax.clear()
        self.ref_ax.clear()
        self._ref_texts = []
        thr = self.snr_spin.value()
        if len(self.flux):
            x = np.arange(len(self.flux))
            disp = self.display_flux()
            self.ax.plot(x, disp, color='tab:blue', lw=0.7)
            if self.log_box.isChecked():
                self.ax.set_yscale('log')
                pos = disp[disp > 0]
                self.ax.set_ylim(max(np.percentile(pos, 5) if len(pos) else 1, 1e-3),
                                 np.nanmax(disp) * 3)
            elif self.norm_box.isChecked():
                self.ax.set_ylim(-0.05, 1.1)
            if len(self.peaks):
                strong = self.peaks[:, 2] >= thr
                for sel, color, size in ((strong, 'tab:red', 5), (~strong, '0.6', 3)):
                    p = self.peaks[sel]
                    self.ax.plot(p[:, 0], np.interp(p[:, 0], x, disp) * 1.08, 'v',
                                 color=color, ms=size, ls='none')
            self.ax.set_xlim(-0.5, len(self.flux) - 0.5)
            for pairs, color in ((self.manual, 'tab:green'), (self.auto, 'tab:cyan')):
                for px, wv in pairs:
                    self.ax.axvline(px, color=color, lw=1.0, alpha=0.9)
                    self.ax.text(px, 1.0, '%.1f' % wv, rotation=90, ha='right', va='top',
                                 transform=self.ax.get_xaxis_transform(), color=color, fontsize=7)
            if self.sel_peak is not None:
                self.ax.axvline(self.sel_peak, color='orange', lw=1.6, ls='--')
            model = self.overlay_model() if self.ov_box.isChecked() and len(self.ref_wave) else None
            if model is not None and np.all(np.diff(model) != 0) and \
                    (np.all(np.diff(model) > 0) or np.all(np.diff(model) < 0)):
                m = self.ref_mask() & (self.ref_wave >= model.min()) & (self.ref_wave <= model.max())
                order = np.argsort(model)
                px = np.interp(self.ref_wave[m], model[order], x[order])
                st = self.ref_strength[m] / max(np.nanmax(self.ref_strength), 1e-12)
                self.ax.vlines(px, 0, np.maximum(st, 0.3) * 0.4, color='tab:orange', lw=0.8,
                               alpha=0.85, transform=self.ax.get_xaxis_transform())
        self.ax.set_xlabel('Pixel (0-based)')
        self.ax.set_ylabel('Normalized flux' if self.norm_box.isChecked() else 'Counts')
        self.ax.set_title((Path(self.lamp_path).name if self.lamp_path else 'Lamp spectrum') +
                          f'   (red ▼: peaks with SNR ≥ {thr:g}; grey: weaker)')

        if len(self.ref_wave):
            strength = self.ref_strength / max(np.nanmax(self.ref_strength), 1e-12)
            m = self.ref_mask()
            self.ref_ax.vlines(self.ref_wave[m], 0, np.maximum(strength[m], 0.15), color='k', lw=0.8)
            self.ref_ax.axhline(0, color='0.5', lw=0.6)
            R = self.res_spin.value()
            lo, hi = self.lo_spin.value(), self.hi_spin.value()
            if R > 0 and hi > lo and m.any():
                fw = (lo + hi) / 2 / R
                grid = np.linspace(lo, hi, int(min(max((hi - lo) / (fw / 6), 2), 40000)))
                sig = fw / 2.355
                prof = np.zeros_like(grid)
                for wv, st in zip(self.ref_wave[m], np.maximum(strength[m], 0.15)):
                    j = np.searchsorted(grid, [wv - 5 * sig, wv + 5 * sig])
                    prof[j[0]:j[1]] += st * np.exp(-0.5 * ((grid[j[0]:j[1]] - wv) / sig) ** 2)
                self.ref_ax.plot(grid, np.minimum(prof, 1.2), color='tab:purple', lw=0.7, alpha=0.8)
            for pairs, color, lw in ((self.manual, 'tab:green', 1.2), (self.auto, 'tab:cyan', 0.8)):
                for _, w in pairs:
                    self.ref_ax.axvline(w, color=color, lw=lw, alpha=0.8)
            if self.sel_ref is not None:
                self.ref_ax.axvline(self.sel_ref, color='orange', lw=1.6, ls='--')
            self.ref_ax.set_xlim(*(view_ref or (self.lo_spin.value() - 20, self.hi_spin.value() + 20)))
            self.ref_ax.set_ylim(0, 1.25)
        self.ref_ax.set_xlabel('Wavelength (Å)')
        self.ref_ax.set_ylabel('Reference (relative)')
        if view_top:
            self.ax.set_xlim(*view_top[0])
            self.ax.set_ylim(*view_top[1])
        self._syncing = False
        self.label_reference_lines()
        self.canvas.draw_idle()

    def label_reference_lines(self):
        """Write wavelengths above the stems when few enough are in view."""
        for t in self._ref_texts:
            t.remove()
        self._ref_texts = []
        if len(self.ref_wave) == 0:
            return
        lo, hi = sorted(self.ref_ax.get_xlim())
        sel = np.where((self.ref_wave >= lo) & (self.ref_wave <= hi) & self.ref_mask())[0]
        if len(sel) > 60:
            return
        strength = self.ref_strength / max(np.nanmax(self.ref_strength), 1e-12)
        for i in sel:
            text = '%.2f' % self.ref_wave[i] + (' ' + self.ref_label[i] if self.ref_label[i] else '')
            self._ref_texts.append(self.ref_ax.text(
                self.ref_wave[i], max(strength[i], 0.15) + 0.02, text, rotation=90, ha='center',
                va='bottom', fontsize=7))

    def sync_reference_view(self):
        """Move the reference panel to the wavelength range predicted for the lamp view."""
        if self._syncing or self.coef is None or not self.sync_box.isChecked() \
                or len(self.ref_wave) == 0 or len(self.flux) == 0:
            return
        x0, x1 = self.ax.get_xlim()
        lo, hi = np.polyval(self.coef, [max(x0, 0), min(x1, len(self.flux) - 1)])
        if lo != hi:
            self._syncing = True
            self.ref_ax.set_xlim(*sorted((lo, hi)))
            self._syncing = False
        self.label_reference_lines()
        self.canvas.draw_idle()

    # ------------------------------------------------------------ events
    @staticmethod
    def _tol(ax):
        return 0.015 * abs(np.diff(ax.get_xlim())[0])

    def on_click(self, event):
        if event.inaxes not in (self.ax, self.ref_ax) or event.xdata is None or self.toolbar.mode:
            return
        top = event.inaxes is self.ax
        if event.button == 3:
            self.remove_nearest(event.xdata, top)
            return
        if event.button != 1:
            return
        if top:
            if len(self.peaks) == 0:
                return
            i = int(np.argmin(np.abs(self.peaks[:, 0] - event.xdata)))
            if abs(self.peaks[i, 0] - event.xdata) > self._tol(self.ax):
                self.statusBar().showMessage('No detected peak near the click; zoom in and retry.')
                return
            self.sel_peak = float(self.peaks[i, 0])
        else:
            if len(self.ref_wave) == 0:
                return
            cand = np.where(self.ref_mask())[0]
            if len(cand) == 0:
                return
            i = int(cand[np.argmin(np.abs(self.ref_wave[cand] - event.xdata))])
            if abs(self.ref_wave[i] - event.xdata) > self._tol(self.ref_ax):
                self.statusBar().showMessage('No reference line near the click; zoom in and retry.')
                return
            self.sel_ref = float(self.ref_wave[i])
        self.update_pick_label()
        if self.sel_peak is not None and self.sel_ref is not None:
            self.add_pair()
        else:
            self.redraw(keep_view=True)

    def remove_nearest(self, xdata, top: bool):
        col = 0 if top else 1
        tol = self._tol(self.ax if top else self.ref_ax)
        best = None
        for lst in (self.manual, self.auto):
            for p in lst:
                d = abs(p[col] - xdata)
                if d <= tol and (best is None or d < best[0]):
                    best = (d, lst, p)
        if best:
            best[1].remove(best[2])
            self.update_model()

    def on_move(self, event):
        if event.xdata is None:
            return
        if event.inaxes is self.ax:
            msg = 'pixel %.1f' % event.xdata
            if self.coef is not None:
                msg += '   λ≈%.2f Å' % np.polyval(self.coef, event.xdata)
            self.statusBar().showMessage(msg)
        elif event.inaxes is self.ref_ax:
            self.statusBar().showMessage('λ = %.2f Å' % event.xdata)

    def on_scroll(self, event):
        if event.inaxes not in (self.ax, self.ref_ax) or event.xdata is None:
            return
        scale = 0.7 if event.button == 'up' else 1 / 0.7
        ax = event.inaxes
        x0, x1 = ax.get_xlim()
        new = (event.xdata - (event.xdata - x0) * scale, event.xdata + (x1 - event.xdata) * scale)
        if ax is self.ax and len(self.flux):
            new = (max(new[0], -0.5), min(new[1], len(self.flux) - 0.5))
        ax.set_xlim(*new)
        if ax is self.ref_ax:
            self.label_reference_lines()
        self.canvas.draw_idle()


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog='specproc-identify',
        description='Compare an extracted lamp spectrum with a reference line list, pick lines '
                    'by hand and auto-identify the rest.')
    parser.add_argument('--lamp', help='extracted 1-D lamp spectrum (SpecProc *_spec.fits or text)')
    parser.add_argument('--linelist', help='reference line list (e.g. calib_data/linelists/fear.dat)')
    parser.add_argument('--anchors', help='existing identified-lines file to continue editing')
    parser.add_argument('--dispersion-axis', choices=['x', 'y'], default='y',
                        help='dispersion axis when --lamp is a 2-D image (default: y)')
    args = parser.parse_args(argv)
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv[:1])
    win = LineIdentifier(args.lamp, args.linelist, args.anchors, args.dispersion_axis)
    win.show()
    return app.exec_()


if __name__ == '__main__':
    sys.exit(main())
