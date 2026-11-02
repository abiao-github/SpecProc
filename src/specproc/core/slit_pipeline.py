"""
Long-slit spectrograph reduction pipeline (e.g. Muztagh-Ata 1.93m slit spectrograph).

Unlike the echelle pipeline there are no orders, no inter-order background and no
blaze function: one wide slit band is recorded on the CCD. Steps:

  1. Overscan removal, master bias subtraction
  2. Pixel-response flat (flat / smoothed flat inside the slit band)
  3. Cosmic-ray removal, flat division
  4. Slit trace on the master flat (+ constant star offset) and sky subtraction
  5. Sum and optimal (Horne) extraction along the trace
  6. Wavelength calibration of the arc extracted along the same trace
"""
# 实现狭缝模式下的平场、提取与波长定标流程。

import os
import warnings
import logging
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from scipy import ndimage
from scipy.signal import savgol_filter

from specproc.core.data_structures import Spectrum, SpectraSet
from specproc.core.line_identification import (
    dispersion_type,
    find_arc_peaks,
    fit_dispersion,
    normalize_poly_type,
    measure_resolution,
    refine_with_catalog,
)

logger = logging.getLogger(__name__)


# 查找可用的配置文件路径。
def _resolve_path(path_str: str) -> Optional[Path]:
    """Resolve a path (absolute, cwd-relative or package-relative)."""
    if not path_str:
        return None
    p = Path(path_str).expanduser()
    if p.is_absolute():
        return p if p.exists() else None
    if (Path.cwd() / p).exists():
        return Path.cwd() / p
    pkg = Path(__file__).resolve().parent.parent / p
    return pkg if pkg.exists() else None


# 管理单孔径狭缝光谱归约流程。
class SlitPipeline:
    """Reduction pipeline for single-aperture long-slit spectra."""

    # 读取狭缝归约参数并初始化状态。
    def __init__(self, config):
        self.config = config
        c = config
        self.out = Path(c.get_output_path())
        self.save_plots = c.get_bool('reduce', 'save_plots', True)
        # Pipeline diagnostics are vector PDFs only.
        self.fig_format = 'pdf'
        self.gain = c.get_float('telescope.detector', 'gain', 1.0)
        self.readnoise = c.get_float('telescope.detector', 'readnoise', 5.0)

        self.overscan_start = c.get_int('data', 'overscan_start_column', -1)
        self.overscan_method = c.get('data', 'overscan_method', 'mean_only').strip().lower()
        self.overscan_smooth_window = c.get_int('data', 'overscan_smooth_window', -1)
        self.overscan_poly_order = c.get_int('data', 'overscan_poly_order', 3)
        self.dispersion_axis = c.get('slit', 'dispersion_axis', 'y').strip().lower()
        self.band_x = (c.get_int('slit', 'band_start', -1), c.get_int('slit', 'band_end', -1))
        self.band_frac = c.get_float('slit', 'band_threshold', 0.05)
        self.flat_smooth = c.get_int('slit', 'flat_smooth_pixels', 21)
        self.flat_edge_margin = c.get_int('slit', 'flat_edge_margin', 40)
        self.trace_order = c.get_int('slit', 'trace_order', 2)
        self.trace_nbins = c.get_int('slit', 'trace_bins', 40)
        self.trace_edge_tol = c.get_float('slit', 'trace_edge_tolerance', 3.0)
        self.aperture = c.get_int('slit', 'extract_half_width', 12)
        self.sky_gap = c.get_int('slit', 'sky_inner_offset', 20)
        self.aperture_mode = c.get('slit', 'aperture_mode', 'fixed').strip().lower()
        self.aperture_fwhm_factor = c.get_float('slit', 'aperture_fwhm_factor', 1.7)
        self.sky_fwhm_factor = c.get_float('slit', 'sky_fwhm_factor', 3.0)
        self.sky_order = c.get_int('slit', 'sky_order', 2)
        self.sky_max_order = max(self.sky_order, c.get_int('slit', 'sky_max_order', 3))
        self.sky_edge_margin = c.get_int('slit', 'sky_edge_margin', 5)
        self.snr_valid_fraction = c.get_float('slit', 'snr_valid_fraction', 0.05)
        self.snr_smooth_window = c.get_int('slit', 'snr_smooth_window', 151)
        self.snr_noise_mode = c.get('slit', 'snr_noise_mode', 'local').strip().lower()
        self._valid_range = None
        self._shift_model = None
        self.sky_match_tol = c.get_float('slit', 'sky_match_tolerance', 3.0)
        self.sky_lowres_list = c.get('slit', 'sky_lowres_linelist', 'calib_data/linelists/lowskylines.dat')
        self.sky_id_list = c.get('slit', 'sky_id_linelist', 'calib_data/linelists/skylines.dat')
        self.sky_id_tol = c.get_float('slit', 'sky_id_tolerance', 1.5)
        self.sky_resolution_split = c.get_float('slit', 'sky_resolution_split', 4000.0)
        self.wave_order = c.get_int('slit', 'wave_poly_order', 5)
        self.wave_poly_type = normalize_poly_type(
            c.get('slit', 'wave_poly_type', 'legendre'))
        self.wave_hint = (c.get_float('slit', 'wave_start', 0.0), c.get_float('slit', 'wave_end', 0.0))
        self.match_tol = c.get_float('slit', 'match_tolerance', 3.0)
        self.peak_snr = c.get_float('slit', 'peak_snr', 10.0)
        self.strong_peak_count = int(c.get_float('slit', 'strong_peak_count', 30))
        self.snr_half_window = c.get_int('slit', 'snr_half_window', 50)
        self.identification_override = None
        self._written_arc_wave = set()
        self.cosmic = c.get_bool('reduce', 'cosmic_enabled', True)

    # ------------------------------------------------------------------ I/O
    # 读取图像并处理方向与过扫描区域。
    def _read(self, path: str) -> Tuple[np.ndarray, fits.Header]:
        with fits.open(path) as hdul:
            data = np.asarray(hdul[0].data, dtype=np.float64)
            header = hdul[0].header.copy()
        self._last_overscan_profile = None
        self._last_raw = data.copy()
        if self.overscan_start > 0 and data.shape[1] >= self.overscan_start:
            rows = np.arange(data.shape[0], dtype=float)
            data_end = self.overscan_start - 1
            data_region = data[:, :data_end]
            sections = []
            for key in ('BIASSEC1', 'BIASSEC2'):
                match = re.match(r'\[(\d+):(\d+),', str(header.get(key, '')))
                if match:
                    sections.append((int(match.group(1)) - 1,
                                     int(match.group(2))))
            if not sections:
                sections = [(self.overscan_start + 3, data.shape[1])]
            sections = [(start, min(end, data.shape[1]))
                        for start, end in sections if start < data.shape[1]]
            profiles = []
            split = data_region.shape[1] // 2
            for index, (start, end) in enumerate(sections):
                overscan = data[:, start:end]
                profile = np.array([
                    sigma_clipped_stats(row, sigma=3.0, maxiters=5)[0]
                    for row in overscan
                ])
                if self.overscan_method == 'mean_savgol' and len(rows) >= 5:
                    window = self.overscan_smooth_window
                    if window < 1:
                        window = max(5, len(rows) // 5)
                    window = min(window, len(rows) if len(rows) % 2 else len(rows) - 1)
                    if window % 2 == 0:
                        window -= 1
                    window = max(5, window)
                    level = savgol_filter(profile, window, min(self.overscan_poly_order, window - 1))
                elif self.overscan_method in ('mean_polynomial', 'polynomial'):
                    order = min(self.overscan_poly_order, max(0, len(rows) - 1))
                    coord = np.linspace(-1.0, 1.0, len(rows))
                    coeff = np.polynomial.legendre.legfit(coord, profile, order)
                    level = np.polynomial.legendre.legval(coord, coeff)
                else:
                    level = profile
                profiles.append((rows, profile, np.asarray(level)))
                if index == 0:
                    data_region[:, :split] -= level[:, None]
                else:
                    data_region[:, split:] -= level[:, None]
            self._last_overscan_profile = profiles
            self._last_overscan_split = split
            data = data_region
        if self.dispersion_axis == 'x':
            data = data.T
        return data, header

    def _plot_overscan_profile(self, path: Path, profile, title: str):
        if not self.save_plots or not profile:
            return
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        path.parent.mkdir(parents=True, exist_ok=True)
        method = self.overscan_method
        if method == 'mean_savgol':
            method_title = f' - mean_savgol (order={self.overscan_poly_order})'
        elif method in ('mean_polynomial', 'polynomial'):
            method_title = (f' - mean_polynomial (legendre, '
                            f'order={self.overscan_poly_order})')
        else:
            method_title = ' - mean_only'
        # Only smoothing/fitting produces a model curve; mean_only is points only
        has_model = method != 'mean_only' and method in (
            'mean_savgol', 'mean_polynomial', 'polynomial')
        if has_model:
            fig, (ax, ax_res) = plt.subplots(
                2, 1, figsize=(10, 8), sharex=True,
                gridspec_kw={'height_ratios': [3, 1]})
        else:
            fig, ax = plt.subplots(figsize=(10, 6))
            ax_res = None
        # Both gates share the same row axis
        for name, colour, (rows, measured, model) in zip(
                ('Left', 'Right'), ('tab:blue', 'tab:orange'), profile):
            ax.plot(rows + 1, measured, '.', color=colour, alpha=0.5, ms=4,
                    label=f'{name} overscan')
            if has_model:
                ax.plot(rows + 1, model, '-', color='tab:red' if name == 'Left' else 'tab:green',
                        lw=2, label=f'{name} model')
                ax_res.plot(rows + 1, measured - model, '.', color=colour,
                            alpha=0.6, ms=3, label=name)
        split = getattr(self, '_last_overscan_split', None)
        if split and len(profile) > 1:
            ax.plot([], [], ' ', label=f'Split at column {split} (Left | Right gates)')
        ax.set_ylabel('Overscan level (ADU)')
        ax.set_title(f'{title}{method_title}')
        ax.legend(loc='best')
        last = ax_res if ax_res is not None else ax
        if ax_res is not None:
            ax_res.axhline(0, color='gray', lw=0.8)
            ax_res.set_ylabel('Residual')
            ax_res.legend(loc='best', fontsize=7)
        last.set_xlabel('Row')
        fig.savefig(Path(str(path) + '.pdf'), bbox_inches='tight')
        plt.close(fig)

    # 保存浮点图像及其头信息。
    def _write(self, path: Path, data: np.ndarray, header: Optional[fits.Header] = None):
        path.parent.mkdir(parents=True, exist_ok=True)
        fits.PrimaryHDU(np.asarray(data, dtype=np.float32), header).writeto(path, overwrite=True)

    # 用红色细线圆圈标出检测到的宇宙线像素。
    def _plot_cosmic_marked(self, path: Path, image: np.ndarray, mask: np.ndarray, name: str):
        if not self.save_plots:
            return
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        fig, ax = plt.subplots(figsize=(8, 8))
        vmin, vmax = np.nanpercentile(image, [1, 99])
        im = ax.imshow(image, origin='lower', aspect='auto', cmap='viridis', vmin=vmin, vmax=vmax)
        fig.colorbar(im, ax=ax).set_label('Counts')
        ys, xs = np.where(mask)
        if xs.size:
            ax.scatter(xs, ys, s=10, facecolors='none', edgecolors='red',
                       linewidths=0.6, alpha=0.9, label=f'Cosmic pixels: {xs.size}')
            ax.legend(loc='upper right')
        ax.set_title(f'{name} detected cosmic rays (marked in red)')
        ax.set_xlabel('spatial pixel')
        ax.set_ylabel('dispersion pixel')
        fig.tight_layout()
        fig.savefig(Path(str(path) + '.' + self.fig_format), dpi=150)
        plt.close(fig)

    # 绘制图像并叠加轨迹标记。
    def _plot_image(self, path: Path, data: np.ndarray, title: str, overlay=None):
        if not self.save_plots:
            return
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        lo, hi = np.nanpercentile(data, [2, 99.5])
        fig, ax = plt.subplots(figsize=(6, 7))
        im = ax.imshow(data, origin='lower', aspect='auto', cmap='viridis', vmin=lo, vmax=hi)
        fig.colorbar(im, ax=ax).set_label('Counts')
        if overlay:
            for xs, ys, style in overlay:
                ax.plot(xs, ys, style, lw=0.8, ms=2)
        ax.set_xlabel('spatial pixel')
        ax.set_ylabel('dispersion pixel')
        ax.set_title(title)
        fig.savefig(Path(str(path) + '.' + self.fig_format), dpi=110, bbox_inches='tight')
        plt.close(fig)

    # --------------------------------------------------------------- stage 1
    # 合并偏置帧并保存主偏置图。
    def _master_bias(self, bias_files: List[str]) -> np.ndarray:
        stack = np.array([self._read(f)[0] for f in bias_files])
        master = np.median(stack, axis=0)
        bias_path = self.out / 'step1_basic' / 'bias_subtracted' / 'master_bias'
        self._write(Path(str(bias_path) + '.fits'), master)
        self._plot_image(bias_path, master, 'Master bias')
        logger.info(f"Master bias from {len(bias_files)} frames: median={np.median(master):.2f}")
        return master

    # 扣除主偏置并保存校正图像。
    def _debias(self, path: str, bias: np.ndarray, write: bool = True):
        """Bias-subtract while keeping the science data on its native pixel grid."""
        data, hdr = self._read(path)
        data = data - bias
        stem = Path(path).stem
        if write:
            out = self.out / 'step1_basic' / 'bias_subtracted' / f'{stem}_bias_corrected'
            self._write(Path(str(out) + '.fits'), data, hdr)
            self._plot_image(out, data, f'{stem} bias-corrected')
        return data, hdr

    # --------------------------------------------------------------- stage 2
    # 确定狭缝光带的空间范围。
    def _detect_band(self, flat: np.ndarray) -> Tuple[int, int]:
        if self.band_x[0] > 0 and self.band_x[1] > 0:
            return self.band_x[0] - 1, self.band_x[1]
        prof = ndimage.uniform_filter1d(np.median(flat, axis=0), 15)
        peak = prof.max()
        base = np.median(prof[:50]) if prof.size > 100 else 0.0
        thr = base + self.band_frac * (peak - base)
        idx = np.where(prof > thr)[0]
        return int(idx.min()), int(idx.max()) + 1

    # 构建像素平场并追踪狭缝区域。
    def _make_flat(self, flat_files: List[str], bias: np.ndarray):
        stack = np.array([self._debias(f, bias)[0] for f in flat_files])
        master = np.median(stack, axis=0)
        x0, x1 = self._detect_band(master)
        logger.info(f"Slit band detected at spatial columns {x0 + 1}-{x1} (1-based)")

        smooth = ndimage.median_filter(master, size=(self.flat_smooth, self.flat_smooth))
        smooth = ndimage.gaussian_filter(smooth, self.flat_smooth / 3.0)
        pix = np.ones_like(master)
        valid = smooth > 0.05 * np.nanmax(smooth[:, x0:x1])
        region = np.zeros_like(valid)
        m = min(self.flat_edge_margin, (x1 - x0) // 4)
        region[:, x0 + m:x1 - m] = True
        use = valid & region
        pix[use] = master[use] / smooth[use]
        pix[(pix < 0.5) | (pix > 1.5) | ~np.isfinite(pix)] = 1.0

        d = self.out / 'step4_flat_corrected'
        self._write(d / 'master_flat.fits', master)
        self._write(d / 'pixel_flat.fits', pix)
        self._plot_image(d / 'master_flat', master, 'Master flat',
                         [([x0, x0], [0, master.shape[0]], 'r-'), ([x1, x1], [0, master.shape[0]], 'r-')])
        self._plot_image(d / 'pixel_flat', pix, 'Pixel flat (normalized)',
                         [([x0, x0], [0, pix.shape[0]], 'r-'), ([x1, x1], [0, pix.shape[0]], 'r-')])
        d2 = self.out / 'step2_trace'
        centre, left, right, (ty, tx) = self._trace_flat(master, (x0, x1))
        yy = np.arange(master.shape[0])
        self._flat_trace = (centre, left, right)
        self._plot_trace(d2 / 'master_flat_trace', master, centre, 'Master flat trace',
                         points=(tx, ty))
        np.savetxt(d2 / 'master_flat_trace.txt', np.column_stack([yy + 1, centre + 1, left + 1, right + 1]),
                   header='dispersion_pixel  centre  left_edge  right_edge  (spatial columns, 1-based)',
                   fmt='%.3f')
        self._flat_trace = (centre, left, right)
        self._master_flat = master
        return pix, (x0, x1)

    # --------------------------------------------------------------- stage 3
    # 狭缝内沿空间方向平滑的天光线模型，供宇宙线探测先行扣除。
    def _line_model(self, data: np.ndarray, size: int = 41) -> np.ndarray:
        """Running median along the slit in the curvature-corrected wavelength frame.

        Only this background model is interpolated; the science pixels remain on their
        native grid. Sky lines fill the slit, so L.A.Cosmic would flag them as cosmic rays.
        """
        model = np.zeros(data.shape, dtype=float)
        if self._band is None:
            return model
        c0, c1 = self._band
        c0, c1 = max(c0 - 10, 0), min(c1 + 10, data.shape[1])
        if c1 - c0 < size:
            return model
        source = np.asarray(data[:, c0:c1], dtype=float)
        if self._shift_model is None:
            model[:, c0:c1] = ndimage.median_filter(
                source, size=(1, size), mode='nearest')
            return model

        aligned, native_rows = self._align_columns(data, c0, c1)
        aligned_model = ndimage.median_filter(aligned, size=(1, size), mode='nearest')
        rows = np.arange(data.shape[0], dtype=float)
        for j, x in enumerate(range(c0, c1)):
            aligned_row = np.interp(rows, native_rows[:, j], rows)
            model[:, x] = ndimage.map_coordinates(
                aligned_model[:, j], [aligned_row], order=3, mode='nearest')
        return model

    # 以模型加邻域残差中值替换被标记像素（L.A.Cosmic 带背景时的清理结果不可用）。
    @staticmethod
    def _fill_masked(data: np.ndarray, mask: np.ndarray, model: np.ndarray, half: int = 2):
        out = np.array(data, dtype=float)
        resid = np.where(mask, np.nan, out - model)
        ny, nx = out.shape
        for y, x in np.argwhere(mask):
            win = resid[max(y - half, 0):y + half + 1, max(x - half, 0):x + half + 1]
            val = np.nanmedian(win) if np.isfinite(win).any() else 0.0
            out[y, x] = model[y, x] + val
        return out

    # 探测并校正图像中的宇宙射线。
    def _cosmic(self, data: np.ndarray) -> np.ndarray:
        self._last_cosmic_mask = None
        self._last_cosmic_input = data
        if self.cosmic:
            try:
                import astroscrappy
                model = self._line_model(data)
                mask, _ = astroscrappy.detect_cosmics(
                    data + 0.0, inbkg=model,
                    gain=self.gain, readnoise=self.readnoise,
                    sigclip=self.config.get_float('reduce', 'cosmic_sigclip', 5.0),
                    objlim=self.config.get_float('reduce', 'cosmic_objlim', 5.0),
                    satlevel=1e9, cleantype='meanmask')
                logger.info(f"Cosmic rays flagged: {int(mask.sum())} px")
                self._last_cosmic_mask = mask
                data = self._fill_masked(data, mask, model)
            except Exception as exc:
                logger.warning(f"Cosmic ray removal failed: {exc}")
        return data

    # --------------------------------------------------------------- stage 4
    # 追踪平场狭缝中心线及边界。
    def _trace_flat(self, master: np.ndarray, band: Tuple[int, int]):
        """Trace the slit band on the master flat (as the echelle orders are).

        In each row bin the left/right band edges are located at 50% of the
        local plateau level; the centre line is their midpoint. Centre and
        edges are fitted with polynomials of order `trace_order`.
        """
        x0, x1 = band
        ny = master.shape[0]
        pad = 40
        a0, a1 = max(x0 - pad, 0), min(x1 + pad, master.shape[1])
        edges = np.linspace(0, ny, self.trace_nbins + 1).astype(int)
        ys, cl, cr, lev = [], [], [], []
        for a, b in zip(edges[:-1], edges[1:]):
            prof = ndimage.gaussian_filter1d(np.median(master[a:b, a0:a1], axis=0), 3)
            plateau = np.percentile(prof, 90)
            base = np.median(np.concatenate([prof[:20], prof[-20:]]))
            if plateau - base <= 0:
                continue
            half = base + 0.5 * (plateau - base)
            idx = np.where(prof > half)[0]
            if idx.size < 50:
                continue
            i0, i1 = idx.min(), idx.max()
            # linear interpolation at the half-level crossings
            xl = i0 - 1 + (half - prof[i0 - 1]) / (prof[i0] - prof[i0 - 1]) if i0 > 0 else i0
            xr = i1 + (prof[i1] - half) / (prof[i1] - prof[i1 + 1]) if i1 < len(prof) - 1 else i1
            ys.append(0.5 * (a + b))
            cl.append(a0 + xl)
            cr.append(a0 + xr)
            lev.append(plateau - base)
        ys, cl, cr, lev = map(np.asarray, (ys, cl, cr, lev))
        good = lev > 0.1 * np.max(lev)
        if good.sum() < self.trace_order + 2:
            raise RuntimeError("Flat trace failed: too few valid bins")
        mid = 0.5 * (cl + cr)
        for _ in range(4):
            coef = np.polyfit(ys[good], mid[good], self.trace_order)
            res = mid - np.polyval(coef, ys)
            sg = 1.4826 * np.median(np.abs(res[good])) + 0.3
            good = good & (np.abs(res) < 3 * sg)
        yy = np.arange(ny)
        centre = np.polyval(np.polyfit(ys[good], mid[good], self.trace_order), yy)
        left = np.polyval(np.polyfit(ys[good], cl[good], self.trace_order), yy)
        right = np.polyval(np.polyfit(ys[good], cr[good], self.trace_order), yy)
        logger.info(f"Flat trace: centre {centre[0] + 1:.1f}->{centre[-1] + 1:.1f}, "
                    f"slit width {np.median(right - left):.1f} px")
        return centre, left, right, (ys[good], mid[good])

    # 独立追迹目标源/灯谱的狭缝中心（逐行）。
    def _trace_source(self, img: np.ndarray, band: Tuple[int, int],
                      initial: np.ndarray, label: str = '') -> np.ndarray:
        """Trace a source independently, using the flat trace only as a window guess."""
        ny = img.shape[0]
        x0, x1 = band
        yy = np.arange(ny, dtype=float)
        edges = np.linspace(0, ny, self.trace_nbins + 1).astype(int)
        ys, centres, flux = [], [], []
        half_window = max(20, int(0.20 * (x1 - x0)))
        for a, b in zip(edges[:-1], edges[1:]):
            if b - a < 3:
                continue
            c = float(initial[(a + b) // 2])
            lo = int(max(c - half_window, x0))
            hi = int(min(c + half_window, x1))
            if hi - lo < 5:
                continue
            prof = np.median(img[a:b, lo:hi], axis=0)
            prof = ndimage.gaussian_filter1d(prof - np.median(prof), 1.5)
            pk = int(np.argmax(prof))
            if not np.isfinite(prof[pk]) or prof[pk] <= 0:
                continue
            l2, h2 = max(pk - 7, 0), min(pk + 8, prof.size)
            w = np.clip(prof[l2:h2], 0, None)
            if w.sum() <= 0:
                continue
            centre = lo + (np.arange(l2, h2) * w).sum() / w.sum()
            ys.append(0.5 * (a + b))
            centres.append(centre)
            flux.append(float(prof[pk]))
        if len(ys) < self.trace_order + 2:
            logger.info("%s trace: too few valid bins; using flat trace", label or 'source')
            return np.asarray(initial, dtype=float)
        ys = np.asarray(ys, dtype=float)
        centres = np.asarray(centres, dtype=float)
        flux = np.asarray(flux, dtype=float)
        # Faint bins (the continuum fades along dispersion) give noise centroids
        # that bend a polynomial; keep only bins well above the brightest ones.
        good = flux > 0.15 * np.percentile(flux, 90)
        order = min(self.trace_order, max(1, len(ys) - 2))
        # Extrapolating a curved fit far beyond the reliable bins is unsafe.
        if good.sum() >= 2 and np.ptp(ys[good]) < 0.7 * np.ptp(ys):
            order = 1
        coef = None
        for _ in range(4):
            if good.sum() < order + 2:
                break
            coef = np.polyfit(ys[good], centres[good], order)
            res = centres - np.polyval(coef, ys)
            sg = 1.4826 * np.median(np.abs(res[good])) + 0.2
            new = good & (np.abs(res) < 3.0 * sg)
            if np.array_equal(new, good):
                break
            good = new
        if coef is None:
            logger.info("%s trace fit failed; using flat trace", label or 'source')
            return np.asarray(initial, dtype=float)
        trace = np.polyval(coef, yy)
        # The star trace must follow the slit-edge midpoint (flat trace) apart
        # from a constant plus a slow linear drift; a bent fit means noisy bins
        # drove it, so rebuild it as the edge trend + constant + linear term.
        init = np.asarray(initial, dtype=float)
        dev = trace - init
        resid = dev - np.polyval(np.polyfit(yy, dev, 1), yy)
        if np.ptp(resid) > self.trace_edge_tol:
            dy = centres - np.interp(ys, yy, init)
            use = good.copy()
            deg = 1 if use.sum() >= 5 and np.ptp(ys[use]) >= 0.5 * np.ptp(ys) else 0
            for _ in range(4):
                if use.sum() < deg + 2:
                    deg = 0
                cf = np.polyfit(ys[use], dy[use], deg)
                r = dy - np.polyval(cf, ys)
                sg = 1.4826 * np.median(np.abs(r[use])) + 0.2
                new_use = use & (np.abs(r) < 3.0 * sg)
                if new_use.sum() < 2 or np.array_equal(new_use, use):
                    break
                use = new_use
            logger.warning("%s trace deviates from the slit-edge trend by %.1f px; "
                           "using edge trend + %s fitted on %d bins",
                           label or 'source', np.ptp(resid),
                           'linear drift' if deg else 'constant offset', use.sum())
            trace = init + np.polyval(cf, yy)
        # Guard against a fit running away from the slit.
        if (np.any(trace < x0) or np.any(trace > x1) or
                np.nanmax(np.abs(trace - initial)) > half_window):
            logger.warning("%s trace fit unstable; using flat trace", label or 'source')
            return np.asarray(initial, dtype=float)
        logger.info("%s trace: %.1f -> %.1f px (%d bins)",
                    label or 'source', trace[0], trace[-1], good.sum())
        return trace

    def _save_trace(self, d: Path, name: str, trace: np.ndarray):
        """Write the per-frame trace centre as text."""
        d.mkdir(parents=True, exist_ok=True)
        np.savetxt(d / f'{name}_trace.txt',
                   np.column_stack([np.arange(len(trace)) + 1, trace + 1]),
                   header='dispersion_pixel  trace_centre (1-based)',
                   fmt='%.3f')

    # 测量星像相对平场中心的偏移。
    def _star_offset(self, img: np.ndarray, centre: np.ndarray, band: Tuple[int, int]) -> float:
        """Constant offset (px) of the star from the flat-derived slit centre line."""
        x0, x1 = band
        ny = img.shape[0]
        edges = np.linspace(0, ny, self.trace_nbins + 1).astype(int)
        offs, fl = [], []
        for a, b in zip(edges[:-1], edges[1:]):
            c = centre[(a + b) // 2]
            lo = int(max(c - 0.4 * (x1 - x0), x0))
            hi = int(min(c + 0.4 * (x1 - x0), x1))
            prof = np.median(img[a:b, lo:hi], axis=0)
            prof = ndimage.gaussian_filter1d(prof - np.median(prof), 1.5)
            pk = int(np.argmax(prof))
            l2, h2 = max(pk - 6, 0), min(pk + 7, prof.size)
            w = np.clip(prof[l2:h2], 0, None)
            if w.sum() <= 0:
                continue
            offs.append(lo + (np.arange(l2, h2) * w).sum() / w.sum() - c)
            fl.append(prof[pk])
        offs, fl = np.asarray(offs), np.asarray(fl)
        if offs.size < 3:
            raise RuntimeError("Star offset measurement failed: too few valid bins")
        sel = fl > 0.3 * np.median(fl)
        med = np.median(offs[sel])
        sel &= np.abs(offs - med) < 3 * (1.4826 * np.median(np.abs(offs[sel] - med)) + 0.3)
        return float(np.median(offs[sel]))

    # 估算星像空间剖面的半高全宽。
    def _measure_fwhm(self, img: np.ndarray, band: Tuple[int, int]) -> Optional[float]:
        """FWHM (px) of the spatial profile, stacked along the whole dispersion axis.

        The profile is aligned on the master-flat centre line plus the star offset,
        collapsed with a median over row bins (robust to cosmic rays), the
        background is removed with a linear fit to its wings and a Gaussian is fitted.
        Returns None if the profile is too faint or the fit is implausible.
        """
        from scipy.optimize import curve_fit
        centre = self._flat_trace[0]
        trace = centre + self._star_offset(img, centre, band)
        ny, nx = img.shape
        half = 40
        k = np.arange(-half, half + 1)
        cols = np.clip(np.rint(trace)[:, None].astype(int) + k[None, :], 0, nx - 1)
        cut = np.take_along_axis(img, cols, axis=1)
        edges = np.linspace(0, ny, self.trace_nbins + 1).astype(int)
        prof = np.mean([np.median(cut[a:b], axis=0) for a, b in zip(edges[:-1], edges[1:])], axis=0)
        wing = np.abs(k) > 25
        bg = np.polyfit(k[wing], prof[wing], 1)
        prof = prof - np.polyval(bg, k)
        noise = 1.4826 * np.median(np.abs(prof[wing])) + 1e-9
        if prof.max() < 8 * noise:
            return None
        def g(x, a, m, sg):
            return a * np.exp(-0.5 * ((x - m) / sg) ** 2)
        try:
            (a, m, sg), _ = curve_fit(g, k, prof, p0=[prof.max(), k[np.argmax(prof)], 3.0], maxfev=2000)
        except Exception:
            return None
        fwhm = 2.3548 * abs(sg)
        return float(fwhm) if 2.0 <= fwhm <= 30.0 and abs(m) < 15 else None

    # 根据星像宽度设置提取孔径。
    def _set_aperture(self, science_files, bias, band):
        if self.aperture_mode != 'auto':
            logger.info(f"Aperture (fixed): half-width {self.aperture} px, sky starts at {self.sky_gap} px")
            return
        fwhms = []
        for f in science_files:
            data, _ = self._debias(f, bias, write=False)
            f_i = self._measure_fwhm(data, band)
            logger.info(f"Spatial FWHM of {Path(f).name}: " + (f"{f_i:.2f} px" if f_i else "not measurable"))
            if f_i:
                fwhms.append(f_i)
        if not fwhms:
            logger.warning(f"FWHM not measurable for any frame; using fixed aperture "
                           f"(half-width {self.aperture} px, sky from {self.sky_gap} px)")
            return
        fwhm = float(np.median(fwhms))
        self.aperture = int(np.ceil(self.aperture_fwhm_factor * fwhm))
        self.sky_gap = int(np.ceil(self.sky_fwhm_factor * fwhm))
        logger.info(f"Auto aperture: median FWHM={fwhm:.2f} px -> extraction half-width "
                    f"{self.aperture} px, sky region starts at {self.sky_gap} px")

    # 拟合并扣除狭缝区域的天空背景。
    # ------------------------------------------------- slit curvature (2D distortion)
    _DIST_BINS = 12      # column bins across the slit in which arc lines are measured
    _DIST_XORDER = 4     # polynomial order along the slit
    _DIST_YORDER = 2     # polynomial order along the dispersion axis (the tilt varies slowly)

    @staticmethod
    def _line_centroid(spec: np.ndarray, px: float, noise: float, half: int = 6,
                       search: float = 3.0):
        from scipy.optimize import curve_fit
        i = int(round(px))
        a, b = max(i - half, 0), min(i + half + 1, len(spec))
        if b - a < 7:
            return None
        xs, ys = np.arange(a, b, dtype=float), spec[a:b]
        amp = float(ys.max() - ys.min())
        if amp < 10 * noise:
            return None

        def g(x, h, mu, sg, base):
            return base + h * np.exp(-0.5 * ((x - mu) / sg) ** 2)
        try:
            popt, _ = curve_fit(g, xs, ys, p0=[amp, px, 1.5, ys.min()],
                                bounds=([0, px - search, 0.7, -np.inf],
                                        [np.inf, px + search, 4.0, np.inf]))
        except Exception:
            return None
        return float(popt[1])

    # 在各列区间测量孤立谱线相对中心列的行位移。
    def _measure_shifts(self, img: np.ndarray, band: Tuple[int, int]):
        x0, x1 = band
        lo, hi = x0 + 25, x1 - 25
        edges = np.linspace(lo, hi, self._DIST_BINS + 1).astype(int)
        xc = 0.5 * (edges[:-1] + edges[1:] - 1)
        spec = np.array([np.nansum(img[:, a:b], axis=1) for a, b in zip(edges[:-1], edges[1:])])
        noise = []
        for sp in spec:
            d = np.diff(sp)
            noise.append(1.4826 * np.median(np.abs(d - np.median(d))) / np.sqrt(2) + 1e-9)
        iref = int(np.argmin(np.abs(xc - 0.5 * (x0 + x1))))
        peaks = self._find_peaks(spec[iref], self.peak_snr)
        cen = np.array([p[0] for p in peaks], dtype=float)
        amp = np.array([p[1] for p in peaks], dtype=float)
        if len(cen) > 1:
            gap = np.minimum(np.r_[np.inf, np.diff(cen)], np.r_[np.diff(cen), np.inf])
            iso = gap > 10
            cen, amp = cen[iso], amp[iso]
        pts = []
        for c in cen:
            y0 = self._line_centroid(spec[iref], c, noise[iref])
            if y0 is None:
                continue
            for b in range(len(xc)):
                cb = y0 if b == iref else self._line_centroid(spec[b], y0, noise[b])
                if cb is not None:
                    pts.append((xc[b], y0, cb - y0))
        return xc, iref, np.asarray(pts, dtype=float).reshape(-1, 3), float(cen[np.argmax(amp)]) if len(cen) else 0.0

    # 拟合位移场 s(x, y)：参考列上行 y 的谱线在列 x 处位于 y + s(x, y)。
    def _fit_shift(self, xc, iref, pts, ny: int, band):
        from numpy.polynomial import legendre as L
        if len(pts) < 60:
            return None
        xref = float(xc[iref])
        model = {'xref': xref, 'xhalf': max(xc.max() - xref, xref - xc.min()),
                 'xmin': float(xc.min()), 'xmax': float(xc.max()), 'ymid': ny / 2.0,
                 'band': (int(band[0]), int(band[1])), 'ny': ny}
        terms = [(i, j) for j in range(1, self._DIST_XORDER + 1) for i in range(self._DIST_YORDER + 1)]

        def design(x, y):
            t = (np.clip(x, model['xmin'], model['xmax']) - xref) / model['xhalf']
            v = (y - model['ymid']) / model['ymid']
            return np.stack([t ** j * L.legval(v, np.eye(i + 1)[i]) for i, j in terms], axis=-1)
        model['design'] = design
        A = design(pts[:, 0], pts[:, 1])
        use = np.ones(len(pts), bool)
        for _ in range(6):
            coef, *_ = np.linalg.lstsq(A[use], pts[use, 2], rcond=None)
            res = pts[:, 2] - A @ coef
            sg = 1.4826 * np.median(np.abs(res[use])) + 1e-6
            new = np.abs(res) < 3 * max(sg, 0.02)
            if (new == use).all() or new.sum() < len(terms) + 5:
                break
            use = new
        model.update(coef=coef, used=use, pts=pts, n=int(use.sum()),
                     rms=float(np.sqrt(np.mean(res[use] ** 2))))
        return model

    @staticmethod
    def _shift(model, x, y):
        x, y = np.broadcast_arrays(np.asarray(x, dtype=float), np.asarray(y, dtype=float))
        return model['design'](x, y) @ model['coef']

    # 仅为弯曲诊断图生成拉平预览，不用于后续科学处理。
    def _rectify(self, img: np.ndarray) -> np.ndarray:
        mdl = self._shift_model
        if mdl is None:
            return img
        ny, nx = img.shape
        c0, c1 = max(mdl['band'][0] - 60, 0), min(mdl['band'][1] + 60, nx)
        yy = np.arange(ny, dtype=float)[:, None]
        xx = np.arange(c0, c1, dtype=float)[None, :]
        shift = self._shift(mdl, xx, yy)
        cols = np.broadcast_to(np.arange(c1 - c0, dtype=float)[None, :], shift.shape)
        out = np.array(img, dtype=float)
        out[:, c0:c1] = ndimage.map_coordinates(
            np.asarray(img[:, c0:c1], dtype=float), [yy + shift, cols], order=3, mode='nearest')
        return out

    def _align_columns(self, img: np.ndarray, c0: int, c1: int):
        """Interpolate only a slit-band view onto common dispersion rows (a row then maps to one wavelength)."""
        ny = img.shape[0]
        rows = np.arange(ny, dtype=float)
        columns = np.arange(c0, c1, dtype=float)
        shift = self._shift(self._shift_model, columns[None, :], rows[:, None])
        native_rows = rows[:, None] + shift
        if np.any(np.diff(native_rows, axis=0) <= 0):
            raise RuntimeError("Slit-curvature mapping is not monotonic in dispersion pixels")
        aligned = np.empty((ny, c1 - c0), dtype=float)
        for j, x in enumerate(range(c0, c1)):
            aligned[:, j] = ndimage.map_coordinates(
                np.asarray(img[:, x], dtype=float), [native_rows[:, j]],
                order=3, mode='nearest')
        return aligned, native_rows

    # 由灯谱建立畸变模型；后续只在背景模型中使用，不重采样处理图像。
    def _build_distortion(self, flat_files, arc_files, bias):
        self._shift_model = None
        if not arc_files:
            logger.warning("No lamp frames: curvature-aware sky modeling is unavailable")
            return
        flat = np.median([self._debias(f, bias, write=False)[0] for f in flat_files], axis=0)
        band = self._detect_band(flat)
        arc = np.median([self._debias(f, bias, write=False)[0] for f in arc_files], axis=0)
        xc, iref, pts, ycut = self._measure_shifts(arc, band)
        model = self._fit_shift(xc, iref, pts, arc.shape[0], band)
        if model is None:
            logger.warning(f"Slit curvature: only {len(pts)} measurements; "
                           "curvature-aware sky modeling is unavailable")
            return
        self._shift_model = model
        rect = self._rectify(arc)
        resid = model['rms']
        span = float(np.ptp(self._shift(model, np.array([model['xmin'], model['xmax']]),
                                        np.full(2, model['ymid']))))
        logger.info(f"Slit curvature model: {model['n']}/{len(pts)} line measurements from "
                    f"{len(arc_files)} lamp frame(s), fit rms={model['rms']:.3f} px, line shift across "
                    f"the slit {span:.2f} px at mid-spectrum, inlier fit residual "
                    f"{resid:.3f} px (rms)")
        self._plot_distortion(self.out / 'step5_extraction' / 'slit_distortion',
                              model, arc, rect, ycut, resid)

    # 绘制谱线弯曲的测量、拟合位移场和仅供诊断的对齐预览。
    def _plot_distortion(self, path, model, before, after, ycut, resid):
        if not self.save_plots:
            return
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        pts, use, ny = model['pts'], model['used'], model['ny']
        fig, ax = plt.subplots(2, 2, figsize=(12, 9))
        a1, a2, a3, a4 = ax.ravel()
        sc = a1.scatter(pts[use, 0], pts[use, 2], c=pts[use, 1], s=8, cmap='viridis')
        if (~use).any():
            a1.scatter(pts[~use, 0], pts[~use, 2], s=14, facecolors='none', edgecolors='0.4',
                       label='rejected')
            a1.legend(fontsize=7)
        fig.colorbar(sc, ax=a1).set_label('line row (dispersion pixel)')
        xs = np.linspace(model['xmin'], model['xmax'], 60)
        line_rows = np.unique(np.round(pts[use, 1], 1))
        cmap = plt.get_cmap('viridis')
        norm = plt.Normalize(float(line_rows.min()), float(line_rows.max())
                             if len(line_rows) > 1 else float(line_rows.min()) + 1.0)
        for yv in line_rows:
            a1.plot(xs, self._shift(model, xs, np.full_like(xs, yv)),
                    lw=0.65, color=cmap(norm(yv)), alpha=0.75)
        a1.axvline(model['xref'], color='0.4', ls=':', lw=0.8,
                   label='reference column (zero shift)')
        n_lines = len(line_rows)
        a1.set_xlabel('spatial pixel')
        a1.set_ylabel('line shift relative to the reference column (rows)')
        a1.set_title(f'{n_lines} measured arc lines; showing each fitted line trace '
                     f'(rms={model["rms"]:.3f} px)')
        a1.legend(fontsize=6, loc='best')
        gy = np.arange(0, ny, 8, dtype=float)
        gx = np.linspace(model['xmin'], model['xmax'], 80)
        grid = self._shift(model, gx[None, :], gy[:, None])
        im = a2.imshow(grid, origin='lower', aspect='auto', cmap='coolwarm',
                       extent=[gx[0], gx[-1], 0, ny])
        fig.colorbar(im, ax=a2).set_label('shift (rows)')
        a2.set_xlabel('spatial pixel')
        a2.set_ylabel('dispersion pixel')
        a2.set_title('Fitted shift field s(x, y)')
        x0, x1 = model['band']
        c0, c1 = max(x0 - 20, 0), min(x1 + 20, before.shape[1])
        r0, r1 = int(max(ycut - 25, 0)), int(min(ycut + 25, ny))
        for a, img, title in ((a3, before, 'Lamp on native pixel grid'),
                              (a4, after, f'Diagnostic-only aligned lamp preview '
                               f'(inlier fit residual {resid:.3f} px rms)')):
            cut = img[r0:r1, c0:c1]
            lo, hi = np.nanpercentile(cut, [2, 99.5])
            a.imshow(cut, origin='lower', aspect='auto', cmap='viridis', vmin=lo, vmax=hi,
                     extent=[c0, c1, r0, r1])
            a.axhline(ycut, color='r', ls='--', lw=0.6)
            a.set_xlabel('spatial pixel')
            a.set_ylabel('dispersion pixel')
            a.set_title(title)
        fig.suptitle('Slit-curvature model (science images remain unrectified)', fontsize=11)
        fig.tight_layout()
        fig.savefig(Path(str(path) + '.' + self.fig_format), dpi=110, bbox_inches='tight')
        plt.close(fig)

    def _sky_bounds(self, ny: int):
        """Per-row slit edges (from the master-flat trace) used for the sky fit and subtraction."""
        _, left, right = self._flat_trace
        return np.asarray(left, dtype=float)[:ny], np.asarray(right, dtype=float)[:ny]

    def _fit_sky_row(self, xn: np.ndarray, values: np.ndarray, use: np.ndarray):
        """Fit one row's cross-slit sky profile with an order chosen by BIC.

        The order runs from sky_order to sky_max_order, so rows with a strong
        sky line (whose profile across the slit is curved) get a higher order
        while faint rows stay linear. Returns Legendre-free power coefficients
        in the normalised column xn, zero-padded to sky_max_order + 1, or None.
        """
        lo, hi = self.sky_order, self.sky_max_order
        best = None
        for degree in range(lo, hi + 1):
            sel = use.copy()
            coef = None
            for _ in range(4):
                if sel.sum() < max(degree + 3, 10):
                    coef = None
                    break
                design = np.vander(xn[sel], N=degree + 1, increasing=True)
                coef, *_ = np.linalg.lstsq(design, values[sel], rcond=None)
                residual = values[sel] - design @ coef
                centre = np.median(residual)
                sigma = 1.4826 * np.median(np.abs(residual - centre))
                if sigma <= 1e-9:
                    break
                keep = np.abs(residual - centre) < 3.0 * sigma
                if keep.all():
                    break
                idx = np.flatnonzero(sel)
                sel = np.zeros_like(sel)
                sel[idx[keep]] = True
            if coef is None:
                continue
            fit_all = np.vander(xn[use], N=degree + 1, increasing=True) @ coef
            res_all = values[use] - fit_all
            # Score on the inliers so clipped blends/stars do not decide the order.
            scale = 1.4826 * np.median(np.abs(res_all - np.median(res_all))) + 1e-9
            inl = np.abs(res_all - np.median(res_all)) < 3.0 * scale
            n = max(int(inl.sum()), degree + 3)
            rss = float(np.sum(res_all[inl] ** 2)) + 1e-12
            bic = n * np.log(rss / n) + (degree + 1) * np.log(n)
            if best is None or bic < best[0]:
                best = (bic, coef)
        if best is None:
            return None
        out = np.zeros(hi + 1)
        out[:len(best[1])] = best[1]
        return out

    def _subtract_sky(self, img: np.ndarray, trace: np.ndarray, band: Tuple[int, int]):
        """Fit the sky per row inside the slit edges and subtract it up to the edges.

        The edges follow the (tilted) slit edges of the master-flat trace, not a fixed
        column range. Pixels within sky_edge_margin of an edge are left out of the fit
        (illumination fall-off) but the sky model is extended to the edge itself.
        """
        if self._shift_model is not None:
            return self._subtract_curved_sky(img, trace)

        ny, nx = img.shape
        left, right = self._sky_bounds(ny)
        xi = np.arange(nx)
        sky = np.zeros_like(img)
        fm = self.sky_edge_margin
        xmid = 0.5 * (nx - 1)
        xscale = max(0.5 * nx, 1.0)
        xn = (xi - xmid) / xscale
        for y in range(ny):
            fit_m = ((xi >= left[y] + fm) & (xi <= right[y] - fm)
                     & (np.abs(xi - trace[y]) > self.sky_gap))
            if fit_m.sum() < 10:
                continue
            c = self._fit_sky_row(xn, img[y], fit_m)
            if c is None:
                continue
            apply = (xi >= np.ceil(left[y])) & (xi <= np.floor(right[y]))
            sky[y, apply] = np.vander(xn[apply], N=len(c), increasing=True) @ c
        return img - sky, sky

    def _subtract_curved_sky(self, img: np.ndarray, trace: np.ndarray):
        """Fit sky at common wavelengths without resampling the science image."""
        ny, nx = img.shape
        left, right = self._sky_bounds(ny)
        c0 = max(0, int(np.floor(np.nanmin(left))))
        c1 = min(nx, int(np.ceil(np.nanmax(right))) + 1)
        rows = np.arange(ny, dtype=float)
        columns = np.arange(c0, c1, dtype=float)
        aligned, native_rows = self._align_columns(img, c0, c1)

        native_index = np.empty((ny, c1 - c0), dtype=int)
        for j in range(c1 - c0):
            native_y = np.interp(rows, native_rows[:, j], rows)
            native_index[:, j] = np.clip(np.rint(native_y).astype(int), 0, ny - 1)

        fit_mask = (
            (columns[None, :] >= left[native_index] + self.sky_edge_margin)
            & (columns[None, :] <= right[native_index] - self.sky_edge_margin)
            & (np.abs(columns[None, :] - trace[native_index]) > self.sky_gap)
            & np.isfinite(aligned)
        )
        xmid = 0.5 * (c0 + c1 - 1)
        xscale = max(0.5 * (c1 - c0), 1.0)
        xn = (columns - xmid) / xscale
        degree = self.sky_max_order
        coefficients = np.full((ny, degree + 1), np.nan, dtype=float)
        for y in range(ny):
            use = fit_mask[y]
            if use.sum() < 10:
                continue
            c = self._fit_sky_row(xn, aligned[y], use)
            if c is not None:
                coefficients[y] = c

        valid = np.flatnonzero(np.isfinite(coefficients[:, 0]))
        if not len(valid):
            raise RuntimeError("Could not fit any rows of the curvature-aware sky model")
        for k in range(degree + 1):
            coefficients[:, k] = np.interp(rows, valid, coefficients[valid, k])

        sky = np.zeros_like(img, dtype=float)
        for j, x in enumerate(range(c0, c1)):
            aligned_row = np.interp(rows, native_rows[:, j], rows)
            x_value = (x - xmid) / xscale
            values = np.zeros(ny, dtype=float)
            for k in range(degree + 1):
                values += ndimage.map_coordinates(
                    coefficients[:, k], [aligned_row], order=3, mode='nearest') * x_value ** k
            apply = (x >= np.ceil(left)) & (x <= np.floor(right))
            sky[apply, x] = values[apply]
        return img - sky, sky

    def _extract_sky(self, img: np.ndarray, trace: np.ndarray, band: Tuple[int, int]) -> np.ndarray:
        """1D sky spectrum measured directly in the sky region of the unsubtracted image.

        With a curvature model, sky samples are aligned to common dispersion rows before averaging,
        then sampled back onto the target trace's native dispersion rows.
        """
        left, right = self._sky_bounds(img.shape[0])
        fm = self.sky_edge_margin
        if self._shift_model is not None:
            ny, nx = img.shape
            c0 = max(0, int(np.floor(np.nanmin(left))))
            c1 = min(nx, int(np.ceil(np.nanmax(right))) + 1)
            rows = np.arange(ny, dtype=float)
            columns = np.arange(c0, c1, dtype=float)
            aligned, native_rows = self._align_columns(img, c0, c1)
            native_index = np.empty((ny, c1 - c0), dtype=int)
            for j in range(c1 - c0):
                native_y = np.interp(rows, native_rows[:, j], rows)
                native_index[:, j] = np.clip(np.rint(native_y).astype(int), 0, ny - 1)
            region = (
                (columns[None, :] >= left[native_index] + fm)
                & (columns[None, :] <= right[native_index] - fm)
                & (np.abs(columns[None, :] - trace[native_index]) > self.sky_gap)
            )
            vals = np.where(region, aligned, np.nan)
            with np.errstate(all='ignore'), warnings.catch_warnings():
                warnings.simplefilter('ignore')
                med = np.nanmedian(vals, axis=1)
                mad = 1.4826 * np.nanmedian(np.abs(vals - med[:, None]), axis=1)
                clipped = np.where(
                    np.abs(vals - med[:, None]) <= 3 * mad[:, None] + 1e-9,
                    vals, np.nan)
                spectrum = np.nanmean(clipped, axis=1)

            star_rows = rows + self._shift(
                self._shift_model, np.asarray(trace, dtype=float), rows)
            if np.any(np.diff(star_rows) <= 0):
                raise RuntimeError("Target-trace curvature mapping is not monotonic")
            target_rows = np.interp(rows, star_rows, rows)
            return ndimage.map_coordinates(
                np.nan_to_num(spectrum, nan=0.0), [target_rows], order=3,
                mode='nearest')

        xi = np.arange(img.shape[1])
        region = ((xi[None, :] >= left[:, None] + fm) & (xi[None, :] <= right[:, None] - fm)
                  & (np.abs(xi[None, :] - trace[:, None]) > self.sky_gap))
        vals = np.where(region, img, np.nan)
        with np.errstate(all='ignore'), warnings.catch_warnings():
            warnings.simplefilter('ignore')
            med = np.nanmedian(vals, axis=1)
            mad = 1.4826 * np.nanmedian(np.abs(vals - med[:, None]), axis=1)
            clipped = np.where(np.abs(vals - med[:, None]) <= 3 * mad[:, None] + 1e-9, vals, np.nan)
            return np.nanmean(clipped, axis=1)

    def _plot_trace(self, path, img, trace, title, aperture=None, sky_gap=None,
                    points=None):
        """Cut-out around the slit band with trace centre, slit edges and, optionally,
        the extraction aperture and sky region boundaries."""
        if not self.save_plots:
            return
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        ny, nx = img.shape
        _, left, right = self._flat_trace
        pad = 25
        lo = max(0, int(np.floor(left.min())) - pad)
        hi = min(nx, int(np.ceil(right.max())) + pad)
        yy = np.arange(ny, dtype=float)
        panel = img[:, lo:hi]
        fig, ax = plt.subplots(figsize=(8, 6))
        vmin, vmax = np.nanpercentile(panel, [2, 99])
        im = ax.imshow(panel, origin='lower', aspect='auto', cmap='viridis',
                       vmin=vmin, vmax=vmax, extent=(lo, hi - 1, 0, ny - 1))
        fig.colorbar(im, ax=ax).set_label('Counts')
        ax.plot(trace, yy, color='tab:red', lw=0.9, label='trace')
        ax.plot(left, yy, color='yellow', lw=0.9, label='slit edge')
        ax.plot(right, yy, color='yellow', lw=0.9)
        if points is not None:
            ax.plot(points[0], points[1], 'c.', ms=2, label='edge measurements')
        if aperture is not None:
            ax.plot(trace + aperture, yy, 'tab:orange', lw=0.8, ls='--', label='aperture')
            ax.plot(trace - aperture, yy, 'tab:orange', lw=0.8, ls='--')
        if sky_gap is not None:
            ax.plot(trace + sky_gap, yy, 'tab:cyan', lw=0.8, ls='--', label='sky starts')
            ax.plot(trace - sky_gap, yy, 'tab:cyan', lw=0.8, ls='--')
        ax.set_xlim(lo, hi - 1)
        ax.set_ylim(0, ny - 1)
        ax.set_xlabel('spatial pixel')
        ax.set_ylabel('dispersion pixel')
        ax.set_title(title)
        ax.legend(loc='upper right', fontsize=7, framealpha=0.85)
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(Path(str(path) + '.pdf'), dpi=130, bbox_inches='tight')
        plt.close(fig)

    def _extract(self, img: np.ndarray, trace: np.ndarray, sky_sub_noise=None,
                 half_width: Optional[int] = None):
        ny, nx = img.shape
        hw = self.aperture if half_width is None else int(half_width)
        offs = np.arange(-hw, hw + 1)
        cols = np.rint(trace)[:, None].astype(int) + offs[None, :]
        cols = np.clip(cols, 0, nx - 1)
        cut = np.take_along_axis(img, cols, axis=1)
        flux_sum = cut.sum(axis=1)

        # Spatial profile in trace-relative coordinates (sub-pixel shifts accounted for)
        rel = cols - trace[:, None]
        tot = np.where(np.abs(flux_sum) > 0, flux_sum, np.nan)
        bright = flux_sum > np.percentile(flux_sum, 50)
        bins = np.arange(-hw - 0.5, hw + 1.5, 1.0)
        prof = np.zeros(len(offs))
        for i in range(len(offs)):
            sel = (rel >= bins[i]) & (rel < bins[i + 1]) & bright[:, None]
            vals = (cut / tot[:, None])[sel]
            prof[i] = np.nanmedian(vals) if vals.size else 0.0
        prof = np.clip(prof, 0, None)
        prof /= prof.sum() if prof.sum() > 0 else 1.0

        var = self.readnoise ** 2 / (self.gain ** 2) + np.clip(cut, 0, None) / max(self.gain, 1e-9)
        # outlier rejection in profile fit
        P = prof[None, :]
        model = flux_sum[:, None] * P
        ok = (cut - model) ** 2 < 25 * var
        w = ok * P / var
        num = (w * cut).sum(axis=1)
        den = (w * P).sum(axis=1)
        flux_opt = np.where(den > 0, num / np.where(den > 0, den, 1), flux_sum)
        err_opt = np.where(den > 0, 1 / np.sqrt(np.where(den > 0, den, 1)), np.nan)
        return flux_sum, flux_opt, err_opt, prof

    # --------------------------------------------------------------- stage 6
    # 检测并定位标定光谱的发射峰。
    @staticmethod
    def _find_peaks(spec: np.ndarray, snr: float = 10.0, top: Optional[int] = None):
        """Peaks above snr x (global clipped-MAD sigma); optionally only the `top`
        highest-S/N ones. Returned as (centre, amplitude), sorted by position."""
        found = find_arc_peaks(spec, snr, with_snr=True, noise_mode='global')
        if top is not None and len(found) > top:
            found = sorted(found, key=lambda p: -p[2])[:top]
        return [(c, a) for c, a, _ in sorted(found)]

    # 读取并排序参考波长目录。
    def _load_catalog(self) -> np.ndarray:
        path = _resolve_path(self.config.get('telescope.linelist', 'full_linelist', ''))
        if path is None:
            logger.warning("Full line list not found (%s); wavelength solution will use "
                           "the anchor lines only",
                           self.config.get('telescope.linelist', 'full_linelist', ''))
            return np.array([])
        w = []
        for line in path.read_text(errors='ignore').splitlines():
            parts = line.split()
            if not parts or parts[0].startswith('#'):
                continue
            try:
                w.append(float(parts[0]))
            except ValueError:
                pass
        return np.array(sorted(w))

    # 读取像素与波长的人工对应表。
    def _anchor_path_for_header(self, header) -> str:
        """Select the R/B anchor file for MSTG data when available."""
        configured = self.config.get('slit', 'identification_file', '')
        tel = str(header.get('TELESCOP', '')).strip().upper()
        ins = str(header.get('INSTRUME', '')).strip().upper()
        if 'MOST' not in tel and 'MUZTAGH' not in tel and 'MSTG' not in tel:
            return configured
        if 'SPECTRUM-B' not in ins and not ins.endswith('-B'):
            return configured
        p = Path(configured)
        if not p.name.endswith('_R.csv'):
            return configured
        candidate = p.with_name(p.name[:-6] + '_B.csv')
        resolved = _resolve_path(str(candidate))
        if resolved is not None:
            logger.info("Using blue-arm anchor file: %s", resolved)
            return str(candidate)
        logger.warning("Blue anchor file %s not found; using %s", candidate, configured)
        return configured

    def _load_id_file(self) -> Optional[np.ndarray]:
        path = _resolve_path(self.identification_override or
                             self.config.get('slit', 'identification_file', ''))
        if path is None:
            return None
        rows = []
        for line in path.read_text().splitlines():
            line = line.split('#')[0].strip()
            if not line:
                continue
            parts = line.replace(',', ' ').split()
            try:
                rows.append((float(parts[0]), float(parts[1])))
            except (ValueError, IndexError):
                continue
        if len(rows) < 2:
            logger.warning(f"Identification file {path} has too few lines")
            return None
        logger.info("Using identification file: %s", path)
        logger.info(f"Using {len(rows)} identified (anchor) lines from {path}")
        return np.array(rows)

    # 依据波长提示匹配谱线与目录。
    def _auto_identify(self, peaks, npix, catalog):
        """Match peaks to catalogue assuming a linear solution from wave_start/wave_end hints."""
        w0, w1 = self.wave_hint
        if w0 <= 0 or w1 <= 0 or len(catalog) == 0:
            return None
        pix = np.array([p[0] for p in peaks])
        domain = (0.0, float(npix - 1))
        fit = fit_dispersion([0.0, float(npix - 1)], [w0, w1], 1,
                             self.wave_poly_type, domain=domain)
        if fit is None:
            return None
        pairs = None
        for it in range(8):
            wl = fit(pix)
            tol = max(self.match_tol / (it + 1), 0.4)
            idx = np.searchsorted(catalog, wl)
            idx = np.clip(idx, 1, len(catalog) - 1)
            near = np.where(np.abs(catalog[idx - 1] - wl) < np.abs(catalog[idx] - wl), idx - 1, idx)
            d = catalog[near] - wl
            sel = np.abs(d) < tol
            if sel.sum() < self.wave_order + 2:
                return None
            pairs = np.column_stack([pix[sel], catalog[near][sel]])
            order = min(self.wave_order, 1 if it < 2 else self.wave_order)
            fit = fit_dispersion(pairs[:, 0], pairs[:, 1], order,
                                 self.wave_poly_type, domain=domain)
            if fit is None:
                return None
        return pairs

    @staticmethod
    def _filter_strong_anchors(anchors, strong_peaks, tolerance: float = 1.5):
        """Keep anchor pairs whose pixel coincides with one of the strongest peaks."""
        if anchors is None or len(anchors) == 0 or not strong_peaks:
            return anchors
        centres = np.asarray([p[0] for p in strong_peaks], dtype=float)
        if not len(centres):
            return anchors
        keep = []
        for x, w in anchors:
            if float(np.min(np.abs(centres - x))) <= max(tolerance, 1.5):
                keep.append((float(x), float(w)))
        return np.asarray(keep, dtype=float).reshape(-1, 2)

    # 拟合并评估标定帧的波长解。
    def _wavelength_solution(self, arc_flux: np.ndarray, d: Path, tag: str = 'arc'):
        npix = len(arc_flux)
        d.mkdir(parents=True, exist_ok=True)
        peaks = self._find_peaks(arc_flux, self.peak_snr)
        pairs = self._load_id_file()
        # Anchor / seed matching uses the highest-S/N peaks (as many as there are anchors);
        # full-catalogue matching and resolution use every detected peak.
        n_strong = len(pairs) if pairs is not None else self.strong_peak_count
        strong_peaks = self._find_peaks(arc_flux, self.peak_snr, top=n_strong)
        logger.info("Arc: %d peaks above %.1f sigma (global MAD), %d strongest used for anchors",
                    len(peaks), self.peak_snr, len(strong_peaks))

        if pairs is not None:
            strong_anchors = self._filter_strong_anchors(pairs, strong_peaks)
            if len(strong_anchors) >= self.wave_order + 2:
                logger.info("Anchor file: %d -> %d pairs coincide with the %d strongest peaks",
                            len(pairs), len(strong_anchors), len(strong_peaks))
                pairs = strong_anchors
            else:
                logger.warning("Only %d anchor pairs coincide with the strongest peaks; "
                               "keeping all %d anchors", len(strong_anchors), len(pairs))
        source = 'identification_file'
        if pairs is not None:
            catalog = self._load_catalog()
            if len(catalog):
                n_anchor = len(pairs)
                pairs, _ = refine_with_catalog(
                    [p[0] for p in peaks], pairs, catalog,
                    self.wave_order, self.match_tol,
                    poly_type=self.wave_poly_type)
                logger.info(f"Anchors {n_anchor} -> {len(pairs)} lines after full-catalogue matching")
                source = 'anchors+catalog'
        if pairs is None:
            catalog = self._load_catalog()
            seed = self._auto_identify(strong_peaks, npix, catalog)
            if seed is not None and len(catalog):
                pairs, _ = refine_with_catalog(
                    [p[0] for p in peaks], seed, catalog,
                    self.wave_order, self.match_tol,
                    poly_type=self.wave_poly_type)
            else:
                pairs = seed
            source = 'auto_match'
        if pairs is None or len(pairs) < self.wave_order + 2:
            logger.warning("No wavelength solution: provide [slit] identification_file or "
                           "wave_start/wave_end hints. Spectra are saved on a pixel axis.")
            return None, None, source

        x, y = pairs[:, 0], pairs[:, 1]
        domain = (0.0, float(npix - 1))
        keep = np.ones(len(x), bool)
        fit = fit_dispersion(x, y, self.wave_order, self.wave_poly_type, domain=domain)
        for _ in range(5):
            if fit is None:
                break
            res = y - fit(x)
            rms = float(np.sqrt(np.mean(res[keep] ** 2)))
            new = np.abs(res) < max(3 * rms, 0.3)
            if new.sum() < self.wave_order + 2 or (new == keep).all():
                break
            keep = new
            fit = fit_dispersion(x[keep], y[keep], self.wave_order,
                                 self.wave_poly_type, domain=domain)
        fit = fit_dispersion(x[keep], y[keep], self.wave_order,
                             self.wave_poly_type, domain=domain)
        if fit is None:
            logger.warning("Not enough wavelength lines for the requested polynomial order")
            return None, None, source
        rms = float(np.sqrt(np.mean((y - fit(x))[keep] ** 2)))
        # Resolution uses every significant detected arc peak, independent of
        # whether it entered the full-catalogue wavelength solution.
        resolution = measure_resolution(
            np.arange(len(arc_flux), dtype=float), arc_flux, peaks, fit)
        if resolution:
            logger.info("Arc resolution: R ~ %.0f (FWHM %.2f A, %d lines, %d rejected)",
                        resolution['resolution'], resolution['fwhm_ang'],
                        resolution['n'], resolution.get('n_rejected', 0))
        logger.info(f"Wavelength solution ({source}, {dispersion_type(fit)} "
                    f"degree {fit.degree()}): {keep.sum()}/{len(x)} lines, rms={rms:.3f} A")
        return fit, (x, y, keep, rms, resolution), source

    # ------------------------------------------------------------------ main
    # 执行狭缝图像校正、提取与定标。
    def run(self, science_files: List[str], bias_files: List[str],
            flat_files: List[str], arc_files: List[str]) -> Dict[str, SpectraSet]:
        logger.info("=" * 60)
        logger.info("STARTING SLIT SPECTRUM REDUCTION PIPELINE")
        logger.info("=" * 60)
        self.out.mkdir(parents=True, exist_ok=True)

        for f in bias_files + flat_files + arc_files + science_files:
            data, hdr = self._read(f)
            d1 = self.out / 'step1_basic' / 'overscan_corrected'
            out = d1 / f'{Path(f).stem}_overscan_corrected'
            self._write(Path(str(out) + '.fits'), data, hdr)
            self._plot_image(d1 / f'{Path(f).stem}_raw', self._last_raw,
                             f'{Path(f).stem} raw')
            self._plot_image(out, data, f'{Path(f).stem} overscan-corrected')
            self._plot_overscan_profile(
                self.out / 'step1_basic' / 'overscan_corrected' /
                f'{Path(f).stem}_overscan_profile',
                self._last_overscan_profile,
                f'{Path(f).stem} overscan profile')

        bias = self._master_bias(bias_files)
        if flat_files:
            pixflat, band = self._make_flat(flat_files, bias)
            self._band = band
        else:
            raise RuntimeError("Slit mode needs flat frames: the trace is derived from the master flat.")

        self._set_aperture(science_files, bias, band)

        # Lamp frames are never combined: each science frame uses the single
        # lamp frame closest in time (DATE-OBS, then file mtime as tie-break).
        arcs = []
        d2 = self.out / 'step2_trace'
        for f in arc_files:
            img_a, hdr_a = self._debias(f, bias)
            arc_name = Path(f).stem
            self._plot_trace(d2 / f'{arc_name}_trace', img_a, self._flat_trace[0],
                             f'{arc_name} trace')
            # Pixel-sensitivity correction only; no sky subtraction for lamp frames.
            img_a = img_a / pixflat
            d4f = self.out / 'step4_flat_corrected'
            self._write(d4f / f'{arc_name}_flat2d_corrected.fits', img_a, hdr_a)
            self._plot_trace(d4f / f'{arc_name}_flat2d_corrected', img_a,
                             self._flat_trace[0], f'{arc_name} flat2d-corrected')
            # Lamp light fills the whole slit: extract along the flat centre line
            # over the full slit width (no point-source trace exists).
            arc_trace = self._flat_trace[0]
            self._save_trace(d2, arc_name, arc_trace)
            arcs.append((self._frame_time(f, hdr_a), os.path.getmtime(f), f,
                         self._extract_calib(img_a, hdr_a, arc_name,
                                             trace=arc_trace, plot=False), hdr_a))
        master_flat_flux = None
        if flat_files:
            master_flat_flux = self._extract_calib(
                self._master_flat, fits.Header(), 'master_flat', plot=False)
        self._valid_range = self._find_valid_range(master_flat_flux, arcs)
        for _, _, f_arc, flux_arc, hdr_arc in arcs:
            self._write_calib_pixel(Path(f_arc).stem, flux_arc, hdr_arc)
        if master_flat_flux is not None:
            self._write_calib_pixel('master_flat', master_flat_flux, fits.Header())

        # Every lamp gets its own wavelength solution, lamp spectrum and diagnostics.
        d5 = self.out / 'step7_wavelength'
        arc_solutions = {}
        for _, _, f_arc, flux_arc, hdr_arc in arcs:
            stem = Path(f_arc).stem
            self.identification_override = self._anchor_path_for_header(hdr_arc)
            fit, info, _src = self._wavelength_solution(flux_arc, d5, stem)
            if fit is None:
                continue
            arc_solutions[f_arc] = (fit, info)
            self._write_lamp_wave(d5, stem, flux_arc, hdr_arc, info, fit)
            self._plot_arc(d5 / f'{stem}_arc', flux_arc, info, fit)
            if master_flat_flux is not None:
                self._write_lamp_wave(d5, 'master_flat', master_flat_flux, fits.Header(),
                                      info, fit, title='master flat')

        # Curvature is needed for sky modeling in extraction, not for step 1,
        # flat construction, or wavelength calibration.
        if arc_files:
            self._build_distortion(flat_files, arc_files, bias)

        results: Dict[str, SpectraSet] = {}
        for path in science_files:
            name = Path(path).stem
            logger.info(f"Reducing {name}")
            img, hdr = self._debias(path, bias)
            self.identification_override = self._anchor_path_for_header(hdr)
            img = self._cosmic(img)
            d1c = self.out / 'step1_basic' / 'cosmic_corrected'
            cosmic_out = d1c / f'{name}_cosmic_corrected'
            self._write(Path(str(cosmic_out) + '.fits'), img, hdr)
            self._plot_image(cosmic_out, img, f'{name} cosmic-corrected')
            if self._last_cosmic_mask is not None:
                self._plot_cosmic_marked(d1c / f'{name}_cosmic_marked',
                                         self._last_cosmic_input, self._last_cosmic_mask, name)

            centre, left, right = self._flat_trace
            offset = self._star_offset(img, centre, band)
            initial = centre + offset
            logger.info(f"Initial star offset from flat trace: {offset:+.2f} px")
            trace = self._trace_source(img, band, initial, label=f'Target {name}')
            self._save_trace(d2, name, trace)
            self._plot_trace(d2 / f'{name}_trace', img, trace, f'{name} trace',
                             aperture=self.aperture, sky_gap=self.sky_gap)

            img = img / pixflat
            d4f = self.out / 'step4_flat_corrected'
            self._write(d4f / f'{name}_flat2d_corrected.fits', img, hdr)
            self._plot_trace(d4f / f'{name}_flat2d_corrected', img, trace,
                             f'{name} flat2d-corrected',
                             aperture=self.aperture, sky_gap=self.sky_gap)

            sub, sky = self._subtract_sky(img, trace, band)
            sky_direct = self._extract_sky(img, trace, band)
            self._write(self.out / 'step5_extraction' / f'{name}_skysub.fits', sub, hdr)

            fsum, fopt, eopt, prof = self._extract(sub, trace)
            d4 = self.out / 'step5_extraction'
            self._plot_trace(d4 / f'{name}_skysub', sub, trace,
                             f'{name} after sky subtraction',
                             aperture=self.aperture, sky_gap=self.sky_gap)
            pixel = np.arange(len(fsum), dtype=float)

            wave = pixel + 1.0
            snr_pix = self._estimate_snr(fopt)
            sky_1d = sky[np.arange(len(trace)), np.clip(np.rint(trace).astype(int), 0, sky.shape[1] - 1)]
            self._plot_sky(d4 / f'{name}_sky_pixel', wave, sky_1d, 'pixel', name, sky_direct)
            self._write_sky(d4 / f'{name}_sky_pixel.fits', wave, sky_1d, 'pixel', hdr, sky_direct)
            self._plot_spectrum(d4 / f'{name}_optimal', wave, fopt, 'pixel', name, snr_pix)
            self._write_spectrum(d4 / f'{name}_optimal.fits', wave, fopt, eopt,
                                 sky_1d, 'pixel', hdr, has_snr=snr_pix)
            snr_wave = []
            wave_fit = None
            wave_rms = None
            wave_resolution = None
            if arcs:
                t_sci = self._frame_time(path, hdr)
                m_sci = os.path.getmtime(path)
                _, _, arc_name, arc_flux, arc_hdr = min(
                    arcs, key=lambda a: (abs(a[0] - t_sci), abs(a[1] - m_sci)))
                logger.info(f"Using lamp frame {Path(arc_name).name} for {name} (closest in time)")
                fit, info = arc_solutions.get(arc_name, (None, None))
                if fit is not None:
                    wave = fit(pixel)
                    wave_fit = fit
                    wave_rms = info[3]
                    wave_resolution = info[4]
                    snr_wave = self._estimate_snr(fopt, wave)
                    self._plot_sky(d5 / f'{name}_sky_wavecal', wave, sky_1d, 'Angstrom', name, sky_direct,
                                   tol=wave_rms, resolution=wave_resolution)
                    self._write_sky(d5 / f'{name}_sky_wavecal.fits', wave, sky_1d,
                                    'Angstrom', hdr, sky_direct, wave_fit=wave_fit,
                                    wave_rms=wave_rms)
                    self._plot_spectrum(d5 / f'{name}_optimal_wavecal', wave, fopt,
                                        'Angstrom', name, snr_wave)
                    self._write_spectrum(d5 / f'{name}_optimal_wavecal.fits', wave, fopt,
                                         eopt, sky_1d, 'Angstrom', hdr,
                                         has_snr=snr_pix + snr_wave,
                                         wave_fit=wave_fit, wave_rms=wave_rms,
                                         wave_resolution=wave_resolution)

            exptime = float(hdr.get('EXPTIME', 1.0) or 1.0)
            spec = Spectrum(aperture=1, order=1, wavelength=wave, flux=fopt, error=eopt,
                            background=sky[np.arange(len(trace)), np.clip(np.rint(trace).astype(int), 0, sky.shape[1] - 1)],
                            exptime=exptime)
            sset = SpectraSet()
            sset.add_spectrum(spec)
            results[path] = sset

        logger.info("Slit reduction finished")
        return results

    # 由主平场（或灯谱）一维谱确定统一的有效像素范围。
    def _find_valid_range(self, master_flat_flux, arcs) -> Optional[Tuple[int, int]]:
        """Pixel range where the master flat (else the lamp) 1D spectrum carries signal."""
        ref = master_flat_flux
        if ref is None and arcs:
            ref = np.nanmedian([a[3] for a in arcs], axis=0)
        if ref is None:
            return None
        ref = np.nan_to_num(np.asarray(ref, dtype=float))
        smooth = ndimage.uniform_filter1d(ref, 31)
        thr = self.snr_valid_fraction * np.percentile(smooth, 95)
        good = np.where(smooth > thr)[0]
        if not len(good):
            return None
        i0, i1 = int(good.min()), int(good.max())
        logger.info(f"Valid pixel range (1D output and SNR): {i0 + 1}-{i1 + 1} "
                    f"(> {self.snr_valid_fraction:.1%} of the 95th-percentile flat/lamp flux)")
        return i0, i1

    # 稳健噪声：1.4826 × MAD，等价于高斯分布下的 sigma。
    @staticmethod
    def _mad_sigma(a: np.ndarray) -> float:
        a = np.asarray(a, dtype=float)
        a = a[np.isfinite(a)]
        return float(1.4826 * np.median(np.abs(a - np.median(a)))) if a.size else 0.0

    # 迭代平滑连续谱并在有效范围内固定位置估算信噪比。
    def _estimate_snr(self, flux: np.ndarray, wave: Optional[np.ndarray] = None):
        """SNR at 30%, 50% and 70% of the valid pixel range.

        A smooth continuum (Savitzky-Golay) is fitted to the valid range, the
        residuals are 3-sigma clipped and the fit is repeated until no further
        outliers are rejected; the rejected regions are widened by 3 pixels so
        line wings are excluded from the noise. SNR = median of the smooth continuum in a small window around
        the position / noise, where the noise is the RMS of the clipped residuals
        (whole range, or the local window when snr_noise_mode = local).
        """
        flux = np.asarray(flux, dtype=float)
        if flux.ndim != 1 or not len(flux):
            return []
        i0, i1 = self._valid_range if self._valid_range else (0, len(flux) - 1)
        i0, i1 = max(i0, 0), min(i1, len(flux) - 1)
        seg = flux[i0:i1 + 1]
        n = len(seg)
        if n < 20:
            return []
        window = min(self.snr_smooth_window, n if n % 2 else n - 1)
        window = max(window - (window % 2 == 0), 7)
        x = np.arange(n, dtype=float)
        keep = np.isfinite(seg)
        smooth = np.zeros(n)
        res = np.zeros(n)
        for _ in range(20):
            filled = np.interp(x, x[keep], seg[keep])
            smooth = savgol_filter(filled, window, 3, mode='interp')
            res = seg - smooth
            sigma = self._mad_sigma(res[keep])
            new_keep = np.isfinite(seg) & (np.abs(res) <= 3 * sigma)
            if (new_keep == keep).all() or new_keep.sum() < 10:
                keep = new_keep if new_keep.sum() >= 10 else keep
                break
            keep = new_keep
        # Widen the rejected regions so the wings of broad lines are excluded too.
        keep = ndimage.binary_erosion(keep, iterations=3, border_value=1)
        rms_all = self._mad_sigma(res[keep])
        half = self.snr_half_window
        out = []
        for frac in (0.3, 0.5, 0.7):
            k = int(round(frac * (n - 1)))
            c = i0 + k
            a, b = max(k - half, 0), min(k + half + 1, n)
            signal = float(np.median(smooth[a:b]))
            noise = rms_all
            if self.snr_noise_mode == 'local':
                loc = res[a:b][keep[a:b]]
                if loc.size > 10:
                    noise = self._mad_sigma(loc)
            # Continuum/line-independent cross-check (DER_SNR style second differences)
            w = seg[a:b]
            noise_der = float(1.482602 / np.sqrt(6.0) * np.nanmedian(
                np.abs(2.0 * w[2:-2] - w[:-4] - w[4:]))) if len(w) > 8 else np.nan
            if np.isfinite(noise_der) and noise_der > 0 and noise > 1.5 * noise_der:
                logger.warning(
                    f"SNR window at {frac:.0%}: residual sigma {noise:.1f} is >1.5x the "
                    f"line-insensitive estimate {noise_der:.1f}; a broad line may be inflating the noise")
            snr = signal / noise if noise > 0 else np.nan
            entry = {'frac': frac, 'index': c, 'signal': signal, 'noise': noise,
                     'noise_der': noise_der, 'snr': float(snr), 'unit': 'pixel',
                     'x': float(wave[c]) if wave is not None else float(c + 1)}
            if wave is not None:
                entry['unit'] = 'Angstrom'
            out.append(entry)
        txt = ', '.join(f"{int(e['frac'] * 100)}%@{e['x']:.1f}: {e['snr']:.1f}" for e in out)
        logger.info(f"SNR at {out[0]['unit']} positions ({self.snr_noise_mode} noise): {txt}")
        return out

    # 读取观测时间并回退到文件时间。
    def _frame_time(self, path: str, hdr) -> float:
        key = self.config.get('data', 'statime_key', 'DATE-OBS')
        try:
            from astropy.time import Time
            return float(Time(hdr[key]).unix)
        except Exception:
            return float(os.path.getmtime(path))

    # 沿主平场轨迹、按整个狭缝宽度提取标定光谱。
    def _extract_calib(self, img: np.ndarray, hdr, name: str,
                       trace=None, plot: bool = True) -> np.ndarray:
        """Extract a flat/lamp frame over the full slit width (no sky subtraction)."""
        if trace is None:
            trace = self._flat_trace[0]
        _, left, right = self._flat_trace
        half_width = int(np.ceil(np.median(right - left) / 2.0))
        fsum, fopt, eopt, _ = self._extract(img, trace, half_width=half_width)
        if plot:
            self._plot_trace(self.out / 'step2_trace' / f'{name}_trace', img, trace,
                             f'{name} trace')
        return fopt

    # 有效范围确定后再输出标定光谱的像素谱。
    def _write_calib_pixel(self, name: str, flux: np.ndarray, hdr):
        d = self.out / 'step5_extraction'
        x = np.arange(len(flux)) + 1.0
        self._write_spectrum(d / f'{name}_optimal.fits', x, flux,
                             np.zeros_like(flux), np.zeros_like(flux), 'pixel', hdr)
        self._plot_spectrum(d / f'{name}_optimal', x, flux, 'pixel', name)

    # 按统一有效像素范围裁剪一维数组（范围之外不输出）。
    def _crop(self, *arrays):
        out = []
        for a in arrays:
            if a is not None and self._valid_range and len(a) > self._valid_range[1]:
                a = np.asarray(a)[self._valid_range[0]:self._valid_range[1] + 1]
            out.append(a)
        return out if len(out) > 1 else out[0]

    # 将光谱、背景和信噪比写入FITS。
    def _write_spectrum(self, path, wave, flux, err, sky, unit, hdr, has_snr=(),
                        wave_fit=None, wave_rms=None, wave_resolution=None):
        path.parent.mkdir(parents=True, exist_ok=True)
        wave, flux, err, sky = self._crop(wave, flux, err, sky)
        cols = [fits.Column(name='WAVELENGTH', format='D', array=wave),
                fits.Column(name='FLUX_OPT', format='D', array=flux),
                fits.Column(name='ERROR', format='D', array=np.nan_to_num(err)),
                fits.Column(name='SKY', format='D', array=sky)]
        t = fits.BinTableHDU.from_columns(cols)
        t.header['WAVEUNIT'] = unit
        if wave_fit is not None:
            t.header['WAVEFIT'] = dispersion_type(wave_fit)
            t.header['WAVEDEG'] = wave_fit.degree()
            if wave_rms is not None:
                t.header['WAVERMS'] = (round(float(wave_rms), 4),
                                       'wavelength fit RMS (A)')
            if wave_resolution:
                t.header['WAVERES'] = (round(float(wave_resolution['resolution']), 1),
                                       'resolving power R = lambda / FWHM')
                t.header['WAVEFWHM'] = (round(float(wave_resolution['fwhm_ang']), 4),
                                        'arc-line FWHM (A)')
        for e in has_snr:
            key = ('SNRP' if e['unit'] == 'pixel' else 'SNRW') + str(int(e['frac'] * 100))
            t.header[key] = (round(e['snr'], 2), f"SNR at {int(e['frac'] * 100)}% of range")
        t.header['EXPTIME'] = float(hdr.get('EXPTIME', 0) or 0)
        t.header['OBJECT'] = str(hdr.get('OBJECT', path.stem.split('_')[0]))
        fits.HDUList([fits.PrimaryHDU(), t]).writeto(path, overwrite=True)

    # 绘制求和与最优提取光谱。
    def _plot_spectrum(self, path, wave, flux, unit, title, snr=()):
        if not self.save_plots:
            return
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        wave, flux = self._crop(wave, flux)
        fig, ax = plt.subplots(figsize=(10, 4))
        ax.plot(wave, flux, lw=0.6, color='b')
        per = 'pixel' if unit == 'pixel' else 'A'
        ymax = np.nanmax(flux)
        for k, e in enumerate(snr):
            ax.axvline(e['x'], color='r', ls='--', lw=0.6)
            ax.text(e['x'], 0.96 - 0.08 * (k % 2), f"SNR={e['snr']:.1f}@{e['x']:.0f} {per}", transform=ax.get_xaxis_transform(),
                    ha='center', va='top', color='r', fontsize=9,
                    bbox=dict(fc='w', ec='r', lw=0.5, alpha=0.9))
        ax.set_ylim(top=ymax * 1.25)
        ax.set_xlabel('Wavelength (A)' if unit == 'Angstrom' else 'Pixel')
        ax.set_ylabel('Counts')
        ax.set_title(title)
        fig.savefig(Path(str(path) + '.pdf'), dpi=110, bbox_inches='tight')
        plt.close(fig)

    # 读取天光线表（兼容 .dat/.txt 后缀）。
    @staticmethod
    def _read_line_list(entry: str):
        path = _resolve_path(entry)
        if path is None:
            alt = str(Path(entry).with_suffix('.txt' if entry.endswith('.dat') else '.dat'))
            path = _resolve_path(alt)
        if path is None:
            logger.warning(f"Sky line list not found: {entry}")
            return []
        rows = []
        for line in path.read_text(errors='ignore').splitlines():
            parts = line.split()
            if not parts or parts[0].startswith('#'):
                continue
            try:
                rows.append((float(parts[0]), ' '.join(parts[1:])))
            except ValueError:
                continue
        rows.sort()
        return rows

    # 低分辨率OH线表用于匹配，skylines.dat 用于确定线的名称。
    def _load_sky_lines(self):
        if getattr(self, '_sky_lines', None) is None:
            low = np.array([r[0] for r in self._read_line_list(self.sky_lowres_list)])
            ident = self._read_line_list(self.sky_id_list)
            self._sky_lines = (low, np.array([r[0] for r in ident]), [r[1] for r in ident])
        return self._sky_lines

    # 绘制天空光谱；波长定标后用天光线表标记发射线。
    def _plot_sky(self, path, x, sky, unit, title, sky_direct=None,
                  tol=None, resolution=None):
        """Night-sky spectrum. In wavelength space the emission lines are matched to
        the sky line lists and marked like the anchor lines of the arc diagnostic."""
        if not self.save_plots:
            return
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        x, sky, sky_direct = self._crop(np.asarray(x, dtype=float), sky, sky_direct)
        fig, ax = plt.subplots(figsize=(10, 4.5))
        ax.plot(x, sky, lw=0.5, color='k', label='model at trace position')
        ref = sky
        if sky_direct is not None:
            ax.plot(x, sky_direct, lw=0.5, color='tab:orange', alpha=0.8,
                    label='extracted in sky region')
            ref = sky_direct
        ref = np.nan_to_num(np.asarray(ref, dtype=float))
        n_match = 0
        if unit == 'Angstrom':
            low, id_w, id_lab = self._load_sky_lines()
            tol = float(tol) if tol and np.isfinite(tol) and tol > 0 else self.sky_match_tol
            r_med = float(resolution['resolution']) if resolution else 0.0
            # tolerance = FWHM/2 (the lamp-fit RMS is far below the line width); RMS only if FWHM unknown
            if resolution and np.isfinite(resolution.get('fwhm_ang', np.nan)):
                tol = 0.5 * float(resolution['fwhm_ang'])
            use_low = r_med < self.sky_resolution_split
            is_oh = np.array([bool(re.match(r'^\d+-\d+', l)) for l in id_lab], dtype=bool)
            peaks = self._find_peaks(ref, self.peak_snr)
            idx = np.arange(len(ref))
            span = max(float(np.nanmax(ref) - np.nanmin(ref)), 1.0)
            for c, _amp in peaks:
                w = float(np.interp(c, idx, x))
                label, colour = f'{w:.1f}', '0.5'
                if use_low:
                    # low/medium resolution: OH blends -> lowskylines; named via skylines.
                    # Non-OH lines (Na, Hg, [OI], ...) are matched separately.
                    j = int(np.argmin(np.abs(low - w))) if len(low) else -1
                    oh_d = abs(low[j] - w) if j >= 0 else np.inf
                    cand = np.where(~is_oh & (np.abs(id_w - w) <= tol))[0] if len(id_w) else []
                    k = int(cand[np.argmin(np.abs(id_w[cand] - w))]) if len(cand) else -1
                    if k >= 0 and abs(id_w[k] - w) < oh_d:
                        label, colour = f'{id_w[k]:.2f} {id_lab[k]}'.strip(), 'tab:blue'
                        n_match += 1
                    elif oh_d <= tol:
                        name = ''
                        if len(id_w):
                            m = int(np.argmin(np.abs(id_w - low[j])))
                            if abs(id_w[m] - low[j]) <= self.sky_id_tol:
                                name = id_lab[m]
                        label, colour = f'{low[j]:.2f} OH {name}'.strip(), 'tab:red'
                        n_match += 1
                    elif len(id_w):
                        # not in the OH list: fall back to the full skylines list
                        m = int(np.argmin(np.abs(id_w - w)))
                        if abs(id_w[m] - w) <= tol:
                            label = f'{id_w[m]:.2f} {id_lab[m]}'.strip()
                            colour = 'tab:red' if is_oh[m] else 'tab:blue'
                            n_match += 1
                elif len(id_w):
                    # higher resolution: compare directly with skylines
                    m = int(np.argmin(np.abs(id_w - w)))
                    if abs(id_w[m] - w) <= tol:
                        label = f'{id_w[m]:.2f} {id_lab[m]}'.strip()
                        colour = 'tab:red' if is_oh[m] else 'tab:blue'
                        n_match += 1
                yv = float(np.interp(w, x, ref))
                ax.plot([w, w], [yv + 0.02 * span, yv + 0.07 * span], color=colour, lw=0.85,
                        solid_capstyle='butt', zorder=4)
                ax.text(w, yv + 0.08 * span, label, rotation=90, ha='center', va='bottom',
                        fontsize=5, color=colour, zorder=5)
            ax.set_ylim(top=float(np.nanmax(ref)) + 0.45 * span)
            title = (f'Night-sky spectrum ({n_match}/{len(peaks)} peaks matched, tol={tol:.2f} A, '
                     f'R~{r_med:.0f}: {"lowskylines + non-OH" if use_low else "skylines"}; red = OH, blue = other)')
        else:
            title = 'Night-sky spectrum'
        ax.legend(loc='upper right', fontsize=7)
        ax.set_xlabel('Wavelength (A)' if unit == 'Angstrom' else 'Pixel')
        ax.set_ylabel('Sky counts / spatial pixel')
        ax.set_title(title)
        fig.savefig(Path(str(path) + '.pdf'), dpi=110, bbox_inches='tight')
        plt.close(fig)

    def _write_sky(self, path, wave, sky, unit, hdr, sky_direct=None, wave_fit=None, wave_rms=None):
        """Write 1D sky spectra: SKY_MODEL (fit at the trace position) and SKY_REGION
        (extracted directly in the sky region), counts per spatial pixel."""
        path.parent.mkdir(parents=True, exist_ok=True)
        wave, sky, sky_direct = self._crop(wave, sky, sky_direct)
        t = fits.BinTableHDU.from_columns([
            fits.Column(name='WAVELENGTH', format='D', array=wave),
            fits.Column(name='SKY_MODEL', format='D', array=sky),
            fits.Column(name='SKY_REGION', format='D',
                        array=np.nan_to_num(sky_direct if sky_direct is not None else sky))])
        t.header['WAVEUNIT'] = unit
        if wave_fit is not None:
            t.header['WAVEFIT'] = dispersion_type(wave_fit)
            t.header['WAVEDEG'] = wave_fit.degree()
            if wave_rms is not None:
                t.header['WAVERMS'] = (round(float(wave_rms), 4), 'wavelength fit RMS (A)')
        t.header['EXPTIME'] = float(hdr.get('EXPTIME', 0) or 0)
        fits.HDUList([fits.PrimaryHDU(), t]).writeto(path, overwrite=True)

    def _write_lamp_wave(self, d: Path, name: str, arc_flux: np.ndarray,
                         hdr, info, fit, title='lamp'):
        """Write the wavelength-calibrated lamp spectrum (once per arc frame)."""
        key = name
        if key in self._written_arc_wave:
            return
        self._written_arc_wave.add(key)
        wave = fit(np.arange(len(arc_flux), dtype=float))
        self._plot_spectrum(d / f'{name}_optimal_wavecal', wave, arc_flux,
                            'Angstrom', f'{name} ({title})')
        self._write_spectrum(d / f'{name}_optimal_wavecal.fits', wave, arc_flux,
                             np.zeros_like(arc_flux), np.zeros_like(arc_flux),
                             'Angstrom', hdr,
                             wave_fit=fit, wave_rms=info[3],
                             wave_resolution=info[4])
        logger.info("Wrote wavelength-calibrated lamp spectrum: %s",
                    d / f'{name}_optimal_wavecal.fits')

    # 绘制标定光谱及波长拟合残差。
    def _plot_arc(self, path, arc, info, fit):
        if not self.save_plots:
            return
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages

        x, y, keep, rms, resolution = info
        wave, arc = self._crop(fit(np.arange(len(arc), dtype=float)), arc)
        wlo, whi = float(np.min(wave)), float(np.max(wave))
        res = y - fit(x)
        fig, (a1, a2) = plt.subplots(
            2, 1, figsize=(10, 7), sharex=True,
            gridspec_kw={'height_ratios': [3.2, 1.2], 'hspace': 0.0})
        a1.plot(wave, arc, color='0.35', lw=0.6, label='arc spectrum')
        arc_span = max(float(np.nanmax(arc) - np.nanmin(arc)), 1.0)
        for w, ok in zip(y, keep):
            if not wlo <= w <= whi:
                continue
            # np.interp needs increasing x; wavelength may fall with pixel.
            order = np.argsort(wave)
            px = float(np.interp(w, wave[order], order.astype(float)))
            centre = int(round(px))
            lo, hi = max(centre - 3, 0), min(centre + 4, len(arc))
            yv = float(np.nanmax(arc[lo:hi]))
            colour = 'tab:red' if ok else '0.55'
            tick_bottom = yv + 0.025 * arc_span
            tick_top = yv + 0.085 * arc_span
            a1.plot([w, w], [tick_bottom, tick_top],
                    color=colour, lw=0.85, solid_capstyle='butt', zorder=4)
            a1.text(w, tick_top + 0.015 * arc_span, f'{w:.4f}', rotation=90,
                    ha='center', va='bottom', fontsize=6.0, color=colour,
                    zorder=5, clip_on=True)
        a1.set_xlim(wlo, whi)
        a1.set_ylim(top=float(np.nanmax(arc)) + 0.30 * arc_span)
        a1.set_ylabel('Counts')
        a1.set_title(f'Arc wavelength solution: {keep.sum()}/{len(y)} matched full-catalog '
                     f'lines, rms={rms:.3f} A '
                     f'({dispersion_type(fit)} degree {fit.degree()})')
        a1.legend(loc='upper right', fontsize=7)
        a2.plot(fit(x)[keep], res[keep], 'o', color='#1565c0', ms=3.5,
                label='inlier')
        if (~keep).any():
            a2.plot(fit(x)[~keep], res[~keep], 'o', mfc='none',
                    mec='#1565c0', ms=4.0, label='rejected')
        a2.axhline(0, color='0.35', lw=0.8)
        a2.axhline(rms, color='0.45', lw=0.7, ls=':')
        a2.axhline(-rms, color='0.45', lw=0.7, ls=':')
        a2.set_ylabel('Residual (A)')
        a2.set_xlabel('Wavelength (A)')
        a2.legend(loc='best', fontsize=7)

        fig2, a3 = plt.subplots(figsize=(10, 5))
        a3.set_ylabel('R = λ/Δλ')
        a3.set_xlabel('Wavelength (A)')
        a3.grid(alpha=0.2)
        if resolution and len(resolution.get('waves', [])):
            if len(resolution.get('rejected_waves', [])):
                a3.plot(resolution['rejected_waves'],
                        resolution['rejected_resolutions'], 'x', ms=3.5,
                        color='0.55', mew=0.8, label='rejected')
            a3.plot(resolution['waves'], resolution['resolutions'], 'o',
                    ms=3.5, color='#1565c0', label='R(λ)')
            tw = resolution.get('trend_waves', [])
            tr = resolution.get('trend_resolutions', [])
            if len(tw) and len(tr):
                a3.plot(tw, tr, color='#1565c0', lw=1.0, ls='--',
                        label='trend')
            else:
                a3.axhline(resolution['resolution'], color='#1565c0',
                           lw=0.8, ls='--', label='median')
            a3.text(0.995, 0.95,
                    f'median R ~ {resolution["resolution"]:.0f} '
                    f'(n={resolution["n"]}, rejected {resolution.get("n_rejected", 0)}, '
                    'all detected peaks)',
                    transform=a3.transAxes, ha='right', va='top', fontsize=7,
                    color='#1565c0')
            a3.legend(loc='best', fontsize=7)
        else:
            a3.text(0.5, 0.5, 'No clean lines for R(lambda)',
                    transform=a3.transAxes, ha='center', va='center',
                    fontsize=8, color='0.45')
        a3.set_xlim(a1.get_xlim())
        out = Path(str(path) + '.pdf')
        out.parent.mkdir(parents=True, exist_ok=True)
        with PdfPages(out) as pdf:
            pdf.savefig(fig, bbox_inches='tight')
            pdf.savefig(fig2, bbox_inches='tight')
        plt.close(fig)
        plt.close(fig2)
