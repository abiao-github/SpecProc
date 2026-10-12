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

import csv
import os
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from astropy.io import fits
from scipy import ndimage

from specproc.core.data_structures import Spectrum, SpectraSet

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
        self.fig_format = c.get('reduce', 'fig_format', 'png')
        self.gain = c.get_float('telescope.detector', 'gain', 1.0)
        self.readnoise = c.get_float('telescope.detector', 'readnoise', 5.0)

        self.overscan_start = c.get_int('data', 'overscan_start_column', -1)
        self.dispersion_axis = c.get('slit', 'dispersion_axis', 'y').strip().lower()
        self.band_x = (c.get_int('slit', 'band_start', -1), c.get_int('slit', 'band_end', -1))
        self.band_frac = c.get_float('slit', 'band_threshold', 0.05)
        self.flat_smooth = c.get_int('slit', 'flat_smooth_pixels', 21)
        self.flat_edge_margin = c.get_int('slit', 'flat_edge_margin', 40)
        self.trace_order = c.get_int('slit', 'trace_order', 2)
        self.trace_nbins = c.get_int('slit', 'trace_bins', 40)
        self.aperture = c.get_int('slit', 'extract_half_width', 12)
        self.sky_gap = c.get_int('slit', 'sky_inner_offset', 20)
        self.aperture_mode = c.get('slit', 'aperture_mode', 'fixed').strip().lower()
        self.aperture_fwhm_factor = c.get_float('slit', 'aperture_fwhm_factor', 1.7)
        self.sky_fwhm_factor = c.get_float('slit', 'sky_fwhm_factor', 3.0)
        self.sky_order = c.get_int('slit', 'sky_order', 1)
        self.wave_order = c.get_int('slit', 'wave_poly_order', 3)
        self.wave_hint = (c.get_float('slit', 'wave_start', 0.0), c.get_float('slit', 'wave_end', 0.0))
        self.match_tol = c.get_float('slit', 'match_tolerance', 3.0)
        self.snr_half_window = c.get_int('slit', 'snr_half_window', 50)
        self.snr_flux_frac = c.get_float('slit', 'snr_effective_fraction', 0.1)
        self.cosmic = c.get_bool('reduce', 'cosmic_enabled', True)

    # ------------------------------------------------------------------ I/O
    # 读取图像并处理方向与过扫描区域。
    def _read(self, path: str) -> Tuple[np.ndarray, fits.Header]:
        with fits.open(path) as hdul:
            data = np.asarray(hdul[0].data, dtype=np.float64)
            header = hdul[0].header.copy()
        if self.overscan_start > 0 and data.shape[1] >= self.overscan_start:
            # Skip the first overscan columns (they contain image bleed)
            ov = data[:, self.overscan_start + 3:]
            level = np.median(ov, axis=1) if ov.size else 0.0
            data = data[:, :self.overscan_start - 1] - np.asarray(level)[:, None]
        if self.dispersion_axis == 'x':
            data = data.T
        return data, header

    # 保存浮点图像及其头信息。
    def _write(self, path: Path, data: np.ndarray, header: Optional[fits.Header] = None):
        path.parent.mkdir(parents=True, exist_ok=True)
        fits.PrimaryHDU(np.asarray(data, dtype=np.float32), header).writeto(path, overwrite=True)

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
        ax.imshow(data, origin='lower', aspect='auto', cmap='gray', vmin=lo, vmax=hi)
        if overlay:
            for xs, ys, style in overlay:
                ax.plot(xs, ys, style, lw=0.8)
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
        self._write(self.out / 'step1_basic' / 'bias_subtracted' / 'master_bias.fits', master)
        logger.info(f"Master bias from {len(bias_files)} frames: median={np.median(master):.2f}")
        return master

    # 扣除主偏置并保存校正图像。
    def _debias(self, path: str, bias: np.ndarray):
        data, hdr = self._read(path)
        data = data - bias
        self._write(self.out / 'step1_basic' / 'bias_subtracted' / f'{Path(path).stem}.fits', data, hdr)
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
        d2 = self.out / 'step2_trace'
        centre, left, right, (ty, tx) = self._trace_flat(master, (x0, x1))
        yy = np.arange(master.shape[0])
        self._plot_image(d2 / 'master_flat_trace', master, 'Master flat trace',
                         [(centre, yy, 'r-'), (left, yy, 'y-'), (right, yy, 'y-'), (tx, ty, 'c.')])
        np.savetxt(d2 / 'master_flat_trace.txt', np.column_stack([yy + 1, centre + 1, left + 1, right + 1]),
                   header='dispersion_pixel  centre  left_edge  right_edge  (spatial columns, 1-based)',
                   fmt='%.3f')
        self._flat_trace = (centre, left, right)
        self._master_flat = master
        return pix, (x0, x1)

    # --------------------------------------------------------------- stage 3
    # 探测并校正图像中的宇宙射线。
    def _cosmic(self, data: np.ndarray) -> np.ndarray:
        if self.cosmic:
            try:
                import astroscrappy
                mask, clean = astroscrappy.detect_cosmics(
                    data + 0.0, gain=self.gain, readnoise=self.readnoise,
                    sigclip=self.config.get_float('reduce', 'cosmic_sigclip', 5.0),
                    objlim=self.config.get_float('reduce', 'cosmic_objlim', 5.0),
                    satlevel=1e9, cleantype='meanmask')
                logger.info(f"Cosmic rays flagged: {int(mask.sum())} px")
                data = clean
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
        g = lambda x, a, m, sg: a * np.exp(-0.5 * ((x - m) / sg) ** 2)
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
            data, _ = self._read(f)
            f_i = self._measure_fwhm(data - bias, band)
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
    def _subtract_sky(self, img: np.ndarray, trace: np.ndarray, band: Tuple[int, int]):
        # Only the flat-corrected region is used (edge margin excluded)
        x0, x1 = band
        m = min(self.flat_edge_margin, (x1 - x0) // 4)
        x0, x1 = x0 + m, x1 - m
        ny, nx = img.shape
        xi = np.arange(nx)
        sky = np.zeros_like(img)
        for y in range(ny):
            m = (xi >= x0) & (xi < x1) & (np.abs(xi - trace[y]) > self.sky_gap)
            if m.sum() < 10:
                continue
            xm, vm = xi[m], img[y, m]
            c = np.polyfit(xm, vm, self.sky_order)
            r = vm - np.polyval(c, xm)
            s = 1.4826 * np.median(np.abs(r)) + 1e-6
            k = np.abs(r) < 3 * s
            if k.sum() > self.sky_order + 2:
                c = np.polyfit(xm[k], vm[k], self.sky_order)
            sky[y, x0:x1] = np.polyval(c, xi[x0:x1])
        return img - sky, sky

    # --------------------------------------------------------------- stage 5
    # 沿狭缝轨迹计算求和与最优光谱。
    def _extract(self, img: np.ndarray, trace: np.ndarray, sky_sub_noise=None):
        ny, nx = img.shape
        hw = self.aperture
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
    def _find_peaks(spec: np.ndarray):
        s = spec - ndimage.minimum_filter1d(ndimage.gaussian_filter1d(spec, 3), 61)
        s = ndimage.gaussian_filter1d(s, 1.0)
        noise = 1.4826 * np.median(np.abs(s - np.median(s))) + 1e-9
        pk = [i for i in range(2, len(s) - 2)
              if s[i] > s[i - 1] and s[i] >= s[i + 1] and s[i] > 8 * noise and s[i] > s[i - 2]]
        out = []
        for i in pk:
            lo, hi = max(i - 3, 0), min(i + 4, len(s))
            w = np.clip(s[lo:hi], 0, None)
            out.append((float((np.arange(lo, hi) * w).sum() / w.sum()), float(s[i])))
        return out

    # 读取并排序参考波长目录。
    def _load_catalog(self) -> np.ndarray:
        path = _resolve_path(self.config.get('telescope.linelist', 'full_linelist', ''))
        if path is None:
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
    def _load_id_file(self) -> Optional[np.ndarray]:
        path = _resolve_path(self.config.get('slit', 'identification_file', ''))
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
        if len(rows) < self.wave_order + 1:
            logger.warning(f"Identification file {path} has too few lines")
            return None
        logger.info(f"Using {len(rows)} identified lines from {path}")
        return np.array(rows)

    # 依据波长提示匹配谱线与目录。
    def _auto_identify(self, peaks, npix, catalog):
        """Match peaks to catalogue assuming a linear solution from wave_start/wave_end hints."""
        w0, w1 = self.wave_hint
        if w0 <= 0 or w1 <= 0 or len(catalog) == 0:
            return None
        pix = np.array([p[0] for p in peaks])
        coef = np.polyfit([0, npix - 1], [w0, w1], 1)
        pairs = None
        for it in range(8):
            wl = np.polyval(coef, pix)
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
            coef = np.polyfit(pairs[:, 0], pairs[:, 1], order)
            if order < self.wave_order:
                coef = np.concatenate([np.zeros(self.wave_order - order), coef])
        return pairs

    # 拟合并评估标定帧的波长解。
    def _wavelength_solution(self, arc_flux: np.ndarray, d: Path, tag: str = 'arc'):
        npix = len(arc_flux)
        d.mkdir(parents=True, exist_ok=True)
        peaks = self._find_peaks(arc_flux)
        with open(d / f'{tag}_peaks.csv', 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['pixel', 'flux'])
            w.writerows(peaks)
        logger.info(f"Arc: {len(peaks)} peaks written to {d / f'{tag}_peaks.csv'}")

        pairs = self._load_id_file()
        source = 'identification_file'
        if pairs is None:
            pairs = self._auto_identify(peaks, npix, self._load_catalog())
            source = 'auto_match'
        if pairs is None or len(pairs) < self.wave_order + 2:
            logger.warning("No wavelength solution: provide [slit] identification_file or "
                           "wave_start/wave_end hints. Spectra are saved on a pixel axis.")
            return None, None, source

        x, y = pairs[:, 0], pairs[:, 1]
        keep = np.ones(len(x), bool)
        for _ in range(5):
            coef = np.polyfit(x[keep], y[keep], self.wave_order)
            res = y - np.polyval(coef, x)
            rms = np.std(res[keep])
            new = np.abs(res) < max(3 * rms, 0.3)
            if new.sum() < self.wave_order + 2 or (new == keep).all():
                break
            keep = new
        coef = np.polyfit(x[keep], y[keep], self.wave_order)
        rms = np.std((y - np.polyval(coef, x))[keep])
        logger.info(f"Wavelength solution ({source}): {keep.sum()}/{len(x)} lines, rms={rms:.3f} A")
        return coef, (x, y, keep, rms), source

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
            self._write(self.out / 'step1_basic' / 'overscan_corrected' / f'{Path(f).stem}.fits', data, hdr)

        bias = self._master_bias(bias_files)
        if flat_files:
            pixflat, band = self._make_flat(flat_files, bias)
        else:
            raise RuntimeError("Slit mode needs flat frames: the trace is derived from the master flat.")

        self._set_aperture(science_files, bias, band)

        # Lamp frames are never combined: each science frame uses the single
        # lamp frame closest in time (DATE-OBS, then file mtime as tie-break).
        arcs = []
        for f in arc_files:
            img_a, hdr_a = self._debias(f, bias)
            arcs.append((self._frame_time(f, hdr_a), os.path.getmtime(f), f,
                         self._extract_calib(img_a, hdr_a, Path(f).stem)))
        if flat_files:
            self._extract_calib(self._master_flat, fits.Header(), 'master_flat')

        results: Dict[str, SpectraSet] = {}
        for path in science_files:
            name = Path(path).stem
            logger.info(f"Reducing {name}")
            img, hdr = self._debias(path, bias)
            img = self._cosmic(img)
            self._write(self.out / 'step1_basic' / 'cosmic_corrected' / f'{name}.fits', img, hdr)

            centre, left, right = self._flat_trace
            offset = self._star_offset(img, centre, band)
            trace = centre + offset
            logger.info(f"Star offset from flat trace: {offset:+.2f} px")
            d2 = self.out / 'step2_trace'
            yy = np.arange(len(trace))
            self._plot_image(d2 / f'{name}_trace', img, f'{name} trace (master flat + {offset:+.1f} px)',
                             [(trace, yy, 'r-'), (left, yy, 'y-'), (right, yy, 'y-')])
            np.savetxt(d2 / f'{name}_trace.txt', np.column_stack([yy + 1, trace + 1]),
                       header=f'dispersion_pixel  spatial_column (1-based); master-flat centre {offset:+.3f} px',
                       fmt='%.3f')

            img = img / pixflat
            d4f = self.out / 'step4_flat_corrected'
            self._write(d4f / f'{name}.fits', img, hdr)

            sub, sky = self._subtract_sky(img, trace, band)
            self._write(self.out / 'step5_extraction' / f'{name}_skysub.fits', sub, hdr)

            fsum, fopt, eopt, prof = self._extract(sub, trace)
            d4 = self.out / 'step5_extraction'
            d5 = self.out / 'step7_wavelength'
            pixel = np.arange(len(fsum), dtype=float)

            wave = pixel + 1.0
            wave_unit = 'pixel'
            snr_pix = self._estimate_snr(fopt)
            sky_1d = sky[np.arange(len(trace)), np.clip(np.rint(trace).astype(int), 0, sky.shape[1] - 1)]
            self._plot_sky(d4 / f'{name}_sky_pixel', wave, sky_1d, 'pixel', name)
            self._plot_spectrum(d4 / f'{name}_spec_pixel', wave, fsum, fopt, 'pixel', name, snr_pix)
            snr_wave = []
            if arcs:
                t_sci = self._frame_time(path, hdr)
                m_sci = os.path.getmtime(path)
                _, _, arc_name, arc_flux = min(
                    arcs, key=lambda a: (abs(a[0] - t_sci), abs(a[1] - m_sci)))
                logger.info(f"Using lamp frame {Path(arc_name).name} for {name} (closest in time)")
                coef, info, src = self._wavelength_solution(arc_flux, d5, Path(arc_name).stem)
                if coef is not None:
                    wave = np.polyval(coef, pixel)
                    wave_unit = 'Angstrom'
                    self._plot_arc(d5 / f'{name}_arc', arc_flux, info, coef)
                    snr_wave = self._estimate_snr(fopt, wave)
                    self._plot_sky(d5 / f'{name}_sky_wave', wave, sky_1d, 'Angstrom', name)
                    self._plot_spectrum(d5 / f'{name}_spec_wave', wave, fsum, fopt,
                                        'Angstrom', name, snr_wave)

            exptime = float(hdr.get('EXPTIME', 1.0) or 1.0)
            spec = Spectrum(aperture=1, order=1, wavelength=wave, flux=fopt, error=eopt,
                            background=sky[np.arange(len(trace)), np.clip(np.rint(trace).astype(int), 0, sky.shape[1] - 1)],
                            exptime=exptime)
            sset = SpectraSet()
            sset.add_spectrum(spec)
            results[path] = sset

            self._write_spectrum(d4 / f'{name}_spec.fits', wave, fsum, fopt, eopt, spec.background, wave_unit, hdr,
                                has_snr=snr_pix + snr_wave)
        logger.info("Slit reduction finished")
        return results

    # 估算指定波段位置的光谱信噪比。
    def _estimate_snr(self, flux: np.ndarray, wave: Optional[np.ndarray] = None):
        """SNR at 25/50/75% of the effective spectral range.

        The effective range is the contiguous span around the flux peak where the smoothed flux exceeds
        snr_effective_fraction of its maximum. In each window (+-snr_half_window
        px) the signal is the median flux and the noise is the standard deviation
        of the residual to a quadratic fit (dimensionless ratio). `x` is the
        position: pixel (1-based) without `wave`, wavelength in Angstrom with it.
        """
        sm = ndimage.uniform_filter1d(np.nan_to_num(flux), 51)
        above = sm > self.snr_flux_frac * sm.max()
        pk = int(np.argmax(sm))
        lo = hi = pk
        while lo > 0 and above[lo - 1]:
            lo -= 1
        while hi < len(sm) - 1 and above[hi + 1]:
            hi += 1
        if hi - lo < 20:
            return []
        out = []
        for frac in (0.25, 0.5, 0.75):
            c = int(round(lo + frac * (hi - lo)))
            a, b = max(c - self.snr_half_window, 0), min(c + self.snr_half_window + 1, len(flux))
            seg = np.asarray(flux[a:b], dtype=float)
            x = np.arange(len(seg))
            keep = np.isfinite(seg)
            for _ in range(3):
                coef = np.polyfit(x[keep], seg[keep], 2)
                res = seg - np.polyval(coef, x)
                sig = 1.4826 * np.median(np.abs(res[keep])) + 1e-12
                keep = np.isfinite(seg) & (np.abs(res) < 4 * sig)
            noise = np.std(res[keep])
            signal = float(np.median(seg[keep]))
            snr = signal / noise if noise > 0 else np.nan
            entry = {'frac': frac, 'index': c, 'signal': signal, 'noise': float(noise),
                     'snr': float(snr), 'unit': 'pixel'}
            entry['x'] = float(wave[c]) if wave is not None else float(c + 1)
            if wave is not None:
                entry['unit'] = 'Angstrom'
            out.append(entry)
        txt = ', '.join(f"{int(e['frac'] * 100)}%@{e['x']:.1f}: {e['snr']:.1f}" for e in out)
        logger.info(f"SNR at {out[0]['unit']} positions: {txt}")
        return out

    # 读取观测时间并回退到文件时间。
    def _frame_time(self, path: str, hdr) -> float:
        key = self.config.get('data', 'statime_key', 'DATE-OBS')
        try:
            from astropy.time import Time
            return float(Time(hdr[key]).unix)
        except Exception:
            return float(os.path.getmtime(path))

    # 沿主平场轨迹提取标定光谱。
    def _extract_calib(self, img: np.ndarray, hdr, name: str) -> np.ndarray:
        """Extract a flat/lamp frame along the master-flat centre line (no sky subtraction)."""
        centre = self._flat_trace[0]
        fsum, fopt, eopt, _ = self._extract(img, centre)
        d = self.out / 'step5_extraction'
        self._write_spectrum(d / f'{name}_spec.fits', np.arange(len(fsum)) + 1.0, fsum, fopt, eopt,
                             np.zeros_like(fsum), 'pixel', hdr)
        self._plot_spectrum(d / f'{name}_spec_pixel', np.arange(len(fsum)) + 1.0, fsum, fopt,
                            'pixel', name)
        return fsum

    # 将光谱、背景和信噪比写入FITS。
    def _write_spectrum(self, path, wave, fsum, fopt, err, sky, unit, hdr, has_snr=()):
        path.parent.mkdir(parents=True, exist_ok=True)
        cols = [fits.Column(name='WAVELENGTH', format='D', array=wave),
                fits.Column(name='FLUX_OPT', format='D', array=fopt),
                fits.Column(name='FLUX_SUM', format='D', array=fsum),
                fits.Column(name='ERROR', format='D', array=np.nan_to_num(err)),
                fits.Column(name='SKY', format='D', array=sky)]
        t = fits.BinTableHDU.from_columns(cols)
        t.header['WAVEUNIT'] = unit
        for e in has_snr:
            key = ('SNRP' if e['unit'] == 'pixel' else 'SNRW') + str(int(e['frac'] * 100))
            t.header[key] = (round(e['snr'], 2), f"SNR at {int(e['frac'] * 100)}% of range")
        t.header['EXPTIME'] = float(hdr.get('EXPTIME', 0) or 0)
        t.header['OBJECT'] = str(hdr.get('OBJECT', path.stem.split('_')[0]))
        fits.HDUList([fits.PrimaryHDU(), t]).writeto(path, overwrite=True)

    # 绘制求和与最优提取光谱。
    def _plot_spectrum(self, path, wave, fsum, fopt, unit, title, snr=()):
        if not self.save_plots:
            return
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        fig, ax = plt.subplots(figsize=(10, 4))
        ax.plot(wave, fsum, lw=0.5, color='0.6', label='sum')
        ax.plot(wave, fopt, lw=0.6, color='b', label='optimal')
        per = 'pixel' if unit == 'pixel' else 'A'
        ymax = np.nanmax(fopt)
        for k, e in enumerate(snr):
            ax.axvline(e['x'], color='r', ls='--', lw=0.6)
            ax.text(e['x'], 0.96 - 0.08 * (k % 2), f"SNR={e['snr']:.1f}@{e['x']:.0f} {per}", transform=ax.get_xaxis_transform(),
                    ha='center', va='top', color='r', fontsize=9,
                    bbox=dict(fc='w', ec='r', lw=0.5, alpha=0.9))
        ax.set_ylim(top=ymax * 1.25)
        ax.set_xlabel('Wavelength (A)' if unit == 'Angstrom' else 'Pixel')
        ax.set_ylabel('Counts')
        ax.set_title(title)
        ax.legend(loc='center left')
        for ext in {self.fig_format, 'pdf'}:
            fig.savefig(Path(str(path) + '.' + ext), dpi=110, bbox_inches='tight')
        plt.close(fig)

    # 绘制天空光谱并标记发射线。
    def _plot_sky(self, path, x, sky, unit, title):
        """Night-sky spectrum (per spatial pixel) with detected emission lines marked."""
        peaks = self._find_peaks(np.nan_to_num(sky))
        idx = np.arange(len(sky))
        px = np.array([np.interp(c, idx, x) for c, _ in peaks])
        pa = np.array([a for _, a in peaks])
        np.savetxt(Path(str(path) + '_lines.csv'), np.column_stack([px, pa]) if len(px) else np.empty((0, 2)),
                   delimiter=',', header=f'position({unit}),amplitude', fmt='%.3f', comments='# ')
        if not self.save_plots:
            return
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(10, 4))
        ax.plot(x, sky, lw=0.5, color='k')
        top = np.argsort(pa)[::-1][:25] if len(pa) else []
        for i in top:
            ax.axvline(px[i], color='r', lw=0.4, alpha=0.6)
            ax.text(px[i], 1.0, f'{px[i]:.1f}', rotation=90, transform=ax.get_xaxis_transform(),
                    ha='center', va='bottom', color='r', fontsize=6)
        ax.set_xlabel('Wavelength (A)' if unit == 'Angstrom' else 'Pixel')
        ax.set_ylabel('Sky counts / spatial pixel')
        ax.set_title(f'{title} night-sky spectrum ({len(peaks)} lines)', pad=28)
        for ext in {self.fig_format, 'pdf'}:
            fig.savefig(Path(str(path) + '.' + ext), dpi=110, bbox_inches='tight')
        plt.close(fig)

    # 绘制标定光谱及波长拟合残差。
    def _plot_arc(self, path, arc, info, coef):
        if not self.save_plots:
            return
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        x, y, keep, rms = info
        fig, (a1, a2) = plt.subplots(2, 1, figsize=(10, 6))
        a1.plot(np.polyval(coef, np.arange(len(arc))), arc, lw=0.5)
        for w in y[keep]:
            a1.axvline(w, color='r', lw=0.3)
        a1.set_xlabel('Wavelength (A)')
        a2.plot(y[keep], (y - np.polyval(coef, x))[keep], 'k.')
        a2.set_ylabel('residual (A)')
        a2.set_title(f'rms={rms:.3f} A')
        fig.savefig(Path(str(path) + '.' + self.fig_format), dpi=110, bbox_inches='tight')
        plt.close(fig)
