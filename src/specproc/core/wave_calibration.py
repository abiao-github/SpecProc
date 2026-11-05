"""
Wavelength calibration stage of spectral reduction pipeline.

Handles ThAr line identification, polynomial fitting,
and wavelength solution computation.
"""
# 匹配参考发射线并拟合二维波长定标模型。

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import cm
import logging
from typing import Tuple, Optional, List, Dict
from scipy.optimize import curve_fit
from pathlib import Path
from specproc.core.data_structures import WaveCalib, Spectrum, SpectraSet
from specproc.core.peak_fitting import fit_peak_center
from specproc.config.config_manager import ConfigManager
from specproc.utils.fits_io import read_fits_image
from specproc.plotting.spectra_plotter import plot_wavelength_calibration

logger = logging.getLogger(__name__)


# 负责建立并应用波长定标模型
class WavelengthCalibrator:
    """Handles wavelength calibration from reference lines."""

    # 初始化对象状态
    def __init__(self):
        """Initialize wavelength calibrator."""
        self.wave_calib = None
        self.reference_lines = {}

    # 读取灯谱参考波长列表
    def load_line_list(self, filepath: str) -> np.ndarray:
        """
        Load reference line list for calibration lamps from a file.
        Assumes wavelength is the first column.
        """
        logger.info(f"Loading line list from {filepath}...")

        try:
            wavelengths = []
            with open(filepath, 'r', encoding='utf-8') as f:
                for line in f:
                    line = line.split('#')[0].strip()
                    if not line:
                        continue
                    parts = line.split()
                    if parts:
                        try:
                            wavelengths.append(float(parts[0]))
                        except ValueError:
                            pass
            wavelengths = np.array(wavelengths)
            wavelengths = wavelengths[~np.isnan(wavelengths)]
            return np.sort(wavelengths)
        except Exception as e:
            logger.error(f"Failed to load linelist from {filepath}: {e}")
            return np.array([])

    # 检测光谱中的发射线峰
    def detect_lines_in_spectrum(self, spectrum: np.ndarray,
                                wavelength: Optional[np.ndarray] = None,
                                threshold: float = 10.0) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Detect emission lines in spectrum.

        Args:
            spectrum: 1D flux array
            wavelength: Optional wavelength array (for initial guess)
            threshold: DetectionThreshold in sigma

        Returns:
            Tuple of (pixel_positions, line_strengths, line_snrs)
        """
        from scipy.signal import find_peaks

        # 修复：使用相邻像素差值的 MAD (中位数绝对偏差) 来稳健估计高斯背景噪声。
        # 发射线光谱中存在极少数超强峰，如果用 np.std(spectrum) 会导致算出的噪声大得离谱，
        # 从而把成百上千条正常强度的发射线全部当成背景过滤掉！
        diffs = np.diff(spectrum)
        # 相邻相减后噪声放大了 sqrt(2) 倍，且 MAD 到 Sigma 的转换系数为 0.6745
        mad_noise = np.nanmedian(np.abs(diffs - np.nanmedian(diffs))) / (np.sqrt(2) * 0.6745)
        
        if mad_noise < 1e-6:
            mad_noise = np.nanstd(spectrum)
            
        # 使用 prominence (相对局部背景的凸起高度) 代替绝对 height
        prom_thresh = max(threshold * mad_noise, 1e-3)
        
        # distance=5 强制要求相邻的 peaks 至少相隔 5 个像素，防止将混叠的“宽肩膀”识别为多根线
        peaks, properties = find_peaks(spectrum, prominence=prom_thresh, distance=5)

        logger.debug(f"Detected {len(peaks)} emission lines")

        grid = np.arange(len(spectrum), dtype=float)
        peaks_int = np.asarray(peaks, dtype=int)
        centres = []
        for k, i in enumerate(peaks_int):
            left_width = right_width = 4
            if k > 0:
                left_width = min(4, max(2, int((i - peaks_int[k - 1]) // 2)))
            if k + 1 < len(peaks_int):
                right_width = min(4, max(2, int((peaks_int[k + 1] - i) // 2)))
            centres.append(fit_peak_center(grid, spectrum, int(i), half_width=4,
                                           left_width=left_width,
                                           right_width=right_width))
        centres = np.asarray(centres, dtype=float)
        snrs = properties['prominences'] / mad_noise if mad_noise > 0 else np.zeros_like(peaks, dtype=float)
        return centres, spectrum[peaks], snrs

    # 读取级次与波长锚点
    def load_anchor_file(self, anchor_file: str) -> Dict[int, np.ndarray]:
        """
        读取灯谱证认的发射线锚点文件。
        预期仅包含两列信息：True_Order(真实级次号), Wavelength(波长)。
        """
        anchors = {}
        try:
            # 使用 genfromtxt 解析 CSV，自动过滤掉带有字符串的 Header 行
            data = np.genfromtxt(anchor_file, delimiter=',')
            data = data[~np.isnan(data).any(axis=1)]
            
            if data.shape[1] >= 2:
                for row in data:
                    m_ref = int(row[0])
                    wave = float(row[1])
                    if m_ref not in anchors:
                        anchors[m_ref] = []
                    anchors[m_ref].append(wave)
                for m in anchors:
                    anchors[m] = np.sort(np.array(anchors[m]))
            else:
                logger.error("Anchor file must contain at least 2 columns: Order and Wavelength.")
                raise ValueError(
                    "提供的锚点文件只有 1 列。请检查是否为 [级次号, 波长] 格式。"
                )
        except Exception as e:
            logger.error(f"Failed to load anchor file {anchor_file}: {e}")
            raise
        return anchors

    # 从文件恢复波长定标模型
    def load_calibration(self, filepath: str) -> WaveCalib:
        """
        从先前保存的 FITS 文件中加载波长定标模型 (Pre-computed Calibration File).
        """
        from astropy.io import fits
        logger.info(f"Loading precomputed wavelength calibration from {filepath}...")
        with fits.open(filepath) as hdul:
            hdr = hdul[0].header
            coefs = hdul[0].data
            wave_calib = WaveCalib(
                poly_coef=coefs,
                xorder=hdr.get('XORDER', 4),
                yorder=hdr.get('YORDER', 4),
                poly_type=hdr.get('POLY_TYP', 'chebyshev'),
                domain_x=(hdr.get('X_MIN', 0.0), hdr.get('X_MAX', 4000.0)),
                domain_y=(hdr.get('Y_MIN', 0.0), hdr.get('Y_MAX', 150.0)),
                line_pixels=np.array([]),
                line_orders=np.array([]),
                line_catalog=np.array([]),
                rms=hdr.get('RMS', 0.0),
                nlines=hdr.get('NLINES', 0),
                calib_type=hdr.get('CALIBTYP', 'Unknown')
            )
            wave_calib.delta_m = hdr.get('DELTAM', 0)
            if 'ORDERFIT' in hdul:
                tab = hdul['ORDERFIT'].data
                wave_calib.order_coefs = {int(m): np.array(c, dtype=float)
                                          for m, c in zip(tab['ORDER'], tab['COEF'])}
        self.wave_calib = wave_calib
        return wave_calib

    # 搜索级次偏移并匹配锚点
    def find_order_offset_and_match(self, detected_peaks: Dict[int, np.ndarray], anchors: Dict[int, np.ndarray],
                                    fixed_offset: Optional[int] = None) -> Tuple[int, List[Tuple[float, int, float]]]:
        """
        寻找最佳偏移 (delta_m = True_Order - Aperture_Number) 并返回成功匹配的锚点 (X_obs, True_Order, wavelength)。
        两阶段盲搜法：
        1. 仅使用少量均匀分布的孔径(Aperture)和强峰，显式区分正反色散方向快速搜索最佳 delta_m。
        2. 利用确定的方向和偏移量，对全量孔径进行秒级匹配。
        """
        all_m_obs = sorted(list(detected_peaks.keys()))
        if not all_m_obs:
            return 0, []

        # ====================================================================
        # Phase 1: 抽取约三分之一均匀分布的孔径进行极速盲搜，确定最佳偏移量
        # ====================================================================
        step = 3 if len(all_m_obs) >= 3 else 1
        search_apertures = all_m_obs[::step]

        best_offset = 0
        best_direction = 1  # 1: Left-to-Right, -1: Right-to-Left
        max_total_matches = -1

        # 动态计算 delta_m 的合理搜索范围
        min_anchor_m = min(anchors.keys())
        
        search_range_start = 0
        search_range_end = min_anchor_m
        if fixed_offset is not None:
            # Blind search is ambiguous when anchors are sparse; trust the configured offset.
            search_range_start = search_range_end = int(fixed_offset)

        logger.info(f"Phase 1: Fast blind matching on {len(search_apertures)} selected apertures. Searching delta_m from {search_range_start} to {search_range_end}...")

        for direction in [1, -1]:
            for delta_m in range(search_range_start, search_range_end + 1):
                total_matches = 0
                for m_obs in search_apertures:
                    X = detected_peaks[m_obs]
                    
                    m_ref = m_obs + delta_m
                    if m_ref not in anchors: continue
                    
                    W = anchors[m_ref]
                    if len(W) < 2 or len(X) < 2: continue

                    best_inliers = []

                    # 遍历所有可能的锚点波长对和观测峰对
                    for k in range(len(W) - 1):
                        for l in range(k + 1, len(W)):
                            dw = W[l] - W[k]
                            for i in range(len(X) - 1):
                                for j in range(i + 1, len(X)):
                                    dx = X[j] - X[i]
                                    if dx == 0: continue

                                    # 色散率 a = direction * |dw/dx|
                                    a_test = direction * abs(dw / dx)
                                    if not (0.002 < abs(a_test) < 0.5): continue
                                    
                                    # 正向色散：波长随 X 增加。所以较小的 X[i] 对应较小的 W[k]
                                    # 反向色散：波长随 X 减小。所以较大的 X[j] 对应较小的 W[k]
                                    b = W[k] - a_test * (X[i] if direction == 1 else X[j])

                                    W_pred = a_test * X + b
                                    diffs = np.abs(W[:, None] - W_pred[None, :])
                                    min_dist_idx = np.argmin(diffs, axis=0)
                                    min_dist = np.min(diffs, axis=0)

                                    valid = min_dist < 1.5
                                    if np.sum(valid) > len(best_inliers):
                                        matched_w = min_dist_idx[valid]
                                        _, unique_indices = np.unique(matched_w, return_index=True)
                                        if len(unique_indices) > len(best_inliers):
                                            valid_x_idx = np.where(valid)[0][unique_indices]
                                            valid_w_idx = matched_w[unique_indices]
                                            best_inliers = [(X[vx], m_ref, W[vw]) for vx, vw in zip(valid_x_idx, valid_w_idx)]
                                            
                    if len(best_inliers) >= max(2, min(3, len(W))):
                        total_matches += len(best_inliers)

                if total_matches > max_total_matches:
                    max_total_matches = total_matches
                    best_offset = delta_m
                    best_direction = direction

        if max_total_matches <= 0:
            logger.warning("盲配模式匹配失败！未找到合适的级次(Order)偏移。请检查峰值提取、锚点文件或放宽色散率限制。")
            return 0, []
            
        dir_str = "Left-to-Right (+)" if best_direction == 1 else "Right-to-Left (-)"
        logger.info(f"Phase 1 complete: Found delta_m={best_offset}, dispersion direction={dir_str}")

        # ====================================================================
        # Phase 2: 使用确定的最佳偏移量和色散方向，全量孔径(Apertures)提取匹配点
        # ====================================================================
        logger.info(f"Phase 2: Extracting matched anchors across all apertures using delta_m={best_offset}...")
        best_global_matched_points = []
        global_total_matches = 0

        for m_obs in all_m_obs:
            X = detected_peaks[m_obs]
            m_ref = m_obs + best_offset
            if m_ref not in anchors: continue
            W = anchors[m_ref]
            if len(W) < 2 or len(X) < 2: continue

            best_inliers = []
            for k in range(len(W) - 1):
                for l in range(k + 1, len(W)):
                    dw = W[l] - W[k]
                    for i in range(len(X) - 1):
                        for j in range(i + 1, len(X)):
                            dx = X[j] - X[i]
                            if dx == 0: continue
                            
                            # 直接使用 Phase 1 确定的最佳色散方向
                            a_test = best_direction * abs(dw / dx)
                            if not (0.002 < abs(a_test) < 0.5): continue
                            
                            b = W[k] - a_test * (X[i] if best_direction == 1 else X[j])
                                
                            W_pred = a_test * X + b
                            diffs = np.abs(W[:, None] - W_pred[None, :])
                            min_dist_idx = np.argmin(diffs, axis=0)
                            min_dist = np.min(diffs, axis=0)
                            
                            valid = min_dist < 2.5  # 放宽初筛容差以包容边缘非线性
                            if np.sum(valid) > len(best_inliers):
                                matched_w = min_dist_idx[valid]
                                _, unique_indices = np.unique(matched_w, return_index=True)
                                
                                # 二次多项式精校准 (针对非线性强烈的蓝端级次)
                                if len(unique_indices) >= 3:
                                    vx = X[np.where(valid)[0][unique_indices]]
                                    vw = W[matched_w[unique_indices]]
                                    try:
                                        p2 = np.polyfit(vx, vw, 2)
                                        W_pred2 = np.polyval(p2, X)
                                        diffs2 = np.abs(W[:, None] - W_pred2[None, :])
                                        min_dist_idx2 = np.argmin(diffs2, axis=0)
                                        min_dist2 = np.min(diffs2, axis=0)
                                        valid2 = min_dist2 < 1.5  # 二次拟合后收紧容差
                                        matched_w2 = min_dist_idx2[valid2]
                                        _, uniq2 = np.unique(matched_w2, return_index=True)
                                        if len(uniq2) >= len(unique_indices):
                                            valid = valid2
                                            matched_w = matched_w2
                                            unique_indices = uniq2
                                    except Exception:
                                        pass

                                if len(unique_indices) > len(best_inliers):
                                    valid_x_idx = np.where(valid)[0][unique_indices]
                                    valid_w_idx = matched_w[unique_indices]
                                    best_inliers = [(X[vx], m_ref, W[vw]) for vx, vw in zip(valid_x_idx, valid_w_idx)]
                                    
            if len(best_inliers) >= max(2, min(3, len(W))):
                global_total_matches += len(best_inliers)
                best_global_matched_points.extend(best_inliers)

        logger.info(f"Blind match fully complete: delta_m={best_offset} with {global_total_matches} lines matched across all apertures.")
            
        return best_offset, best_global_matched_points

    # 利用锚点表中的 X 坐标匹配：先找全局 X 平移和级次偏移，再逐级次微调
    def match_anchor_positions(self, scoring_peaks: Dict[int, np.ndarray], all_peaks: Dict[int, np.ndarray],
                               ref_table: np.ndarray, fixed_offset: Optional[int] = None,
                               max_shift: float = 150.0, tol: float = 3.0,
                               refine_shift: float = 6.0) -> Tuple[int, List[Tuple[float, int, float]]]:
        """Match an identified-line table (order, X, wavelength) to detected peaks by position.

        The lamp frame is shifted by a small rigid amount relative to the table, so the shift
        (and the order offset delta_m = True_Order - Aperture) that makes the strongest peaks
        coincide with the table X positions is searched first; each aperture is then refined
        within a few pixels and every table line takes its nearest unused peak within ``tol``.
        """
        table = {int(m): (ref_table[ref_table[:, 0].astype(int) == m, 1],
                          ref_table[ref_table[:, 0].astype(int) == m, 2])
                 for m in np.unique(ref_table[:, 0].astype(int))}

        def count(peaks, xref, shift):
            if len(peaks) == 0 or len(xref) == 0:
                return 0
            d = np.abs(peaks[None, :] - (xref[:, None] + shift))
            return int(np.sum(d.min(axis=1) < tol))

        apertures = sorted(scoring_peaks)
        deltas = [int(fixed_offset)] if fixed_offset is not None else range(0, max(table) + 1)
        shifts = np.arange(-max_shift, max_shift + 0.5, 1.0)
        best = (0, 0, 0.0)
        for dm in deltas:
            for sh in shifts:
                score = sum(count(np.asarray(scoring_peaks[a]), table[a + dm][0], sh)
                            for a in apertures if a + dm in table)
                if score > best[0]:
                    best = (score, int(dm), float(sh))
        score, delta_m, shift0 = best
        if score <= 0:
            logger.warning("Anchor position matching failed: no table line coincides with a detected peak.")
            return 0, []
        logger.info(f"Anchor position matching: delta_m={delta_m}, global X shift {shift0:+.0f} px "
                    f"({score} strong peaks coincide with table lines).")

        matched: List[Tuple[float, int, float]] = []
        for a in sorted(all_peaks):
            m_ref = a + delta_m
            if m_ref not in table:
                continue
            xref, wref = table[m_ref]
            peaks = np.asarray(all_peaks[a], dtype=float)
            if len(peaks) == 0:
                continue
            offsets = shift0 + np.arange(-refine_shift, refine_shift + 0.5, 0.5)
            sh = max(offsets, key=lambda o: (count(peaks, xref, o), -abs(o - shift0)))
            pairs = sorted((abs(peaks[j] - (xref[i] + sh)), i, j)
                           for i in range(len(xref)) for j in range(len(peaks)))
            used_i, used_j = set(), set()
            for dist, i, j in pairs:
                if dist >= tol:
                    break
                if i in used_i or j in used_j:
                    continue
                used_i.add(i)
                used_j.add(j)
                matched.append((float(peaks[j]), m_ref, float(wref[i])))
        logger.info(f"Anchor position matching complete: {len(matched)} lines matched across all apertures.")
        return delta_m, matched

    # 用已定标的一维参考灯谱确定级次偏移并生成粗定标锚点
    def match_reference_spectrum(self, spectra: Dict[int, object], reference_file: str,
                                 min_corr: float = 0.5, max_lag: int = 80,
                                 sample_step: int = 16) -> Tuple[int, List[Tuple[float, int, float]]]:
        """Match extracted lamp orders to a calibrated, flux-normalised reference lamp spectrum.

        The reference file has columns ``order, wavelength, flux`` (one block per
        true order, pixel-ordered). Each extracted order is cross-correlated with every
        reference order; the true-order offset ``delta_m`` maximising the summed
        correlation is selected and (pixel, order, wavelength) anchors are sampled from
        the reference with the per-order pixel lag applied.
        """
        from scipy.ndimage import median_filter

        ref = np.loadtxt(reference_file)
        ref_orders = sorted(int(o) for o in np.unique(ref[:, 0]))
        ref_blocks = {o: ref[ref[:, 0] == o] for o in ref_orders}

        def normalise(y):
            y = np.nan_to_num(np.asarray(y, dtype=float))
            y = y - median_filter(y, 201)
            return y / (y.std() + 1e-12)

        corr = {}
        for m_obs, spec in spectra.items():
            a = normalise(spec.flux if hasattr(spec, 'flux') else spec)
            fa = np.fft.rfft(a)
            for o in ref_orders:
                b = normalise(np.pad(ref_blocks[o][:, 2], (0, max(0, len(a) - len(ref_blocks[o])))))[:len(a)]
                c = np.fft.irfft(fa * np.conj(np.fft.rfft(b)), n=len(a)) / len(a)
                window = np.concatenate([c[-max_lag:], c[:max_lag + 1]])
                k = int(np.argmax(window))
                corr[(m_obs, o)] = (float(window[k]), k - max_lag)

        m_list = sorted(spectra.keys())
        best_delta, best_score = 0, -1.0
        for delta in range(ref_orders[0] - m_list[-1], ref_orders[-1] - m_list[0] + 1):
            score = sum(max(corr[(m, m + delta)][0], 0.0) for m in m_list if (m, m + delta) in corr)
            if score > best_score:
                best_delta, best_score = delta, score

        from scipy.signal import find_peaks
        anchors = []
        self.last_peak_anchors = []
        used = 0
        for m in m_list:
            key = (m, m + best_delta)
            if key not in corr or corr[key][0] < min_corr:
                continue
            lag = corr[key][1]
            block = ref_blocks[m + best_delta]
            n = len(spectra[m].flux if hasattr(spectra[m], 'flux') else spectra[m])
            for j in range(0, len(block), sample_step):
                x = j + lag
                if 0 <= x < n:
                    anchors.append((float(x), m + best_delta, float(block[j, 1])))
            # 参考谱中的显著峰吸附到观测谱的局部极大值, 用于诊断图标注
            obs = np.nan_to_num(np.asarray(spectra[m].flux if hasattr(spectra[m], 'flux') else spectra[m], float))
            d = np.diff(obs)
            noise = np.median(np.abs(d - np.median(d))) / (np.sqrt(2) * 0.6745) + 1e-12
            ref_flux = np.nan_to_num(block[:, 2])
            pk, props = find_peaks(ref_flux, prominence=max(5 * ref_flux.std() * 0.2, 1e-12))
            for j in pk[np.argsort(props['prominences'])[::-1][:15]]:
                c = int(round(j + lag))
                lo, hi = max(c - 4, 1), min(c + 5, n - 1)
                if hi - lo < 3:
                    continue
                k = lo + int(np.argmax(obs[lo:hi]))
                if obs[k] >= obs[k - 1] and obs[k] >= obs[k + 1] and obs[k] - np.median(obs) > 8 * noise:
                    self.last_peak_anchors.append((float(k), m + best_delta, float(block[j, 1])))
            used += 1
        self.last_match_orders = used
        self.last_match_score = best_score
        logger.info(f"Reference-spectrum match: delta_m={best_delta}, {used}/{len(m_list)} orders "
                    f"correlate (r >= {min_corr}), {len(anchors)} anchor points.")
        return best_delta, anchors

    # 用粗定标匹配完整谱线目录
    def match_full_catalog(self, detected_peaks: Dict[int, np.ndarray], rough_calib: WaveCalib, delta_m: int, full_linelist: np.ndarray, tolerance: float = 0.5,
                           tolerance_fitted: Optional[float] = None) -> Tuple[np.ndarray, np.ndarray]:
        """
        利用拟合的 2D 粗校准模型，预测所有 detected_peaks 的波长并匹配全量线库。
        若模型含独立级次拟合，这些级次可用更严格的 tolerance_fitted。
        """
        pix_pos = []
        matched_wave = []

        if len(full_linelist) == 0:
            logger.error("全量线库为空，跳过全量匹配阶段。")
            return np.array([]), np.array([])

        for m_obs, peaks_x in detected_peaks.items():
            if len(peaks_x) == 0:
                continue
            m_true = m_obs + delta_m
            
            wave_pred = rough_calib.apply_to_pixel(peaks_x, np.full_like(peaks_x, m_true))
            
            tol = tolerance
            if tolerance_fitted is not None and rough_calib.order_coefs and m_true in rough_calib.order_coefs:
                tol = tolerance_fitted

            # 全库匹配寻找最近线
            for x_val, w_pred in zip(peaks_x, wave_pred):
                dist = np.abs(full_linelist - w_pred)
                idx = np.argmin(dist)
                if dist[idx] < tol:
                    pix_pos.append([x_val, m_true])
                    matched_wave.append(full_linelist[idx])
                    
        return np.array(pix_pos), np.array(matched_wave)

    # 拟合像素与级次的二维波长关系
    def fit_wavelength_polynomial(self, pixel_positions: np.ndarray,
                                 wavelengths: np.ndarray,
                                 xorder: int = 4, yorder: int = 4,
                                 poly_type: str = 'chebyshev',
                                 max_residual: float = 0.5) -> WaveCalib:
        """
        Fit 2D wavelength polynomial.

        Args:
            pixel_positions: Array of (x, y) pixel coordinates
            wavelengths: Corresponding wavelengths
            xorder: Polynomial order in X
            yorder: Polynomial order in Y
            poly_type: Polynomial basis ('chebyshev', 'legendre' or 'polynomial')
            max_residual: Absolute maximum residual (in Angstroms) to keep a line

        Returns:
            WaveCalib object with polynomial solution
        """
        logger.info(f"Fitting wavelength 2D {poly_type} polynomial (orders: {xorder}, {yorder})...")

        if len(pixel_positions) != len(wavelengths):
            raise ValueError("Pixel positions and wavelengths must have same length")

        # Extract x, y coordinates
        if pixel_positions.ndim == 2:
            x_pix = pixel_positions[:, 0]
            y_pix = pixel_positions[:, 1]
        else:
            x_pix = pixel_positions
            y_pix = np.zeros_like(x_pix)

        # 改为拟合 m * lambda
        m_lambda = wavelengths * y_pix

        # 获取参与定标数据的边界值，用于归一化
        domain_x = (float(np.min(x_pix)), float(np.max(x_pix)))
        domain_y = (float(np.min(y_pix)), float(np.max(y_pix)))

        min_required = (xorder + 1) * (yorder + 1)
        
        # Auto-downgrade polynomial order if not enough points
        orig_xorder, orig_yorder = xorder, yorder
        while len(m_lambda) < min_required and (xorder > 1 or yorder > 1):
            if xorder > yorder:
                xorder -= 1
            elif yorder > xorder:
                yorder -= 1
            else:
                xorder -= 1
                yorder -= 1
            min_required = (xorder + 1) * (yorder + 1)
            
        if len(m_lambda) < min_required:
            raise RuntimeError(f"Not enough matched lines ({len(m_lambda)}) to fit even a 1x1 2D polynomial. "
                               f"Requires at least {min_required} lines. Please check your data and anchor files.")
                               
        if xorder != orig_xorder or yorder != orig_yorder:
            logger.warning(f"Auto-downgraded polynomial order from {orig_xorder}x{orig_yorder} "
                           f"to {xorder}x{yorder} due to limited points ({len(m_lambda)}).")

        # 将像素和级次归一化到 [-1, 1] 区间
        x_norm = 2.0 * (x_pix - domain_x[0]) / max(domain_x[1] - domain_x[0], 1e-6) - 1.0
        y_norm = 2.0 * (y_pix - domain_y[0]) / max(domain_y[1] - domain_y[0], 1e-6) - 1.0

        # 利用 numpy.polynomial 的 Vander 函数快速构建 2D 正交基设计矩阵
        if poly_type == 'chebyshev':
            from numpy.polynomial.chebyshev import chebvander2d
            A = chebvander2d(x_norm, y_norm, [xorder, yorder])
        elif poly_type == 'legendre':
            from numpy.polynomial.legendre import legvander2d
            A = legvander2d(x_norm, y_norm, [xorder, yorder])
        else:
            from numpy.polynomial.polynomial import polyvander2d
            A = polyvander2d(x_norm, y_norm, [xorder, yorder])

        # Solve least squares
        try:
            keep = np.ones(len(m_lambda), dtype=bool)
            coeffs_1d = None
            for _ in range(5):
                if np.sum(keep) < (xorder + 1) * (yorder + 1):
                    logger.warning("Too few points left for 2D poly fit after clipping.")
                    break
                    
                A_fit = A[keep]
                m_lambda_fit = m_lambda[keep]
                
                coeffs_1d, residuals_lstsq, rank, s = np.linalg.lstsq(A_fit, m_lambda_fit, rcond=None)
                
                pred_m_lambda = A @ coeffs_1d
                pred_lambda = pred_m_lambda / y_pix
                residuals = pred_lambda - wavelengths
                
                rms = np.sqrt(np.mean(residuals[keep] ** 2))
                
                # 修复：设置 3-sigma 剔除的下限 (例如 0.1 埃)。
                # 防止刚性模型（如低阶多项式）在早期迭代时因为核心点拟合得太好导致 rms 极小，
                # 进而产生“滚雪球”式的过度剔除 (Over-clipping)。
                clip_threshold = max(3.0 * rms, 0.1)
                
                # Clip by both robust 3-sigma and absolute maximum residual limit
                new_keep = (np.abs(residuals) < clip_threshold) & (np.abs(residuals) <= max_residual)
                if np.array_equal(new_keep, keep):
                    break
                keep = new_keep
                
            n_rejected = len(keep) - np.sum(keep)
            if n_rejected > 0:
                logger.info(f"Sigma clipping rejected {n_rejected} outliers.")

            if coeffs_1d is None:
                raise RuntimeError(f"Failed to fit 2D polynomial (xorder={xorder}, yorder={yorder}). "
                                   f"Too many outliers rejected. Not enough valid points ({np.sum(keep)}) remaining.")

            pred_m_lambda = A @ coeffs_1d
            pred_lambda = pred_m_lambda / y_pix
            rms = np.sqrt(np.mean((pred_lambda[keep] - wavelengths[keep]) ** 2))

            logger.info(f"Wavelength fit RMS: {rms:.4f} Angstrom (on {np.sum(keep)} lines)")

            # Reshape coefficients into 2D array
            poly_coef = coeffs_1d.reshape((xorder + 1, yorder + 1))

            # Create WaveCalib object (only storing valid inliers in standard fields)
            wave_calib = WaveCalib(
                poly_coef=poly_coef,
                xorder=xorder,
                yorder=yorder,
                poly_type=poly_type,
                domain_x=domain_x,
                domain_y=domain_y,
                line_pixels=x_pix[keep],
                line_orders=y_pix[keep],
                line_catalog=wavelengths[keep],
                rms=rms,
                nlines=int(np.sum(keep)),
                calib_type='ThAr'
            )
            
            # Attach full point arrays for diagnostic plotting
            wave_calib.all_line_pixels = x_pix
            wave_calib.all_line_orders = y_pix
            wave_calib.all_line_catalog = wavelengths
            wave_calib.keep_mask = keep

            self.wave_calib = wave_calib
            return wave_calib

        except Exception as e:
            logger.error(f"Error fitting wavelength polynomial: {e}")
            raise

    # 在全局二维拟合基础上，对每个级次独立拟合一维多项式 lambda(x)
    def fit_per_order(self, wave_calib: WaveCalib, degree: int = 4,
                      min_lines: int = 8, min_span: float = 0.6,
                      max_residual: float = 0.5, max_order_rms: float = 0.0) -> WaveCalib:
        """
        Refit every well-constrained order independently with a 1D Chebyshev polynomial.

        The global 2D fit supplies the initial outlier rejection and remains the fallback
        for orders with too few lines or too little coverage of the detector.

        Args:
            wave_calib: Result of ``fit_wavelength_polynomial`` (inliers in ``keep_mask``).
            degree: Maximum polynomial degree per order (reduced when lines are scarce).
            min_lines: Minimum number of inlier lines required for an order to get its own fit.
            min_span: Minimum fraction of the global x-domain the order's lines must cover.
            max_residual: Absolute residual limit (Angstrom) for lines kept in an order fit.
            max_order_rms: If > 0, an order whose own fit RMS exceeds this falls back to the 2D model.
        """
        from numpy.polynomial.chebyshev import chebfit, chebval

        x_all = wave_calib.all_line_pixels
        m_all = wave_calib.all_line_orders
        w_all = wave_calib.all_line_catalog
        keep_global = wave_calib.keep_mask
        x0, x1 = wave_calib.domain_x
        order_coefs: Dict[int, np.ndarray] = {}
        keep = keep_global.copy()
        order_rms: Dict[int, float] = {}
        # A rigid global model may have rejected genuine lines of orders that deviate from it,
        # so each order starts from every line within a generous window of the global fit.
        res_global = wave_calib.apply_to_pixel(x_all, m_all) - w_all
        candidate = np.abs(res_global) < 3.0 * max_residual

        for m in np.unique(np.rint(m_all)).astype(int):
            sel = (np.rint(m_all) == m) & candidate
            n = int(sel.sum())
            if n < min_lines:
                continue
            xs, ws = x_all[sel], w_all[sel]
            if (xs.max() - xs.min()) < min_span * (x1 - x0):
                continue
            xn = 2.0 * (xs - x0) / max(x1 - x0, 1e-6) - 1.0
            deg = int(max(1, min(degree, (n - 2) // 3)))
            inl = np.ones(n, dtype=bool)
            coef = None
            for _ in range(5):
                if inl.sum() < deg + 2:
                    coef = None
                    break
                coef = chebfit(xn[inl], ws[inl], deg)
                res = chebval(xn, coef) - ws
                rms = np.sqrt(np.mean(res[inl] ** 2))
                new_inl = np.abs(res) < min(max(3.0 * rms, 0.05), max_residual)
                if np.array_equal(new_inl, inl):
                    break
                inl = new_inl
            if coef is None:
                continue
            if max_order_rms > 0 and np.sqrt(np.mean((chebval(xn, coef) - ws)[inl] ** 2)) > max_order_rms:
                logger.info(f"Order {m}: independent fit RMS above {max_order_rms:.2f} A, using the 2D model.")
                continue
            order_coefs[int(m)] = coef
            idx = np.where(sel)[0]
            keep[np.rint(m_all) == m] = False
            keep[idx[inl]] = True
            order_rms[int(m)] = float(np.sqrt(np.mean((chebval(xn, coef) - ws)[inl] ** 2)))

        wave_calib.order_coefs = order_coefs
        wave_calib.keep_mask = keep
        wave_calib.line_pixels = x_all[keep]
        wave_calib.line_orders = m_all[keep]
        wave_calib.line_catalog = w_all[keep]
        wave_calib.nlines = int(keep.sum())
        res_all = wave_calib.apply_to_pixel(x_all, m_all) - w_all
        wave_calib.rms = float(np.sqrt(np.mean(res_all[keep] ** 2)))
        rms_txt = ', '.join(f"{m}:{r:.3f}" for m, r in sorted(order_rms.items()))
        logger.info(f"Per-order fit: {len(order_coefs)} orders fitted independently "
                    f"(degree<= {degree}); other orders use the 2D model. "
                    f"Final RMS {wave_calib.rms:.4f} A on {wave_calib.nlines} lines. "
                    f"Per-order RMS [A]: {rms_txt}")
        return wave_calib

    # 将波长定标应用于一维光谱
    def apply_wavelength_calibration(self, spectrum: np.ndarray,
                                    pixel_array: np.ndarray,
                                    aperture_y: float = 0.0) -> np.ndarray:
        """
        Apply wavelength calibration to spectrum pixels.

        Args:
            spectrum: 1D spectrum array
            pixel_array: Pixel coordinate array
            aperture_y: Y coordinate (True Order number) of aperture

        Returns:
            Wavelength array
        """
        if self.wave_calib is None:
            raise RuntimeError("No wavelength calibration available")

        # Apply polynomial
        wavelength = self.wave_calib.apply_to_pixel(pixel_array, aperture_y * np.ones_like(pixel_array))

        return wavelength

    # 保存波长定标系数和元数据
    def save_calibration(self, output_path: str):
        """Save wavelength calibration to FITS file."""
        if self.wave_calib is None:
            raise RuntimeError("No wavelength calibration to save")

        from astropy.io import fits

        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Build a proper fits.Header (PrimaryHDU rejects plain dicts)
        hdr = fits.Header()
        hdr['XORDER']    = self.wave_calib.xorder
        hdr['YORDER']    = self.wave_calib.yorder
        hdr['RMS']       = self.wave_calib.rms
        hdr['NLINES']    = self.wave_calib.nlines
        hdr['CALIBTYP']  = self.wave_calib.calib_type
        hdr['DELTAM']    = getattr(self.wave_calib, 'delta_m', 0)
        hdr['POLY_TYP']  = getattr(self.wave_calib, 'poly_type', 'chebyshev')
        if hasattr(self.wave_calib, 'domain_x'):
            hdr['X_MIN'], hdr['X_MAX'] = self.wave_calib.domain_x
            hdr['Y_MIN'], hdr['Y_MAX'] = self.wave_calib.domain_y

        hdul = fits.HDUList([
            fits.PrimaryHDU(data=self.wave_calib.poly_coef, header=hdr),
        ])
        order_coefs = getattr(self.wave_calib, 'order_coefs', None)
        if order_coefs:
            width = max(len(c) for c in order_coefs.values())
            orders = sorted(order_coefs)
            table = np.zeros((len(orders), width))
            for i, m in enumerate(orders):
                table[i, :len(order_coefs[m])] = order_coefs[m]
            hdul.append(fits.BinTableHDU.from_columns([
                fits.Column(name='ORDER', format='J', array=np.array(orders)),
                fits.Column(name='COEF', format=f'{width}D', array=table),
            ], name='ORDERFIT'))

        hdul.writeto(str(output_path), overwrite=True)
        logger.info(f"Saved wavelength calibration to {output_path}")


# 短波（蓝）到长波（深红）的渐变色
def _wavelength_cmap():
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list(
        'blue_to_darkred', ['#1f3bff', '#00a8ff', '#00b050', '#f0c800', '#ff5a00', '#8b0000'])


# 绘制定标散点与残差诊断图
def _plot_calib_diagnostic(wave_calib: WaveCalib, plot_file: str, max_residual: float = 0.5):
    """生成波长定标诊断图（散点图及拟合残差），颜色随波长由蓝渐变到深红。"""
    from matplotlib.collections import LineCollection
    from matplotlib.colors import Normalize
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True, gridspec_kw={'height_ratios': [3, 1]})

    x_pix = getattr(wave_calib, 'all_line_pixels', wave_calib.line_pixels)
    m_orders = getattr(wave_calib, 'all_line_orders', wave_calib.line_orders)
    wave_ref = getattr(wave_calib, 'all_line_catalog', wave_calib.line_catalog)
    keep_mask = getattr(wave_calib, 'keep_mask', np.ones(len(x_pix), dtype=bool))

    n_total = len(x_pix)
    n_valid = wave_calib.nlines

    wave_pred = wave_calib.apply_to_pixel(x_pix, m_orders)
    residuals = wave_pred - wave_ref

    cmap = _wavelength_cmap()
    norm = Normalize(float(np.nanmin(wave_ref)), float(np.nanmax(wave_ref)))
    sc = None
    for m in np.unique(m_orders):
        mask = (m_orders == m)
        inlier_mask = mask & keep_mask
        if np.any(inlier_mask):
            sc = ax1.scatter(x_pix[inlier_mask], wave_ref[inlier_mask], c=wave_ref[inlier_mask],
                             cmap=cmap, norm=norm, s=15, alpha=0.9)
            ax2.scatter(x_pix[inlier_mask], residuals[inlier_mask], c=wave_ref[inlier_mask],
                        cmap=cmap, norm=norm, s=15, alpha=0.9)

        outlier_mask = mask & ~keep_mask
        if np.any(outlier_mask):
            ax1.scatter(x_pix[outlier_mask], wave_ref[outlier_mask], color='gray', marker='x', s=15, alpha=0.6)
            ax2.scatter(x_pix[outlier_mask], residuals[outlier_mask], color='gray', marker='x', s=15, alpha=0.6)

        x_smooth = np.linspace(x_pix[mask].min(), x_pix[mask].max(), 100)
        wave_smooth = wave_calib.apply_to_pixel(x_smooth, np.full_like(x_smooth, m))
        pts = np.column_stack([x_smooth, wave_smooth]).reshape(-1, 1, 2)
        lc = LineCollection(np.concatenate([pts[:-1], pts[1:]], axis=1), cmap=cmap, norm=norm, alpha=0.6)
        lc.set_array(0.5 * (wave_smooth[:-1] + wave_smooth[1:]))
        ax1.add_collection(lc)
    ax1.autoscale_view()

    ax1.set_ylabel(r'Wavelength ($\AA$)')
    ax1.set_title(fr'Wavelength Calibration (RMS = {wave_calib.rms:.4f} $\AA$, N = {n_valid}/{n_total})')
    ax1.grid(True, alpha=0.3)

    ax2.axhline(0, color='k', linestyle='--', alpha=0.5)
    ax2.axhline(max_residual, color='r', linestyle=':', alpha=0.6)
    ax2.axhline(-max_residual, color='r', linestyle=':', alpha=0.6)
    ax2.set_xlabel('Pixel (X)')
    ax2.set_ylabel(r'Residual ($\AA$)')

    # 动态对称缩放 Y 轴，包含异常值
    max_data_res = np.max(np.abs(residuals)) if len(residuals) > 0 else 0.1
    y_lim = max(max_data_res * 1.1, max_residual * 1.5)
    ax2.set_ylim(-y_lim, y_lim)
    ax2.grid(True, alpha=0.3)
    fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=[ax1, ax2], label=r'Wavelength ($\AA$)', pad=0.02)

    plt.savefig(plot_file, dpi=150, bbox_inches='tight')
    plt.close(fig)

# 绘制二维定标曲面与残差图
def _plot_calib_surface_diagnostic(wave_calib: WaveCalib, plot_file: str, max_residual: float = 0.5):
    """生成二维波长定标拟合曲面和残差诊断图。"""
    from matplotlib.colors import Normalize
    fig = plt.figure(figsize=(14, 6))
    
    ax1 = fig.add_subplot(121, projection='3d')
    x_pix = getattr(wave_calib, 'all_line_pixels', wave_calib.line_pixels)
    m_orders = getattr(wave_calib, 'all_line_orders', wave_calib.line_orders)
    wave_ref = getattr(wave_calib, 'all_line_catalog', wave_calib.line_catalog)
    keep_mask = getattr(wave_calib, 'keep_mask', np.ones(len(x_pix), dtype=bool))
    
    n_total = len(x_pix)
    n_valid = wave_calib.nlines
    
    cmap = _wavelength_cmap()
    norm = Normalize(float(np.nanmin(wave_ref)), float(np.nanmax(wave_ref)))

    # 生成用于绘制平滑曲面的网格
    x_grid = np.linspace(x_pix.min(), x_pix.max(), 50)
    m_grid = np.linspace(m_orders.min(), m_orders.max(), 50)
    X, M = np.meshgrid(x_grid, m_grid)
    WAVE = wave_calib.apply_to_pixel(X, M)
    
    # 绘制预测曲面和实际匹配散点
    ax1.plot_surface(X, M, WAVE, cmap=cmap, norm=norm, alpha=0.6, edgecolor='none')
    
    ax1.scatter(x_pix[keep_mask], m_orders[keep_mask], wave_ref[keep_mask], color='k', s=10, alpha=0.8, label='Matched Lines')
    if np.any(~keep_mask):
        ax1.scatter(x_pix[~keep_mask], m_orders[~keep_mask], wave_ref[~keep_mask], color='gray', marker='x', s=10, alpha=0.6, label='Rejected Lines')
        
    ax1.set_xlabel('Pixel (X)')
    ax1.set_ylabel('Order (m)')
    ax1.set_zlabel(r'Wavelength ($\AA$)')
    ax1.set_title(f'2D Wavelength Solution Surface\n({wave_calib.poly_type.capitalize()} Poly, {wave_calib.xorder}x{wave_calib.yorder})')
    
    # 2. 拟合残差图 (二维散点，颜色表示波长)
    ax2 = fig.add_subplot(122)
    residuals = wave_calib.apply_to_pixel(x_pix, m_orders) - wave_ref
    
    scatter = ax2.scatter(x_pix[keep_mask], residuals[keep_mask], c=wave_ref[keep_mask], cmap=cmap, norm=norm, s=15, alpha=0.9)
    if np.any(~keep_mask):
        ax2.scatter(x_pix[~keep_mask], residuals[~keep_mask], color='gray', marker='x', s=15, alpha=0.6)
        
    ax2.axhline(0, color='k', linestyle='--', alpha=0.5)
    ax2.axhline(max_residual, color='r', linestyle=':', alpha=0.8, label=rf'Limit ($\pm${max_residual:.2f})')
    ax2.axhline(-max_residual, color='r', linestyle=':', alpha=0.8)
    ax2.legend(loc='upper right', fontsize='small')
    ax2.set_xlabel('Pixel (X)')
    ax2.set_ylabel(r'Residual ($\AA$)')
    ax2.set_title(fr'Fitting Residuals (RMS = {wave_calib.rms:.4f} $\AA$, N = {n_valid}/{n_total})')
    
    # 动态对称缩放 Y 轴，包含异常值
    max_data_res = np.max(np.abs(residuals)) if len(residuals) > 0 else 0.1
    y_lim = max(max_data_res * 1.1, max_residual * 1.5)
    ax2.set_ylim(-y_lim, y_lim)
    ax2.grid(True, alpha=0.3)
    fig.colorbar(scatter, ax=ax2, label=r'Wavelength ($\AA$)')
    
    plt.tight_layout()
    plt.savefig(plot_file, dpi=150)
    plt.close(fig)

# 绘制各级次锚点匹配并保存为 PDF
def _rematch_anchor_table(detected_peaks: Dict[int, np.ndarray], wave_calib: WaveCalib, delta_m: int,
                          anchor_table: Dict[int, np.ndarray], tol: float = 0.5) -> List[Tuple[float, int, float]]:
    """用最终定标把锚点线表匹配到实测 peaks (每个 peak 至多一条锚点线, 取最近者)。"""
    out: List[Tuple[float, int, float]] = []
    for m_obs, peaks in detected_peaks.items():
        m_true = m_obs + delta_m
        waves = anchor_table.get(m_true)
        if waves is None or len(peaks) == 0:
            continue
        pk = np.asarray(peaks, dtype=float)
        pk_wave = wave_calib.apply_to_pixel(pk, np.full_like(pk, m_true))
        for w in waves:
            j = int(np.argmin(np.abs(pk_wave - w)))
            if abs(pk_wave[j] - w) <= tol:
                out.append((float(pk[j]), m_true, float(w)))
    return out


def _plot_matched_anchors_pdf(lamp_spectra, detected_peaks: Dict[int, np.ndarray], 
                              matched_anchors: List[Tuple[float, int, float]], 
                              delta_m: int, 
                              pdf_path: str,
                              snr_used: Optional[Dict[int, float]] = None,
                              max_labels: int = 5,
                              max_peaks: int = 10,
                              wave_calib: Optional[WaveCalib] = None):
    """将每个级次的光谱及匹配到的锚点画入一个多页 PDF 文件（每级次只标注最强的 max_labels 个）。"""
    with PdfPages(pdf_path) as pdf:
        anchor_dict = {}
        for x_obs, m_true, wave in matched_anchors:
            if m_true not in anchor_dict:
                anchor_dict[m_true] = []
            anchor_dict[m_true].append((x_obs, wave))
            
        spectra_dict = getattr(lamp_spectra, 'spectra', lamp_spectra) if not isinstance(lamp_spectra, dict) else lamp_spectra
        
        for m_obs in sorted(spectra_dict.keys()):
            m_true = m_obs + delta_m
            spec = spectra_dict[m_obs]
            flux = getattr(spec, 'flux', spec)
            
            flux = np.asarray(flux, dtype=float)
            pix = np.arange(len(flux), dtype=float)
            if wave_calib is not None:
                xs = wave_calib.apply_to_pixel(pix, np.full_like(pix, m_true))
            else:
                xs = pix

            fig, ax = plt.subplots(figsize=(10, 4))
            ax.plot(xs, flux, label='Lamp Flux', color='k', linewidth=0.7)
            
            # 计算并绘制 SNR=50 的近似阈值基准线
            if len(flux) > 1:
                diffs = np.diff(flux)
                mad_noise = np.nanmedian(np.abs(diffs - np.nanmedian(diffs))) / (np.sqrt(2) * 0.6745)
                if mad_noise < 1e-6:
                    mad_noise = np.nanstd(flux)
                thr_snr = (snr_used or {}).get(m_obs, 50.0)
                level = np.nanmedian(flux) + thr_snr * mad_noise
                ax.axhline(y=level, color='r', linestyle=':', linewidth=1.2, alpha=0.6,
                           label=f'SNR={thr_snr:g} Threshold')

            if m_obs in detected_peaks and len(detected_peaks[m_obs]) > 0:
                p_x = detected_peaks[m_obs].astype(int)
                p_x = p_x[(p_x >= 0) & (p_x < len(flux))]
                p_x = p_x[np.argsort(flux[p_x])[::-1][:max_peaks]]
                ax.plot(xs[p_x], flux[p_x], 'r+', markersize=6, label=f'Top {max_peaks} Detected Peaks')
                
            # 绘制锚点匹配的标记 (向下短箭头 + 双行波长文本)
            if m_true in anchor_dict:
                max_flux = np.nanmax(flux) if len(flux) > 0 and np.nanmax(flux) > 0 else 1.0
                def _height(a):
                    i = int(np.round(a[0]))
                    return flux[i] if 0 <= i < len(flux) and np.isfinite(flux[i]) else -np.inf
                strongest = sorted(anchor_dict[m_true], key=_height, reverse=True)[:max_labels]
                gap = max_flux * 0.03
                line_len = max_flux * 0.08
                for x_obs, wave in strongest:
                    # 以峰附近局部最大值作为峰顶, 竖线与文字都放在峰正上方
                    c = int(np.round(x_obs))
                    lo, hi = max(c - 3, 0), min(c + 4, len(flux))
                    if hi <= lo:
                        continue
                    k = lo + int(np.nanargmax(flux[lo:hi]))
                    top = flux[k]
                    xk = float(xs[k])
                    ax.plot([xk, xk], [top + gap, top + gap + line_len], color='b', lw=1.0)
                    ax.text(xk, top + gap + line_len * 1.1, f"{wave:.4f}", color='b',
                            rotation=90, ha='center', va='bottom', fontsize=8)

            ax.set_title(f"True Order {m_true} (Aperture: {m_obs}, Offset: {delta_m})")
            ax.set_xlabel("Wavelength (Angstrom)" if wave_calib is not None else "Pixel (X)")
            ax.set_ylabel("Flux")
            
            # 留出顶部空间写波长文本，避免箭头和文字出界
            max_val = np.nanmax(flux) if len(flux) > 0 and np.nanmax(flux) > 0 else 1.0
            ax.set_ylim(bottom=0, top=max_val * 1.6)
            ax.legend(loc='upper right')
            
            plt.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

# 查找配置数据文件的实际路径
def _resolve_data_path(filepath: str) -> str:
    """Resolve calibration data path (supports running from arbitrary working directories)."""
    if not filepath:
        return filepath
        
    path = Path(filepath)
    # 如果是绝对路径或者在当前工作目录下能直接找到，直接返回
    if path.is_absolute() or path.exists():
        return str(path)
        
    # 如果在当前目录找不到，回退到源代码根目录寻找
    # __file__ 当前位于 src/specproc/core/wave_calibration.py
    pkg_dir = Path(__file__).parent.parent
    project_root = pkg_dir.parent.parent
    
    if (project_root / filepath).exists():
        return str(project_root / filepath)
        
    if (pkg_dir / filepath).exists():
        return str(pkg_dir / filepath)
        
    return filepath

# 完成灯谱识别、拟合与定标验证
def process_wavelength_stage(lamp_spectra: SpectraSet,
                             config: ConfigManager,
                             output_dir_base: str,
                             anchor_file: Optional[str] = None,
                             lamp_type: str = 'ThAr',
                             save_plots: bool = True,
                             fig_format: str = 'png',
                             lamp_name: str = 'wavelength_calibration') -> WaveCalib:
    """
    执行波长定标 (Step 7).

    Args:
        lamp_spectra: 1D 提取的灯谱 SpectraSet
        config: ConfigManager 实例
        output_dir_base: 基础输出目录
        anchor_file: 灯谱证认锚点 CSV 文件路径。为空则从 config 读取默认路径

    Returns:
        WaveCalib object
    """
    logger.info("Starting wavelength calibration stage...")
    calibrator = WavelengthCalibrator()

    # 1. 载入全量线库
    full_linelist_path = config.get('telescope.linelist', 'full_linelist', 'calib_data/linelists/thar-noao.dat')
    full_linelist = calibrator.load_line_list(_resolve_data_path(full_linelist_path))

    # 2. 从 1D 灯谱中提取观测 peaks
    detected_peaks = {}
    detected_snrs = {}
    spectra_dict = getattr(lamp_spectra, 'spectra', lamp_spectra) if not isinstance(lamp_spectra, dict) else lamp_spectra
    
    for m_obs, spectrum in spectra_dict.items():
        flux = spectrum.flux if hasattr(spectrum, 'flux') else spectrum
        # 统一使用恒定信噪比阈值 10.0 (过滤掉容易引起错配的微弱/混叠线)
        peaks, _, peak_snrs = calibrator.detect_lines_in_spectrum(flux, threshold=10.0)
        detected_peaks[m_obs] = peaks
        detected_snrs[m_obs] = peak_snrs

    use_precomputed = config.get_bool('telescope.linelist', 'use_precomputed_calibration', False)
    precomputed_file = config.get('telescope.linelist', 'precomputed_calib_file', '')
    
    prev_calib = None
    if use_precomputed and precomputed_file:
        precomputed_path = _resolve_data_path(precomputed_file)
        if Path(precomputed_path).exists():
            try:
                prev_calib = calibrator.load_calibration(precomputed_path)
                logger.info(f"Successfully loaded historical calibration. delta_m = {prev_calib.delta_m}.")
            except Exception as e:
                logger.warning(f"Failed to load precomputed calibration: {e}")
        else:
            logger.warning(f"Precomputed calibration file not found: {precomputed_path}")

    anchor_snr_used: Dict[int, float] = {}
    valid_detected_peaks: Dict[int, np.ndarray] = {}
    plot_anchors: List[Tuple[float, int, float]] = []
    plot_anchor_table: Optional[Dict[int, np.ndarray]] = None
    if prev_calib is not None:
        # 使用历史定标直接预测并全库匹配，完全跳过盲搜
        delta_m = prev_calib.delta_m
        pix_pos = []
        matched_wave = []
        wave_tol = config.get_float('reduce.wlcalib', 'rms_threshold', 0.5) * 4.0 # 预测容差可适当放宽
        
        for m_obs, peaks_x in detected_peaks.items():
            if len(peaks_x) == 0: continue
            m_true = m_obs + delta_m
            wave_pred = prev_calib.apply_to_pixel(peaks_x, np.full_like(peaks_x, m_true))
            
            for x_val, w_pred in zip(peaks_x, wave_pred):
                dist = np.abs(full_linelist - w_pred)
                idx = np.argmin(dist)
                if dist[idx] < wave_tol:
                    pix_pos.append([x_val, m_true])
                    matched_wave.append(full_linelist[idx])
                    
        pix_pos = np.array(pix_pos)
        matched_wave = np.array(matched_wave)
        logger.info(f"Matched {len(matched_wave)} lines using precomputed historical model.")
        matched_anchors = [] # 占位
        fit_peaks = detected_peaks
        
    else:
        # 3. 载入锚点文件；参考谱文件若为 (级次, X, 波长) 三列的证认表（无流量），则按锚点线表使用
        reference_spectrum = config.get('telescope.linelist', 'reference_spectrum', '').strip()
        anchors = {}
        anchor_xtable = None
        if reference_spectrum:
            ref_table = np.loadtxt(_resolve_data_path(reference_spectrum), ndmin=2)
            if ref_table.shape[1] >= 3 and ref_table[:, 2].max() > 100.0:
                logger.info(f"Reference file {reference_spectrum} lists identified lines "
                            "(order, X, wavelength); using it as the anchor table.")
                anchor_xtable = ref_table[:, :3]
                for m_ref in np.unique(ref_table[:, 0].astype(int)):
                    anchors[int(m_ref)] = np.sort(ref_table[ref_table[:, 0].astype(int) == m_ref, 2])
                reference_spectrum = ''
        if not reference_spectrum and not anchors:
            if anchor_file is None:
                anchor_file = config.get('telescope.linelist', 'anchor_file', 'calib_data/telescopes/xinglong216hrs/xinglong216_thar_anchors.csv')
            anchors = calibrator.load_anchor_file(_resolve_data_path(anchor_file))

        # 读取剔除最高/最低级次配置
        discard_top = config.get_int('reduce.wlcalib', 'discard_highest_apertures', 0)
        discard_bottom = config.get_int('reduce.wlcalib', 'discard_lowest_apertures', 0)
        all_m_obs = sorted(list(detected_peaks.keys()))
        
        start_idx = discard_bottom
        end_idx = len(all_m_obs) - discard_top
        if start_idx < end_idx:
            valid_m_obs = all_m_obs[start_idx:end_idx]
            if discard_bottom > 0 or discard_top > 0:
                logger.info(f"Discarding {discard_bottom} bottom and {discard_top} top apertures. Using {len(valid_m_obs)} for anchors.")
        else:
            valid_m_obs = all_m_obs
            logger.warning("Discard parameters too large, ignoring discard configuration.")

        # 锚点盲配：每个级次取最亮的前 10 个峰。优先要求 SNR > anchor_snr；若某级次此类峰
        # 少于 anchor_min_lines（暗级次），则放宽到 SNR > anchor_snr_floor，避免整级次被丢弃。
        anchor_snr = config.get_float('reduce.wlcalib', 'anchor_snr', 50.0)
        anchor_floor = config.get_float('reduce.wlcalib', 'anchor_snr_floor', 15.0)
        anchor_min = config.get_int('reduce.wlcalib', 'anchor_min_lines', 5)
        logger.info(f"Using top 10 detected peaks per aperture for anchor matching "
                    f"(SNR > {anchor_snr:g}; relaxed to > {anchor_floor:g} when fewer than {anchor_min} lines).")
        peaks_for_matching = {}
        all_peaks_for_matching = {}
        for m_obs in valid_m_obs:
            if m_obs in detected_peaks:
                snrs = detected_snrs[m_obs]
                thr = anchor_snr if np.sum(snrs > anchor_snr) >= anchor_min else min(anchor_snr, anchor_floor)
                mask = snrs > thr
                if np.sum(mask) > 0:
                    valid_peaks = detected_peaks[m_obs][mask]
                    valid_snrs = snrs[mask]
                    top_indices = np.sort(np.argsort(valid_snrs)[-min(10, len(valid_peaks)):])
                    peaks_for_matching[m_obs] = valid_peaks[top_indices]
                    all_peaks_for_matching[m_obs] = valid_peaks
                    anchor_snr_used[m_obs] = thr
                    if thr < anchor_snr:
                        logger.info(f"Aperture {m_obs}: fewer than {anchor_min} peaks above SNR {anchor_snr:g}; "
                                    f"using SNR > {thr:g} ({len(top_indices)} peaks).")

        # 4. 匹配锚点寻找最佳级次偏移 delta_m
        if reference_spectrum:
            delta_m, matched_anchors = calibrator.match_reference_spectrum(
                {m: spectra_dict[m] for m in valid_m_obs}, _resolve_data_path(reference_spectrum))
        else:
            order_offset = config.get_int('reduce.wlcalib', 'order_offset', -1)
            if anchor_xtable is not None:
                # The table gives X positions: match by position instead of wavelength patterns.
                delta_m, matched_anchors = calibrator.match_anchor_positions(
                    peaks_for_matching, all_peaks_for_matching, anchor_xtable,
                    fixed_offset=order_offset if order_offset >= 0 else None,
                    tol=config.get_float('reduce.wlcalib', 'anchor_match_tol', 3.0))
            else:
                delta_m, matched_anchors = calibrator.find_order_offset_and_match(
                    peaks_for_matching, anchors, fixed_offset=order_offset if order_offset >= 0 else None)

        # 5. 基于锚点拟合 2D 粗定标模型
        anchor_pix = np.array([[pt[0], pt[1]] for pt in matched_anchors])
        anchor_wave = np.array([pt[2] for pt in matched_anchors])
        poly_type = config.get('reduce.wlcalib', 'poly_type', 'chebyshev')
        
        # Dense reference-derived anchors support the full polynomial order; sparse ones need a low-order rough fit.
        rough_cap = 99 if reference_spectrum else 3
        rough_xorder = min(config.get_int('reduce.wlcalib', 'xorder', 4), rough_cap)
        rough_yorder = min(config.get_int('reduce.wlcalib', 'yorder', 4), rough_cap)
        max_res = config.get_float('reduce.wlcalib', 'rms_threshold', 0.5)
        
        logger.info(f"Fitting initial rough 2D polynomial (orders: {rough_xorder}, {rough_yorder}) using {len(anchor_wave)} anchor points...")
        rough_calib = calibrator.fit_wavelength_polynomial(
            anchor_pix, anchor_wave, 
            xorder=rough_xorder, yorder=rough_yorder, poly_type=poly_type,
            max_residual=max_res * 5.0  # Allow larger residuals for rough fit
        )

        # 6. 使用粗定标模型预测波长，并进行全库匹配
        match_tol = config.get_float('reduce.wlcalib', 'match_tolerance', 2.0)
        # 全库匹配：放开手脚，使用信噪比 > 5 的所有检测线
        valid_detected_peaks = {m: detected_peaks[m] for m in valid_m_obs if m in detected_peaks}
        fit_peaks = valid_detected_peaks
        pix_pos, matched_wave = calibrator.match_full_catalog( 
            valid_detected_peaks, rough_calib, delta_m, full_linelist, tolerance=match_tol
        )
        logger.info(f"Successfully matched {len(matched_wave)} lines from the full catalog (tolerance={match_tol:.2f} A).")

        plot_anchors = getattr(calibrator, 'last_peak_anchors', None) or matched_anchors
        if not reference_spectrum and anchors:
            plot_anchor_table = anchors

    # 7. 全局精确 2D 拟合
    x_order = config.get_int('reduce.wlcalib', 'xorder', 4)
    y_order = config.get_int('reduce.wlcalib', 'yorder', 4)
    min_required = (x_order + 1) * (y_order + 1)

    if len(matched_wave) < min_required:
        logger.warning(f"Only {len(matched_wave)} lines matched from full catalog (need {min_required}). "
                       "Falling back to using only matched anchor points.")
        if len(matched_anchors) >= min_required:
            pix_pos = np.array([[pt[0], pt[1]] for pt in matched_anchors])
            matched_wave = np.array([pt[2] for pt in matched_anchors])
            logger.info(f"Using {len(matched_wave)} anchor points for final 2D fit.")
        else:
            raise RuntimeError(f"Not enough lines matched from either catalog ({len(matched_wave)}) or anchors ({len(matched_anchors)}) "
                               f"for 2D polynomial fitting (requires {min_required}).")

    poly_type = config.get('reduce.wlcalib', 'poly_type', 'chebyshev')
    wave_calib = calibrator.fit_wavelength_polynomial(
        pix_pos, matched_wave, xorder=x_order, yorder=y_order, poly_type=poly_type, max_residual=max_res
    )
    wave_calib.delta_m = delta_m

    fit_mode = config.get('reduce.wlcalib', 'fit_mode', 'global').strip().lower()
    if fit_mode == 'per_order':
        order_kwargs = dict(degree=config.get_int('reduce.wlcalib', 'order_degree', 4),
                            min_lines=config.get_int('reduce.wlcalib', 'order_min_lines', 8),
                            min_span=config.get_float('reduce.wlcalib', 'order_min_span', 0.6),
                            max_order_rms=config.get_float('reduce.wlcalib', 'order_max_rms', 0.0),
                            max_residual=max_res)
        calibrator.fit_per_order(wave_calib, **order_kwargs)

        # Independently fitted orders predict wavelengths far better than the first rough model,
        # so rematch the catalog with a tight tolerance to drop wrong matches in crowded orders.
        refine_tol = config.get_float('reduce.wlcalib', 'refine_tolerance', 0.1)
        for it in range(config.get_int('reduce.wlcalib', 'refine_iterations', 2)):
            if not wave_calib.order_coefs:
                break
            new_pix, new_wave = calibrator.match_full_catalog(
                fit_peaks, wave_calib, delta_m, full_linelist,
                tolerance=config.get_float('reduce.wlcalib', 'match_tolerance', 2.0),
                tolerance_fitted=refine_tol)
            if len(new_wave) < min_required:
                break
            logger.info(f"Refinement {it + 1}: {len(new_wave)} lines matched (tolerance {refine_tol:.2f} A "
                        "for independently fitted orders).")
            wave_calib = calibrator.fit_wavelength_polynomial(
                new_pix, new_wave, xorder=x_order, yorder=y_order, poly_type=poly_type, max_residual=max_res)
            wave_calib.delta_m = delta_m
            calibrator.fit_per_order(wave_calib, **order_kwargs)

    # 保存校准结果
    calib_file = Path(output_dir_base) / 'step7_wavelength' / f'wcal_{lamp_name}.fits'
    calib_file.parent.mkdir(parents=True, exist_ok=True)
    calibrator.save_calibration(str(calib_file))

    # 锚点匹配诊断图: 横坐标用最终定标后的波长
    if plot_anchor_table:
        # 盲匹配的锚点可能含误配: 用最终定标重新把锚点线表匹配到实测 peaks
        plot_anchors = _rematch_anchor_table(valid_detected_peaks, wave_calib, delta_m, plot_anchor_table,
                                             tol=config.get_float('reduce.wlcalib', 'anchor_match_tol_wave', 0.5))
    if save_plots and (plot_anchors or valid_detected_peaks):
        pdf_path = calib_file.parent / f'wcal_{lamp_name}_matched_anchors.pdf'
        try:
            _plot_matched_anchors_pdf(lamp_spectra, valid_detected_peaks, plot_anchors, delta_m,
                                      str(pdf_path), snr_used=anchor_snr_used, wave_calib=wave_calib)
            logger.info(f"Saved matched anchors diagnostic PDF to {pdf_path}")
        except Exception as e:
            logger.warning(f"Could not generate matched anchors plot: {e}")

    # 绘制验证图并保存
    if save_plots and getattr(wave_calib, 'line_pixels', None) is not None:
        out_dir = calib_file.parent
        plot_file = out_dir / f'wcal_{lamp_name}.{fig_format}'
        surf_plot_file = out_dir / f'wcal_{lamp_name}_surface.{fig_format}'
        try:
            _plot_calib_diagnostic(wave_calib, str(plot_file), max_res)
            _plot_calib_surface_diagnostic(wave_calib, str(surf_plot_file), max_res)
        except Exception as e:
            logger.warning(f"Could not generate wavelength calibration plot: {e}")

    return wave_calib
