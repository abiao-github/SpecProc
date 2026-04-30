"""
Wavelength calibration stage of spectral reduction pipeline.

Handles ThAr line identification, polynomial fitting,
and wavelength solution computation.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import cm
import logging
from typing import Tuple, Optional, List, Dict
from scipy.optimize import curve_fit
from pathlib import Path
from src.core.data_structures import WaveCalib, Spectrum, SpectraSet
from src.config.config_manager import ConfigManager
from src.utils.fits_io import read_fits_image
from src.plotting.spectra_plotter import plot_wavelength_calibration

logger = logging.getLogger(__name__)


class WavelengthCalibrator:
    """Handles wavelength calibration from reference lines."""

    def __init__(self):
        """Initialize wavelength calibrator."""
        self.wave_calib = None
        self.reference_lines = {}

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

    def detect_lines_in_spectrum(self, spectrum: np.ndarray,
                                wavelength: Optional[np.ndarray] = None,
                                threshold: float = 3.0) -> Tuple[np.ndarray, np.ndarray]:
        """
        Detect emission lines in spectrum.

        Args:
            spectrum: 1D flux array
            wavelength: Optional wavelength array (for initial guess)
            threshold: DetectionThreshold in sigma

        Returns:
            Tuple of (pixel_positions, line_strengths)
        """
        # Find peaks in spectrum
        from scipy.signal import find_peaks

        mean = np.mean(spectrum)
        std = np.std(spectrum)
        peaks, properties = find_peaks(spectrum, height=threshold * std, distance=5)

        logger.debug(f"Detected {len(peaks)} emission lines")

        return peaks, spectrum[peaks]

    def load_anchor_file(self, anchor_file: str) -> Dict[int, np.ndarray]:
        """
        读取灯谱证认的发射线锚点文件。
        预期仅包含两列信息：Order(级次号), Wavelength(波长)。
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

    def find_order_offset_and_match(self, detected_peaks: Dict[int, np.ndarray], anchors: Dict[int, np.ndarray]) -> Tuple[int, List[Tuple[float, int, float]]]:
        """
        寻找最佳级次偏移 (delta_m = m_ref - m_obs) 并返回成功匹配的锚点 (X_obs, m_true, wavelength)。
        两阶段盲搜法：
        1. 仅使用中间级次，显式区分正反色散方向搜索最佳 delta_m。
        2. 利用确定的方向和偏移量，对全量级次进行秒级匹配。
        """
        all_m_obs = sorted(list(detected_peaks.keys()))
        if not all_m_obs:
            return 0, []

        # ====================================================================
        # Phase 1: 选取中间级次进行快速盲搜，确定最佳偏移量和色散方向
        # ====================================================================
        n_orders = len(all_m_obs)
        start_idx = n_orders // 4          # 25% 位置
        end_idx = n_orders * 3 // 4        # 75% 位置
        if end_idx - start_idx < 3:        # 如果总级次太少，则启动兜底保护
            start_idx, end_idx = 0, n_orders
        middle_orders = all_m_obs[start_idx:end_idx]

        best_offset = 0
        best_direction = 1  # 1: Left-to-Right, -1: Right-to-Left
        max_total_matches = -1

        logger.info(f"Phase 1: Blind matching on middle {len(middle_orders)} orders to find delta_m and dispersion direction...")

        for direction in [1, -1]:
            for delta_m in range(-40, 41):
                total_matches = 0
                for m_obs in middle_orders:
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
                                    if not (0.005 < abs(a_test) < 0.15): continue
                                    
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
            logger.warning("盲配模式匹配失败！未找到合适的级次偏移。请检查峰值提取、锚点文件或放宽色散率限制。")
            return 0, []
            
        dir_str = "Left-to-Right (+)" if best_direction == 1 else "Right-to-Left (-)"
        logger.info(f"Phase 1 complete: Found delta_m={best_offset}, dispersion direction={dir_str}")

        # ====================================================================
        # Phase 2: 使用确定的最佳偏移量和色散方向，全量级次提取匹配点
        # ====================================================================
        logger.info(f"Phase 2: Extracting matched anchors across all orders using delta_m={best_offset}...")
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
                            if not (0.005 < abs(a_test) < 0.15): continue
                            
                            b = W[k] - a_test * (X[i] if best_direction == 1 else X[j])
                                
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
                global_total_matches += len(best_inliers)
                best_global_matched_points.extend(best_inliers)

        logger.info(f"Blind match fully complete: delta_m={best_offset} with {global_total_matches} lines matched across all orders.")
            
        return best_offset, best_global_matched_points

    def rough_calibration(self, matched_points: List[Tuple[float, int, float]]) -> np.ndarray:
        """
        粗定标拟合全局方程: m * lambda = a * X^2 + b * X + c
        """
        if len(matched_points) < 3:
            raise ValueError(f"匹配的锚点数量 ({len(matched_points)}) 太少，无法执行粗校准(最少3点)")
        
        X_arr = np.array([pt[0] for pt in matched_points])
        m_arr = np.array([pt[1] for pt in matched_points])
        wave_arr = np.array([pt[2] for pt in matched_points])
        
        m_lambda = m_arr * wave_arr
        # 拟合 2 阶多项式
        coeffs = np.polyfit(X_arr, m_lambda, 2)
        logger.info(f"Rough calibration m*lambda = aX^2 + bX + c coefficients: {coeffs}")
        return coeffs

    def match_full_catalog(self, detected_peaks: Dict[int, np.ndarray], rough_coeffs: np.ndarray, delta_m: int, full_linelist: np.ndarray, tolerance: float = 0.5) -> Tuple[np.ndarray, np.ndarray]:
        """
        利用拟合的全局粗校准方程，预测所有 detected_peaks 的波长并匹配全量线库。
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
            
            # 利用粗定标预测波长: lambda = (aX^2 + bX + c) / m
            m_lambda_pred = np.polyval(rough_coeffs, peaks_x)
            wave_pred = m_lambda_pred / m_true
            
            # 全库匹配寻找最近线
            for x_val, w_pred in zip(peaks_x, wave_pred):
                dist = np.abs(full_linelist - w_pred)
                idx = np.argmin(dist)
                if dist[idx] < tolerance:
                    pix_pos.append([x_val, m_true])
                    matched_wave.append(full_linelist[idx])
                    
        return np.array(pix_pos), np.array(matched_wave)

    def fit_wavelength_polynomial(self, pixel_positions: np.ndarray,
                                 wavelengths: np.ndarray,
                                 xorder: int = 4, yorder: int = 4,
                                 poly_type: str = 'chebyshev') -> WaveCalib:
        """
        Fit 2D wavelength polynomial.

        Args:
            pixel_positions: Array of (x, y) pixel coordinates
            wavelengths: Corresponding wavelengths
            xorder: Polynomial order in X
            yorder: Polynomial order in Y
            poly_type: Polynomial basis ('chebyshev', 'legendre' or 'polynomial')

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
                
                new_keep = np.abs(residuals) < 3.0 * rms
                if np.array_equal(new_keep, keep):
                    break
                keep = new_keep
                
            n_rejected = len(keep) - np.sum(keep)
            if n_rejected > 0:
                logger.info(f"Sigma clipping rejected {n_rejected} outliers.")

            pred_m_lambda = A @ coeffs_1d
            pred_lambda = pred_m_lambda / y_pix
            rms = np.sqrt(np.mean((pred_lambda[keep] - wavelengths[keep]) ** 2))

            logger.info(f"Wavelength fit RMS: {rms:.4f} Angstrom (on {np.sum(keep)} lines)")

            # Reshape coefficients into 2D array
            poly_coef = coeffs_1d.reshape((xorder + 1, yorder + 1))

            # Create WaveCalib object
            wave_calib = WaveCalib(
                poly_coef=poly_coef,
                xorder=xorder,
                yorder=yorder,
                poly_type=poly_type,
                domain_x=domain_x,
                domain_y=domain_y,
                line_pixels=x_pix,
                line_orders=y_pix,
                line_catalog=wavelengths,
                rms=rms,
                nlines=len(wavelengths),
                calib_type='ThAr'
            )

            self.wave_calib = wave_calib
            return wave_calib

        except Exception as e:
            logger.error(f"Error fitting wavelength polynomial: {e}")
            raise

    def apply_wavelength_calibration(self, spectrum: np.ndarray,
                                    pixel_array: np.ndarray,
                                    aperture_y: float = 0.0) -> np.ndarray:
        """
        Apply wavelength calibration to spectrum pixels.

        Args:
            spectrum: 1D spectrum array
            pixel_array: Pixel coordinate array
            aperture_y: Y coordinate (order number) of aperture

        Returns:
            Wavelength array
        """
        if self.wave_calib is None:
            raise RuntimeError("No wavelength calibration available")

        # Apply polynomial
        wavelength = self.wave_calib.apply_to_pixel(pixel_array, aperture_y * np.ones_like(pixel_array))

        return wavelength

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

        hdul.writeto(str(output_path), overwrite=True)
        logger.info(f"Saved wavelength calibration to {output_path}")


def _plot_calib_diagnostic(wave_calib: WaveCalib, plot_file: str):
    """生成波长定标诊断图（散点图及拟合残差）。"""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True, gridspec_kw={'height_ratios': [3, 1]})
    
    x_pix = wave_calib.line_pixels
    m_orders = wave_calib.line_orders
    wave_ref = wave_calib.line_catalog
    wave_pred = wave_calib.apply_to_pixel(x_pix, m_orders)
    residuals = wave_pred - wave_ref
    
    unique_orders = np.unique(m_orders)
    cmap = cm.get_cmap('turbo')
    
    for i, m in enumerate(unique_orders):
        mask = (m_orders == m)
        color = cmap(i / max(1, len(unique_orders) - 1))
        
        # 上图：拟合曲线与参考点
        label_str = f'Order {int(m)}' if len(unique_orders) <= 15 else None
        ax1.scatter(x_pix[mask], wave_ref[mask], color=color, s=15, alpha=0.8, label=label_str)
        x_smooth = np.linspace(x_pix[mask].min(), x_pix[mask].max(), 100)
        wave_smooth = wave_calib.apply_to_pixel(x_smooth, np.full_like(x_smooth, m))
        ax1.plot(x_smooth, wave_smooth, color=color, alpha=0.5)
        
        # 下图：残差
        ax2.scatter(x_pix[mask], residuals[mask], color=color, s=15, alpha=0.8)
        
    ax1.set_ylabel(r'Wavelength ($\AA$)')
    ax1.set_title(fr'Wavelength Calibration (RMS = {wave_calib.rms:.4f} $\AA$, N = {wave_calib.nlines})')
    if len(unique_orders) <= 15:
        ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize='small')
    ax1.grid(True, alpha=0.3)
    
    ax2.axhline(0, color='k', linestyle='--', alpha=0.5)
    ax2.set_xlabel('Pixel (X)')
    ax2.set_ylabel(r'Residual ($\AA$)')
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(plot_file, dpi=150)
    plt.close(fig)

def _plot_calib_surface_diagnostic(wave_calib: WaveCalib, plot_file: str):
    """生成二维波长定标拟合曲面和残差诊断图。"""
    fig = plt.figure(figsize=(14, 6))
    
    # 1. 3D 拟合曲面 (Wavelength Surface)
    ax1 = fig.add_subplot(121, projection='3d')
    x_pix = wave_calib.line_pixels
    m_orders = wave_calib.line_orders
    wave_ref = wave_calib.line_catalog
    
    # 生成用于绘制平滑曲面的网格
    x_grid = np.linspace(x_pix.min(), x_pix.max(), 50)
    m_grid = np.linspace(m_orders.min(), m_orders.max(), 50)
    X, M = np.meshgrid(x_grid, m_grid)
    WAVE = wave_calib.apply_to_pixel(X, M)
    
    # 绘制预测曲面和实际匹配散点
    ax1.plot_surface(X, M, WAVE, cmap='viridis', alpha=0.6, edgecolor='none')
    ax1.scatter(x_pix, m_orders, wave_ref, color='r', s=10, alpha=0.8, label='Matched Lines')
    ax1.set_xlabel('Pixel (X)')
    ax1.set_ylabel('Order (m)')
    ax1.set_zlabel(r'Wavelength ($\AA$)')
    ax1.set_title(f'2D Wavelength Solution Surface\n({wave_calib.poly_type.capitalize()} Poly, {wave_calib.xorder}x{wave_calib.yorder})')
    
    # 2. 拟合残差图 (二维散点，颜色区分级次)
    ax2 = fig.add_subplot(122)
    residuals = wave_calib.apply_to_pixel(x_pix, m_orders) - wave_ref
    scatter = ax2.scatter(x_pix, residuals, c=m_orders, cmap='turbo', s=15, alpha=0.8)
    ax2.axhline(0, color='k', linestyle='--', alpha=0.5)
    ax2.set_xlabel('Pixel (X)')
    ax2.set_ylabel(r'Residual ($\AA$)')
    ax2.set_title(fr'Fitting Residuals (RMS = {wave_calib.rms:.4f} $\AA$)')
    ax2.grid(True, alpha=0.3)
    fig.colorbar(scatter, ax=ax2, label='Order (m)')
    
    plt.tight_layout()
    plt.savefig(plot_file, dpi=150)
    plt.close(fig)

def _plot_matched_anchors_pdf(lamp_spectra, detected_peaks: Dict[int, np.ndarray], 
                              matched_anchors: List[Tuple[float, int, float]], delta_m: int, 
                              pdf_path: str):
    """将每个级次的光谱及匹配到的锚点画入一个多页 PDF 文件。"""
    with PdfPages(pdf_path) as pdf:
        matched_dict = {}
        for x_obs, m_true, wave in matched_anchors:
            if m_true not in matched_dict:
                matched_dict[m_true] = []
            matched_dict[m_true].append((x_obs, wave))
            
        spectra_dict = getattr(lamp_spectra, 'spectra', lamp_spectra) if not isinstance(lamp_spectra, dict) else lamp_spectra
        
        for m_obs in sorted(spectra_dict.keys()):
            m_true = m_obs + delta_m
            spec = spectra_dict[m_obs]
            flux = getattr(spec, 'flux', spec)
            
            fig, ax = plt.subplots(figsize=(10, 4))
            ax.plot(flux, label='Lamp Flux', color='k', linewidth=0.7)
            
            if m_obs in detected_peaks and len(detected_peaks[m_obs]) > 0:
                p_x = detected_peaks[m_obs].astype(int)
                ax.plot(p_x, flux[p_x], 'r+', markersize=6, label='Detected Peaks')
                
            if m_true in matched_dict:
                max_flux = np.nanmax(flux) if len(flux) > 0 and np.nanmax(flux) > 0 else 1.0
                for x_obs, wave in matched_dict[m_true]:
                    ax.axvline(x=x_obs, color='b', linestyle='--', linewidth=0.8, alpha=0.7)
                    ax.text(x_obs, max_flux * 1.05, f"{wave:.3f}", color='b', 
                            rotation=90, va='bottom', ha='center', fontsize=8)
                            
            ax.set_title(f"Order {m_true} (Original: {m_obs}, Offset: {delta_m})")
            ax.set_xlabel("Pixel (X)")
            ax.set_ylabel("Flux")
            
            # 留出顶部 30% 空间写波长文本
            max_val = np.nanmax(flux) if len(flux) > 0 and np.nanmax(flux) > 0 else 1.0
            ax.set_ylim(bottom=0, top=max_val * 1.3)
            ax.legend(loc='upper right')
            
            plt.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

def process_wavelength_stage(lamp_spectra: SpectraSet,
                             config: ConfigManager,
                             output_dir_base: str,
                             anchor_file: Optional[str] = None,
                             lamp_type: str = 'ThAr',
                             save_plots: bool = True,
                             fig_format: str = 'png') -> WaveCalib:
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
    full_linelist = calibrator.load_line_list(full_linelist_path)

    # 2. 从 1D 灯谱中提取观测 peaks
    detected_peaks = {}
    detected_fluxes = {}
    spectra_dict = getattr(lamp_spectra, 'spectra', lamp_spectra) if not isinstance(lamp_spectra, dict) else lamp_spectra
    clipping_threshold = config.get_float('reduce.wlcalib', 'clipping', 3.0)
    
    for m_obs, spectrum in spectra_dict.items():
        flux = spectrum.flux if hasattr(spectrum, 'flux') else spectrum
        peaks, peak_fluxes = calibrator.detect_lines_in_spectrum(flux, threshold=clipping_threshold)
        detected_peaks[m_obs] = peaks
        detected_fluxes[m_obs] = peak_fluxes

    # 3. 载入锚点文件
    if anchor_file is None:
        anchor_file = config.get('telescope.linelist', 'anchor_file', 'calib_data/telescopes/xinglong216hrs/xinglong_thar_lines.csv')
        
    anchors = calibrator.load_anchor_file(anchor_file)

    # 提取最明显的 peaks 参与锚点盲配，加快速度并降低弱线干扰
    # 保留数量设为单级次内最大锚点数的 1.5 倍，且至少保留 15 个
    max_anchors = max([len(pts) for pts in anchors.values()]) if anchors else 10
    top_n = max(int(max_anchors * 1.5), 15)
    
    prominent_peaks = {}
    for m_obs, peaks_x in detected_peaks.items():
        fluxes = detected_fluxes[m_obs]
        if len(peaks_x) > top_n:
            top_indices = np.argsort(fluxes)[-top_n:]
            prominent_peaks[m_obs] = peaks_x[top_indices]
        else:
            prominent_peaks[m_obs] = peaks_x

    # 4. 匹配锚点寻找最佳级次偏移 delta_m (仅使用过滤后的强峰)
    delta_m, matched_anchors = calibrator.find_order_offset_and_match(prominent_peaks, anchors)
    logger.info(f"Determined order offset (delta_m = m_ref - m_obs): {delta_m}")

    # 输出匹配锚点诊断图 PDF
    if save_plots:
        pdf_path = Path(output_dir_base) / 'step7_wavelength' / 'matched_anchors_diagnostic.pdf'
        pdf_path.parent.mkdir(parents=True, exist_ok=True)
        _plot_matched_anchors_pdf(lamp_spectra, prominent_peaks, matched_anchors, delta_m, str(pdf_path))
        logger.info(f"Saved matched anchors diagnostic PDF to {pdf_path}")

    # 5. 基于锚点拟合粗定标方程 (m * lambda = a*X^2 + b*X + c)
    rough_coeffs = calibrator.rough_calibration(matched_anchors)

    # 6. 使用粗定标预测波长，并进行全库匹配
    wave_tol = config.get_float('reduce.wlcalib', 'rms_threshold', 0.5)
    pix_pos, matched_wave = calibrator.match_full_catalog(
        detected_peaks, rough_coeffs, delta_m, full_linelist, tolerance=wave_tol
    )
    logger.info(f"Successfully matched {len(matched_wave)} lines from the full catalog.")

    if len(matched_wave) < 15:
        logger.warning("Very few lines matched to the full catalog. Falling back to using only matched anchor points.")
        if len(matched_anchors) >= 15:
            pix_pos = np.array([[pt[0], pt[1]] for pt in matched_anchors])
            matched_wave = np.array([pt[2] for pt in matched_anchors])
            logger.info(f"Using {len(matched_wave)} anchor points for final 2D fit.")
        else:
            raise RuntimeError("Not enough lines matched (from either catalog or anchors) for 2D polynomial fitting.")

    # 7. 全局精确 2D 拟合
    x_order = config.get_int('reduce.wlcalib', 'xorder', 4)
    y_order = config.get_int('reduce.wlcalib', 'yorder', 4)
    poly_type = config.get('reduce.wlcalib', 'poly_type', 'chebyshev')
    wave_calib = calibrator.fit_wavelength_polynomial(pix_pos, matched_wave, xorder=x_order, yorder=y_order, poly_type=poly_type)
    wave_calib.delta_m = delta_m

    # 保存校准结果
    calib_file = Path(output_dir_base) / 'step7_wavelength' / 'wavelength_calibration.fits'
    calib_file.parent.mkdir(parents=True, exist_ok=True)
    calibrator.save_calibration(str(calib_file))

    # 绘制验证图并保存
    if save_plots and getattr(wave_calib, 'line_pixels', None) is not None:
        out_dir = calib_file.parent
        plot_file = out_dir / f'wavelength_calibration.{fig_format}'
        surf_plot_file = out_dir / f'wavelength_calibration_surface.{fig_format}'
        try:
            _plot_calib_diagnostic(wave_calib, str(plot_file))
            _plot_calib_surface_diagnostic(wave_calib, str(surf_plot_file))
        except Exception as e:
            logger.warning(f"Could not generate wavelength calibration plot: {e}")

    return wave_calib
