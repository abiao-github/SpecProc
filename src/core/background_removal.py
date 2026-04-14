"""
Background subtraction stage of spectral reduction pipeline.

Handles inter-order background (stray light) estimation and removal.
"""

import numpy as np
import logging
from pathlib import Path
from typing import Tuple, Optional
from scipy.ndimage import median_filter
import matplotlib.pyplot as plt
from src.utils.image_processing import estimate_background_2d
from src.config.config_manager import ConfigManager
from src.core.data_structures import ApertureSet
from src.plotting.spectra_plotter import plot_background_residuals
from src.utils.fits_io import write_fits_image

logger = logging.getLogger(__name__)


class BackgroundRemover:
    """Handles background/stray light estimation and correction."""

    def __init__(self, config: ConfigManager):
        """
        Initialize background remover.

        Args:
            config: Configuration manager
        """
        self.config = config
        self.background_2d = None

    def _build_aperture_mask(self, image_shape: Tuple[int, int], apertures: ApertureSet,
                             margin_pixels: int = 3) -> np.ndarray:
        """Build boolean mask using traced lower/upper edges plus an outer margin."""
        rows, cols = image_shape
        mask = np.zeros((rows, cols), dtype=bool)
        x_coords = np.arange(cols)

        for aperture in apertures.apertures.values():
            center = aperture.get_position(x_coords)
            lower = aperture.get_lower(x_coords)
            upper = aperture.get_upper(x_coords)

            width_guess = float(aperture.width if aperture.width is not None else 8.0)
            half_w_max = max(4.0, min(0.12 * rows, 2.2 * width_guess))

            for x in x_coords:
                c = float(center[x])
                lo = float(lower[x])
                hi = float(upper[x])

                if (not np.isfinite(c)) or c < -margin_pixels or c > (rows - 1 + margin_pixels):
                    continue
                if (not np.isfinite(lo)) or (not np.isfinite(hi)):
                    continue

                y_lo = min(lo, hi) - margin_pixels
                y_hi = max(lo, hi) + margin_pixels

                # Guard against pathological edge fits that produce unrealistically
                # wide masks and wipe out inter-order regions.
                if (y_hi - y_lo) > 2.0 * half_w_max:
                    y_lo = c - half_w_max
                    y_hi = c + half_w_max

                y1 = max(0, int(np.floor(y_lo)))
                y2 = min(rows, int(np.ceil(y_hi + 1)))
                if y2 > y1:
                    mask[y1:y2, x] = True

        return mask

    def _estimate_background_split_halves(self, image: np.ndarray, order_mask: np.ndarray,
                                          order: int, method: str,
                                          smooth_sigma: float, sigma_clip: float,
                                          maxiters: int, bspline_smooth: float
                                          ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Estimate background separately for upper/lower detector halves."""
        rows, cols = image.shape
        split = rows // 2
        background = np.zeros_like(image, dtype=np.float64)
        model_lower = np.zeros((split, cols), dtype=np.float64)
        model_upper = np.zeros((rows - split, cols), dtype=np.float64)

        halves = [
            (0, split, "lower"),
            (split, rows, "upper"),
        ]
        for y0, y1, label in halves:
            sub = image[y0:y1, :].astype(np.float64)
            sub_mask = order_mask[y0:y1, :]
            valid = ~sub_mask

            sampled = sub.copy()
            sampled[sub_mask] = np.nan

            if np.mean(valid) < 0.03:
                # Fallback in very dense order region: smooth heavily on unmasked
                # values while leaving masked zones as NaN constraints.
                logger.warning(
                    f"Inter-order fraction too small in {label} half ({np.mean(valid):.4f}); "
                    "using robust median-smoothed fallback background."
                )
                fill = np.nanmedian(sampled)
                if not np.isfinite(fill):
                    fill = float(np.nanmedian(sub))
                tmp = np.where(np.isfinite(sampled), sampled, fill)
                sub_bg = median_filter(tmp, size=31).astype(np.float64)
            else:
                sub_bg = estimate_background_2d(
                    sampled,
                    order=order,
                    valid_mask=valid,
                    method=method,
                    smooth_sigma=smooth_sigma,
                    sigma_clip=sigma_clip,
                    maxiters=maxiters,
                    bspline_smooth=bspline_smooth,
                ).astype(np.float64)

            background[y0:y1, :] = sub_bg
            if label == "lower":
                model_lower[:, :] = sub_bg
            else:
                model_upper[:, :] = sub_bg

        return background, model_lower, model_upper

    def estimate_background(self, image: np.ndarray, order: int = 2) -> np.ndarray:
        """
        Estimate 2D background using polynomial fitting.

        Args:
            image: 2D image array
            order: Polynomial order

        Returns:
            Background model
        """
        logger.info(f"Estimating 2D background (polyorder={order})...")

        background = estimate_background_2d(image, order=order)
        self.background_2d = background

        logger.info(f"Background estimated: min={np.min(background):.1f}, "
                   f"max={np.max(background):.1f}")

        return background

    def estimate_background_interorder(self, image: np.ndarray, apertures: ApertureSet,
                                       order: int = 2,
                                       mask_margin_pixels: int = 3,
                                       method: str = 'chebyshev',
                                       smooth_sigma: float = 20.0,
                                       sigma_clip: float = 3.0,
                                       maxiters: int = 4,
                                       bspline_smooth: float = 1.0) -> Tuple[np.ndarray, np.ndarray]:
        """
        Estimate 2D background only from inter-order pixels.

        Steps:
        1) Build aperture mask from traced orders.
        2) Keep only inter-order pixels as fit constraints.
        3) Fit smooth 2D background surface.
        """
        logger.info(
            f"Estimating inter-order background (method={method}, order={order}, margin_px={mask_margin_pixels})..."
        )

        order_mask = self._build_aperture_mask(image.shape, apertures, margin_pixels=mask_margin_pixels)
        inter_order = ~order_mask

        sampled = image.astype(np.float64).copy()
        sampled[order_mask] = np.nan

        # Use inter-order pixels directly as constraints with robust clipping.
        background = estimate_background_2d(
            sampled,
            order=order,
            valid_mask=inter_order,
            method=method,
            smooth_sigma=smooth_sigma,
            sigma_clip=sigma_clip,
            maxiters=maxiters,
            bspline_smooth=bspline_smooth,
        )
        self.background_2d = background

        logger.info(
            f"Inter-order background estimated ({method}): min={np.min(background):.1f}, max={np.max(background):.1f}, "
            f"inter-order fraction={np.mean(inter_order):.3f}"
        )
        return background, order_mask

    def apply_background_correction(self, image: np.ndarray,
                                   background: Optional[np.ndarray] = None) -> np.ndarray:
        """
        Apply background correction to image.

        Args:
            image: Original science image
            background: Background model (uses self.background_2d if None)

        Returns:
            Background-corrected image
        """
        bkg = background if background is not None else self.background_2d

        if bkg is None:
            logger.warning("No background model available, returning uncorrected image")
            return image

        if image.shape != bkg.shape:
            raise ValueError(f"Image shape {image.shape} does not match "
                           f"background shape {bkg.shape}")

        corrected = image.astype(float) - bkg.astype(float)
        logger.debug(f"Background corrected: min={np.min(corrected):.1f}, "
                    f"max={np.max(corrected):.1f}")

        return corrected

    def extract_interorder_background(self, image: np.ndarray,
                                     aperture_positions: np.ndarray) -> np.ndarray:
        """
        Extract inter-order background by masking order regions.

        Args:
            image: 2D science image
            aperture_positions: Aperture center positions

        Returns:
            Inter-order background map
        """
        logger.info("Extracting inter-order background...")

        # Create aperture mask
        rows, cols = image.shape
        aperture_mask = np.zeros((rows, cols), dtype=bool)

        # Mark pixels within orders
        for x in range(cols):
            for aperture in aperture_positions:
                center = aperture['center'] + aperture.get('center_slope', 0) * x
                width = aperture.get('width', 20)
                y_min = max(0, int(center - width/2))
                y_max = min(rows, int(center + width/2))
                aperture_mask[y_min:y_max, x] = True

        # Image outside orders is inter-order light
        interorder_pixels = image[~aperture_mask]

        if len(interorder_pixels) > 0:
            median_bkg = np.median(interorder_pixels)
            logger.info(f"Median inter-order background: {median_bkg:.1f}")

            return median_bkg
        else:
            return 0.0

    def save_background(self, output_path: str):
        """Save background model to FITS file."""
        if self.background_2d is None:
            raise RuntimeError("No background model to save")

        from astropy.io import fits
        from pathlib import Path

        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        hdu = fits.PrimaryHDU(data=self.background_2d.astype(np.float32))
        hdu.header['EXTNAME'] = '2D BACKGROUND'
        hdu.header['ORIGIN'] = 'SpecProc'
        hdul = fits.HDUList([hdu])
        hdul.writeto(str(output_path), overwrite=True)

        logger.info(f"Saved background model to {output_path}")


def process_background_stage(config: ConfigManager, science_image: np.ndarray,
                            midpath: str = None,
                            output_subdir: str = 'step2_scatterlight',
                            output_tag: str = '',
                            apertures: Optional[ApertureSet] = None,
                            mask_margin_scale: float = 1.2) -> np.ndarray:
    """
    Execute background subtraction stage.

    Args:
        config: Configuration manager
        science_image: 2D science image
        midpath: Intermediate output directory

    Returns:
        Background model
    """
    # Use default midpath if not provided
    if midpath is None:
        base_output_path = config.get_output_path()
        midpath = str(Path(base_output_path) / 'midpath')

    remover = BackgroundRemover(config)

    # Get polynomial order from config
    poly_order = config.get_int('reduce.background', 'poly_order', 3)
    bg_method = config.get('reduce.background', 'method', 'chebyshev')
    if bg_method == '2d_poly':
        bg_method = 'chebyshev'
    smooth_sigma = config.get_float('reduce.background', 'smooth_sigma', 20.0)
    sigma_clip = config.get_float('reduce.background', 'sigma_clip', 3.0)
    maxiters = config.get_int('reduce.background', 'sigma_clip_maxiters', 4)
    bspline_smooth = config.get_float('reduce.background', 'bspline_smooth', 1.0)
    mask_margin_pixels = config.get_int('reduce.background', 'mask_margin_pixels', 3)
    mask_margin_pixels = max(1, int(round(mask_margin_pixels * float(mask_margin_scale))))

    # Estimate background (inter-order only when apertures are available)
    order_mask = None
    sampled_background = None
    model_lower = None
    model_upper = None
    if apertures is not None and apertures.norders > 0:
        logger.info("Step 3: fitting scattered-light background in split halves (lower/upper CCD)...")
        order_mask = remover._build_aperture_mask(science_image.shape, apertures,
                                                  margin_pixels=mask_margin_pixels)
        sampled_background = science_image.astype(np.float64).copy()
        sampled_background[order_mask] = np.nan

        background, model_lower, model_upper = remover._estimate_background_split_halves(
            science_image,
            order_mask,
            order=poly_order,
            method=bg_method,
            smooth_sigma=smooth_sigma,
            sigma_clip=sigma_clip,
            maxiters=maxiters,
            bspline_smooth=bspline_smooth,
        )
        remover.background_2d = background

        logger.info(
            f"Split-half background estimated ({bg_method}): "
            f"min={np.nanmin(background):.1f}, max={np.nanmax(background):.1f}, "
            f"inter-order fraction={np.mean(~order_mask):.3f}"
        )
    else:
        logger.warning("No apertures available in Step 3; using split-half background fit without order mask.")
        order_mask = np.zeros_like(science_image, dtype=bool)
        sampled_background = science_image.astype(np.float64).copy()
        background, model_lower, model_upper = remover._estimate_background_split_halves(
            science_image,
            order_mask,
            order=poly_order,
            method=bg_method,
            smooth_sigma=smooth_sigma,
            sigma_clip=sigma_clip,
            maxiters=maxiters,
            bspline_smooth=bspline_smooth,
        )
        remover.background_2d = background

    # Save background
    base_output_path = config.get_output_path()
    tag = f"{output_tag}" if output_tag else ''
    background_file = Path(base_output_path) / output_subdir / f'{tag}background_model.fits'
    background_file.parent.mkdir(parents=True, exist_ok=True)
    remover.save_background(str(background_file))

    if order_mask is not None:
        out_dir = Path(base_output_path) / output_subdir
        mask_file = out_dir / f'{tag}order_mask.fits'
        write_fits_image(str(mask_file), order_mask.astype(np.int16), dtype='int16')

        # 1) 级次mask之后的图
        masked_image = science_image.astype(np.float32).copy()
        masked_image[order_mask] = np.nan
        write_fits_image(str(out_dir / f'{tag}order_masked_image.fits'), masked_image, dtype='float32')

        # 2) 扣除级次后剩下的背景图（仅保留道间像素，级次区域设0）
        interorder_only = np.where(order_mask, 0.0, science_image).astype(np.float32)
        write_fits_image(str(out_dir / f'{tag}background_observed.fits'), interorder_only, dtype='float32')

        # 3) 分上下两门的拟合/平滑背景模型
        if model_lower is not None and model_upper is not None:
            rows, cols = science_image.shape
            split = rows // 2
            bg_lower_full = np.zeros((rows, cols), dtype=np.float32)
            bg_upper_full = np.zeros((rows, cols), dtype=np.float32)
            bg_lower_full[:split, :] = model_lower.astype(np.float32)
            bg_upper_full[split:, :] = model_upper.astype(np.float32)
            write_fits_image(str(out_dir / f'{tag}background_model_lower.fits'), bg_lower_full, dtype='float32')
            write_fits_image(str(out_dir / f'{tag}background_model_upper.fits'), bg_upper_full, dtype='float32')

        # 4) background - background_model 残差图（在道间区域评估）
        if sampled_background is not None:
            residual = sampled_background.astype(np.float64) - background.astype(np.float64)
            residual[order_mask] = np.nan
            write_fits_image(str(out_dir / f'{tag}background_minus_model_residual.fits'),
                             residual.astype(np.float32), dtype='float32')

    # 5) 科学图像减去background_model 的图（最终Step 3结果）
    corrected = (science_image.astype(np.float64) - background.astype(np.float64)).astype(np.float32)
    write_fits_image(
        str(Path(base_output_path) / output_subdir / f'{tag}science_minus_backgroundmodel.fits'),
        corrected,
        dtype='float32'
    )

    # Save diagnostic plot if enabled
    save_plots = config.get_bool('reduce', 'save_plots', True)
    if save_plots:
        out_dir = background_file.parent
        fig_format = config.get('reduce', 'fig_format', 'png')
        plot_file = out_dir / f'{tag}background_residuals.{fig_format}'
        plot_background_residuals(science_image, background, str(plot_file))

        # Additional Step 3 diagnostics requested by user.
        def _save_img(arr: np.ndarray, out_name: str, title: str, cmap: str = 'viridis'):
            plt.figure(figsize=(10, 8))
            show = np.array(arr, dtype=float, copy=True)
            finite = np.isfinite(show)
            if np.any(finite):
                vmin = np.nanpercentile(show[finite], 1)
                vmax = np.nanpercentile(show[finite], 99)
            else:
                vmin, vmax = 0.0, 1.0
            plt.imshow(show, origin='lower', aspect='auto', cmap=cmap, vmin=vmin, vmax=vmax)
            plt.colorbar()
            plt.title(title)
            plt.tight_layout()
            plt.savefig(str(out_dir / f'{tag}{out_name}.{fig_format}'), dpi=150)
            plt.close()

        if order_mask is not None:
            masked_image = science_image.astype(float).copy()
            masked_image[order_mask] = np.nan
            _save_img(masked_image, 'order_masked_image', 'Step 3: Order-masked Image (orders removed)')

            interorder_only = np.where(order_mask, np.nan, science_image.astype(float))
            _save_img(interorder_only, 'background_observed', 'Step 3: Inter-order Observed Background')

        if model_lower is not None and model_upper is not None:
            rows, cols = science_image.shape
            split = rows // 2
            bg_lower_full = np.full((rows, cols), np.nan, dtype=float)
            bg_upper_full = np.full((rows, cols), np.nan, dtype=float)
            bg_lower_full[:split, :] = model_lower
            bg_upper_full[split:, :] = model_upper
            _save_img(bg_lower_full, 'background_model_lower', 'Step 3: Background Model (Lower Half)')
            _save_img(bg_upper_full, 'background_model_upper', 'Step 3: Background Model (Upper Half)')

        corrected_full = science_image.astype(float) - background.astype(float)
        if sampled_background is not None:
            residual_bg = sampled_background.astype(float) - background.astype(float)
            if order_mask is not None:
                residual_bg[order_mask] = np.nan
            _save_img(residual_bg, 'background_minus_model_residual', 'Step 3: Background(observed) - Background Model Residual')

        _save_img(corrected_full, 'science_minus_backgroundmodel', 'Step 3: Science - Background Model')

    return background
