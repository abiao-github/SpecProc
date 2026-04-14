"""
Cosmic ray detection and correction stage.
"""

import numpy as np
import logging
import matplotlib.pyplot as plt
from pathlib import Path
from typing import List, Optional
from src.utils.fits_io import read_fits_image, write_fits_image
from src.config.config_manager import ConfigManager
from src.plotting.spectra_plotter import plot_2d_image_to_file

logger = logging.getLogger(__name__)


def _detect_cosmics(image: np.ndarray, sigma: float = 5.0, window: int = 5,
                    fine_sigma: float = 2.0,
                    line_sigma: float = 1.5,
                    grow_sigma: float = 2.5,
                    maxsize: int = 8) -> np.ndarray:
    """Detect cosmic rays using local noise, coherence rejection, and size filtering.

    Approach
    --------
    1. Subtract a small-window local median to reveal sharp deviations.
    2. Build a per-pixel noise estimate from a large-window background model
       (noise ∝ √signal – Poisson dominated) so the threshold adapts to the
       local brightness of each spectral order instead of using a single
       global sigma that is dominated by the dark inter-order background.
     3. Reject coherent spectral structure along BOTH detector axes:
         - cross-dispersion protection: ±1-row neighbours elevated
         - dispersion-direction protection: ±1-column neighbours elevated
         Real orders/lines are spatially coherent; cosmic rays are compact.
     4. Keep only compact connected components and optionally grow them by a
         lower threshold to capture adjacent CR-hit pixels.

    Parameters
    ----------
    image : 2-D array
        Raw or bias-corrected science image (counts/ADU).
    sigma : float
        Detection threshold in units of per-pixel noise (default 5.0).
    window : int
        Kernel size for the fine-structure local median (default 5).
    fine_sigma : float
        SNR threshold for the cross-dispersion PSF test (default 2.0).
        Increase toward 3 to be more aggressive at protecting real lines;
        decrease toward 1.5 to allow sparser PSFs.
    line_sigma : float
        SNR threshold for along-dispersion line continuity protection.
    grow_sigma : float
        Lower SNR threshold used when growing compact cosmic-ray seeds.
    maxsize : int
        Maximum connected-component size kept as a cosmic-ray cluster.
    """
    from scipy import ndimage
    from scipy.ndimage import median_filter

    # --- (1) Fine-structure deviation ---
    fine_med = median_filter(image, size=window)
    deviation = image - fine_med          # positive = locally brighter

    # --- (2) Per-pixel noise model via large-window background ---
    # Background window spans several order widths so it averages over
    # bright orders and dark inter-order gaps, giving a spatially smooth
    # background level used to estimate the local Poisson noise.
    bg_window = max(window * 6, 31)
    if bg_window % 2 == 0:
        bg_window += 1
    bg = median_filter(image, size=bg_window)

    # noise ~ sqrt(max(bg, 0) + noise_floor^2)
    # noise_floor catches regions where bg ≈ 0 (dark inter-order).
    pos_pixels = image[image > 0]
    noise_floor = float(np.percentile(pos_pixels, 5)) if pos_pixels.size > 0 else 1.0
    noise_floor = max(noise_floor, 1.0)
    noise = np.sqrt(np.maximum(bg, 0.0) + noise_floor ** 2)
    noise = np.maximum(noise, 1.0)

    # --- (3) SNR image and initial candidate mask ---
    snr = deviation / noise
    candidates = snr > sigma

    # --- (4) Directional coherence tests ---
    # Rows = cross-dispersion direction; columns = dispersion direction.
    snr_up = np.empty_like(snr)
    snr_dn = np.empty_like(snr)
    snr_lf = np.empty_like(snr)
    snr_rt = np.empty_like(snr)
    snr_up[1:] = snr[:-1];  snr_up[0] = 0.0
    snr_dn[:-1] = snr[1:];  snr_dn[-1] = 0.0
    snr_lf[:, 1:] = snr[:, :-1]; snr_lf[:, 0] = 0.0
    snr_rt[:, :-1] = snr[:, 1:]; snr_rt[:, -1] = 0.0

    psf_y = (snr_up >= fine_sigma) & (snr_dn >= fine_sigma)
    psf_x = (snr_lf >= line_sigma) & (snr_rt >= line_sigma)

    # Require the candidate to stand out above its neighbours, not just above
    # the local noise. This suppresses bright but smooth spectral structure.
    neighbour_max = np.maximum.reduce([snr_up, snr_dn, snr_lf, snr_rt, np.zeros_like(snr)])
    peak_excess = snr - neighbour_max

    seeds = candidates & (~psf_y) & (~psf_x) & (peak_excess > 1.0)

    # --- (5) Compact-component filtering and limited growth ---
    structure = np.ones((3, 3), dtype=int)
    labels, nlabel = ndimage.label(seeds, structure=structure)
    cosmic_mask = np.zeros_like(seeds, dtype=bool)

    protected = psf_y | psf_x
    growable = (snr > grow_sigma) & (~protected)

    for label_id in range(1, nlabel + 1):
        comp = labels == label_id
        if not np.any(comp):
            continue

        yy, xx = np.where(comp)
        size = yy.size
        span_y = int(yy.max() - yy.min() + 1)
        span_x = int(xx.max() - xx.min() + 1)

        # Real spectral features tend to be extended along dispersion. Keep only
        # small, compact clusters as cosmic-ray candidates.
        if size > maxsize or span_x > 3 or span_y > 3 or (span_x * span_y) > 12:
            continue

        grown = ndimage.binary_propagation(comp, structure=structure, mask=(growable | comp))
        gyy, gxx = np.where(grown)
        gsize = gyy.size
        if gsize > max(2 * maxsize, 12):
            continue

        gspan_y = int(gyy.max() - gyy.min() + 1)
        gspan_x = int(gxx.max() - gxx.min() + 1)
        if gspan_x > 4 or gspan_y > 4:
            continue

        cosmic_mask |= grown

    return cosmic_mask


def _fix_cosmics(image: np.ndarray, cosmic_mask: np.ndarray, window: int = 7) -> np.ndarray:
    """Replace cosmic pixels with local median."""
    from scipy.ndimage import median_filter

    corrected = image.copy()
    local_med = median_filter(image, size=window)
    corrected[cosmic_mask] = local_med[cosmic_mask]

    return corrected


def process_cosmic_stage(config: ConfigManager, image_filenames: List[str],
                        midpath: Optional[str] = None,
                        output_subdir: str = 'step1_basic/cosmic_corrected') -> List[str]:
    """
    Execute cosmic ray correction stage for a list of images.

    Args:
        config: Configuration manager
        image_filenames: List of FITS file paths
        midpath: Intermediate output directory

    Returns:
        List of corrected image file paths
    """
    logger.info("=" * 60)
    logger.info("STEP 1: COSMIC RAY REMOVAL")
    logger.info("=" * 60)

    if midpath is None:
        base_output_path = config.get_output_path()

    cosmic_path = Path(base_output_path) / output_subdir
    cosmic_path.mkdir(parents=True, exist_ok=True)

    sigma = config.get_float('reduce', 'cosmic_sigma', 5.0)
    window = config.get_int('reduce', 'cosmic_window', 5)
    fine_sigma = config.get_float('reduce', 'cosmic_fine_sigma', 2.0)
    line_sigma = config.get_float('reduce', 'cosmic_line_sigma', 1.5)
    grow_sigma = config.get_float('reduce', 'cosmic_grow_sigma', 2.5)
    maxsize = config.get_int('reduce', 'cosmic_maxsize', 8)
    save_plots = config.get_bool('reduce', 'save_plots', True)
    fig_format = config.get('reduce', 'fig_format', 'png')

    corrected_files = []

    for filename in image_filenames:
        try:
            input_path = Path(filename)
            output_file = cosmic_path / input_path.name

            image, header = read_fits_image(str(filename))

            mask = _detect_cosmics(
                image,
                sigma=sigma,
                window=window,
                fine_sigma=fine_sigma,
                line_sigma=line_sigma,
                grow_sigma=grow_sigma,
                maxsize=maxsize,
            )
            corrected = _fix_cosmics(image, mask, window=window)

            # Add processing marker to header
            header['COSMICR'] = (True, 'Cosmic ray removal completed')

            write_fits_image(str(output_file), corrected, header=header, overwrite=True, dtype='float32')

            if save_plots:
                plot_file = cosmic_path / f"{input_path.stem}_cosmic_corrected.{fig_format}"
                plot_2d_image_to_file(corrected, str(plot_file), "Cosmic Corrected Image")

                # Mark removed cosmic pixels on top of the original image.
                marked_file = cosmic_path / f"{input_path.stem}_cosmic_marked.{fig_format}"
                fig, ax = plt.subplots(figsize=(10, 8))
                vmin = np.percentile(image, 1)
                vmax = np.percentile(image, 99)
                ax.imshow(image, cmap='viridis', origin='lower', vmin=vmin, vmax=vmax)

                ys, xs = np.where(mask)
                if xs.size > 0:
                    ax.scatter(
                        xs,
                        ys,
                        s=10,
                        facecolors='none',
                        edgecolors='red',
                        linewidths=0.6,
                        alpha=0.9,
                        label=f'Cosmic pixels: {xs.size}'
                    )
                    ax.legend(loc='upper right')

                ax.set_title('Detected Cosmic Rays (marked in red)')
                ax.set_xlabel('Column')
                ax.set_ylabel('Row')
                fig.tight_layout()
                fig.savefig(str(marked_file), dpi=150)
                plt.close(fig)

            corrected_files.append(str(output_file))
            logger.info(f"Cosmic correction: {filename} -> {output_file} ({np.sum(mask)} pixels fixed)")

        except Exception as e:
            logger.error(f"Error cosmic-correcting {filename}: {e}")
            corrected_files.append(str(filename))

    logger.info(f"✓ Cosmic correction complete: {len(corrected_files)} images processed")
    return corrected_files
