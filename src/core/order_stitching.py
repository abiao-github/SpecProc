"""
Order stitching stage.

Merge overlapping neighboring echelle orders into one continuous 1D spectrum.
"""

import logging
from pathlib import Path
from typing import Optional, Tuple

import numpy as np

from src.config.config_manager import ConfigManager
from src.core.data_structures import SpectraSet
from src.plotting.spectra_plotter import plot_spectrum_to_file

logger = logging.getLogger(__name__)


def _clean_order_arrays(wavelength: np.ndarray, flux: np.ndarray,
                        error: Optional[np.ndarray] = None) -> Tuple[np.ndarray, np.ndarray, Optional[np.ndarray]]:
    """Clean and sort one order by wavelength."""
    good = np.isfinite(wavelength) & np.isfinite(flux)
    if error is not None:
        good &= np.isfinite(error)

    wl = wavelength[good]
    fl = flux[good]
    er = error[good] if error is not None else None

    if wl.size == 0:
        return wl, fl, er

    idx = np.argsort(wl)
    wl = wl[idx]
    fl = fl[idx]
    if er is not None:
        er = er[idx]
    return wl, fl, er


def stitch_orders(spectra_set: SpectraSet,
                 min_overlap_points: int = 20) -> Tuple[np.ndarray, np.ndarray, Optional[np.ndarray]]:
    """
    Stitch all orders into one continuous spectrum.

    The algorithm scales neighboring orders in overlapping wavelength regions,
    then merges overlap by weighted averaging.
    """
    if spectra_set is None or len(spectra_set.spectra) == 0:
        raise RuntimeError("No spectra available for order stitching")

    items = []
    for spec in spectra_set.spectra.values():
        wl, fl, er = _clean_order_arrays(spec.wavelength, spec.flux, spec.error)
        if wl.size > 0:
            items.append((np.nanmedian(wl), wl, fl, er, spec.aperture))

    if not items:
        raise RuntimeError("No valid wavelength/flux points for stitching")

    items.sort(key=lambda t: t[0])

    stitched_wl = items[0][1].copy()
    stitched_fl = items[0][2].copy()
    stitched_er = items[0][3].copy() if items[0][3] is not None else None

    for _, next_wl, next_fl, next_er, aperture in items[1:]:
        w_lo = max(stitched_wl.min(), next_wl.min())
        w_hi = min(stitched_wl.max(), next_wl.max())

        if w_hi > w_lo:
            prev_overlap = (stitched_wl >= w_lo) & (stitched_wl <= w_hi)
            next_overlap = (next_wl >= w_lo) & (next_wl <= w_hi)

            if np.sum(prev_overlap) >= min_overlap_points and np.sum(next_overlap) >= min_overlap_points:
                grid_n = int(min(np.sum(prev_overlap), np.sum(next_overlap)))
                grid = np.linspace(w_lo, w_hi, max(grid_n, min_overlap_points))

                prev_interp = np.interp(grid, stitched_wl, stitched_fl)
                next_interp = np.interp(grid, next_wl, next_fl)

                good = np.isfinite(prev_interp) & np.isfinite(next_interp) & (np.abs(next_interp) > 1e-12)
                if np.any(good):
                    scale = np.nanmedian(prev_interp[good] / next_interp[good])
                    if np.isfinite(scale) and scale > 0:
                        next_fl = next_fl * scale
                        if next_er is not None:
                            next_er = next_er * scale
                        logger.info(f"Order stitching: scaled aperture {aperture} by factor {scale:.5f}")

                prev_interp = np.interp(grid, stitched_wl, stitched_fl)
                next_interp = np.interp(grid, next_wl, next_fl)

                if stitched_er is not None and next_er is not None:
                    prev_ei = np.interp(grid, stitched_wl, stitched_er)
                    next_ei = np.interp(grid, next_wl, next_er)
                    w1 = 1.0 / np.clip(prev_ei, 1e-12, None) ** 2
                    w2 = 1.0 / np.clip(next_ei, 1e-12, None) ** 2
                    overlap_flux = (prev_interp * w1 + next_interp * w2) / (w1 + w2 + 1e-12)
                    overlap_err = np.sqrt(1.0 / (w1 + w2 + 1e-12))
                else:
                    overlap_flux = 0.5 * (prev_interp + next_interp)
                    overlap_err = None

                prev_keep = stitched_wl < w_lo
                next_keep = next_wl > w_hi

                new_wl = np.concatenate([stitched_wl[prev_keep], grid, next_wl[next_keep]])
                new_fl = np.concatenate([stitched_fl[prev_keep], overlap_flux, next_fl[next_keep]])

                if stitched_er is not None and next_er is not None and overlap_err is not None:
                    new_er = np.concatenate([stitched_er[prev_keep], overlap_err, next_er[next_keep]])
                else:
                    new_er = None

                stitched_wl, stitched_fl, stitched_er = new_wl, new_fl, new_er
                continue

        stitched_wl = np.concatenate([stitched_wl, next_wl])
        stitched_fl = np.concatenate([stitched_fl, next_fl])
        if stitched_er is not None and next_er is not None:
            stitched_er = np.concatenate([stitched_er, next_er])
        else:
            stitched_er = None

    idx = np.argsort(stitched_wl)
    stitched_wl = stitched_wl[idx]
    stitched_fl = stitched_fl[idx]
    if stitched_er is not None:
        stitched_er = stitched_er[idx]

    return stitched_wl, stitched_fl, stitched_er


def save_stitched_spectrum(output_path: str, wavelength: np.ndarray,
                          flux: np.ndarray, error: Optional[np.ndarray] = None):
    """Save stitched 1D spectrum as FITS table."""
    from astropy.io import fits

    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    cols = [
        fits.Column(name='WAVELENGTH', format='D', array=wavelength),
        fits.Column(name='FLUX', format='D', array=flux),
    ]
    if error is not None:
        cols.append(fits.Column(name='ERROR', format='D', array=error))

    table = fits.BinTableHDU.from_columns(cols)
    pri = fits.PrimaryHDU()
    pri.header['NSPEC'] = (len(wavelength), 'Number of stitched pixels')
    pri.header['CONTENT'] = ('Stitched 1D spectrum', 'Order stitching output')

    fits.HDUList([pri, table]).writeto(str(out), overwrite=True)
    logger.info(f"Saved stitched spectrum to {out}")


def process_order_stitching_stage(spectra_set: SpectraSet,
                                 output_dir_base: str,
                                 output_subdir: str = 'step8_stitching',
                                 save_plots: bool = True,
                                 fig_format: str = 'png') -> Tuple[np.ndarray, np.ndarray, Optional[np.ndarray]]:
    """Run order stitching stage and save outputs."""
    wl, fl, er = stitch_orders(spectra_set)

    out_dir = Path(output_dir_base) / output_subdir
    out_dir.mkdir(parents=True, exist_ok=True)

    save_stitched_spectrum(str(out_dir / 'stitched_spectrum.fits'), wl, fl, er)

    if save_plots:
        plot_spectrum_to_file(
            wl,
            fl,
            str(out_dir / f'stitched_spectrum.{fig_format}'),
            er,
            'Stitched 1D Spectrum'
        )

    return wl, fl, er
