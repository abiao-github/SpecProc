"""
2D flat-field correction stage — Step 4 of the spectral reduction pipeline.

Applies the pixel-to-pixel flat correction map (built in Step 2 / order_tracing)
to a science image, and persists all flat-model artefacts to the step4 output
directory.

Main entry point: process_flat_correction_stage()
"""

import logging
import pickle
from pathlib import Path
from typing import Optional

import numpy as np



from scipy import ndimage
from scipy.signal import savgol_filter
from scipy.interpolate import UnivariateSpline
from scipy.ndimage import gaussian_filter1d
from typing import Dict, Tuple, List
from src.core.data_structures import FlatField, ApertureSet
from src.utils.fits_io import write_fits_image
from src.plotting.spectra_plotter import plot_2d_image_to_file

logger = logging.getLogger(__name__)


def save_flat_correction_products(flat_field: FlatField,
                                  apertures: Optional[ApertureSet],
                                  out_dir: Path,
                                  processor=None,
                                  fig_format: str = 'png',
                                  save_plots: bool = True) -> None:
    """
    Persist Step-4 flat-model artefacts (blaze profiles, sensitivity map,
    pixel flat, 2D correction map, response map, diagnostic plots).

    Parameters
    ----------
    processor : FlatCorrectionModelBuilder, optional
        If provided, re-use its already-computed response_map and diagnostics
        instead of rebuilding from scratch (saves ~5 s).
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    if flat_field.blaze_profiles:
        with open(out_dir / 'blaze_profiles.pkl', 'wb') as f:
            pickle.dump(flat_field.blaze_profiles, f)

    if flat_field.flat_sens is not None:
        write_fits_image(
            str(out_dir / 'flat_sensitivity.fits'),
            flat_field.flat_sens.astype(np.float32),
            dtype='float32',
        )

    if flat_field.smoothed_model is not None:
        write_fits_image(
            str(out_dir / 'flat_smoothed_model.fits'),
            flat_field.smoothed_model.astype(np.float32),
            dtype='float32',
        )

    if flat_field.pixel_flat is not None:
        write_fits_image(
            str(out_dir / 'pixel_to_pixel_flat.fits'),
            flat_field.pixel_flat.astype(np.float32),
            dtype='float32',
        )

    if flat_field.flat_corr_2d is not None:
        write_fits_image(
            str(out_dir / 'flat_correction_2d.fits'),
            flat_field.flat_corr_2d.astype(np.float32),
            dtype='float32',
        )

    # If processor is provided, save its diagnostics directly.
    if processor is not None:
        if processor.response_map is not None:
            write_fits_image(
                str(out_dir / 'flat_response_map.fits'),
                processor.response_map.astype(np.float32),
                dtype='float32',
            )

        # Compute clean_flat for diagnostics
        clean_flat = flat_field.flat_data
        if flat_field.scattered_light is not None:
            clean_flat = np.clip(
                flat_field.flat_data.astype(np.float32) - flat_field.scattered_light.astype(np.float32),
                1e-6,
                None,
            )

        clean_flat_corrected = None
        if flat_field.flat_corr_2d is not None:
            safe_flat = flat_field.flat_corr_2d.astype(np.float32)
            bad = (~np.isfinite(safe_flat)) | (safe_flat <= 0.05)
            if np.any(bad):
                safe_flat = safe_flat.copy()
                safe_flat[bad] = 1.0
            clean_flat_corrected = clean_flat.astype(np.float32) / safe_flat

        processor.save_step4_diagnostics(
            out_dir,
            clean_flat=clean_flat.astype(np.float32),
            clean_flat_corrected=clean_flat_corrected,
        )



class FlatCorrectionModelBuilder:
    def __init__(self):
        self.flat_data = None
        self.flat_mask = None
        self.response_map = None
        self.smoothed_model = None
        self.pixel_flat = None
        self.illumination_flat = None
        self.flat_corr_2d = None
        self.flat_sens = None
        self.order_diagnostics = {}

    def _smooth_order_response(self, flux_1d: np.ndarray, method: str = 'savgol',
                               window: int = 21, bspline_smooth: float = 0.5) -> np.ndarray:
        """Smooth per-order 1D blaze extraction with selectable method."""
        x = np.arange(flux_1d.size)
        y = np.asarray(flux_1d, dtype=float)

        valid = np.isfinite(y) & (y > 0)
        if np.sum(valid) < max(20, int(0.2 * y.size)):
            med = np.nanmedian(y[valid]) if np.any(valid) else 1.0
            return np.full_like(y, med if np.isfinite(med) and med > 0 else 1.0)

        y_interp = np.interp(x, x[valid], y[valid])
        win = max(5, int(window))
        win = min(win, y.size - 1 if y.size % 2 == 0 else y.size)
        if win < 7:
            win = 7
        if win % 2 == 0:
            win += 1

        m = (method or 'savgol').strip().lower()
        if m == 'median':
            smooth = ndimage.median_filter(y_interp, size=win, mode='nearest')
            smooth = savgol_filter(smooth, window_length=win, polyorder=min(3, win - 2), mode='interp')
        elif m == 'bspline':
            k = min(3, max(1, len(x) - 1))
            s = max(0.0, float(bspline_smooth)) * len(x)
            smooth = UnivariateSpline(x, y_interp, k=k, s=s)(x)
        else:
            smooth = savgol_filter(y_interp, window_length=win, polyorder=min(3, win - 2), mode='interp')

        # Iterative clipping similar to apnormalize-style robust fitting.
        fit = smooth.copy()
        keep = np.ones_like(y_interp, dtype=bool)
        for _ in range(3):
            resid = y_interp - fit
            sigma = np.std(resid[keep]) if np.any(keep) else 0.0
            if sigma <= 0:
                break
            new_keep = np.abs(resid) < 3.0 * sigma
            if np.array_equal(new_keep, keep):
                break
            keep = new_keep
            if np.sum(keep) < 10:
                break
            fit = np.interp(x, x[keep], y_interp[keep])
            fit = savgol_filter(fit, window_length=win, polyorder=3, mode='interp')

        fit = np.where(fit > 0, fit, np.nanmedian(fit[fit > 0]) if np.any(fit > 0) else 1.0)
        return fit

    def _smooth_blaze_sparse_bspline(self, flux_1d: np.ndarray,
                                      knot_spacing: int = 500,
                                      edge_nknots: int = 6) -> np.ndarray:
        """Sparse-knots B-spline blaze extraction (PypeIt-style).

        Interior knots are placed at wide spacing (e.g. 500 px) so the
        spline is rigid enough to crush fringe ripples, while extra knots
        are packed near both edges to faithfully track the steep blaze
        roll-off.

        Parameters
        ----------
        flux_1d : array
            Raw 1-D aperture-summed blaze.
        knot_spacing : int
            Spacing between interior knots in the flat central region (px).
        edge_nknots : int
            Number of extra knots placed in each edge transition zone.
        """
        from scipy.interpolate import splrep, splev

        x = np.arange(flux_1d.size)
        y = np.asarray(flux_1d, dtype=float)

        valid = np.isfinite(y) & (y > 0)
        n_valid = int(np.sum(valid))
        if n_valid < max(20, int(0.2 * y.size)):
            med = np.nanmedian(y[valid]) if np.any(valid) else 1.0
            return np.full_like(y, med if np.isfinite(med) and med > 0 else 1.0)

        # Light pre-clean: narrow median kills single-pixel outliers / cosmic rays
        y_clean = ndimage.median_filter(
            np.interp(x, x[valid], y[valid]), size=5, mode='nearest')

        # Determine the valid (signal) region: > 5% of peak
        peak = np.max(y_clean)
        sig_mask = y_clean > (peak * 0.05)
        sig_idx = np.where(sig_mask)[0]
        if sig_idx.size < 20:
            return y_clean.copy()
        start, end = int(sig_idx[0]), int(sig_idx[-1])

        # --- Build knot vector ---
        # Edge transition zone ≈ 10% of order width, at least 50 px
        span = end - start
        edge_zone = max(50, int(span * 0.10))

        # Interior: sparse knots in the flat central region
        interior_start = start + edge_zone
        interior_end = end - edge_zone
        spacing = max(50, int(knot_spacing))

        if interior_end > interior_start + spacing:
            mid_knots = np.arange(interior_start, interior_end, spacing, dtype=float)
        else:
            # Very short order: just a few uniform knots
            mid_knots = np.linspace(interior_start, interior_end,
                                    max(3, span // 100), dtype=float)

        # Edges: dense knots (logarithmic-ish spacing toward the boundary)
        n_ek = max(2, int(edge_nknots))
        left_offsets = np.geomspace(10, edge_zone, n_ek)
        right_offsets = np.geomspace(10, edge_zone, n_ek)
        left_knots = start + left_offsets
        right_knots = end - right_offsets

        all_knots = np.sort(np.unique(np.concatenate([
            left_knots, mid_knots, right_knots
        ])))

        # Ensure knots are strictly within data range (required by splrep)
        x_fit = x[sig_mask].astype(float)
        y_fit = y_clean[sig_mask]
        all_knots = all_knots[(all_knots > x_fit[0]) & (all_knots < x_fit[-1])]
        if all_knots.size < 2:
            # Fallback: uniform sparse knots
            all_knots = np.linspace(x_fit[0] + 10, x_fit[-1] - 10,
                                    max(4, span // spacing))

        # Cubic B-spline fit with prescribed knots
        try:
            tck = splrep(x_fit, y_fit, t=all_knots, k=3, task=-1)
            blaze_smooth = splev(x, tck)
        except Exception:
            # Fallback: smoothing spline if custom knots fail
            tck = splrep(x_fit, y_fit, k=3, s=n_valid)
            blaze_smooth = splev(x, tck)

        # Outside valid region: copy raw values
        blaze_smooth[~sig_mask] = y_clean[~sig_mask]

        # Safety: no negative or zero values
        blaze_smooth = np.where(
            blaze_smooth > 0, blaze_smooth,
            np.nanmedian(blaze_smooth[blaze_smooth > 0])
            if np.any(blaze_smooth > 0) else 1.0)
        return blaze_smooth

    def build_order_response_map(self, apertures: ApertureSet,
                                 source_image: Optional[np.ndarray] = None,
                                 blaze_knot_spacing: int = 500,
                                 blaze_edge_nknots: int = 6,
                                 width_smooth_window: int = 41,
                                 profile_bin_step: float = 0.01,
                                 n_profile_segments: int = 100,
                                 profile_smooth_sigma: float = 6.0,
                                 pixel_flat_min: float = 0.5,
                                 pixel_flat_max: float = 1.5
                                 ) -> Tuple[np.ndarray, Dict[int, np.ndarray], Dict[int, np.ndarray], np.ndarray, np.ndarray, Optional[np.ndarray]]:
        """
        Build per-order response model and convert it to a 2D sensitivity map.

        This follows the echelle flat-field philosophy used in IRAF/PypeIt/REDUCE:
        estimate order-by-order blaze/response first, then construct a 2D correction map.

        Improvements over a single global profile:
        1. **Segmented cross-dispersion profiles**: The profile is built in
           ``n_profile_segments`` bins along X, then linearly interpolated at
           each column, so the model tracks PSF shape changes (aberrations,
           defocus, vignetting).
        2. **Sparse-knots B-spline blaze fitting**: Interior knots are placed
           at wide spacing (``blaze_knot_spacing``) to crush fringe ripples,
           while dense knots near edges track the steep blaze roll-off
           — so the 2D model is fringe-free and the resulting pixel_flat
           absorbs fringe structure for correction.
        """
        if self.flat_data is None and source_image is None:
            raise RuntimeError("No flat data")

        image = source_image if source_image is not None else self.flat_data
        image = np.asarray(image, dtype=np.float64)
        height, width = image.shape
        x_coords = np.arange(width)

        # Clamp segments to a sane range
        n_profile_segments = max(1, min(n_profile_segments, 200))

        response_map = np.ones_like(image, dtype=np.float32)
        coverage = np.zeros_like(image, dtype=bool)

        blaze_profiles: Dict[int, np.ndarray] = {}
        cross_profiles: Dict[int, np.ndarray] = {}
        order_diag: Dict[int, Dict[str, np.ndarray]] = {}

        for aperture_id, aperture in apertures.apertures.items():
            center = aperture.get_position(x_coords)
            lower = aperture.get_lower(x_coords)
            upper = aperture.get_upper(x_coords)

            # 1D aperture sum extraction inside traced boundaries.
            raw_blaze = np.full(width, np.nan, dtype=float)
            width_raw = np.full(width, np.nan, dtype=float)
            per_col_bounds = []

            for x in x_coords:
                y_lo = min(lower[x], upper[x])
                y_hi = max(lower[x], upper[x])
                y1 = max(0, int(np.floor(y_lo)))
                y2 = min(height, int(np.ceil(y_hi + 1)))
                if y2 - y1 < 3:
                    per_col_bounds.append(None)
                    continue

                stripe = image[y1:y2, x]
                yy = np.arange(y1, y2)
                if self.flat_mask is not None:
                    local_bad = self.flat_mask[y1:y2, x] > 0
                    if np.all(local_bad):
                        per_col_bounds.append(None)
                        continue
                    stripe = stripe[~local_bad]
                    yy = yy[~local_bad]

                if stripe.size < 3:
                    per_col_bounds.append(None)
                    continue

                raw_blaze[x] = np.sum(stripe)
                width_raw[x] = max(1.0, y_hi - y_lo)
                per_col_bounds.append((yy, stripe))

            # --- Smooth blaze (sparse-knots B-spline: rigid interior, flexible edges) ---
            blaze_1d = self._smooth_blaze_sparse_bspline(
                raw_blaze,
                knot_spacing=blaze_knot_spacing,
                edge_nknots=blaze_edge_nknots,
            )

            width_smooth = self._smooth_order_response(
                width_raw,
                method='savgol',
                window=width_smooth_window,
                bspline_smooth=0.0,
            )
            width_smooth = np.clip(width_smooth, 1.0, None)

            norm = np.nanmedian(blaze_1d[np.isfinite(blaze_1d) & (blaze_1d > 0)])
            if not np.isfinite(norm) or norm <= 0:
                norm = 1.0
            blaze = blaze_1d / norm
            blaze = np.where(np.isfinite(blaze) & (blaze > 1e-6), blaze, 1.0)
            blaze_profiles[aperture_id] = blaze.astype(np.float32)

            valid_cols = np.where(np.isfinite(raw_blaze))[0]
            if valid_cols.size == 0:
                logger.warning(f"Order {aperture_id}: no valid columns — skipped")
                continue

            sample_cols = set()
            if valid_cols.size > 0:
                sample_cols = set(np.unique(np.linspace(valid_cols.min(), valid_cols.max(), num=min(6, valid_cols.size), dtype=int)))

            # ==============================================================
            # Build segmented cross-dispersion profiles along X
            # ==============================================================
            step = max(0.002, float(profile_bin_step))
            n_seg = n_profile_segments

            # Determine column ranges for each segment
            x_lo_seg = valid_cols.min()
            x_hi_seg = valid_cols.max()
            seg_edges = np.linspace(x_lo_seg, x_hi_seg + 1, n_seg + 1).astype(int)
            seg_centers = 0.5 * (seg_edges[:-1] + seg_edges[1:])  # X center of each segment

            # Build one profile per segment
            seg_profiles_x = []   # list of (prof_x_i,) arrays — normalized y coordinates
            seg_profiles_y = []   # list of (prof_y_i,) arrays — profile values
            common_yn_grid = None  # will be unified after first pass

            for si in range(n_seg):
                xs0, xs1 = int(seg_edges[si]), int(seg_edges[si + 1])
                yn_seg = []
                fn_seg = []
                for x in range(xs0, min(xs1, width)):
                    if per_col_bounds[x] is None:
                        continue
                    yy, stripe = per_col_bounds[x]
                    bval = blaze_1d[x]
                    wval = width_smooth[x]
                    if (not np.isfinite(bval)) or bval <= 0 or (not np.isfinite(wval)) or wval <= 0:
                        continue
                    yn = (yy.astype(float) - center[x]) / wval
                    fn = stripe / bval
                    good = np.isfinite(yn) & np.isfinite(fn)
                    if np.any(good):
                        yn_seg.append(yn[good])
                        fn_seg.append(fn[good])

                if not yn_seg:
                    seg_profiles_x.append(None)
                    seg_profiles_y.append(None)
                    continue

                yn_cat = np.concatenate(yn_seg)
                fn_cat = np.concatenate(fn_seg)

                yn_min = max(-1.2, float(np.nanpercentile(yn_cat, 1)))
                yn_max_v = min(1.2, float(np.nanpercentile(yn_cat, 99)))
                if yn_max_v - yn_min < 0.2:
                    yn_min, yn_max_v = -0.8, 0.8
                bins = np.arange(yn_min, yn_max_v + step, step)
                if bins.size < 8:
                    bins = np.linspace(-0.8, 0.8, 81)

                idx = np.digitize(yn_cat, bins) - 1
                px, py = [], []
                for i in range(len(bins) - 1):
                    pick = idx == i
                    if np.sum(pick) < 6:
                        continue
                    px.append(0.5 * (bins[i] + bins[i + 1]))
                    py.append(np.nanmedian(fn_cat[pick]))

                if len(px) < 8:
                    seg_profiles_x.append(None)
                    seg_profiles_y.append(None)
                    continue

                px = np.asarray(px, dtype=float)
                py = np.asarray(py, dtype=float)
                py = np.where(np.isfinite(py), py, 1.0)
                py = np.clip(py, 1e-6, None)

                # Smooth the binned profile
                py = self._smooth_order_response(
                    py,
                    method='savgol',
                    window=max(11, int(11 * max(1.0, step / 0.01))),
                    bspline_smooth=0.0,
                )

                seg_profiles_x.append(px)
                seg_profiles_y.append(py)

            # Count valid segments
            valid_seg_idx = [i for i in range(n_seg) if seg_profiles_x[i] is not None]

            if not valid_seg_idx:
                logger.warning(f"Order {aperture_id}: no valid profile segments — skipped")
                continue

            # Build a common y_norm grid that covers all segments
            all_yn_min = min(seg_profiles_x[i].min() for i in valid_seg_idx)
            all_yn_max = max(seg_profiles_x[i].max() for i in valid_seg_idx)
            common_yn = np.arange(all_yn_min, all_yn_max + step, step)

            # Resample each segment profile onto the common grid
            seg_on_grid = np.full((n_seg, common_yn.size), np.nan, dtype=float)
            for i in valid_seg_idx:
                seg_on_grid[i] = np.interp(common_yn, seg_profiles_x[i], seg_profiles_y[i],
                                           left=np.nan, right=np.nan)

            # Fill missing segments by nearest-neighbour copy
            for i in range(n_seg):
                if i in valid_seg_idx:
                    continue
                # Find nearest valid segment
                dists = [abs(i - vi) for vi in valid_seg_idx]
                nearest = valid_seg_idx[int(np.argmin(dists))]
                seg_on_grid[i] = seg_on_grid[nearest]

            # ----------------------------------------------------------
            # Smooth profiles along dispersion direction (axis 0 = segment index)
            # so that CSP varies continuously from segment to segment.
            # ----------------------------------------------------------
            if n_seg >= 3 and profile_smooth_sigma > 0:
                for j in range(seg_on_grid.shape[1]):
                    col = seg_on_grid[:, j]
                    ok = np.isfinite(col)
                    if np.sum(ok) < 3:
                        continue
                    # Interpolate NaN gaps before smoothing
                    idx_ok = np.where(ok)[0]
                    col_filled = np.interp(np.arange(n_seg), idx_ok, col[idx_ok])
                    seg_on_grid[:, j] = gaussian_filter1d(
                        col_filled, sigma=profile_smooth_sigma, mode='nearest')

            # Store the composite profile for diagnostics (average of all segments)
            with np.errstate(all='ignore'):
                mean_prof_y = np.nanmean(seg_on_grid, axis=0)
            mean_prof_y = np.where(np.isfinite(mean_prof_y), mean_prof_y, 1e-6)
            csum = np.nansum(mean_prof_y)
            if np.isfinite(csum) and csum > 0:
                cross_profiles[aperture_id] = (mean_prof_y / csum).astype(np.float32)

            sampled_realspace = []

            # ==============================================================
            # Reconstruct 2D model: Model(X,Y) = P(y_norm; X) × Blaze(X)
            # with P interpolated between segment profiles
            # ==============================================================
            for x in x_coords:
                y_lo = min(lower[x], upper[x])
                y_hi = max(lower[x], upper[x])
                y1 = max(0, int(np.floor(y_lo)))
                y2 = min(height, int(np.ceil(y_hi + 1)))
                if y2 - y1 < 3:
                    continue

                bval = blaze_1d[x]
                wval = width_smooth[x]
                if (not np.isfinite(bval)) or bval <= 0 or (not np.isfinite(wval)) or wval <= 0:
                    continue

                yy = np.arange(y1, y2, dtype=float)
                yreq = (yy - center[x]) / wval

                # Interpolate profile between the two nearest segments
                if n_seg == 1:
                    pval = np.interp(yreq, common_yn, seg_on_grid[0], left=1.0, right=1.0)
                else:
                    # Find the fractional segment index for this column
                    frac = np.clip((float(x) - seg_centers[0]) / max(1.0, seg_centers[-1] - seg_centers[0]),
                                   0.0, 1.0)
                    seg_f = frac * (n_seg - 1)
                    si_lo = max(0, min(int(np.floor(seg_f)), n_seg - 2))
                    si_hi = si_lo + 1
                    w_hi = seg_f - si_lo
                    w_lo = 1.0 - w_hi

                    p_lo = np.interp(yreq, common_yn, seg_on_grid[si_lo], left=1.0, right=1.0)
                    p_hi = np.interp(yreq, common_yn, seg_on_grid[si_hi], left=1.0, right=1.0)
                    pval = w_lo * p_lo + w_hi * p_hi

                model_col = pval * bval

                if x in sample_cols:
                    sampled_realspace.append({
                        'x': int(x),
                        'y': yy.astype(np.float32),
                        'ynorm': yreq.astype(np.float32),
                        'profile_val': pval.astype(np.float32),
                        'model_flux': model_col.astype(np.float32),
                    })

                response_map[y1:y2, x] = model_col.astype(np.float32)
                coverage[y1:y2, x] = True

            order_diag[aperture_id] = {
                'raw_blaze': raw_blaze.astype(np.float32),
                'blaze_fit': blaze_1d.astype(np.float32),
                'width_smooth': width_smooth.astype(np.float32),
                'profile_x': common_yn.astype(np.float32),
                'profile_y': mean_prof_y.astype(np.float32),
                'seg_on_grid': seg_on_grid.astype(np.float32),
                'seg_centers': seg_centers.astype(np.float32),
                'sampled_realspace': sampled_realspace,
                'n_profile_segments': n_seg,
                'n_valid_segments': len(valid_seg_idx),
            }

        # Smoothed model keeps blaze + cross-order envelope, without pixel noise.
        smoothed_model = response_map.astype(np.float32)

        # Pixel-to-pixel flat: high-frequency component around unity.

        pixel_flat = np.ones_like(image, dtype=np.float32)
        pixel_flat[coverage] = image[coverage] / (smoothed_model[coverage] + 1e-10)
        bad_cov = coverage & (~np.isfinite(pixel_flat) | (pixel_flat < pixel_flat_min) | (pixel_flat > pixel_flat_max))
        if np.any(bad_cov):
            pixel_flat[bad_cov] = 1.0
        if np.any(coverage):
            pf_med = np.nanmedian(pixel_flat[coverage][np.isfinite(pixel_flat[coverage])])
            if np.isfinite(pf_med) and pf_med > 0:
                pixel_flat[coverage] = pixel_flat[coverage] / pf_med
        pixel_flat = np.where(np.isfinite(pixel_flat), pixel_flat, 1.0)
        pixel_flat[~coverage] = 1.0
        # 进一步消除级次边界突变：将 coverage 区域内 <0.7 或 >1.3 的像元直接赋值为 1
        edge_mask = coverage & ((pixel_flat < 0.7) | (pixel_flat > 1.3))
        if np.any(edge_mask):
            pixel_flat[edge_mask] = 1.0
        pixel_flat = np.clip(pixel_flat, 0.1, 10.0).astype(np.float32)

        # Illumination flat is intentionally not used/output in pixel-to-pixel-only workflow.
        illumination_flat = None

        # Final 2D correction map equals pixel flat in this workflow.
        flat_corr_2d = np.clip(pixel_flat, 0.1, 10.0).astype(np.float32)

        sens = flat_corr_2d.copy()

        # Normalize sensitivity map to unity on valid pixels.
        if np.any(coverage):
            med = np.median(sens[coverage][np.isfinite(sens[coverage])])
            if np.isfinite(med) and med > 0:
                sens[coverage] /= med

        sens = np.clip(sens, 0.25, 4.0).astype(np.float32)

        self.response_map = response_map
        self.smoothed_model = smoothed_model
        self.pixel_flat = pixel_flat
        self.illumination_flat = illumination_flat
        self.flat_corr_2d = flat_corr_2d
        self.flat_sens = sens
        self.order_diagnostics = order_diag
        logger.info("Built Step4 flat model via 1D blaze + normalized-profile reconstruction")
        return sens, blaze_profiles, cross_profiles, smoothed_model, pixel_flat, illumination_flat

    def save_step4_diagnostics(self, output_dir: Path,
                               clean_flat: Optional[np.ndarray] = None,
                               clean_flat_corrected: Optional[np.ndarray] = None,
                               fig_format: str = 'png',
                               save_plots: bool = True):
        """Save step4 diagnostic artifacts: blaze/profile/inverse-map/model/pixel-flat and corrected flat."""
        out_dir = Path(output_dir)
        out_dir.mkdir(parents=True, exist_ok=True)

        # fig_format and save_plots are now passed via method parameters
        # (kept as instance defaults for backward compat with save_step4_diagnostics signature)

        # ---- FITS outputs ----
        if self.response_map is not None:
            write_fits_image(str(out_dir / 'model_flat_2d.fits'), self.response_map.astype(np.float32), dtype='float32')
        if self.pixel_flat is not None:
            write_fits_image(str(out_dir / 'pixel_flat_2d.fits'), self.pixel_flat.astype(np.float32), dtype='float32')
        if self.flat_sens is not None:
            write_fits_image(str(out_dir / 'sensitivity_map_2d.fits'), self.flat_sens.astype(np.float32), dtype='float32')
        if clean_flat is not None:
            write_fits_image(str(out_dir / 'clean_master_flat.fits'), clean_flat.astype(np.float32), dtype='float32')
        if clean_flat_corrected is not None:
            write_fits_image(str(out_dir / 'clean_master_flat_flat2d_corrected.fits'), clean_flat_corrected.astype(np.float32), dtype='float32')

        # ---- 2D overview PNG plots ----
        if save_plots:
            if self.response_map is not None:
                plot_2d_image_to_file(self.response_map, str(out_dir / f'model_flat_2d.{fig_format}'), 'Reconstructed 2D Flat Model')
            if self.pixel_flat is not None:
                plot_2d_image_to_file(self.pixel_flat, str(out_dir / f'pixel_flat_2d.{fig_format}'),
                                      'Pixel-to-Pixel Flat', vmin=0.85, vmax=1.15)
            if clean_flat is not None:
                plot_2d_image_to_file(clean_flat, str(out_dir / f'clean_master_flat.{fig_format}'), 'Clean Master Flat')
            if clean_flat_corrected is not None:
                plot_2d_image_to_file(clean_flat_corrected, str(out_dir / f'clean_master_flat_flat2d_corrected.{fig_format}'),
                                      'Clean Master Flat After 2D Flat Correction')

        if not self.order_diagnostics:
            return

        import matplotlib.pyplot as plt
        from matplotlib.cm import ScalarMappable
        from matplotlib.colors import Normalize

        blaze_dir = out_dir / 'blaze_per_order'
        profile_dir = out_dir / 'profile_per_order'
        realspace_dir = out_dir / 'realspace_profile_per_order'
        blaze_dir.mkdir(parents=True, exist_ok=True)
        profile_dir.mkdir(parents=True, exist_ok=True)
        realspace_dir.mkdir(parents=True, exist_ok=True)

        # Collect data for combined overview plots.
        all_order_ids = sorted(self.order_diagnostics.keys())

        # ---- Per-order diagnostic loop ----
        for aperture_id, diag in self.order_diagnostics.items():
            raw_blaze = diag.get('raw_blaze')
            blaze_fit = diag.get('blaze_fit')
            profile_x = diag.get('profile_x')
            profile_y = diag.get('profile_y')
            seg_on_grid = diag.get('seg_on_grid')
            seg_centers = diag.get('seg_centers')
            sampled_realspace = diag.get('sampled_realspace', [])
            n_seg = diag.get('n_profile_segments', 1)
            n_valid_seg = diag.get('n_valid_segments', 0)

            # Save numerical data
            np.savez(
                blaze_dir / f'order_{aperture_id:03d}_blaze.npz',
                raw_blaze=raw_blaze,
                blaze_fit=blaze_fit,
                width_smooth=diag.get('width_smooth'),
            )
            np.savez(
                profile_dir / f'order_{aperture_id:03d}_profile.npz',
                y_norm=profile_x,
                profile_mean=profile_y,
                seg_on_grid=seg_on_grid,
                seg_centers=seg_centers,
            )

            if not save_plots:
                continue

            # --- Per-order blaze plot ---
            if raw_blaze is not None and blaze_fit is not None:
                plt.figure(figsize=(10, 4))
                x = np.arange(raw_blaze.size)
                plt.plot(x, raw_blaze, color='0.6', linewidth=0.8, label='Raw Blaze 1D')
                plt.plot(x, blaze_fit, color='tab:red', linewidth=1.4, label='B-spline Fit')
                plt.xlabel('X (column)')
                plt.ylabel('Flux')
                plt.title(f'Order {aperture_id} — Blaze (segments={n_seg}, valid={n_valid_seg})')
                plt.legend()
                plt.grid(alpha=0.2)
                plt.tight_layout()
                plt.savefig(blaze_dir / f'order_{aperture_id:03d}_blaze.{fig_format}', dpi=150)
                plt.close()

            # --- Per-order profile plot (2D heatmap + mean) ---
            if profile_x is not None and profile_y is not None:
                fig = plt.figure(figsize=(10, 4.5))
                
                if seg_on_grid is not None and seg_centers is not None and n_seg > 1:
                    ax1 = plt.subplot(121)
                    # seg_on_grid is shape (n_seg, n_y). We transpose to (n_y, n_seg)
                    img_data = seg_on_grid.T
                    # Origin='lower' so y goes up.
                    extent = [seg_centers[0], seg_centers[-1], profile_x[0], profile_x[-1]]
                    im = ax1.imshow(img_data, aspect='auto', origin='lower', extent=extent, cmap='viridis')
                    ax1.set_xlabel('Dispersion X (pixel)')
                    ax1.set_ylabel('Spatial y_norm')
                    ax1.set_title('Profile Variation along Dispersion')
                    fig.colorbar(im, ax=ax1, label='Normalized Flux')
                    
                    ax2 = plt.subplot(122)
                else:
                    ax2 = plt.subplot(111)
                    
                ax2.plot(profile_x, profile_y, color='black', linewidth=1.8, label='Mean profile')
                ax2.set_xlabel('Spatial y_norm')
                ax2.set_ylabel('Normalized Flux')
                ax2.set_title(f'Order {aperture_id} — Mean Profile')
                ax2.grid(alpha=0.2)
                ax2.legend()
                
                plt.suptitle(f'Order {aperture_id} — Cross-Dispersion Profile ({n_seg} segments)')
                plt.tight_layout()
                plt.savefig(profile_dir / f'order_{aperture_id:03d}_profile.{fig_format}', dpi=150)
                plt.close()

            # --- Per-order real-space profile plot ---
            if sampled_realspace:
                plt.figure(figsize=(8, 5))
                for prof in sampled_realspace:
                    y = prof['y']
                    m = prof['model_flux']
                    xcol = prof['x']
                    plt.plot(y, m, linewidth=1.0, label=f'X={xcol}')
                plt.xlabel('Y (physical pixel)')
                plt.ylabel('Model Flux')
                plt.title(f'Order {aperture_id} — Inverse-Mapped Real-Space Profiles')
                if len(sampled_realspace) <= 8:
                    plt.legend(fontsize=8)
                plt.grid(alpha=0.2)
                plt.tight_layout()
                plt.savefig(realspace_dir / f'order_{aperture_id:03d}_realspace_profiles.{fig_format}', dpi=150)
                plt.close()

        # ================================================================
        # Combined overview: all orders' blaze in one figure
        # ================================================================
        if save_plots and all_order_ids:
            n_orders = len(all_order_ids)
            cmap = plt.cm.turbo
            # Reverse: low order IDs (bottom, long wavelength) → red,
            #          high order IDs (top, short wavelength) → blue.
            norm_c = Normalize(vmin=0, vmax=max(1, n_orders - 1))

            def _order_color(i):
                return cmap(norm_c(n_orders - 1 - i))

            # --- All-orders blaze overview ---
            fig, axes = plt.subplots(2, 1, figsize=(14, 8), sharex=True,
                                     gridspec_kw={'height_ratios': [3, 1]},
                                     constrained_layout=True)
            for i, oid in enumerate(all_order_ids):
                diag = self.order_diagnostics[oid]
                bfit = diag.get('blaze_fit')
                if bfit is None:
                    continue
                x = np.arange(bfit.size)
                clr = _order_color(i)
                axes[0].plot(x, bfit, color=clr, linewidth=0.6, alpha=0.8)
                # Ratio: raw / fit (shows fringe residual)
                raw = diag.get('raw_blaze')
                if raw is not None:
                    ratio = np.where((np.isfinite(bfit)) & (bfit > 0),
                                     raw / bfit, np.nan)
                    axes[1].plot(x, ratio, color=clr, linewidth=0.4, alpha=0.6)
            axes[0].set_ylabel('Heavy-Smooth Blaze')
            axes[0].set_title(f'All Orders Blaze Overview ({n_orders} orders)')
            axes[0].grid(alpha=0.15)
            axes[1].set_xlabel('X (column)')
            axes[1].set_ylabel('Raw / Heavy')
            axes[1].axhline(1.0, color='k', linewidth=0.5, linestyle='--')
            axes[1].set_ylim(0.90, 1.10)
            axes[1].grid(alpha=0.15)
            cmap_r = plt.cm.turbo_r
            sm = ScalarMappable(cmap=cmap_r, norm=norm_c)
            sm.set_array([])
            cbar = fig.colorbar(sm, ax=axes, shrink=0.6, pad=0.02)
            order_ids_arr = np.array(all_order_ids)
            cbar.set_ticks(np.linspace(0, n_orders - 1, min(10, n_orders)))
            cbar.set_ticklabels([str(order_ids_arr[int(t)]) for t in np.linspace(0, n_orders - 1, min(10, n_orders))])
            cbar.set_label('Order ID')
            fig.savefig(out_dir / f'all_orders_blaze_overview.{fig_format}', dpi=150)
            plt.close(fig)

            # --- All-orders profile overview ---
            fig, ax = plt.subplots(figsize=(10, 7))
            for i, oid in enumerate(all_order_ids):
                diag = self.order_diagnostics[oid]
                px = diag.get('profile_x')
                py = diag.get('profile_y')
                if px is None or py is None:
                    continue
                clr = _order_color(i)
                ax.plot(px, py, color=clr, linewidth=0.7, alpha=0.8)
            ax.set_xlabel('y_norm')
            ax.set_ylabel('Normalized Flux')
            ax.set_title(f'All Orders Cross-Dispersion Profile Overview ({n_orders} orders)')
            ax.grid(alpha=0.15)
            sm = ScalarMappable(cmap=cmap_r, norm=norm_c)
            sm.set_array([])
            cbar = fig.colorbar(sm, ax=ax, shrink=0.7, pad=0.02)
            cbar.set_ticks(np.linspace(0, n_orders - 1, min(10, n_orders)))
            cbar.set_ticklabels([str(order_ids_arr[int(t)]) for t in np.linspace(0, n_orders - 1, min(10, n_orders))])
            cbar.set_label('Order ID')
            fig.tight_layout()
            fig.savefig(out_dir / f'all_orders_profile_overview.{fig_format}', dpi=150)
            plt.close(fig)

            logger.info(f"Step 4 diagnostics saved to {out_dir} "
                        f"({n_orders} orders, format={fig_format})")



def process_flat_correction_stage(science_image: np.ndarray,
                                   flat_field: FlatField,
                                   output_dir_base: str,
                                   apertures: Optional[ApertureSet] = None,
                                   science_name: str = 'science',
                                   # Blaze / profile parameters for on-demand model build
                                   blaze_knot_spacing: int = 500,
                                   blaze_edge_nknots: int = 6,
                                   width_smooth_window: int = 41,
                                   profile_bin_step: float = 0.01,
                                   n_profile_segments: int = 100,
                                   profile_smooth_sigma: float = 6.0,
                                   pixel_flat_min: float = 0.5,
                                   pixel_flat_max: float = 1.5,
                                   save_plots: bool = True,
                                   fig_format: str = 'png') -> np.ndarray:
    """
    Step 4: apply the 2D pixel-flat correction map to a science image.
    """

    processor = None
    if flat_field.pixel_flat is None:
        if apertures is None:
            raise RuntimeError("Apertures are required to build the flat model, but none were provided.")
        logger.info("  pixel_flat not found. Building 2D flat model (blaze × profile → pixel flat)...")
        processor = FlatCorrectionModelBuilder()
        processor.flat_data = flat_field.flat_data
        processor.flat_mask = flat_field.flat_mask

        clean_flat = flat_field.flat_data
        if flat_field.scattered_light is not None:
            clean_flat = np.clip(
                flat_field.flat_data.astype(np.float32) - flat_field.scattered_light.astype(np.float32),
                1e-6,
                None,
            )

        flat_sens, blaze_profiles, cross_profiles, smoothed_model, pixel_flat, illum = (
            processor.build_order_response_map(
                apertures, source_image=clean_flat,
                blaze_knot_spacing=blaze_knot_spacing,
                blaze_edge_nknots=blaze_edge_nknots,
                width_smooth_window=width_smooth_window,
                profile_bin_step=profile_bin_step,
                n_profile_segments=n_profile_segments,
                profile_smooth_sigma=profile_smooth_sigma,
                pixel_flat_min=pixel_flat_min,
                pixel_flat_max=pixel_flat_max)
        )
        flat_field.flat_sens = flat_sens
        flat_field.blaze_profiles = blaze_profiles
        flat_field.cross_profiles = cross_profiles
        flat_field.smoothed_model = smoothed_model
        flat_field.pixel_flat = pixel_flat
        flat_field.illumination_flat = illum

        # 保存 pixel_flat_2d.fits
        pf_path = Path(output_dir_base) / 'step4_flat_corrected' / 'pixel_flat_2d.fits'
        pf_path.parent.mkdir(parents=True, exist_ok=True)
        from src.utils.fits_io import write_fits_image
        write_fits_image(str(pf_path), pixel_flat.astype(np.float32), dtype='float32')
        logger.info(f"  Built and saved pixel flat to {pf_path.name}")

        # 保存 Step 4 诊断文件（blaze 曲线、交叉色散轮廓等）

        processor.save_step4_diagnostics(
            pf_path.parent,
            clean_flat=clean_flat,
            clean_flat_corrected=None,
            fig_format=fig_format,
            save_plots=save_plots,
        )

    flat_corr = flat_field.pixel_flat
    if flat_corr is None:
        flat_corr = flat_field.flat_corr_2d
    if flat_corr is None:
        raise RuntimeError("No flat correction map available. Run Step 2 first.")

    if flat_corr.shape != science_image.shape:
        raise RuntimeError(
            f"Science image shape {science_image.shape} does not match "
            f"flat correction shape {flat_corr.shape}"
        )

    safe_flat = flat_corr.astype(np.float32)
    bad = (~np.isfinite(safe_flat)) | (safe_flat <= 0.05)
    if np.any(bad):
        safe_flat = safe_flat.copy()
        safe_flat[bad] = 1.0

    corrected = science_image.astype(np.float32) / safe_flat

    out_dir = Path(output_dir_base) / 'step4_flat_corrected'
    out_dir.mkdir(parents=True, exist_ok=True)

    # Save all flat-model artefacts (passing processor if it was freshly built here to save its diagnostics).
    save_flat_correction_products(flat_field, apertures, out_dir, processor=processor,
                                  fig_format=fig_format, save_plots=save_plots)

    write_fits_image(
        str(out_dir / f'{science_name}_science_flat2d_corrected.fits'),
        corrected,
        dtype='float32',
    )

    if save_plots:
        plot_2d_image_to_file(
            corrected,
            str(out_dir / f'{science_name}_science_flat2d_corrected.{fig_format}'),
            'Science After 2D Flat Correction',
        )
        plot_2d_image_to_file(
            safe_flat,
            str(out_dir / f'flat_correction_used.{fig_format}'),
            '2D Flat Correction Used',
        )
        write_fits_image(
            str(out_dir / 'flat_correction_used.fits'),
            safe_flat.astype(np.float32),
            dtype='float32',
        )

    logger.info(f"✓ 2D flat correction complete: {science_name}")
    return corrected
