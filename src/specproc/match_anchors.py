"""
``specproc-match-anchors``: build an anchor line file from a lamp spectrum.

The tool detects the strong emission lines in an extracted 1-D lamp spectrum,
matches them against a high-resolution reference lamp spectrum (for example
the NOIRLAB/KPNO FeAr atlas shipped as ``calib_data/linelists/fear.fits``),
fits a low-order pixel-to-wavelength solution and writes a compact anchor CSV
plus a diagnostic figure.

Typical use for the Muztagh-Ata 1.93 m slit spectrograph::

    specproc-match-anchors \
        --lamp output/step5_extraction/RR2000V_FeAr_0.8silt_0.2s_1_spec.fits \
        --reference src/specproc/calib_data/linelists/fear.fits \
        --linelist src/specproc/calib_data/linelists/fear.txt \
        --config MSTG193_slit \
        --output src/specproc/calib_data/telescopes/mstg193slit/mstg193slit_R_ids.csv \
        --plot docs/mstg193slit_R_anchor_match.png

The first two columns of the output (0-based pixel and wavelength in Angstrom)
are exactly what ``[slit] identification_file`` expects.  A few strong,
isolated lines are enough to seed the full-catalogue wavelength calibration
performed later by the slit pipeline.
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Optional, Sequence, Tuple

import numpy as np

from specproc.core.line_identification import load_linelist
from specproc.core.anchor_matching import (
    Reference,
    build_catalog,
    detect_peaks,
    load_observed,
    load_reference,
    plot_diagnostic,
    select_anchors,
    solve_dispersion,
    write_anchor_file,
)

logger = logging.getLogger('specproc.match_anchors')


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog='specproc-match-anchors',
        description='Match strong lamp lines and write an anchor line file.')
    p.add_argument('--lamp', required=True,
                   help='extracted 1-D lamp spectrum (FITS table or text)')
    p.add_argument('--reference', required=True,
                   help='high-resolution reference lamp spectrum (e.g. fear.fits) or line list')
    p.add_argument('--linelist',
                   help='optional line catalogue used for ion labels (e.g. fear.txt)')
    p.add_argument('--wave-range', nargs=2, type=float, metavar=('W0', 'W1'),
                   help='approximate wavelength coverage in Angstrom (first/last '
                        'pixel).  It selects the reference lines used for the match '
                        'and seeds the curvature scan, so it should be within a few '
                        'percent of the true range.  Read from --config or the FITS '
                        'header when omitted; without any hint a blind pattern scan '
                        'is tried.')
    p.add_argument('--config',
                   help='optional SpecProc config/profile providing [slit] wave_start/wave_end '
                        '(e.g. MSTG193_slit or a path to a .cfg file)')
    p.add_argument('--output',
                   help='output anchor CSV (default: <lamp>_anchors.csv)')
    p.add_argument('--plot',
                   help='output diagnostic figure (default: <lamp>_anchors.png)')
    p.add_argument('--poly-order', type=int, default=3, choices=(1, 2, 3, 4),
                   help='order of the preliminary dispersion polynomial (default 3; '
                        'the final pipeline fits its own solution)')
    p.add_argument('--resolution', type=float, default=2000.0,
                   help='approximate resolving power, used to reject blended anchors '
                        '(default 2000; 0 disables the check)')
    p.add_argument('--min-snr', type=float, default=8.0,
                   help='minimum S/N of a lamp peak used as an anchor (default 8)')
    p.add_argument('--match-tol', type=float, default=2.0,
                   help='initial matching tolerance in Angstrom (default 2.0)')
    p.add_argument('--max-residual', type=float, default=0.6,
                   help='maximum residual of an anchor against the preliminary fit '
                        '(default 0.6 A; the wavelength itself comes from the catalogue)')
    p.add_argument('--min-separation', type=float, default=4.0,
                   help='minimum separation from another lamp peak in pixels (default 4)')
    p.add_argument('--max-anchors', type=int, default=30,
                   help='maximum number of anchors (default 30)')
    p.add_argument('--catalog-snr', type=float, default=20.0,
                   help='peak-detection S/N in the reference spectrum (default 20)')
    p.add_argument('--catalog-frac', type=float, default=0.01,
                   help='keep reference lines above this fraction of the strongest '
                        'line inside the requested range (default 0.01)')
    p.add_argument('--catalog-min-snr', type=float, default=50.0,
                   help='absolute minimum reference-line S/N (default 50)')
    p.add_argument('--zoom', nargs=2, type=float, metavar=('W0', 'W1'),
                   help='wavelength range for the zoom panel')
    p.add_argument('--log-level', default='INFO',
                   choices=('DEBUG', 'INFO', 'WARNING', 'ERROR'))
    return p


def _range_from_config(config_path: Optional[str]) -> Optional[Tuple[float, float]]:
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


def _observed_span(peaks, wave_range, n_pixels):
    """Return the wavelength span actually covered by the detected peaks."""
    if not peaks:
        return None
    px = np.array([p.pixel for p in peaks], dtype=float)
    coef = np.polyfit([0, n_pixels - 1], list(wave_range), 1)
    return tuple(sorted((float(np.polyval(coef, np.percentile(px, 2))),
                         float(np.polyval(coef, np.percentile(px, 98))))))


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level.upper()),
                        format='%(asctime)s - %(levelname)s - %(message)s')

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
        logger.info('Approximate wavelength coverage: %.1f-%.1f A', wave_range[0], wave_range[1])
    else:
        logger.info('No wavelength hint; running the blind line-pattern scan.')

    # Build the *strong-line* catalogue from the whole reference atlas.  Do
    # not cut it to --wave-range: that range is extrapolated and approximate,
    # so a hard cut can hide the true strong lines and force a red alias.
    # solve_dispersion treats the range as a soft score bonus instead.
    try:
        catalog = build_catalog(
            ref, wave_range=None, snr=args.catalog_snr,
            strength_frac=args.catalog_frac, min_strength=args.catalog_min_snr,
            linelist_wave=label_wave, linelist_labels=label_text)
    except Exception as exc:
        logger.error('Cannot build the reference catalogue: %s', exc)
        return 2
    if len(catalog) < 5:
        logger.error('Too few strong reference lines (%d); lower --catalog-min-snr',
                     len(catalog))
        return 2

    # Detect the lamp peaks with a slightly lower threshold than the anchor
    # cut so that the pattern search has enough lines to work with.
    detect_snr = min(args.min_snr, 8.0)
    peaks = detect_peaks(flux, snr=detect_snr, min_width=0.5, max_width=10.0,
                         distance=4, baseline_size=61)
    logger.info('Detected %d lamp peaks above S/N=%.1f', len(peaks), detect_snr)
    if len(peaks) < 4:
        logger.error('Too few lamp peaks; lower --min-snr or check the spectrum')
        return 2

    solution = solve_dispersion(
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

    lamp_path = Path(args.lamp)
    out = Path(args.output) if args.output else lamp_path.with_name(lamp_path.stem + '_anchors.csv')
    fig = Path(args.plot) if args.plot else lamp_path.with_name(lamp_path.stem + '_anchors.png')
    note = (f'seeded with {wave_range[0]:.1f}-{wave_range[1]:.1f} A'
            if wave_range else 'blind pattern match')
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


if __name__ == '__main__':  # pragma: no cover
    sys.exit(main())
