"""
SpecProc: Main application entry point.

PyQt-based GUI for echelle spectrograph FITS data reduction.
"""

import sys
import argparse
import logging
import shutil
import time
from logging.handlers import RotatingFileHandler
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from astropy.io import fits

from specproc.config.config_manager import ConfigManager
from specproc.core.processing_pipeline import ProcessingPipeline


COPY_DEFAULT_CONFIG_SENTINEL = '__COPY_DEFAULT_CONFIG__'


def setup_logging(log_file: str = 'specproc.log'):
    """Setup application logging."""
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)  # Set to INFO to reduce debug noise

    # Reduce matplotlib logging - set more specific loggers
    logging.getLogger('matplotlib').setLevel(logging.WARNING)
    logging.getLogger('matplotlib.font_manager').setLevel(logging.WARNING)
    logging.getLogger('matplotlib.pyplot').setLevel(logging.WARNING)
    logging.getLogger('PIL').setLevel(logging.WARNING)
    
    # Also disable font manager's verbose output
    import matplotlib
    matplotlib.rcParams['font.family'] = 'sans-serif'
    matplotlib.rcParams['font.sans-serif'] = ['DejaVu Sans']
    matplotlib.rcParams['axes.unicode_minus'] = False

    # File handler
    file_handler = RotatingFileHandler(
        log_file, maxBytes=10*1024*1024, backupCount=5
    )
    file_handler.setLevel(logging.INFO)  # File also at INFO level

    # Console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)

    # Formatter
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    file_handler.setFormatter(formatter)
    console_handler.setFormatter(formatter)

    # Add handlers
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)

    logger.info("=" * 60)
    logger.info("SpecProc Application Started")
    logger.info("=" * 60)


def _collect_input_fits_files(inputs: List[str]) -> List[str]:
    """Collect FITS files from positional CLI inputs (files and/or directories)."""
    paths = inputs if inputs else ['.']
    discovered: List[Path] = []

    for item in paths:
        p = Path(item).expanduser()
        if p.is_dir():
            # Recursively scan directory and subdirectories for FITS files.
            for pattern in ('*.fits', '*.fit', '*.FITS', '*.FIT'):
                discovered.extend(sorted(p.rglob(pattern)))
        elif p.is_file() and p.suffix.lower() in ('.fits', '.fit'):
            discovered.append(p)

    # Remove duplicates while preserving deterministic order.
    unique_files: List[str] = []
    seen = set()
    for file_path in sorted(discovered):
        key = str(file_path.resolve())
        if key not in seen:
            seen.add(key)
            unique_files.append(str(file_path.resolve()))

    return unique_files


def _read_frame_type_hint(fits_path: str) -> Tuple[str, Optional[float]]:
    """Read frame-type hints and exposure time from FITS headers."""
    collected: List[str] = []
    exposure_value: Optional[float] = None

    with fits.open(fits_path) as hdul:
        for idx in range(len(hdul)):
            header = getattr(hdul[idx], 'header', None)
            if header is None:
                continue

            # Prioritize dedicated type keywords.
            for key in ('OBSTYPE', 'IMAGETYP', 'OBJECT'):
                value = str(header.get(key, '')).strip()
                if value:
                    collected.append(value)

            if exposure_value is None:
                raw_exptime = header.get('EXPTIME', header.get('EXPOSURE', None))
                try:
                    if raw_exptime is not None:
                        exposure_value = float(raw_exptime)
                except (TypeError, ValueError):
                    pass

            # Some observatories store frame type only in COMMENT.
            comments = header.get('COMMENT', [])
            if isinstance(comments, str):
                if comments.strip():
                    collected.append(comments.strip())
            else:
                for comment in comments:
                    c = str(comment).strip()
                    if c:
                        collected.append(c)

    # Deduplicate while preserving order.
    seen = set()
    unique_values: List[str] = []
    for value in collected:
        key = value.lower()
        if key not in seen:
            seen.add(key)
            unique_values.append(value)

    return ' | '.join(unique_values), exposure_value


def _classify_obstype(obstype: str) -> Optional[str]:
    """Map header type hints to one of: bias, flat, calib, science."""
    normalized = obstype.strip().lower().replace('_', ' ')
    if not normalized:
        return None

    bias_tokens = {'bias', 'zero'}
    flat_tokens = {'flat', 'quartz', 'halogen', 'continuum', 'tungsten'}
    calib_tokens = {'thar', 'lamp', 'arc', 'comparison', 'comp', 'fear', 'ne', 'ar', 'he'}
    science_tokens = {'science', 'object', 'target', 'star', 'source'}

    if any(tok in normalized for tok in bias_tokens):
        return 'bias'
    if any(tok in normalized for tok in flat_tokens):
        return 'flat'
    if any(tok in normalized for tok in calib_tokens):
        return 'calib'
    if any(tok in normalized for tok in science_tokens):
        return 'science'

    return None


def run_cli_mode(config_path: Optional[str], inputs: List[str]) -> Dict:
    """Run full pipeline in CLI mode using OBSTYPE-based FITS classification."""
    logger = logging.getLogger(__name__)
    started_at = time.perf_counter()
    fits_files = _collect_input_fits_files(inputs)
    if not fits_files:
        raise FileNotFoundError(
            'No FITS files found. Provide FITS files or directories, e.g. '\
            'specproc --mode cli 20241102_hrs'
        )

    grouped: Dict[str, List[str]] = {
        'bias': [],
        'flat': [],
        'calib': [],
        'science': [],
    }
    unknown: List[Tuple[str, str]] = []

    for file_path in fits_files:
        try:
            obstype, exposure_value = _read_frame_type_hint(file_path)
            category = _classify_obstype(obstype)
            if category is None and exposure_value is not None and exposure_value > 0:
                category = 'science'

            if category is None:
                unknown.append((file_path, obstype))
                continue
            grouped[category].append(file_path)
        except Exception as exc:
            logger.warning(f"Skip unreadable FITS file {file_path}: {exc}")

    if unknown:
        for file_path, obstype in unknown:
            logger.warning(
                f"Unrecognized frame type for {Path(file_path).name}: '{obstype}'. File is skipped."
            )

    missing = [name for name in ('bias', 'flat', 'calib', 'science') if not grouped[name]]
    if missing:
        counts = ', '.join([f"{key}={len(val)}" for key, val in grouped.items()])
        raise RuntimeError(
            f"OBSTYPE classification incomplete ({counts}). Missing: {', '.join(missing)}."
        )

    logger.info(
        "CLI classification result: "
        f"bias={len(grouped['bias'])}, flat={len(grouped['flat'])}, "
        f"calib={len(grouped['calib'])}, science={len(grouped['science'])}"
    )

    processed_files = grouped['bias'] + grouped['flat'] + grouped['calib'] + grouped['science']
    total_bytes = sum(Path(file_path).stat().st_size for file_path in processed_files if Path(file_path).exists())
    total_mb = total_bytes / (1024.0 * 1024.0)

    config = ConfigManager(config_path)
    pipeline = ProcessingPipeline(config)
    result = pipeline.run_full_pipeline(
        raw_image_paths=grouped['science'],
        bias_filenames=grouped['bias'],
        flat_filenames=grouped['flat'],
        calib_filenames=grouped['calib'],
    )

    elapsed_minutes = (time.perf_counter() - started_at) / 60.0
    speed = total_mb / elapsed_minutes if elapsed_minutes > 0 else 0.0
    logger.info(
        f"CLI processing summary: elapsed={elapsed_minutes:.2f} min, "
        f"input_size={total_mb:.2f} MB, throughput={speed:.2f} MB/min"
    )

    return result


def copy_default_config_to_cwd(target_filename: str = 'specproc.cfg') -> Path:
    """Copy packaged default config into current working directory."""
    default_config = Path(__file__).parent / 'config' / 'default_config.cfg'
    if not default_config.exists():
        raise FileNotFoundError(f'Default config file not found: {default_config}')

    target_path = Path.cwd() / target_filename
    if target_path.exists():
        print(f'Config already exists: {target_path}')
        return target_path

    shutil.copyfile(default_config, target_path)
    print(f'Default config copied to: {target_path}')
    return target_path


def main(config_path: Optional[str] = None, mode: Optional[str] = None):
    """Main application entry point.
    
    Args:
        config_path: Path to configuration file (optional)
        mode: Run mode ('gui' or 'cli')
    """
    # Parse command line arguments first (before QApplication)
    parser = argparse.ArgumentParser(
        description='SpecProc - Spectral Processing Tool for Echelle Spectrographs',
        epilog='Examples:\n  specproc                        # Start GUI mode\n  specproc 20241102_hrs           # Auto CLI mode with directory input\n  specproc --config               # Copy default config to current directory\n  specproc --mode cli 20241102_hrs\n  specproc --mode cli --config custom.cfg ./raw_data',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--mode', choices=['gui', 'cli'], default=None,
                        help='Run mode: gui or cli (auto: cli when inputs are provided, otherwise gui)')
    parser.add_argument('--config', nargs='?', const=COPY_DEFAULT_CONFIG_SENTINEL, default=None,
                        help='Path to configuration file; pass without value to copy default config')
    parser.add_argument('inputs', nargs='*',
                        help='FITS files or directories (CLI mode). If omitted, current directory is used.')
    parser.add_argument('--version', action='version', version='SpecProc 1.0.0')
    args = parser.parse_args()

    # Override with function arguments if provided
    if config_path is not None:
        args.config = config_path
    if mode is not None:
        args.mode = mode

    if args.config == COPY_DEFAULT_CONFIG_SENTINEL:
        copy_default_config_to_cwd()
        return

    # Auto mode selection:
    # - If positional inputs are provided, default to CLI
    # - Otherwise default to GUI
    if args.mode is None:
        args.mode = 'cli' if args.inputs else 'gui'

    # Setup logging
    setup_logging()

    if args.mode == 'cli':
        run_cli_mode(args.config, args.inputs)
    else:
        # GUI mode
        from PyQt5.QtWidgets import QApplication
        from specproc.gui.main_window import MainWindow

        config = ConfigManager(args.config)
        app = QApplication(sys.argv)
        window = MainWindow(config)
        window.show()
        sys.exit(app.exec_())


if __name__ == '__main__':
    main()
