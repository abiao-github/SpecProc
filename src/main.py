"""
SpecProc: Main application entry point.

PyQt-based GUI for echelle spectrograph FITS data reduction.
"""

import sys
import argparse
import logging
from logging.handlers import RotatingFileHandler
from pathlib import Path

from PyQt5.QtWidgets import QApplication

from src.config.config_manager import ConfigManager
from src.gui.main_window import MainWindow


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


def main(config_path: str = None, mode: str = 'gui'):
    """Main application entry point.
    
    Args:
        config_path: Path to configuration file (optional)
        mode: Run mode ('gui' or 'cli')
    """
    # Parse command line arguments first (before QApplication)
    parser = argparse.ArgumentParser(
        description='SpecProc - Spectral Processing Tool for Echelle Spectrographs',
        epilog='Examples:\n  specproc              # Start GUI mode (default)\n  specproc --mode cli   # Start CLI mode\n  specproc --config custom.cfg  # Use custom config file',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--mode', choices=['gui', 'cli'], default='gui',
                        help='Run mode: gui (default) or cli')
    parser.add_argument('--config', type=str, default=None,
                        help='Path to configuration file')
    parser.add_argument('--version', action='version', version='SpecProc 1.0.0')
    args = parser.parse_args()

    # Override with function arguments if provided
    if config_path is not None:
        args.config = config_path
    if mode is not None:
        args.mode = mode

    # Setup logging
    setup_logging()

    # Load or create configuration
    config = ConfigManager(args.config)

    if args.mode == 'cli':
        # CLI mode - import and run from run.py
        import run
        run.run_cli(args.config)
    else:
        # GUI mode
        app = QApplication(sys.argv)
        window = MainWindow(config)
        window.show()
        sys.exit(app.exec_())


if __name__ == '__main__':
    main()
