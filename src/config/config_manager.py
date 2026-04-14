"""
Configuration management for SpecProc.

Handles reading and writing INI-format configuration files,
with support for multiple telescope profiles and processing parameters.
"""

import configparser
from pathlib import Path
from typing import Dict, Any, Optional
import logging
import os

logger = logging.getLogger(__name__)


class ConfigManager:
    """Management of telescope and processing configuration."""

    def __init__(self, config_path: Optional[str] = None):
        """
        Initialize configuration manager.

        Args:
            config_path: Path to config file. If None, auto-detects latest .cfg file.
        """
        self.config = configparser.ConfigParser()
        self.config_path = Path(config_path) if config_path else None

        if config_path and Path(config_path).exists():
            self.load(config_path)
        else:
            # Auto-detect latest .cfg file in current directory
            detected_config = self._find_latest_config()
            if detected_config:
                self.load(detected_config)
                self.config_path = detected_config
            else:
                # No config found: load defaults and create a new config in current directory
                self._load_defaults()
                default_user_config_path = Path.cwd() / 'specproc.cfg'
                self.config_path = default_user_config_path

                # If not exists, create empty file and write defaults
                if not default_user_config_path.exists():
                    try:
                        self.save(default_user_config_path)
                    except Exception as e:
                        logger.warning(f"Failed to create default user config at {default_user_config_path}: {e}")

    def _find_latest_config(self) -> Optional[Path]:
        """
        Find the latest .cfg file in current directory (excluding default_config.cfg).

        Returns:
            Path to latest config file, or None if not found.
        """
        current_dir = Path.cwd()
        cfg_files = [f for f in current_dir.glob("*.cfg") if f.name != "default_config.cfg"]

        if not cfg_files:
            return None

        # Return the most recently modified .cfg file
        return max(cfg_files, key=lambda f: f.stat().st_mtime)

    def _load_defaults(self):
        """Load default configuration values."""
        # Load from default_config.cfg in the package directory
        default_config_path = Path(__file__).parent.parent / "default_config.cfg"
        if default_config_path.exists():
            # Load default config without setting config_path
            self.config.read(default_config_path)
        else:
            # Fallback to hardcoded defaults if file not found
            self._load_hardcoded_defaults()

    def _load_hardcoded_defaults(self):
        """Load hardcoded default configuration values."""
        # Data section
        self.config['data'] = {
            'statime_key': 'DATE-OBS',
            'exptime_key': 'EXPTIME',
            'direction': 'xr-',
            'overscan_start_column': '-1',  # Default: no overscan correction
            'overscan_method': 'mean_savgol',
            'trim_x_start': '-1',
            'trim_x_end': '-1',
            'trim_y_start': '-1',
            'trim_y_end': '-1',
        }

        # Reduce section
        self.config['reduce'] = {
            'output_path': './output',
            'mode': 'normal',
            'fig_format': 'png',
            'oned_suffix': '_ods',
            'ncores': 'max',
            'auto_process_all_images': 'yes',
            'cosmic_enabled': 'yes',
            'cosmic_sigma': '5.0',
            'cosmic_window': '5',
            'cosmic_fine_sigma': '2.0',
            'cosmic_line_sigma': '1.5',
            'cosmic_grow_sigma': '2.5',
            'cosmic_maxsize': '8',
        }

        # Bias section
        self.config['reduce.bias'] = {
            'combine_method': 'median',
            'combine_sigma': '3.0',
        }

        # Flat section
        self.config['reduce.flat'] = {
            'combine_method': 'median',
            'q_threshold': '0.5',
            'mosaic_maxcount': '65535',
<<<<<<< HEAD
            'blaze_smooth_method': 'savgol',
            'blaze_smooth_window': '21',
            'blaze_bspline_smooth': '0.5',
=======
            'blaze_knot_spacing': '500',
            'blaze_edge_nknots': '6',
>>>>>>> cef6f04 (	modified:   README.md)
            'width_smooth_window': '41',
            'profile_bin_step': '0.01',
            'pixel_flat_min': '0.5',
            'pixel_flat_max': '1.5',
        }

        # Trace section
        self.config['reduce.trace'] = {
            'separation': '30.0',
<<<<<<< HEAD
            'seed_threshold': '0.30',
=======
            'snr_threshold': '5.0',
>>>>>>> cef6f04 (	modified:   README.md)
            'prominence_scale': '0.50',
            'search_half_scale': '0.45',
            'step_denominator': '220',
            'fill_missing_orders': 'yes',
            'gap_fill_factor': '1.6',
<<<<<<< HEAD
            'fit_method': 'polynomial',
            'bspline_smooth': '0.2',
            'edge_degree': '3',
            'aperture_root_fraction': '0.03',
            'aperture_noise_floor_sigma': '3.0',
            'filling': '0.3',
            'degree': '3',
=======
            'bspline_smooth': '0.2',
            'edge_degree': '3',
            'aperture_boundary_snr': '3.0',
            'n_extend_below': '0',
            'n_extend_above': '0',
            'n_mask_below': '2',
            'n_mask_above': '1',
>>>>>>> cef6f04 (	modified:   README.md)
        }

        # Extract section
        self.config['reduce.extract'] = {
            'method': 'sum',  # 'sum' or 'optimal'
        }

        # Wavelength calibration section
        self.config['reduce.wlcalib'] = {
            'linelist': 'ThAr',
            'search_database': 'yes',
            'use_prev_fitpar': 'yes',
            'xorder': '4',
            'yorder': '4',
            'window_size': '10',
            'clipping': '3.0',
            'auto_selection': 'yes',
            'rms_threshold': '0.1',
            'time_diff': '600',  # seconds
        }

        # Background section
        self.config['reduce.background'] = {
            'method': 'chebyshev',
            'poly_order': '3',
            'smooth_sigma': '20.0',
            'sigma_clip': '3.0',
            'sigma_clip_maxiters': '4',
<<<<<<< HEAD
            'mask_margin_pixels': '3',
=======
            'mask_margin_pixels': '1',
>>>>>>> cef6f04 (	modified:   README.md)
            'bspline_smooth': '1.0',
            'thar_pollution': 'no',
        }

        # Telescope section (generic echelle)
        self.config['telescope'] = {
            'name': 'generic',
            'focal_length': '0.0',  # mm
            'pixel_size': '0.015',  # mm
        }

        # Instrument section
        self.config['instrument'] = {
            'detector': 'generic_ccd',
            'gain': '1.0',
            'readnoise': '5.0',
        }

    def load(self, config_path: str):
        """Load configuration from file."""
        try:
            self.config.read(config_path)
            self.config_path = Path(config_path)
            logger.info(f"Loaded configuration from {config_path}")
        except Exception as e:
            logger.error(f"Error loading config from {config_path}: {e}")
            raise

    def save(self, config_path: Optional[str] = None):
        """Save configuration to file."""
        path = config_path or self.config_path
        if not path:
            # Auto-detect or create default config file
            detected_config = self._find_latest_config()
            if detected_config:
                path = detected_config
            else:
                path = Path.cwd() / "specproc.cfg"

        try:
            with open(path, 'w') as f:
                self.config.write(f)
            logger.info(f"Saved configuration to {path}")
            self.config_path = path  # Update the stored path
        except Exception as e:
            logger.error(f"Error saving config to {path}: {e}")
            raise

    def get_section(self, section: str) -> Dict[str, str]:
        """Get all values in a section."""
        if section not in self.config:
            logger.warning(f"Section [{section}] not found, returning empty dict")
            return {}
        return dict(self.config[section])

    def get(self, section: str, key: str, fallback: Any = None) -> Any:
        """
        Get configuration value with type conversion.

        Args:
            section: Configuration section
            key: Configuration key
            fallback: Default value if not found

        Returns:
            Configuration value
        """
        try:
            if self.config.has_option(section, key):
                return self.config.get(section, key)
            else:
                return fallback
        except Exception:
            return fallback

    def get_int(self, section: str, key: str, fallback: int = 0) -> int:
        """Get integer configuration value."""
        try:
            if self.config.has_option(section, key):
                return self.config.getint(section, key)
            else:
                return fallback
        except Exception:
            return fallback

    def get_float(self, section: str, key: str, fallback: float = 0.0) -> float:
        """Get float configuration value."""
        try:
            if self.config.has_option(section, key):
                return self.config.getfloat(section, key)
            else:
                return fallback
        except Exception:
            return fallback

    def get_bool(self, section: str, key: str, fallback: bool = False) -> bool:
        """Get boolean configuration value."""
        try:
            if self.config.has_option(section, key):
                return self.config.getboolean(section, key)
            else:
                return fallback
        except Exception:
            return fallback

    def set(self, section: str, key: str, value: Any):
        """Set configuration value."""
        if section not in self.config:
            self.config[section] = {}
        self.config[section][key] = str(value)

    def get_rawdata_path(self) -> str:
        """
        Get raw data directory path with proper expansion.

        Returns:
            Expanded path (absolute or relative to current working directory)
        """
        rawdata_path = self.get('data', 'rawdata_path', './rawdata')
        # Expand ~ to user's home directory
        path = Path(rawdata_path).expanduser()
        # If path is relative (not absolute), make it relative to current working directory
        if not path.is_absolute():
            # Resolve to get absolute path relative to current working directory
            path = Path.cwd() / path
        return path.as_posix()

    def get_output_path(self) -> str:
        """
        Get output directory path with proper expansion.

        Returns:
            Expanded path (absolute or relative to current working directory)
        """
        output_path = self.get('reduce', 'output_path', './output')
        path = Path(output_path).expanduser()
        if not path.is_absolute():
            path = Path.cwd() / path
        return path.as_posix()

    def create_directories(self):
        """Create output directories specified in config."""
        output_path = self.get_output_path()

        dirs = [
            output_path,
        ]

        for dir_path in dirs:
            if dir_path:
                Path(dir_path).mkdir(parents=True, exist_ok=True)
                logger.info(f"Created/verified directory: {dir_path}")

    def __repr__(self) -> str:
        """String representation."""
        return f"ConfigManager(path={self.config_path})"


def create_default_config(output_path: str):
    """Create a default configuration file."""
    manager = ConfigManager()
    manager.save(output_path)
    logger.info(f"Created default config file: {output_path}")
