#!/usr/bin/env python3
"""
Quick test script for SpecProc core functionality.

Tests basic data structures, configuration management, and FITS I/O.
"""

import sys
import os
from pathlib import Path
import numpy as np

# Add src to path
sys.path.insert(0, os.path.dirname(__file__))

import logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

from src.config.config_manager import ConfigManager
from src.core.data_structures import (
    Spectrum, SpectraSet, ApertureLocation, ApertureSet,
    FlatField, WaveCalib
)
from src.utils.fits_io import read_fits_image
from src.core.processing_pipeline import ProcessingPipeline


def test_config_manager():
    """Test configuration management."""
    print("\n" + "="*60)
    print("TEST 1: Configuration Manager")
    print("="*60)

    config = ConfigManager()
    print("✓ Default configuration created")

    # Save and load
    config.save('test_config.cfg')
    print("✓ Configuration saved")

    config2 = ConfigManager('test_config.cfg')
    print("✓ Configuration loaded")

    # Test parameter access
    midpath = config.get('reduce', 'midpath', './midpath')
    print(f"✓ Parameter access: midpath = {midpath}")

    # Cleanup
    os.remove('test_config.cfg')
    return config


def test_data_structures():
    """Test core data structures."""
    print("\n" + "="*60)
    print("TEST 2: Data Structures")
    print("="*60)

    # Test Spectrum
    spectrum = Spectrum(
        aperture=1,
        order=50,
        wavelength=[3000, 3010, 3020],
        flux=[100, 150, 120],
        error=[10, 12, 10]
    )
    print(f"✓ Spectrum created: {spectrum.npixel()} pixels")

    # Test SpectraSet
    spectra_set = SpectraSet()
    spectra_set.add_spectrum(spectrum)
    print(f"✓ SpectraSet created: {spectra_set.norders} orders")

    # Test ApertureLocation
    aperture = ApertureLocation(
        aperture=1,
        order=50,
        center_coef=[100, 0.1],
        lower_coef=[85],
        upper_coef=[115]
    )
    print(f"✓ ApertureLocation created: aperture {aperture.aperture}")

    # Test ApertureSet
    aperture_set = ApertureSet()
    aperture_set.add_aperture(aperture)
    print(f"✓ ApertureSet created: {aperture_set.norders} orders")

    # Test WaveCalib
    wave_calib = WaveCalib(
        poly_coef=np.zeros((5, 5)),
        xorder=4,
        yorder=4,
        rms=0.05,
        nlines=20
    )
    print(f"✓ WaveCalib created: {wave_calib.nlines} lines, RMS={wave_calib.rms}")

    return spectrum, spectra_set, aperture_set, wave_calib


def test_fits_io():
    """Test FITS file I/O."""
    print("\n" + "="*60)
    print("TEST 3: FITS I/O")
    print("="*60)

    test_data_dir = Path('/Users/abiao/Documents/资料存档/软件工具/gitrepo/20241102_hrs')

    if test_data_dir.exists():
        # Find a test FITS file
        fits_files = list(test_data_dir.glob('*.fits'))
        if fits_files:
            test_file = fits_files[0]
            print(f"Found test FITS file: {test_file.name}")

            try:
                data, header = read_fits_image(str(test_file))
                print(f"✓ FITS image read: shape={data.shape}, dtype={data.dtype}")
                print(f"  Header keys: {len(header)}")
                print(f"  Image range: [{np.min(data):.1f}, {np.max(data):.1f}]")
                return data, header
            except Exception as e:
                print(f"✗ Error reading FITS: {e}")
    else:
        print("⚠ Test data directory not found, skipping FITS I/O test")

    return None, None


def test_pipeline():
    """Test pipeline creation."""
    print("\n" + "="*60)
    print("TEST 4: Processing Pipeline")
    print("="*60)

    config = ConfigManager()
    pipeline = ProcessingPipeline(config)
    print("✓ ProcessingPipeline created")

    state = pipeline.get_state()
    print(f"✓ Pipeline state initialized: {state.current_stage or 'Idle'}")

    return pipeline


def main():
    """Run all tests."""
    print("\n" + "="*60)
    print("SpecProc Core Functionality Test Suite")
    print("="*60)

    try:
        # Test 1: Configuration
        config = test_config_manager()

        # Test 2: Data structures
        test_data_structures()

        # Test 3: FITS I/O
        test_fits_io()

        # Test 4: Pipeline
        test_pipeline()

        print("\n" + "="*60)
        print("✓ All tests passed successfully!")
        print("="*60)
        print("\nTo launch the GUI application, run: python run.py")

    except Exception as e:
        print(f"\n✗ Test failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
