#!/usr/bin/env python3
"""
Test script for overscan correction functionality.
"""

import sys
import os
from pathlib import Path
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))

import logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

from src.utils.overscan import OverscanCorrector, process_overscan_correction
from src.utils.fits_io import read_fits_image


def test_overscan_correction():
    """Test overscan correction with real FITS data."""
    print("\n" + "="*60)
    print("TEST: Overscan Correction")
    print("="*60)

    test_data_dir = Path('/Users/abiao/Documents/资料存档/软件工具/gitrepo/20241102_hrs')
    fits_files = sorted(test_data_dir.glob('*.fits'))

    if not fits_files:
        print("✗ No test files found")
        return False

    try:
        # Read first test image
        test_file = fits_files[0]
        logger.info(f"Reading test image: {test_file.name}")

        raw_image, header = read_fits_image(str(test_file))
        logger.info(f"✓ Image loaded: shape={raw_image.shape}")

        # Create overscan corrector with test configuration
        config = {
            'right': (4110, 4160),  # Right overscan for generic 4160x4136 CCD
        }

        corrector = OverscanCorrector(config)
        logger.info(f"✓ OverscanCorrector created with config: {config}")

        # Extract overscan regions
        regions = corrector.extract_overscan_regions(raw_image)
        logger.info(f"✓ Overscan regions extracted: {list(regions.keys())}")

        if 'right' in regions:
            logger.info(f"  Right overscan shape: {regions['right'].shape}")

        # Estimate overscan bias
        bias = corrector.estimate_overscan_bias(raw_image, method='mean_savgol')
        logger.info(f"✓ Overscan bias estimated: shape={bias.shape}")

        # Apply correction
        corrected = corrector.apply_overscan_correction(raw_image, bias)
        logger.info(f"✓ Overscan correction applied")
        logger.info(f"  Original range: [{np.min(raw_image)}, {np.max(raw_image)}]")
        logger.info(f"  Corrected range: [{np.min(corrected):.1f}, {np.max(corrected):.1f}]")

        # Trim overscan regions
        trimmed = corrector.trim_image(corrected)
        logger.info(f"✓ Image trimmed: {raw_image.shape} → {trimmed.shape}")

        return True

    except Exception as e:
        logger.error(f"✗ Overscan test error: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_overscan_stage():
    """Test overscan correction stage."""
    print("\n" + "="*60)
    print("TEST: Overscan Correction Stage")
    print("="*60)

    from src.config.config_manager import ConfigManager
    from src.core.basic_reduction import OverscanCorrectionStage

    test_data_dir = Path('/Users/abiao/Documents/资料存档/软件工具/gitrepo/20241102_hrs')
    fits_files = sorted(test_data_dir.glob('*.fits'))[:2]

    if not fits_files:
        print("✗ No test files found")
        return False

    try:
        # Create configuration
        config = ConfigManager()
        logger.info("✓ Configuration created")

        # Create stage
        stage = OverscanCorrectionStage(
            overscan_start_column=config.get_int('data', 'overscan_start_column', -1),
            overscan_method=config.get('data', 'overscan_method', 'mean_only'),
            overscan_smooth_window=config.get_int('data', 'overscan_smooth_window', -1),
            overscan_poly_order=config.get_int('data', 'overscan_poly_order', 3),
            overscan_poly_type=config.get('data', 'overscan_poly_type', 'legendre'),
            trim_x_start=config.get_int('data', 'trim_x_start', -1),
            trim_x_end=config.get_int('data', 'trim_x_end', -1),
            trim_y_start=config.get_int('data', 'trim_y_start', -1),
            trim_y_end=config.get_int('data', 'trim_y_end', -1),
            save_plots=config.get_bool('reduce', 'save_plots', True),
            fig_format=config.get('reduce', 'fig_format', 'png')
        )
        logger.info("✓ OverscanCorrectionStage created")

        # Test single file correction
        test_file = str(fits_files[0])
        logger.info(f"Testing file correction: {Path(test_file).name}")

        raw_image, _ = read_fits_image(test_file)
        corrected = stage.correct_image(raw_image, trim=False)

        logger.info(f"✓ Image corrected successfully: shape={corrected.shape}")

        return True

    except Exception as e:
        logger.error(f"✗ Overscan stage test error: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    """Run all overscan tests."""
    print("\n" + "="*60)
    print("SpecProc Overscan Correction Tests")
    print("="*60)

    results = []

    # Test 1: Basic overscan correction
    try:
        results.append(("Overscan Correction", test_overscan_correction()))
    except Exception as e:
        logger.error(f"Test exception: {e}")
        results.append(("Overscan Correction", False))

    # Test 2: Overscan stage
    try:
        results.append(("Overscan Stage", test_overscan_stage()))
    except Exception as e:
        logger.error(f"Test exception: {e}")
        results.append(("Overscan Stage", False))

    # Summary
    print("\n" + "="*60)
    print("TEST SUMMARY")
    print("="*60)

    for name, result in results:
        status = "✓ PASS" if result else "✗ FAIL"
        print(f"{name:30s}: {status}")

    all_pass = all(r for _, r in results)
    if all_pass:
        print("\n✓ All overscan tests passed!")
        return 0
    else:
        print("\n✗ Some tests failed")
        return 1


if __name__ == '__main__':
    sys.exit(main())
