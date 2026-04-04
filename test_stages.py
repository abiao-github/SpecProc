#!/usr/bin/env python3
"""
Quick test script for SpecProc - simple single-stage test.

Tests basic processing with minimal data.
"""

import sys
import os
from pathlib import Path
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))

import logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

from src.config.config_manager import ConfigManager
from src.utils.fits_io import read_fits_image
from src.core.bias_correction import BiasCorrector
from src.core.flat_fielding import FlatFieldProcessor


def test_bias_stage():
    """Test bias correction stage."""
    print("\n" + "="*60)
    print("TEST: Bias Correction Stage")
    print("="*60)

    test_data_dir = Path('/Users/abiao/Documents/资料存档/软件工具/gitrepo/20241102_hrs')
    fits_files = sorted(test_data_dir.glob('*.fits'))[:5]

    if not fits_files:
        print("✗ No test files found")
        return False

    try:
        config = ConfigManager()
        corrector = BiasCorrector(config)

        bias_files = [str(f) for f in fits_files[:3]]
        logger.info(f"Combining {len(bias_files)} bias frames...")

        master_bias, uncertainty = corrector.combine_bias_frames(bias_files)

        if master_bias is not None:
            logger.info(f"✓ Master bias created: shape={master_bias.shape}")
            return True
        else:
            logger.error("✗ Failed to create master bias")
            return False

    except Exception as e:
        logger.error(f"✗ Bias stage error: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_flat_stage():
    """Test flat fielding stage."""
    print("\n" + "="*60)
    print("TEST: Flat Fielding Stage")
    print("="*60)

    test_data_dir = Path('/Users/abiao/Documents/资料存档/软件工具/gitrepo/20241102_hrs')
    fits_files = sorted(test_data_dir.glob('*.fits'))[3:8]

    if not fits_files:
        print("✗ No test files found")
        return False

    try:
        config = ConfigManager()
        processor = FlatFieldProcessor(config)

        flat_files = [str(f) for f in fits_files[:3]]
        logger.info(f"Combining {len(flat_files)} flat frames...")

        flat_data, flat_mask = processor.combine_flat_frames(flat_files)

        if flat_data is not None:
            logger.info(f"✓ Master flat created: shape={flat_data.shape}")

            # Normalize
            flat_norm = processor.normalize_flat()
            logger.info(f"✓ Flat normalized")

            # Extract sensitivity
            flat_sens = processor.extract_sensitivity_map()
            logger.info(f"✓ Sensitivity map extracted")

            # Detect orders
            apertures = processor.detect_orders(threshold=0.3)
            logger.info(f"✓ Orders detected: {apertures.norders} orders")

            return True
        else:
            logger.error("✗ Failed to create master flat")
            return False

    except Exception as e:
        logger.error(f"✗ Flat stage error: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    """Run stage tests."""
    print("\n" + "="*60)
    print("SpecProc Pipeline Stage Tests")
    print("="*60)

    results = []

    # Test bias stage
    try:
        results.append(("Bias Correction", test_bias_stage()))
    except Exception as e:
        logger.error(f"Bias test exception: {e}")
        results.append(("Bias Correction", False))

    # Test flat stage
    try:
        results.append(("Flat Fielding", test_flat_stage()))
    except Exception as e:
        logger.error(f"Flat test exception: {e}")
        results.append(("Flat Fielding", False))

    # Summary
    print("\n" + "="*60)
    print("TEST SUMMARY")
    print("="*60)

    for name, result in results:
        status = "✓ PASS" if result else "✗ FAIL"
        print(f"{name:30s}: {status}")

    all_pass = all(r for _, r in results)
    if all_pass:
        print("\n✓ All tests passed!")
        return 0
    else:
        print("\n✗ Some tests failed")
        return 1


if __name__ == '__main__':
    sys.exit(main())
