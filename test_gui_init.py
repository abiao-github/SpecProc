#!/usr/bin/env python3
"""
Quick test to verify GUI launches correctly.
"""

import sys
import os

sys.path.insert(0, os.path.dirname(__file__))

from PyQt5.QtWidgets import QApplication
from src.gui.main_window import MainWindow
from src.config.config_manager import ConfigManager
import logging

logging.basicConfig(level=logging.INFO)

def test_gui():
    """Test GUI initialization."""
    print("Testing GUI initialization...")

    # Create config
    config = ConfigManager()
    print("✓ Configuration created")

    # Create app
    app = QApplication(sys.argv)
    print("✓ QApplication created")

    # Create main window
    window = MainWindow(config)
    print("✓ MainWindow created")

    # Verify widgets exist
    assert hasattr(window, 'bias_list'), "bias_list not found"
    assert hasattr(window, 'flat_list'), "flat_list not found"
    assert hasattr(window, 'calib_list'), "calib_list not found"
    assert hasattr(window, 'raw_list'), "raw_list not found"
    print("✓ All list widgets exist")

    # Verify methods exist
    assert hasattr(window, '_remove_bias_file'), "_remove_bias_file not found"
    assert hasattr(window, '_remove_flat_file'), "_remove_flat_file not found"
    assert hasattr(window, '_remove_calib_file'), "_remove_calib_file not found"
    assert hasattr(window, '_remove_raw_file'), "_remove_raw_file not found"
    print("✓ All remove methods exist")

    # Check window is visible
    window.show()
    print("✓ MainWindow shown")

    # Get window geometry
    geom = window.geometry()
    print(f"✓ Window size: {geom.width()}x{geom.height()}")

    print("\n✅ GUI initialization test PASSED!")
    print("\n💡 To run the application, execute: python run.py")

    # Don't block - just verify it can be created
    return 0

if __name__ == '__main__':
    try:
        sys.exit(test_gui())
    except Exception as e:
        print(f"\n❌ GUI test FAILED: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
