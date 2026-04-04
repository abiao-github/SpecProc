#!/usr/bin/env python3
"""
Test GUI with processing steps selection.
"""

import sys
import os

sys.path.insert(0, os.path.dirname(__file__))

from PyQt5.QtWidgets import QApplication
from src.gui.main_window import MainWindow
from src.config.config_manager import ConfigManager
import logging

logging.basicConfig(level=logging.INFO)

def test_gui_with_steps():
    """Test GUI with processing steps selection."""
    print("Testing GUI with processing steps selection...")

    # Create config
    config = ConfigManager()
    print("✓ Configuration created")

    # Create app
    app = QApplication(sys.argv)
    print("✓ QApplication created")

    # Create main window
    window = MainWindow(config)
    print("✓ MainWindow created")

    # Verify stage checkboxes exist
    assert hasattr(window, 'stage_checkboxes'), "stage_checkboxes not found"
    print(f"✓ Stage checkboxes exists: {len(window.stage_checkboxes)} stages")

    # Verify each stage checkbox
    expected_stages = ['stage_0', 'stage_1', 'stage_2', 'stage_3', 'stage_4', 'stage_5']
    for stage_id in expected_stages:
        assert stage_id in window.stage_checkboxes, f"Stage {stage_id} not found"
        checkbox = window.stage_checkboxes[stage_id]
        assert checkbox.isChecked(), f"Stage {stage_id} should be checked by default"
    print(f"✓ All {len(expected_stages)} stage checkboxes verified and checked by default")

    # Verify new buttons exist
    assert hasattr(window, 'run_selected_btn'), "run_selected_btn not found"
    print("✓ 'Run Selected Steps' button exists")

    # Verify new methods exist
    assert hasattr(window, '_run_selected_pipeline'), "_run_selected_pipeline not found"
    assert hasattr(window, '_execute_selected_stages'), "_execute_selected_stages not found"
    assert hasattr(window, '_select_all_stages'), "_select_all_stages not found"
    assert hasattr(window, '_clear_all_stages'), "_clear_all_stages not found"
    print("✓ All new methods exist")

    # Test select all
    window._clear_all_stages()
    for checkbox in window.stage_checkboxes.values():
        assert not checkbox.isChecked(), "All checkboxes should be unchecked after clear"
    print("✓ 'Clear All' functionality works")

    window._select_all_stages()
    for checkbox in window.stage_checkboxes.values():
        assert checkbox.isChecked(), "All checkboxes should be checked after select all"
    print("✓ 'Select All' functionality works")

    # Test individual stage control
    window.stage_checkboxes['stage_0'].setChecked(False)
    assert not window.stage_checkboxes['stage_0'].isChecked(), "Stage 0 should be unchecked"
    assert window.stage_checkboxes['stage_1'].isChecked(), "Other stages should remain checked"
    print("✓ Individual stage control works")

    # Show window
    window.show()
    print("✓ MainWindow shown")

    print("\n✅ GUI processing steps selection test PASSED!")

    return 0

if __name__ == '__main__':
    try:
        sys.exit(test_gui_with_steps())
    except Exception as e:
        print(f"\n❌ GUI test FAILED: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
