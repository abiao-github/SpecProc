"""
Main PyQt GUI window for SpecProc application.

Provides the primary user interface with file management, process control,
and result visualization.
"""

import sys
import logging
from pathlib import Path
from typing import List, Optional
import json

from PyQt5.QtWidgets import (
    QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QLabel, QFileDialog, QListWidget, QListWidgetItem,
    QProgressBar, QTextEdit, QTabWidget, QSplitter, QMessageBox,
    QMenuBar, QMenu, QStatusBar, QCheckBox, QGroupBox, QDialog
)
from PyQt5.QtCore import Qt, QThread, pyqtSignal, QTimer
from PyQt5.QtGui import QIcon, QFont

from src.config.config_manager import ConfigManager
from src.core.processing_pipeline import ProcessingPipeline

logger = logging.getLogger(__name__)


class ProcessingWorker(QThread):
    """Worker thread for pipeline execution."""

    progress_updated = pyqtSignal(float, str)
    execution_complete = pyqtSignal(bool, str)

    def __init__(self, pipeline: ProcessingPipeline, raw_image: str,
                 biases: List[str], flats: List[str], calib: str):
        """
        Initialize worker thread.

        Args:
            pipeline: ProcessingPipeline instance
            raw_image: Science image path
            biases: List of bias frame paths
            flats: List of flat frame paths
            calib: Wavelength calibration frame path
        """
        super().__init__()
        self.pipeline = pipeline
        self.raw_image = raw_image
        self.biases = biases
        self.flats = flats
        self.calib = calib

    def run(self):
        """Execute pipeline in worker thread."""
        try:
            # Set progress callback
            self.pipeline.set_progress_callback(
                lambda p, s: self.progress_updated.emit(p, s)
            )

            # Run pipeline
            self.pipeline.run_full_pipeline(
                self.raw_image, self.biases, self.flats, self.calib
            )

            self.execution_complete.emit(True, "Pipeline execution completed successfully")
        except Exception as e:
            logger.error(f"Pipeline execution error: {e}")
            self.execution_complete.emit(False, f"Error: {str(e)}")


class MainWindow(QMainWindow):
    """Main application window."""

    def __init__(self, config: Optional[ConfigManager] = None):
        """
        Initialize main window.

        Args:
            config: Optional configuration manager
        """
        super().__init__()
        self.config = config or ConfigManager()
        self.pipeline = ProcessingPipeline(self.config)
        self.worker_thread = None

        # File lists
        self.raw_files = []
        self.bias_files = []
        self.flat_files = []
        self.calib_file = []

        self.init_ui()
        self.setup_logging()

    def init_ui(self):
        """Initialize user interface."""
        self.setWindowTitle("SpecProc - Spectral Data Reduction")
        self.setGeometry(100, 100, 1400, 900)

        # Create main widget
        main_widget = QWidget()
        self.setCentralWidget(main_widget)

        # Main layout
        layout = QHBoxLayout()
        main_widget.setLayout(layout)

        # Left panel: file management
        left_panel = self._create_file_panel()
        layout.addWidget(left_panel)

        # Right panel: process control & results
        right_panel = self._create_process_panel()
        layout.addWidget(right_panel)

        # Create menu bar
        self._create_menu_bar()

        # Create status bar
        self.statusBar = QStatusBar()
        self.setStatusBar(self.statusBar)
        self.statusBar.showMessage("Ready")

    def _create_file_panel(self) -> QWidget:
        """Create file management panel."""
        panel = QWidget()
        layout = QVBoxLayout()
        layout.setContentsMargins(5, 5, 5, 5)
        layout.setSpacing(8)
        panel.setLayout(layout)

        # ===== BIAS FRAMES =====
        layout.addWidget(QLabel("📊 Bias Frames:"))

        bias_btn_layout = QHBoxLayout()
        bias_add_btn = QPushButton("➕ Add Bias Frames...")
        bias_add_btn.clicked.connect(self._select_bias_files)
        bias_del_btn = QPushButton("❌ Remove Selected")
        bias_del_btn.clicked.connect(self._remove_bias_file)
        bias_btn_layout.addWidget(bias_add_btn)
        bias_btn_layout.addWidget(bias_del_btn)
        layout.addLayout(bias_btn_layout)

        self.bias_list = QListWidget()
        self.bias_list.setMinimumHeight(40)
        self.bias_list.setSelectionMode(QListWidget.ExtendedSelection)
        layout.addWidget(self.bias_list, 1)

        # ===== FLAT FRAMES =====
        layout.addWidget(QLabel("📹 Flat Frames:"))

        flat_btn_layout = QHBoxLayout()
        flat_add_btn = QPushButton("➕ Add Flat Frames...")
        flat_add_btn.clicked.connect(self._select_flat_files)
        flat_del_btn = QPushButton("❌ Remove Selected")
        flat_del_btn.clicked.connect(self._remove_flat_file)
        flat_btn_layout.addWidget(flat_add_btn)
        flat_btn_layout.addWidget(flat_del_btn)
        layout.addLayout(flat_btn_layout)

        self.flat_list = QListWidget()
        self.flat_list.setMinimumHeight(40)
        self.flat_list.setSelectionMode(QListWidget.ExtendedSelection)
        layout.addWidget(self.flat_list, 1)

        # ===== CALIBRATION FRAMES =====
        layout.addWidget(QLabel("🔬 Wavelength Calibration Frames (ThAr/Ar/Ne):"))

        calib_btn_layout = QHBoxLayout()
        calib_add_btn = QPushButton("➕ Add Calibration Frame(s)...")
        calib_add_btn.clicked.connect(self._select_calib_files)
        calib_del_btn = QPushButton("❌ Remove Selected")
        calib_del_btn.clicked.connect(self._remove_calib_file)
        calib_btn_layout.addWidget(calib_add_btn)
        calib_btn_layout.addWidget(calib_del_btn)
        layout.addLayout(calib_btn_layout)

        self.calib_list = QListWidget()
        self.calib_list.setMinimumHeight(40)
        self.calib_list.setSelectionMode(QListWidget.ExtendedSelection)
        layout.addWidget(self.calib_list, 1)

        # ===== SCIENCE IMAGES =====
        layout.addWidget(QLabel("🌟 Science Image(s):"))

        raw_btn_layout = QHBoxLayout()
        raw_add_btn = QPushButton("➕ Add Science Image(s)...")
        raw_add_btn.clicked.connect(self._select_raw_files)
        raw_del_btn = QPushButton("❌ Remove Selected")
        raw_del_btn.clicked.connect(self._remove_raw_file)
        raw_btn_layout.addWidget(raw_add_btn)
        raw_btn_layout.addWidget(raw_del_btn)
        layout.addLayout(raw_btn_layout)

        self.raw_list = QListWidget()
        self.raw_list.setMinimumHeight(40)
        self.raw_list.setSelectionMode(QListWidget.ExtendedSelection)
        layout.addWidget(self.raw_list, 1)

        layout.addStretch()

        return panel

    def _create_process_panel(self) -> QWidget:
        """Create processing control panel."""
        panel = QWidget()
        layout = QVBoxLayout()
        panel.setLayout(layout)

        # ===== PROCESS CONTROL BUTTONS =====
        button_layout = QHBoxLayout()

        self.settings_btn = QPushButton("⚙ Settings")
        self.settings_btn.clicked.connect(self._show_settings)
        button_layout.addWidget(self.settings_btn)

        self.run_all_btn = QPushButton("▶️ Run All Steps")
        self.run_all_btn.clicked.connect(self._run_full_pipeline)
        button_layout.addWidget(self.run_all_btn)

        self.stop_btn = QPushButton("⏹ Stop")
        self.stop_btn.clicked.connect(self._stop_processing)
        self.stop_btn.setEnabled(False)
        button_layout.addWidget(self.stop_btn)

        layout.addLayout(button_layout)

        # ===== PROGRESS INDICATOR =====
        layout.addWidget(QLabel("Progress:"))
        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 100)
        self.progress_bar.setValue(0)
        layout.addWidget(self.progress_bar)

        # Stage indicator
        self.stage_label = QLabel("Idle")
        stage_font = QFont()
        stage_font.setBold(True)
        self.stage_label.setFont(stage_font)
        layout.addWidget(self.stage_label)

        # ===== PROCESSING STEPS SELECTION =====
        steps_group = QGroupBox("📋 Processing Steps")
        steps_layout = QVBoxLayout()

        # Checkboxes for each stage
        self.stage_checkboxes = {}
        stages = [
            ("stage_0", "Stage 0: Overscan Correction", True),
            ("stage_1", "Stage 1: Bias Correction", True),
            ("stage_2", "Stage 2: Flat Fielding & Order Tracing", True),
            ("stage_3", "Stage 3: Wavelength Calibration", True),
            ("stage_4", "Stage 4: Background Subtraction", True),
            ("stage_5", "Stage 5: Spectrum Extraction", True),
        ]

        for stage_id, stage_name, default_checked in stages:
            checkbox = QCheckBox(stage_name)
            checkbox.setChecked(default_checked)
            self.stage_checkboxes[stage_id] = checkbox
            steps_layout.addWidget(checkbox)

        # Button layout for steps
        steps_button_layout = QHBoxLayout()

        self.run_selected_btn = QPushButton("▶️ Run Selected Steps")
        self.run_selected_btn.clicked.connect(self._run_selected_pipeline)
        steps_button_layout.addWidget(self.run_selected_btn)

        select_all_btn = QPushButton("☑️ Select All")
        select_all_btn.clicked.connect(self._select_all_stages)
        steps_button_layout.addWidget(select_all_btn)

        clear_all_btn = QPushButton("☐ Clear All")
        clear_all_btn.clicked.connect(self._clear_all_stages)
        steps_button_layout.addWidget(clear_all_btn)

        steps_layout.addLayout(steps_button_layout)
        steps_group.setLayout(steps_layout)
        layout.addWidget(steps_group)

        # ===== TABS FOR RESULTS =====
        self.tabs = QTabWidget()

        # Log tab
        self.log_text = QTextEdit()
        self.log_text.setReadOnly(True)
        self.tabs.addTab(self.log_text, "📋 Log")

        # Results tab (placeholder)
        results_widget = QTextEdit()
        results_widget.setText("Results will appear here after processing")
        results_widget.setReadOnly(True)
        self.tabs.addTab(results_widget, "📊 Results")

        layout.addWidget(self.tabs)

        return panel

    def _create_menu_bar(self):
        """Create menu bar."""
        menubar = self.menuBar()

        # File menu
        file_menu = menubar.addMenu("File")
        open_config_action = file_menu.addAction("Open Configuration...")
        open_config_action.triggered.connect(self._open_config)
        settings_action = file_menu.addAction("Settings...")
        settings_action.triggered.connect(self._show_settings)
        file_menu.addSeparator()
        exit_action = file_menu.addAction("Exit")
        exit_action.triggered.connect(self.close)

        # Help menu
        help_menu = menubar.addMenu("Help")
        about_action = help_menu.addAction("About")
        about_action.triggered.connect(self._show_about)

    def _select_bias_files(self):
        """Select bias frame files."""
        files, _ = QFileDialog.getOpenFileNames(
            self, "Select Bias Frames", "", "FITS Files (*.fits *.fit)"
        )
        if files:
            for f in files:
                self.bias_files.append(f)
                self.bias_list.addItem(Path(f).name)
            self.statusBar.showMessage(f"Added {len(files)} bias frames")

    def _remove_bias_file(self):
        """Remove selected bias file from list."""
        for item in self.bias_list.selectedItems():
            try:
                idx = self.bias_list.row(item)
                if idx < len(self.bias_files):
                    self.bias_files.pop(idx)
                    self.bias_list.takeItem(idx)
                self.statusBar.showMessage(f"Bias frames: {len(self.bias_files)} remaining")
            except Exception as e:
                self.statusBar.showMessage(f"Error removing bias file: {e}")

    def _select_flat_files(self):
        """Select flat frame files."""
        files, _ = QFileDialog.getOpenFileNames(
            self, "Select Flat Frames", "", "FITS Files (*.fits *.fit)"
        )
        if files:
            for f in files:
                self.flat_files.append(f)
                self.flat_list.addItem(Path(f).name)
            self.statusBar.showMessage(f"Added {len(files)} flat frames")

    def _remove_flat_file(self):
        """Remove selected flat file from list."""
        for item in self.flat_list.selectedItems():
            try:
                idx = self.flat_list.row(item)
                if idx < len(self.flat_files):
                    self.flat_files.pop(idx)
                    self.flat_list.takeItem(idx)
                self.statusBar.showMessage(f"Flat frames: {len(self.flat_files)} remaining")
            except Exception as e:
                self.statusBar.showMessage(f"Error removing flat file: {e}")

    def _select_calib_files(self):
        """Select wavelength calibration frame(s)."""
        files, _ = QFileDialog.getOpenFileNames(
            self, "Select Wavelength Calibration Frame(s)", "", "FITS Files (*.fits *.fit)"
        )
        if files:
            # Initialize as list if not already
            if not isinstance(self.calib_file, list):
                self.calib_file = []
            # Append new files to existing list
            for f in files:
                self.calib_file.append(f)
                self.calib_list.addItem(Path(f).name)
            self.statusBar.showMessage(f"Added {len(files)} calibration frame(s), total: {len(self.calib_file)}")

    def _remove_calib_file(self):
        """Remove selected calibration file from list."""
        for item in self.calib_list.selectedItems():
            try:
                idx = self.calib_list.row(item)
                if isinstance(self.calib_file, list) and idx < len(self.calib_file):
                    self.calib_file.pop(idx)
                    self.calib_list.takeItem(idx)
                self.statusBar.showMessage(f"Calibration frames: {len(self.calib_file)} remaining")
            except Exception as e:
                self.statusBar.showMessage(f"Error removing calibration file: {e}")

    def _select_raw_files(self):
        """Select science image file(s)."""
        files, _ = QFileDialog.getOpenFileNames(
            self, "Select Science Image(s)", "", "FITS Files (*.fits *.fit)"
        )
        if files:
            # Append new files to existing list
            for f in files:
                self.raw_files.append(f)
                self.raw_list.addItem(Path(f).name)
            self.statusBar.showMessage(f"Added {len(files)} science image(s), total: {len(self.raw_files)}")

    def _remove_raw_file(self):
        """Remove selected science image from list."""
        for item in self.raw_list.selectedItems():
            try:
                idx = self.raw_list.row(item)
                if idx < len(self.raw_files):
                    self.raw_files.pop(idx)
                    self.raw_list.takeItem(idx)
                self.statusBar.showMessage(f"Science images: {len(self.raw_files)} remaining")
            except Exception as e:
                self.statusBar.showMessage(f"Error removing science image: {e}")

    def _run_full_pipeline(self):
        """Execute full processing pipeline."""
        # Validate inputs
        if not self.raw_files:
            QMessageBox.warning(self, "Error", "Please select at least one science image")
            return
        if not self.bias_files:
            QMessageBox.warning(self, "Error", "Please select bias frames")
            return
        if not self.flat_files:
            QMessageBox.warning(self, "Error", "Please select flat frames")
            return

        # Handle calib files (can be list)
        calib_file = None
        if isinstance(self.calib_file, list):
            if self.calib_file:
                # Use first calib file (or most recent by modification time)
                calib_file = sorted(self.calib_file, key=lambda x: Path(x).stat().st_mtime)[-1]
        else:
            calib_file = self.calib_file

        if not calib_file:
            QMessageBox.warning(self, "Error", "Please select a wavelength calibration frame")
            return

        # Ask if user wants to process all or just first
        if len(self.raw_files) > 1:
            auto_process = self.config.get_bool('reduce', 'auto_process_all_images', False)
            if not auto_process:
                reply = QMessageBox.question(
                    self, "Multiple Science Images",
                    f"Found {len(self.raw_files)} science image(s).\n"
                    "Process all of them sequentially?",
                    QMessageBox.Yes | QMessageBox.No
                )
                if reply == QMessageBox.No:
                    # Use only first
                    science_files = self.raw_files[:1]
                else:
                    science_files = self.raw_files
            else:
                science_files = self.raw_files
        else:
            science_files = self.raw_files

        # Disable controls
        self.run_all_btn.setEnabled(False)
        self.stop_btn.setEnabled(True)

        # Process each science image
        self.log_text.append("=" * 60)
        self.log_text.append(f"STARTING PIPELINE - Processing {len(science_files)} image(s)")
        self.log_text.append("=" * 60)

        for idx, science_file in enumerate(science_files):
            self.log_text.append(f"\n[Image {idx+1}/{len(science_files)}] {Path(science_file).name}")

            # Create worker thread
            self.worker_thread = ProcessingWorker(
                self.pipeline,
                science_file,
                self.bias_files,
                self.flat_files,
                calib_file
            )
            self.worker_thread.progress_updated.connect(self._on_progress)
            self.worker_thread.execution_complete.connect(
                lambda success, msg, img_idx=idx: self._on_execution_complete(success, msg, img_idx, len(science_files))
            )
            self.worker_thread.start()
            self.worker_thread.wait()  # Wait for this image to complete before next

        self.run_all_btn.setEnabled(True)
        self.stop_btn.setEnabled(False)

    def _run_selected_pipeline(self):
        """Execute only the selected processing steps."""
        # Check which steps are selected
        selected_stages = []
        for stage_id, checkbox in self.stage_checkboxes.items():
            if checkbox.isChecked():
                selected_stages.append(stage_id)

        if not selected_stages:
            QMessageBox.warning(self, "Error", "Please select at least one processing step")
            return

        # Validate inputs based on selected steps
        if not self.raw_files:
            QMessageBox.warning(self, "Error", "Please select at least one science image")
            return

        # Log selected steps
        self.log_text.append("\n" + "=" * 60)
        self.log_text.append("SELECTED STEPS:")
        stage_names = {
            "stage_0": "Stage 0: Overscan Correction",
            "stage_1": "Stage 1: Bias Correction",
            "stage_2": "Stage 2: Flat Fielding",
            "stage_3": "Stage 3: Wavelength Calibration",
            "stage_4": "Stage 4: Background Subtraction",
            "stage_5": "Stage 5: Spectrum Extraction",
        }
        for stage_id in selected_stages:
            self.log_text.append(f"  ✓ {stage_names.get(stage_id, stage_id)}")
        self.log_text.append("=" * 60)

        # Validate required inputs
        if "stage_1" in selected_stages and not self.bias_files:
            QMessageBox.warning(self, "Error", "Stage 1 (Bias) requires bias frames")
            return
        if "stage_2" in selected_stages and not self.flat_files:
            QMessageBox.warning(self, "Error", "Stage 2 (Flat) requires flat frames")
            return
        if "stage_3" in selected_stages:
            if not isinstance(self.calib_file, list) or not self.calib_file:
                QMessageBox.warning(self, "Error", "Stage 3 (Wavelength) requires calibration frames")
                return

        # Ask about multiple science images
        if len(self.raw_files) > 1:
            auto_process = self.config.get_bool('reduce', 'auto_process_all_images', False)
            if not auto_process:
                reply = QMessageBox.question(
                    self, "Multiple Science Images",
                    f"Found {len(self.raw_files)} science image(s).\n"
                    "Process all of them sequentially?",
                    QMessageBox.Yes | QMessageBox.No
                )
                science_files = self.raw_files if reply == QMessageBox.Yes else self.raw_files[:1]
            else:
                science_files = self.raw_files
        else:
            science_files = self.raw_files

        # Disable controls
        self.run_all_btn.setEnabled(False)
        self.run_selected_btn.setEnabled(False)
        self.stop_btn.setEnabled(True)

        # Get calib file
        calib_file = None
        if isinstance(self.calib_file, list) and self.calib_file:
            calib_file = sorted(self.calib_file, key=lambda x: Path(x).stat().st_mtime)[-1]

        # Create a modified pipeline for selected steps
        self.log_text.append(f"\nProcessing {len(science_files)} science image(s)...")

        for idx, science_file in enumerate(science_files):
            self.log_text.append(f"\n[Image {idx+1}/{len(science_files)}] {Path(science_file).name}")

            # Here we would execute only selected stages
            # For now, we'll pass the selected stages info to a modified processing
            try:
                # Create a custom processing based on selected stages
                self._execute_selected_stages(science_file, selected_stages)
                self._on_execution_complete(True, f"Image {idx+1} processed successfully", idx, len(science_files))
            except Exception as e:
                self._on_execution_complete(False, f"Error processing image {idx+1}: {str(e)}", idx, len(science_files))
                break

        self.run_all_btn.setEnabled(True)
        self.run_selected_btn.setEnabled(True)
        self.stop_btn.setEnabled(False)

    def _execute_selected_stages(self, science_file: str, selected_stages: List[str]):
        """Execute only the selected processing stages."""
        from src.core.overscan_correction import process_overscan_stage
        from src.core.bias_correction import process_bias_stage
        from src.core.flat_fielding import process_flat_stage
        from src.core.wave_calibration import process_wavelength_stage
        from src.core.background_removal import process_background_stage
        from src.core.extraction import process_extraction_stage

        self.log_text.append(f"  Executing {len(selected_stages)} selected stage(s)...")

        # For overscan correction, process all image types
        if "stage_0" in selected_stages:
            self.log_text.append("  ✓ Overscan Correction")
            try:
                all_raw_files = [science_file] + self.bias_files + self.flat_files + self.calib_file
                corrected_files = process_overscan_stage(self.config, all_raw_files, None)
                # Update file paths for subsequent stages
                science_file = corrected_files[0]
                self.bias_files = corrected_files[1:1+len(self.bias_files)]
                self.flat_files = corrected_files[1+len(self.bias_files):1+len(self.bias_files)+len(self.flat_files)]
                self.calib_file = corrected_files[-1] if isinstance(self.calib_file, list) else corrected_files[-1]
            except Exception as e:
                self.log_text.append(f"  ✗ Overscan Correction failed: {e}")

        # For other stages, process only science image
        if "stage_1" in selected_stages:
            self.log_text.append("  ✓ Bias Correction")
            try:
                master_bias = process_bias_stage(self.config, self.bias_files)
            except Exception as e:
                self.log_text.append(f"  ✗ Bias Correction failed: {e}")

        if "stage_2" in selected_stages:
            self.log_text.append("  ✓ Flat Fielding")
            try:
                flat_field, apertures = process_flat_stage(self.config, self.flat_files)
            except Exception as e:
                self.log_text.append(f"  ✗ Flat Fielding failed: {e}")

        if "stage_3" in selected_stages:
            self.log_text.append("  ✓ Wavelength Calibration")
            try:
                wave_calib = process_wavelength_stage(self.config, self.calib_file)
            except Exception as e:
                self.log_text.append(f"  ✗ Wavelength Calibration failed: {e}")

        if "stage_4" in selected_stages:
            self.log_text.append("  ✓ Background Subtraction")
            try:
                process_background_stage(self.config, science_file)
            except Exception as e:
                self.log_text.append(f"  ✗ Background Subtraction failed: {e}")

        if "stage_5" in selected_stages:
            self.log_text.append("  ✓ Spectrum Extraction")
            try:
                process_extraction_stage(self.config, science_file)
            except Exception as e:
                self.log_text.append(f"  ✗ Spectrum Extraction failed: {e}")

    def _select_all_stages(self):
        """Select all processing stages."""
        for checkbox in self.stage_checkboxes.values():
            checkbox.setChecked(True)
        self.statusBar.showMessage("All processing stages selected")

    def _clear_all_stages(self):
        """Clear all processing stage selections."""
        for checkbox in self.stage_checkboxes.values():
            checkbox.setChecked(False)
        self.statusBar.showMessage("All processing stages cleared")

    def _stop_processing(self):
        """Stop processing."""
        if self.worker_thread and self.worker_thread.isRunning():
            self.worker_thread.terminate()
            self.worker_thread.wait()
            self.log_text.append("\n[STOPPED] Pipeline execution halted by user")

        self.run_all_btn.setEnabled(True)
        self.stop_btn.setEnabled(False)

    def _on_progress(self, progress: float, stage: str):
        """Handle progress update."""
        self.progress_bar.setValue(int(progress * 100))
        self.stage_label.setText(f"Current: {stage} ({progress*100:.1f}%)")
        self.log_text.append(f"[{progress*100:5.1f}%] {stage}")

    def _on_execution_complete(self, success: bool, message: str,
                             img_idx: int = 0, total_imgs: int = 1):
        """Handle pipeline completion."""
        self.log_text.append("\n" + "=" * 60)
        self.log_text.append(message)
        self.log_text.append("=" * 60)

        if success:
            if img_idx == total_imgs - 1:
                # Last image processed
                QMessageBox.information(self, "Success",
                                       f"All {total_imgs} image(s) processed successfully!")
                self.progress_bar.setValue(100)
                self.run_all_btn.setEnabled(True)
                self.stop_btn.setEnabled(False)
        else:
            QMessageBox.critical(self, "Error", message)
            self.progress_bar.setValue(0)
            self.run_all_btn.setEnabled(True)
            self.stop_btn.setEnabled(False)

    def _open_config(self):
        """Open configuration file."""
        file, _ = QFileDialog.getOpenFileName(
            self, "Open Configuration", "", "Config Files (*.cfg);;All Files (*)"
        )
        if file:
            try:
                self.config.load(file)
                self.statusBar.showMessage(f"Loaded config: {Path(file).name}")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to load config: {e}")

    def _show_settings(self):
        """Show settings dialog."""
        from src.gui.settings_dialog import SettingsDialog
        dialog = SettingsDialog(self.config, self)
        if dialog.exec_() == QDialog.Accepted:
            self.statusBar.showMessage("Settings saved")

    def _show_about(self):
        """Show about dialog."""
        QMessageBox.information(
            self,
            "About SpecProc",
            "SpecProc v0.1.0\n\nPyQt-based spectral reduction for echelle spectrographs\n\n"
            "Inspired by the gamse package"
        )

    def setup_logging(self):
        """Setup logging to GUI."""
        class QTextEditHandler(logging.Handler):
            def __init__(self, text_edit):
                super().__init__()
                self.text_edit = text_edit

            def emit(self, record):
                msg = self.format(record)
                # Limit log output to prevent GUI slowdown
                if self.text_edit.document().lineCount() > 1000:
                    self.text_edit.clear()
                self.text_edit.append(msg)

        handler = QTextEditHandler(self.log_text)
        formatter = logging.Formatter('[%(levelname)s] %(message)s')
        handler.setFormatter(formatter)
        logging.getLogger().addHandler(handler)

    def closeEvent(self, event):
        """Handle window close."""
        if self.worker_thread and self.worker_thread.isRunning():
            self.worker_thread.terminate()
            self.worker_thread.wait()
        super().closeEvent(event)
