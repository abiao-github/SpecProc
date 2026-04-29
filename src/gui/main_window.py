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
import numpy as np

from PyQt5.QtWidgets import (
    QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QLabel, QFileDialog, QListWidget, QListWidgetItem,
    QProgressBar, QTextEdit, QTabWidget, QSplitter, QMessageBox,
    QMenuBar, QMenu, QStatusBar, QCheckBox, QGroupBox, QDialog, QFrame
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

    def __init__(self, pipeline: ProcessingPipeline, raw_images: List[str],
                 biases: List[str], flats: List[str], calibs: List[str]):
        """
        Initialize worker thread.

        Args:
            pipeline: ProcessingPipeline instance
            raw_images: List of science image paths
            biases: List of bias frame paths
            flats: List of flat frame paths
            calibs: List of wavelength calibration frame paths
        """
        super().__init__()
        self.pipeline = pipeline
        self.raw_images = raw_images
        self.biases = biases
        self.flats = flats
        self.calibs = calibs

    def run(self):
        """Execute pipeline in worker thread."""
        try:
            # Set progress callback
            self.pipeline.set_progress_callback(
                lambda p, s: self.progress_updated.emit(p, s)
            )

            # Run pipeline
            self.pipeline.run_full_pipeline(
                self.raw_images, self.biases, self.flats, self.calibs
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
        self.last_selected_dir = str(Path.cwd())
        
        # Store original file paths to handle deleted processed files
        # Format: {filename: original_path}
        self.original_file_paths = {}
        
        # Processing state tracking
        self.processing_state = {
            'stage_0_completed': False,  # Overscan correction
            'stage_1_completed': False,  # Bias correction
            'stage_2_completed': False,  # Flat fielding
            'stage_3_completed': False,  # Wavelength calibration
            'stage_4_completed': False,  # Background subtraction
            'stage_5_completed': False,  # Spectrum extraction
            
            # Processed file paths
            'overscan_corrected_files': [],  # Stage 0 output
            'master_bias_path': None,        # Stage 1 output
            'master_flat_path': None,        # Stage 2 output
            'wavelength_calib_path': None,   # Stage 3 output
            'background_model_path': None,   # Stage 4 output
            'extracted_spectra_path': None,  # Stage 5 output
        }

        self.init_ui()
        self.setup_logging()
        
        # Check overscan configuration and update UI accordingly
        self._update_overscan_checkbox_from_config()

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
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(5)
        panel.setLayout(layout)

        # ===== BIAS FRAMES =====
        bias_group = QGroupBox("📊 Bias Frames")
        bias_layout = QVBoxLayout()

        bias_btn_layout = QHBoxLayout()
        bias_add_btn = QPushButton("➕ Add Bias Frames...")
        bias_add_btn.clicked.connect(self._select_bias_files)
        bias_del_btn = QPushButton("❌ Remove Selected")
        bias_del_btn.clicked.connect(self._remove_bias_file)
        bias_btn_layout.addWidget(bias_add_btn)
        bias_btn_layout.addWidget(bias_del_btn)
        bias_layout.addLayout(bias_btn_layout)

        self.bias_list = QListWidget()
        self.bias_list.setMinimumHeight(40)
        self.bias_list.setMaximumHeight(120)
        self.bias_list.setSelectionMode(QListWidget.ExtendedSelection)
        bias_layout.addWidget(self.bias_list, 1)

        bias_group.setLayout(bias_layout)
        layout.addWidget(bias_group)

        # Add separator with spacing
        layout.addSpacing(8)
        separator1 = QFrame()
        separator1.setFrameShape(QFrame.HLine)
        separator1.setFrameShadow(QFrame.Sunken)
        separator1.setContentsMargins(10, 0, 10, 0)
        layout.addWidget(separator1)
        layout.addSpacing(8)

        # ===== FLAT FRAMES =====
        flat_group = QGroupBox("📹 Flat Frames")
        flat_layout = QVBoxLayout()

        flat_btn_layout = QHBoxLayout()
        flat_add_btn = QPushButton("➕ Add Flat Frames...")
        flat_add_btn.clicked.connect(self._select_flat_files)
        flat_del_btn = QPushButton("❌ Remove Selected")
        flat_del_btn.clicked.connect(self._remove_flat_file)
        flat_btn_layout.addWidget(flat_add_btn)
        flat_btn_layout.addWidget(flat_del_btn)
        flat_layout.addLayout(flat_btn_layout)

        self.flat_list = QListWidget()
        self.flat_list.setMinimumHeight(40)
        self.flat_list.setMaximumHeight(120)
        self.flat_list.setSelectionMode(QListWidget.ExtendedSelection)
        flat_layout.addWidget(self.flat_list, 1)

        flat_group.setLayout(flat_layout)
        layout.addWidget(flat_group)

        # Add separator with spacing
        layout.addSpacing(8)
        separator2 = QFrame()
        separator2.setFrameShape(QFrame.HLine)
        separator2.setFrameShadow(QFrame.Sunken)
        separator2.setContentsMargins(10, 0, 10, 0)
        layout.addWidget(separator2)
        layout.addSpacing(8)

        # ===== CALIBRATION FRAMES =====
        calib_group = QGroupBox("🔬 Wavelength Calibration Frames")
        calib_layout = QVBoxLayout()

        calib_btn_layout = QHBoxLayout()
        calib_add_btn = QPushButton("➕ Add Calibration Frame(s)...")
        calib_add_btn.clicked.connect(self._select_calib_files)
        calib_del_btn = QPushButton("❌ Remove Selected")
        calib_del_btn.clicked.connect(self._remove_calib_file)
        calib_btn_layout.addWidget(calib_add_btn)
        calib_btn_layout.addWidget(calib_del_btn)
        calib_layout.addLayout(calib_btn_layout)

        self.calib_list = QListWidget()
        self.calib_list.setMinimumHeight(40)
        self.calib_list.setMaximumHeight(80)
        self.calib_list.setSelectionMode(QListWidget.ExtendedSelection)
        calib_layout.addWidget(self.calib_list, 1)

        calib_group.setLayout(calib_layout)
        layout.addWidget(calib_group)

        # Add separator with spacing
        layout.addSpacing(8)
        separator3 = QFrame()
        separator3.setFrameShape(QFrame.HLine)
        separator3.setFrameShadow(QFrame.Sunken)
        separator3.setContentsMargins(10, 0, 10, 0)
        layout.addWidget(separator3)
        layout.addSpacing(8)

        # ===== SCIENCE IMAGES =====
        science_group = QGroupBox("🌟 Science Images")
        science_layout = QVBoxLayout()

        raw_btn_layout = QHBoxLayout()
        raw_add_btn = QPushButton("➕ Add Science Image(s)...")
        raw_add_btn.clicked.connect(self._select_raw_files)
        raw_del_btn = QPushButton("❌ Remove Selected")
        raw_del_btn.clicked.connect(self._remove_raw_file)
        raw_btn_layout.addWidget(raw_add_btn)
        raw_btn_layout.addWidget(raw_del_btn)
        science_layout.addLayout(raw_btn_layout)

        self.raw_list = QListWidget()
        self.raw_list.setMinimumHeight(40)
        self.raw_list.setSelectionMode(QListWidget.ExtendedSelection)
        science_layout.addWidget(self.raw_list, 1)

        science_group.setLayout(science_layout)
        layout.addWidget(science_group, 2)

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
            ("stage_0", "Step 1: Basic Pre-processing", True),
            ("stage_1", "Step 2: Order Tracing", True),
            ("stage_2", "Step 3: Scattered Light Subtraction", True),
            ("stage_3", "Step 4: 2D Flat-Field Correction", True),
            ("stage_4", "Step 5: 1D Extraction", True),
            ("stage_5", "Step 6: De-blazing", True),
            ("stage_6", "Step 7: Wavelength Calibration", True),
        ]

        # Step1 sub-step toggles: overscan / bias / cosmic
        self.step1_substep_checkboxes = {}
        step1_sub_layout = QHBoxLayout()
        step1_sub_layout.setContentsMargins(22, 0, 0, 0)

        for stage_id, stage_name, default_checked in stages:
            checkbox = QCheckBox(stage_name)
            checkbox.setChecked(default_checked)
            self.stage_checkboxes[stage_id] = checkbox
            steps_layout.addWidget(checkbox)

            # Keep Step1 options immediately after Step1 in the UI.
            if stage_id == "stage_0":
                step1_sub_layout.addWidget(QLabel("("))
                for key, text in [
                    ("overscan", "Overscan"),
                    ("bias", "Bias"),
                    ("cosmic", "Cosmic"),
                ]:
                    cb = QCheckBox(text)
                    cb.setChecked(True)
                    cb.stateChanged.connect(self._refresh_step1_checkbox_text)
                    self.step1_substep_checkboxes[key] = cb
                    step1_sub_layout.addWidget(cb)
                step1_sub_layout.addWidget(QLabel(")"))
                step1_sub_layout.addStretch()
                steps_layout.addLayout(step1_sub_layout)
        
        # Special handling for stage_0 (overscan) checkbox
        if 'stage_0' in self.stage_checkboxes:
            stage_0_checkbox = self.stage_checkboxes['stage_0']
            stage_0_checkbox.stateChanged.connect(self._on_overscan_checkbox_changed)

        self._refresh_step1_checkbox_text()

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
            self, "Select Bias Frames", self.last_selected_dir, "FITS Files (*.fits *.fit)"
        )
        if files:
            self.last_selected_dir = str(Path(files[0]).parent)
            for f in files:
                self.bias_files.append(f)
                self.original_file_paths[Path(f).name] = f  # Save original path
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
            self, "Select Flat Frames", self.last_selected_dir, "FITS Files (*.fits *.fit)"
        )
        if files:
            self.last_selected_dir = str(Path(files[0]).parent)
            for f in files:
                self.flat_files.append(f)
                self.original_file_paths[Path(f).name] = f  # Save original path
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
            self, "Select Wavelength Calibration Frame(s)", self.last_selected_dir, "FITS Files (*.fits *.fit)"
        )
        if files:
            self.last_selected_dir = str(Path(files[0]).parent)
            # Initialize as list if not already
            if not isinstance(self.calib_file, list):
                self.calib_file = []
            # Append new files to existing list
            for f in files:
                self.calib_file.append(f)
                self.original_file_paths[Path(f).name] = f  # Save original path
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
            self, "Select Science Image(s)", self.last_selected_dir, "FITS Files (*.fits *.fit)"
        )
        if files:
            self.last_selected_dir = str(Path(files[0]).parent)
            # Reset processing state when new files are selected
            self._reset_processing_state()
            # Append new files to existing list
            for f in files:
                self.raw_files.append(f)
                self.original_file_paths[Path(f).name] = f  # Save original path
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
    
    def _reset_processing_state(self):
        """Reset processing state when new files are selected or processing is restarted."""
        # Keep original file paths - don't reset them as they point to raw data
        # Reset processed file paths in the lists
        self.processing_state = {
            'stage_0_completed': False,  # Overscan correction
            'stage_1_completed': False,  # Bias correction
            'stage_2_completed': False,  # Flat fielding
            'stage_3_completed': False,  # Wavelength calibration
            'stage_4_completed': False,  # Background subtraction
            'stage_5_completed': False,  # Spectrum extraction
            
            # Processed file paths
            'overscan_corrected_files': [],  # Stage 0 output
            'master_bias_path': None,        # Stage 1 output
            'master_flat_path': None,        # Stage 2 output
            'wavelength_calib_path': None,   # Stage 3 output
            'background_model_path': None,   # Stage 4 output
            'extracted_spectra_path': None,  # Stage 5 output
        }
        self.statusBar.showMessage("Processing state reset")
    
    def _get_overscan_file_status(self) -> tuple:
        """
        Get the status of overscan processed files.
        
        Returns:
            Tuple of (existing_files: list, missing_files: list)
            - existing_files: list of filenames that already have processed files
            - missing_files: list of filenames that need processing
        """
        from pathlib import Path
        
        # Get all files that should be processed
        all_files = []
        
        # Bias files
        for bias_file in self.bias_files:
            all_files.append(Path(bias_file).name)
        
        # Flat files
        for flat_file in self.flat_files:
            all_files.append(Path(flat_file).name)
        
        # Calibration files
        calib_list = self.calib_file if isinstance(self.calib_file, list) else [self.calib_file]
        for calib_file in calib_list:
            if calib_file:
                all_files.append(Path(calib_file).name)
        
        # Science files
        for science_file in self.raw_files:
            all_files.append(Path(science_file).name)
        
        # Check each file in output directory
        output_dir = self.config.get_output_path()
        overscan_dir = Path(output_dir) / 'step1_basic' / 'overscan_corrected'
        
        existing_files = []
        missing_files = []
        
        for filename in all_files:
            processed_file = overscan_dir / filename
            if processed_file.exists():
                existing_files.append(filename)
            else:
                missing_files.append(filename)
        
        return existing_files, missing_files
    
    def _check_processed_files_exist(self, stage_id: str) -> bool:
        """
        Check if processed files for a given stage still exist.
        
        Args:
            stage_id: Stage identifier (e.g., 'stage_0')
        
        Returns:
            True if all processed files exist, False otherwise
        """
        from pathlib import Path
        
        if stage_id == 'stage_0':
            # Check if Stage 0 was marked as completed
            if not self.processing_state.get('stage_0_completed', False):
                return False  # Not completed yet, need to process
            
            # Check all files that should have been processed
            # Get all input files
            all_files = []
            
            # Bias files
            for bias_file in self.bias_files:
                all_files.append(Path(bias_file).name)
            
            # Flat files
            for flat_file in self.flat_files:
                all_files.append(Path(flat_file).name)
            
            # Calibration files
            calib_list = self.calib_file if isinstance(self.calib_file, list) else [self.calib_file]
            for calib_file in calib_list:
                if calib_file:
                    all_files.append(Path(calib_file).name)
            
            # Science files
            for science_file in self.raw_files:
                all_files.append(Path(science_file).name)
            
            # Check each file in output directory
            output_dir = self.config.get_output_path()
            overscan_dir = Path(output_dir) / 'step1_basic' / 'overscan_corrected'
            
            if not overscan_dir.exists():
                return False  # Directory doesn't exist
            
            for filename in all_files:
                if not (overscan_dir / filename).exists():
                    return False  # File missing
            
            return True  # All files exist
        
        elif stage_id == 'stage_1' and self.processing_state['stage_1_completed']:
            return Path(self.processing_state['master_bias_path']).exists() if self.processing_state['master_bias_path'] else False
        
        elif stage_id == 'stage_2' and self.processing_state['stage_2_completed']:
            return Path(self.processing_state['master_flat_path']).exists() if self.processing_state['master_flat_path'] else False
        
        elif stage_id == 'stage_3' and self.processing_state['stage_3_completed']:
            return Path(self.processing_state['wavelength_calib_path']).exists() if self.processing_state['wavelength_calib_path'] else False
        
        elif stage_id == 'stage_4' and self.processing_state['stage_4_completed']:
            return Path(self.processing_state['background_model_path']).exists() if self.processing_state['background_model_path'] else False
        
        elif stage_id == 'stage_5' and self.processing_state['stage_5_completed']:
            return Path(self.processing_state['extracted_spectra_path']).exists() if self.processing_state['extracted_spectra_path'] else False
        
        return False  # Stage not completed or no path stored
    
    def _should_reprocess_stage(self, stage_id: str) -> tuple:
        """
        Determine if a stage should be reprocessed.
        Note: stage_0 (overscan) has its own dialog logic in the processing section.
        
        Args:
            stage_id: Stage identifier
        
        Returns:
            Tuple of (should_reprocess: bool, user_choice_made: bool, force_all: bool)
            - should_reprocess: True if stage should be reprocessed
            - user_choice_made: True if user was asked and made a choice (Yes/No)
            - force_all: True if user chose to reprocess all files (not just missing ones)
        """
        # stage_0 (overscan) has its own comprehensive dialog logic in the processing section
        if stage_id == 'stage_0':
            # Always return True to let the processing section handle the dialogs
            return True, False, False
        
        # If stage not marked as completed, definitely need to process
        if not self.processing_state.get(f'{stage_id}_completed', False):
            return True, False, False  # Need to reprocess, no user choice needed, not force_all
        
        # Check if processed files still exist
        if not self._check_processed_files_exist(stage_id):
            # Files missing, need to reprocess, no user choice needed
            return True, False, False
        
        # All files exist - always reprocess (overwrite) without asking
        return True, False, True  # Reprocess, no user choice made, force_all

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

        # Handle calib files
        calib_files = self.calib_file if isinstance(self.calib_file, list) else [self.calib_file]
        calib_files = [f for f in calib_files if f]

        if not calib_files:
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

        # Calculate total images to process
        # Calib file count
        calib_count = 0
        if self.calib_file:
            if isinstance(self.calib_file, list):
                calib_count = len(self.calib_file)
            else:
                calib_count = 1  # Single file as path

        # Total images that will be processed (all images go through overscan)
        total_imgs = len(science_files) + len(self.bias_files) + len(self.flat_files) + calib_count

        # Process each science image
        self.log_text.append("=" * 60)
        self.log_text.append(f"STARTING PIPELINE - Processing {len(science_files)} science image(s)")
        self.log_text.append(f"Total images to be processed: {total_imgs} (bias: {len(self.bias_files)}, flat: {len(self.flat_files)}, calib: {calib_count}, science: {len(science_files)})")
        self.log_text.append("=" * 60)

        self.worker_thread = ProcessingWorker(
            self.pipeline,
            science_files,
            self.bias_files,
            self.flat_files,
            calib_files
        )
        self.worker_thread.progress_updated.connect(self._on_progress)
        self.worker_thread.execution_complete.connect(
            lambda success, msg: self._on_execution_complete(success, msg, len(science_files)-1, len(science_files), total_imgs)
        )
        self.worker_thread.start()

    def _run_selected_pipeline(self):
        """Execute only the selected processing steps."""
        # Reset run-scoped cache to avoid reusing stale paths from previous runs.
        self._flat_files_after_bias = None
        self._flat_files_after_step1 = None
        self._step3_mask = None

        # Check which steps are selected
        selected_stages = []
        for stage_id, checkbox in self.stage_checkboxes.items():
            if checkbox.isChecked():
                selected_stages.append(stage_id)

        step1_selected = "stage_0" in selected_stages
        step1_do_overscan = step1_selected and self.step1_substep_checkboxes["overscan"].isChecked()
        step1_do_bias = step1_selected and self.step1_substep_checkboxes["bias"].isChecked()
        step1_do_cosmic = step1_selected and self.step1_substep_checkboxes["cosmic"].isChecked()

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
            "stage_0": "Step 1: Basic Pre-processing (Overscan/Bias/Cosmic)",
            "stage_1": "Step 2: Order Tracing",
            "stage_2": "Step 3: Scattered Light Subtraction",
            "stage_3": "Step 4: 2D Flat-Field Correction",
            "stage_4": "Step 5: 1D Extraction",
            "stage_5": "Step 6: De-blazing",
            "stage_6": "Step 7: Wavelength Calibration",
        }
        for stage_id in selected_stages:
            self.log_text.append(f"  ✓ {stage_names.get(stage_id, stage_id)}")
        self.log_text.append("=" * 60)

        # Validate required inputs
        if step1_do_bias and not self.bias_files:
            QMessageBox.warning(self, "Error", "Step 1 (Basic Pre-processing) requires bias frames")
            return
        if "stage_1" in selected_stages and not self.flat_files:
            QMessageBox.warning(self, "Error", "Step 2 (Order Tracing) requires flat frames")
            return
        if "stage_5" in selected_stages:
            # Check if we have any calibration files
            has_calib = False
            if isinstance(self.calib_file, list):
                has_calib = len(self.calib_file) > 0
            else:
                has_calib = bool(self.calib_file)  # Could be string or None
            
            if not has_calib:
                QMessageBox.warning(self, "Error", "Step 6 (Wavelength Calibration) requires calibration frames")
                return
        
        # Check if any selected stages need reprocessing
        need_reprocess = False
        user_chose_no = False
        force_reprocess_stages = set()  # Track which stages should force reprocess all files
        
        for stage_id in selected_stages:
            should_reprocess, user_choice_made, force_all = self._should_reprocess_stage(stage_id)
            
            if user_choice_made and not should_reprocess:
                # User was asked and chose NOT to reprocess this stage
                user_chose_no = True
                # User chose NOT to reprocess, cancel the entire operation
                self.statusBar.showMessage("Processing cancelled by user")
                self.log_text.append("\n" + "=" * 60)
                self.log_text.append("PROCESSING CANCELLED")
                self.log_text.append("User chose not to reprocess existing results")
                self.log_text.append("=" * 60)
                self.run_all_btn.setEnabled(True)
                self.run_selected_btn.setEnabled(True)
                self.stop_btn.setEnabled(False)
                return
            
            if should_reprocess:
                need_reprocess = True
                # Reset the stage completion flag
                self.processing_state[f'{stage_id}_completed'] = False
                
                # If user chose to reprocess all files (not just missing ones)
                if force_all:
                    force_reprocess_stages.add(stage_id)
        
        # If no stages need reprocessing (and user wasn't asked or said yes to all)
        if not need_reprocess:
            # This means all stages are already processed and files exist
            # Show a brief status message but don't pop up a dialog
            self.statusBar.showMessage("All selected steps are already processed")
            self.log_text.append("\n" + "=" * 60)
            self.log_text.append("PROCESSING SKIPPED")
            self.log_text.append("All selected steps are already processed")
            self.log_text.append("=" * 60)
            self.run_all_btn.setEnabled(True)
            self.run_selected_btn.setEnabled(True)
            self.stop_btn.setEnabled(False)
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

        # Calculate total images to process
        # Calib file count
        calib_count = 0
        if self.calib_file:
            if isinstance(self.calib_file, list):
                calib_count = len(self.calib_file)
            else:
                calib_count = 1  # Single file as path

        # Total images that will be processed (all images go through overscan if stage_0 is selected)
        # Otherwise, only science images are processed
        if step1_do_overscan:
            total_imgs = len(science_files) + len(self.bias_files) + len(self.flat_files) + calib_count
        else:
            total_imgs = len(science_files)

        # Disable controls
        self.run_all_btn.setEnabled(False)
        self.run_selected_btn.setEnabled(False)
        self.stop_btn.setEnabled(True)

        # --- Data Flow Management for Selected Steps ---
        # Initialize file lists for each type. These will be updated as steps run.
        current_science_files = list(science_files)
        current_flat_files = list(self.flat_files)
        current_bias_files = list(self.bias_files)
        current_calib_files = list(self.calib_file) if isinstance(self.calib_file, list) else [self.calib_file]

        # --- Smart Input Finder for Step 1 ---
        # Before running any Step 1 sub-step, determine the correct input files by
        # checking for existing processed outputs from previous sub-steps.
        output_dir = Path(self.config.get_output_path())
        overscan_dir = output_dir / 'step1_basic' / 'overscan_corrected'
        bias_dir = output_dir / 'step1_basic' / 'bias_subtracted'
        cosmic_dir = output_dir / 'step1_basic' / 'cosmic_corrected'

        def find_latest_input(raw_file_path, check_cosmic=False, check_bias=False, check_overscan=False):
            """Find the most processed version of a file."""
            raw_name = Path(raw_file_path).name
            if check_cosmic:
                cosmic_path = cosmic_dir / raw_name
                if cosmic_path.exists():
                    return str(cosmic_path)
            if check_bias:
                bias_path = bias_dir / raw_name
                if bias_path.exists():
                    return str(bias_path)
            if check_overscan:
                overscan_path = overscan_dir / raw_name
                if overscan_path.exists():
                    return str(overscan_path)
            return raw_file_path

        # Before running any stage, find the correct input files.
        # If a Step 1 sub-step is selected, its inputs are the raw files.
        # Otherwise, find the latest processed version.
        if not step1_do_overscan:
            # If not doing overscan, inputs for bias/cosmic/later stages might be overscan-corrected.
            current_science_files = [find_latest_input(f, check_overscan=True) for f in science_files]
            current_flat_files = [find_latest_input(f, check_overscan=True) for f in self.flat_files]
            current_bias_files = [find_latest_input(f, check_overscan=True) for f in self.bias_files]
            current_calib_files = [find_latest_input(f, check_overscan=True) for f in current_calib_files]
        
        if not step1_do_bias:
            # If not doing bias, inputs for cosmic/later stages might be bias-corrected.
            current_science_files = [find_latest_input(f, check_overscan=True, check_bias=True) for f in science_files]
            current_flat_files = [find_latest_input(f, check_overscan=True, check_bias=True) for f in self.flat_files]
            current_calib_files = [find_latest_input(f, check_overscan=True, check_bias=True) for f in current_calib_files]
        
        if not step1_do_cosmic:
            # If not doing cosmic, inputs for later stages might be cosmic-corrected.
            current_science_files = [find_latest_input(f, check_overscan=True, check_bias=True, check_cosmic=True) for f in science_files]


        if "stage_0" in selected_stages and step1_do_overscan:
            self.log_text.append("\n" + "=" * 60)
            self.log_text.append("STEP 1.1: OVERSCAN CORRECTION")
            self.log_text.append("Processing all image types: Bias, Flat, Calibration, Science")
            self.log_text.append("=" * 60)
            
            try:
                from src.core.basic_reduction import process_overscan_stage
                
                # Collect all files to process in the correct order
                all_files_to_process = []
                all_files_to_process.extend(current_bias_files)
                all_files_to_process.extend(current_flat_files)
                all_files_to_process.extend(current_calib_files)
                all_files_to_process.extend(current_science_files)
                all_files_to_process = sorted(list(set(all_files_to_process))) # Unique list

                overscan_kwargs = {
                    'overscan_start_column': self.config.get_int('data', 'overscan_start_column', -1),
                    'overscan_method': self.config.get('data', 'overscan_method', 'mean_only'),
                    'overscan_smooth_window': self.config.get_int('data', 'overscan_smooth_window', -1),
                    'overscan_poly_order': self.config.get_int('data', 'overscan_poly_order', 3),
                    'overscan_poly_type': self.config.get('data', 'overscan_poly_type', 'legendre'),
                }
                processed_files = process_overscan_stage(
                    all_files_to_process, self.config.get_output_path(), **overscan_kwargs
                )
                
                # Create a map of original filename to processed path
                processed_map = {Path(p).name: p for p in processed_files}

                # Update current file lists to point to the new overscan-corrected files
                current_science_files = [processed_map.get(Path(f).name, f) for f in current_science_files]
                current_flat_files = [processed_map.get(Path(f).name, f) for f in current_flat_files]
                current_bias_files = [processed_map.get(Path(f).name, f) for f in current_bias_files]
                current_calib_files = [processed_map.get(Path(f).name, f) for f in current_calib_files]
                
                self.log_text.append(f"  ✓ Overscan correction complete: {len(processed_files)} images processed")

            except Exception as e:
                self.log_text.append(f"  ✗ Overscan correction failed: {e}")
                raise
        elif "stage_0" in selected_stages:
            self.log_text.append("\n" + "=" * 60)
            self.log_text.append("STEP 1.1: OVERSCAN CORRECTION")
            self.log_text.append("  Skipped (option unchecked)")
            self.log_text.append("=" * 60)

        # --- Step 1.2: Bias Correction ---
        master_bias = None
        if "stage_0" in selected_stages and step1_do_bias and self.bias_files:
            try:
                self.log_text.append("\n============================================================")
                self.log_text.append("STEP 1.2: BIAS SUBTRACTION")
                from src.core.basic_reduction import process_bias_stage, BiasCorrector, read_fits_image, write_fits_image
                
                bias_kwargs = {
                    'combine_method': self.config.get('reduce.bias', 'combine_method', 'median'),
                    'combine_sigma': self.config.get_float('reduce.bias', 'combine_sigma', 3.0),
                    'save_plots': self.config.get_bool('reduce', 'save_plots', True),
                    'fig_format': self.config.get('reduce', 'fig_format', 'png'),
                }
                master_bias, master_bias_path = process_bias_stage(
                    current_bias_files, self.config.get_output_path(), **bias_kwargs
                )
                self._master_bias = master_bias # Cache for later use

                # Apply bias to all other relevant file types
                corrector = BiasCorrector(**bias_kwargs)
                corrector.master_bias = master_bias

                # Apply to science, flat, and calib files
                files_to_correct = {
                    "science": current_science_files,
                    "flat": current_flat_files,
                    "calib": current_calib_files
                }
                output_dir = Path(self.config.get_output_path()) / 'step1_basic' / 'bias_subtracted'
                output_dir.mkdir(parents=True, exist_ok=True)

                for file_type, file_list in files_to_correct.items():
                    if not file_list: continue
                    self.log_text.append(f"  Applying bias correction to {len(file_list)} {file_type} frames...")
                    corrected_list = []
                    for f in file_list:
                        img, hdr = read_fits_image(f)
                        corrected_img = corrector.apply_bias_correction(img)
                        out_path = output_dir / Path(f).name
                        hdr['BIASCOR'] = (True, 'Bias correction applied')

                        # Manually save diagnostic plot for each corrected image
                        if bias_kwargs['save_plots']:
                            from src.plotting.spectra_plotter import plot_2d_image_to_file
                            plot_file = output_dir / f'{Path(f).stem}_bias_corrected.{bias_kwargs["fig_format"]}'
                            plot_2d_image_to_file(corrected_img, str(plot_file),
                                                  f"Bias Corrected {file_type.capitalize()} Image")

                        write_fits_image(str(out_path), corrected_img, header=hdr, dtype='float32')
                        corrected_list.append(str(out_path))
                    
                    # Update the current file list for this type
                    if file_type == "science": current_science_files = corrected_list
                    elif file_type == "flat": current_flat_files = corrected_list
                    elif file_type == "calib": current_calib_files = corrected_list

                self.log_text.append(f"  ✓ Bias correction application complete.")

            except Exception as e:
                self.log_text.append(f"  ✗ Bias correction failed: {e}")
                raise
        elif "stage_0" in selected_stages:
                self.log_text.append("\n============================================================")
                self.log_text.append("STEP 1.2: BIAS SUBTRACTION")
                self.log_text.append("  Skipped (option unchecked or no bias files)")
                self.log_text.append("============================================================")

            # Apply cosmic-ray removal in Step1 basic preprocessing (science only).
        if "stage_0" in selected_stages and step1_do_cosmic:
                if self.config.get_bool('reduce', 'cosmic_enabled', True):
                    from src.core.basic_reduction import process_cosmic_stage

                    self.log_text.append("\n============================================================")
                    self.log_text.append("STEP 1.3: COSMIC RAY REMOVAL")
                    self.log_text.append("============================================================")
                    self.log_text.append("  Applying cosmic ray removal to SCIENCE frames only...")

                    cosmic_kwargs = {
                        'cosmic_sigclip': self.config.get_float('reduce', 'cosmic_sigclip', 5.0),
                        'cosmic_objlim': self.config.get_float('reduce', 'cosmic_objlim', 5.0),
                        'cosmic_gain': self.config.get_float('reduce', 'cosmic_gain', 1.0),
                        'cosmic_readnoise': self.config.get_float('reduce', 'cosmic_readnoise', 5.0),
                    }

                    if current_science_files:
                        science_cosmic_outputs = process_cosmic_stage(
                            current_science_files,
                            output_dir_base=self.config.get_output_path(),
                            **cosmic_kwargs
                        )
                        current_science_files = science_cosmic_outputs
                        self.log_text.append("  ✓ Cosmic ray removal complete.")
                    else:
                        self.log_text.append("  No science files to process for cosmic rays.")
                else:
                    self.log_text.append("\n============================================================")
                    self.log_text.append("STEP 1.3: COSMIC RAY REMOVAL")
                    self.log_text.append("============================================================")
                    self.log_text.append("  Cosmic ray removal disabled by config")
        elif "stage_0" in selected_stages:
            self.log_text.append("\n============================================================")
            self.log_text.append("STEP 1.3: COSMIC RAY REMOVAL")
            self.log_text.append("============================================================")
            self.log_text.append("  Skipped (option unchecked)")

        # Update processing state for Step 1
        if "stage_0" in selected_stages:
            self.processing_state['stage_1_completed'] = True
            if master_bias is not None:
                self.processing_state['master_bias_path'] = str(master_bias_path)

            self.log_text.append("\n============================================================")
            self.log_text.append("Step 1 pre-processing complete.")
            self.log_text.append(f"  Output directory: {Path(self.config.get_output_path()) / 'step1_basic'}")
            self.log_text.append(f"  Science files ready for next steps: {len(current_science_files)}")
            self.log_text.append(f"  Flat files ready for next steps: {len(current_flat_files)}")
            self.log_text.append(f"  Calibration files ready for next steps: {len(current_calib_files)}")
            self.log_text.append("============================================================")

        if "stage_1" in selected_stages and self.flat_files:
            self.log_text.append("\n============================================================")
            self.log_text.append("STEP 2: ORDER TRACING")
            from src.core.order_tracing import process_order_tracing_stage
            self.log_text.append(f"  Using {len(current_flat_files)} flat frames as input...")

            trace_kwargs = {
                'output_dir_base': self.config.get_output_path(),
                'combine_method': self.config.get('reduce.flat', 'combine_method', 'median'),
                'combine_sigma': self.config.get_float('reduce.bias', 'combine_sigma', 3.0),
                'mosaic_maxcount': self.config.get_float('reduce.flat', 'mosaic_maxcount', 65535),
                'snr_threshold': self.config.get_float('reduce.trace', 'snr_threshold', 5.0),
                'gap_fill_factor': self.config.get_float('reduce.trace', 'gap_fill_factor', 1.6),
                'min_trace_coverage': self.config.get_float('reduce.trace', 'min_trace_coverage', 0.20),
                'trace_degree': self.config.get_int('reduce.trace', 'degree', 4),
                'width_cheb_degree': self.config.get_int('reduce.trace', 'width_cheb_degree', 3),
                'aperture_boundary_snr': self.config.get_float('reduce.trace', 'aperture_boundary_snr', 3.0),
                'save_plots': self.config.get_bool('reduce', 'save_plots', True),
                'fig_format': self.config.get('reduce', 'fig_format', 'png'),
            }

            flat_field, apertures = process_order_tracing_stage(
                current_flat_files, **trace_kwargs
            )
            self.log_text.append(f"  ✓ Order tracing complete: {apertures.norders if apertures else 0} orders detected.")
            self._flat_field = flat_field
            self._apertures = apertures
            self.processing_state['stage_2_completed'] = True

        # Attempt to load flat_field and apertures if needed by downstream stages
        needs_flat = any(s in selected_stages for s in ["stage_2", "stage_3", "stage_4", "stage_5", "stage_6"])
        if needs_flat:
            flat_field = getattr(self, '_flat_field', getattr(self.pipeline.state, 'flat_field', None))
            apertures = getattr(self, '_apertures', getattr(self.pipeline.state, 'apertures', None))
            
            if flat_field is None or apertures is None:
                self.log_text.append("  Loading MasterFlat and Apertures from disk...")
                try:
                    from src.utils.fits_io import read_fits_image
                    from src.core.data_structures import ApertureSet, ApertureLocation, FlatField
                    
                    output_dir = Path(self.config.get_output_path())
                    flat_path_step3 = output_dir / 'step3_scatterlight' / 'MasterFlat.fits'
                    flat_path_step2 = output_dir / 'step2_trace' / 'MasterFlat.fits'
                    flat_path = flat_path_step3 if flat_path_step3.exists() else flat_path_step2
                    coefs_path = output_dir / 'step2_trace' / 'Orders_trace_coefs.json'
                    
                    if flat_path.exists() and coefs_path.exists():
                        flat_img, _ = read_fits_image(str(flat_path))
                        
                        with open(coefs_path, 'r') as f:
                            coefs_data = json.load(f)
                            
                        loaded_apertures = ApertureSet()
                        width = flat_img.shape[1]
                        for ap_id_str, ap_data in coefs_data.get('orders', {}).items():
                            ap = ApertureLocation(
                                aperture=int(ap_id_str), order=int(ap_id_str),
                                is_chebyshev=True, domain=(0.0, float(width - 1))
                            )
                            if 'center_arr' in ap_data:
                                ap.center_arr = np.array(ap_data['center_arr'])
                            if 'w_up_cheb' in ap_data:
                                ap.w_up_cheb_coef = np.array(ap_data['w_up_cheb'])
                            if 'w_low_cheb' in ap_data:
                                ap.w_low_cheb_coef = np.array(ap_data['w_low_cheb'])
                            loaded_apertures.add_aperture(ap)
                            
                        flat_field = FlatField(
                            flat_data=flat_img, flat_mask=None, flat_norm=None, 
                            flat_sens=None, scattered_light=None, smoothed_model=None, 
                            pixel_flat=None, illumination_flat=None, aperture_set=loaded_apertures
                        )
                        apertures = loaded_apertures
                        
                        # Load blaze profiles from Step 4
                        blaze_pkl = output_dir / 'step4_flat_corrected' / 'blaze_profiles.pkl'
                        if blaze_pkl.exists():
                            try:
                                import pickle
                                with open(blaze_pkl, 'rb') as fb:
                                    flat_field.blaze_profiles = pickle.load(fb)
                            except Exception as e:
                                self.log_text.append(f"  ! Failed to load blaze_profiles.pkl: {e}")

                        # Load 2D smoothed model from Step 4 for True Real-Space Optimal Extraction
                        model_2d_path = output_dir / 'step4_flat_corrected' / 'model_flat_2d.fits'
                        if model_2d_path.exists():
                            try:
                                model_img, _ = read_fits_image(str(model_2d_path))
                                flat_field.smoothed_model = model_img
                            except Exception as e:
                                self.log_text.append(f"  ! Failed to load model_flat_2d.fits: {e}")

                        self._flat_field = flat_field
                        self._apertures = apertures
                        self.pipeline.state.flat_field = flat_field
                        self.pipeline.state.apertures = apertures
                        self.log_text.append(f"  ✓ Successfully loaded {flat_path.name} and Apertures from disk.")
                    else:
                        self.log_text.append("  ! Could not find MasterFlat.fits or Orders_trace_coefs.json on disk. Run Step 2 first.")
                except Exception as e:
                    self.log_text.append(f"  ✗ Failed to load from disk: {e}")

        if "stage_2" in selected_stages:
            flat_field = getattr(self, '_flat_field', getattr(self.pipeline.state, 'flat_field', None))
            apertures = getattr(self, '_apertures', getattr(self.pipeline.state, 'apertures', None))
            if flat_field is not None and apertures is not None:
                self.log_text.append("\n============================================================")
                self.log_text.append("STEP 3: SCATTERED LIGHT SUBTRACTION (MASTER FLAT)")
                try:
                    from src.core.scattered_light import create_widened_mask
                    from src.core.scattered_light import process_background_stage
                    from src.utils.fits_io import write_fits_image
                    
                    out_dir = Path(self.config.get_output_path()) / 'step3_scatterlight'
                    out_dir.mkdir(parents=True, exist_ok=True)
                    
                    mask_margin_pixels = self.config.get_int('reduce.background', 'mask_margin_pixels', 1)
                    n_mask_below = self.config.get_int('reduce.trace', 'n_mask_below', 4)
                    n_mask_above = self.config.get_int('reduce.trace', 'n_mask_above', 4)
                    h, w = flat_field.flat_data.shape
                    
                    step3_mask, lo_traces, hi_traces, full_ids = create_widened_mask(
                        apertures, flat_field.flat_data, mask_margin_pixels,
                        n_mask_below, n_mask_above
                    )
                    
                    write_fits_image(str(out_dir / 'Scattered_Light_Mask.fits'), step3_mask, dtype='uint8')
                    self._step3_mask = step3_mask
                    
                    if self.config.get_bool('reduce', 'save_plots', True):
                        import matplotlib.pyplot as plt
                        fig_format = self.config.get('reduce', 'fig_format', 'png')
                        plt.figure(figsize=(12, 10))
                        plt.imshow(step3_mask, aspect='auto', origin='lower', cmap='gray', vmin=0, vmax=1)
                        
                        x_anno = w - 1
                        for i, ap_id in enumerate(full_ids):
                            yc = 0.5 * (lo_traces[i, x_anno] + hi_traces[i, x_anno])
                            if np.isfinite(yc) and 0 <= yc < h:
                                plt.text(x_anno + 5, yc, str(ap_id), color='red', fontsize=3.5, ha='left', va='center', clip_on=False)
                        
                        plt.title(f'Step 3 Widened Mask (margin={mask_margin_pixels})\nWhite=Masked (Orders + Virtual), Black=Background')
                        plt.xlabel('Pixel (X)')
                        plt.ylabel('Pixel (Y)')
                        plt.xlim(0, w - 1)
                        plt.ylim(0, h - 1)
                        plot_file = out_dir / f'scattered_light_mask.{fig_format}'
                        plt.savefig(str(plot_file), dpi=150, bbox_inches='tight')
                        plt.close()
                        self.log_text.append(f"  ✓ Saved widened order mask and plot (margin={mask_margin_pixels})")

                    bg_kwargs = {
                        'output_subdir': 'step3_scatterlight',
                        'output_tag': 'MasterFlat_',
                        'mask_margin_scale': 1.0,
                        'mask_margin_pixels': 0, # Already widened
                        'order_mask': step3_mask,
                        'poly_order': self.config.get_int('reduce.background', 'poly_order', 3),
                        'bg_method': self.config.get('reduce.background', 'method', 'convolution'),
                        'sigma_clip_val': self.config.get_float('reduce.background', 'sigma_clip', 3.0),
                        'maxiters': self.config.get_int('reduce.background', 'sigma_clip_maxiters', 4),
                        'bspline_smooth': self.config.get_float('reduce.background', 'bspline_smooth', 1.0),
                        'n_mask_below': self.config.get_int('reduce.trace', 'n_mask_below', 4),
                        'n_mask_above': self.config.get_int('reduce.trace', 'n_mask_above', 4),
                        'clip_mode': self.config.get('reduce.background', 'sigma_clip_mode', 'upper'),
                        'split_row': self.config.get_int('data', 'detector_split_row', 2068),
                        'kernel_sigma_x': self.config.get_float('reduce.background', 'kernel_sigma_x', 13.0),
                        'kernel_sigma_y': self.config.get_float('reduce.background', 'kernel_sigma_y', 13.0),
                        'spline_smooth_factor': self.config.get_float('reduce.background', 'spline_smooth_factor', 1.0),
                        'spline_post_smooth_x': self.config.get_float('reduce.background', 'spline_post_smooth_x', 5.0),
                        'save_plots': self.config.get_bool('reduce', 'save_plots', True),
                        'fig_format': self.config.get('reduce', 'fig_format', 'png')
                    }
                    
                    self.log_text.append(f"  Subtracting scattered light from Master Flat...")
                    flat_background = process_background_stage(
                        science_image=flat_field.flat_data,
                        output_dir_base=self.config.get_output_path(),
                        apertures=apertures,
                        **bg_kwargs
                    )
                    flat_clean = np.clip(
                        flat_field.flat_data.astype(np.float32) - flat_background.astype(np.float32),
                        1e-6,
                        None
                    )
                    flat_field.scattered_light = flat_background
                    self._flat_field = flat_field
                    self.pipeline.state.flat_field = flat_field
                    
                    master_flat_path = Path(self.config.get_output_path()) / 'step2_trace' / 'MasterFlat.fits'
                    from src.utils.fits_io import read_fits_image
                    _, flat_header = read_fits_image(str(master_flat_path)) if master_flat_path.exists() else (None, None)
                    if flat_header is not None:
                        flat_header['BKGSCAT'] = (True, 'Scattered light background subtracted')
                        
                    write_fits_image(str(out_dir / 'MasterFlat.fits'), flat_clean.astype('float32'), header=flat_header, dtype='float32')
                    self.log_text.append(f"  ✓ Master Flat scattered light subtraction complete.")
                    self.processing_state['stage_3_completed'] = True
                    
                except Exception as e:
                    self.log_text.append(f"  ✗ Master Flat scattered light subtraction failed: {e}")
                    raise
            else:
                self.log_text.append("\n  ! Step 3 (Master Flat): No flat field or apertures available. Run Step 2 first.")

        total_stages = len(selected_stages)
        current_stage = 0
        if "stage_0" in selected_stages: current_stage += 1
        if "stage_1" in selected_stages: current_stage += 1

        extracted_spectra_dict = {}
        deblazed_spectra_dict = {}
        calibrated_spectra_dict = {}

        apertures = getattr(self, '_apertures', getattr(self.pipeline.state, 'apertures', None))
        flat_field = getattr(self, '_flat_field', getattr(self.pipeline.state, 'flat_field', None))

        # STAGE 2: Scattered Light Subtraction (Science Images)
        if "stage_2" in selected_stages:
            current_stage += 1
            progress = (current_stage / total_stages) * 100
            self.progress_bar.setValue(int(progress))
            self.stage_label.setText("Scattered Light Subtraction")
            
            self.log_text.append("\n" + "=" * 60)
            self.log_text.append("STEP 3: SCATTERED LIGHT SUBTRACTION (Science Images)")
            self.log_text.append("=" * 60)
            
            from src.core.scattered_light import process_background_stage
            from src.utils.fits_io import read_fits_image, write_fits_image
            
            new_science_files = []
            for sci_file in current_science_files:
                try:
                    sci_img, sci_header = read_fits_image(sci_file)
                    sci_name = Path(sci_file).stem
                    order_mask = getattr(self, '_step3_mask', None)
                    if order_mask is None:
                        mask_path = Path(self.config.get_output_path()) / 'step3_scatterlight' / 'Scattered_Light_Mask.fits'
                        if mask_path.exists():
                            mask_img, _ = read_fits_image(str(mask_path))
                            order_mask = mask_img
                            
                    background = process_background_stage(
                        science_image=sci_img,
                        output_dir_base=self.config.get_output_path(),
                        output_subdir='step3_scatterlight',
                        output_tag=f'{sci_name}_',
                        apertures=apertures,
                        mask_margin_scale=1.0,
                        sci_header=sci_header,
                        order_mask=order_mask,
                        poly_order=self.config.get_int('reduce.background', 'poly_order', 3),
                        bg_method=self.config.get('reduce.background', 'method', 'convolution'),
                        sigma_clip_val=self.config.get_float('reduce.background', 'sigma_clip', 3.0),
                        maxiters=self.config.get_int('reduce.background', 'sigma_clip_maxiters', 4),
                        bspline_smooth=self.config.get_float('reduce.background', 'bspline_smooth', 1.0),
                        mask_margin_pixels=0,
                        n_mask_below=self.config.get_int('reduce.trace', 'n_mask_below', 4),
                        n_mask_above=self.config.get_int('reduce.trace', 'n_mask_above', 4),
                        clip_mode=self.config.get('reduce.background', 'sigma_clip_mode', 'upper'),
                        split_row=self.config.get_int('data', 'detector_split_row', 2068),
                        kernel_sigma_x=self.config.get_float('reduce.background', 'kernel_sigma_x', 13.0),
                        kernel_sigma_y=self.config.get_float('reduce.background', 'kernel_sigma_y', 13.0),
                        spline_smooth_factor=self.config.get_float('reduce.background', 'spline_smooth_factor', 1.0),
                        spline_post_smooth_x=self.config.get_float('reduce.background', 'spline_post_smooth_x', 5.0),
                        save_plots=self.config.get_bool('reduce', 'save_plots', True),
                        fig_format=self.config.get('reduce', 'fig_format', 'png')
                    )
                    corrected = sci_img.astype('float32') - background.astype('float32')
                    out_dir = Path(self.config.get_output_path()) / 'step3_scatterlight'
                    out_dir.mkdir(parents=True, exist_ok=True)
                    corrected_path = out_dir / Path(sci_file).name
                    if sci_header is not None:
                        sci_header['BKGSCAT'] = (True, 'Scattered light background subtracted')
                    write_fits_image(str(corrected_path), corrected, header=sci_header, dtype='float32')
                    self.log_text.append(f"  ✓ Step 3 output for {sci_name}: {corrected_path.name}")
                    new_science_files.append(str(corrected_path))
                except Exception as e:
                    self.log_text.append(f"  ✗ Step 3 failed for {Path(sci_file).name}: {e}")
                    new_science_files.append(sci_file)
            current_science_files = new_science_files

        # STAGE 3: 2D Flat-Field Correction
        if "stage_3" in selected_stages:
            current_stage += 1
            progress = (current_stage / total_stages) * 100
            self.progress_bar.setValue(int(progress))
            self.stage_label.setText("2D Flat-Field Correction")
            
            self.log_text.append("\n" + "=" * 60)
            self.log_text.append("STEP 4: 2D FLAT-FIELD CORRECTION")
            self.log_text.append("=" * 60)
            
            from src.core.flat_correction import process_flat_correction_stage
            from src.utils.fits_io import read_fits_image, write_fits_image

            calib_files = self.calib_file if isinstance(self.calib_file, list) else [self.calib_file]
            new_calib_files = []
            if calib_files and flat_field is not None:
                for calib_f in calib_files:
                    if not calib_f: continue
                    try:
                        calib_corr_path = self.pipeline.stage_flat_fielding_calib_2d(calib_f)
                        self.log_text.append(f"  ✓ Step 4 calib output: {Path(calib_corr_path).name}")
                        new_calib_files.append(calib_corr_path)
                    except Exception as e:
                        self.log_text.append(f"  ! Step 4 calib correction skipped/failed: {e}")
                        new_calib_files.append(calib_f)
                current_calib_files = new_calib_files
            
            if flat_field is None:
                self.log_text.append("  ! Step 4: no MasterFlat available – Step 2 must run first")
            else:
                new_science_files = []
                for sci_file in current_science_files:
                    try:
                        sci_img, sci_header = read_fits_image(sci_file)
                        sci_name = Path(sci_file).stem
                        flat_corr_kwargs = {
                            'blaze_smooth_factor': self.config.get_float('reduce.flat', 'blaze_smooth_factor', 1.0),
                            'width_smooth_window': self.config.get_int('reduce.flat', 'width_smooth_window', 41),
                            'profile_bin_step': self.config.get_float('reduce.flat', 'profile_bin_step', 0.01),
                            'n_profile_segments': self.config.get_int('reduce.flat', 'n_profile_segments', 100),
                            'profile_smooth_sigma': self.config.get_float('reduce.flat', 'profile_smooth_sigma', 6.0),
                            'pixel_flat_min': self.config.get_float('reduce.flat', 'pixel_flat_min', 0.5),
                            'pixel_flat_max': self.config.get_float('reduce.flat', 'pixel_flat_max', 1.5),
                            'fringe_orders': self.config.get_int('reduce.flat', 'fringe_orders', 20),
                            'save_plots': self.config.get_bool('reduce', 'save_plots', True),
                            'fig_format': self.config.get('reduce', 'fig_format', 'png'),
                        }
                        corrected = process_flat_correction_stage(
                            science_image=sci_img,
                            flat_field=flat_field,
                            output_dir_base=self.config.get_output_path(),
                            apertures=apertures,
                            science_name=sci_name,
                            **flat_corr_kwargs
                        )
                        out_dir = Path(self.config.get_output_path()) / 'step4_flat_corrected'
                        corrected_path = out_dir / Path(sci_file).name
                        if sci_header is not None:
                            sci_header['FLATCOR'] = (True, '2D pixel flat correction applied')
                            write_fits_image(str(corrected_path), corrected, header=sci_header, dtype='float32')
                        self.log_text.append(f"  ✓ Step 4 output for {sci_name}: {corrected_path.name}")
                        new_science_files.append(str(corrected_path))
                    except Exception as e:
                        self.log_text.append(f"  ✗ Step 4 failed for {Path(sci_file).name}: {e}")
                        new_science_files.append(sci_file)
                current_science_files = new_science_files
                self._flat_field = flat_field
                self.pipeline.state.flat_field = flat_field

        # STAGE 4: 1D Extraction
        if "stage_4" in selected_stages:
            current_stage += 1
            progress = (current_stage / total_stages) * 100
            self.progress_bar.setValue(int(progress))
            self.stage_label.setText("1D Extraction")
            
            if apertures is None:
                self.log_text.append("  ! Step 5: no apertures available – Step 2 must run first")
            else:
                self.log_text.append("\n" + "=" * 60)
                self.log_text.append("STEP 5: 1D EXTRACTION")
                self.log_text.append("=" * 60)
                
                from src.core.extraction import process_extraction_stage
                from src.utils.fits_io import read_fits_image
                extract_kwargs = {
                    'output_dir_base': self.config.get_output_path(),
                    'optimal_sigma': self.config.get_float('reduce.extract', 'optimal_sigma', 3.0),
                    'save_plots': self.config.get_bool('reduce', 'save_plots', True),
                    'fig_format': self.config.get('reduce', 'fig_format', 'png'),
                }

                # Extract MasterFlat
                mf_path = Path(self.config.get_output_path()) / 'step4_flat_corrected' / 'MasterFlat.fits'
                if mf_path.exists():
                    try:
                        mf_img, _ = read_fits_image(str(mf_path))
                        mf_name = "MasterFlat"
                        process_extraction_stage(mf_img, apertures, wavelength_calib=None, flat_field=flat_field, method_override='sum', output_filename=f'{mf_name}_1D_sum.fits', plot_prefix=f'{mf_name}_1D_sum', **extract_kwargs)
                        process_extraction_stage(mf_img, apertures, wavelength_calib=None, flat_field=flat_field, method_override='optimal', output_filename=f'{mf_name}_1D_optimal.fits', plot_prefix=f'{mf_name}_1D_optimal', **extract_kwargs)
                        self.log_text.append(f"  ✓ Extracted MasterFlat 1D spectra (sum & optimal)")
                    except Exception as e:
                        self.log_text.append(f"  ! Failed to extract MasterFlat: {e}")

                # Extract Calib
                for calib_f in current_calib_files:
                    if not calib_f: continue
                    best_calib = self.pipeline.get_best_calib_file(calib_f)
                    if best_calib and Path(best_calib).exists():
                        try:
                            c_img, _ = read_fits_image(best_calib)
                            c_name = Path(best_calib).stem
                            process_extraction_stage(c_img, apertures, wavelength_calib=None, flat_field=flat_field, method_override='sum', output_filename=f'{c_name}_1D_sum.fits', plot_prefix=f'{c_name}_1D_sum', **extract_kwargs)
                            process_extraction_stage(c_img, apertures, wavelength_calib=None, flat_field=flat_field, method_override='optimal', output_filename=f'{c_name}_1D_optimal.fits', plot_prefix=f'{c_name}_1D_optimal', **extract_kwargs)
                            self.log_text.append(f"  ✓ Extracted Calibration 1D spectra for {c_name} (sum & optimal)")
                        except Exception as e:
                            self.log_text.append(f"  ! Failed to extract Calibration {c_name}: {e}")

                # Extract Science
                for sci_file in current_science_files:
                    try:
                        sci_img, _ = read_fits_image(sci_file)
                        sci_name = Path(sci_file).stem
                        process_extraction_stage(sci_img, apertures, wavelength_calib=None, flat_field=flat_field, method_override='sum', output_filename=f'{sci_name}_1D_sum.fits', plot_prefix=f'{sci_name}_1D_sum', **extract_kwargs)
                        spectra_opt = process_extraction_stage(sci_img, apertures, wavelength_calib=None, flat_field=flat_field, method_override='optimal', output_filename=f'{sci_name}_1D_optimal.fits', plot_prefix=f'{sci_name}_1D_optimal', **extract_kwargs)
                        extracted_spectra_dict[sci_name] = spectra_opt
                        self.log_text.append(f"  ✓ Extracted Science 1D spectra for {sci_name} (sum & optimal)")
                    except Exception as e:
                        self.log_text.append(f"  ! Failed to extract Science {Path(sci_file).name}: {e}")

        # STAGE 5: De-blazing
        if "stage_5" in selected_stages:
            current_stage += 1
            progress = (current_stage / total_stages) * 100
            self.progress_bar.setValue(int(progress))
            self.stage_label.setText("De-blazing")
            
            if flat_field is None or not flat_field.blaze_profiles:
                self.log_text.append("  ! Step 6: no blaze function available – Step 4 must run first")
            else:
                self.log_text.append("\n" + "=" * 60)
                self.log_text.append("STEP 6: DE-BLAZING")
                self.log_text.append("=" * 60)
                
                from src.core.de_blazing import process_de_blazing_stage
                from src.core.extraction import load_extracted_spectra
                
                def do_deblaze(target_name, spectra_obj=None):
                    deblazed_opt = None
                    for method in ['sum', 'optimal']:
                        spectra = spectra_obj if (method == 'optimal' and spectra_obj is not None) else None
                        if spectra is None:
                            spectra_path = Path(self.config.get_output_path()) / 'step5_extraction' / f'{target_name}_1D_{method}.fits'
                            if spectra_path.exists():
                                try: spectra = load_extracted_spectra(str(spectra_path))
                                except Exception: pass
                        if spectra is not None:
                            try:
                                db_spec = process_de_blazing_stage(
                                    spectra,
                                    output_dir_base=self.config.get_output_path(),
                                    flat_field=flat_field,
                                    save_deblaze=self.config.get_bool('reduce.save_intermediate', 'save_deblaze', True),
                                    output_filename=f'{target_name}_1D_{method}_Deblaze.fits',
                                    plot_prefix=f'{target_name}_1D_{method}_Deblaze',
                                    save_plots=self.config.get_bool('reduce', 'save_plots', True),
                                    fig_format=self.config.get('reduce', 'fig_format', 'png')
                                )
                                self.log_text.append(f"  ✓ Step 6 De-blazing complete for {target_name} ({method})")
                                if method == 'optimal': deblazed_opt = db_spec
                            except Exception as e:
                                self.log_text.append(f"  ! Failed to de-blaze {target_name} ({method}): {e}")
                    return deblazed_opt

                do_deblaze("MasterFlat")
                for calib_f in current_calib_files:
                    if not calib_f: continue
                    calib_name = Path(self.pipeline.get_best_calib_file(calib_f)).stem
                    do_deblaze(calib_name)
                    
                for sci_file in current_science_files:
                    sci_name = Path(sci_file).stem
                    db_opt = do_deblaze(sci_name, extracted_spectra_dict.get(sci_name))
                    if db_opt is not None:
                        deblazed_spectra_dict[sci_name] = db_opt

        # STAGE 6: Wavelength Calibration
        if "stage_6" in selected_stages:
            current_stage += 1
            progress = (current_stage / total_stages) * 100
            self.progress_bar.setValue(int(progress))
            self.stage_label.setText("Wavelength Calibration")
            
            if not current_calib_files or not current_calib_files[0]:
                self.log_text.append("  ! Step 7: no calibration file provided")
            else:
                 self.log_text.append("\n" + "=" * 60)
                 self.log_text.append("STEP 7: WAVELENGTH CALIBRATION")
                 self.log_text.append("=" * 60)
                 
                 from src.core.wave_calibration import process_wavelength_stage, WavelengthCalibrator
                 from src.utils.fits_io import read_fits_image
                 from src.core.extraction import process_extraction_stage, load_extracted_spectra
                 import os
                 from astropy.time import Time
                 
                 time_key = self.config.get('data', 'statime_key', 'DATE-OBS')
                 def get_file_time(filepath):
                     if not filepath or not Path(filepath).exists(): return 0
                     try:
                         _, hdr = read_fits_image(filepath)
                         if hdr and time_key in hdr:
                             return Time(hdr[time_key]).unix
                     except Exception: pass
                     return os.path.getmtime(filepath)

                 calibrations = []
                 for calib_file in current_calib_files:
                     if not calib_file: continue
                     best_calib = self.pipeline.get_best_calib_file(calib_file)
                     calib_name = Path(best_calib).stem
                     lamp_spectra = None
                     
                     lamp_db_path = Path(self.config.get_output_path()) / 'step6_deblazing' / f'{calib_name}_1D_optimal_Deblaze.fits'
                     if lamp_db_path.exists():
                         try: lamp_spectra = load_extracted_spectra(str(lamp_db_path))
                         except Exception: pass
                         
                     if lamp_spectra is None:
                         lamp_ex_path = Path(self.config.get_output_path()) / 'step5_extraction' / f'{calib_name}_1D_optimal.fits'
                         if lamp_ex_path.exists():
                             try: lamp_spectra = load_extracted_spectra(str(lamp_ex_path))
                             except Exception: pass
                             
                     if lamp_spectra is None:
                         self.log_text.append(f"  ! Step 7: Calibration lamp 1D optimal spectrum not found for {calib_name}. Please run Step 5 or 6 first.")
                         continue
                         
                     wave_calib = process_wavelength_stage(
                         lamp_spectra=lamp_spectra, config=self.config, output_dir_base=self.config.get_output_path(),
                         lamp_type=self.config.get('telescope.linelist', 'linelist_type', 'ThAr'),
                         save_plots=self.config.get_bool('reduce', 'save_plots', True),
                         fig_format=self.config.get('reduce', 'fig_format', 'png'),
                     )
                     c_time = get_file_time(best_calib)
                     calibrations.append((c_time, wave_calib, calib_name))
                 
                 if not calibrations:
                     self.log_text.append("  ! Step 7: No valid wavelength calibrations could be generated.")
                 else:
                     targets = [("MasterFlat", None)]
                     for c_file in current_calib_files:
                         if not c_file: continue
                         best_calib = self.pipeline.get_best_calib_file(c_file)
                         targets.append((Path(best_calib).stem, best_calib))
                     for s_file in current_science_files:
                         targets.append((Path(s_file).stem, s_file))
                         
                     for target_name, target_file in targets:
                         target_time = get_file_time(target_file)
                         closest_calib = min(calibrations, key=lambda x: abs(x[0] - target_time))
                         wave_calib = closest_calib[1]
                         self.log_text.append(f"  Using calibration {closest_calib[2]} for {target_name} (closest in time)")
                         
                         calib_opt = None
                         for method in ['sum', 'optimal']:
                             spectra = None
                             deblaze_path = Path(self.config.get_output_path()) / 'step6_deblazing' / f'{target_name}_1D_{method}_Deblaze.fits'
                             if deblaze_path.exists():
                                 try: spectra = load_extracted_spectra(str(deblaze_path))
                                 except Exception: pass
                             
                             if spectra is None:
                                 extract_path = Path(self.config.get_output_path()) / 'step5_extraction' / f'{target_name}_1D_{method}.fits'
                                 if extract_path.exists():
                                     try: spectra = load_extracted_spectra(str(extract_path))
                                     except Exception: pass
                             
                             if spectra is None:
                                 self.log_text.append(f"  ! Step 7: no {method} spectra to calibrate for {target_name}")
                                 continue
                                 
                             delta_m = getattr(wave_calib, 'delta_m', 0)
                             if delta_m != 0:
                                 spectra.shift_orders(delta_m)
                                 
                             calibrator = WavelengthCalibrator()
                             calibrator.wave_calib = wave_calib
                             
                             from src.core.data_structures import SpectraSet
                             calibrated_spectra = SpectraSet()
                             
                             for spectrum in spectra.spectra.values():
                                 n_pixels = len(spectrum.flux)
                                 pixel_array = np.arange(n_pixels)
                                 wavelengths = calibrator.apply_wavelength_calibration(
                                     spectrum.flux, pixel_array, aperture_y=float(spectrum.aperture)
                                 )
                                 calibrated_spectrum = spectrum.copy()
                                 calibrated_spectrum.wavelength = wavelengths
                                 calibrated_spectra.add_spectrum(calibrated_spectrum)
                                 
                             out_dir = Path(self.config.get_output_path()) / 'step7_wavelength'
                             out_dir.mkdir(parents=True, exist_ok=True)
                             from src.core.de_blazing import save_deblazed_spectra
                             save_deblazed_spectra(str(out_dir / f'{target_name}_1D_{method}_calibrated.fits'), calibrated_spectra)
                             
                             self.log_text.append(f"  ✓ Step 7 Wavelength calibration applied for {target_name} ({method})")
                             if method == 'optimal':
                                 calib_opt = calibrated_spectra
                                 
                         if calib_opt is not None and target_name in [Path(s).stem for s in current_science_files]:
                             calibrated_spectra_dict[target_name] = calib_opt

        # STAGE 7: Order Stitching
        if "stage_7" in selected_stages:
            current_stage += 1
            progress = (current_stage / total_stages) * 100
            self.progress_bar.setValue(int(progress))
            self.stage_label.setText("Order Stitching")
            
            self.log_text.append("\n" + "=" * 60)
            self.log_text.append("STEP 8: ORDER STITCHING")
            self.log_text.append("=" * 60)
            
            from src.core.order_stitching import process_order_stitching_stage
            from src.core.extraction import load_extracted_spectra
            
            for sci_file in current_science_files:
                sci_name = Path(sci_file).stem
                spectra = calibrated_spectra_dict.get(sci_name)
                if spectra is None:
                    calibrated_path = Path(self.config.get_output_path()) / 'step7_wavelength' / f'{sci_name}_1D_calibrated.fits'
                    if calibrated_path.exists():
                        try: spectra = load_extracted_spectra(str(calibrated_path))
                        except Exception: pass
                        
                if spectra is None:
                    self.log_text.append(f"  ! Step 8: no calibrated spectra available for {sci_name}")
                else:
                    try:
                        stitched = process_order_stitching_stage(
                            spectra,
                            output_dir_base=self.config.get_output_path(),
                            output_subdir=f'step8_stitching/{sci_name}',
                            save_plots=self.config.get_bool('reduce', 'save_plots', True),
                            fig_format=self.config.get('reduce', 'fig_format', 'png'),
                        )
                        self.log_text.append(f"  ✓ Step 8 Order Stitching complete for {sci_name}")
                    except Exception as e:
                        self.log_text.append(f"  ✗ Step 8 Order Stitching failed for {sci_name}: {e}")

        self.log_text.append("\n============================================================")
        self.log_text.append(f"✓ {len(current_science_files)} science images processed through selected steps.")
        self.log_text.append("============================================================")
        
        self.progress_bar.setValue(100)
        self.stage_label.setText("Completed")
        self.run_all_btn.setEnabled(True)
        self.run_selected_btn.setEnabled(True)
        self.stop_btn.setEnabled(False)

    def _select_all_stages(self):
        """Select all processing stages."""
        for checkbox in self.stage_checkboxes.values():
            checkbox.setChecked(True)
        self.statusBar.showMessage("All processing steps selected")

    def _clear_all_stages(self):
        """Clear all processing stage selections."""
        for checkbox in self.stage_checkboxes.values():
            checkbox.setChecked(False)
        self.statusBar.showMessage("All processing steps cleared")

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
                             img_idx: int = 0, total_science: int = 1, total_imgs: int = 1):
        """Handle pipeline completion.

        Args:
            success: Whether processing was successful
            message: Status message
            img_idx: Index of the current science image
            total_science: Total number of science images
            total_imgs: Total number of all images processed (bias + flat + calib + science)
        """
        self.log_text.append("\n" + "=" * 60)
        self.log_text.append(message)
        self.log_text.append("=" * 60)

        if success:
            if img_idx == total_science - 1:
                # Last science image processed
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
            # Update overscan checkbox state after settings are saved
            self._update_overscan_checkbox_from_config()
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
        # 由于日志实时显示在GUI中难以实现，我们只保留终端和文件日志
        # 不再添加GUI日志处理器，避免性能问题
        pass

    def closeEvent(self, event):
        """Handle window close."""
        if self.worker_thread and self.worker_thread.isRunning():
            self.worker_thread.terminate()
            self.worker_thread.wait()
        super().closeEvent(event)

    def _refresh_step1_checkbox_text(self):
        """Update Step1 display text so options are shown in parentheses."""
        stage_0_checkbox = self.stage_checkboxes.get('stage_0')
        if not stage_0_checkbox:
            return

        if not hasattr(self, 'step1_substep_checkboxes') or not self.step1_substep_checkboxes:
            stage_0_checkbox.setText("Step 1: Basic Pre-processing")
            return

        selected = []
        for key, label in [('overscan', 'Overscan'), ('bias', 'Bias'), ('cosmic', 'Cosmic')]:
            cb = self.step1_substep_checkboxes.get(key)
            if cb and cb.isChecked():
                selected.append(label)

        options_text = '/'.join(selected) if selected else 'None'
        stage_0_checkbox.setText(f"Step 1: Basic Pre-processing ({options_text})")
    
    def _update_overscan_checkbox_from_config(self):
        """Update overscan checkbox based on configuration."""
        overscan_start_col = self.config.get_int('data', 'overscan_start_column', -1)
        stage_0_checkbox = self.stage_checkboxes.get('stage_0')
        overscan_cb = self.step1_substep_checkboxes.get('overscan')
        
        if stage_0_checkbox and overscan_cb:
            # DO NOT set stage_0_checkbox.setChecked() here. This is the user's explicit choice.
            # Only update the 'overscan' sub-checkbox based on config.

            # Temporarily block signals to prevent triggering the change handler
            overscan_cb.blockSignals(True)
            
            # If overscan_start_column is -1 (disabled), uncheck the checkbox
            overscan_cb.setChecked(overscan_start_col != -1)
            
            # Restore signals for the sub-checkbox
            overscan_cb.blockSignals(False)

            # Keep Step1 substep toggles disabled when Step1 itself is unchecked.
            enabled = stage_0_checkbox.isChecked()
            for cb in self.step1_substep_checkboxes.values():
                cb.setEnabled(enabled)

            self._refresh_step1_checkbox_text()
    
    def _on_overscan_checkbox_changed(self, state):
        """Handle overscan checkbox state change."""
        stage_0_checkbox = self.stage_checkboxes.get('stage_0')
        if not stage_0_checkbox:
            return

        enabled = (state == Qt.Checked)
        for cb in self.step1_substep_checkboxes.values():
            cb.setEnabled(enabled)
        self._refresh_step1_checkbox_text()
            
        overscan_start_col = self.config.get_int('data', 'overscan_start_column', -1)
        overscan_cb = self.step1_substep_checkboxes.get('overscan')
        
        if state == Qt.Checked and overscan_cb and overscan_cb.isChecked():
            # Check if overscan_start_column is set (not -1)
            if overscan_start_col == -1:
                # Ask user to set the overscan configuration
                QMessageBox.information(
                    self, "Overscan Configuration Required",
                    "Overscan start column is not configured. Please set the overscan "
                    "configuration in the Settings dialog."
                )
                
                # Uncheck only the overscan substep, keep Step 1 checked
                overscan_cb.blockSignals(True)
                overscan_cb.setChecked(False)
                overscan_cb.blockSignals(False)
                self._refresh_step1_checkbox_text()
                
                # Open settings dialog to the Data Reduction tab (index 1) which contains overscan settings
                from src.gui.settings_dialog import SettingsDialog
                dialog = SettingsDialog(self.config, self, initial_tab=1)
                result = dialog.exec_()
                
                # After dialog closes, update checkbox state based on configuration
                self._update_overscan_checkbox_from_config()
                
                if result == QDialog.Accepted:
                    overscan_start_col = self.config.get_int('data', 'overscan_start_column', -1)
                    if overscan_start_col != -1:
                        overscan_cb.blockSignals(True)
                        overscan_cb.setChecked(True)
                        overscan_cb.blockSignals(False)
                        self._refresh_step1_checkbox_text()
                        self.statusBar.showMessage(f"Overscan configuration updated: start column = {overscan_start_col}")
                    else:
                        self.statusBar.showMessage("Overscan configuration saved (disabled)")
                else:
                    self.statusBar.showMessage("Settings changes cancelled")
            else:
                # Configuration is already set, just update status message
                self.statusBar.showMessage(f"Overscan correction enabled (start column = {overscan_start_col})")
        else:
            # If unchecked, just update status message - don't modify configuration
            self.statusBar.showMessage("Step 1 toggled (configuration preserved)")
