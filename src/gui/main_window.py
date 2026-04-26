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
            ("stage_5", "Step 6: Wavelength Calibration", True),
            ("stage_6", "Step 7: Blaze Correction", True),
            ("stage_7", "Step 8: Order Stitching", True),
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

        # Handle calib files (can be list or string)
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
                lambda success, msg, img_idx=idx, total_science=len(science_files), total_imgs=total_imgs:
                    self._on_execution_complete(success, msg, img_idx, total_science, total_imgs)
            )
            self.worker_thread.start()
            self.worker_thread.wait()  # Wait for this image to complete before next

        self.run_all_btn.setEnabled(True)
        self.stop_btn.setEnabled(False)

    def _run_selected_pipeline(self):
        """Execute only the selected processing steps."""
        # Reset run-scoped cache to avoid reusing stale paths from previous runs.
        self._flat_files_after_bias = None
        self._flat_files_after_step1 = None

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
            "stage_5": "Step 6: Wavelength Calibration",
            "stage_6": "Step 7: Blaze Correction",
            "stage_7": "Step 8: Order Stitching",
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

        # Get calib file
        calib_file = None
        if isinstance(self.calib_file, list) and self.calib_file:
            calib_file = sorted(self.calib_file, key=lambda x: Path(x).stat().st_mtime)[-1]

        # Create a modified pipeline for selected steps
        # Process all images with overscan correction if stage_0 is selected
        if "stage_0" in selected_stages:
            self.log_text.append("\n" + "=" * 60)
            self.log_text.append("STEP 1.1: OVERSCAN CORRECTION")
            self.log_text.append("Processing all image types in order: Bias → Flat → Calibration → Science")
            self.log_text.append("=" * 60)
            
            try:
                from src.core.overscan_correction import process_overscan_stage
                
                # Collect all files to process in the correct order
                all_files_to_process = []
                file_groups = []
                
                # 1. Bias files
                if self.bias_files:
                    for bias_file in self.bias_files:
                        file_name = Path(bias_file).name
                        all_files_to_process.append(bias_file)
                        file_groups.append(("bias", bias_file))
                
                # 2. Flat files
                if self.flat_files:
                    for flat_file in self.flat_files:
                        file_name = Path(flat_file).name
                        all_files_to_process.append(flat_file)
                        file_groups.append(("flat", flat_file))
                
                # 3. Calibration files (灯谱)
                calib_list = self.calib_file if isinstance(self.calib_file, list) else [self.calib_file]
                if calib_list and any(calib_list):
                    for calib_file in calib_list:
                        if calib_file:  # Skip if None
                            file_name = Path(calib_file).name
                            all_files_to_process.append(calib_file)
                            file_groups.append(("calib", calib_file))
                
                # 4. Science files
                science_files_to_process = []
                for science_file in science_files:
                    file_name = Path(science_file).name
                    all_files_to_process.append(science_file)
                    file_groups.append(("science", science_file))
                    science_files_to_process.append(science_file)
                
                # Check which files actually need processing
                output_dir = self.config.get_output_path()
                overscan_dir = Path(output_dir) / 'step1_basic' / 'overscan_corrected'
                
                # Get file status
                existing_files, missing_files = self._get_overscan_file_status()
                
                # Determine what action to take based on scenarios
                force_all = "stage_0" in force_reprocess_stages
                should_process = False
                skip_overscan = False  # Flag to skip overscan processing
                
                if force_all:
                    # User explicitly chose to reprocess all from previous dialog
                    should_process = True
                    self.log_text.append("\n  Reprocessing all files (user requested)")
                else:
                    # No dialog - determine action based on file status
                    total_files = len(existing_files) + len(missing_files)
                    
                    if len(missing_files) == 0:
                        # All files exist - reprocess all
                        should_process = True
                        force_all = True
                        self.log_text.append(f"\n  Reprocessing all {len(existing_files)} files")
                    elif len(existing_files) == 0:
                        # All missing - process all
                        should_process = True
                        self.log_text.append(f"\n  All {total_files} files need processing")
                    else:
                        # Some exist, some missing - process only missing
                        should_process = True
                        self.log_text.append(f"\n  Processing only missing files: {len(missing_files)}")
                
                # Build file lists based on user choice
                files_to_process = []
                
                if should_process:
                    for file_type, file_path in file_groups:
                        file_name = Path(file_path).name
                        
                        # Check if file_path exists, if not, use the original raw path
                        if not Path(file_path).exists():
                            # Try to find original raw file
                            if file_name in self.original_file_paths:
                                original_path = self.original_file_paths[file_name]
                                if Path(original_path).exists():
                                    file_path = original_path
                                else:
                                    continue  # Skip if original not found
                            else:
                                continue  # Skip if no original tracked
                        
                        processed_file = overscan_dir / file_name
                        
                        # Determine if this file should be processed
                        if force_all:
                            # Reprocessing all files - use original raw file
                            if file_name in self.original_file_paths:
                                raw_file = self.original_file_paths[file_name]
                                if Path(raw_file).exists():
                                    files_to_process.append(raw_file)
                                else:
                                    files_to_process.append(file_path)
                            else:
                                files_to_process.append(file_path)
                        else:
                            # Only process missing files
                            if not processed_file.exists():
                                files_to_process.append(file_path)
                            else:
                                # Update path to processed file
                                if file_type == "bias":
                                    idx = self.bias_files.index(file_path)
                                    self.bias_files[idx] = str(processed_file)
                                elif file_type == "flat":
                                    idx = self.flat_files.index(file_path)
                                    self.flat_files[idx] = str(processed_file)
                                elif file_type == "calib":
                                    if isinstance(self.calib_file, list):
                                        idx = self.calib_file.index(file_path)
                                        self.calib_file[idx] = str(processed_file)
                                    else:
                                        self.calib_file = str(processed_file)
                                elif file_type == "science":
                                    idx = science_files.index(file_path)
                                    science_files[idx] = str(processed_file)
                
                # Mark force_reprocess_stages for consistency
                if force_all:
                    force_reprocess_stages.add("stage_0")
                
                # Process files if any need processing
                if files_to_process:
                    self.log_text.append("  Processing images with overscan correction...")
                    
                    # Process in batches to show progress
                    processed_files = process_overscan_stage(self.config, files_to_process, None)
                    
                    # Update file paths for processed science files
                    for processed_file in processed_files:
                        processed_name = Path(processed_file).name
                        # Find which science file this corresponds to
                        for i, science_file in enumerate(science_files_to_process):
                            if Path(science_file).name == processed_name:
                                # Make sure index is valid for science_files list
                                if i < len(science_files):
                                    science_files[i] = processed_file
                                else:
                                    logger.warning(f"Index {i} out of bounds for science_files list")
                                break
                    
                    self.log_text.append(f"  ✓ Overscan correction complete: {len(processed_files)} images processed")
                    
                    # Update processing state
                    self.processing_state['stage_0_completed'] = True
                    # Collect all overscan-corrected file paths
                    overscan_corrected_files = []
                    for file_type, file_path in file_groups:
                        file_name = Path(file_path).name
                        processed_file = overscan_dir / file_name
                        if processed_file.exists():
                            overscan_corrected_files.append(str(processed_file))
                    self.processing_state['overscan_corrected_files'] = overscan_corrected_files
                else:
                    self.log_text.append("\n  All files already overscan-corrected")
                    # User chose not to reprocess - don't mark stage_0 as completed
                    # but still collect overscan file paths for subsequent stages to use
                    overscan_corrected_files = []
                    for file_type, file_path in file_groups:
                        file_name = Path(file_path).name
                        processed_file = overscan_dir / file_name
                        if processed_file.exists():
                            overscan_corrected_files.append(str(processed_file))
                    self.processing_state['overscan_corrected_files'] = overscan_corrected_files
                    
            except Exception as e:
                self.log_text.append(f"  ✗ Overscan correction failed: {e}")
                # Continue anyway, the error will be caught when processing further stages

        # Continue with selected-step execution.
        # Step 1 includes both overscan and bias preprocessing.
        other_stages_selected = True
        
        if other_stages_selected:
            # Pre-combine bias frames once before processing science images
            master_bias_path = None
            output_dir = self.config.get_output_path()
            overscan_dir = Path(output_dir) / 'step1_basic' / 'overscan_corrected'
            
            # Check if overscan-corrected files exist for science, bias and flat
            overscan_science_files = []
            overscan_bias_files = []
            overscan_flat_files = []
            
            if overscan_dir.exists():
                # Check science files
                for sci_file in science_files:
                    overscan_sci_file = overscan_dir / Path(sci_file).name
                    if overscan_sci_file.exists():
                        overscan_science_files.append(str(overscan_sci_file))
                
                # Check bias files - use overscan version if exists, otherwise fall back to original
            all_bias_overscan = True  # Track if all bias files are overscan-corrected
            for bias_file in self.bias_files:
                bias_name = Path(bias_file).name
                overscan_bias_file = overscan_dir / bias_name
                if overscan_bias_file.exists():
                    overscan_bias_files.append(str(overscan_bias_file))
                elif bias_name in self.original_file_paths:
                    original = self.original_file_paths[bias_name]
                    if Path(original).exists():
                        overscan_bias_files.append(original)
                        all_bias_overscan = False
                        logger.info(f"Using original bias file (overscan version missing): {bias_name}")

            # Check flat files - use overscan version if exists, otherwise fall back to original
            all_flat_overscan = True  # Track if all flat files are overscan-corrected
            for flat_file in self.flat_files:
                flat_name = Path(flat_file).name
                overscan_flat_file = overscan_dir / flat_name
                if overscan_flat_file.exists():
                    overscan_flat_files.append(str(overscan_flat_file))
                elif flat_name in self.original_file_paths:
                    original = self.original_file_paths[flat_name]
                    if Path(original).exists():
                        overscan_flat_files.append(original)
                        all_flat_overscan = False
                        logger.info(f"Using original flat file (overscan version missing): {flat_name}")
            
            # Determine which version to use for science/flat inputs.
            # For selected-step execution, Step 2+ must depend ONLY on the substeps
            # the user checked in Step 1, rather than silently reusing old outputs.
            use_overscan_science = step1_do_overscan
            use_overscan_bias = (len(overscan_bias_files) == len(self.bias_files)) if self.bias_files else False
            
            # Must be consistent: both science and bias use overscan OR both use raw.
            # If overscan was not explicitly selected, later steps start from raw inputs.
            if use_overscan_science and use_overscan_bias and all_bias_overscan:
                self._use_overscan_science = True
                # Update science_files to use overscan-corrected versions
                for i, sci_file in enumerate(science_files):
                    sci_name = Path(sci_file).name
                    overscan_sci_file = overscan_dir / sci_name
                    if overscan_sci_file.exists():
                        science_files[i] = str(overscan_sci_file)
                    elif sci_name in self.original_file_paths:
                        # Fall back to original raw file if overscan version missing
                        original = self.original_file_paths[sci_name]
                        if Path(original).exists():
                            science_files[i] = original
            else:
                self._use_overscan_science = False
                # Log the reason for falling back to raw files
                if step1_do_overscan and not all_bias_overscan:
                    logger.info("Some bias files missing overscan correction - using raw files for both science and bias")
                elif step1_do_overscan and not use_overscan_bias:
                    logger.info("Bias overscan files incomplete - using raw files for both science and bias")
                elif step1_do_overscan and not use_overscan_science:
                    logger.info("Science overscan files incomplete - using raw files for both science and bias")
                # Use original files for science if bias uses originals
                for i, sci_file in enumerate(science_files):
                    sci_name = Path(sci_file).name
                    if sci_name in self.original_file_paths:
                        original = self.original_file_paths[sci_name]
                        if Path(original).exists():
                            science_files[i] = original
            
            # Step 1 chaining for Step 2: start from raw flats, then advance only
            # through the substeps explicitly selected by the user.
            current_flat_files = [f for f in self.flat_files if Path(f).exists()]
            if step1_do_overscan:
                overscan_selected_flats = [f for f in overscan_flat_files if Path(f).exists()]
                if len(overscan_selected_flats) == len(current_flat_files):
                    current_flat_files = overscan_selected_flats

            if step1_do_bias and self.bias_files:
                self.log_text.append("\n============================================================")
                self.log_text.append("STEP 1.2: BIAS SUBTRACTION")
                from src.core.bias_correction import BiasCorrector
                
                # Check which science files already have bias-corrected versions
                output_dir = self.config.get_output_path()
                bias_output_dir = Path(output_dir) / 'step1_basic' / 'bias_subtracted'
                
                existing_bias_sci = []  # Files that already have bias correction
                missing_bias_sci = []   # Files that need bias correction
                science_files_bias_paths = {}  # Map: source_file -> bias_corrected_path
                
                # Build list of science file names
                sci_file_names = [Path(f).name for f in science_files]
                
                for sci_name in sci_file_names:
                    bias_sci_path = bias_output_dir / sci_name
                    if bias_sci_path.exists():
                        existing_bias_sci.append(sci_name)
                        science_files_bias_paths[sci_name] = str(bias_sci_path)
                    else:
                        missing_bias_sci.append(sci_name)
                
                # Determine which science files need bias correction
                if len(missing_bias_sci) == 0:
                    # All exist - reprocess all
                    reprocess_all = True
                    self.log_text.append(f"  All {len(science_files)} files already bias-corrected, reprocessing all")
                else:
                    # Some missing - only process missing ones
                    reprocess_all = False
                    self.log_text.append(f"  Already bias-corrected: {len(existing_bias_sci)}, missing: {len(missing_bias_sci)}")
                
                # Update science_files to use bias-corrected paths for existing files
                for i, sci_file in enumerate(science_files):
                    sci_name = Path(sci_file).name
                    if sci_name in science_files_bias_paths:
                        science_files[i] = science_files_bias_paths[sci_name]
                
                # Get list of files to actually process
                science_files_to_bias = []
                if reprocess_all:
                    # Reprocess all - use source files consistent with bias files
                    for sci_file in science_files:
                        sci_name = Path(sci_file).name
                        if self._use_overscan_science:
                            # Use overscan-corrected version
                            overscan_file = overscan_dir / sci_name
                            if overscan_file.exists():
                                science_files_to_bias.append(str(overscan_file))
                        else:
                            # Use raw original version
                            if sci_name in self.original_file_paths:
                                source_file = self.original_file_paths[sci_name]
                                if Path(source_file).exists():
                                    science_files_to_bias.append(source_file)
                else:
                    # Only process missing ones
                    for sci_name in missing_bias_sci:
                        # Find the source file path (overscan-corrected or raw)
                        for sci_file in science_files:
                            if Path(sci_file).name == sci_name:
                                science_files_to_bias.append(sci_file)
                                break
                
                try:
                    # Determine which bias files to use
                    if self._use_overscan_science:
                        bias_files_to_use = overscan_bias_files
                        self.log_text.append(f"  Using bias files (with overscan correction)")
                    else:
                        # Use original bias files
                        bias_files_to_use = []
                        for bias_file in self.bias_files:
                            bias_name = Path(bias_file).name
                            if bias_name in self.original_file_paths:
                                original = self.original_file_paths[bias_name]
                                if Path(original).exists():
                                    bias_files_to_use.append(original)
                        self.log_text.append(f"  Using raw bias files")
                    
                    self.log_text.append(f"  Combining {len(bias_files_to_use)} bias frames...")
                    corrector = BiasCorrector(self.config)
                    master_bias, _ = corrector.combine_bias_frames(bias_files_to_use)
                    
                    # Save master bias
                    bias_output_dir.mkdir(parents=True, exist_ok=True)
                    master_bias_path = bias_output_dir / 'master_bias.fits'
                    corrector.save_master_bias(str(master_bias_path))
                    self.log_text.append(f"  Master bias saved: {master_bias_path.name}")
                    self._master_bias = master_bias
                except Exception as e:
                    self.log_text.append(f"  ✗ Failed to combine bias frames: {e}")
                    raise

                # Step1 basic preprocessing applies to both science and flat frames.
                need_flat_preprocess = True
                flat_sources_for_bias = []
                if need_flat_preprocess:
                    flat_sources_for_bias = current_flat_files
                    flat_sources_for_bias = [f for f in flat_sources_for_bias if Path(f).exists()]
                
                # Apply bias correction only to files that need it
                if science_files_to_bias:
                    if reprocess_all:
                        self.log_text.append(f"  Applying bias correction to all files: {len(science_files_to_bias)} image(s)")
                    else:
                        self.log_text.append(f"  Applying bias correction to missing files: {len(science_files_to_bias)} image(s)")
                    from src.utils.fits_io import read_fits_image, write_fits_image
                    from src.plotting.spectra_plotter import plot_2d_image_to_file
                    fig_format = self.config.get('reduce', 'fig_format', 'png')
                    
                    for sci_file in science_files_to_bias:
                        try:
                            sci_img, sci_header = read_fits_image(sci_file)
                            corrector = BiasCorrector(self.config)
                            corrector.master_bias = master_bias
                            bias_corrected = corrector.apply_bias_correction(sci_img, master_bias)
                            
                            output_science_path = bias_output_dir / Path(sci_file).name
                            sci_header['BIASCOR'] = (True, 'Bias correction applied')
                            write_fits_image(str(output_science_path), bias_corrected, header=sci_header, dtype='float32')
                            
                            # Save diagnostic plot
                            save_plots = self.config.get_bool('reduce', 'save_plots', True)
                            if save_plots:
                                plot_file = bias_output_dir / f'{Path(sci_file).stem}_bias_corrected.{fig_format}'
                                plot_2d_image_to_file(bias_corrected, str(plot_file), "Bias Corrected Science Image")
                            
                            # Update science_files list to use bias-corrected version (by filename)
                            sci_name = Path(sci_file).name
                            for i, sf in enumerate(science_files):
                                if Path(sf).name == sci_name:
                                    science_files[i] = str(output_science_path)
                                    break
                        except Exception as e:
                            self.log_text.append(f"    ✗ {Path(sci_file).name}: {e}")
                    
                    self.log_text.append(f"  ✓ Bias correction complete")

                # Apply bias correction to flat frames too (required for correct step 2 input)
                bias_corrected_flat_files = []

                if need_flat_preprocess and flat_sources_for_bias:
                    self.log_text.append(f"  Applying bias correction to flat frames: {len(flat_sources_for_bias)} image(s)")
                    from src.utils.fits_io import read_fits_image, write_fits_image
                    from src.plotting.spectra_plotter import plot_2d_image_to_file
                    save_plots = self.config.get_bool('reduce', 'save_plots', True)
                    fig_format = self.config.get('reduce', 'fig_format', 'png')

                    for flat_file in flat_sources_for_bias:
                        try:
                            flat_img, flat_header = read_fits_image(flat_file)
                            corrector = BiasCorrector(self.config)
                            corrector.master_bias = master_bias
                            flat_corrected = corrector.apply_bias_correction(flat_img, master_bias)

                            output_flat_path = bias_output_dir / Path(flat_file).name
                            flat_header['BIASCOR'] = (True, 'Bias correction applied')
                            write_fits_image(str(output_flat_path), flat_corrected, header=flat_header, dtype='float32')

                            if save_plots:
                                plot_file = bias_output_dir / f'{Path(flat_file).stem}_bias_corrected.{fig_format}'
                                plot_2d_image_to_file(flat_corrected, str(plot_file), "Bias Corrected Flat Image")

                            bias_corrected_flat_files.append(str(output_flat_path))
                        except Exception as e:
                            self.log_text.append(f"    ✗ flat {Path(flat_file).name}: {e}")

                    self.log_text.append(f"  ✓ Flat bias correction complete: {len(bias_corrected_flat_files)} saved")
                elif need_flat_preprocess:
                    self.log_text.append("  ! No flat frames available for bias correction")
                else:
                    self.log_text.append("  ! No flat frames available for bias correction")

                self._flat_files_after_bias = bias_corrected_flat_files if bias_corrected_flat_files else flat_sources_for_bias
                current_flat_files = list(self._flat_files_after_bias)

            # Apply cosmic-ray removal in Step1 basic preprocessing (science only).
            if step1_do_cosmic:
                if self.config.get_bool('reduce', 'cosmic_enabled', True):
                    from src.core.cosmic_removal import process_cosmic_stage

                    self.log_text.append("\n============================================================")
                    self.log_text.append("STEP 1.3: COSMIC RAY REMOVAL")
                    self.log_text.append("============================================================")
                    self.log_text.append(
                        "  Applying cosmic ray removal to science and flat frames "
                        "(using the latest selected Step 1 output as input)..."
                    )

                    # If bias wasn't selected, the current science_files list already reflects
                    # either overscan outputs (if selected) or raw inputs (if not selected).

                    # Cosmic-correct science files
                    if science_files:
                        science_cosmic_outputs = process_cosmic_stage(
                            self.config,
                            science_files,
                            None,
                            output_subdir='step1_basic/cosmic_corrected'
                        )

                        science_cosmic_map = {Path(p).name: p for p in science_cosmic_outputs}
                        for i, sci_file in enumerate(science_files):
                            sci_name = Path(sci_file).name
                            if sci_name in science_cosmic_map:
                                science_files[i] = science_cosmic_map[sci_name]

                    # Cosmic-correct flat files too so Step 2 starts from the final Step 1 output.
                    if current_flat_files:
                        flat_cosmic_outputs = process_cosmic_stage(
                            self.config,
                            current_flat_files,
                            None,
                            output_subdir='step1_basic/cosmic_corrected'
                        )
                        current_flat_files = [p for p in flat_cosmic_outputs if Path(p).exists()]

                    self._flat_files_after_step1 = list(current_flat_files)

                    self.log_text.append("  ✓ Cosmic ray removal complete (saved in step1_basic/cosmic_corrected)")
                else:
                    self.log_text.append("\n============================================================")
                    self.log_text.append("STEP 1.3: COSMIC RAY REMOVAL")
                    self.log_text.append("============================================================")
                    self.log_text.append("  Cosmic ray removal disabled by config")
            elif step1_selected:
                self.log_text.append("\n============================================================")
                self.log_text.append("STEP 1.3: COSMIC RAY REMOVAL")
                self.log_text.append("============================================================")
                self.log_text.append("  Skipped (option unchecked)")
                
                # Update processing state
                self.processing_state['stage_1_completed'] = True
                self.processing_state['master_bias_path'] = str(master_bias_path)

            # If Step 1 ended without cosmic, the latest Step 1 flat input is either
            # raw, overscan-corrected, or bias-corrected depending on selections.
            if self._flat_files_after_step1 is None:
                self._flat_files_after_step1 = list(current_flat_files)
            
            # Only run individual processing loop if there are stages beyond step 1
            stages_beyond_1 = [s for s in selected_stages if s not in ("stage_0", "stage_1")]
            
            if "stage_1" in selected_stages and self.flat_files:
                self.log_text.append("\n============================================================")
                self.log_text.append("STEP 2: Order Tracing")
                from src.core.flat_fielding import process_flat_stage
                
                # Step 2 must start from the final selected Step 1 output:
                # raw -> overscan -> bias -> cosmic, whichever is the last one executed.
                flat_files_to_use = getattr(self, '_flat_files_after_step1', None)
                if not flat_files_to_use:
                    flat_files_to_use = [f for f in self.flat_files if Path(f).exists()]

                self.log_text.append(f"  Combining {len(flat_files_to_use)} flat frames...")
                flat_field, apertures = process_flat_stage(self.config, flat_files_to_use)
                self.log_text.append(f"  ✓ Flat field created: {len(apertures)} orders detected")
                self._flat_field = flat_field
                self._apertures = apertures
                
                # Update processing state
                self.processing_state['stage_2_completed'] = True
            
            if stages_beyond_1:
                self.log_text.append("\n============================================================")
                self.log_text.append(
                    f"Processing {len(science_files)} science image(s) from the final Step 1 output..."
                )
                for i, sci_file in enumerate(science_files, start=1):
                    self.log_text.append(f"[{i}/{len(science_files)}] {Path(sci_file).name}")
                self.log_text.append("============================================================")

                had_error = False
                for idx, science_file in enumerate(science_files):
                    try:
                        # Create a custom processing based on selected stages
                        self._execute_selected_stages(science_file, selected_stages, idx, len(science_files))
                    except Exception as e:
                        had_error = True
                        self._on_execution_complete(False, f"Error: {str(e)}", idx, len(science_files), total_imgs)
                        break

                if not had_error:
                    self.log_text.append(f"\n✓ {len(science_files)} science images processed")
                else:
                    self.log_text.append("\n✗ Processing stopped due to error")
        else:
            # Only stage_0 is selected
            # Summary is already printed in the processing section above
            self.statusBar.showMessage("Overscan correction completed")

        self.run_all_btn.setEnabled(True)
        self.run_selected_btn.setEnabled(True)
        self.stop_btn.setEnabled(False)

    def _execute_selected_stages(self, science_file: str, selected_stages: List[str], idx: int = 0, total: int = 1):
        """Execute only the selected processing stages."""
        from src.core.overscan_correction import process_overscan_stage
        from src.core.flat_fielding import process_flat_stage
        from src.core.wave_calibration import process_wavelength_stage
        from src.core.scattered_light import process_background_stage
        from src.core.extraction import process_extraction_stage
        from src.utils.fits_io import read_fits_image, write_fits_image

        self.progress_bar.setValue(10)
        self.stage_label.setText(f"Processing [{idx+1}/{total}]...")

        total_stages = len(selected_stages)
        current_stage = 0

        # For overscan correction, check if already processed
        if "stage_0" in selected_stages:
            current_stage += 1
            progress = (current_stage / total_stages) * 100
            self.progress_bar.setValue(int(progress))
            self.stage_label.setText(f"[{idx+1}/{total}] Overscan Correction")

        # For other stages, check if overscan-corrected science file should be used
        if hasattr(self, '_use_overscan_science') and self._use_overscan_science:
            output_dir = self.config.get_output_path()
            overscan_dir = Path(output_dir) / 'step1_basic' / 'overscan_corrected'
            overscan_science_file = overscan_dir / Path(science_file).name
            if overscan_science_file.exists():
                science_file = str(overscan_science_file)
        
        if "stage_1" in selected_stages:
            current_stage += 1
            progress = (current_stage / total_stages) * 100
            self.progress_bar.setValue(int(progress))
            self.stage_label.setText(f"[{idx+1}/{total}] Bias Correction")
            # Bias correction is already handled before the loop - science_file is already bias-corrected

        if "stage_2" in selected_stages:
            current_stage += 1
            progress = (current_stage / total_stages) * 100
            self.progress_bar.setValue(int(progress))
            self.stage_label.setText(f"[{idx+1}/{total}] Scattered Light Subtraction")

            # Step 3: build/subtract scattered-light background for each science image.
            sci_img, sci_header = read_fits_image(science_file)
            sci_name = Path(science_file).stem
            apertures = getattr(self, '_apertures', None)

            background = process_background_stage(
                self.config,
                sci_img,
                None,
                output_subdir='step3_scatterlight',
                output_tag=f'{sci_name}_',
                apertures=apertures,
                mask_margin_scale=1.2,
            )
            corrected = sci_img.astype('float32') - background.astype('float32')

            out_dir = Path(self.config.get_output_path()) / 'step3_scatterlight'
            out_dir.mkdir(parents=True, exist_ok=True)
            corrected_path = out_dir / f'{sci_name}_science_background_subtracted.fits'

            if sci_header is not None:
                sci_header['BKGSCAT'] = (True, 'Scattered light background subtracted')
            write_fits_image(str(corrected_path), corrected, header=sci_header, dtype='float32')

            self.log_text.append(f"    ✓ Step 3 output: {corrected_path.name}")
            science_file = str(corrected_path)

        if "stage_3" in selected_stages:
            current_stage += 1
            progress = (current_stage / total_stages) * 100
            self.progress_bar.setValue(int(progress))
            self.stage_label.setText(f"[{idx+1}/{total}] 2D Flat-Field Correction")

            # Step 4: divide science image by the pixel flat from Step 2.
            flat_field = getattr(self, '_flat_field', None)
            pixel_flat = flat_field.pixel_flat if flat_field is not None else None
            if pixel_flat is None:
                self.log_text.append("  ! Step 4: no pixel_flat available – Step 2 must run first")
            else:
                sci_img, sci_header = read_fits_image(science_file)
                safe_flat = np.where((pixel_flat > 0.1) & np.isfinite(pixel_flat), pixel_flat, 1.0)
                corrected = sci_img.astype('float32') / safe_flat.astype('float32')
                out_dir = Path(self.config.get_output_path()) / 'step4_flat_corrected'
                out_dir.mkdir(parents=True, exist_ok=True)
                sci_name = Path(science_file).stem
                corrected_path = out_dir / f'{sci_name}_flat_corrected.fits'
                if sci_header is not None:
                    sci_header['FLATCOR'] = (True, '2D pixel flat correction applied')
                write_fits_image(str(corrected_path), corrected, header=sci_header, dtype='float32')
                self.log_text.append(f"  ✓ Step 4 output: {corrected_path.name}")
                science_file = str(corrected_path)

        if "stage_4" in selected_stages:
            current_stage += 1
            progress = (current_stage / total_stages) * 100
            self.progress_bar.setValue(int(progress))
            self.stage_label.setText(f"[{idx+1}/{total}] Spectrum Extraction")

            apertures = getattr(self, '_apertures', None)
            if apertures is None:
                self.log_text.append("  ! Step 5: no apertures available – Step 2 must run first")
            else:
                sci_img, _ = read_fits_image(science_file)
                sci_name = Path(science_file).stem
                spectra = process_extraction_stage(
                    self.config, sci_img, apertures,
                    output_filename=f'{sci_name}_extracted_spectra.fits',
                    plot_prefix=sci_name,
                )
                self.log_text.append(f"  ✓ Step 5 extraction complete: {len(spectra.spectra)} orders")

        if "stage_5" in selected_stages:
            current_stage += 1
            progress = (current_stage / total_stages) * 100
            self.progress_bar.setValue(int(progress))
            self.stage_label.setText(f"[{idx+1}/{total}] Wavelength Calibration")
            calib_for_wave = self.calib_file
            if isinstance(self.calib_file, list) and self.calib_file:
                calib_for_wave = sorted(self.calib_file, key=lambda x: Path(x).stat().st_mtime)[-1]
            wave_calib = process_wavelength_stage(self.config, calib_for_wave)

        # Final progress update
        self.progress_bar.setValue(100)
        self.stage_label.setText("Completed")

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
        
        if stage_0_checkbox:
            # Temporarily block signals to prevent triggering the change handler
            stage_0_checkbox.blockSignals(True)
            
            # If overscan_start_column is -1 (disabled), uncheck the checkbox
            if overscan_start_col == -1:
                stage_0_checkbox.setChecked(False)
            else:
                # If overscan_start_column is set to a valid value, check the checkbox
                stage_0_checkbox.setChecked(True)
            
            # Restore signals
            stage_0_checkbox.blockSignals(False)

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
        
        if state == Qt.Checked:
            # Check if overscan_start_column is set (not -1)
            if overscan_start_col == -1:
                # Ask user to set the overscan configuration
                QMessageBox.information(
                    self, "Overscan Configuration Required",
                    "Overscan start column is not configured. Please set the overscan "
                    "configuration in the Settings dialog."
                )
                
                # Temporarily block signals to prevent recursion
                stage_0_checkbox.blockSignals(True)
                stage_0_checkbox.setChecked(False)
                stage_0_checkbox.blockSignals(False)
                
                # Open settings dialog to the Data Reduction tab (index 1) which contains overscan settings
                from src.gui.settings_dialog import SettingsDialog
                dialog = SettingsDialog(self.config, self, initial_tab=1)
                result = dialog.exec_()
                
                # After dialog closes, update checkbox state based on configuration
                self._update_overscan_checkbox_from_config()
                
                if result == QDialog.Accepted:
                    overscan_start_col = self.config.get_int('data', 'overscan_start_column', -1)
                    if overscan_start_col != -1:
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
            self.statusBar.showMessage("Overscan correction disabled (configuration preserved)")
