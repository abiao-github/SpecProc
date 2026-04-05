"""
Settings dialog for SpecProc configuration.

Provides a GUI for editing configuration parameters organized by tabs.
"""

from PyQt5.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit,
    QSpinBox, QDoubleSpinBox, QComboBox, QCheckBox,
    QPushButton, QTabWidget, QWidget, QGroupBox, QFormLayout, QScrollArea
)
from PyQt5.QtCore import Qt

from src.config.config_manager import ConfigManager


class SettingsDialog(QDialog):
    """Settings dialog for configuration editing."""

    def __init__(self, config: ConfigManager, parent=None):
        """
        Initialize settings dialog.

        Args:
            config: Configuration manager instance
            parent: Parent widget
        """
        super().__init__(parent)
        self.config = config
        self.setWindowTitle("SpecProc Settings")
        self.setModal(True)
        self.resize(800, 600)

        self.init_ui()
        self.load_current_values()

    def init_ui(self):
        """Initialize user interface."""
        layout = QVBoxLayout()

        # Tab widget with different configuration categories
        tab_widget = QTabWidget()

        # Create tabs merging related parameters
        tab_widget.addTab(self._create_data_tab(), "Data & Calibration")
        tab_widget.addTab(self._create_correction_tab(), "Bias & Flat")
        tab_widget.addTab(self._create_wlcalib_tab(), "Wavelength Calibration")
        tab_widget.addTab(self._create_extraction_tab(), "Trace & Extraction")
        tab_widget.addTab(self._create_processing_tab(), "Processing & Output")

        layout.addWidget(tab_widget)

        # Buttons
        button_layout = QHBoxLayout()
        save_btn = QPushButton("Save")
        save_btn.clicked.connect(self._save_settings)
        cancel_btn = QPushButton("Cancel")
        cancel_btn.clicked.connect(self.reject)
        button_layout.addStretch()
        button_layout.addWidget(save_btn)
        button_layout.addWidget(cancel_btn)
        layout.addLayout(button_layout)

        self.setLayout(layout)

    def _create_data_tab(self) -> QWidget:
        """Create data input, telescope and overscan settings tab."""
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        widget = QWidget()
        layout = QFormLayout()

        # Data input group
        data_group = QGroupBox("Data Input")
        data_layout = QFormLayout()

        self.rawdata_path_edit = QLineEdit()
        data_layout.addRow("Raw Data Path:", self.rawdata_path_edit)

        self.statime_key_edit = QLineEdit()
        data_layout.addRow("Start Time Key (DATE-OBS):", self.statime_key_edit)

        self.exptime_key_edit = QLineEdit()
        data_layout.addRow("Exposure Time Key (EXPTIME):", self.exptime_key_edit)

        self.direction_combo = QComboBox()
        self.direction_combo.addItems(['xr-', 'xl-', 'yr-', 'yl-'])
        data_layout.addRow("Dispersion Direction:", self.direction_combo)

        data_group.setLayout(data_layout)
        layout.addRow(data_group)

        # Overscan configuration
        overscan_group = QGroupBox("Overscan Configuration")
        overscan_layout = QFormLayout()

        self.overscan_start_column_spin = QSpinBox()
        self.overscan_start_column_spin.setRange(1, 10000)
        overscan_layout.addRow("Overscan Start Column (1-based):", self.overscan_start_column_spin)

        self.overscan_method_combo = QComboBox()
        self.overscan_method_combo.addItems(['median', 'polynomial'])
        overscan_layout.addRow("Overscan Method:", self.overscan_method_combo)

        overscan_group.setLayout(overscan_layout)
        layout.addRow(overscan_group)

        # Telescope configuration
        telescope_group = QGroupBox("Telescope Configuration")
        telescope_layout = QFormLayout()

        self.telescope_name_edit = QLineEdit()
        telescope_layout.addRow("Telescope Name:", self.telescope_name_edit)

        self.instrument_edit = QLineEdit()
        telescope_layout.addRow("Instrument Name:", self.instrument_edit)

        telescope_group.setLayout(telescope_layout)
        layout.addRow(telescope_group)

        # Calibration configuration (merged from Telescope & Calibration tab)
        calib_group = QGroupBox("Calibration Configuration")
        calib_layout = QFormLayout()

        self.linelist_combo = QComboBox()
        self.linelist_combo.addItems(['ThAr', 'Ar', 'Ne', 'He', 'Fe', 'custom'])
        calib_layout.addRow("Calibration Lamp Type:", self.linelist_combo)

        self.linelist_path_edit = QLineEdit()
        calib_layout.addRow("Linelist Path:", self.linelist_path_edit)

        self.linelist_file_edit = QLineEdit()
        calib_layout.addRow("Linelist File:", self.linelist_file_edit)

        self.use_precomputed_calib_check = QCheckBox("Use pre-identified calibration files")
        calib_layout.addRow(self.use_precomputed_calib_check)

        self.calibration_path_edit = QLineEdit()
        calib_layout.addRow("Calibration Path:", self.calibration_path_edit)

        self.calibration_file_edit = QLineEdit()
        calib_layout.addRow("Calibration File:", self.calibration_file_edit)

        calib_group.setLayout(calib_layout)
        layout.addRow(calib_group)

        # Detector parameters
        detector_group = QGroupBox("Detector Parameters")
        detector_layout = QFormLayout()

        self.detector_nrows_spin = QSpinBox()
        self.detector_nrows_spin.setRange(1, 10000)
        detector_layout.addRow("Number of Rows:", self.detector_nrows_spin)

        self.detector_ncols_spin = QSpinBox()
        self.detector_ncols_spin.setRange(1, 10000)
        detector_layout.addRow("Number of Columns:", self.detector_ncols_spin)

        self.detector_gain_spin = QDoubleSpinBox()
        self.detector_gain_spin.setRange(0.1, 100.0)
        self.detector_gain_spin.setSingleStep(0.1)
        detector_layout.addRow("CCD Gain:", self.detector_gain_spin)

        self.detector_readnoise_spin = QDoubleSpinBox()
        self.detector_readnoise_spin.setRange(0.0, 100.0)
        self.detector_readnoise_spin.setSingleStep(0.5)
        detector_layout.addRow("Read Noise:", self.detector_readnoise_spin)

        detector_group.setLayout(detector_layout)
        layout.addRow(detector_group)

        widget.setLayout(layout)
        scroll.setWidget(widget)
        return scroll

    def _create_correction_tab(self) -> QWidget:
        """Create bias, flat and background correction settings tab."""
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        widget = QWidget()
        layout = QFormLayout()

        # Bias correction
        bias_group = QGroupBox("Bias Correction")
        bias_layout = QFormLayout()

        self.bias_method_combo = QComboBox()
        self.bias_method_combo.addItems(['combine', 'load'])
        bias_layout.addRow("Bias Method:", self.bias_method_combo)

        self.bias_combine_method_combo = QComboBox()
        self.bias_combine_method_combo.addItems(['mean', 'median', 'sigma_clip'])
        bias_layout.addRow("Combine Method:", self.bias_combine_method_combo)

        self.bias_combine_sigma_spin = QDoubleSpinBox()
        self.bias_combine_sigma_spin.setRange(1.0, 10.0)
        self.bias_combine_sigma_spin.setSingleStep(0.1)
        bias_layout.addRow("Sigma Clipping Threshold:", self.bias_combine_sigma_spin)

        bias_group.setLayout(bias_layout)
        layout.addRow(bias_group)

        # Flat fielding
        flat_group = QGroupBox("Flat Fielding")
        flat_layout = QFormLayout()

        self.flat_method_combo = QComboBox()
        self.flat_method_combo.addItems(['combine', 'load'])
        flat_layout.addRow("Flat Method:", self.flat_method_combo)

        self.flat_combine_method_combo = QComboBox()
        self.flat_combine_method_combo.addItems(['mean', 'median'])
        flat_layout.addRow("Combine Method:", self.flat_combine_method_combo)

        self.flat_q_threshold_spin = QDoubleSpinBox()
        self.flat_q_threshold_spin.setRange(0.0, 1.0)
        self.flat_q_threshold_spin.setSingleStep(0.05)
        flat_layout.addRow("Quality Threshold:", self.flat_q_threshold_spin)

        self.flat_mosaic_maxcount_spin = QDoubleSpinBox()
        self.flat_mosaic_maxcount_spin.setRange(0.0, 100000.0)
        flat_layout.addRow("Mosaic Max Count:", self.flat_mosaic_maxcount_spin)

        flat_group.setLayout(flat_layout)
        layout.addRow(flat_group)

        # Background subtraction
        bg_group = QGroupBox("Background Subtraction")
        bg_layout = QFormLayout()

        self.bg_method_combo = QComboBox()
        self.bg_method_combo.addItems(['2d_poly', 'median'])
        bg_layout.addRow("Background Method:", self.bg_method_combo)

        self.bg_poly_order_spin = QSpinBox()
        self.bg_poly_order_spin.setRange(1, 10)
        bg_layout.addRow("Polynomial Order:", self.bg_poly_order_spin)

        bg_group.setLayout(bg_layout)
        layout.addRow(bg_group)

        widget.setLayout(layout)
        scroll.setWidget(widget)
        return scroll

    def _create_wlcalib_tab(self) -> QWidget:
        """Create wavelength calibration settings tab."""
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        widget = QWidget()
        layout = QFormLayout()

        # Main calibration parameters (merged lamp type with Data & Calibration tab)
        self.wlcalib_search_database_check = QCheckBox("Search calibration database")
        layout.addRow(self.wlcalib_search_database_check)

        self.wlcalib_use_prev_fitpar_check = QCheckBox("Use previous fitted parameters as initial guess")
        layout.addRow(self.wlcalib_use_prev_fitpar_check)

        # Polynomial solution
        poly_group = QGroupBox("Polynomial Solution")
        poly_layout = QFormLayout()

        self.wlcalib_xorder_spin = QSpinBox()
        self.wlcalib_xorder_spin.setRange(1, 10)
        poly_layout.addRow("X Order (dispersion):", self.wlcalib_xorder_spin)

        self.wlcalib_yorder_spin = QSpinBox()
        self.wlcalib_yorder_spin.setRange(1, 10)
        poly_layout.addRow("Y Order (spatial):", self.wlcalib_yorder_spin)

        poly_group.setLayout(poly_layout)
        layout.addRow(poly_group)

        # Quality control
        quality_group = QGroupBox("Quality Control")
        quality_layout = QFormLayout()

        self.wlcalib_rms_threshold_spin = QDoubleSpinBox()
        self.wlcalib_rms_threshold_spin.setRange(0.01, 1.0)
        self.wlcalib_rms_threshold_spin.setSingleStep(0.01)
        quality_layout.addRow("RMS Threshold (Angstrom):", self.wlcalib_rms_threshold_spin)

        quality_group.setLayout(quality_layout)
        layout.addRow(quality_group)

        # Advanced parameters
        adv_group = QGroupBox("Advanced Parameters")
        adv_layout = QFormLayout()

        self.wlcalib_time_diff_spin = QSpinBox()
        self.wlcalib_time_diff_spin.setRange(0, 3600)
        adv_layout.addRow("Time Window (seconds):", self.wlcalib_time_diff_spin)

        self.wlcalib_auto_selection_check = QCheckBox("Auto-select calibration reference")
        adv_layout.addRow(self.wlcalib_auto_selection_check)

        self.wlcalib_window_size_spin = QSpinBox()
        self.wlcalib_window_size_spin.setRange(1, 50)
        adv_layout.addRow("Fitting Window Size:", self.wlcalib_window_size_spin)

        self.wlcalib_clipping_spin = QDoubleSpinBox()
        self.wlcalib_clipping_spin.setRange(1.0, 10.0)
        self.wlcalib_clipping_spin.setSingleStep(0.1)
        adv_layout.addRow("Sigma Clipping Threshold:", self.wlcalib_clipping_spin)

        self.wlcalib_q_threshold_spin = QDoubleSpinBox()
        self.wlcalib_q_threshold_spin.setRange(0.0, 1.0)
        self.wlcalib_q_threshold_spin.setSingleStep(0.05)
        adv_layout.addRow("Line Intensity Threshold:", self.wlcalib_q_threshold_spin)

        adv_group.setLayout(adv_layout)
        layout.addRow(adv_group)

        widget.setLayout(layout)
        scroll.setWidget(widget)
        return scroll

    def _create_extraction_tab(self) -> QWidget:
        """Create order tracing and spectrum extraction settings tab."""
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        widget = QWidget()
        layout = QFormLayout()

        # Order tracing
        trace_group = QGroupBox("Order Tracing")
        trace_layout = QFormLayout()

        self.trace_scan_step_spin = QSpinBox()
        self.trace_scan_step_spin.setRange(1, 100)
        trace_layout.addRow("Scan Step (pixels):", self.trace_scan_step_spin)

        self.trace_minimum_spin = QDoubleSpinBox()
        self.trace_minimum_spin.setRange(0.0, 1000.0)
        trace_layout.addRow("Minimum Detection Strength:", self.trace_minimum_spin)

        self.trace_separation_spin = QDoubleSpinBox()
        self.trace_separation_spin.setRange(0.0, 100.0)
        trace_layout.addRow("Order Separation (pixels):", self.trace_separation_spin)

        self.trace_filling_spin = QDoubleSpinBox()
        self.trace_filling_spin.setRange(0.0, 1.0)
        self.trace_filling_spin.setSingleStep(0.05)
        trace_layout.addRow("Filling Factor Threshold:", self.trace_filling_spin)

        self.trace_degree_spin = QSpinBox()
        self.trace_degree_spin.setRange(1, 10)
        trace_layout.addRow("Polynomial Degree:", self.trace_degree_spin)

        trace_group.setLayout(trace_layout)
        layout.addRow(trace_group)

        # Spectrum extraction
        extract_group = QGroupBox("Spectrum Extraction")
        extract_layout = QFormLayout()

        self.extract_method_combo = QComboBox()
        self.extract_method_combo.addItems(['sum', 'optimal'])
        extract_layout.addRow("Extraction Method:", self.extract_method_combo)

        self.extract_lower_limit_spin = QDoubleSpinBox()
        self.extract_lower_limit_spin.setRange(-100.0, 0.0)
        self.extract_lower_limit_spin.setSingleStep(0.5)
        extract_layout.addRow("Lower Limit (pixels):", self.extract_lower_limit_spin)

        self.extract_upper_limit_spin = QDoubleSpinBox()
        self.extract_upper_limit_spin.setRange(0.0, 100.0)
        self.extract_upper_limit_spin.setSingleStep(0.5)
        extract_layout.addRow("Upper Limit (pixels):", self.extract_upper_limit_spin)

        extract_group.setLayout(extract_layout)
        layout.addRow(extract_group)

        widget.setLayout(layout)
        scroll.setWidget(widget)
        return scroll

    def _create_processing_tab(self) -> QWidget:
        """Create general processing and output settings tab."""
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        widget = QWidget()
        layout = QFormLayout()

        # Processing mode
        self.mode_combo = QComboBox()
        self.mode_combo.addItems(['normal', 'debug'])
        layout.addRow("Processing Mode:", self.mode_combo)

        # Output path
        self.output_path_edit = QLineEdit()
        layout.addRow("Output Path:", self.output_path_edit)

        # Figure format
        self.fig_format_combo = QComboBox()
        self.fig_format_combo.addItems(['png', 'pdf', 'jpg'])
        layout.addRow("Figure Format:", self.fig_format_combo)

        # Output suffix
        self.oned_suffix_edit = QLineEdit()
        layout.addRow("Output Suffix:", self.oned_suffix_edit)

        # Processing cores
        self.ncores_edit = QLineEdit()
        layout.addRow("Number of Cores:", self.ncores_edit)

        # Auto process
        self.auto_process_check = QCheckBox("Auto process all science images without confirmation")
        layout.addRow(self.auto_process_check)

        # Cosmic ray correction
        cosmic_group = QGroupBox("Cosmic Ray Correction")
        cosmic_layout = QFormLayout()

        self.cosmic_enabled_check = QCheckBox("Enable cosmic ray correction")
        cosmic_layout.addRow(self.cosmic_enabled_check)

        self.cosmic_sigma_spin = QDoubleSpinBox()
        self.cosmic_sigma_spin.setRange(0.1, 20.0)
        self.cosmic_sigma_spin.setSingleStep(0.1)
        cosmic_layout.addRow("Sigma Threshold:", self.cosmic_sigma_spin)

        self.cosmic_window_spin = QSpinBox()
        self.cosmic_window_spin.setRange(1, 50)
        cosmic_layout.addRow("Window Size:", self.cosmic_window_spin)

        cosmic_group.setLayout(cosmic_layout)
        layout.addRow(cosmic_group)

        # Output options
        output_group = QGroupBox("Output Options")
        output_layout = QFormLayout()

        self.save_plots_check = QCheckBox("Save diagnostic plots")
        output_layout.addRow(self.save_plots_check)

        self.save_overscan_check = QCheckBox("Save overscan correction results")
        output_layout.addRow(self.save_overscan_check)

        self.save_bias_check = QCheckBox("Save bias correction results")
        output_layout.addRow(self.save_bias_check)

        self.save_flat_check = QCheckBox("Save flat fielding results")
        output_layout.addRow(self.save_flat_check)

        self.save_background_check = QCheckBox("Save background subtraction results")
        output_layout.addRow(self.save_background_check)

        self.save_cosmic_check = QCheckBox("Save cosmic ray correction results")
        output_layout.addRow(self.save_cosmic_check)

        self.save_extraction_check = QCheckBox("Save spectrum extraction results")
        output_layout.addRow(self.save_extraction_check)

        self.save_wlcalib_check = QCheckBox("Save wavelength calibration results")
        output_layout.addRow(self.save_wlcalib_check)

        self.save_deblaze_check = QCheckBox("Save de-blazing results")
        output_layout.addRow(self.save_deblaze_check)

        output_group.setLayout(output_layout)
        layout.addRow(output_group)

        widget.setLayout(layout)
        scroll.setWidget(widget)
        return scroll

    def load_current_values(self):
        """Load current configuration values into UI."""
        # Data & Calibration tab
        self.rawdata_path_edit.setText(self.config.get('data', 'rawdata_path', ''))
        self.statime_key_edit.setText(self.config.get('data', 'statime_key', 'DATE-OBS'))
        self.exptime_key_edit.setText(self.config.get('data', 'exptime_key', 'EXPTIME'))
        self.direction_combo.setCurrentText(self.config.get('data', 'direction', 'xr-'))
        self.overscan_start_column_spin.setValue(self.config.get_int('data', 'overscan_start_column', 4097))
        self.overscan_method_combo.setCurrentText(self.config.get('data', 'overscan_method', 'median'))
        self.telescope_name_edit.setText(self.config.get('telescope', 'name', 'xinglong216hrs'))
        self.instrument_edit.setText(self.config.get('telescope', 'instrument', 'hrs'))
        self.linelist_combo.setCurrentText(self.config.get('telescope.linelist', 'linelist_type', 'ThAr'))
        self.linelist_path_edit.setText(self.config.get('telescope.linelist', 'linelist_path', 'calib_data/linelists/'))
        self.linelist_file_edit.setText(self.config.get('telescope.linelist', 'linelist_file', 'thar-noao.dat'))
        self.use_precomputed_calib_check.setChecked(self.config.get_bool('telescope.linelist', 'use_precomputed_calibration', False))
        self.calibration_path_edit.setText(self.config.get('telescope.linelist', 'calibration_path', 'calib_data/telescopes/xinglong216hrs/'))
        self.calibration_file_edit.setText(self.config.get('telescope.linelist', 'calibration_file', 'wlcalib_20211123011_A.fits'))
        self.detector_nrows_spin.setValue(self.config.get_int('telescope.detector', 'nrows', 2048))
        self.detector_ncols_spin.setValue(self.config.get_int('telescope.detector', 'ncols', 4096))
        self.detector_gain_spin.setValue(self.config.get_float('telescope.detector', 'gain', 1.0))
        self.detector_readnoise_spin.setValue(self.config.get_float('telescope.detector', 'readnoise', 5.0))

        # Bias & Flat & Background tab
        self.bias_method_combo.setCurrentText(self.config.get('reduce.bias', 'method', 'combine'))
        self.bias_combine_method_combo.setCurrentText(self.config.get('reduce.bias', 'combine_method', 'median'))
        self.bias_combine_sigma_spin.setValue(self.config.get_float('reduce.bias', 'combine_sigma', 3.0))
        self.flat_method_combo.setCurrentText(self.config.get('reduce.flat', 'method', 'combine'))
        self.flat_combine_method_combo.setCurrentText(self.config.get('reduce.flat', 'combine_method', 'median'))
        self.flat_q_threshold_spin.setValue(self.config.get_float('reduce.flat', 'q_threshold', 0.5))
        self.flat_mosaic_maxcount_spin.setValue(self.config.get_float('reduce.flat', 'mosaic_maxcount', 65535.0))
        self.bg_method_combo.setCurrentText(self.config.get('reduce.background', 'method', '2d_poly'))
        self.bg_poly_order_spin.setValue(self.config.get_int('reduce.background', 'poly_order', 2))

        # Wavelength Calibration tab
        self.wlcalib_search_database_check.setChecked(self.config.get_bool('reduce.wlcalib', 'search_database', True))
        self.wlcalib_use_prev_fitpar_check.setChecked(self.config.get_bool('reduce.wlcalib', 'use_prev_fitpar', True))
        self.wlcalib_xorder_spin.setValue(self.config.get_int('reduce.wlcalib', 'xorder', 4))
        self.wlcalib_yorder_spin.setValue(self.config.get_int('reduce.wlcalib', 'yorder', 4))
        self.wlcalib_rms_threshold_spin.setValue(self.config.get_float('reduce.wlcalib', 'rms_threshold', 0.1))
        self.wlcalib_time_diff_spin.setValue(self.config.get_int('reduce.wlcalib', 'time_diff', 600))
        self.wlcalib_auto_selection_check.setChecked(self.config.get_bool('reduce.wlcalib', 'auto_selection', True))
        self.wlcalib_window_size_spin.setValue(self.config.get_int('reduce.wlcalib', 'window_size', 10))
        self.wlcalib_clipping_spin.setValue(self.config.get_float('reduce.wlcalib', 'clipping', 3.0))
        self.wlcalib_q_threshold_spin.setValue(self.config.get_float('reduce.wlcalib', 'q_threshold', 0.5))

        # Trace & Extraction tab
        self.trace_scan_step_spin.setValue(self.config.get_int('reduce.trace', 'scan_step', 10))
        self.trace_minimum_spin.setValue(self.config.get_float('reduce.trace', 'minimum', 50.0))
        self.trace_separation_spin.setValue(self.config.get_float('reduce.trace', 'separation', 30.0))
        self.trace_filling_spin.setValue(self.config.get_float('reduce.trace', 'filling', 0.3))
        self.trace_degree_spin.setValue(self.config.get_int('reduce.trace', 'degree', 3))
        self.extract_method_combo.setCurrentText(self.config.get('reduce.extract', 'method', 'sum'))
        self.extract_lower_limit_spin.setValue(self.config.get_float('reduce.extract', 'lower_limit', -5.0))
        self.extract_upper_limit_spin.setValue(self.config.get_float('reduce.extract', 'upper_limit', 5.0))

        # Processing & Output tab
        self.mode_combo.setCurrentText(self.config.get('reduce', 'mode', 'normal'))
        self.output_path_edit.setText(self.config.get('reduce', 'output_path', 'output'))
        self.fig_format_combo.setCurrentText(self.config.get('reduce', 'fig_format', 'png'))
        self.oned_suffix_edit.setText(self.config.get('reduce', 'oned_suffix', '_ods'))
        self.ncores_edit.setText(self.config.get('reduce', 'ncores', 'max'))
        self.auto_process_check.setChecked(self.config.get_bool('reduce', 'auto_process_all_images', False))
        self.cosmic_enabled_check.setChecked(self.config.get_bool('reduce', 'cosmic_enabled', True))
        self.cosmic_sigma_spin.setValue(self.config.get_float('reduce', 'cosmic_sigma', 5.0))
        self.cosmic_window_spin.setValue(self.config.get_int('reduce', 'cosmic_window', 7))
        self.save_plots_check.setChecked(self.config.get_bool('reduce', 'save_plots', True))
        self.save_overscan_check.setChecked(self.config.get_bool('reduce.save_intermediate', 'save_overscan', True))
        self.save_bias_check.setChecked(self.config.get_bool('reduce.save_intermediate', 'save_bias', True))
        self.save_flat_check.setChecked(self.config.get_bool('reduce.save_intermediate', 'save_flat', True))
        self.save_background_check.setChecked(self.config.get_bool('reduce.save_intermediate', 'save_background', True))
        self.save_cosmic_check.setChecked(self.config.get_bool('reduce.save_intermediate', 'save_cosmic', True))
        self.save_extraction_check.setChecked(self.config.get_bool('reduce.save_intermediate', 'save_extraction', True))
        self.save_wlcalib_check.setChecked(self.config.get_bool('reduce.save_intermediate', 'save_wlcalib', True))
        self.save_deblaze_check.setChecked(self.config.get_bool('reduce.save_intermediate', 'save_deblaze', True))

    def _save_settings(self):
        """Save settings to configuration."""
        try:
            # Data & Calibration tab
            self.config.set('data', 'rawdata_path', self.rawdata_path_edit.text())
            self.config.set('data', 'statime_key', self.statime_key_edit.text())
            self.config.set('data', 'exptime_key', self.exptime_key_edit.text())
            self.config.set('data', 'direction', self.direction_combo.currentText())
            self.config.set('data', 'overscan_start_column', str(self.overscan_start_column_spin.value()))
            self.config.set('data', 'overscan_method', self.overscan_method_combo.currentText())
            self.config.set('telescope', 'name', self.telescope_name_edit.text())
            self.config.set('telescope', 'instrument', self.instrument_edit.text())
            self.config.set('telescope.linelist', 'linelist_type', self.linelist_combo.currentText())
            self.config.set('telescope.linelist', 'linelist_path', self.linelist_path_edit.text())
            self.config.set('telescope.linelist', 'linelist_file', self.linelist_file_edit.text())
            self.config.set('telescope.linelist', 'use_precomputed_calibration',
                           'yes' if self.use_precomputed_calib_check.isChecked() else 'no')
            self.config.set('telescope.linelist', 'calibration_path', self.calibration_path_edit.text())
            self.config.set('telescope.linelist', 'calibration_file', self.calibration_file_edit.text())
            self.config.set('telescope.detector', 'nrows', str(self.detector_nrows_spin.value()))
            self.config.set('telescope.detector', 'ncols', str(self.detector_ncols_spin.value()))
            self.config.set('telescope.detector', 'gain', str(self.detector_gain_spin.value()))
            self.config.set('telescope.detector', 'readnoise', str(self.detector_readnoise_spin.value()))

            # Bias & Flat & Background tab
            self.config.set('reduce.bias', 'method', self.bias_method_combo.currentText())
            self.config.set('reduce.bias', 'combine_method', self.bias_combine_method_combo.currentText())
            self.config.set('reduce.bias', 'combine_sigma', str(self.bias_combine_sigma_spin.value()))
            self.config.set('reduce.flat', 'method', self.flat_method_combo.currentText())
            self.config.set('reduce.flat', 'combine_method', self.flat_combine_method_combo.currentText())
            self.config.set('reduce.flat', 'q_threshold', str(self.flat_q_threshold_spin.value()))
            self.config.set('reduce.flat', 'mosaic_maxcount', str(self.flat_mosaic_maxcount_spin.value()))
            self.config.set('reduce.background', 'method', self.bg_method_combo.currentText())
            self.config.set('reduce.background', 'poly_order', str(self.bg_poly_order_spin.value()))

            # Wavelength Calibration tab
            self.config.set('reduce.wlcalib', 'search_database',
                           'yes' if self.wlcalib_search_database_check.isChecked() else 'no')
            self.config.set('reduce.wlcalib', 'use_prev_fitpar',
                           'yes' if self.wlcalib_use_prev_fitpar_check.isChecked() else 'no')
            self.config.set('reduce.wlcalib', 'xorder', str(self.wlcalib_xorder_spin.value()))
            self.config.set('reduce.wlcalib', 'yorder', str(self.wlcalib_yorder_spin.value()))
            self.config.set('reduce.wlcalib', 'rms_threshold', str(self.wlcalib_rms_threshold_spin.value()))
            self.config.set('reduce.wlcalib', 'time_diff', str(self.wlcalib_time_diff_spin.value()))
            self.config.set('reduce.wlcalib', 'auto_selection',
                           'yes' if self.wlcalib_auto_selection_check.isChecked() else 'no')
            self.config.set('reduce.wlcalib', 'window_size', str(self.wlcalib_window_size_spin.value()))
            self.config.set('reduce.wlcalib', 'clipping', str(self.wlcalib_clipping_spin.value()))
            self.config.set('reduce.wlcalib', 'q_threshold', str(self.wlcalib_q_threshold_spin.value()))

            # Trace & Extraction tab
            self.config.set('reduce.trace', 'scan_step', str(self.trace_scan_step_spin.value()))
            self.config.set('reduce.trace', 'minimum', str(self.trace_minimum_spin.value()))
            self.config.set('reduce.trace', 'separation', str(self.trace_separation_spin.value()))
            self.config.set('reduce.trace', 'filling', str(self.trace_filling_spin.value()))
            self.config.set('reduce.trace', 'degree', str(self.trace_degree_spin.value()))
            self.config.set('reduce.extract', 'method', self.extract_method_combo.currentText())
            self.config.set('reduce.extract', 'lower_limit', str(self.extract_lower_limit_spin.value()))
            self.config.set('reduce.extract', 'upper_limit', str(self.extract_upper_limit_spin.value()))

            # Processing & Output tab
            self.config.set('reduce', 'mode', self.mode_combo.currentText())
            self.config.set('reduce', 'output_path', self.output_path_edit.text())
            self.config.set('reduce', 'fig_format', self.fig_format_combo.currentText())
            self.config.set('reduce', 'oned_suffix', self.oned_suffix_edit.text())
            self.config.set('reduce', 'ncores', self.ncores_edit.text())
            self.config.set('reduce', 'auto_process_all_images', 'yes' if self.auto_process_check.isChecked() else 'no')
            self.config.set('reduce', 'cosmic_enabled', 'yes' if self.cosmic_enabled_check.isChecked() else 'no')
            self.config.set('reduce', 'cosmic_sigma', str(self.cosmic_sigma_spin.value()))
            self.config.set('reduce', 'cosmic_window', str(self.cosmic_window_spin.value()))
            self.config.set('reduce', 'save_plots', 'yes' if self.save_plots_check.isChecked() else 'no')
            self.config.set('reduce.save_intermediate', 'save_overscan',
                           'yes' if self.save_overscan_check.isChecked() else 'no')
            self.config.set('reduce.save_intermediate', 'save_bias',
                           'yes' if self.save_bias_check.isChecked() else 'no')
            self.config.set('reduce.save_intermediate', 'save_flat',
                           'yes' if self.save_flat_check.isChecked() else 'no')
            self.config.set('reduce.save_intermediate', 'save_background',
                           'yes' if self.save_background_check.isChecked() else 'no')
            self.config.set('reduce.save_intermediate', 'save_cosmic',
                           'yes' if self.save_cosmic_check.isChecked() else 'no')
            self.config.set('reduce.save_intermediate', 'save_extraction',
                           'yes' if self.save_extraction_check.isChecked() else 'no')
            self.config.set('reduce.save_intermediate', 'save_wlcalib',
                           'yes' if self.save_wlcalib_check.isChecked() else 'no')
            self.config.set('reduce.save_intermediate', 'save_deblaze',
                           'yes' if self.save_deblaze_check.isChecked() else 'no')

            # Save to file
            self.config.save()

            self.accept()

        except Exception as e:
            from PyQt5.QtWidgets import QMessageBox
            QMessageBox.critical(self, "Error", f"Failed to save settings: {e}")