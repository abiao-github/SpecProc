
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

    def __init__(self, config: ConfigManager, parent=None, initial_tab: int = 0):
        """
        Initialize settings dialog.

        Args:
            config: Configuration manager instance
            parent: Parent widget
            initial_tab: Initial tab index to show (0: Data & Instruments, 1: Data Reduction, etc.)
        """
        super().__init__(parent)
        self.config = config
        self.setWindowTitle("SpecProc Settings")
        self.setModal(True)
        self.resize(800, 600)
        
        self.initial_tab = initial_tab
        self.tab_widget = None

        self.init_ui()
        self.load_current_values()
        
        # Set initial tab
        if self.tab_widget:
            self.tab_widget.setCurrentIndex(initial_tab)

    def init_ui(self):
        """Initialize user interface."""
        layout = QVBoxLayout()

        # Tab widget with different configuration categories
        self.tab_widget = QTabWidget()

        # Create tabs following data reduction pipeline order
        self.tab_widget.addTab(self._create_data_tab(), "Data & Instruments")
        self.tab_widget.addTab(self._create_reduction_tab(), "Data Reduction")
<<<<<<< HEAD
        self.tab_widget.addTab(self._create_extraction_tab(), "Trace & Extraction")
=======
        self.tab_widget.addTab(self._create_extraction_tab(), "Orders Tracing")
>>>>>>> cef6f04 (	modified:   README.md)
        self.tab_widget.addTab(self._create_wlcalib_tab(), "Wavelength Calibration")
        self.tab_widget.addTab(self._create_processing_tab(), "Processing & Output")

        layout.addWidget(self.tab_widget)

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
        """Create data input and instrument settings tab."""
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

        self.output_path_edit = QLineEdit()
        data_layout.addRow("Output Path:", self.output_path_edit)

        data_group.setLayout(data_layout)
        layout.addRow(data_group)

        # Telescope configuration
        telescope_group = QGroupBox("Telescope & Instrument")
        telescope_layout = QFormLayout()

        self.telescope_name_edit = QLineEdit()
        telescope_layout.addRow("Telescope Name:", self.telescope_name_edit)

        self.instrument_edit = QLineEdit()
        telescope_layout.addRow("Instrument Name:", self.instrument_edit)

        telescope_group.setLayout(telescope_layout)
        layout.addRow(telescope_group)

        # Detector parameters
        detector_group = QGroupBox("Detector Parameters")
        detector_layout = QFormLayout()

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
<<<<<<< HEAD
=======
        widget.setMinimumWidth(400)
>>>>>>> cef6f04 (	modified:   README.md)
        scroll.setWidget(widget)
        return scroll

    def _create_reduction_tab(self) -> QWidget:
        """Create bias, flat, background and cosmic ray correction settings tab."""
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        widget = QWidget()
        layout = QFormLayout()

        # Overscan configuration
        overscan_group = QGroupBox("Overscan Correction")
        overscan_layout = QFormLayout()

        self.overscan_start_column_spin = QSpinBox()
        self.overscan_start_column_spin.setRange(-1, 10000)
        self.overscan_start_column_spin.setSpecialValueText("None")
        overscan_layout.addRow("Overscan Start Column (1-based):", self.overscan_start_column_spin)

        self.overscan_method_combo = QComboBox()
        self.overscan_method_combo.addItems(['mean_only', 'mean_savgol', 'mean_polynomial'])
        self.overscan_method_combo.setToolTip(
            "mean_only: Sigma-clipped mean only, no smoothing\n"
            "mean_savgol: Sigma-clipped mean with Savitzky-Golay smoothing\n"
            "mean_polynomial: Sigma-clipped mean with polynomial fitting (IRAF style)"
        )
        overscan_layout.addRow("Overscan Method:", self.overscan_method_combo)

        self.overscan_smooth_window_spin = QSpinBox()
        self.overscan_smooth_window_spin.setRange(-1, 10000)
        self.overscan_smooth_window_spin.setSpecialValueText("Auto (1/5 height)")
        self.overscan_smooth_window_spin.setToolTip("Savitzky-Golay smoothing window size. -1 for auto (1/5 of image height)")
        overscan_layout.addRow("Smooth Window:", self.overscan_smooth_window_spin)

        self.overscan_poly_type_combo = QComboBox()
        self.overscan_poly_type_combo.addItems(['legendre', 'chebyshev', 'polynomial'])
        self.overscan_poly_type_combo.setToolTip(
            "Polynomial type for mean_polynomial method:\n"
            "legendre: Legendre polynomial (IRAF default)\n"
            "chebyshev: Chebyshev polynomial\n"
            "polynomial: Standard power series"
        )
        overscan_layout.addRow("Poly Type:", self.overscan_poly_type_combo)

        self.overscan_poly_order_spin = QSpinBox()
        self.overscan_poly_order_spin.setRange(1, 10)
        self.overscan_poly_order_spin.setValue(3)
        self.overscan_poly_order_spin.setToolTip("Polynomial order (degree) for mean_polynomial method")
        overscan_layout.addRow("Poly Order:", self.overscan_poly_order_spin)

<<<<<<< HEAD
=======
        self.detector_split_row_spin = QSpinBox()
        self.detector_split_row_spin.setRange(-1, 20000)
        self.detector_split_row_spin.setSpecialValueText("Auto (rows/2)")
        self.detector_split_row_spin.setToolTip(
            "Detector amplifier split row (1-based).\n"
            "Divides the CCD into lower and upper readout halves.\n"
            "Used by overscan correction and step 3 background fitting.\n"
            "For Xinglong 216 HRS 4136-row CCD: 2068. Set -1 for auto."
        )
        overscan_layout.addRow("Amplifier Split Row:", self.detector_split_row_spin)

>>>>>>> cef6f04 (	modified:   README.md)
        overscan_group.setLayout(overscan_layout)
        layout.addRow(overscan_group)

        # Connect overscan method change to update parameter states
        self.overscan_method_combo.currentTextChanged.connect(self._update_overscan_params_enabled)

        # Bias correction
        bias_group = QGroupBox("Bias Correction")
        bias_layout = QFormLayout()

        self.bias_combine_method_combo = QComboBox()
        self.bias_combine_method_combo.addItems(['mean', 'median', 'sigma_clip'])
        bias_layout.addRow("Combine Method:", self.bias_combine_method_combo)

        self.bias_combine_sigma_spin = QDoubleSpinBox()
        self.bias_combine_sigma_spin.setRange(1.0, 10.0)
        self.bias_combine_sigma_spin.setSingleStep(0.1)
        bias_layout.addRow("Sigma Clipping Threshold:", self.bias_combine_sigma_spin)

        bias_group.setLayout(bias_layout)
        layout.addRow(bias_group)

<<<<<<< HEAD
        # Flat fielding
        flat_group = QGroupBox("Flat Fielding")
        flat_layout = QFormLayout()

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

        self.flat_blaze_smooth_method_combo = QComboBox()
        self.flat_blaze_smooth_method_combo.addItems(['median', 'savgol', 'bspline'])
        flat_layout.addRow("Blaze Smooth Method:", self.flat_blaze_smooth_method_combo)

        self.flat_blaze_smooth_window_spin = QSpinBox()
        self.flat_blaze_smooth_window_spin.setRange(5, 401)
        self.flat_blaze_smooth_window_spin.setSingleStep(2)
        flat_layout.addRow("Blaze Smooth Window:", self.flat_blaze_smooth_window_spin)

        self.flat_blaze_bspline_smooth_spin = QDoubleSpinBox()
        self.flat_blaze_bspline_smooth_spin.setRange(0.0, 20.0)
        self.flat_blaze_bspline_smooth_spin.setSingleStep(0.1)
        flat_layout.addRow("Blaze B-spline Smooth:", self.flat_blaze_bspline_smooth_spin)

        self.flat_width_smooth_window_spin = QSpinBox()
        self.flat_width_smooth_window_spin.setRange(5, 401)
        self.flat_width_smooth_window_spin.setSingleStep(2)
        flat_layout.addRow("Width Smooth Window:", self.flat_width_smooth_window_spin)

        self.flat_profile_bin_step_spin = QDoubleSpinBox()
        self.flat_profile_bin_step_spin.setRange(0.002, 0.2)
        self.flat_profile_bin_step_spin.setDecimals(3)
        self.flat_profile_bin_step_spin.setSingleStep(0.001)
        flat_layout.addRow("y_norm Bin Step:", self.flat_profile_bin_step_spin)

        self.flat_pixel_min_spin = QDoubleSpinBox()
        self.flat_pixel_min_spin.setRange(0.01, 2.0)
        self.flat_pixel_min_spin.setSingleStep(0.01)
        flat_layout.addRow("Pixel Flat Min:", self.flat_pixel_min_spin)

        self.flat_pixel_max_spin = QDoubleSpinBox()
        self.flat_pixel_max_spin.setRange(0.1, 5.0)
        self.flat_pixel_max_spin.setSingleStep(0.01)
        flat_layout.addRow("Pixel Flat Max:", self.flat_pixel_max_spin)

        flat_group.setLayout(flat_layout)
        layout.addRow(flat_group)

        # Background subtraction
=======
        # Cosmic ray correction (Step 1)
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

        # Background subtraction (Step 3)
>>>>>>> cef6f04 (	modified:   README.md)
        bg_group = QGroupBox("Background Subtraction")
        bg_layout = QFormLayout()

        self.bg_method_combo = QComboBox()
<<<<<<< HEAD
        self.bg_method_combo.addItems(['chebyshev', 'bspline', 'smooth'])
        bg_layout.addRow("Background Method:", self.bg_method_combo)
=======
        self.bg_method_combo.addItems(['convolution', 'column_spline', 'chebyshev', 'bspline', 'gaussian_smooth'])
        self.bg_method_combo.setToolTip(
            "'convolution': astropy NaN-aware 2D Gaussian convolution (recommended,\n"
            "  reproduces local ripples/fringes);\n"
            "'column_spline': per-column 1D smoothing spline (follows ripple peaks/troughs,\n"
            "  with light horizontal smoothing);\n"
            "'chebyshev': global 2D Chebyshev polynomial fit;\n"
            "'bspline': bivariate spline fit;\n"
            "'gaussian_smooth': iterative Gaussian smoothing with valid-pixel weighting.")
        bg_layout.addRow("Background Method:", self.bg_method_combo)
        self.bg_method_combo.currentTextChanged.connect(self._update_bg_params_enabled)

        self.bg_kernel_sigma_x_spin = QDoubleSpinBox()
        self.bg_kernel_sigma_x_spin.setRange(5.0, 100.0)
        self.bg_kernel_sigma_x_spin.setSingleStep(1.0)
        self.bg_kernel_sigma_x_spin.setValue(13.0)
        self.bg_kernel_sigma_x_spin.setToolTip(
            "Gaussian kernel σ_x (色散方向, 水平, X) for convolution method.\n"
            "建议 10–20 px。\n"
            "σ_x 应极细以还原干涉波峰（如 60-110 px），较小的值（13）能完美刻画横向条纹。")
        bg_layout.addRow("Kernel σ_x (px):", self.bg_kernel_sigma_x_spin)

        self.bg_kernel_sigma_y_spin = QDoubleSpinBox()
        self.bg_kernel_sigma_y_spin.setRange(5.0, 80.0)
        self.bg_kernel_sigma_y_spin.setSingleStep(1.0)
        self.bg_kernel_sigma_y_spin.setValue(13.0)
        self.bg_kernel_sigma_y_spin.setToolTip(
            "Gaussian kernel σ_y (交叉色散方向, 垂直, Y) for convolution method.\n"
            "建议 12–20 px。\n"
            "σ_y (13) 刚好能跨越 mask 宽度（10-12 px），同时保留纵向高频变化。")
        bg_layout.addRow("Kernel σ_y (px):", self.bg_kernel_sigma_y_spin)

        self.bg_spline_smooth_factor_spin = QDoubleSpinBox()
        self.bg_spline_smooth_factor_spin.setRange(0.01, 10.0)
        self.bg_spline_smooth_factor_spin.setSingleStep(0.05)
        self.bg_spline_smooth_factor_spin.setDecimals(3)
        self.bg_spline_smooth_factor_spin.setToolTip(
            "Smoothing factor multiplier for column_spline method.\n"
            "s = factor × N_bg_pixels.  Smaller values let the spline follow\n"
            "ripple peaks/troughs more closely.  Typical 0.1–1.0.")
        bg_layout.addRow("Spline Smooth Factor:", self.bg_spline_smooth_factor_spin)

        self.bg_spline_post_smooth_x_spin = QDoubleSpinBox()
        self.bg_spline_post_smooth_x_spin.setRange(0.0, 50.0)
        self.bg_spline_post_smooth_x_spin.setSingleStep(1.0)
        self.bg_spline_post_smooth_x_spin.setToolTip(
            "Gaussian σ (pixels) for light horizontal smoothing after column-wise\n"
            "spline assembly.  Removes column-to-column discontinuities.\n"
            "Set to 0 to disable.  Typical 3–10 px.")
        bg_layout.addRow("Spline Post-smooth σ_x:", self.bg_spline_post_smooth_x_spin)
>>>>>>> cef6f04 (	modified:   README.md)

        self.bg_poly_order_spin = QSpinBox()
        self.bg_poly_order_spin.setRange(1, 10)
        bg_layout.addRow("Fit Order:", self.bg_poly_order_spin)

        self.bg_smooth_sigma_spin = QDoubleSpinBox()
        self.bg_smooth_sigma_spin.setRange(1.0, 200.0)
        self.bg_smooth_sigma_spin.setSingleStep(1.0)
        bg_layout.addRow("Smooth Sigma:", self.bg_smooth_sigma_spin)

        self.bg_sigma_clip_spin = QDoubleSpinBox()
        self.bg_sigma_clip_spin.setRange(0.0, 10.0)
        self.bg_sigma_clip_spin.setSingleStep(0.1)
        bg_layout.addRow("Sigma Clip:", self.bg_sigma_clip_spin)

        self.bg_sigma_clip_iter_spin = QSpinBox()
        self.bg_sigma_clip_iter_spin.setRange(1, 20)
        bg_layout.addRow("Clip Iterations:", self.bg_sigma_clip_iter_spin)

<<<<<<< HEAD
=======
        self.bg_clip_mode_combo = QComboBox()
        self.bg_clip_mode_combo.addItems(['upper', 'both', 'lower'])
        self.bg_clip_mode_combo.setToolTip(
            "'upper': reject only bright outliers (recommended for flats with "
            "unmasked faint orders);\n"
            "'both': symmetric rejection;\n"
            "'lower': reject only faint outliers.")
        bg_layout.addRow("Clip Mode:", self.bg_clip_mode_combo)

>>>>>>> cef6f04 (	modified:   README.md)
        self.bg_mask_margin_spin = QSpinBox()
        self.bg_mask_margin_spin.setRange(0, 20)
        bg_layout.addRow("Mask Margin (px):", self.bg_mask_margin_spin)

        self.bg_bspline_smooth_spin = QDoubleSpinBox()
        self.bg_bspline_smooth_spin.setRange(0.0, 100.0)
        self.bg_bspline_smooth_spin.setSingleStep(0.1)
        bg_layout.addRow("B-spline Smooth:", self.bg_bspline_smooth_spin)

        bg_group.setLayout(bg_layout)
        layout.addRow(bg_group)

<<<<<<< HEAD
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
=======
        # Flat fielding (Step 4)
        flat_group = QGroupBox("Flat Fielding")
        flat_layout = QFormLayout()

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

        self.flat_blaze_knot_spacing_spin = QSpinBox()
        self.flat_blaze_knot_spacing_spin.setRange(50, 2000)
        self.flat_blaze_knot_spacing_spin.setSingleStep(50)
        self.flat_blaze_knot_spacing_spin.setToolTip(
            "Interior B-spline knot spacing (pixels) for blaze fitting.\n"
            "Wider = more rigid (better fringe rejection).\n"
            "Must exceed the fringe period (typically >200 px).\n"
            "Recommended 400-600 for 4k detectors.")
        flat_layout.addRow("Blaze Knot Spacing:", self.flat_blaze_knot_spacing_spin)

        self.flat_blaze_edge_nknots_spin = QSpinBox()
        self.flat_blaze_edge_nknots_spin.setRange(2, 20)
        self.flat_blaze_edge_nknots_spin.setSingleStep(1)
        self.flat_blaze_edge_nknots_spin.setToolTip(
            "Number of extra dense knots at each edge of the blaze.\n"
            "More knots = closer fit to the steep edge roll-off.\n"
            "4-8 recommended.")
        flat_layout.addRow("Blaze Edge Knots:", self.flat_blaze_edge_nknots_spin)

        self.flat_width_smooth_window_spin = QSpinBox()
        self.flat_width_smooth_window_spin.setRange(5, 401)
        self.flat_width_smooth_window_spin.setSingleStep(2)
        flat_layout.addRow("Width Smooth Window:", self.flat_width_smooth_window_spin)

        self.flat_n_profile_segments_spin = QSpinBox()
        self.flat_n_profile_segments_spin.setRange(1, 200)
        self.flat_n_profile_segments_spin.setSingleStep(10)
        self.flat_n_profile_segments_spin.setToolTip(
            "Number of segments along X (dispersion) for the\n"
            "cross-dispersion profile. Each segment builds a median\n"
            "binned profile, then all segments are Gaussian-smoothed\n"
            "along X to produce a continuously varying P(y_norm; x).\n"
            "50-100 recommended for 4k detectors.")
        flat_layout.addRow("Profile Segments:", self.flat_n_profile_segments_spin)

        self.flat_profile_smooth_sigma_spin = QDoubleSpinBox()
        self.flat_profile_smooth_sigma_spin.setRange(0.0, 50.0)
        self.flat_profile_smooth_sigma_spin.setSingleStep(1.0)
        self.flat_profile_smooth_sigma_spin.setToolTip(
            "Gaussian sigma (in segment units) for smoothing the\n"
            "cross-dispersion profiles along the dispersion direction.\n"
            "0 = no smoothing. 4-8 typical for 100 segments.")
        flat_layout.addRow("Profile Smooth Sigma:", self.flat_profile_smooth_sigma_spin)

        self.flat_profile_bin_step_spin = QDoubleSpinBox()
        self.flat_profile_bin_step_spin.setRange(0.002, 0.2)
        self.flat_profile_bin_step_spin.setDecimals(3)
        self.flat_profile_bin_step_spin.setSingleStep(0.001)
        flat_layout.addRow("y_norm Bin Step:", self.flat_profile_bin_step_spin)

        self.flat_pixel_min_spin = QDoubleSpinBox()
        self.flat_pixel_min_spin.setRange(0.01, 2.0)
        self.flat_pixel_min_spin.setSingleStep(0.01)
        flat_layout.addRow("Pixel Flat Min:", self.flat_pixel_min_spin)

        self.flat_pixel_max_spin = QDoubleSpinBox()
        self.flat_pixel_max_spin.setRange(0.1, 5.0)
        self.flat_pixel_max_spin.setSingleStep(0.01)
        flat_layout.addRow("Pixel Flat Max:", self.flat_pixel_max_spin)

        flat_group.setLayout(flat_layout)
        layout.addRow(flat_group)
>>>>>>> cef6f04 (	modified:   README.md)

        widget.setLayout(layout)
        scroll.setWidget(widget)
        return scroll

    def _update_overscan_params_enabled(self):
        """Enable/disable overscan parameters based on selected method."""
        method = self.overscan_method_combo.currentText()
        
        if method == 'mean_only':
            # No additional parameters needed
            self.overscan_smooth_window_spin.setEnabled(False)
            self.overscan_poly_type_combo.setEnabled(False)
            self.overscan_poly_order_spin.setEnabled(False)
        elif method == 'mean_savgol':
            # Only smooth window needed
            self.overscan_smooth_window_spin.setEnabled(True)
            self.overscan_poly_type_combo.setEnabled(False)
            self.overscan_poly_order_spin.setEnabled(False)
        elif method == 'mean_polynomial':
            # Poly type and order needed
            self.overscan_smooth_window_spin.setEnabled(False)
            self.overscan_poly_type_combo.setEnabled(True)
            self.overscan_poly_order_spin.setEnabled(True)

<<<<<<< HEAD
=======
    def _update_bg_params_enabled(self):
        """Enable/disable background parameters based on selected method."""
        method = self.bg_method_combo.currentText()

        # Method-specific params
        is_conv = (method == 'convolution')
        is_spline = (method == 'column_spline')
        is_cheby = (method == 'chebyshev')
        is_bspline = (method == 'bspline')
        is_gauss = (method == 'gaussian_smooth')

        self.bg_kernel_sigma_x_spin.setEnabled(is_conv)
        self.bg_kernel_sigma_y_spin.setEnabled(is_conv)
        self.bg_spline_smooth_factor_spin.setEnabled(is_spline)
        self.bg_spline_post_smooth_x_spin.setEnabled(is_spline)
        self.bg_poly_order_spin.setEnabled(is_cheby or is_bspline)
        self.bg_smooth_sigma_spin.setEnabled(is_gauss)
        self.bg_bspline_smooth_spin.setEnabled(is_bspline)

        # Common params always enabled: sigma_clip, clip_iterations, clip_mode, mask_margin

>>>>>>> cef6f04 (	modified:   README.md)
    def _create_wlcalib_tab(self) -> QWidget:
        """Create wavelength calibration settings tab."""
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        widget = QWidget()
        layout = QFormLayout()

        # Calibration lamp type and files
        calib_group = QGroupBox("Calibration Files")
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

        # Main calibration parameters
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
        trace_group = QGroupBox("Orders Tracing")
        trace_layout = QFormLayout()

        # Order Spacing Tolerance
        self.trace_spacing_tol_spin = QDoubleSpinBox()
        self.trace_spacing_tol_spin.setRange(0.05, 1.0)
        self.trace_spacing_tol_spin.setSingleStep(0.01)
        self.trace_spacing_tol_spin.setDecimals(2)
        self.trace_spacing_tol_spin.setToolTip(
            "Order spacing tolerance: maximum allowed deviation (fraction) from predicted spacing for order acceptance.\n"
            "Default 0.30 = ±30%. Orders with spacing deviating more than this from prediction are rejected.")
        trace_layout.addRow("Order Spacing Tolerance:", self.trace_spacing_tol_spin)

        self.trace_step_denominator_spin = QSpinBox()
        self.trace_step_denominator_spin.setRange(5, 2000)
        self.trace_step_denominator_spin.setToolTip(
            "Column step size for trace refinement (e.g., 20 means 1 sample every 20 columns).")
        trace_layout.addRow("Trace Step (columns):", self.trace_step_denominator_spin)

        self.trace_snr_threshold_spin = QDoubleSpinBox()
        self.trace_snr_threshold_spin.setRange(1.0, 100000.0)
        self.trace_snr_threshold_spin.setSingleStep(0.5)
        self.trace_snr_threshold_spin.setToolTip(
            "Order detection SNR threshold (multiples of sigma above background).\n"
            "Typical: 3.0–10.0. Lower to detect very faint orders.")
        trace_layout.addRow("Detection SNR Threshold:", self.trace_snr_threshold_spin)

<<<<<<< HEAD
        self.trace_seed_threshold_spin = QDoubleSpinBox()
        self.trace_seed_threshold_spin.setRange(0.01, 1.0)
        self.trace_seed_threshold_spin.setSingleStep(0.01)
        trace_layout.addRow("Seed Threshold (0-1):", self.trace_seed_threshold_spin)

        self.trace_prominence_scale_spin = QDoubleSpinBox()
        self.trace_prominence_scale_spin.setRange(0.01, 5.0)
        self.trace_prominence_scale_spin.setSingleStep(0.05)
        trace_layout.addRow("Prominence Scale:", self.trace_prominence_scale_spin)

        self.trace_search_half_scale_spin = QDoubleSpinBox()
        self.trace_search_half_scale_spin.setRange(0.1, 2.0)
        self.trace_search_half_scale_spin.setSingleStep(0.05)
        trace_layout.addRow("Search Half Scale:", self.trace_search_half_scale_spin)

        self.trace_step_denominator_spin = QSpinBox()
        self.trace_step_denominator_spin.setRange(40, 2000)
        trace_layout.addRow("Trace Step Denominator:", self.trace_step_denominator_spin)

        self.trace_fill_missing_check = QCheckBox("Fill missing weak orders")
        trace_layout.addRow(self.trace_fill_missing_check)

        self.trace_gap_fill_factor_spin = QDoubleSpinBox()
        self.trace_gap_fill_factor_spin.setRange(1.0, 4.0)
        self.trace_gap_fill_factor_spin.setSingleStep(0.1)
        trace_layout.addRow("Gap Fill Factor:", self.trace_gap_fill_factor_spin)

        self.trace_fit_method_combo = QComboBox()
        self.trace_fit_method_combo.addItems(['polynomial', 'chebyshev', 'bspline'])
        trace_layout.addRow("Trace Fit Method (Center+Edges):", self.trace_fit_method_combo)

        self.trace_bspline_smooth_spin = QDoubleSpinBox()
        self.trace_bspline_smooth_spin.setRange(0.0, 5.0)
        self.trace_bspline_smooth_spin.setSingleStep(0.05)
        trace_layout.addRow("B-spline Smooth:", self.trace_bspline_smooth_spin)

        self.trace_edge_degree_spin = QSpinBox()
        self.trace_edge_degree_spin.setRange(1, 8)
        trace_layout.addRow("Edge Fit Degree:", self.trace_edge_degree_spin)

        self.trace_aperture_root_fraction_spin = QDoubleSpinBox()
        self.trace_aperture_root_fraction_spin.setRange(0.0, 0.5)
        self.trace_aperture_root_fraction_spin.setSingleStep(0.005)
        self.trace_aperture_root_fraction_spin.setDecimals(3)
        trace_layout.addRow("Aperture Root Fraction:", self.trace_aperture_root_fraction_spin)

        self.trace_aperture_noise_floor_sigma_spin = QDoubleSpinBox()
        self.trace_aperture_noise_floor_sigma_spin.setRange(0.0, 20.0)
        self.trace_aperture_noise_floor_sigma_spin.setSingleStep(0.1)
        self.trace_aperture_noise_floor_sigma_spin.setDecimals(2)
        trace_layout.addRow("Aperture Noise Floor Sigma:", self.trace_aperture_noise_floor_sigma_spin)

        self.trace_filling_spin = QDoubleSpinBox()
        self.trace_filling_spin.setRange(0.0, 1.0)
        self.trace_filling_spin.setSingleStep(0.05)
        trace_layout.addRow("Filling Factor Threshold:", self.trace_filling_spin)
=======
        self.trace_min_coverage_spin = QDoubleSpinBox()
        self.trace_min_coverage_spin.setRange(0.05, 0.80)
        self.trace_min_coverage_spin.setSingleStep(0.05)
        self.trace_min_coverage_spin.setDecimals(2)
        self.trace_min_coverage_spin.setToolTip(
            "Minimum fraction of detector width an order must span.\n"
            "Lower to accept faint partial orders; raise to reject short false traces.")
        trace_layout.addRow("Min Trace Coverage:", self.trace_min_coverage_spin)
>>>>>>> cef6f04 (	modified:   README.md)

        self.trace_degree_spin = QSpinBox()
        self.trace_degree_spin.setRange(1, 10)
        self.trace_degree_spin.setToolTip(
            "Polynomial degree for Chebyshev fit of the order trace (applies to both center and edges).")
        trace_layout.addRow("Trace Degree:", self.trace_degree_spin)

        self.trace_width_cheb_degree_spin = QSpinBox()
        self.trace_width_cheb_degree_spin.setRange(1, 10)
        self.trace_width_cheb_degree_spin.setToolTip("Chebyshev polynomial degree for fitting the order width.")
        trace_layout.addRow("Width Cheb Degree:", self.trace_width_cheb_degree_spin)

        self.trace_boundary_frac_spin = QDoubleSpinBox()
        self.trace_boundary_frac_spin.setRange(0.001, 0.2)
        self.trace_boundary_frac_spin.setSingleStep(0.001)
        self.trace_boundary_frac_spin.setDecimals(3)
        self.trace_boundary_frac_spin.setToolTip(
            "Initial boundary threshold: fraction of peak for preliminary order edge (default 0.02).\n"
            "Lower = wider initial region for Moffat fit.")
        trace_layout.addRow("Boundary Fraction (Initial):", self.trace_boundary_frac_spin)

        self.trace_fwhm_scale_spin = QDoubleSpinBox()
        self.trace_fwhm_scale_spin.setRange(0.5, 3.0)
        self.trace_fwhm_scale_spin.setSingleStep(0.05)
        self.trace_fwhm_scale_spin.setDecimals(2)
        self.trace_fwhm_scale_spin.setToolTip(
            "Final boundary = center ± FWHM × scale (default 1.5).\n"
            "Controls the width of the extracted order region.")
        trace_layout.addRow("FWHM Scale (Final):", self.trace_fwhm_scale_spin)

        # --- Gap-fill parameters (only relevant when checkbox is checked) ---

        self.trace_gap_fill_factor_spin = QDoubleSpinBox()
        self.trace_gap_fill_factor_spin.setRange(1.0, 4.0)
        self.trace_gap_fill_factor_spin.setSingleStep(0.1)
        trace_layout.addRow("Gap Fill Factor:", self.trace_gap_fill_factor_spin)

        self.trace_gap_fill_factor_interp_spin = QDoubleSpinBox()
        self.trace_gap_fill_factor_interp_spin.setRange(1.0, 4.0)
        self.trace_gap_fill_factor_interp_spin.setSingleStep(0.05)
        self.trace_gap_fill_factor_interp_spin.setToolTip("Gap factor to trigger missing-order insertion during interpolation.")
        trace_layout.addRow("Gap Fill Factor (Interp):", self.trace_gap_fill_factor_interp_spin)

        self.trace_gap_fill_snr_spin = QDoubleSpinBox()
        self.trace_gap_fill_snr_spin.setRange(0.1, 5.0)
        self.trace_gap_fill_snr_spin.setSingleStep(0.1)
        self.trace_gap_fill_snr_spin.setToolTip(
            "Gap-fill peak acceptance threshold (multiples of robust noise).\n"
            "Lower values recover fainter orders; raise to avoid false insertions.")
        trace_layout.addRow("Gap Fill SNR:", self.trace_gap_fill_snr_spin)

        self.trace_n_extend_below_spin = QSpinBox()
        self.trace_n_extend_below_spin.setRange(0, 10)
        trace_layout.addRow("Extend Below (Step2):", self.trace_n_extend_below_spin)

        self.trace_n_extend_above_spin = QSpinBox()
        self.trace_n_extend_above_spin.setRange(0, 10)
        trace_layout.addRow("Extend Above (Step2):", self.trace_n_extend_above_spin)

        self.trace_n_mask_below_spin = QSpinBox()
        self.trace_n_mask_below_spin.setRange(0, 10)
        trace_layout.addRow("Mask Below (Step3):", self.trace_n_mask_below_spin)

        self.trace_n_mask_above_spin = QSpinBox()
        self.trace_n_mask_above_spin.setRange(0, 10)
        trace_layout.addRow("Mask Above (Step3):", self.trace_n_mask_above_spin)

        trace_group.setLayout(trace_layout)

        # Add Orders tracing参数组到主layout
        layout.addRow(trace_group)

        extract_group = QGroupBox("Spectrum Extraction")
        extract_layout = QFormLayout()
        self.extract_method_combo = QComboBox()
        self.extract_method_combo.addItems(['sum', 'optimal'])
        extract_layout.addRow("Extraction Method:", self.extract_method_combo)

        extract_group.setLayout(extract_layout)
        layout.addRow(extract_group)

        widget.setLayout(layout)
        scroll.setWidget(widget)
        # Initialize UI state
        self._update_trace_param_interactions()
        return scroll

    def _update_trace_param_interactions(self):
        """Enable/disable parameters (Placeholder, logic simplified)."""
        pass

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

        # Figure format
        self.fig_format_combo = QComboBox()
        self.fig_format_combo.addItems(['png', 'pdf', 'jpg'])
        layout.addRow("Figure Format:", self.fig_format_combo)

        # Output suffix
        self.oned_suffix_edit = QLineEdit()
        layout.addRow("Output Suffix:", self.oned_suffix_edit)

        # Auto process
        self.auto_process_check = QCheckBox("Auto process all science images without confirmation")
        layout.addRow(self.auto_process_check)

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
        # Data & Instruments tab
        self.rawdata_path_edit.setText(self.config.get('data', 'rawdata_path', './rawdata'))
        self.statime_key_edit.setText(self.config.get('data', 'statime_key', 'DATE-OBS'))
        self.exptime_key_edit.setText(self.config.get('data', 'exptime_key', 'EXPTIME'))
        self.direction_combo.setCurrentText(self.config.get('data', 'direction', 'xr-'))
        self.output_path_edit.setText(self.config.get('reduce', 'output_path', 'output'))
        self.telescope_name_edit.setText(self.config.get('telescope', 'name', 'xinglong216hrs'))
        self.instrument_edit.setText(self.config.get('telescope', 'instrument', 'hrs'))
        self.detector_gain_spin.setValue(self.config.get_float('telescope.detector', 'gain', 1.0))
        self.detector_readnoise_spin.setValue(self.config.get_float('telescope.detector', 'readnoise', 5.0))

        # Data Reduction tab
        overscan_val = self.config.get_int('data', 'overscan_start_column', -1)
        self.overscan_start_column_spin.setValue(overscan_val)
        self.overscan_method_combo.setCurrentText(self.config.get('data', 'overscan_method', 'mean_only'))
        self.overscan_smooth_window_spin.setValue(self.config.get_int('data', 'overscan_smooth_window', -1))
        self.overscan_poly_type_combo.setCurrentText(self.config.get('data', 'overscan_poly_type', 'legendre'))
        self.overscan_poly_order_spin.setValue(self.config.get_int('data', 'overscan_poly_order', 3))
<<<<<<< HEAD
=======
        self.detector_split_row_spin.setValue(self.config.get_int('data', 'detector_split_row', 2068))
>>>>>>> cef6f04 (	modified:   README.md)
        # Update enabled state of overscan params based on selected method
        self._update_overscan_params_enabled()
        self.bias_combine_method_combo.setCurrentText(self.config.get('reduce.bias', 'combine_method', 'median'))
        self.bias_combine_sigma_spin.setValue(self.config.get_float('reduce.bias', 'combine_sigma', 3.0))
        self.flat_combine_method_combo.setCurrentText(self.config.get('reduce.flat', 'combine_method', 'median'))
        self.flat_q_threshold_spin.setValue(self.config.get_float('reduce.flat', 'q_threshold', 0.5))
        self.flat_mosaic_maxcount_spin.setValue(self.config.get_float('reduce.flat', 'mosaic_maxcount', 65535.0))
<<<<<<< HEAD
        self.flat_blaze_smooth_method_combo.setCurrentText(self.config.get('reduce.flat', 'blaze_smooth_method', 'savgol'))
        self.flat_blaze_smooth_window_spin.setValue(self.config.get_int('reduce.flat', 'blaze_smooth_window', 21))
        self.flat_blaze_bspline_smooth_spin.setValue(self.config.get_float('reduce.flat', 'blaze_bspline_smooth', 0.5))
        self.flat_width_smooth_window_spin.setValue(self.config.get_int('reduce.flat', 'width_smooth_window', 41))
        self.flat_profile_bin_step_spin.setValue(self.config.get_float('reduce.flat', 'profile_bin_step', 0.01))
        self.flat_pixel_min_spin.setValue(self.config.get_float('reduce.flat', 'pixel_flat_min', 0.5))
        self.flat_pixel_max_spin.setValue(self.config.get_float('reduce.flat', 'pixel_flat_max', 1.5))
        bg_method = self.config.get('reduce.background', 'method', 'chebyshev')
        if bg_method == '2d_poly':
            bg_method = 'chebyshev'
        self.bg_method_combo.setCurrentText(bg_method)
=======
        self.flat_blaze_knot_spacing_spin.setValue(self.config.get_int('reduce.flat', 'blaze_knot_spacing', 500))
        self.flat_blaze_edge_nknots_spin.setValue(self.config.get_int('reduce.flat', 'blaze_edge_nknots', 6))
        self.flat_width_smooth_window_spin.setValue(self.config.get_int('reduce.flat', 'width_smooth_window', 41))
        self.flat_n_profile_segments_spin.setValue(self.config.get_int('reduce.flat', 'n_profile_segments', 100))
        self.flat_profile_smooth_sigma_spin.setValue(self.config.get_float('reduce.flat', 'profile_smooth_sigma', 6.0))
        self.flat_profile_bin_step_spin.setValue(self.config.get_float('reduce.flat', 'profile_bin_step', 0.01))
        self.flat_pixel_min_spin.setValue(self.config.get_float('reduce.flat', 'pixel_flat_min', 0.5))
        self.flat_pixel_max_spin.setValue(self.config.get_float('reduce.flat', 'pixel_flat_max', 1.5))
        bg_method = self.config.get('reduce.background', 'method', 'convolution')
        if bg_method == '2d_poly':
            bg_method = 'chebyshev'
        elif bg_method == 'smooth':
            bg_method = 'gaussian_smooth'
        self.bg_method_combo.setCurrentText(bg_method)
        self.bg_kernel_sigma_x_spin.setValue(self.config.get_float('reduce.background', 'kernel_sigma_x', 13.0))
        self.bg_kernel_sigma_y_spin.setValue(self.config.get_float('reduce.background', 'kernel_sigma_y', 13.0))
        self.bg_spline_smooth_factor_spin.setValue(self.config.get_float('reduce.background', 'spline_smooth_factor', 1.0))
        self.bg_spline_post_smooth_x_spin.setValue(self.config.get_float('reduce.background', 'spline_post_smooth_x', 5.0))
>>>>>>> cef6f04 (	modified:   README.md)
        self.bg_poly_order_spin.setValue(self.config.get_int('reduce.background', 'poly_order', 3))
        self.bg_smooth_sigma_spin.setValue(self.config.get_float('reduce.background', 'smooth_sigma', 20.0))
        self.bg_sigma_clip_spin.setValue(self.config.get_float('reduce.background', 'sigma_clip', 3.0))
        self.bg_sigma_clip_iter_spin.setValue(self.config.get_int('reduce.background', 'sigma_clip_maxiters', 4))
<<<<<<< HEAD
        self.bg_mask_margin_spin.setValue(self.config.get_int('reduce.background', 'mask_margin_pixels', 3))
        self.bg_bspline_smooth_spin.setValue(self.config.get_float('reduce.background', 'bspline_smooth', 1.0))
=======
        clip_mode = self.config.get('reduce.background', 'sigma_clip_mode', 'upper')
        self.bg_clip_mode_combo.setCurrentText(clip_mode)
        self.bg_mask_margin_spin.setValue(self.config.get_int('reduce.background', 'mask_margin_pixels', 1))
        self.bg_bspline_smooth_spin.setValue(self.config.get_float('reduce.background', 'bspline_smooth', 1.0))
        self._update_bg_params_enabled()
>>>>>>> cef6f04 (	modified:   README.md)
        self.cosmic_enabled_check.setChecked(self.config.get_bool('reduce', 'cosmic_enabled', True))
        self.cosmic_sigma_spin.setValue(self.config.get_float('reduce', 'cosmic_sigma', 5.0))
        self.cosmic_window_spin.setValue(self.config.get_int('reduce', 'cosmic_window', 5))

        # Wavelength Calibration tab
        self.linelist_combo.setCurrentText(self.config.get('telescope.linelist', 'linelist_type', 'ThAr'))
        self.linelist_path_edit.setText(self.config.get('telescope.linelist', 'linelist_path', 'calib_data/linelists/'))
        self.linelist_file_edit.setText(self.config.get('telescope.linelist', 'linelist_file', 'thar-noao.dat'))
        self.use_precomputed_calib_check.setChecked(self.config.get_bool('telescope.linelist', 'use_precomputed_calibration', False))
        self.calibration_path_edit.setText(self.config.get('telescope.linelist', 'calibration_path', 'calib_data/telescopes/xinglong216hrs/'))
        self.calibration_file_edit.setText(self.config.get('telescope.linelist', 'calibration_file', 'wlcalib_20211123011_A.fits'))
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

<<<<<<< HEAD
        # Trace & Extraction tab
        self.trace_scan_step_spin.setValue(self.config.get_int('reduce.trace', 'scan_step', 10))
        self.trace_minimum_spin.setValue(self.config.get_float('reduce.trace', 'minimum', 50.0))
        self.trace_separation_spin.setValue(self.config.get_float('reduce.trace', 'separation', 30.0))
        self.trace_seed_threshold_spin.setValue(self.config.get_float('reduce.trace', 'seed_threshold', 0.30))
        self.trace_prominence_scale_spin.setValue(self.config.get_float('reduce.trace', 'prominence_scale', 0.50))
        self.trace_search_half_scale_spin.setValue(self.config.get_float('reduce.trace', 'search_half_scale', 0.45))
        self.trace_step_denominator_spin.setValue(self.config.get_int('reduce.trace', 'step_denominator', 220))
        self.trace_fill_missing_check.setChecked(self.config.get_bool('reduce.trace', 'fill_missing_orders', True))
        self.trace_gap_fill_factor_spin.setValue(self.config.get_float('reduce.trace', 'gap_fill_factor', 1.6))
        self.trace_fit_method_combo.setCurrentText(self.config.get('reduce.trace', 'fit_method', 'polynomial'))
        self.trace_bspline_smooth_spin.setValue(self.config.get_float('reduce.trace', 'bspline_smooth', 0.2))
        self.trace_edge_degree_spin.setValue(self.config.get_int('reduce.trace', 'edge_degree', 3))
        self.trace_aperture_root_fraction_spin.setValue(self.config.get_float('reduce.trace', 'aperture_root_fraction', 0.03))
        self.trace_aperture_noise_floor_sigma_spin.setValue(self.config.get_float('reduce.trace', 'aperture_noise_floor_sigma', 3.0))
        self.trace_filling_spin.setValue(self.config.get_float('reduce.trace', 'filling', 0.3))
        self.trace_degree_spin.setValue(self.config.get_int('reduce.trace', 'degree', 3))
        self.extract_method_combo.setCurrentText(self.config.get('reduce.extract', 'method', 'sum'))
        self.extract_lower_limit_spin.setValue(self.config.get_float('reduce.extract', 'lower_limit', -5.0))
        self.extract_upper_limit_spin.setValue(self.config.get_float('reduce.extract', 'upper_limit', 5.0))
=======
        # Orders Tracing tab
        self.trace_spacing_tol_spin.setValue(self.config.get_float('reduce.trace', 'spacing_tol', 0.3))
        self.trace_step_denominator_spin.setValue(self.config.get_int('reduce.trace', 'step_denominator', 20))
        self.trace_snr_threshold_spin.setValue(self.config.get_float('reduce.trace', 'snr_threshold', 5.0))
        self.trace_gap_fill_factor_spin.setValue(self.config.get_float('reduce.trace', 'gap_fill_factor', 1.6))
        self.trace_gap_fill_factor_interp_spin.setValue(self.config.get_float('reduce.trace', 'gap_fill_factor_interp', 1.35))
        self.trace_gap_fill_snr_spin.setValue(self.config.get_float('reduce.trace', 'gap_fill_snr', 2.5))
        self.trace_min_coverage_spin.setValue(self.config.get_float('reduce.trace', 'min_trace_coverage', 0.20))
        self.trace_degree_spin.setValue(self.config.get_int('reduce.trace', 'degree', 4))
        self.trace_width_cheb_degree_spin.setValue(self.config.get_int('reduce.trace', 'width_cheb_degree', 3))
        self.trace_n_extend_below_spin.setValue(self.config.get_int('reduce.trace', 'n_extend_below', 0))
        self.trace_n_extend_above_spin.setValue(self.config.get_int('reduce.trace', 'n_extend_above', 0))
        self.trace_n_mask_below_spin.setValue(self.config.get_int('reduce.trace', 'n_mask_below', 2))
        self.trace_n_mask_above_spin.setValue(self.config.get_int('reduce.trace', 'n_mask_above', 1))
        self.trace_boundary_frac_spin.setValue(self.config.get_float('reduce.trace', 'boundary_frac', 0.02))
        self.trace_fwhm_scale_spin.setValue(self.config.get_float('reduce.trace', 'fwhm_scale', 1.5))
        self.extract_method_combo.setCurrentText(self.config.get('reduce.extract', 'method', 'optimal'))
>>>>>>> cef6f04 (	modified:   README.md)

        # Processing & Output tab
        self.mode_combo.setCurrentText(self.config.get('reduce', 'mode', 'normal'))
        self.fig_format_combo.setCurrentText(self.config.get('reduce', 'fig_format', 'png'))
        self.oned_suffix_edit.setText(self.config.get('reduce', 'oned_suffix', '_ods'))
        self.auto_process_check.setChecked(self.config.get_bool('reduce', 'auto_process_all_images', False))
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
            # Data & Instruments tab
            self.config.set('data', 'rawdata_path', self.rawdata_path_edit.text())
            self.config.set('data', 'statime_key', self.statime_key_edit.text())
            self.config.set('data', 'exptime_key', self.exptime_key_edit.text())
            self.config.set('data', 'direction', self.direction_combo.currentText())
            self.config.set('reduce', 'output_path', self.output_path_edit.text())
            self.config.set('telescope', 'name', self.telescope_name_edit.text())
            self.config.set('telescope', 'instrument', self.instrument_edit.text())
            self.config.set('telescope.detector', 'gain', str(self.detector_gain_spin.value()))
            self.config.set('telescope.detector', 'readnoise', str(self.detector_readnoise_spin.value()))

            # Data Reduction tab
            self.config.set('data', 'overscan_start_column', str(self.overscan_start_column_spin.value()))
            self.config.set('data', 'overscan_method', self.overscan_method_combo.currentText())
            self.config.set('data', 'overscan_smooth_window', str(self.overscan_smooth_window_spin.value()))
            self.config.set('data', 'overscan_poly_type', self.overscan_poly_type_combo.currentText())
            self.config.set('data', 'overscan_poly_order', str(self.overscan_poly_order_spin.value()))
<<<<<<< HEAD
=======
            self.config.set('data', 'detector_split_row', str(self.detector_split_row_spin.value()))
>>>>>>> cef6f04 (	modified:   README.md)
            self.config.set('reduce.bias', 'combine_method', self.bias_combine_method_combo.currentText())
            self.config.set('reduce.bias', 'combine_sigma', str(self.bias_combine_sigma_spin.value()))
            self.config.set('reduce.flat', 'combine_method', self.flat_combine_method_combo.currentText())
            self.config.set('reduce.flat', 'q_threshold', str(self.flat_q_threshold_spin.value()))
            self.config.set('reduce.flat', 'mosaic_maxcount', str(self.flat_mosaic_maxcount_spin.value()))
<<<<<<< HEAD
            self.config.set('reduce.flat', 'blaze_smooth_method', self.flat_blaze_smooth_method_combo.currentText())
            self.config.set('reduce.flat', 'blaze_smooth_window', str(self.flat_blaze_smooth_window_spin.value()))
            self.config.set('reduce.flat', 'blaze_bspline_smooth', str(self.flat_blaze_bspline_smooth_spin.value()))
            self.config.set('reduce.flat', 'width_smooth_window', str(self.flat_width_smooth_window_spin.value()))
=======
            self.config.set('reduce.flat', 'blaze_knot_spacing', str(self.flat_blaze_knot_spacing_spin.value()))
            self.config.set('reduce.flat', 'blaze_edge_nknots', str(self.flat_blaze_edge_nknots_spin.value()))
            self.config.set('reduce.flat', 'width_smooth_window', str(self.flat_width_smooth_window_spin.value()))
            self.config.set('reduce.flat', 'n_profile_segments', str(self.flat_n_profile_segments_spin.value()))
            self.config.set('reduce.flat', 'profile_smooth_sigma', str(self.flat_profile_smooth_sigma_spin.value()))
>>>>>>> cef6f04 (	modified:   README.md)
            self.config.set('reduce.flat', 'profile_bin_step', str(self.flat_profile_bin_step_spin.value()))
            self.config.set('reduce.flat', 'pixel_flat_min', str(self.flat_pixel_min_spin.value()))
            self.config.set('reduce.flat', 'pixel_flat_max', str(self.flat_pixel_max_spin.value()))
            self.config.set('reduce.background', 'method', self.bg_method_combo.currentText())
<<<<<<< HEAD
=======
            self.config.set('reduce.background', 'kernel_sigma_x', str(self.bg_kernel_sigma_x_spin.value()))
            self.config.set('reduce.background', 'kernel_sigma_y', str(self.bg_kernel_sigma_y_spin.value()))
            self.config.set('reduce.background', 'spline_smooth_factor', str(self.bg_spline_smooth_factor_spin.value()))
            self.config.set('reduce.background', 'spline_post_smooth_x', str(self.bg_spline_post_smooth_x_spin.value()))
>>>>>>> cef6f04 (	modified:   README.md)
            self.config.set('reduce.background', 'poly_order', str(self.bg_poly_order_spin.value()))
            self.config.set('reduce.background', 'smooth_sigma', str(self.bg_smooth_sigma_spin.value()))
            self.config.set('reduce.background', 'sigma_clip', str(self.bg_sigma_clip_spin.value()))
            self.config.set('reduce.background', 'sigma_clip_maxiters', str(self.bg_sigma_clip_iter_spin.value()))
<<<<<<< HEAD
=======
            self.config.set('reduce.background', 'sigma_clip_mode', self.bg_clip_mode_combo.currentText())
>>>>>>> cef6f04 (	modified:   README.md)
            self.config.set('reduce.background', 'mask_margin_pixels', str(self.bg_mask_margin_spin.value()))
            self.config.set('reduce.background', 'bspline_smooth', str(self.bg_bspline_smooth_spin.value()))
            self.config.set('reduce', 'cosmic_enabled', 'yes' if self.cosmic_enabled_check.isChecked() else 'no')
            self.config.set('reduce', 'cosmic_sigma', str(self.cosmic_sigma_spin.value()))
            self.config.set('reduce', 'cosmic_window', str(self.cosmic_window_spin.value()))
<<<<<<< HEAD
            self.config.set('reduce.bias', 'combine_sigma', str(self.bias_combine_sigma_spin.value()))
            self.config.set('reduce.flat', 'combine_method', self.flat_combine_method_combo.currentText())
            self.config.set('reduce.flat', 'q_threshold', str(self.flat_q_threshold_spin.value()))
            self.config.set('reduce.flat', 'mosaic_maxcount', str(self.flat_mosaic_maxcount_spin.value()))
            self.config.set('reduce.flat', 'blaze_smooth_method', self.flat_blaze_smooth_method_combo.currentText())
            self.config.set('reduce.flat', 'blaze_smooth_window', str(self.flat_blaze_smooth_window_spin.value()))
            self.config.set('reduce.flat', 'blaze_bspline_smooth', str(self.flat_blaze_bspline_smooth_spin.value()))
            self.config.set('reduce.flat', 'width_smooth_window', str(self.flat_width_smooth_window_spin.value()))
            self.config.set('reduce.flat', 'profile_bin_step', str(self.flat_profile_bin_step_spin.value()))
            self.config.set('reduce.flat', 'pixel_flat_min', str(self.flat_pixel_min_spin.value()))
            self.config.set('reduce.flat', 'pixel_flat_max', str(self.flat_pixel_max_spin.value()))
            self.config.set('reduce.background', 'method', self.bg_method_combo.currentText())
            self.config.set('reduce.background', 'poly_order', str(self.bg_poly_order_spin.value()))
            self.config.set('reduce.background', 'smooth_sigma', str(self.bg_smooth_sigma_spin.value()))
            self.config.set('reduce.background', 'sigma_clip', str(self.bg_sigma_clip_spin.value()))
            self.config.set('reduce.background', 'sigma_clip_maxiters', str(self.bg_sigma_clip_iter_spin.value()))
            self.config.set('reduce.background', 'mask_margin_pixels', str(self.bg_mask_margin_spin.value()))
            self.config.set('reduce.background', 'bspline_smooth', str(self.bg_bspline_smooth_spin.value()))
=======
>>>>>>> cef6f04 (	modified:   README.md)

            # Wavelength Calibration tab
            self.config.set('telescope.linelist', 'linelist_type', self.linelist_combo.currentText())
            self.config.set('telescope.linelist', 'linelist_path', self.linelist_path_edit.text())
            self.config.set('telescope.linelist', 'linelist_file', self.linelist_file_edit.text())
            self.config.set('telescope.linelist', 'use_precomputed_calibration',
                           'yes' if self.use_precomputed_calib_check.isChecked() else 'no')
            self.config.set('telescope.linelist', 'calibration_path', self.calibration_path_edit.text())
            self.config.set('telescope.linelist', 'calibration_file', self.calibration_file_edit.text())
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

<<<<<<< HEAD
            # Trace & Extraction tab
            self.config.set('reduce.trace', 'scan_step', str(self.trace_scan_step_spin.value()))
            self.config.set('reduce.trace', 'minimum', str(self.trace_minimum_spin.value()))
            self.config.set('reduce.trace', 'separation', str(self.trace_separation_spin.value()))
            self.config.set('reduce.trace', 'seed_threshold', str(self.trace_seed_threshold_spin.value()))
            self.config.set('reduce.trace', 'prominence_scale', str(self.trace_prominence_scale_spin.value()))
            self.config.set('reduce.trace', 'search_half_scale', str(self.trace_search_half_scale_spin.value()))
            self.config.set('reduce.trace', 'step_denominator', str(self.trace_step_denominator_spin.value()))
            self.config.set('reduce.trace', 'fill_missing_orders',
                           'yes' if self.trace_fill_missing_check.isChecked() else 'no')
            self.config.set('reduce.trace', 'gap_fill_factor', str(self.trace_gap_fill_factor_spin.value()))
            self.config.set('reduce.trace', 'fit_method', self.trace_fit_method_combo.currentText())
            self.config.set('reduce.trace', 'bspline_smooth', str(self.trace_bspline_smooth_spin.value()))
            self.config.set('reduce.trace', 'edge_degree', str(self.trace_edge_degree_spin.value()))
            self.config.set('reduce.trace', 'aperture_root_fraction', str(self.trace_aperture_root_fraction_spin.value()))
            self.config.set('reduce.trace', 'aperture_noise_floor_sigma', str(self.trace_aperture_noise_floor_sigma_spin.value()))
            self.config.set('reduce.trace', 'filling', str(self.trace_filling_spin.value()))
=======
            # Orders Tracing tab
            self.config.set('reduce.trace', 'step_denominator', str(self.trace_step_denominator_spin.value()))
            self.config.set('reduce.trace', 'snr_threshold', str(self.trace_snr_threshold_spin.value()))
            self.config.set('reduce.trace', 'gap_fill_factor', str(self.trace_gap_fill_factor_spin.value()))
            self.config.set('reduce.trace', 'gap_fill_factor_interp', str(self.trace_gap_fill_factor_interp_spin.value()))
            self.config.set('reduce.trace', 'gap_fill_snr', str(self.trace_gap_fill_snr_spin.value()))
            self.config.set('reduce.trace', 'min_trace_coverage', str(self.trace_min_coverage_spin.value()))
>>>>>>> cef6f04 (	modified:   README.md)
            self.config.set('reduce.trace', 'degree', str(self.trace_degree_spin.value()))
            self.config.set('reduce.trace', 'width_cheb_degree', str(self.trace_width_cheb_degree_spin.value()))
            self.config.set('reduce.trace', 'n_extend_below', str(self.trace_n_extend_below_spin.value()))
            self.config.set('reduce.trace', 'n_extend_above', str(self.trace_n_extend_above_spin.value()))
            self.config.set('reduce.trace', 'n_mask_below', str(self.trace_n_mask_below_spin.value()))
            self.config.set('reduce.trace', 'n_mask_above', str(self.trace_n_mask_above_spin.value()))
            self.config.set('reduce.trace', 'boundary_frac', str(self.trace_boundary_frac_spin.value()))
            self.config.set('reduce.trace', 'fwhm_scale', str(self.trace_fwhm_scale_spin.value()))
            self.config.set('reduce.trace', 'spacing_tol', str(self.trace_spacing_tol_spin.value()))
            self.config.set('reduce.extract', 'method', self.extract_method_combo.currentText())

            # Processing & Output tab
            self.config.set('reduce', 'mode', self.mode_combo.currentText())
            self.config.set('reduce', 'fig_format', self.fig_format_combo.currentText())
            self.config.set('reduce', 'oned_suffix', self.oned_suffix_edit.text())
            self.config.set('reduce', 'auto_process_all_images', 'yes' if self.auto_process_check.isChecked() else 'no')
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