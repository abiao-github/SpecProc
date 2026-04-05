"""
Settings dialog for SpecProc configuration.

Provides a GUI for editing configuration parameters.
"""

from PyQt5.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit,
    QSpinBox, QDoubleSpinBox, QComboBox, QCheckBox,
    QPushButton, QTabWidget, QWidget, QGroupBox, QFormLayout
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
        self.resize(600, 500)

        self.init_ui()
        self.load_current_values()

    def init_ui(self):
        """Initialize user interface."""
        layout = QVBoxLayout()

        # Single settings tab containing all options
        all_tab = self._create_all_settings_tab()
        layout.addWidget(all_tab)

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

    def _create_all_settings_tab(self) -> QWidget:
        """Create combined settings tab containing all sections."""
        widget = QWidget()
        layout = QVBoxLayout()

        data_section = self._create_data_tab()
        reduce_section = self._create_reduce_tab()
        overscan_section = self._create_overscan_tab()

        layout.addWidget(data_section)
        layout.addWidget(reduce_section)
        layout.addWidget(overscan_section)

        widget.setLayout(layout)
        return widget

    def _create_data_tab(self) -> QWidget:
        """Create data settings section."""
        widget = QWidget()
        layout = QFormLayout()

        self.statime_key_edit = QLineEdit()
        layout.addRow("Start Time Key:", self.statime_key_edit)

        self.exptime_key_edit = QLineEdit()
        layout.addRow("Exposure Time Key:", self.exptime_key_edit)

        # Direction
        self.direction_combo = QComboBox()
        self.direction_combo.addItems(['xr-', 'xl-', 'yr-', 'yl-'])
        layout.addRow("Dispersion Direction:", self.direction_combo)

        widget.setLayout(layout)
        return widget

    def _create_reduce_tab(self) -> QWidget:
        """Create processing settings tab."""
        widget = QWidget()
        layout = QFormLayout()

        # Paths
        self.output_path_edit = QLineEdit()
        layout.addRow("Output Path:", self.output_path_edit)

        # Processing options
        self.mode_combo = QComboBox()
        self.mode_combo.addItems(['normal', 'debug'])
        layout.addRow("Processing Mode:", self.mode_combo)

        self.fig_format_combo = QComboBox()
        self.fig_format_combo.addItems(['png', 'pdf', 'jpg'])
        layout.addRow("Figure Format:", self.fig_format_combo)

        self.oned_suffix_edit = QLineEdit()
        layout.addRow("Output Suffix:", self.oned_suffix_edit)

        # Auto process all images
        self.auto_process_check = QCheckBox("Auto process all science images")
        layout.addRow(self.auto_process_check)

        widget.setLayout(layout)
        return widget

    def _create_overscan_tab(self) -> QWidget:
        """Create overscan settings tab."""
        widget = QWidget()
        layout = QVBoxLayout()

        # Overscan method
        method_layout = QHBoxLayout()
        method_layout.addWidget(QLabel("Overscan Method:"))
        self.overscan_method_combo = QComboBox()
        self.overscan_method_combo.addItems(['median', 'polynomial'])
        method_layout.addWidget(self.overscan_method_combo)
        method_layout.addStretch()
        layout.addLayout(method_layout)

        # Simplified overscan configuration
        overscan_group = QGroupBox("Overscan Region (Simplified)")
        overscan_layout = QFormLayout()

        self.overscan_start_column = QSpinBox()
        self.overscan_start_column.setRange(1, 10000)
        self.overscan_start_column.setValue(4097)
        overscan_layout.addRow("Overscan Start Column (1-based):", self.overscan_start_column)

        overscan_group.setLayout(overscan_layout)
        layout.addWidget(overscan_group)

        layout.addStretch()
        widget.setLayout(layout)
        return widget

    def load_current_values(self):
        """Load current configuration values into UI."""
        # Data tab
        self.statime_key_edit.setText(self.config.get('data', 'statime_key', 'DATE-OBS'))
        self.exptime_key_edit.setText(self.config.get('data', 'exptime_key', 'EXPTIME'))
        self.direction_combo.setCurrentText(self.config.get('data', 'direction', 'xr-'))

        # Reduce tab
        self.output_path_edit.setText(self.config.get('reduce', 'output_path', './output'))
        self.mode_combo.setCurrentText(self.config.get('reduce', 'mode', 'normal'))
        self.fig_format_combo.setCurrentText(self.config.get('reduce', 'fig_format', 'png'))
        self.oned_suffix_edit.setText(self.config.get('reduce', 'oned_suffix', '_ods'))
        self.auto_process_check.setChecked(self.config.get_bool('reduce', 'auto_process_all_images', False))

        # Overscan tab
        self.overscan_method_combo.setCurrentText(self.config.get('data', 'overscan_method', 'median'))

        # Simplified overscan
        self.overscan_start_column.setValue(self.config.get_int('data', 'overscan_start_column', 4097))


    def _save_settings(self):
        """Save settings to configuration."""
        try:
            # Data tab
            self.config.set('data', 'statime_key', self.statime_key_edit.text())
            self.config.set('data', 'exptime_key', self.exptime_key_edit.text())
            self.config.set('data', 'direction', self.direction_combo.currentText())

            # Reduce tab
            self.config.set('reduce', 'output_path', self.output_path_edit.text())
            self.config.set('reduce', 'mode', self.mode_combo.currentText())
            self.config.set('reduce', 'fig_format', self.fig_format_combo.currentText())
            self.config.set('reduce', 'oned_suffix', self.oned_suffix_edit.text())
            self.config.set('reduce', 'auto_process_all_images', 'yes' if self.auto_process_check.isChecked() else 'no')

            # Overscan tab
            self.config.set('data', 'overscan_method', self.overscan_method_combo.currentText())

            # Simplified overscan
            self.config.set('data', 'overscan_start_column', str(self.overscan_start_column.value()))

            # Save to file
            self.config.save()

            self.accept()

        except Exception as e:
            from PyQt5.QtWidgets import QMessageBox
            QMessageBox.critical(self, "Error", f"Failed to save settings: {e}")