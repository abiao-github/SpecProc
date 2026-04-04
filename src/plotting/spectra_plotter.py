"""
Plotting and visualization utilities for spectral data.

Provides functions for creating diagnostic plots and result visualization.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from typing import Optional, List
import logging

logger = logging.getLogger(__name__)


class MatplotlibCanvas(FigureCanvas):
    """Matplotlib canvas widget for PyQt5."""

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        """Initialize canvas."""
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111)
        super().__init__(self.fig)
        self.setParent(parent)

    def plot_spectrum(self, wavelength: np.ndarray, flux: np.ndarray,
                     error: Optional[np.ndarray] = None, title: str = ""):
        """Plot 1D spectrum."""
        self.axes.clear()
        self.axes.plot(wavelength, flux, 'b-', linewidth=0.7, label="Flux")

        if error is not None:
            self.axes.fill_between(wavelength, flux - error, flux + error,
                                   alpha=0.3, color='blue', label="Error")

        self.axes.set_xlabel("Wavelength (Angstrom)")
        self.axes.set_ylabel("Flux")
        self.axes.set_title(title)
        self.axes.legend()
        self.axes.grid(True, alpha=0.3)
        self.fig.tight_layout()
        self.draw()

    def plot_2d_image(self, image: np.ndarray, title: str = "",
                     vmin: Optional[float] = None, vmax: Optional[float] = None):
        """Plot 2D image with colorbar."""
        self.axes.clear()

        if vmin is None:
            vmin = np.percentile(image, 1)
        if vmax is None:
            vmax = np.percentile(image, 99)

        im = self.axes.imshow(image, cmap='viridis', origin='lower',
                             vmin=vmin, vmax=vmax)
        self.axes.set_title(title)
        self.axes.set_xlabel("Column")
        self.axes.set_ylabel("Row")

        cbar = self.fig.colorbar(im, ax=self.axes)
        cbar.set_label("Counts")

        self.fig.tight_layout()
        self.draw()

    def plot_aperture_traces(self, image: np.ndarray, apertures: dict,
                            title: str = ""):
        """Plot 2D image with aperture traces overlaid."""
        self.axes.clear()

        # Plot image
        im = self.axes.imshow(image, cmap='gray', origin='lower')

        # Overlay apertures
        cols = np.arange(image.shape[1])
        colors = plt.cm.rainbow(np.linspace(0, 1, len(apertures)))

        for (aperture_id, aperture), color in zip(apertures.items(), colors):
            center_pos = aperture.get_position(cols)
            self.axes.plot(cols, center_pos, color=color, linewidth=1,
                          label=f"Order {aperture_id}")

            # Add extraction limits
            upper = aperture.get_upper(cols)
            lower = aperture.get_lower(cols)
            self.axes.fill_between(cols, lower, upper, alpha=0.2, color=color)

        self.axes.set_title(title)
        self.axes.set_xlabel("Column")
        self.axes.set_ylabel("Row")
        self.axes.legend(fontsize=8, loc='upper right')

        self.fig.tight_layout()
        self.draw()


def plot_spectrum_to_file(wavelength: np.ndarray, flux: np.ndarray,
                         output_path: str, error: Optional[np.ndarray] = None,
                         title: str = ""):
    """Save spectrum plot to file."""
    fig, ax = plt.subplots(figsize=(12, 6))

    ax.plot(wavelength, flux, 'b-', linewidth=0.7, label="Flux")

    if error is not None:
        ax.fill_between(wavelength, flux - error, flux + error,
                        alpha=0.3, color='blue', label="Error")

    ax.set_xlabel("Wavelength (Angstrom)")
    ax.set_ylabel("Flux")
    ax.set_title(title)
    ax.legend()
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)

    logger.info(f"Saved spectrum plot to {output_path}")


def plot_2d_image_to_file(image: np.ndarray, output_path: str,
                         title: str = "", cmap: str = 'viridis'):
    """Save 2D image plot to file."""
    fig, ax = plt.subplots(figsize=(10, 8))

    vmin = np.percentile(image, 1)
    vmax = np.percentile(image, 99)

    im = ax.imshow(image, cmap=cmap, origin='lower', vmin=vmin, vmax=vmax)
    ax.set_title(title)
    ax.set_xlabel("Column")
    ax.set_ylabel("Row")

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Counts")

    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)

    logger.info(f"Saved 2D image plot to {output_path}")


def plot_background_residuals(science_image: np.ndarray, background: np.ndarray,
                             output_path: str):
    """Plot background estimation and residuals."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    # Original
    im0 = axes[0].imshow(science_image, cmap='viridis', origin='lower')
    axes[0].set_title("Original Image")
    fig.colorbar(im0, ax=axes[0])

    # Background
    im1 = axes[1].imshow(background, cmap='viridis', origin='lower')
    axes[1].set_title("Background Model")
    fig.colorbar(im1, ax=axes[1])

    # Residual
    residual = science_image - background
    im2 = axes[2].imshow(residual, cmap='viridis', origin='lower')
    axes[2].set_title("Residual (after background subtraction)")
    fig.colorbar(im2, ax=axes[2])

    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)

    logger.info(f"Saved background residuals plot to {output_path}")


def plot_wavelength_calibration(line_pixels: np.ndarray, line_wavelengths: np.ndarray,
                               fitted_wavelengths: np.ndarray,
                               output_path: str):
    """Plot wavelength calibration fit."""
    fig, axes = plt.subplots(2, 1, figsize=(12, 8))

    # Wavelength vs pixel
    axes[0].plot(line_pixels, line_wavelengths, 'bo', label="Detected lines")
    axes[0].plot(line_pixels, fitted_wavelengths, 'r-', linewidth=2,
                label="Polynomial fit")
    axes[0].set_xlabel("Pixel")
    axes[0].set_ylabel("Wavelength (Angstrom)")
    axes[0].set_title("Wavelength Calibration")
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    # Residuals
    residuals = line_wavelengths - fitted_wavelengths
    axes[1].plot(line_pixels, residuals, 'go')
    axes[1].axhline(0, color='r', linestyle='--')
    axes[1].set_xlabel("Pixel")
    axes[1].set_ylabel("Residual (Angstrom)")
    axes[1].set_title(f"Calibration Residuals (RMS: {np.std(residuals):.4f} Å)")
    axes[1].grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)

    logger.info(f"Saved wavelength calibration plot to {output_path}")
