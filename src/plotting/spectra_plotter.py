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


def plot_overscan_corrected_image(image: np.ndarray, output_path: str,
                                 title: str = "", cmap: str = 'viridis',
                                 overscan_start_col: int = -1,
                                 original_width: int = None):
    """Save overscan-corrected image plot (without marking overscan region).
    
    Args:
        image: The corrected image (may have overscan region trimmed)
        output_path: Path to save the plot
        title: Plot title
        cmap: Colormap for the image
        overscan_start_col: 1-based column where overscan starts
        original_width: Original image width before trimming (optional)
    """
    fig, ax = plt.subplots(figsize=(12, 8))

    vmin = np.percentile(image, 1)
    vmax = np.percentile(image, 99)

    im = ax.imshow(image, cmap=cmap, origin='lower', vmin=vmin, vmax=vmax)
    
    height, current_width = image.shape
    
    # Set title (no overscan region marking as requested)
    ax.set_title(f"{title}")
    
    ax.set_xlabel("Column (0-based)")
    ax.set_ylabel("Row (0-based)")
    
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Counts")

    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)

    logger.info(f"Saved overscan-corrected image plot to {output_path}")


def plot_overscan_profile(overscan_profile: dict, output_path: str,
                         title: str = "Overscan Profile"):
    """Plot overscan profile (smoothed or fitted curve).
    
    Args:
        overscan_profile: Dictionary containing overscan profile data
            For non-split profiles: 'raw', 'smoothed' or 'fitted', 'method', 
                          'window_length', 'polyorder', etc.
            For split profiles (split=True): 'upper' and 'lower' dicts with same structure,
                          plus 'split': True and 'method'.
        output_path: Path to save the plot
        title: Plot title
    """
    method = overscan_profile.get('method', 'unknown')
    is_split = overscan_profile.get('split', False)

    has_model = False
    if is_split:
        bottom = overscan_profile.get('bottom', {})
        top = overscan_profile.get('top', {})
        has_model = (
            ('fitted' in bottom or 'smoothed' in bottom) or
            ('fitted' in top or 'smoothed' in top)
        )
    else:
        has_model = ('fitted' in overscan_profile) or ('smoothed' in overscan_profile)

    if has_model:
        fig, (ax, ax_res) = plt.subplots(
            2, 1,
            figsize=(10, 8),
            sharex=True,
            gridspec_kw={'height_ratios': [3, 1]}
        )
    else:
        fig, ax = plt.subplots(figsize=(10, 6))
        ax_res = None

    def _get_method_title(profile: dict, method_name: str) -> str:
        display_method = 'mean_savgol' if method_name == 'iraf_savgol' else method_name
        if display_method == 'mean_polynomial':
            poly_type = profile.get('poly_type', 'polynomial')
            poly_order = profile.get('poly_order', 'N/A')
            return f" - {display_method} ({poly_type}, order={poly_order})"
        if display_method == 'mean_savgol':
            poly_order = profile.get('polyorder', 'N/A')
            return f" - {display_method} (order={poly_order})"
        if display_method == 'mean_only':
            return f" - {display_method}"
        return f" - {display_method}"

    method_title = _get_method_title(overscan_profile, method)
    
    if is_split:
        # Handle split profile (top and bottom halves)
        # numpy: image[0,:] = ds9 row 1, image[4135,:] = ds9 row 4136
        # ds9: x=0 (left) = small row, x=total_height-1 (right) = large row
        # 'bottom' = image[:mid_row,:] = small rows (ds9 bottom half)
        # 'top' = image[mid_row:,:] = large rows (ds9 top half)
        bottom = overscan_profile.get('bottom', {})  # image[:mid_row,:] = small rows
        top = overscan_profile.get('top', {})  # image[mid_row:,:] = large rows
        
        n_bottom = len(bottom.get('raw', []))
        n_top = len(top.get('raw', []))
        total_height = n_bottom + n_top
        mid_row = n_bottom
        
        # Bottom half ('bottom' = image[:mid_row,:] = small rows) -> plot on left
        if 'raw' in bottom:
            # x=0 = ds9 row 1, x=n_bottom-1 = ds9 row mid_row
            row_indices = np.arange(n_bottom)
            ax.plot(row_indices, bottom['raw'], 'b.', alpha=0.5, markersize=4)
            
            bottom_model = None
            if 'fitted' in bottom:
                bottom_model = bottom['fitted']
                ax.plot(row_indices, bottom_model, 'g-', linewidth=2)
            elif 'smoothed' in bottom:
                bottom_model = bottom['smoothed']
                ax.plot(row_indices, bottom_model, 'r-', linewidth=2)

            if ax_res is not None and bottom_model is not None:
                residual_bottom = np.asarray(bottom['raw']) - np.asarray(bottom_model)
                ax_res.plot(row_indices, residual_bottom, 'k.', alpha=0.6, markersize=3)
        
        # Top half ('top' = image[mid_row:,:] = large rows) -> plot on right
        if 'raw' in top:
            # x=n_bottom = ds9 row mid_row+1, x=total_height-1 = ds9 row total_height
            row_indices = np.arange(n_bottom, total_height)
            ax.plot(row_indices, top['raw'], 'c.', alpha=0.5, markersize=4)
            
            top_model = None
            if 'fitted' in top:
                top_model = top['fitted']
                ax.plot(row_indices, top_model, 'orange', linewidth=2)
            elif 'smoothed' in top:
                top_model = top['smoothed']
                ax.plot(row_indices, top_model, 'm-', linewidth=2)

            if ax_res is not None and top_model is not None:
                residual_top = np.asarray(top['raw']) - np.asarray(top_model)
                ax_res.plot(row_indices, residual_top, 'k.', alpha=0.6, markersize=3)
        
        # Add vertical line at split boundary
        ax.axvline(
            x=n_bottom,
            color='gray',
            linestyle='--',
            linewidth=1,
            alpha=0.7,
            label=f'Split at row {n_bottom}'
        )
        if ax_res is not None:
            ax_res.axvline(x=n_bottom, color='gray', linestyle='--', linewidth=1, alpha=0.7)

        # For split mode, enrich title details using one half profile if available.
        detail_source = bottom if bottom else top
        method_title = _get_method_title(detail_source, method)
        
    else:
        # Non-split profile (single region)
        n = len(overscan_profile['raw'])
        # ds9: row 1 at bottom (left), row n at top (right)
        row_indices = np.arange(n)  # 0 to n-1 corresponds to ds9 row 1 to row n
        
        # Plot raw overscan values
        ax.plot(row_indices, overscan_profile['raw'], 'b.', alpha=0.5, label='Raw (sigma-clipped mean)')
        
        # Plot smoothed curve if available
        model = None
        if 'smoothed' in overscan_profile:
            model = overscan_profile['smoothed']
            ax.plot(row_indices, overscan_profile['smoothed'], 'r-', linewidth=2, 
                    label=f"Savitzky-Golay smoothed (window={overscan_profile.get('window_length', 'N/A')}, order={overscan_profile.get('polyorder', 'N/A')})")
            method_title = _get_method_title(overscan_profile, method)
        
        # Plot polynomial fitted curve if available
        if 'fitted' in overscan_profile:
            model = overscan_profile['fitted']
            ax.plot(row_indices, overscan_profile['fitted'], 'g-', linewidth=2,
                    label=f"Polynomial fit (order={overscan_profile.get('poly_order', 'N/A')})")
            method_title = _get_method_title(overscan_profile, method)
            
            # If coefficients are available, display them
            if 'coefficients' in overscan_profile:
                coeffs = overscan_profile['coefficients']
                coeff_text = "Coefficients: "
                for i, coeff in enumerate(coeffs[::-1]):  # Reverse for standard polynomial notation
                    if i == 0:
                        coeff_text += f"{coeff:.3f}"
                    elif i == 1:
                        coeff_text += f" + {coeff:.3f}x"
                    else:
                        coeff_text += f" + {coeff:.3f}x^{i}"
                
                ax.text(0.02, 0.02, coeff_text,
                        transform=ax.transAxes, fontsize=9,
                        verticalalignment='bottom',
                        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))

        if ax_res is not None and model is not None:
            residual = np.asarray(overscan_profile['raw']) - np.asarray(model)
            ax_res.plot(row_indices, residual, 'k.', alpha=0.6, markersize=3)
    
    ax.set_xlabel("ds9 Row Number (1=bottom, increasing upward)")
    ax.set_ylabel("Overscan Level (counts)")
    ax.set_title(f"{title}{method_title}")
    
    # Add grid
    ax.grid(True, alpha=0.3)
    
    # Add legend
    ax.legend(loc='best')

    if ax_res is not None:
        ax_res.axhline(y=0.0, color='gray', linestyle='--', linewidth=1, alpha=0.7)
        ax_res.set_ylabel("Residual")
        ax_res.set_xlabel("ds9 Row Number (1=bottom, increasing upward)")
        ax_res.grid(True, alpha=0.3)
    else:
        ax.set_xlabel("ds9 Row Number (1=bottom, increasing upward)")
    
    # Tight layout
    plt.tight_layout()
    
    # Save the plot
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    
    logger.info(f"Saved overscan profile plot to {output_path}")


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
