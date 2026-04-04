# SpecProc: PyQt GUI for Echelle Spectrograph FITS Data Reduction

A complete PyQt-based graphical interface for reducing echelle spectrograph FITS data. Inspired by the [gamse](https://github.com/wangleon/gamse) package.

## Features

- **Complete Reduction Pipeline**: 8-stage automated spectral reduction
  1. Overscan correction - Overscan改正 (所有图像)
  2. Bias subtraction - 本底减除 (均值/中值合并)
  3. Flat fielding & order tracing - 平场改正与阶序追踪
  4. Background subtraction - 背景扣除
  5. Cosmic ray correction - 宇宙线去除 (仅科学图像)
  6. 1D spectrum extraction - 一维谱提取
  7. Wavelength calibration - 波长定标 (对提取的 1D 光谱应用)
  8. De-blazing - Blaze 函数改正

- **Interactive GUI**: PyQt5-based user interface with:
  - File management for bias, flat, and science frames
  - Real-time progress tracking
  - Processing logs and diagnostics
  - One-click pipeline execution or step-by-step processing

- **Configuration-Driven**: INI-format configuration files for easy parameter adjustment

- **Generic Spectrograph Support**: Configurable for different echelle spectrographs

## Installation

### Quick Install (Recommended)

```bash
# Use automated installation script
cd SpecProc
./install.sh

# Activate environment and run
conda activate specproc
python run.py
```

### Manual Installation

#### Prerequisites

- Python 3.8+
- Anaconda or Miniconda (recommended)

#### Step 1: Create Conda Environment

```bash
conda create -n specproc python=3.8
conda activate specproc
```

#### Step 2: Install Dependencies

```bash
# Option 1: Using conda (recommended)
conda install -c conda-forge pyqt5 numpy scipy astropy matplotlib

# Option 2: Using pip
pip install -r requirements.txt
```

#### Step 3: Install SpecProc

```bash
cd SpecProc
pip install -e .
```

#### Step 4: Verify Installation

```bash
specproc --help
python run.py --help
```

### Running SpecProc

After installation, you can use the `specproc` command from any directory.

#### GUI Mode (Graphical Interface) - Default

```bash
specproc              # GUI mode (default)
# or
specproc --mode gui
# or
specproc --config custom.cfg  # with custom config
```

#### CLI Mode (Command Line)

```bash
specproc --mode cli
# or
specproc --mode cli --config custom.cfg
```

**Note**: The `specproc` command works from any directory. You don't need to be in the SpecProc directory to run it.

For detailed installation instructions, see:
- `INSTALLATION_GUIDE_CN.md` (中文安装指南)
- `INSTALLATION_GUIDE.md` (English installation guide)

## Project Structure

```
SpecProc/
├── src/
│   ├── gui/                       # GUI components
│   │   ├── main_window.py         # Main application window
│   │   └── widgets.py             # Custom Qt widgets
│   ├── core/                      # Processing pipeline
│   │   ├── data_structures.py     # Core data classes
│   │   ├── bias_correction.py     # Stage 1
│   │   ├── flat_fielding.py       # Stage 2
│   │   ├── wave_calibration.py    # Stage 3
│   │   ├── background_removal.py  # Stage 4
│   │   ├── extraction.py          # Stage 5
│   │   └── processing_pipeline.py # Orchestrator
│   ├── config/
│   │   └── config_manager.py      # Configuration management
│   ├── utils/
│   │   ├── fits_io.py             # FITS file I/O
│   │   └── image_processing.py    # Image utilities
│   ├── plotting/
│   │   └── spectra_plotter.py     # Visualization utilities
│   └── main.py                    # Application entry point
├── data/
│   ├── telescope/                 # Telescope configurations
│   ├── calibration/               # Calibration data
│   └── linelists/                 # Spectral line lists
├── config_default.cfg             # Default configuration
└── run.py                         # Launcher script
```

## Usage

### Starting the Application

```bash
python run.py
```

### Processing Workflow

1. **Select Input Files**:
   - Add bias frames (multiple)
   - Add flat frames (multiple)
   - Select wavelength calibration (ThAr/Ar/Ne) frame
   - Select science image to reduce

2. **Configure Processing** (optional):
   - Open configuration file to adjust parameters
   - Modify extraction method, polynomial orders, etc.

3. **Run Pipeline**:
   - Click "Run All" to execute complete reduction
   - View progress and logs in real-time

4. **Review Results**:
   - Check processing log for diagnostic information
   - Output spectra saved to `output/` directory

## Configuration

Edit `config_default.cfg` or create custom configuration files to adjustparameters:

- **Bias correction**: Combination method, sigma clipping
- **Flat fielding**: Order detection thresholds
- **Wavelength calibration**: Polynomial orders, line list
- **Extraction**: Sum vs. optimal extraction, aperture limits
- **Background**: 2D polynomial fitting order

## Core Processing Pipeline

### 1. Overscan Correction
- Removes CCD readout bias from overscan regions
- **Applied to ALL image types**: bias, flat, ThAr, science
- Must be the first step in the pipeline

### 2. Bias Correction
- Combines multiple bias frames (median or mean)
- **Note**: Bias frames are short exposure (0s), no cosmic ray correction needed
- Applies bias subtraction to all science frames

### 3. Flat Fielding & Order Tracing
- Combines multiple flat field images
- **Note**: Cosmic rays removed during combination, no separate cosmic ray correction needed
- Detects echelle orders using cross-correlation
- Extracts spatial profiles and sensitivity map

### 4. Background Subtraction
- Estimates 2D inter-order background
- Fits polynomial model to background pixels
- Subtracts background from science image

### 5. Cosmic Ray Correction
- **Applied to science frames ONLY** (long exposure, single frame)
- **Why only science?**
  - Bias: 0s exposure, no cosmic rays
  - Flat: multiple frames combined, cosmic rays removed during combination
  - ThAr: lamp spectra, typically no cosmic ray correction needed
  - Science: long exposure (minutes to hours), single frame, cosmic rays must be removed
- Uses sigma-threshold detection and median filter interpolation
- **Placement**: After background subtraction, before spectrum extraction
- Configurable via `cosmic_sigma` and `cosmic_window` parameters

### 6. Spectrum Extraction
- Extracts 1D spectra using:
  - **Sum extraction**: Simple aperture sum
  - **Optimal extraction**: Gaussian-weighted (Horne 1986)
- Computes extraction errors
- **Note**: Extraction is done in pixel space first, producing 1D spectra in pixel units

### 7. Wavelength Calibration
- **First**: Calibrates ThAr/Ar/Ne frame to create wavelength solution
  - Extracts 1D spectrum from calibration frame (ThAr)
  - Identifies emission lines
  - Fits 2D wavelength polynomial (xorder, yorder): λ(x,y) = Σ p_ij·x^i·y^j
- **Then**: Applies wavelength calibration to extracted science 1D spectra
  - Converts pixel coordinates to wavelength units
  - **Critical**: Must be done AFTER spectrum extraction
  - Cannot be done on 2D images - only works on extracted 1D spectra

### 8. De-blazing
- Divides wavelength-calibrated spectra by blaze function from flat field
- Removes wavelength-dependent instrument response
- Applied after wavelength calibration

## Output Files

Processed data stored in:
- `output/`: Final 1D spectra (FITS)
- `midpath/`: Intermediate processing results
- `figures/`: Diagnostic plots (PNG/PDF)

## Key Data Structures

### ApertureSet
Collection of aperture (order) locations with polynomial traces

### Spectrum
1D spectrum with wavelength, flux, error arrays

### FlatField
Flat field calibration including sensitivity map and order information

### WaveCalib
Wavelength solution with polynomial coefficients and line identifications

## Testing

Test FITS data available in `/Users/abiao/Documents/资料存档/软件工具/gitrepo/20241102_hrs/`

## Development Notes

- Processing modules designed for easy extension to new stages
- All intermediate results saved for debugging and reproducibility
- Configuration-driven architecture enables instrument customization
- Multi-threaded GUI prevents freezing during long operations

## References

-gamse package: https://github.com/wangleon/gamse
- Horne, K. 1986, PASP 98, 609 (optimal extraction method)

## License

[Specify License]

## Authors

SpecProc Development Team
