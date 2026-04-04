# SpecProc Quick Reference

## Commands

### Installation

```bash
# Quick install
./install.sh

# Manual install
conda create -n specproc python=3.8
conda activate specproc
pip install -e .
```

### Running

```bash
# GUI mode (default)
specproc

# CLI mode
specproc --mode cli

# With custom config
specproc --config custom.cfg
```

### Help

```bash
specproc --help
specproc --version
```

## Configuration Files

### Default Config
- Location: `default_config.cfg`
- Contains all default settings

### User Config
- Location: `specproc.cfg` (optional)
- Overrides default settings

### Key Config Sections

```ini
[data]
rawpath = ./rawdata

[telescope]
name = xinglong216hrs
instrument = hrs

[telescope.linelist]
linelist_type = ThAr
linelist_file = thar-noao.dat
use_precomputed_calibration = yes
calibration_file = wlcalib_20211123011_A.fits
```

## Data Directory Structure

```
SpecProc/
├── rawdata/              # Input FITS files
├── output/
│   ├── midpath/         # Intermediate results
│   ├── spectra/         # Final 1D spectra
│   └── figures/         # Diagnostic plots
├── calib_data/
│   ├── linelists/       # Lamp line lists
│   └── telescopes/      # Telescope calibrations
└── specproc.cfg          # User config (optional)
```

## Processing Stages

| Stage | Description | Command (CLI) |
|--------|-------------|----------------|
| 0 | Overscan correction | Select step 0 |
| 1 | Bias subtraction | Select step 1 |
| 2 | Flat fielding & order tracing | Select step 2 |
| 3 | Background subtraction | Select step 3 |
| 4 | Cosmic ray correction | Select step 4 |
| 5 | Spectrum extraction | Select step 5 |
| 6 | Wavelength calibration | Auto (after extraction) |
| 7 | De-blazing | Auto (after calibration) |

## GUI Quick Start

1. Launch: `specproc`
2. Select files:
   - Bias files
   - Flat files
   - Calibration files (ThAr)
   - Science files
3. Click "Run Full Pipeline"

## CLI Quick Start

```bash
# Activate environment
conda activate specproc

# Run CLI
specproc --mode cli

# Follow prompts:
# - Select bias files (e.g., 1-3)
# - Select flat files (e.g., 4-6)
# - Select science files (e.g., 8-9)
# - Select calibration file (e.g., 7)
# - Select stages (Enter for all)
```

## Common Issues

### ImportError: No module named 'PyQt5'
```bash
pip install PyQt5
```

### specproc command not found
```bash
conda activate specproc
pip install -e .
```

### Configuration file not found
```bash
# Copy default config
cp default_config.cfg specproc.cfg
```

## File Naming Conventions

### Calibration Files
- Line lists: `<lamp>_lines_<min>_<max>.csv`
- Calibration: `<instrument>_wlcalib_<YYYYMMDDHH>_<suffix>.fits`

### Output Files
- 1D spectra: `<science_filename>_ods.fits`
- Diagnostic figures: `<step>_<description>.png`

## Documentation

- [README.md](README.md) - Main documentation
- [INSTALLATION_GUIDE.md](INSTALLATION_GUIDE.md) - Installation guide
- [QUICK_START.md](QUICK_START.md) - Quick start
- [PIPELINE_FLOWCHART.md](PIPELINE_FLOWCHART.md) - Workflow
- [calib_data/README.md](calib_data/README_CN.md) - Calibration data
