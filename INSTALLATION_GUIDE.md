# SpecProc Installation and Usage Guide

## System Requirements

- Python 3.8 or higher
- Anaconda or Miniconda (recommended)

## Installation

### Quick Install (Recommended)

```bash
# Navigate to SpecProc directory
cd SpecProc

# Run installation script
./install.sh

# Activate environment
conda activate specproc

# Launch application
specproc
```

### Manual Installation

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
pip install -e .
```

#### Step 4: Verify Installation

```bash
specproc --help
specproc --version
```

## Usage

### GUI Mode (Default)

```bash
specproc
# or
specproc --mode gui
```

### CLI Mode

```bash
specproc --mode cli
```

### With Configuration File

```bash
specproc --config /path/to/config.cfg
```

## Data Preparation

### Working Directory Setup

**Important**: SpecProc should be run in your working directory, NOT in the source code directory.

### Correct Workflow

```bash
# 1. Create your working directory (e.g., for an observation project)
mkdir -p /my/project/2024/obs1
cd /my/project/2024/obs1

# 2. Create subdirectories for data processing
mkdir -p rawdata output/midpath output/spectra output/figures

# 3. Copy/move FITS data files to rawdata directory
cp /somewhere/bias_*.fits ./rawdata/
cp /somewhere/flat_*.fits ./rawdata/
cp /somewhere/science_*.fits ./rawdata/

# 4. Create user config file (optional)
# Copy default config from SpecProc source directory
cp /path/to/SpecProc/default_config.cfg ./specproc.cfg

# 5. Run SpecProc in your working directory
specproc --config ./specproc.cfg
```

### Directory Structure

**Working directory structure** (where you process data):
```
/my/project/2024/obs1/       # Your working directory
├── rawdata/                  # Input FITS files
│   ├── bias_*.fits
│   ├── flat_*.fits
│   ├── thar_*.fits
│   └── science_*.fits
├── output/                    # Processing results (auto-created)
│   ├── midpath/              # Intermediate results
│   ├── spectra/              # Final 1D spectra
│   └── figures/              # Diagnostic plots
├── specproc.cfg              # User config file (optional)
└── ...
```

**Note**:
- ❌ Do NOT run `specproc` in SpecProc source directory (`/path/to/SpecProc`)
- ✅ Run `specproc` in your working directory
- ✅ `rawdata` and `output` will be created in your working directory

## Configuration

### Default Configuration

Default config is in `default_config.cfg`.

### Custom Configuration

```bash
cp default_config.cfg specproc.cfg
nano specproc.cfg
```

### Key Settings

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

## Processing Pipeline

### 8 Stages

1. Overscan correction
2. Bias subtraction
3. Flat fielding & order tracing
4. Background subtraction
5. Cosmic ray correction (science only)
6. 1D spectrum extraction
7. Wavelength calibration
8. De-blazing

### Output

```
output/
├── midpath/         # Intermediate results
├── spectra/         # Final 1D spectra
└── figures/         # Diagnostic plots
```

## Troubleshooting

### Command not found

```bash
conda activate specproc
pip install -e .
```

### Configuration issues

Ensure config file path is correct.

### Data path issues

Set correct path in config file:
```ini
[data]
rawpath = /absolute/path/to/rawdata
```

## More Information

- [README.md](README.md) - Main documentation
- [QUICK_START.md](QUICK_START.md) - Quick start guide
- [PIPELINE_FLOWCHART.md](PIPELINE_FLOWCHART.md) - Processing workflow
- [calib_data/README.md](calib_data/README_CN.md) - Calibration data
