# Calibration Data Directory

This directory contains all calibration data files for SpecProc.

## Directory Structure

```
calib_data/
├── linelists/              # Lamp emission line catalogs
│   ├── thar-noao.dat      # ThAr lamp lines (Xinglong 2.16m HRS)
│   ├── thar.dat           # Standard ThAr lamp lines
│   ├── FeAr.dat           # FeAr lamp lines
│   └── ...
└── telescopes/             # Telescope-specific calibration files
    ├── generic/           # Generic configuration template
    └── xinglong216hrs/    # Xinglong 2.16m telescope
        ├── wlcalib_20141103049.fits
        ├── wlcalib_20171202012.fits
        ├── wlcalib_20190905028_A.fits
        └── wlcalib_20211123011_A.fits
```

## Lamp Line Lists

Lamp linelist files are stored in `linelists/` directory.

**Available linelist files**:
- `thar-noao.dat` - ThAr lamp lines (Xinglong 2.16m HRS recommended)
- `thar.dat` - Standard ThAr lamp lines
- `FeAr.dat` - FeAr lamp lines

**Supported lamp types**:
- `ThAr` - Thorium-Argon (most common for echelle spectrographs)
- `FeAr` - Iron-Argon
- `Ar` - Argon
- `Ne` - Neon
- `He` - Helium
- `Fe` - Iron

## Telescope Calibrations

Telescope-specific calibration files are stored in `telescopes/<telescope_name>/` directory.

**Available calibration files for Xinglong 2.16m HRS**:
- `wlcalib_20141103049.fits` - 2014-11-03 04:50
- `wlcalib_20171202012.fits` - 2017-12-02 01:20
- `wlcalib_20190905028_A.fits` - 2019-09-05 02:50 (version A)
- `wlcalib_20211123011_A.fits` - 2021-11-23 01:10 (version A) - **Latest**

## Configuration

Telescope and instrument configuration is set in `default_config.cfg`:

```ini
[telescope]
# Telescope name for calibration lookup
name = xinglong216hrs

# Spectrograph instrument name
instrument = hrs

[telescope.linelist]
# Lamp type to use
linelist_type = ThAr

# Path to linelist files
linelist_path = calib_data/linelists/

# Specific linelist file to use
linelist_file = thar-noao.dat

# Use pre-identified calibration files
use_precomputed_calibration = yes
calibration_path = calib_data/telescopes/xinglong216hrs/

# Specific calibration file to use
calibration_file = wlcalib_20211123011_A.fits
```

## Usage

### Using Pre-computed Calibration (Recommended)

```ini
use_precomputed_calibration = yes
calibration_file = wlcalib_20211123011_A.fits
```

### Re-fitting Wavelength Calibration

```ini
use_precomputed_calibration = no
linelist_file = thar-noao.dat
```

## Adding Custom Calibrations

### Adding a New Lamp Linelist:

1. Create file in `linelists/` directory
2. Follow file format (wavelength, intensity, comment)
3. Configure in `default_config.cfg`: `linelist_file = <filename>`

### Adding a New Telescope:

1. Create directory: `calib_data/telescopes/<telescope_name>/`
2. Place calibration files with proper naming
3. Configure in `default_config.cfg`:
   ```ini
   [telescope]
   name = <telescope_name>
   instrument = <instrument_name>
   ```

## Recommended Configuration for Xinglong 2.16m HRS

```ini
[telescope]
name = xinglong216hrs
instrument = hrs

[telescope.linelist]
linelist_type = ThAr
linelist_file = thar-noao.dat
use_precomputed_calibration = yes
calibration_path = calib_data/telescopes/xinglong216hrs/
calibration_file = wlcalib_20211123011_A.fits
```
