# Calibration Data Directory

This directory contains all calibration data files for SpecProc.

## Directory Structure

```
calib_data/
├── linelists/              # Lamp emission line catalogs
│   ├── thar-noao.dat      # ThAr lamp lines (Xinglong 2.16m HRS)
│   ├── thar.dat           # Standard ThAr lamp lines
│   ├── fear.txt           # FeAr lamp lines (Willmarth)
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

**Available line-list files**:
- `thar-noao.dat` - Th-Ar line list; recommended for Xinglong 2.16m HRS.
- `thar.dat` - Th-Ar line list (3180-9665 Å).
- `vacthar.dat` - `thar.dat` converted to vacuum wavelengths.
- `fear.txt` - Fe-Ar hollow-cathode line list (3000-10,000 Å, Willmarth).
- `fear.fits` - NOIRLAB/KPNO FeAr atlas spectrum (0.2 Å sampling, linear wavelength WCS) used by `specproc-match-anchors`.
- `cuar.dat` - Cu-Ar line list.
- `henear.dat`, `henearhres.dat`, `idhenear.dat` - He-Ne-Ar line lists, including a high-resolution list and an extended list.
- `vacidhenear.dat` - `idhenear.dat` converted to vacuum wavelengths.
- `ctiohenear.dat`, `ctiohear.dat` - He-Ne-Ar and He-Ar lists used at CTIO.
- `argon.dat`, `krypton.dat`, `xenon.dat` - Infrared emission lines of Ar, Kr, and Xe.
- `skylines.dat`, `lowskylines.dat`, `ohlines.dat` - Night-sky and OH emission-line lists.
- `thorium.dat`, `vacthorium.dat` - Obsolete aliases; use `thar.dat` and `vacthar.dat`, respectively.

The catalogs are plain-text files and their columns vary by list. `airvac.cl` is a wavelength air/vacuum conversion utility, not a line list. Select a catalog appropriate to the lamp and wavelength range; the files are not CSV catalogs.

### Air/vacuum wavelength conversion

SpecProc includes a command-line tool based on the IRAF `airvac.cl` formula (Allen's *Astrophysical Quantities*, 1973). Wavelengths are in Angstroms. From the repository root, run:

```bash
specproc-air2vac input.dat output.dat
specproc-vac2air input.dat output.dat
```

Install or update SpecProc to add both commands to the active Python environment (for a local editable install, run `python -m pip install -e .`). The defaults are `-t 15` °C, `-p 760` mmHg, and `-f 4` mmHg water-vapour pressure; override them with `-t`, `-p`, and `-f`. The tool converts the first numeric column and preserves comments and remaining text. Write to a different output file to keep the original list unchanged.

When running directly from the repository without installing, use `PYTHONPATH=src python -m specproc.utils.airvac input.dat output.dat --type air2vac` or `--type vac2air`.

### Wavelength convention of `thar-noao.dat`

The file has no header stating whether its wavelengths are in air or vacuum. A comparison with `thar.dat` supports the conclusion that `thar-noao.dat` is an air-wavelength list: within the `thar.dat` range (3184.9486-9664.6983 Å), 2,948 wavelength values match `thar.dat` exactly. The IRAF line-list documentation identifies `vacthar.dat` as the vacuum-wavelength conversion of `thar.dat`, so the exact overlap is consistent with both source lists using air wavelengths. This is strong evidence, not an explicit declaration in the `thar-noao.dat` file; verify the original source metadata when available.

### Manual line identification (long-slit spectra)

For the Muztagh-Ata 1.93m slit spectrograph there is no pre-computed calibration. Identify a few obvious lamp lines by hand:

```bash
specproc-identify --lamp output/step5_extraction/<lamp>_spec.fits --linelist calib_data/linelists/fear.dat
```

(`python -m specproc.gui.line_identifier` works without installation; `--lamp` and `--linelist` are optional and can also be loaded from the window.) The window shows the extracted lamp spectrum (top, versus pixel) and the reference list drawn as delta functions on a zero baseline (bottom, versus wavelength). You can click an obvious lamp peak and its matching reference line to add a manual anchor; right-click removes a pair. Set the **low/upper wavelength range** and **Peak SNR threshold**, then press **Auto identify**. Manual anchors are optional: with none, the tool matches the pattern of line spacings (ratios of spacings among triples of the strongest lamp peaks versus catalogue lines, in either dispersion direction) to find a starting solution, discarding solutions whose wavelength coverage falls outside the selected range, then iteratively matches lamp peaks above the SNR threshold. Any manual anchors are retained. You can also tick **Overlay reference lines**: the catalogue is drawn as vertical lines on the lamp spectrum using the current model (or the low/upper range, with **Flip** for decreasing wavelength). Align them by hand with **Shift**/**Stretch** (or the arrow keys: ←/→ shift, ↑/↓ stretch), then press **Fit overlay** to match the peaks and fit. Save the result and set `identification_file` in `[slit]` to it (set `peak_snr` to the same threshold): the pipeline refits using these lines, matches against `full_linelist` and refits the final polynomial. Use a catalogue in the same medium (air or vacuum) as the identified wavelengths.

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
