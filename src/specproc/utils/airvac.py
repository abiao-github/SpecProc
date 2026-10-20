"""Air <-> vacuum wavelength conversion for line-list files.

Port of the IRAF ``airvac.cl`` script. The refractive index formula is from
Allen's Astrophysical Quantities (1973). Wavelengths are in Angstroms.

Usage:
    specproc-air2vac INPUT OUTPUT [-t T] [-p P] [-f F]
    specproc-vac2air INPUT OUTPUT [-t T] [-p P] [-f F]
    python -m specproc.utils.airvac INPUT OUTPUT --type air2vac|vac2air
"""

import argparse
import sys

import numpy as np


def refractive_index(wavelength, t=15.0, p=760.0, f=4.0):
    """Refractive index of air at the given vacuum-like wavelength (Angstrom).

    t: temperature (C), p: pressure (mmHg), f: water vapour pressure (mmHg).
    """
    w = np.asarray(wavelength, dtype=float)
    b = (10000.0 / w) ** 2
    a = 64.328 + 29498.1 / (146.0 - b) + 255.4 / (41.0 - b)
    a = a * p * (1.0 + (1.049 - 0.0157 * t) * 1e-6 * p) / (720.883 * (1.0 + 0.003661 * t))
    a = a - (0.0624 - 0.000680 * b) / (1.0 + 0.003661 * t) * f
    return 1.0 + a / 1e6


def air_to_vac(wavelength, t=15.0, p=760.0, f=4.0):
    """Convert air wavelengths to vacuum wavelengths (same as IRAF airvac)."""
    return refractive_index(wavelength, t, p, f) * np.asarray(wavelength, dtype=float)


def vac_to_air(wavelength, t=15.0, p=760.0, f=4.0):
    """Convert vacuum wavelengths to air wavelengths (same as IRAF airvac)."""
    return np.asarray(wavelength, dtype=float) / refractive_index(wavelength, t, p, f)


def convert_line(line, conv_type="air2vac", t=15.0, p=760.0, f=4.0):
    """Convert the first column of one line; keep the remaining text as is.

    Blank lines, comments and lines not starting with a number are unchanged.
    """
    parts = line.split(None, 1)
    if not parts:
        return line
    try:
        w = float(parts[0])
    except ValueError:
        return line
    func = air_to_vac if conv_type == "air2vac" else vac_to_air
    new_w = float(func(w, t, p, f))
    label = parts[1].strip() if len(parts) > 1 else ""
    return "%.4f %s" % (new_w, label) if label else "%.4f" % new_w


def convert_file(input_path, output_path, conv_type="air2vac", t=15.0, p=760.0, f=4.0):
    """Convert a line-list file; returns the number of converted lines."""
    if conv_type not in ("air2vac", "vac2air"):
        raise ValueError("conv_type must be 'air2vac' or 'vac2air'")
    count = 0
    with open(input_path, "r") as fin, open(output_path, "w") as fout:
        fout.write("# airvac (%s, %s, type=%s, t=%.6g, p=%.6g, f=%.6g)\n"
                   % (input_path, output_path, conv_type, t, p, f))
        for raw in fin:
            line = raw.rstrip("\r\n")
            new = convert_line(line, conv_type, t, p, f)
            if new is not line and new != line:
                count += 1
            fout.write(new + "\n")
    return count


def main(argv=None, conv_type=None):
    prog = "python -m specproc.utils.airvac" if conv_type is None else "specproc-" + conv_type
    direction = {
        "air2vac": ("air", "vacuum"),
        "vac2air": ("vacuum", "air"),
    }
    parser = argparse.ArgumentParser(
        prog=prog,
        description=("Convert a line-list file between air and vacuum wavelengths "
                     "(Angstrom). Comments and text after the first column are preserved.")
        if conv_type is None else
        "Convert line-list wavelengths from %s to %s (Angstrom). "
        "Comments and text after the first column are preserved."
        % direction[conv_type])
    parser.add_argument("input", help="input line-list file")
    parser.add_argument("output", help="output line-list file")
    if conv_type is None:
        parser.add_argument("--type", dest="conv_type", default="air2vac",
                            choices=["air2vac", "vac2air"],
                            help="conversion (default: air2vac)")
    parser.add_argument("-t", type=float, default=15.0, help="temperature in C (default: 15)")
    parser.add_argument("-p", type=float, default=760.0, help="pressure in mmHg (default: 760)")
    parser.add_argument("-f", type=float, default=4.0,
                        help="water vapour pressure in mmHg (default: 4)")
    args = parser.parse_args(argv)

    try:
        conversion = conv_type or args.conv_type
        n = convert_file(args.input, args.output, conversion, args.t, args.p, args.f)
    except (OSError, ValueError) as exc:
        print("Error: %s" % exc, file=sys.stderr)
        return 1
    print("Converted %d wavelengths (%s): %s" % (n, conversion, args.output))
    return 0


def air2vac_main():
    """Command-line entry point for air-to-vacuum conversion."""
    return main(conv_type="air2vac")


def vac2air_main():
    """Command-line entry point for vacuum-to-air conversion."""
    return main(conv_type="vac2air")


if __name__ == "__main__":
    sys.exit(main())
