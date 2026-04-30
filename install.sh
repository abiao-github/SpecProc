#!/bin/bash
# SpecProc Installation Script
# Installs SpecProc into the currently active conda/Python environment.

set -e  # Exit on error

echo "==================================="
echo "SpecProc Installation Script"
echo "==================================="
echo ""

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo "Error: conda is not found. Please install Anaconda or Miniconda first."
    echo "Download from: https://www.anaconda.com/products/distribution"
    exit 1
fi

# Detect the currently active conda environment
CURRENT_ENV="${CONDA_DEFAULT_ENV:-base}"
echo "Detected active conda environment: $CURRENT_ENV"
echo ""

if [ "$CURRENT_ENV" = "base" ]; then
    echo "⚠  You are in the 'base' environment."
    read -p "   Install into 'base'? [y/N] " -n 1 -r
    echo ""
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "Aborted. Please activate your target environment first:"
        echo "  conda activate <your_env>"
        exit 0
    fi
fi

echo "Step 1: Installing dependencies into '$CURRENT_ENV' ..."

# Install PyQt5 via conda (pip version can be problematic)
echo "Installing PyQt5..."
conda install -c conda-forge pyqt5 -y

# Install scientific packages
echo "Installing scientific packages..."
conda install numpy scipy astropy matplotlib -y
echo "✓ Dependencies installed"
echo ""

echo "Step 2: Installing SpecProc (editable mode)..."
echo "Current directory: $(pwd)"
pip install -e .
echo "✓ SpecProc installed in development mode"
echo ""

echo "Step 3: Verifying installation..."
if command -v specproc &> /dev/null; then
    echo "✓ SpecProc command is available"
else
    echo "Warning: specproc command not found in PATH"
    echo "You can still run: python run.py"
fi

echo ""
echo "==================================="
echo "Installation Complete!"
echo "==================================="
echo ""
echo "SpecProc is installed in environment: $CURRENT_ENV"
echo ""
echo "Usage:"
echo ""
echo "  Run GUI mode (default):"
echo "    specproc"
echo "    or"
echo "    specproc --mode gui"
echo ""
echo "  Run CLI mode:"
echo "    specproc --mode cli"
echo ""
echo "  Specify config file:"
echo "    specproc --config /path/to/config.cfg"
echo ""
echo "Note: The 'specproc' command works from any directory!"
echo ""
echo "For more information, see:"
echo "  - README_CN.md (中文)"
echo "  - README.md (English)"
echo ""
