#!/bin/bash
# SpecProc Installation Script
# This script will set up SpecProc in a new conda environment

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

# Environment name
ENV_NAME="specproc"
PYTHON_VERSION="3.8"

echo "Step 1: Creating conda environment..."
conda create -n $ENV_NAME python=$PYTHON_VERSION -y
echo "✓ Conda environment '$ENV_NAME' created with Python $PYTHON_VERSION"
echo ""

echo "Step 2: Installing dependencies..."
source $(conda info --base)/etc/profile.d/conda.sh
conda activate $ENV_NAME

# Install PyQt5
echo "Installing PyQt5..."
conda install -c conda-forge pyqt5 -y

# Install scientific packages
echo "Installing scientific packages..."
conda install numpy scipy astropy matplotlib -y
echo "✓ Dependencies installed"
echo ""

echo "Step 3: Installing SpecProc..."
echo "Current directory: $(pwd)"
pip install -e .
echo "✓ SpecProc installed in development mode"
echo ""

echo "Step 4: Verifying installation..."
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
echo "To use SpecProc:"
echo ""
echo "1. Activate the environment:"
echo "   conda activate $ENV_NAME"
echo ""
echo "2. Run GUI mode (default):"
echo "   specproc"
echo "   or"
echo "   specproc --mode gui"
echo ""
echo "3. Run CLI mode:"
echo "   specproc --mode cli"
echo ""
echo "4. Specify config file:"
echo "   specproc --config /path/to/config.cfg"
echo ""
echo "Note: The 'specproc' command works from any directory!"
echo ""
echo "For more information, see:"
echo "  - INSTALLATION_GUIDE_CN.md (中文)"
echo "  - INSTALLATION_GUIDE.md (English)"
echo "  - README.md"
echo ""
