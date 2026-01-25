#!/bin/bash
set -e

# BUILD_AND_PREP_PYPI - Helper script to build and prepare HERA for PyPI.
#
# Syntax:
#   ./deploy/build_and_prep_pypi.sh
#
# Description:
#   This script automates the process of building the HERA Python package
#   locally and preparing the distribution files (.whl, .tar.gz) for
#   upload to a GitHub Release.
#
# Workflow:
#   1.  MATLAB Build: Triggers 'build_HERA_python' if MATLAB is available.
#   2.  Verification: Checks if the package was generated correctly.
#   3.  Distribution: Runs 'python3 -m build' to generate distribution artifacts.
#
# Inputs:
#   None (Relies on 'build_HERA_python.m' in the same directory)
#
# Outputs:
#   Distribution files are placed in 'deploy/dist'.
#
# Author: Lukas von Erdmannsdorff
#

# Define paths
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
BUILD_OUTPUT_DIR="$PROJECT_ROOT/deploy/output/python"
DIST_DIR="$PROJECT_ROOT/deploy/dist"

echo "======================"
echo "HERA PyPI Build Helper"
echo "======================"

# 1. MATLAB Build
# Check if MATLAB is available to run the build automatically
if ! command -v matlab &> /dev/null; then
    echo "Warning: 'matlab' command not found. Assuming you have already run the MATLAB build manually."
    echo "If not, please run 'build_HERA_python()' in MATLAB first."
else
    echo "Starting MATLAB build (this might take a while)..."
     matlab -batch "addpath('$SCRIPT_DIR'); try, build_HERA_python(); catch e, disp(e.message); exit(1); end; exit(0);"
fi

# 2. Verification
# Verify output exists (look for setup.py which indicates the package root)
if [ ! -f "$BUILD_OUTPUT_DIR/setup.py" ]; then
    echo "Error: 'setup.py' not found in $BUILD_OUTPUT_DIR"
    echo "Did the MATLAB build complete successfully?"
    exit 1
fi

PKG_DIR="$BUILD_OUTPUT_DIR"

echo "Found package root at: $PKG_DIR"

# 2b. Preparation & Metadata Injection
# Run the preparation script to correct metadata and inject runtime checks
echo "Running pre-PyPI preparation script..."
if [ -f "$PROJECT_ROOT/.github/scripts/prepare_pypi.py" ]; then
    python3 "$PROJECT_ROOT/.github/scripts/prepare_pypi.py" "$PKG_DIR"
else
    echo "Warning: Preparation script not found at .github/scripts/prepare_pypi.py"
fi

# 3. Distribution
# Build Python Distribution
echo "Building Python Distribution (sdist and wheel)..."
rm -rf "$DIST_DIR"
mkdir -p "$DIST_DIR"

cd "$PKG_DIR"

# Set up a virtual environment to avoid PEP 668 "externally-managed-environment" errors
VENV_DIR="$SCRIPT_DIR/.venv_build"
echo "Creating/Using virtual environment at $VENV_DIR..."
python3 -m venv "$VENV_DIR"
source "$VENV_DIR/bin/activate"

# Ensure build tool is installed in the venv
echo "Installing 'build' tool in virtual environment..."
pip3 install --upgrade build

# Build
echo "Building package..."
python3 -m build --outdir "$DIST_DIR"

# Cleanup (optional, but good to deactivate)
deactivate

echo "========================================"
echo "Build Complete!"
echo "Artifacts are ready in: $DIST_DIR"
echo "========================================"
echo "Next Steps:"
echo "1. Create a new Release on GitHub (with a new tag, e.g., v1.0.0)."
echo "2. Upload the files from '$DIST_DIR' to that release."
echo "3. Go to GitHub Actions -> 'Publish to PyPI' -> Run Workflow."

