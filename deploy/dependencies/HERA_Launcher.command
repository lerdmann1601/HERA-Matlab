#!/bin/bash
# HERA Launcher for macOS
# This script launches HERA_Runtime in a Terminal window.

# Get the directory where this script is located
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

echo "====================="
echo "HERA Launcher (macOS)"
echo "====================="
echo ""

# ---------------------------------------------------------
# 1. Locate the HERA application execution script
# ---------------------------------------------------------
echo "Searching for HERA application..."

# Define potential paths where the app might be installed
# The installer typically installs to /Applications/HERA_Runtime/application
POSSIBLE_APP_PATHS=(
    "$DIR/run_HERA_Runtime.sh"                                         # Portable / Local
    "/Applications/HERA_Runtime/application/run_HERA_Runtime.sh"       # Default System Install
    "$HOME/Applications/HERA_Runtime/application/run_HERA_Runtime.sh"  # Default User Install
)

APP_SCRIPT=""

for script_path in "${POSSIBLE_APP_PATHS[@]}"; do
    if [ -f "$script_path" ]; then
        APP_SCRIPT="$script_path"
        break
    fi
done

if [ -z "$APP_SCRIPT" ]; then
    echo " -> Error: Could not find 'run_HERA_Runtime.sh'."
    echo "    Make sure HERA is installed or this launcher is in the application folder."
    echo ""
    read -p "Press Enter to exit..."
    exit 1
else
    echo " -> Found App: $APP_SCRIPT"
fi
echo ""

# ---------------------------------------------------------
# 2. Locate the MATLAB Runtime
# ---------------------------------------------------------
echo "Searching for MATLAB Runtime..."

# Define potential Runtime paths
POSSIBLE_RUNTIME_PATHS=(
    "/Applications/MATLAB_R2025b.app"
    "/Applications/MATLAB/MATLAB_Runtime/R2025b"
    "/Applications/MATLAB/MATLAB_Runtime/v918" 
)

RUNTIME_PATH=""

# Check for existing MATLAB intallation 
for path in "${POSSIBLE_RUNTIME_PATHS[@]}"; do
    if [ -d "$path" ]; then
        RUNTIME_PATH="$path"
        break
    fi
done

# If not found, try to find any MATLAB app
if [ -z "$RUNTIME_PATH" ]; then
    FOUND_MATLAB=$(find /Applications -maxdepth 1 -name "MATLAB*.app" 2>/dev/null | head -n 1)
    if [ -n "$FOUND_MATLAB" ]; then
        RUNTIME_PATH="$FOUND_MATLAB"
    fi
fi

# If still not found, ask user
if [ -z "$RUNTIME_PATH" ]; then
    echo " -> Could not automatically find MATLAB or MATLAB Runtime."
    echo ""
    echo "Please drag and drop your MATLAB/Runtime installation folder here and press Enter:"
    read USER_PATH
    RUNTIME_PATH="$USER_PATH"
else
    echo " -> Found Runtime: $RUNTIME_PATH"
fi

echo ""
echo "Launching HERA..."
echo ""

# ---------------------------------------------------------
# 3. Execute
# ---------------------------------------------------------

"$APP_SCRIPT" "$RUNTIME_PATH"

# Keep terminal open if it crashes immediately
if [ $? -ne 0 ]; then
    echo ""
    echo "HERA exited with an error."
    read -p "Press Enter to close..."
fi
