"""
Runtime health checker and installer utility for HERA-Matlab.

This module provides functionality to verify if the required MATLAB Runtime
(R2025b / v25.2) is correctly installed and configured on the system path.
If the runtime is missing, it provides platform-specific instructions for installation.

Typical usage example:

    python3 -m hera_matlab.install_runtime
"""

import os
import platform
import sys


def check_runtime() -> None:
    """
    Checks if the correct MATLAB Runtime is installed and accessible.
    
    This function performs the following steps:
    1. Detects the current operating system and architecture.
    2. Constructs platform-specific download URLs for the MATLAB Runtime R2025b.
    3. Scans the system's dynamic library path (PATH, LD_LIBRARY_PATH, DYLD_LIBRARY_PATH)
       to verify if the runtime's core library (`mclmcrrt`) is loadable.
    4. Prints a status report to the console.

    If the runtime is missing, it prints error details and a direct download link
    to resolve the issue.
    """
    # -------------------------------------------------------------------------
    # 1. Platform Detection
    # -------------------------------------------------------------------------
    system = platform.system()
    arch = ""
    
    if system == 'Windows':
        bit_str = platform.architecture()[0]
        if bit_str == '64bit':
            arch = 'win64'
        else:
            print(f"Error: {bit_str} Windows is not supported.")
            return
    elif system == 'Linux':
        arch = 'glnxa64'
    elif system == 'Darwin': # macOS
        if platform.mac_ver()[-1] == 'arm64':
             arch = 'maca64'
        else:
             arch = 'maci64'
    else:
        print(f"Error: Operating system {system} is not supported.")
        return

    # -------------------------------------------------------------------------
    # 2. Runtime Verification
    # -------------------------------------------------------------------------
    
    # Official MathWorks download base URL for R2025b
    base_url = "https://ssd.mathworks.com/supportfiles/downloads/R2025b/Release/0/deployment_files/installer/complete"
    
    download_urls = {
        'win64': f"{base_url}/win64/MATLAB_Runtime_R2025b_win64.zip",
        'glnxa64': f"{base_url}/glnxa64/MATLAB_Runtime_R2025b_glnxa64.zip",
        'maci64': f"{base_url}/maci64/MATLAB_Runtime_R2025b_maci64.dmg",
        'maca64': f"{base_url}/maca64/MATLAB_Runtime_R2025b_maca64.dmg.zip"
    }

    # Determine environment variable and library name conventions based on OS
    path_var = 'PATH'
    lib_prefix = ''
    ext = 'dll'
    
    if system == 'Linux':
        path_var = 'LD_LIBRARY_PATH'
        lib_prefix = 'libmw'
        ext = 'so'
    elif system == 'Darwin':
        path_var = 'DYLD_LIBRARY_PATH'
        lib_prefix = 'libmw'
        ext = 'dylib'
    
    # Target Runtime Version: R2025b (v25.2)
    runtime_version_underscore = '25_2'
    runtime_version_dot = '25.2'
    
    # Construct the filename of the library to search for
    file_to_find = ""
    if system == 'Windows':
        file_to_find = f"{lib_prefix}mclmcrrt{runtime_version_underscore}.{ext}"
    elif system == 'Linux':
         file_to_find = f"{lib_prefix}mclmcrrt.{ext}.{runtime_version_dot}"
    elif system == 'Darwin':
         file_to_find = f"{lib_prefix}mclmcrrt.{runtime_version_dot}.{ext}"
         
    # Scan the path for the library
    path_val = os.environ.get(path_var, '')
    found = False
    found_path = ""
    
    if path_val:
        for path_elem in path_val.split(os.pathsep):
            potential_path = os.path.join(path_elem, file_to_find)
            if os.path.isfile(potential_path):
                found = True
                found_path = potential_path
                break
    
    # -------------------------------------------------------------------------
    # 3. User Reporting
    # -------------------------------------------------------------------------
    if found:
        print(f"MATLAB Runtime R2025b found.")
        print(f"Platform: {system} ({arch})")
        print(f"Location: {found_path}")
    else:
        # Provide clear, actionable instructions
        print(f"Error: MATLAB Runtime R2025b not found on '{path_var}'.")
        print(f"Platform: {system} ({arch})")
        print("")
        print(f"Please install the MATLAB Runtime R2025b (v25.2) from:")
        print(f"{download_urls.get(arch, 'https://www.mathworks.com/products/compiler/matlab-runtime.html')}")
        
        if system == 'Darwin':
             print("")
             print(f"Note for macOS: SIP might clear DYLD_LIBRARY_PATH. Try running with 'mwpython'.")
        
        sys.exit(1)

if __name__ == "__main__":
    check_runtime()
