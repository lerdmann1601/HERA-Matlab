========================
HERA Runtime Application
========================

Thank you for downloading HERA (Hierarchical-Compensatory, Effect-Size-Driven Ranking Algorithm).
This package contains the standalone application and installer.

NOTE: NO MATLAB LICENSE REQUIRED!
Neither a commercial MATLAB license nor a manual Java (JRE) installation is needed.
HERA runs entirely on the free, self-contained MATLAB Runtime.

1. Installation
===============
1. Run "HERA_Runtime_Installer" (.app on macOS, .exe on Windows, .install on Linux).
2. Follow the on-screen setup CLI. The installer will automatically download
   and install the free MATLAB Runtime (R2025b) if needed.
3. Note (macOS): At the end of the installation, MathWorks displays a generic
   informational dialog suggesting manual configuration ("DYLD_LIBRARY_PATH").
   You can safely dismiss this message. Manual path configuration is not required,
   as the provided "HERA_Launcher.command" configures all runtime paths automatically.

2. macOS Security Authorization (Gatekeeper)
============================================
HERA is open-source academic software. While the binaries are ad-hoc signed,
they are distributed without a commercial Apple Developer ID certificate and are
not notarized by Apple. Consequently, macOS Gatekeeper intercepts execution by default
with security dialogs ("Unidentified Developer" or "Apple cannot check for malware").

To run the application, you need to authorize the software TWICE: once for the
Installer and once for the Launcher ("HERA_Launcher.command"). Note: A simple
right-click alone is often insufficient on modern macOS versions. You can choose
between a graphical workaround (Option 1) or a permanent terminal solution (Option 2).

Option 1: Authorization via System Settings
-------------------------------------------
1. Right-click the file (Installer or Launcher) and select "Open".
2. If macOS still blocks it ("Unidentified Developer" or "App is damaged"):
   a. Go to "System Settings" -> "Privacy & Security".
   b. Scroll down to the "Security" section.
   c. Click "Open Anyway" (or "Allow") for the blocked file.
   d. Confirm by clicking "Open" (enter user password or Touch ID).
Repeat these steps for both the Installer and the Launcher when prompted.

Option 2: Permanent Solution via Terminal
------------------------------------------
1. Open the Terminal app (via Spotlight or Applications -> Utilities).
2. Run: sudo xattr -cr /path/to/extracted/HERA_folder
   (Tip: You can type "sudo xattr -cr " with a trailing space, then drag and drop
   the extracted folder from Finder into the Terminal window).
3. Press Enter, enter your administrator password, and press Enter.
4. Both the Installer and Launcher can now be opened directly without Gatekeeper prompts.

DISCLAIMER & SECURITY WARNING:
This software and associated instructions are provided "as is", without warranty
of any kind. The authors disclaim all liability for any damages, errors, or security
risks. Modifying file attributes skips Apple's security checks; proceed at your own risk.

3. How to Run
=============
macOS:
------
1. Double-click "HERA_Launcher.command".
2. A Terminal window will open, automatically locate the MATLAB Runtime, and
   launch HERA in interactive mode. Follow the on-screen prompts.

Windows:
--------
1. Double-click "HERA_Runtime.exe".
2. A Command Prompt window will open and launch the application.

Linux:
------
1. Open a terminal.
2. Run: ./run_HERA_Runtime.sh <Path_to_MATLAB_Runtime>
   (e.g., ./run_HERA_Runtime.sh /usr/local/MATLAB/MATLAB_Runtime/R2025b)

4. Command-Line & Non-Interactive Usage (Terminal / Scripts / Pipelines)
========================================================================
To run HERA directly from the command line, automated scripts, or headless
environments without interactive prompts, invoke the runtime script with your
MATLAB Runtime path and a JSON configuration file:

macOS:
/Applications/HERA_Runtime/application/run_HERA_Runtime.sh /Applications/MATLAB/MATLAB_Runtime/R2025b configFile "/path/to/config.json"

Linux:
./run_HERA_Runtime.sh /usr/local/MATLAB/MATLAB_Runtime/R2025b configFile "/path/to/config.json"

Windows:
HERA_Runtime.exe configFile "C:\path\to\config.json"

To run unit tests:
macOS: /Applications/HERA_Runtime/application/run_HERA_Runtime.sh /Applications/MATLAB/MATLAB_Runtime/R2025b runtest true
Linux: ./run_HERA_Runtime.sh <Path_to_Runtime> runtest true
Windows: HERA_Runtime.exe runtest true

For complete documentation, guides, and tutorials, visit:
https://lerdmann1601.github.io/HERA-Matlab/Standalone_Runtime

================================
Built by Lukas von Erdmannsdorff
================================
