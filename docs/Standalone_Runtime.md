# Standalone Runtime

HERA can be compiled into a standalone application for macOS, Linux, and
Windows. The build process generates an **installer** that automatically
downloads the required MATLAB Runtime, making it easy to distribute.

> **Download:** A pre-built installer for macOS (Apple Silicon) is available as
> a ZIP archive in the [Releases](https://github.com/lerdmann1601/HERA-Matlab/releases)
> section.

## Building the Installer

To build the installer, you need a MATLAB installation with the **MATLAB
Compiler** toolbox.

1. Open MATLAB and navigate to the project root.
2. Run the build script:

   ```matlab
   cd deploy
   build_HERA_matlab
   ```

3. The artifacts (Installer + ZIP) will be generated in `deploy/output/matlab`.

## Installation and Usage

The generated installer handles the dependency setup for you.

1. **Run the Installer**:
   * **General**: Download and extract the ZIP archive from the release.
   * **Windows**: Double-click `HERA_Runtime_Installer.exe`.
   * **macOS**: Double-click `HERA_Runtime_Installer.app`.
   * **Linux**: Run the installer executable from the terminal.
2. **Follow the Prompts**: The installer will automatically download and
   install the correct MATLAB Runtime if it's missing.
3. **Run HERA**:
   * **Windows**: Launch `HERA_Runtime` from the installation directory.
   * **macOS**: Double-click `HERA_Launcher.command` (provided with the release).
   * **Linux**: Run `./run_HERA_Runtime.sh <RuntimePath>` from the
     terminal.
