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

## 2. Installation

The generated installer handles the dependency setup for you.

1. **Run the Installer**:
   * **General**: Download and extract the ZIP archive from the release.
   * **Windows**: Double-click `HERA_Runtime_Installer.exe`.
   * **macOS**: Double-click `HERA_Runtime_Installer.app`.
   * **Linux**: Run the installer executable from the terminal.
2. **Finish Installation**: Follow the on-screen prompts. The installer will automatically download and install the correct MATLAB Runtime if it's missing on your system.

## 3. Usage

### GUI Launcher (Interactive)

To start the application in the standard interactive mode:

* **Windows**: Launch `HERA_Runtime` from the installation directory.
* **macOS**: Double-click the `HERA_Launcher.command` script provided with the release.
* **Linux**: Run `./run_HERA_Runtime.sh <RuntimePath>` from the terminal.

### Command Line Interface (macOS)

For advanced usage, you can run the application directly from the terminal using the `HERA_Launcher.command` script. This allows you to pass arguments for batch processing, testing, or analysis.

**Note:** Ensure you are in the directory containing the launcher script.

#### 1. Batch Processing (Non-Interactive)

Run a full ranking analysis using a configuration file, skipping the UI.
For more details on configuration parameters, see [Configuration & Parameters](Configuration_&_Parameters.md).

```bash
./HERA_Launcher.command configFile "/absolute/path/to/config.json"
```

#### 2. Run Unit Tests

Execute the internal test suite to verify the integrity of the installation.
For more details on the test suite, see the [Testing](../README.md#testing) section in the main documentation.

```bash
./HERA_Launcher.command runtest true
```

#### 3. Run Convergence Analysis

Perform the robust convergence verification study. For more details, refer to the [Convergence Analysis Documentation](Convergence_Analysis.md).

```bash
./HERA_Launcher.command convergence true
```
