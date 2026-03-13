# Standalone Runtime

HERA can be compiled into a standalone application for macOS, Linux, and
Windows. The build process generates an **installer** that automatically
downloads the required MATLAB Runtime, making it easy to distribute.

> [!IMPORTANT]
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

> [!IMPORTANT]
> **Important Note for macOS Users:**
> Due to macOS security protocols for software distributed outside the App Store, you may need to manually authorize  the
> application **twice**: once for the Installer and once for the Launcher.
>
> **Option 1**
> If macOS prevents the file from opening ("App is damaged" or "Unidentified Developer"):
>
> 1. **Right-click** the file (Installer or Launcher) and select **Open**.
> 2. Attempt to open it. If it still fails:
> 3. Go to **System Settings** -> **Privacy & Security**.
> 4. Scroll down to the **Security** section.
> 5. Click **Open Anyway** (or "Allow") for the blocked application.
> 6. Confirm by clicking **Open**.
> Repeat these steps for both the Installer and the Launcher if prompted (might be every time you want to use it).
>
> **Option 2: Permanent Fix**
> This removes all security warnings for HERA permanently by clearing the quarantine attributes.
>
> 1. Open the **Terminal** app (Press `Cmd + Space`, type "Terminal", and press Enter).
> 2. Type `sudo xattr -cr` (make sure there is a **space** at the end).
> 3. **Drag and drop** the entire extracted `HERA_Runtime_...` folder from Finder into the Terminal window or copy and paste the folder path.
> 4. Press **Enter**.
> 5. Type your login password if prompted (typing will be invisible) and press **Enter**.
> 6. Open the **Terminal** app (Press `Cmd + Space`, type "Terminal", and press Enter).
> 7. Type `sudo xattr -cr` (make sure there is a **space** at the end).
> 8. **Drag and drop** the entire extracted `HERA_Runtime_...` folder from Finder into the Terminal window or copy and paste the folder path.
> 9. Press **Enter**.
> 10. Type your login password if prompted (typing will be invisible) and press **Enter**.
> 11. Now you can double-click the Installer and Launcher without issues.
>
> *Disclaimer: This software and the associated instructions are provided "as is", without warranty of any kind. The author expressly disclaim all liability for any damages, errors, or security risks arising from the use of this software or the execution of system commands. Modifying file attributes skips Apple's security checks; proceed at your own risk.*

## 3. Usage

### GUI Launcher (Interactive)

To start the application in the standard interactive mode:

* **Windows**: Launch `HERA_Runtime` from the installation directory.
* **macOS**: Double-click the `HERA_Launcher.command` script provided with the release.
* **Linux**: Run `./run_HERA_Runtime.sh <RuntimePath>` from the terminal.

### Command Line Interface (macOS)

For advanced usage, you can run the application directly from the terminal using the `HERA_Launcher.command` script. This allows you to pass arguments for batch processing, testing, or analysis.

> [!NOTE]
> Ensure you are in the directory containing the launcher script.

#### 1. Batch Processing (Non-Interactive)

Run a full ranking analysis using a configuration file, skipping the UI.
For more details on configuration parameters, see [Configuration & Parameters](Configuration_&_Parameters.md).

```bash
./HERA_Launcher.command configFile "/absolute/path/to/config.json"
```

#### 2. Run Unit Tests

Execute the internal test suite to verify the integrity of the installation.
For more details on the test suite, see the [Testing](https://lerdmann1601.github.io/HERA-Matlab/#testing) section in the main documentation.

```bash
./HERA_Launcher.command runtest true
```

#### 3. Run Convergence Analysis

Perform the robust convergence verification study. For more details, refer to the [Convergence Analysis Documentation](Convergence_Analysis.md).

```bash
./HERA_Launcher.command convergence true
```
