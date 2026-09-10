# Standalone Runtime

HERA can be compiled into a standalone application for macOS, Linux, and Windows. It allows researchers to use the complete feature set of HERA without requiring a MATLAB installation or license.

> [!NOTE]
> **Key Facts for End Users:**
>
> * **Zero License Requirement**: You do **not** need a MATLAB license or commercial MATLAB installation.
> * **No Manual JRE Installation**: You do **not** need to install a Java Runtime Environment (JRE) manually; all internal runtime dependencies are handled automatically by the installer.
> * **Automatic Setup**: The installer automatically downloads and configures the free, self-contained MATLAB Runtime.
> * **Pre-Built Releases**: A pre-built installer is currently provided for **macOS (Apple Silicon)** in the [Releases](https://github.com/lerdmann1601/HERA-Matlab/releases) section. For **Windows and Linux**, non-MATLAB users should use the native [Python Integration](Python_Integration.md) (`hera-matlab`), which is pre-compiled, runs completely license-free on all operating systems, and requires no build steps. *(Compiling the standalone application from source on Windows or Linux requires a local MATLAB installation with the MATLAB Compiler toolbox).*

---

## 1. Installation Guide (End Users)

### Step 1: Download & Extract

1. Navigate to the [GitHub Releases](https://github.com/lerdmann1601/HERA-Matlab/releases) page.
2. Download the pre-built ZIP archive for macOS: `HERA_Runtime_<version>_maca64.zip`.
3. Extract the ZIP archive (double-click in Finder). The extracted folder contains:
   * `HERA_Runtime_Installer_<version>.app` (The setup CLI)
   * `HERA_Launcher.command` (The application launcher for macOS)
   * `readme.txt` (Quick reference guide)

> [!NOTE]
> Users on Windows or Linux who do not have MATLAB should refer to the [Python Integration Guide](Python_Integration.md).
> Developers with a MATLAB Compiler license can build standalone executables for Windows and Linux as described in the [Developer Build section](#4-building-the-installer-from-source-for-developers)).

---

### Step 2: Run the Installer

1. Double-click the installer for your system:
   * **macOS**: `HERA_Runtime_Installer_<version>.app`
   * **Windows**: `HERA_Runtime_Installer_<version>.exe`
   * **Linux**: `./HERA_Runtime_Installer_<version>.install`
2. Follow the on-screen setup CLI. The installer will automatically download and install the correct MATLAB Runtime (R2025b) if it is missing on your system.

> [!IMPORTANT]
> **Note on the Post-Installation MathWorks Notice (macOS):**
>
> When the installer finishes, MathWorks automatically displays a standard informational dialog suggesting manual path configuration:
>
> ```text
> Your installation may require additional configuration steps...
> If the environment variable DYLD_LIBRARY_PATH is undefined, set it to...
> ```
>
> You can safely dismiss this dialog. Manual environment configuration is not necessary because the included `HERA_Launcher.command` script resolves all paths and runtime libraries automatically.

---

### Step 3: macOS Security Authorization (Gatekeeper)

HERA is open-source academic software. While the binaries are ad-hoc signed, they are distributed without a commercial Apple Developer ID certificate and are not notarized by Apple. Consequently, macOS Gatekeeper intercepts execution by default with a security dialog (e.g., *"App is damaged"*, *"Unidentified Developer"*, or *"macOS cannot verify that this app is free from malware"*).

To run the application, macOS requires you to authorize the software **twice**: once for the **Installer** and once for the **Launcher** (`HERA_Launcher.command`). You can choose between a **graphical workaround** (Option 1: via System Settings for individual sessions) or a **permanent terminal solution** (Option 2: clearing quarantine attributes):

> [!IMPORTANT]
> On recent macOS versions (such as Sequoia, Sonoma, or Ventura), Gatekeeper may require explicit authorization through System Settings. If macOS blocks execution, follow Option 1 or Option 2 below.

---

#### Option 1: Manual Authorization via System Settings

If macOS prevents the file from opening:

1. **Right-click** (or Control-click) the application (Installer or Launcher) and select **Open**.
2. If macOS displays a dialog preventing execution, navigate to **System Settings** -> **Privacy & Security**.
3. Scroll down to the **Security** section.
4. Click **Open Anyway** (or **Allow**) next to the blocked application.
5. Confirm by clicking **Open** (using your administrator password or Touch ID).

Repeat these steps for both the **Installer** and the **Launcher** if prompted by macOS.

#### Option 2: Permanent Solution via Terminal (Clearing Quarantine Attributes)

This removes Gatekeeper warnings for HERA by clearing macOS quarantine attributes for all extracted files at once:

1. Open the **Terminal** app (via Spotlight search or `Applications` -> `Utilities`).
2. Run the following command (substituting the path to your extracted folder):

   ```bash
   sudo xattr -cr /path/to/extracted/HERA_folder
   ```

3. Press **Enter**, enter your administrator password when prompted (input characters are not displayed on screen), and confirm.
4. Both the Installer and Launcher can now be opened directly without Gatekeeper prompts.

> [!TIP]
> You can type `sudo xattr -cr` with a trailing space, then drag and drop the extracted HERA folder from Finder into the Terminal window to automatically insert its path.

---

> [!WARNING]
> **Disclaimer & Security Warning:**  
> *This software and the associated instructions are provided "as is", without warranty of any kind. The author expressly disclaims all liability for any damages, errors, or security risks arising from the use of this software or the execution of system commands. Modifying file attributes skips Apple's security checks; proceed at your own risk.*

---

## 2. Running HERA (End Users)

### Standard Interactive Mode

For standard use with the interactive, guided command-line interface:

* **macOS**: **Double-click `HERA_Launcher.command`** in Finder.
  * A Terminal window will open automatically.
  * The launcher auto-detects your MATLAB Runtime installation, sets all internal dynamic library paths, and starts HERA.
  * The guided CLI prompt (`HERA.start_ranking()`) will guide you step-by-step through data loading, metric selection, and statistical parameterization.
* **Windows**: **Double-click `HERA_Runtime.exe`**.
  * A Command Prompt window will open and launch the interactive application.
* **Linux**: Open a terminal and run the execution script with the path to your MATLAB Runtime:

  ```bash
  ./run_HERA_Runtime.sh /usr/local/MATLAB/MATLAB_Runtime/R2025b
  ```

> [!TIP]
> On macOS, keep `HERA_Launcher.command` in the folder where it was extracted. 
> If you installed the MATLAB Runtime in a non-standard location and the launcher cannot locate it automatically, it will prompt you to simply drag and drop your MATLAB Runtime folder into the terminal window.

---

### Command-Line & Terminal Usage (CLI / Non-Interactive)

For automated pipelines, scripted execution, or headless terminal sessions, HERA can be run non-interactively using a JSON configuration file. On both macOS and Linux, this is invoked via the runtime execution script (`run_HERA_Runtime.sh`) with the path to the MATLAB Runtime:

> [!TIP]
> On macOS, double-clicking `HERA_Launcher.command` in Finder is the recommended way to start the interactive CLI. 
> When running from an existing terminal window, automated script, or pipeline, execute `run_HERA_Runtime.sh` directly with your MATLAB Runtime path and `configFile`.

---

#### 1. Analysis with JSON Configuration (Single Run, Automated Pipeline, or Batch)

Run a full ranking analysis using a JSON configuration file, bypassing all interactive UI prompts (ideal for single runs from the terminal, automated scripts, or batch processing):

* **macOS**:

  ```bash
  /Applications/HERA_Runtime/application/run_HERA_Runtime.sh /Applications/MATLAB/MATLAB_Runtime/R2025b configFile "/absolute/path/to/config.json"
  ```

* **Linux**:

  ```bash
  ./run_HERA_Runtime.sh /usr/local/MATLAB/MATLAB_Runtime/R2025b configFile "/absolute/path/to/config.json"
  ```

* **Windows**:

  ```cmd
  HERA_Runtime.exe configFile "C:\absolute\path\to\config.json"
  ```

For configuration parameter specifications and templates, see [Configuration & Parameters](Configuration_&_Parameters.md).

#### 2. Run Verification Unit Tests

Execute the comprehensive 46-test validation suite to verify algorithmic integrity on your system:

* **macOS**: `/Applications/HERA_Runtime/application/run_HERA_Runtime.sh /Applications/MATLAB/MATLAB_Runtime/R2025b runtest true`
* **Linux**: `./run_HERA_Runtime.sh <Path_to_Runtime> runtest true`
* **Windows**: `HERA_Runtime.exe runtest true`

#### 3. Run Convergence Analysis

Perform the robust convergence verification study:

* **macOS**: `/Applications/HERA_Runtime/application/run_HERA_Runtime.sh /Applications/MATLAB/MATLAB_Runtime/R2025b convergence true`
* **Linux**: `./run_HERA_Runtime.sh <Path_to_Runtime> convergence true`
* **Windows**: `HERA_Runtime.exe convergence true`

For further details, refer to the [Convergence Analysis Documentation](Convergence_Analysis.md).

---

## 3. Troubleshooting & FAQ (End Users)

### Q1: What causes dynamic library errors (`dyld: Library not loaded`) on macOS?

* **Explanation**: On macOS, compiled MATLAB applications require dynamic library references to be linked against the
MATLAB Runtime directory. If the execution script `run_HERA_Runtime.sh` is invoked directly in a terminal without passing the MATLAB Runtime path as its first argument, the dynamic linker cannot find the required libraries.
* **Solution**:
  * For **interactive mode**, double-click `HERA_Launcher.command` in Finder (it automatically discovers your MATLAB Runtime installation and supplies all paths).
  * For **command-line or script usage**, always provide the runtime directory path as the first argument:

    ```bash
    /Applications/HERA_Runtime/application/run_HERA_Runtime.sh /Applications/MATLAB/MATLAB_Runtime/R2025b
    ```

### Q2: How should I run HERA in terminal scripts or non-interactive environments?

* **Explanation**: The interactive prompt (`HERA.start_ranking()`) is designed for an active terminal session with live user input. In automated scripts, remote shells, or non-interactive pipelines where standard input cannot be assigned interactively, HERA should be run via configuration files.
* **Solution**:
  * **Interactive CLI**: Double-click `HERA_Launcher.command` in Finder (macOS) or `HERA_Runtime.exe` (Windows) to run the guided step-by-step prompt in a dedicated terminal window.
  * **Terminal & Automated Pipelines**: Pass your analysis parameters via a JSON configuration file using the `configFile` parameter. This executes headlessly without prompting and writes all results directly to disk:
    * **macOS**: `/Applications/HERA_Runtime/application/run_HERA_Runtime.sh <Path_to_Runtime> configFile "/path/to/config.json"`
    * **Linux**: `./run_HERA_Runtime.sh <Path_to_Runtime> configFile "/path/to/config.json"`
    * **Windows**: `HERA_Runtime.exe configFile "C:\path\to\config.json"`

### Q3: Does HERA require an active internet connection or license?

* **Installation**: An internet connection is only needed during initial setup if the installer downloads the free MATLAB Runtime.
* **Execution**: Once installed, HERA operates completely offline. No internet connection, user account, or license check is ever required.

### Q4: What should I do if macOS reports that the app is damaged or blocks execution?

* **Explanation**: Because HERA is open-source academic software distributed without a commercial Apple Developer ID certificate, macOS Gatekeeper prevents execution by default.
* **Solution**: Follow the quick steps in [Step 3: macOS Security Authorization](#step-3-macos-security-authorization-gatekeeper) to approve the application in System Settings (Option 1) or remove quarantine flags via Terminal (Option 2).

---

## 4. Building the Installer from Source (For Developers)

> [!IMPORTANT]
> The instructions below are **only for developers and maintainers** compiling HERA from source. End users do **not** need to build the application.

### Developer Requirements

* **MATLAB** (R2020a or later, R2025b recommended)
* **MATLAB Compiler** toolbox (`compiler.build` and `compiler.package`)
* **Statistics and Machine Learning Toolbox**
* **Parallel Computing Toolbox**

> [!NOTE]
> **Automated Builds:** As documented in [Automated Build (GitHub Actions)](Automated_Build.md), automated cloud builds via GitHub Actions cannot be run by default because the repository maintainer does not possess a cloud MATLAB license for GitHub runners. Standalone installers are therefore built locally in MATLAB using the script below.

### Build Procedure

1. Open MATLAB and navigate to the project root directory.
2. Run the build script:

   ```matlab
   cd deploy
   build_HERA_matlab
   ```

3. The script executes the following automated pipeline:
   * Compiles `+HERA/start_ranking.m` into a standalone application (`compiler.build.standaloneApplication`).
   * Packages the executable into a web installer (`compiler.package.installer` with `'RuntimeDelivery', 'web'`).
   * On macOS, bundles `HERA_Launcher.command` and `readme.txt`, and applies ad-hoc code signatures (`codesign --force --deep -s -`) to minimize Gatekeeper friction.
   * Compresses the output into a distribution ZIP archive in `deploy/output/matlab/`.
