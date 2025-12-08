## HERA v1.0.1

HERA is a robust MATLAB toolbox designed for the scientific benchmarking of clustered data. It offers hierarchical ranking logic, advanced effect-size statistics, and automated reporting capabilities.

### Platform Support

* **macOS (Apple Silicon)**

### Installation

This release includes a standalone installer that manages all necessary dependencies. **No MATLAB license is required.**

1. Download and extract the provided ZIP file.
2. Run the **Installer** (`HERA_Runtime_Installer_...`) to set up the application and automatically download the MATLAB Runtime.
3. Start the application using the **HERA_Launcher**.

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
> 6. Now you can double-click the Installer and Launcher without issues.
>
> *Disclaimer: This software and the associated instructions are provided "as is", without warranty of any kind. The authors expressly disclaim all liability for any damages, errors, or security risks arising from the use of this software or the execution of system commands. Modifying file attributes skips Apple's security checks; proceed at your own risk.*

### Key Features

* **Hierarchical Ranking:** Facilitates the comparison of methods using multiple weighted metrics.
* **Statistical Rigor:** Includes built-in support for Wilcoxon signed-rank tests, Cliff's Delta, Relative Difference, and Bootstrapping.
* **Automated Reporting:** Generates comprehensive PDF reports and visualizations, including Sankey diagrams and Win-Loss matrices.
* **Reproducibility:** Ensures consistent results through configuration-based workflows with fixed seeds.

**Full Changelog**: <https://github.com/lerdmann1601/HERA-Matlab/commits/v1.0.1>
