# Automated Build (GitHub Actions)

> **Note:** The automated build workflow requires a valid MATLAB license to be configured
> as a secret(`MATLAB_LICENSE`) in the repository settings.
> I sadly can not provide a license for this respository.
> Student licenses may not support this feature.

To enable GitHub Actions for building and testing, you need to provide a valid
MATLAB license:

1. Go to your repository's **Settings** > **Secrets and variables** > **Actions**.
2. Click **New repository secret**.
3. Name the secret `MATLAB_LICENSE` and paste the contents of your license file.

## Running the Build

This workflow is set to **manual execution** (`workflow_dispatch`) to save
resources but will be triggered once if a new version tag is pushed to the repository.

1. Navigate to the **Actions** tab in the repository.
2. Select **Build HERA Runtime** from the sidebar.
3. Click the **Run workflow** button.

The workflow performs the following steps:

1. **Unit Testing**: Runs the full test suite (`HERA.run_unit_test`) to ensure
   code integrity.
2. **Compilation**: Builds the standalone application for the target operating
   system (macOS/Linux/Windows) using the `deploy/build_HERA_matlab.m` script.
3. **Toolbox Packaging**: Packages the code as a MATLAB Toolbox (`.mltbx`) using
`package_HERA_toolbox.m`.
4. **Python Build**: Compiles the Python interface using `build_HERA_python.m`.
5. **Artifact Upload**: Uploads the compiled installer, toolbox, and Python
package as build artifacts, which can be downloaded from the GitHub Actions
run page.
