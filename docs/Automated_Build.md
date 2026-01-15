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
resources. It is currently configured **not** to trigger automatically on new tags because I, the repository owner, do not have a MATLAB license for the GitHub runner.

> **Note for Forks:** If you have a valid MATLAB license configured, you can re-enable automatic builds by uncommenting the `release` trigger in `.github/workflows/build_release.yml`.

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

## Publishing to PyPI

To publish the Python package to PyPI, follow this **manual workflow** (since the GitHub Runner cannot build the package due to licensing):

### Prerequisites

1. **PyPI Trusted Publishing**: Configure GitHub Actions as a trusted publisher in your PyPI account settings.
   - Owner: `lerdmann1601`
   - Repository: `HERA-Matlab`
   - Workflow name: `publish_pypi.yml`
   - Environment: `pypi`
2. **Local Environment**: Ensure you have MATLAB and `python`, `pip` installed locally.

### Release Steps

1. **Prepare Artifacts**:
   Run the helper script locally to build the package and generate the distribution files:

   ```bash
   ./deploy/build_and_prep_pypi.sh
   ```

   This will create a `deploy/dist` folder containing `.whl` and `.tar.gz` files.

2. **Create Release**:
   - Go to GitHub -> Releases -> **Draft a new release**.
   - Choose a text tag (e.g., `v1.2.0`).
   - **Upload** the files from `deploy/dist/` to the release.
   - Publish the release.

3. **Publish**:
   - Go to the **Actions** tab in GitHub.
   - Select **Publish to PyPI**.
   - Click **Run workflow**.
   - The action will download the files from your release and upload them to PyPI.
