# Python Integration

HERA provides a compiled Python interface that allows seamless integration into
Python-based data science pipelines. This package wraps the MATLAB functions and
provides them as native Python objects.

> **Note:** The package utilizes the **MATLAB Runtime**. It must be installed separately
> as described below.

## 1. Installation (For End Users)

The easiest way to install the package is via `pip` from PyPI.

### Step 1: Install Package

```bash
pip3 install hera-matlab
```

### Step 2: Install MATLAB Runtime

HERA requires the **MATLAB Runtime R2025b (v25.2)**.
After installing the package, run the following command to check if you have the correct runtime installed or to get the direct download link:

```bash
python3 -m hera_matlab.install_runtime
```

Follow the instructions provided by this command to download and install the runtime if it is missing.

### Step 3: macOS Specifics (Critical)

On macOS, you cannot use the standard `python` interpreter to import the package, as it will likely fail with a library loading error. Instead, you must use the `mwpython` wrapper provided by the MATLAB Runtime.

**Location of mwpython:**
Typically found at:
`/Applications/MATLAB/MATLAB_Runtime/R2025b/bin/mwpython`

**Usage:**
Instead of running `python3 script.py`, run:

```bash
/Applications/MATLAB/MATLAB_Runtime/R2025b/bin/mwpython script.py
```

You may want to add this to your PATH variable for easier access.

> [!WARNING]
> **Python Version Compatibility:**
> HERA (via MATLAB Runtime) currently supports **Python 3.9, 3.10, 3.11, and 3.12**.
> **Python 3.13 is NOT supported.**
>
> If you have Python 3.13 installed, `mwpython` will fail to load your packages.
> **Solution:** Create and activate a Virtual Environment with Python 3.12 before running `mwpython`.
>
> ```bash
> # Example using brew to install python 3.12
> brew install python@3.12
> /opt/homebrew/bin/python3.12 -m venv .venv_hera
> source .venv_hera/bin/activate
> pip3 install hera-matlab
> # Now mwpython will automatically use this environment
> mwpython script.py
> ```

## 2. Usage Modes

### A. Standard Pipeline (File-Based)

This mode replicates the MATLAB batch processing workflow.
It runs the complete analysis based on a JSON configuration file and
automatically generates all PDF reports and plots on disk.

> **Note:** The interactive command-line interface (CLI) is **not supported**
in the Python package. You must use a configuration file.

```python
import hera_matlab

# Initialize Runtime
hera = hera_matlab.initialize()

# Run with JSON configuration
# Outputs (PDFs, Images) will be saved to the 'output_dir' defined in the config
hera.start_ranking('configFile', 'analysis_config.json', nargout=0)

hera.terminate()
```

### B. Direct Data Integration (NumPy/Pandas)

This mode allows you to use HERA as a computational engine within your Python
scripts (e.g., Jupyter Notebooks). You can pass data directly from NumPy/Pandas
and receive the ranking results as a Python dictionary, enabling seamless
integration into larger data science pipelines.

```python
import hera_matlab
import matlab

# Initialize
hera = hera_matlab.initialize()

# Prepare Data (Convert NumPy arrays to MATLAB types)
# Example: 2 Subjects x 2 Methods
data_m1 = matlab.double([[0.1, 0.5], [0.2, 0.4]])
data_m2 = matlab.double([[1.0, 3.0], [1.2, 2.9]])

# Configure Analysis
config = {
    'custom_data': [data_m1, data_m2],
    'metric_names': ['Runtime', 'Accuracy'],
    'dataset_names': ['Method A', 'Method B'],
    'ranking_mode': 'M1_M2',
    'output_dir': 'my_hera_results' # Optional: Specify output folder
}

# Execute Ranking and retrieve Dictionary
results = hera.run_ranking(config, nargout=1)

# Access Results
print(f"Final Ranks: {results['final_rank']}")
print(f"Effect Sizes (Cliff's Delta): {results['d_vals_all']}")
print(f"Effect Sizes (Rel Diff): {results['rel_vals_all']}")
print(f"P-Values: {results['p_vals_all']}")

hera.terminate()
```

See [Results Structure Reference](https://github.com/lerdmann1601/HERA-Matlab/blob/main/docs/Results_Structure_Reference.md) for a complete list of available fields

## 3. Build Instructions (For Maintainers)

To generate the installer and Python package from source (requires MATLAB Compiler SDK):

1. **Run the Build Helper:**

   ```bash
   ./deploy/build_and_prep_pypi.sh
   ```

   This script compiles the MATLAB code, injects the runtime checks, and prepares
   the distribution artifacts (`.whl`, `.tar.gz`) in `deploy/dist`.

2. **Distribution:**
   Upload the generated artifacts from `deploy/dist` to PyPI or GitHub Releases.

## 4. Running the Test Suite

You can execute the internal HERA verification and test suite (Scientific, System, and Unit tests) directly from Python to ensure the installation is valid and mathematically correct.

```python
import hera_matlab

# Initialize
hera = hera_matlab.initialize()

# Launch Test Mode
# Arguments: 'runtest', 'true'
# This will execute all tests and print the results to standard output/log files.
hera.start_ranking('runtest', 'true', nargout=0)

hera.terminate()
```

This is useful for verifying deployments on new machines (e.g., CI/CD).

## 5. Running Convergence Analysis

You can also trigger the robust convergence analysis directly from Python. This is identical to the MATLAB `HERA.start_ranking('convergence', 'true')` command. For more details on the analysis and its parameters, see [Convergence Analysis](Convergence_Analysis.md).

```python
import hera_matlab

# Initialize
hera = hera_matlab.initialize()

# Run Convergence Analysis
# Arguments: 'convergence', 'true'
# Optional: 'sims', 50 (for higher precision)
# Optional: 'logPath', '/path/to/log' (or 'interactive')
hera.start_ranking('convergence', 'true', nargout=0)

hera.terminate()
```
