# Python Integration

HERA provides a compiled Python interface that allows seamless integration into
Python-based data science pipelines. This package wraps the MATLAB functions and
provides them as native Python objects.

> **Note:** The package utilizes the **MATLAB Runtime**. It must be installed separately
> as described below.

## 1. Installation (For End Users)

The easiest way to install the package is via `pip` after installing the
MATLAB Runtime.

### Step 1: Install MATLAB Runtime

Download and install the **MATLAB Runtime (R2025b)** for your operating system
from the [MathWorks Website](https://www.mathworks.com/products/compiler/matlab-runtime.html).

### Step 2: Install Package

Navigate to the `hera_matlab` folder (provided in the release or build output)
and run:

```bash
pip install .
```

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

To generate the installer and Python package from source
(requires MATLAB Compiler SDK):

```matlab
cd deploy
build_HERA_python
```

The output will be generated in `deploy/output/python`.
