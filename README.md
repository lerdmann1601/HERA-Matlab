<!-- markdownlint-disable MD013 MD041 -->
<!-- markdownlint-disable MD033 -->
<div align="center">

<img src="assets/hera_logo.svg" alt="HERA Logo" width="300"/>

# HERA: Hierarchical-Compensatory, Effect-Size-Driven and Non-Parametric Ranking Algorithm

[![MATLAB](https://img.shields.io/badge/MATLAB-R2024a%2B-orange.svg)](https://www.mathworks.com/products/matlab.html)
[![Statistics Toolbox](https://img.shields.io/badge/Toolbox-Statistics_and_Machine_Learning-blue.svg)](https://www.mathworks.com/products/statistics.html)
[![Parallel Computing Toolbox](https://img.shields.io/badge/Toolbox-Parallel_Computing-blue.svg)](https://www.mathworks.com/products/parallel-computing.html)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![View on GitHub](https://img.shields.io/badge/GitHub-View_on_GitHub-blue?logo=github)](https://github.com/lerdmann1601/HERA-Matlab)
[![Issues](https://img.shields.io/github/issues/lerdmann1601/HERA-Matlab)](https://github.com/lerdmann1601/HERA-Matlab/issues)
<!-- markdownlint-disable-next-line MD013 -->
[![ORCID](https://img.shields.io/badge/ORCID-0009--0009--3758--7363-green.svg)](https://orcid.org/0009-0009-3758-7363)

<!-- markdownlint-disable-next-line MD036 -->
**Made for Scientific Benchmarking**

[Key Features](#key-features) • [Installation](#installation) • [Quick Start](#quick-start) • [Documentation](#documentation) • [Citation](#citation)

</div>
<!-- markdownlint-enable MD033 MD013 -->

---

## Overview

**HERA** is a MATLAB toolbox designed to automate the objective comparison of
algorithms, experimental conditions, or datasets across multiple quality
metrics. Unlike traditional ranking methods that rely solely on mean values or
p-values, HERA employs a **hierarchical-compensatory logic** that integrates:

* **Significance Testing**: Wilcoxon signed-rank tests for paired data.
* **Effect Sizes**: Cliff's Delta and Relative Mean Difference for practical relevance.
* **Bootstrapping**: Data-driven thresholds and BCa confidence intervals.

This ensures that a "win" is only counted if it is both **statistically
significant** and **practically relevant**, providing a robust and nuanced
ranking system.

---

## Key Features

* **Hierarchical Logic**: Define primary and secondary metrics. Secondary
  metrics can act as tie-breakers or rank correctors (e.g., `M1_M2`,
  `M1_M2_M3`).
* **Data-Driven Thresholds**: Automatically calculates adaptive effect size
  thresholds using Percentile Bootstrapping.
* **Robustness**: Utilizes Bias-Corrected and Accelerated (BCa) confidence
  intervals and Cluster Bootstrapping for rank stability.
* **Automated Reporting**: Generates PDF reports, Win-Loss Matrices, Sankey
  Diagrams, and machine-readable JSON/CSV exports.
* **Reproducibility**: Supports fixed-seed execution and configuration
  file-based workflows.

---

## Installation

### Requirements

* **MATLAB** (R2024a or later Required)
* **Statistics and Machine Learning Toolbox** (Required)
* **Parallel Computing Toolbox** (Required for performance)

### Setup

#### Option A: MATLAB Toolbox (Recommended)

1. Download the latest `HERA_vX.Y.Z.mltbx` from the
   [Releases](https://github.com/lerdmann1601/HERA-Matlab/releases) page.
2. Double-click the file to install it.
3. Done! HERA is now available as a command (`HERA.start_ranking`) in MATLAB.

#### Option B: Git Clone (for Developers)

1. **Clone the repository:**

    ```bash
    git clone https://github.com/lerdmann1601/HERA-Matlab.git
    ```

2. **Install/Configure Path:**

    Navigate to the repository folder and run the setup script to add HERA to
    your MATLAB path.

    ```matlab
    cd HERA-Matlab
    setup_HERA
    ```

<!-- markdownlint-disable MD033 -->
<details>
<summary><strong>Standalone Runtime</strong></summary>

HERA can be compiled into a standalone application for macOS, Linux, and
Windows. The build process generates an **installer** that automatically
downloads the required MATLAB Runtime, making it easy to distribute.

> **Download:** A pre-built installer for macOS (Apple Silicon) is available as
> a ZIP archive in the [Releases](https://github.com/lerdmann1601/HERA-Matlab/releases)
> section.

#### Building the Installer

To build the installer, you need a MATLAB installation with the **MATLAB
Compiler** toolbox.

1. Open MATLAB and navigate to the project root.
2. Run the build script:

   ```matlab
   cd deploy
   build_hera
   ```

3. The artifacts (Installer + ZIP) will be generated in `deploy/output`.

#### Installation and Usage

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

</details>
<!-- markdownlint-enable MD033 -->

<!-- markdownlint-disable MD033 -->
<details>
<summary><strong>Automated Build (GitHub Actions)</strong></summary>

> **Note:** The automated build workflow requires a valid MATLAB license to be configured
> as a secret(`MATLAB_LICENSE`) in the repository settings.
> I sadly can not provide a license for this respository.
> Student licenses may not support this feature.

To enable GitHub Actions for building and testing, you need to provide a valid
MATLAB license:

1. Go to your repository's **Settings** > **Secrets and variables** > **Actions**.
2. Click **New repository secret**.
3. Name the secret `MATLAB_LICENSE` and paste the contents of your license file.

**Running the Build:**

This workflow is set to **manual execution** (`workflow_dispatch`) to save
resources.

1. Navigate to the **Actions** tab in the repository.
2. Select **Build HERA Runtime** from the sidebar.
3. Click the **Run workflow** button.

The workflow performs the following steps:

1. **Unit Testing**: Runs the full test suite (`HERA.run_unit_test`) to ensure
   code integrity.
2. **Compilation**: Builds the standalone application for the target operating
   system (macOS/Linux/Windows) using the `deploy/build_hera.m` script.
3. **Artifact Upload**: Uploads the compiled installer and application as build
   artifacts, which can be downloaded from the GitHub Actions run page.

</details>
<!-- markdownlint-enable MD033 -->

<!-- markdownlint-disable MD033 -->
</details>
<!-- markdownlint-enable MD033 -->

<!-- markdownlint-disable MD033 -->
<details>
<summary><strong>Python Integration</strong></summary>

HERA provides a compiled Python interface that allows seamless integration into
Python-based data science pipelines. This package wraps the MATLAB functions and
provides them as native Python objects.

> **Note:** The package utilizes the **MATLAB Runtime**. The provided installer
> automates the setup of this dependency.

#### 1. Installation (For End Users)

The easiest way to install the package is via `pip` after installing the
MATLAB Runtime.

#### Step 1: Install MATLAB Runtime

Download and install the **MATLAB Runtime (R2025b)** for your operating system
from the [MathWorks Website](https://www.mathworks.com/products/compiler/matlab-runtime.html).

#### Step 2: Install Package

Navigate to the `hera_matlab` folder (provided in the release or build output)
and run:

```bash
pip install .
```

> **Note:** If an automated installer (`hera_matlab_Installer`) is provided for
your specific OS, you can run it to handle both steps automatically.

#### 2. Usage Modes

#### A. Standard Pipeline (File-Based)

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

#### B. Direct Data Integration (NumPy/Pandas)

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
# See "Results Structure Reference" below for a complete list of available fields
print(f"Final Ranks: {results['final_rank']}")
print(f"Effect Sizes (Cliff's Delta): {results['d_vals_all']}")
print(f"Effect Sizes (Rel Diff): {results['rel_vals_all']}")
print(f"P-Values: {results['p_vals_all']}")

hera.terminate()
```

#### 3. Build Instructions (For Maintainers)

To generate the installer and Python package from source (requires MATLAB Compiler SDK):

```matlab
cd deploy
build_hera_python
```

The output will be generated in `deploy/output/python`.

</details>
<!-- markdownlint-enable MD033 -->

---

## Quick Start

### 1. Interactive Mode (Recommended for Beginners)

The interactive command-line interface guides you through every step of the configuration,
from data selection to statistical parameters.
If you are new to HERA, this is the recommended mode.
At any point, you can exit the interface by typing `exit` or `quit` or `q`.

```matlab
HERA.start_ranking()
```

### 2. Batch Mode (Reproducible / Server)

For automated analysis or reproducible research, use a JSON configuration file.

```matlab
HERA.start_ranking('configFile', 'config.json')
```

### 3. Unit Test Mode

Run the built-in validation suite to ensure HERA is working correctly on your system.

```matlab
% Run tests and save log to default location
HERA.start_ranking('runtest', 'true')

% Run tests and save log to a specific folder
HERA.start_ranking('runtest', 'true', 'logPath', '/path/to/logs')
```

> **Note:** Example use cases with synthetic datasets and results are
> provided in the `data/examples` directory. See `data/README.md` for a
> walkthrough of the example use cases and visual examples of the ranking
> outputs.
>
> **Note:** HERA is designed for high-performance scientific computing, featuring
> **fully parallelized bootstrap procedures** and **automatic memory management**
> to optimize efficiency. However, specifically due to the extensive use of
> bootstrapping, it remains a **CPU-intensive application**. Please ensure you
> have access to enough CPU cores for reasonable performance.
---

## Documentation

<!-- markdownlint-disable MD033 -->
<details>
<summary><strong>Repository Structure</strong></summary>

The codebase is organized as a MATLAB package (`+HERA`) to ensure namespace isolation.

```text
HERA-Matlab/
├── +HERA/                     % Main Package Namespace
│   ├── +output/               % Report Generation (PDF, JSON, CSV)
│   ├── +plot/                 % Visualization (Sankey, Heatmaps)
│   ├── +run/                  % Execution Logic
│   ├── +start/                % CLI & Configuration Logic
│   ├── +stats/                % Statistical Core (Cliff's Delta, Convergence Check)
│   ├── +test/                 % Unit Test Suite
│   ├── language/              % Localization Files
│   ├── bootstrap_ranking.m    % Cluster Bootstrap Analysis
│   ├── borda_ranking.m        % Consensus Ranking
│   ├── calculate_bca_ci.m     % BCa Confidence Intervals
│   ├── calculate_ranking.m    % Core Ranking Logic
│   ├── calculate_thresholds.m % Threshold Calculation
│   ├── default.m              % Global Defaults
│   ├── design.m               % Style & Design Definitions
│   ├── generate_output.m      % Output Generation Controller
│   ├── generate_plots.m       % Plot Generation Controller
│   ├── get_language.m         % Language Loader
│   ├── get_version.m          % Version Retrieval
│   ├── load_data.m            % Data Import & Validation
│   ├── power_analysis.m       % Power Analysis
│   ├── run_ranking.m          % Core Function (Developer API)
│   ├── start_ranking.m        % Main Entry Point (User API)
│   └── run_unit_test.m        % Test Runner
├── assets/                    % Images & Logos
├── data/                      % Data Directory
│   ├── examples/              % Synthetic Example Datasets
│   └── README.md              % Data Documentation
├── deploy/                    % Build Scripts (Standalone App)
├── paper/                     % Paper Resources
├── tests/                     % Unit Test & Analysis Reports
├── setup_HERA.m               % Path Setup Script
├── CITATION.cff               % Citation Metadata
├── CODE_OF_CONDUCT.md         % Community Standards
├── CONTRIBUTING.md            % Contribution Guidelines
├── LICENSE                    % License File
└── README.md                  % Global Documentation
```

</details>
<!-- markdownlint-enable MD033 -->

<!-- markdownlint-disable MD033 MD013 MD060 -->
<details>
<summary><strong>Ranking Modes Explained</strong></summary>

| Mode | Behavior | Use Case |
| :--- | :--- | :--- |
| `M1` | Ranks strictly by Metric 1. | Single-metric evaluation. |
| `M1_M2` | Metric 1 is primary. Metric 2 can correct (swap) ranks if a lower-ranked method significantly outperforms a higher-ranked one in Metric 2. | Balancing performance vs. cost. |
| `M1_M3A` | Metric 1 is primary. Metric 2 acts strictly as a tie-breaker. | Tie-breaking without overriding primary results. |
| `M1_M2_M3` | Full hierarchy. M1 is primary. M2 corrects M1 (iterative). M3 applies two sub-logics: (1) **One-time correction** if M2 is neutral, and (2) **Iterative tie-breaking** if both M1 and M2 are neutral. | Complex multi-objective benchmarking. |

</details>
<!-- markdownlint-enable MD033 MD013 MD060 -->

<!-- markdownlint-disable MD033 -->
<details>
<summary><strong>Input Data Specification</strong></summary>

* **Format**: CSV or Excel (`.xlsx`).
* **Organization**: One file per metric.
* **Filename**: The filename (excluding extension) must strictly match the
    corresponding entry in `metric_names` (e.g., `Accuracy.csv` for metric
    `Accuracy`).
* **Dimensions**: Rows = Observations (*n*, e.g., Subjects), Columns =
    Datasets (*N*, e.g., Methods).
* **Consistency**: All files must have identical dimensions. Uneven sample sizes
    (missing data) are handled by automatic `NaN` padding (empty cells) to
    ensure a uniform matrix size, and pairwise deletion is applied during
    analysis.

</details>
<!-- markdownlint-enable MD033 -->

<!-- markdownlint-disable MD033 -->
<details>
<summary><strong>Methodological Guidelines & Limitations</strong></summary>

The statistical rigor of HERA (e.g., Holm-Bonferroni correction, Bootstrapping)
imposes practical limits on the number of datasets (*N*) and sample size (*n*).
Therefore the following guidelines are provided as theoretical considerations
but should not be taken as strict requirements.

**Number of Datasets (*N*):**
Increasing *N* quadratically increases the number of pairwise comparisons (*m* =
*N*(*N*-1)/2), which reduces statistical power due to strict corrections.

* **Minimum (*N* = 3)**: Required for a meaningful ranking. (*N* = 2 is just a
    simple comparison).
* **Optimal (*N* ≈ 8–10)**: Balances ranking depth with statistical
    power (28–45 comparisons).
* **Upper Limit (*N* ≈ 15)**: Not generally recommended. The loss of power
    from FWER corrections makes detecting true differences unlikely. However,
    it is possible to use HERA with *N* > 15 and you can just give it a try.

> **Visual Limit (*N* ≤ 20)**: While HERA should technically compute rankings
> for any *N* (exported to CSV/JSON), the generated plots (e.g. Win-Loss Matrix,
> Final Summary) visually degrade beyond *N* = 20. For *N* > 20, I recommend
> relying on the machine-readable and text-based outputs. You can disable plots
> (`create_reports: false`) to save runtime.
>
> **Recommendation:**
> If you have a large pool of candidates (*N* >> 15), it could be a good idea to
> apply a global screening method (e.g., **Friedman Test** followed by Nemenyi post-hoc)
> to identify the top tier of algorithms. Ranking the entire set with HERA may be
> overly strict; instead, select the top performing subset (e.g., the best 10-15)
> and use HERA for the final ranking of the best candidates.

**Sample Size (*n*):**
A balance between statistical stability and computational feasibility is
required.

* **Minimum (*n* = 16)**: Required for the Wilcoxon test to use the Normal
    Approximation in Matlab.
* **Robust Min (*n* ≈ 25–30)**: Necessary for stable BCa confidence
    intervals and Jackknife estimates (Although it automatically switches
    to Percentil Bootstrap if Bias or Jackknife estimates become unstable).
* **Optimal (*n* ≈ 50–300)**: Best balance of power, stability, and
    runtime.
* **Upper Limit (*n* ≈ 1,000–5,000)**: Higher *n* improves statistics
    but linearly scales runtime. *n* ≫ 5,000 may be computationally
    impractical due to extensive bootstrapping.

> **Recommendation:** Perform an *a priori* power analysis to estimate the
> required *n* for your chosen *N*.

**Missing Data Handling (NaN):**
HERA is robust against missing data (`NaN`) but handling it comes with trade-offs:

* **Pairwise Deletion**: HERA employs pairwise deletion to maximize data
    usage without requiring complex imputation. While this assumes data is missing
    completely at random (MCAR), it remains methodologically robust: By relying
    on discrete, independent pairwise comparisons, the algorithm avoids the
    mathematical inconsistencies (e.g., non-positive definite matrices) that
    typically compromise pairwise deletion in global multivariate statistics.
* **Robust Bootstrapping**: If `NaN`s are detected, HERA automatically switches
    to a "Robust Path". This dynamically filters invalid data for each bootstrap
    resample to ensure correctness, which **significantly increases runtime**
    especially for large sample sizes (*n*).
* **Automatic Warning**: A warning is issued if valid data drops below 80% for
    any comparison however it is not a strict requirement.

> **Recommendation**: Minimize `NaN`s to preserve statistical power and
> performance. For critical analyses with substantial data loss, use
> validated imputation methods (e.g., MICE) *before* running HERA.

**Number of Metrics (*M* ≤ 3):**
HERA is designed for a maximum of 3 hierarchical metrics to maintain methodological
robustness. This limit is inherent to the hierarchical-compensatory design and is
based on the following methodological considerations:

* **Loss of Interpretability**: With every additional hierarchical level, the causal
    chain of the ranking decision becomes opaque and increasingly difficult to trace.
    Limiting the depth to 3 levels ensures that the ranking logic remains transparent
    and empirically justifiable.
* **Increased Risk of Collinearity**: Adding more metrics increases the probability
    of introducing redundant criteria (e.g., two metrics measuring valid features
    of the same underlying property). In a sequential logic, these correlates
    would be falsely treated as independent evidence, distorting the ranking.
* **Functional Saturation**: The hierarchical-compensatory logic is fully saturated
    by three levels (Sorting, Correction, Finalization). Adding a fourth metric
    yields diminishing margins of utility, as the probability of meaningful rank
    adjustments approaches zero, while the complexity of the decision model
    increases disproportionately.

> **Recommendation:**
> If you want to consider more than 3 metrics and use HERA you could first
> perform a check for collinearity (e.g., using a correlation matrix).
> Strongly correlated metrics could be aggregated into a common factor
> (e.g., via Principal Component Analysis (PCA)) before running the HERA
> analysis with up to 3 metrics.
>
> If your study design requires the **simultaneous integration** of a large number
> of metrics ($M \gg 3$) HERA is not feasible and compensatory or outranking MCDA
> methods are methodologically more appropriate. In this case, approaches like
> **TOPSIS** or **PROMETHEE** might be a better choice.

</details>
<!-- markdownlint-enable MD033 -->

<!-- markdownlint-disable MD033 -->
<details>
<summary><strong>Configuration Reference (Full)</strong></summary>

To run HERA in **Batch Mode**, create a `.json` file (e.g.,
`analysis_config.json`). Below is the complete list of available parameters.

**Minimal Configuration:**

```json
{
  "userInput": {
    "folderPath": "/Users/Name/Data",
    "fileType": ".csv",
    "metric_names": ["Accuracy", "Runtime", "Memory"],
    "ranking_mode": "M1_M2_M3",
    "output_dir": "/Users/Name/Results"
  }
}
```

<details>
<summary><strong>Full Configuration Example (JSON)</strong></summary>

This example shows **all** possible parameters with their default values.
Parameters inside `system` and `bootstrap_*` must be nested correctly as shown.

```json
{
  "userInput": {
    "folderPath": "/Path/To/Data",
    "fileType": ".csv",
    "metric_names": ["Metric1", "Metric2", "Metric3"],
    "output_dir": "/Path/To/Results",
    "language": "en",
    "ranking_mode": "M1_M2_M3",
    
    "reproducible": true,
    "seed": 123,
    "num_workers": "auto",
    "create_reports": true,
    "plot_theme": "light",
    
    "ci_level": 0.95,
    "alphas": [0.05, 0.05, 0.05],
    "run_sensitivity_analysis": true,
    "run_power_analysis": true,
    "power_simulations": 10000,
    "min_data_completeness": 0.80,

    "manual_B_thr": 2000,
    "manual_B_ci": 10000,
    "manual_B_rank": 500,

    "system": {
      "target_memory": "auto",
      "jack_parfor_thr": 300,
      "jack_vec_limit": 150,
      "delta_mat_limit": 30000,
      "min_batch_size": 100
    },

    "bootstrap_thresholds": {
      "B_start": 100,
      "B_step": 100,
      "B_end": 10000,
      "n_trials": 25,
      "convergence_tolerance": 0.01,
      "smoothing_window": 3,
      "convergence_streak_needed": 3,
      "min_steps_for_convergence_check": 1
    },

    "bootstrap_ci": {
      "B_start": 100,
      "B_step": 200,
      "B_end": 20000,
      "n_trials": 30,
      "convergence_tolerance": 0.03,
      "smoothing_window": 3,
      "convergence_streak_needed": 3,
      "min_steps_for_convergence_check": 1
    },

    "bootstrap_ranks": {
      "B_start": 50,
      "B_step": 25,
      "B_end": 2500,
      "n_trials": 15,
      "convergence_tolerance": 0.005,
      "smoothing_window": 3,
      "convergence_streak_needed": 3,
      "min_steps_for_convergence_check": 1
    }
  }
}
```

</details>

</details>

<details>
<summary><strong>Parameter Dictionary</strong></summary>

<!-- markdownlint-disable MD013 MD060 -->
| Category | Parameter | Type | Default | Description |
| :--- | :--- | :--- | :--- | :--- |
| **Input/Output** | `folderPath` | string | - | Absolute path to data folder. |
| | `fileType` | string | - | `.csv` or `.xlsx`. |
| | `metric_names` | array | - | List of filenames (metrics) in hierarchical order. |
| | `output_dir` | string | - | Path to save results. |
| | `language` | string | `"en"` | Output language code. |
| **Logic** | `ranking_mode` | string | `"M1_M2_M3"` | Logic mode (`M1`, `M1_M2`, `M1_M3A`, `M1_M2_M3`). |
| | `run_sensitivity_analysis` | bool | `true` | Run ranking for all metric permutations. |
| | `run_power_analysis` | bool | `true` | Run post-hoc power analysis. |
| | `min_data_completeness` | double | `0.80` | Min fraction of valid pairs required (0.8 = 80%). |
| **Statistics** | `ci_level` | double | `0.95` | Confidence interval level (e.g., 0.95 for 95%). |
| | `alphas` | array | `[0.05, ...]` | Significance level for each metric. |
| | `power_simulations` | int | `10000` | Number of simulations for power analysis. |
| **System** | `reproducible` | bool | `true` | Use fixed RNG seed. |
| | `seed` | int | `123` | RNG seed value. |
| | `num_workers` | int/str | `"auto"` | Number of parallel workers. `"auto"` uses `parcluster('local').NumWorkers`. |
| | `system.target_memory` | int | `"auto"` | Target memory per chunk (MB). Automatically calculated based on available RAM, but can be manually defined via JSON config file. |
| | `system.jack_parfor_thr` | int | `300` | Min N to trigger parallel execution. |
| | `system.jack_vec_limit` | int | `150` | Max N for vectorized Jackknife calculations. |
| | `system.delta_mat_limit` | int | `30000` | Max N*M product for matrix-based Cliff's Delta. |
| | `system.min_batch_size` | int | `100` | Min batch size for parallel processing. |
| **Graphics** | `create_reports` | bool | `true` | Generate PDF reports and high-res plots. If `false`, only essential convergence and diagnostics plots are saved. |
| | `plot_theme` | string | `"light"` | `"light"` or `"dark"`. |
| **Bootstrap (Manual)** | `manual_B_thr` | int | `2000` | Iterations for Thresholds (empty = auto). |
| | `manual_B_ci` | int | `5000` | Iterations for CIs (empty = auto). |
| | `manual_B_rank` | int | `500` | Iterations for Rank Stability (empty = auto). |
| **Bootstrap (Auto)** | `bootstrap_thresholds` | struct | (See Below) | Config for Threshold convergence. |
| | `bootstrap_ci` | struct | (See Below) | Config for CI convergence. |
| | `bootstrap_ranks` | struct | (See Below) | Config for Rank Stability convergence. |

</details>

<details>
<summary><strong>Bootstrap Configuration (Auto-Convergence)</strong></summary>

The `bootstrap_*` structs (e.g., `bootstrap_thresholds`) support the following
nested parameters to control the convergence algorithm:

| Field | Type | Default (Thr/CI/Rank) | Description |
| :--- | :--- | :--- | :--- |
| `B_start` | int | `100` / `100` / `50` | Initial number of bootstrap iterations. |
| `B_step` | int | `100` / `200` / `25` | Iterations to add in each step. |
| `B_end` | int | `10000` / `20000` / `2500` | Maximum number of iterations. |
| `n_trials` | int | `25` / `30` / `15` | Number of independent trials per step to check stability. |
| `convergence_tolerance` | double | `0.01` / `0.03` / `0.005` | Max allowed variation (e.g., 0.005 = 0.5%). |
| `smoothing_window` | int | `3` / `3` / `3` | Window size for moving average smoothing. |
| `convergence_streak_needed` | int | `3` / `3` / `3` | Consecutive steps required to pass tolerance. |
| `min_steps_for_convergence_check` | int | `1` | Minimum steps before checking convergence. |
<!-- markdownlint-enable MD013 MD060 -->

<details>
<summary><strong>Convergence Modes</strong></summary>

* **Simple Convergence**: Used when `smoothing_window` is empty. Checks if the
    value changes less than `convergence_tolerance` between steps.
* **Robust Convergence**: Used when `smoothing_window` is set (default). Uses
    a moving average to smooth fluctuations and requires
    `convergence_streak_needed` consecutive stable steps.
* **Fallback Convergence** (Elbow Method): Triggered if the primary
    plateau criterion is not met. It uses a heuristic to detect the point of
    diminishing returns. **Warning:** This is a diagnostic aid, not a guarantee!
    Results should be treated as a suggestion to be verified by visual inspection
    and for further analysis (for more details see Troubleshooting Convergence).

> **Note:** The different parameters used in the **Robust Mode** for
> convergence checking have been validated in a empirical Convergence Analysis.
> The default parameters were found to be sufficient for the majority of scenarios.
> However it is not a guarantee!
> Convergence failures using the elbow method were rare (0.3% of cases,
> exclusively in Ranking). Simple Convergence Mode was not assessed in this analysis!
> For detailed results, see
> [tests/Convergence_Analysis.md](tests/Convergence_Analysis.md).

</details>

<details>
<summary><strong>Troubleshooting Convergence</strong></summary>

While the automated check should work for most datasets, "difficult" data with high
variance or flat likelihood landscapes may fail to converge within
*B*<sub>max</sub>. In this case, you can try the following:

1. **Check the Elbow**: Inspect the generated stability plots. If you see a
    clear "elbow" where the curve flattens but fluctuates slightly above the
    strict tolerance *ε*, the convergence parameters might be too strict
    for your data's noise level.
2. **Adjust Parameters**: You can relax `convergence_tolerance` (e.g., to
    0.02) or increase `n_trials` and/ or `smoothing_window` in the configuration.
3. **Use Simple Convergence**: Select Simple Convergence in the CLI or set
    `smoothing_window` to empty to use simple convergence. Be aware that
    this might not be the most robust option! Choose a higher `min_steps_for_convergence_check`
    e.g. 3 to ensure that the convergence check will not be influenced by high
    initial fluctuations of stability measures.
4. **Manual Override**: If no clear convergence is found (no elbow), or for
    theoretical guarantees of large numbers, you can just
    use fixed high *B* values (e.g., `manual_B_ci = 15000`, `manual_B_thr = 2000`)
    as per literature recommendations.

> **Reproducibility Note:** Visual inspection of convergence plots is strongly
> recommended for final reporting! If the automated check finds convergence too
> early (e.g., in robust mode), you can make the parameters stricter (e.g.,
> decrease `convergence_tolerance` or increase `B_step` and/or `B_start`). All
> procedures should use the same fixed random seed for full reproducibility.

</details>

</details>
<!-- markdownlint-enable MD033 -->

<!-- markdownlint-disable MD033 -->
<details>
<summary><strong>Advanced Usage (Developer Mode)</strong></summary>

Developers can call `HERA.run_ranking` directly with data matrices, bypassing
file I/O. This is useful for integration into other pipelines or simulation
studies.

**Syntax:**

```matlab
results = HERA.run_ranking(userInput);
```

**Input Structure (`userInput`):**
Instead of `folderPath`, provide `custom_data`:

```matlab
% 1. Prepare Data
% Cell array of (N_Subjects x N_Methods) matrices
data_m1 = randn(50, 5); 
data_m2 = randn(50, 5);
custom_data = {data_m1, data_m2};

% 2. Configure User Input
userInput = struct();
userInput.custom_data = custom_data;
userInput.metric_names = {'Accuracy', 'Runtime'}; % Must match data count
userInput.dataset_names = {'Method A', 'Method B', 'Method C', 'Method D',
    'Method E'};
userInput.ranking_mode = 'M1_M2';
userInput.output_dir = pwd;

% 3. Run
results = HERA.run_ranking(userInput);
```

**Output (`results` struct):**

* `.final_rank`: Final ranking vector.
* `.d_vals_all`: Effect sizes (Cliff's Delta).
* `.p_vals_all`: Raw p-values.
* `.ci_d_all`: Confidence intervals.

</details>

<!-- markdownlint-disable MD033 -->
<details>
<summary><strong>Results Structure Reference</strong></summary>

When running `results = HERA.run_ranking(...)`, the returned structure contains:

<!-- markdownlint-disable MD013 MD060 -->
| Field | Dimensions | Description |
|---|---|---|
| **Final Results** | | |
| `final_rank` | `[N x 1]` | Final rank for each dataset (1 = Best). |
| `final_order` | `[1 x N]` | Indices of datasets sorted by rank. |
| `final_bootstrap_ranks` | `[N x B]` | Bootstrapped rank distribution for stability analysis. |
| `ci_lower_rank` | `[N x 1]` | Lower bound of Rank Confidence Interval. |
| `ci_upper_rank` | `[N x 1]` | Upper bound of Rank Confidence Interval. |
| `thresholds` | struct | Thresholds for Cliff's Delta and RelDiff. |
| **Statistics** | | |
| `d_vals_all` | `[Pairs x M]` | Cliff's Delta effect sizes for all pairs/metrics. |
| `rel_vals_all` | `[Pairs x M]` | Relative Mean Differences. |
| `ci_d_all` | `[Pairs x 2 x M]` | BCa Confidence Intervals for Delta. |
| `ci_r_all` | `[Pairs x 2 x M]` | BCa Confidence Intervals for RelDiff. |
| `all_p_value_matrices` | `{1 x M}` | Cell array of raw p-values (Wilcoxon). |
| `all_alpha_matrices` | `{1 x M}` | Cell array of Holm-Bonferroni corrected alphas. |
| `all_sig_matrices` | `{1 x M}` | Logical matrices indicating significant wins. |
| **Diagnostics** | | |
| `swap_details` | struct | Log of logic-based rank swaps (M1 vs M2). |
| `intermediate_orders` | struct | Rankings after each hierarchical stage. |
| `borda_results` | struct | Consensus ranking from sensitivity analysis. |
| `power_results` | struct | Post-hoc power analysis data. |
| `all_permutation_ranks` | `[N x Perms]` | Ranks for every metric permutation tested. |
| `selected_permutations` | `[Perms x M]` | Indices of metrics for each permutation. |
<!-- markdownlint-enable MD013 -->

</details>

---

## Outputs

HERA generates a timestamped directory containing:

<!-- markdownlint-disable MD013 -->
```text
Ranking_<Timestamp>/
├── Output/
│   ├── results_*.csv                 % Final ranking table (Mean ± SD of metrics and rank CI)
│   ├── data_*.json                   % Complete analysis record (Inputs, Config, Stats, Results)
│   ├── log_*.csv                     % Detailed log of pairwise comparisons and logic
│   ├── sensitivity_details_*.csv     % Results of the Borda sensitivity analysis
│   ├── BCa_Correction_Factors_*.csv  % Correction factors (Bias/Skewness) for BCa CIs
│   └── bootstrap_rank_*.csv          % Complete distribution of bootstrapped ranks
├── Graphics/                         % High-res PNGs organized in subfolders
│   ├── Ranking/
│   ├── Detail_Comparison/
│   ├── CI_Histograms/
│   └── Threshold_Analysis/
├── PDF/                              % Specialized reports
│   ├── Ranking_Report.pdf
│   ├── Convergence_Report.pdf
│   └── Bootstrap_Report.pdf
├── Final_Ranking_*.png               % Summary graphic of ranking result
├── Final_Report_*.pdf                % Consolidated graphical report of the main results
├── Ranking_*.txt                     % Complete console log of the session
└── configuration.json                % Reusable configuration file to reproduce the run
```
<!-- markdownlint-enable MD013 -->

---

## Testing

HERA includes a comprehensive validation framework (`run_unit_test.m`)
comprising **46 test cases** organized into four suites:

1. **Unit Tests** (19 cases): Checks individual components, helper functions, and
    execution logic (Run/Start packages) to ensure specific parts of the code work
    correctly.
2. **Statistical Tests** (5 cases): Verifies the core mathematical functions
    (e.g., Jackknife, Cliff's Delta) and ensures the performance optimizations
    (hybrid switching) work as intended.
3. **Scientific Tests** (19 cases): Comprehensive validation of ranking logic,
    statistical accuracy, and robustness against edge cases (e.g., zero
    variance, outliers).
4. **System Tests** (3 cases): Runs the entire HERA pipeline from start to
    finish to ensure that the JSON configuration (batch), Developer API and NaN
    Data handling are working correctly.

### Running Tests

You can run the test suite in three ways:

1. **Auto-Log Mode (Default)**
    Automatically finds a writable folder (e.g., Documents) to save the log
    file.

    ```matlab
    HERA.run_unit_test()
    ```

2. **Interactive Mode**
    Opens a dialog to select where to save the log file.

    ```matlab
    HERA.run_unit_test('interactive')
    ```

3. **Custom Path Mode**
    Saves the log file to a specific directory.

    ```matlab
    HERA.run_unit_test('/path/to/my/logs')
    ```

### GitHub Actions (Cloud Testing)

For reviewers or users without a local MATLAB license, you can run the test
suite directly on GitHub:

1. Go to the **Actions** tab in this repository.
2. Select **Testing HERA** from the left sidebar.
3. Click **Run workflow**.

---

## Contributing

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for details.

1. Fork the repository.
2. Create a feature branch.
3. Commit your changes.
4. Open a Pull Request.

---

## Citation

If you use HERA in your research, please cite:

```bibtex
@software{HERA_Matlab,
  author = {von Erdmannsdorff, Lukas},
  title = {HERA: A Hierarchical-Compensatory, Effect-Size Driven and Non-parametric
  Ranking Algorithm using Data-Driven Thresholds and Bootstrap Validation},
  url = {https://github.com/lerdmann1601/HERA-Matlab},
  version = {1.0.3},
  year = {2026}
}
```

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file
for details.

---

<!-- markdownlint-disable MD033 -->
<div align="center">
  <sub>Built by Lukas von Erdmannsdorff</sub>
</div>
<!-- markdownlint-enable MD033 -->
