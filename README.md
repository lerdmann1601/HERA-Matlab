
<div align="center">

<img src="https://raw.githubusercontent.com/lerdmann1601/HERA-Matlab/main/assets/hera_logo.svg" alt="HERA Logo" width="300"/>

# HERA: Hierarchical-Compensatory, Effect-Size-Driven and Non-Parametric Ranking Algorithm

[![MATLAB](https://img.shields.io/badge/MATLAB-R2024a%2B-orange.svg)](https://www.mathworks.com/products/matlab.html)
[![Statistics Toolbox](https://img.shields.io/badge/Toolbox-Statistics_and_Machine_Learning-blue.svg)](https://www.mathworks.com/products/statistics.html)
[![Parallel Computing Toolbox](https://img.shields.io/badge/Toolbox-Parallel_Computing-blue.svg)](https://www.mathworks.com/products/parallel-computing.html)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![View on GitHub](https://img.shields.io/badge/GitHub-View_on_GitHub-blue?logo=github)](https://github.com/lerdmann1601/HERA-Matlab)
[![Issues](https://img.shields.io/github/issues/lerdmann1601/HERA-Matlab)](https://github.com/lerdmann1601/HERA-Matlab/issues)
[![PyPI Downloads](https://static.pepy.tech/personalized-badge/hera-matlab?period=total&units=INTERNATIONAL_SYSTEM&left_color=BLACK&right_color=GREEN&left_text=PyPI%20Downloads)](https://pepy.tech/projects/hera-matlab)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.18274871-blue.svg)](https://doi.org/10.5281/zenodo.18274871)
[![ORCID](https://img.shields.io/badge/ORCID-0009--0009--3758--7363-green.svg)](https://orcid.org/0009-0009-3758-7363)

**Made for Scientific Benchmarking**

[Key Features](#key-features) • [Installation](#installation) • [Quick Start](#quick-start) • [Documentation](#documentation) • [Citation](#citation)

</div>

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

👉 [Standalone Runtime](https://github.com/lerdmann1601/HERA-Matlab/blob/main/docs/Standalone_Runtime.md)

👉 [Python Integration](https://github.com/lerdmann1601/HERA-Matlab/blob/main/docs/Python_Integration.md)

👉 [Automated Build (GitHub Actions)](https://github.com/lerdmann1601/HERA-Matlab/blob/main/docs/Automated_Build.md)

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
> provided in the `data/examples` directory. See [docs/Example_Analysis.md](https://github.com/lerdmann1601/HERA-Matlab/blob/main/docs/Example_Analysis.md) for a
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

👉 [Repository Structure](https://github.com/lerdmann1601/HERA-Matlab/blob/main/docs/Repository_Structure.md)

👉 [Theoretical Background](https://github.com/lerdmann1601/HERA-Matlab/blob/main/paper/paper.pdf)

👉 [Ranking Modes Explained](https://github.com/lerdmann1601/HERA-Matlab/blob/main/docs/Ranking_Modes_Explained.md)

👉 [Input Data Specification](https://github.com/lerdmann1601/HERA-Matlab/blob/main/docs/Input_Data_Specification.md)

👉 [Example Analysis](https://github.com/lerdmann1601/HERA-Matlab/blob/main/docs/Example_Analysis.md)

👉 [Methodological Guidelines & Limitations](https://github.com/lerdmann1601/HERA-Matlab/blob/main/docs/Methodological_Guidelines_&_Limitations.md)

👉 [Configuration & Parameters](https://github.com/lerdmann1601/HERA-Matlab/blob/main/docs/Configuration_&_Parameters.md)

👉 [Bootstrap Configuration](https://github.com/lerdmann1601/HERA-Matlab/blob/main/docs/Bootstrap_Configuration.md)

👉 [Convergence Modes](https://github.com/lerdmann1601/HERA-Matlab/blob/main/docs/Convergence_Modes.md)

👉 [Convergence Analysis](https://github.com/lerdmann1601/HERA-Matlab/blob/main/docs/Convergence_Analysis.md)

👉 [Advanced Usage (Developer Mode)](https://github.com/lerdmann1601/HERA-Matlab/blob/main/docs/Advanced_Usage.md)

👉 [Results Structure Reference](https://github.com/lerdmann1601/HERA-Matlab/blob/main/docs/Results_Structure_Reference.md)

---

## Outputs

HERA generates a timestamped directory containing:

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

We welcome contributions! Please see [CONTRIBUTING.md](https://github.com/lerdmann1601/HERA-Matlab/blob/main/CONTRIBUTING.md) for details.

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
  version = {1.1.2},
  doi = {10.5281/zenodo.18274871},
  year = {2026}
}
```

---

## License

This project is licensed under the MIT License - see the [LICENSE](https://github.com/lerdmann1601/HERA-Matlab/blob/main/license.txt) file
for details.

---
