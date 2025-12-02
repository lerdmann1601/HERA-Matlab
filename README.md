
<div align="center">

<img src="assets/hera_logo.svg" alt="HERA Logo" width="300"/>

# HERA: Hierarchical-Compensatory, Effect-Size-Driven and Non-Parametric Ranking Algorithm 

[![MATLAB](https://img.shields.io/badge/MATLAB-R2021a%2B-orange.svg)](https://www.mathworks.com/products/matlab.html)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![View on GitHub](https://img.shields.io/badge/GitHub-View_on_GitHub-blue?logo=github)](https://github.com/lerdmann1601/HERA-Matlab)
[![Issues](https://img.shields.io/github/issues/lerdmann1601/HERA-Matlab)](https://github.com/lerdmann1601/HERA-Matlab/issues)

**Made for Scientific Benchmarking**

[Key Features](#key-features) • [Installation](#installation) • [Quick Start](#quick-start) • [Documentation](#documentation) • [Citation](#citation)

</div>

---

## Overview

**HERA** is a MATLAB toolbox designed to automate the objective comparison of algorithms, experimental conditions, or datasets across multiple quality metrics. Unlike traditional ranking methods that rely solely on mean values or p-values, HERA employs a **hierarchical-compensatory logic** that integrates:

* **Significance Testing**: Wilcoxon signed-rank tests for paired data.
* **Effect Sizes**: Cliff's Delta and Relative Mean Difference for practical relevance.
* **Bootstrapping**: Data-driven thresholds and BCa confidence intervals.

This ensures that a "win" is only counted if it is both **statistically significant** and **practically relevant**, providing a robust and nuanced ranking system.

---

## Key Features

* **Hierarchical Logic**: Define primary and secondary metrics. Secondary metrics can act as tie-breakers or rank correctors (e.g., `M1_M2`, `M1_M2_M3`).
* **Data-Driven Thresholds**: Automatically calculates adaptive effect size thresholds using Percentile Bootstrapping.
* **Robustness**: Utilizes Bias-Corrected and Accelerated (BCa) confidence intervals and Cluster Bootstrapping for rank stability.
* **Automated Reporting**: Generates PDF reports, Win-Loss Matrices, Sankey Diagrams, and machine-readable JSON/CSV exports.
* **Reproducibility**: Supports fixed-seed execution and configuration file-based workflows.

---

## Quick Start

### 1. Interactive Mode (Recommended for Beginners)

The interactive wizard guides you through every step of the configuration, from data selection to statistical parameters.

```matlab
import HERA.start_ranking
HERA.start_ranking()
```

### 2. Batch Mode (Reproducible / Server)

For automated analysis or reproducible research, use a JSON configuration file.

```matlab
import HERA.start_ranking
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

---

## Installation

### Requirements

* **MATLAB** (R2021a or later recommended)
* **Statistics and Machine Learning Toolbox** (Required)
* **Parallel Computing Toolbox** (Recommended for performance)

### Setup

1. **Clone the repository:**

    ```bash
    git clone https://github.com/lerdmann1601/HERA-Matlab.git
    ```

2. **Add to MATLAB Path:**

    ```matlab
    addpath(genpath('/path/to/HERA-Matlab'));
    savepath;
    ```

---

## Documentation

<details>
<summary><strong>Ranking Modes Explained</strong></summary>

| Mode | Behavior | Use Case |
|------|----------|----------|
| `M1` | Ranks strictly by Metric 1. | Single-metric evaluation. |
| `M1_M2` | Metric 1 is primary. Metric 2 can correct (swap) ranks if a lower-ranked method significantly outperforms a higher-ranked one in Metric 2. | Balancing performance vs. cost. |
| `M1_M3A` | Metric 1 is primary. Metric 2 acts strictly as a tie-breaker. | Tie-breaking without overriding primary results. |
| `M1_M2_M3` | Full hierarchy. M1 is primary. M2 corrects M1 (iterative). M3 applies two sub-logics: (1) **One-time correction** if M2 is neutral, and (2) **Iterative tie-breaking** if both M1 and M2 are neutral. | Complex multi-objective benchmarking. |

</details>

<details>
<summary><strong>Input Data Specification</strong></summary>

* **Format**: CSV or Excel (`.xlsx`).
* **Organization**: One file per metric.
* **Filename**: Must match `metric_names` (e.g., `Accuracy.csv`).
* **Format**: CSV or Excel (`.xlsx`).
* **Organization**: One file per metric.
* **Filename**: Must match `metric_names` (e.g., `Accuracy.csv`).
* **Dimensions**: Rows = Subjects ($N$), Columns = Methods ($M$).
* **Consistency**: All files must have identical dimensions.

</details>

<details>
<summary><strong>Configuration Reference (Full)</strong></summary>

To run HERA in **Batch Mode**, create a `.json` file (e.g., `analysis_config.json`). Below is the complete list of available parameters.

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

**Parameter Dictionary:**

| Category | Parameter | Type | Default | Description |
|---|---|---|---|---|
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
| **Graphics** | `create_reports` | bool | `true` | Generate PDF reports and high-res plots. |
| | `plot_theme` | string | `"light"` | `"light"` or `"dark"`. |
| **Bootstrap (Manual)** | `manual_B_thr` | int | `2000` | Iterations for Thresholds (empty = auto). |
| | `manual_B_ci` | int | `5000` | Iterations for CIs (empty = auto). |
| | `manual_B_rank` | int | `500` | Iterations for Rank Stability (empty = auto). |
| **Bootstrap (Auto)** | `bootstrap_thresholds` | struct | (See Below) | Config for Threshold convergence. |
| | `bootstrap_ci` | struct | (See Below) | Config for CI convergence. |
| | `bootstrap_ranks` | struct | (See Below) | Config for Rank Stability convergence. |

**Detailed Bootstrap Configuration (Auto-Convergence):**

The `bootstrap_*` structs (e.g., `bootstrap_thresholds`) support the following nested parameters to control the convergence algorithm:

| Field | Type | Default (Thr/CI/Rank) | Description |
|---|---|---|---|
| `B_start` | int | `100` / `100` / `50` | Initial number of bootstrap iterations. |
| `B_step` | int | `100` / `200` / `10` | Iterations to add in each step. |
| `B_end` | int | `10000` / `20000` / `1500` | Maximum number of iterations. |
| `n_trials` | int | `25` / `25` / `15` | Number of independent trials per step to check stability. |
| `convergence_tolerance` | double | `0.005` / `0.01` / `0.005` | Max allowed variation (e.g., 0.005 = 0.5%). |
| `smoothing_window` | int | `3` / `4` / `3` | Window size for moving average smoothing. |
| `convergence_streak_needed` | int | `3` / `3` / `3` | Consecutive steps required to pass tolerance. |
| `min_steps_for_convergence_check` | int | `1` | Minimum steps before checking convergence. |

**Convergence Modes:**

* **Simple Convergence**: Used when `smoothing_window` is empty. Checks if the value changes less than `convergence_tolerance` between steps.
* **Robust Convergence**: Used when `smoothing_window` is set (default). Uses a moving average to smooth fluctuations and requires `convergence_streak_needed` consecutive stable steps.

</details>

<details>
<summary><strong>Advanced Usage (Developer Mode)</strong></summary>

Developers can call `HERA.run_ranking` directly with data matrices, bypassing file I/O. This is useful for integration into other pipelines or simulation studies.

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
userInput.dataset_names = {'Method A', 'Method B', 'Method C', 'Method D', 'Method E'};
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

<details>
<summary><strong>Results Structure Reference</strong></summary>

When running `results = HERA.run_ranking(...)`, the returned structure contains:

| Field | Dimensions | Description |
|---|---|---|
| **Final Results** | | |
| `final_rank` | $[N \times 1]$ | Final rank for each dataset (1 = Best). |
| `final_order` | $[1 \times N]$ | Indices of datasets sorted by rank. |
| `final_bootstrap_ranks` | $[N \times B]$ | Bootstrapped rank distribution for stability analysis. |
| **Statistics** | | |
| `d_vals_all` | $[Pairs \times M]$ | Cliff's Delta effect sizes for all pairs/metrics. |
| `rel_vals_all` | $[Pairs \times M]$ | Relative Mean Differences. |
| `ci_d_all` | $[Pairs \times 2 \times M]$ | BCa Confidence Intervals for Delta. |
| `ci_r_all` | $[Pairs \times 2 \times M]$ | BCa Confidence Intervals for RelDiff. |
| `all_p_value_matrices` | $\{1 \times M\}$ | Cell array of raw p-values (Wilcoxon). |
| `all_alpha_matrices` | $\{1 \times M\}$ | Cell array of Holm-Bonferroni corrected alphas. |
| `all_sig_matrices` | $\{1 \times M\}$ | Logical matrices indicating significant wins. |
| **Diagnostics** | | |
| `swap_details` | struct | Log of logic-based rank swaps (M1 vs M2). |
| `intermediate_orders` | struct | Rankings after each hierarchical stage. |
| `borda_results` | struct | Consensus ranking from sensitivity analysis. |
| `power_results` | struct | Post-hoc power analysis data. |
| `all_permutation_ranks` | $[N \times Perms]$ | Ranks for every metric permutation tested. |
| `selected_permutations` | $[Perms \times M]$ | Indices of metrics for each permutation. |

</details>

---

## Outputs

HERA generates a timestamped directory containing:

* **`Output/results.csv`**: Final ranking table (Mean ± SD).
* **`Output/analysis_data.json`**: Complete analysis record (Inputs, Config, Stats, Results) for matlab independent processing.
* **`Output/log.csv`**: Detailed log of pairwise comparisons and logic.
* **`Graphics/`**: High-res PNGs (Win-Loss Matrices, Sankey Diagrams).
* **`PDF/`**: Compiled reports (Ranking Report, Convergence Report).

---

## Testing

HERA includes a comprehensive scientific validation suite (`run_unit_test.m`) comprising **18 test cases**. These tests verify statistical accuracy (e.g., Wilcoxon exact p-values), logical robustness (e.g., outlier handling), and algorithm stability.

### Running Tests

You can run the test suite in three ways:

1. **Auto-Log Mode (Default)**
    Automatically finds a writable folder (e.g., Documents) to save the log file.

    ```matlab
    import HERA.run_unit_test
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

For reviewers or users without a local MATLAB license, you can run the test suite directly on GitHub:

1. Go to the **Actions** tab in this repository.
2. Select **Manual Test Suite** from the left sidebar.
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
  title = {HERA: Hierarchical-Compensatory Ranking Algorithm (MATLAB)},
  url = {https://github.com/lerdmann1601/HERA-Matlab},
  version = {1.0.0},
  year = {2025}
}
```

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

<div align="center">
  <sub>Built by Lukas von Erdmannsdorff</sub>
</div>
