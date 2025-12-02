
<div align="center">

<img src="assets/hera_logo.svg" alt="HERA Logo" width="300"/>

# HERA: Hierarchical-Compensatory, Effect-Size-Driven and Non-Parametric Ranking Algorithm 

[![MATLAB](https://img.shields.io/badge/MATLAB-R2021a%2B-orange.svg)](https://www.mathworks.com/products/matlab.html)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![View on GitHub](https://img.shields.io/badge/GitHub-View_on_GitHub-blue?logo=github)](https://github.com/lerdmann1601/HERA-Matlab)
[![Issues](https://img.shields.io/github/issues/lerdmann1601/HERA-Matlab)](https://github.com/lerdmann1601/HERA-Matlab/issues)

**HERA is made for Scientific Benchmarking**

[Key Features](#key-features) • [Installation](#installation) • [Quick Start](#quick-start) • [Documentation](#documentation) • [Citation](#citation)

</div>

---

## Overview

**HERA** is a MATLAB toolbox designed to automate the objective comparison of algorithms, experimental conditions, or datasets across multiple quality metrics. Unlike traditional ranking methods that rely solely on mean values or p-values, HERA employs a **hierarchical-compensatory logic** that integrates:

*   **Significance Testing**: Wilcoxon signed-rank tests for paired data.
*   **Effect Sizes**: Cliff's Delta and Relative Mean Difference for practical relevance.
*   **Bootstrapping**: Data-driven thresholds and BCa confidence intervals.

This ensures that a "win" is only counted if it is both **statistically significant** and **practically relevant**, providing a robust and nuanced ranking system.

---

## Key Features

*   **Hierarchical Logic**: Define primary and secondary metrics. Secondary metrics can act as tie-breakers or rank correctors (e.g., `M1_M2`, `M1_M2_M3`).
*   **Data-Driven Thresholds**: Automatically calculates adaptive effect size thresholds using Percentile Bootstrapping.
*   **Robustness**: Utilizes Bias-Corrected and Accelerated (BCa) confidence intervals and Cluster Bootstrapping for rank stability.
*   **Automated Reporting**: Generates PDF reports, Win-Loss Matrices, Sankey Diagrams, and machine-readable JSON/CSV exports.
*   **Reproducibility**: Supports fixed-seed execution and configuration file-based workflows.

---

## Quick Start

### Interactive Mode
For first-time users, the interactive wizard is the easiest way to start. It guides you through data selection and parameter configuration.

```matlab
import HERA.start_ranking
HERA.start_ranking()
```

### Batch Mode
For reproducible studies, run HERA with a configuration file or structure.

```matlab
% Example: Run with a config file
import HERA.start_ranking
HERA.start_ranking('configFile', 'config.json')
```

---

## Installation

### Requirements
*   **MATLAB** (R2021a or later recommended)
*   **Statistics and Machine Learning Toolbox** (Required)
*   **Parallel Computing Toolbox** (Recommended for performance)

### Setup
1.  **Clone the repository:**
    ```bash
    git clone https://github.com/lerdmann1601/HERA-Matlab.git
    ```
2.  **Add to MATLAB Path:**
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

*   **Format**: CSV or Excel (`.xlsx`).
*   **Organization**: One file per metric.
*   **Filename**: Must match `metric_names` (e.g., `Accuracy.csv`).
*   **Dimensions**: Rows = Subjects ($N$), Columns = Methods ($M$).
*   **Consistency**: All files must have identical dimensions.

</details>

<details>
<summary><strong>Configuration Reference</strong></summary>

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `folderPath` | string | - | Path to input data. |
| `metric_names` | cell | - | Names of metrics (filenames). |
| `ranking_mode` | string | `M1` | Logic mode. |
| `reproducible` | bool | `true` | Fixes RNG seed. |
| `create_reports` | bool | `true` | Generates PDF reports. |

</details>

---

## Outputs

HERA generates a timestamped directory containing:

*   **`Output/results.csv`**: Final ranking table (Mean ± SD).
*   **`Output/log.csv`**: Detailed log of pairwise comparisons and logic.
*   **`Graphics/`**: High-res PNGs (Win-Loss Matrices, Sankey Diagrams).
*   **`PDF/`**: Compiled reports (Ranking Report, Convergence Report).

---

## Testing

Validate your installation with the included unit test suite (18 tests covering edge cases).

```matlab
HERA.start_ranking('runtest', true)
```

### GitHub Actions (Cloud Testing)
For reviewers or users without a local MATLAB license, you can run the test suite directly on GitHub:
1.  Go to the **Actions** tab in this repository.
2.  Select **Manual Test Suite** from the left sidebar.
3.  Click **Run workflow**.


---

## Contributing

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for details.

1.  Fork the repository.
2.  Create a feature branch.
3.  Commit your changes.
4.  Open a Pull Request.

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
