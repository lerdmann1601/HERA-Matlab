# HERA: Hierarchical-Compensatory Ranking Algorithm (MATLAB)

HERA is a MATLAB toolbox designed for the objective, non-parametric, and effect-size-driven ranking of multiple datasets (e.g., comparison of algorithms, medical imaging, or experimental conditions). It employs a hierarchical logic that prioritizes primary metrics while allowing secondary metrics to resolve ties or correct rankings based on practical significance.

---

## Summary

In many scientific disciplines, there is a need to objectively compare the performance of different methods across multiple quality metrics. Traditional ranking often relies on simple mean values or p-values, which may fail to capture practical significance or robustness against outliers.

HERA addresses this by implementing a sequential ranking algorithm based on:

- **Significance Tests**: Wilcoxon signed-rank tests for paired data.
- **Effect Sizes**: Cliff's Delta for stochastic dominance and Relative Mean Difference for practical relevance.
- **Bootstrapping**: Data-driven determination of thresholds and confidence intervals (BCa) to ensure statistical validity without arbitrary cutoffs.

---

## Statement of Need

Objective benchmarking is critical in research. However, manually comparing methods across multiple metrics often leads to subjective weighting schemes (e.g., "Method A is better because metric X is slightly higher, even though metric Y is lower"). Furthermore, relying solely on p-values can be misleading in large datasets, where trivial differences become significant, or in small datasets, where relevant effects are missed.

HERA fills this gap by providing a formalized, hierarchical-compensatory ranking system. It ensures that a "win" is only counted if it is both statistically significant (Holm-Bonferroni corrected) and practically relevant (exceeding data-driven effect size thresholds). This allows researchers to obtain a nuanced and robust ranking of methods, facilitating better decision-making and fairer comparisons in benchmarking studies.

---

## Key Features

- **Hierarchical Logic**: Supports flexible ranking modes (e.g., M1, M1_M2, M1_M2_M3) where secondary metrics can act as tie-breakers or correctors.
- **Non-Parametric**: Does not assume normal distribution of data; uses robust statistics.
- **Robustness**: Includes Bias-Corrected and Accelerated (BCa) bootstrap confidence intervals and Cluster Bootstrapping for rank stability analysis.
- **Data-Driven Thresholds**: Uses bootstrapping to automatically determine what constitutes a "relevant" difference based on the data's variability.
- **Automated Reporting**: Generates comprehensive PDF reports, high-resolution plots (Win-Loss Matrices, Sankey Diagrams), and machine-readable JSON/CSV exports.
- **Reproducibility**: Supports fixed-seed execution for consistent results.

---

## Installation

### Requirements

- MATLAB (R2021a or later recommended for tiledlayout support).
- Statistics and Machine Learning Toolbox (Required).
- Parallel Computing Toolbox (Recommended for faster bootstrapping).

### Setup

1. Clone the repository:

```bash
git clone https://github.com/lerdmann1601/HERA-Matlab.git
```

2. Add the repository folder (and subfolders) to your MATLAB path:

```matlab
addpath(genpath('path/to/HERA-Matlab'));
savepath;
```

---

## Usage

HERA supports **three operation modes** to accommodate different use cases: **Interactive Mode** (CLI-guided), **Batch/Server Mode** (config file-based), and **Unit Test Mode** (validation).

### 1. Interactive Mode (Quick Start)

Simply run the start command and follow the on-screen instructions to select your data and configure the analysis:

```matlab
import HERA.start_ranking
HERA.start_ranking()
```

**What happens:**
- You'll be guided through all settings via CLI prompts (reproducibility, parallel processing, data selection, ranking logic, etc.).
- Choose between **Standard Mode** (uses defaults) or **Manual Mode** (customize all parameters).
- You can also load a previously saved configuration file.

---

### 2. Batch/Server Mode (Reproducible Research)

For scientific publications or automated workflows, use a configuration file. This ensures your analysis is exactly reproducible.

#### Method A: Create a JSON Configuration File Manually (Non-MATLAB Users)

Create a plain text file with a `.json` extension (e.g., `myanalysis.json`) using any text editor:

```json
{
  "userInput": {
    "folderPath": "C:\\MyData\\Study1",
    "fileType": ".csv",
    "metric_names": ["DiceCoefficient", "HausdorffDistance"],
    "ranking_mode": "M1_M2",
    "output_dir": "C:\\MyResults",
    "reproducible": true,
    "seed": 123,
    "num_workers": 4,
    "plot_theme": "dark"
  }
}
```

**Important Notes:**
- Windows paths must use double backslashes (`\\`). macOS/Linux paths use single forward slashes (`/`).
- `metric_names` defines the hierarchical order of analysis.
- `ranking_mode` must match the number of metrics (see **Ranking Modes** below).

**Run the analysis:**

```matlab
import HERA.start_ranking
HERA.start_ranking('configFile', 'C:\MyData\myanalysis.json')
```

#### Method B: Create a JSON Configuration File Programmatically (MATLAB Users)

```matlab
% 1. Create configuration structure
userInput = struct();

% 2. Define Data Source
userInput.folderPath   = 'C:\MyData\Study1'; 
userInput.fileType     = '.csv';

% 3. Define Metric Hierarchy
userInput.metric_names = {'DiceCoefficient', 'HausdorffDistance'}; 

% 4. Select Ranking Logic
userInput.ranking_mode = 'M1_M2'; 

% 5. Output and Settings
userInput.output_dir   = 'C:\MyResults';
userInput.reproducible = true;
userInput.seed         = 123;
userInput.num_workers  = 4;

% 6. Save configuration as JSON
dataToSave = struct('userInput', userInput);
jsonText = jsonencode(dataToSave, 'PrettyPrint', true);
fid = fopen('C:\MyData\myanalysis.json', 'w');
fprintf(fid, '%s', jsonText);
fclose(fid);

% 7. Run analysis
import HERA.start_ranking
HERA.start_ranking('configFile', 'C:\MyData\myanalysis.json')
```

---

### 3. Unit Test Mode

To validate the installation and ensure the algorithm performs correctly on your system, run the included unit test suite:

```matlab
import HERA.start_ranking
HERA.start_ranking('runtest', true)
```

**Optional:** Specify a custom log path:

```matlab
HERA.start_ranking('runtest', true, 'logPath', 'C:\MyLogs')
```

This suite covers 18 statistical and logical test cases, including edge cases like zero-variance data, Efron's Dice (cycles), and missing data handling.

---

## 📊 Ranking Modes Explained

HERA supports **four ranking modes** that define how metrics interact in the hierarchical decision process:

| Mode | Description | Use Case |
|------|-------------|----------|
| **M1** | Ranks based **only** on Metric 1 (primary metric). Secondary metrics are ignored. | Single-metric benchmarking (e.g., accuracy only). |
| **M1_M2** | Metric 1 is primary. Metric 2 **corrects/swaps** ranks if differences are statistically and practically significant. | Dual-metric evaluation where the secondary metric can override the primary (e.g., accuracy vs. computational cost). |
| **M1_M3A** | Metric 1 is primary. Metric 2 acts as a **tie-breaker only** when Metric 1 shows no significant difference. | Dual-metric evaluation where the secondary metric refines ties without overriding the primary. |
| **M1_M2_M3** | Full 3-metric hierarchy. Metric 1 is primary, Metric 2 corrects, and Metric 3 acts as a tie-breaker. | Complex benchmarking with multiple quality dimensions (e.g., medical imaging: Dice, Hausdorff, runtime). |

**Selection Logic:**
- **1 metric selected** → `M1` (automatically set)
- **2 metrics selected** → Choose between `M1_M2` or `M1_M3A`
- **3 metrics selected** → `M1_M2_M3` (automatically set)

---

## 📁 Input Data Format

Data should be provided as **CSV** or **Excel** files.

**Structure:**
- **Rows**: Subjects/Samples (Observations).
- **Columns**: Methods/Algorithms/Datasets to be compared.
- **Files**: One file per metric.
- **Naming**: Filenames must match `metric_names` in the config (e.g., `DiceCoefficient.csv` and `HausdorffDistance.csv`).
- **Consistency**: All files must have the same number of rows (subjects) and columns (methods).

**Example:**

If you have 3 algorithms (A, B, C) and 50 test cases, each CSV file should be 50 rows × 3 columns.

---

## 📦 Outputs

HERA creates a timestamped folder in your output directory containing:

### Folder Structure

```
Ranking_20231201_143022/
├── PDF/
│   ├── RankingReport.pdf
│   └── ConvergenceDiagnostics.pdf
├── Graphics/
│   ├── WinLossMatrix.png
│   ├── Ranking.png
│   ├── SankeyDiagram.png
│   └── ConvergencePlots.png
├── Output/
│   ├── results.csv          (Final ranking table with Mean ± SD)
│   ├── log.csv              (Detailed pairwise comparison log)
│   └── analysis_data.json   (Complete analysis data)
└── configuration.json       (Your input configuration)
```

### Key Files

- **`results.csv`**: Final ranking table with Mean ± SD for all metrics.
- **`log.csv`**: Detailed log of every pairwise comparison and the decision logic used.
- **`analysis_data.json`**: Full analysis data for programmatic post-processing (effect sizes, p-values, confidence intervals, etc.).
- **PDF Reports** (if enabled): Comprehensive visual summaries with Win-Loss Matrices, Sankey Diagrams, and ranking stability plots.

---

## Configuration Options Reference

Below is a comprehensive list of configuration parameters you can set in your `userInput` struct or JSON file:

### Essential Parameters

| Parameter | Type | Description | Example |
|-----------|------|-------------|---------|
| `folderPath` | string | Path to folder containing data files | `'C:\MyData\Study1'` |
| `fileType` | string | File extension (`.csv` or `.xlsx`) | `'.csv'` |
| `metric_names` | cell array | Metric names in hierarchical order | `{'Dice', 'Hausdorff'}` |
| `ranking_mode` | string | Ranking logic (`M1`, `M1_M2`, `M1_M3A`, `M1_M2_M3`) | `'M1_M2'` |
| `output_dir` | string | Path to save results | `'C:\MyResults'` |

### Reproducibility & Performance

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `reproducible` | boolean | `true` | Use fixed random seed for reproducibility |
| `seed` | integer | `42` | Random seed value |
| `num_workers` | integer or `'auto'` | `'auto'` | Number of parallel CPU cores |

### Statistical Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `alphas` | array | `[0.05]` | Significance level for each metric (e.g., `[0.05, 0.05, 0.05]`) |
| `ci_level` | float | `0.95` | Confidence interval level (0-1) |

### Bootstrap Configuration

For each of the three bootstrap analyses (thresholds, BCa CI, rank stability), you can choose:
- **Manual**: Fixed number of bootstrap samples (e.g., `manualBthr = 5000`)
- **Automatic**: Convergence search (`manualBthr = []` triggers automatic mode)

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `manualBthr` | integer or `[]` | `[]` | Manual B-value for threshold calculation (empty = auto) |
| `manualBci` | integer or `[]` | `[]` | Manual B-value for BCa CI (empty = auto) |
| `manualBrank` | integer or `[]` | `[]` | Manual B-value for rank stability (empty = auto) |

### Advanced Analyses

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `run_sensitivity_analysis` | boolean | `true` | Test alternative metric hierarchies (only for 2+ metrics) |
| `selected_permutations` | matrix | all | Which metric orders to test (auto-generated if empty) |
| `run_power_analysis` | boolean | `true` | Perform post-hoc power analysis |
| `power_simulations` | integer | `1000` | Number of power simulations |

### Graphics & Reporting

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `plot_theme` | string | `'light'` | Plot theme (`'light'` or `'dark'`) |
| `create_reports` | boolean | `true` | Generate PDF reports and full graphics (set to `false` for diagnostic mode) |

---

## Advanced Features

### Sensitivity Analysis

When enabled (automatically for 2+ metrics), HERA tests **all permutations** of your metric hierarchy to assess how sensitive the final ranking is to the order of metrics.

**Example:**
- Primary hierarchy: `Dice → Hausdorff → Runtime`
- Tested permutations: All 6 possible orders

The **Borda ranking** aggregates results across all permutations to produce a consensus ranking.

**How to disable:**

```matlab
userInput.run_sensitivity_analysis = false;
```

### Power Analysis

HERA performs a **post-hoc, non-parametric power analysis** to estimate the statistical power of your study. This helps you understand whether your sample size was sufficient to detect meaningful differences.

**How to configure:**

```matlab
userInput.run_power_analysis = true;
userInput.power_simulations = 1000; % Number of Monte Carlo simulations
```

### Diagnostic Batch Mode

For large-scale automated experiments, you can disable PDF generation and heavy plotting to speed up execution:

```matlab
userInput.create_reports = false;
```

**What's still saved:**
- All CSV and JSON data
- Convergence diagnostic plots
- Essential Win-Loss matrices

---

## Testing

To validate the installation and ensure the algorithm performs correctly on your system, run the included unit test suite:

```matlab
import HERA.start_ranking
HERA.start_ranking('runtest', true)
```

This suite covers 18 statistical and logical test cases, including edge cases like zero-variance data, Efron's Dice (cycles), and missing data handling.

---

## Example Workflows

### Example 1: Quick Start with Default Settings

```matlab
import HERA.start_ranking
HERA.start_ranking()
% Follow the CLI prompts and select "Standard Mode"
```

### Example 2: Reproducible Research with Custom Config

```matlab
userInput = struct();
userInput.folderPath = '/data/medical_imaging/study1';
userInput.fileType = '.csv';
userInput.metric_names = {'DiceCoefficient', 'HausdorffDistance', 'Runtime'};
userInput.ranking_mode = 'M1_M2_M3';
userInput.output_dir = '/results/study1';
userInput.reproducible = true;
userInput.seed = 42;
userInput.num_workers = 8;

% Save and run
dataToSave = struct('userInput', userInput);
jsonText = jsonencode(dataToSave, 'PrettyPrint', true);
fid = fopen('/data/study1_config.json', 'w');
fprintf(fid, '%s', jsonText);
fclose(fid);

import HERA.start_ranking
HERA.start_ranking('configFile', '/data/study1_config.json')
```

### Example 3: Batch Mode for High-Throughput Analysis

```matlab
userInput = struct();
userInput.folderPath = '/data/benchmark/algorithms';
userInput.fileType = '.csv';
userInput.metric_names = {'Accuracy', 'F1Score'};
userInput.ranking_mode = 'M1_M2';
userInput.output_dir = '/results/benchmark';
userInput.create_reports = false; % Diagnostic mode
userInput.num_workers = 16;

import HERA.start_ranking
HERA.start_ranking(userInput) % Direct struct input
```

---

## Contributing

Contributions are welcome! Please open an issue or submit a pull request on GitHub.

---

## License

This project is licensed under the MIT License - see the LICENSE file for details.

---

## Contact

**Author**: Lukas von Erdmannsdorff  
**GitHub**: [https://github.com/lerdmann1601/HERA-Matlab](https://github.com/lerdmann1601/HERA-Matlab)  
**Version**: 1.0
