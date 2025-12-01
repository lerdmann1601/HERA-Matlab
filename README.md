# HERA: Hierarchical-Compensatory Ranking Algorithm (MATLAB)

## Summary

HERA is a MATLAB toolbox designed for the objective, non-parametric, and effect-size-driven ranking of multiple datasets (e.g., comparison of algorithms, medical treatments, or experimental conditions). It employs a hierarchical logic that prioritizes primary metrics while allowing secondary metrics to resolve ties or correct rankings based on practical significance.

In many scientific disciplines, there is a need to objectively compare the performance of different methods across multiple quality metrics. Traditional ranking often relies on simple mean values or p-values, which may fail to capture practical significance or robustness against outliers.

HERA addresses this by implementing a sequential ranking algorithm based on:

- **Significance Tests**: Wilcoxon signed-rank tests for paired data.
- **Effect Sizes**: Cliff's Delta for stochastic dominance and Relative Mean Difference for practical relevance.
- **Bootstrapping**: Data-driven determination of thresholds and confidence intervals (BCa) to ensure statistical validity without arbitrary cutoffs.

## Statement of Need

Objective benchmarking is critical in research. However, manually comparing methods across multiple metrics often leads to subjective weighting schemes (e.g., "Method A is better because metric X is slightly higher, even though metric Y is lower"). Furthermore, relying solely on p-values can be misleading in large datasets, where trivial differences become significant, or in small datasets, where relevant effects are missed.

HERA fills this gap by providing a formalized, hierarchical-compensatory ranking system. It ensures that a "win" is only counted if it is both statistically significant (Holm-Bonferroni corrected) and practically relevant (exceeding data-driven effect size thresholds). This allows researchers to obtain a nuanced and robust ranking of methods, facilitating better decision-making and fairer comparisons in benchmarking studies.

## Key Features

- **Hierarchical Logic**: Supports flexible ranking modes (e.g., M1, M1_M2, M1_M2_M3) where secondary metrics can act as tie-breakers or correctors.
- **Non-Parametric**: Does not assume normal distribution of data; uses robust statistics.
- **Robustness**: Includes Bias-Corrected and Accelerated (BCa) bootstrap confidence intervals and Cluster Bootstrapping for rank stability analysis.
- **Data-Driven Thresholds**: Uses bootstrapping to automatically determine what constitutes a "relevant" difference based on the data's variability.
- **Automated Reporting**: Generates comprehensive PDF reports, high-resolution plots (Win-Loss Matrices, Sankey Diagrams), and machine-readable JSON/CSV exports.
- **Reproducibility**: Supports fixed-seed execution for consistent results.

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

## Usage

HERA can be used in two ways: Interactive Mode (GUI-like prompts in the terminal) or Script/Batch Mode (for reproducibility).

### 1. Interactive Mode (Quick Start)

Simply run the start command and follow the on-screen instructions to select your data and configure the analysis:

```matlab
import HERA.start_ranking
HERA.start_ranking()
```

### 2. Batch Mode (Reproducible Research)

For scientific publications, it is recommended to use a configuration struct. This ensures that your analysis is exactly reproducible.

```matlab
% 1. Create configuration structure
userInput = struct();

% 2. Define Data Source
% Path to the folder containing your .csv or .xlsx files
userInput.folderPath   = 'C:\MyData\Study1'; 
userInput.fileType     = '.csv';

% 3. Define Metric Hierarchy
% The order here defines the hierarchy (Metric 1 -> Metric 2 -> Metric 3)
% Filenames must match these names (e.g., 'DiceCoefficient.csv')
userInput.metric_names = {'DiceCoefficient', 'HausdorffDistance'}; 

% 4. Select Ranking Logic
% 'M1': Sort by Metric 1 only.
% 'M1_M2': Metric 1 is primary, Metric 2 corrects/swaps ranks if differences are significant.
% 'M1_M3A': Metric 2 acts as a tie-breaker only if Metric 1 is neutral.
userInput.ranking_mode = 'M1_M2'; 

% 5. Output and Settings
userInput.output_dir   = 'C:\MyResults';
userInput.reproducible = true;
userInput.seed         = 123; % Fixed seed for exact reproducibility
userInput.num_workers  = 4;   % Number of parallel workers

% 6. Start analysis
import HERA.start_ranking
HERA.start_ranking(userInput);
```

### Input Data Format

Data should be provided as CSV or Excel files.

**Structure:**
- Rows: Subjects/Samples (Observations).
- Columns: Methods/Algorithms/Datasets to be compared.
- Files: You need one file per metric.
- Naming: The filenames must match the metric_names provided in the config (e.g., DiceCoefficient.csv and HausdorffDistance.csv).
- Consistency: All files must have the same number of rows (subjects) and columns (methods).

## Outputs

HERA creates a timestamped folder in your output directory containing:

- **PDF/**: Comprehensive reports (Final Report, Convergence Diagnostics).
- **Graphics/**: High-resolution PNGs of all plots (Win-Loss Matrix, Ranking Bars, Sankey Diagrams).
- **Output/**:
  - `_results.csv`: Final ranking table with Mean ± SD.
  - `_log.csv`: Detailed log of every pairwise comparison and the decision logic used.
  - `_data.json`: Full analysis data for programmatic post-processing.

## Testing

To validate the installation and ensure the algorithm performs correctly on your system, run the included unit test suite. This suite covers 18 statistical and logical test cases, including edge cases like zero-variance data, Efron's Dice (cycles), and missing data handling.

```matlab
run('tests/run_tests.m')
```
