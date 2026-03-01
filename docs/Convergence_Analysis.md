# Convergence Analysis for Robust Mode

This document provides an overview of the validation performed for the **Robust Convergence Mode** in HERA.

## Overview

The purpose of this analysis was to ensure that the default parameters used in the robust convergence mode are effective across a variety of data scenarios. Rather than a formal theoretical study, this was a practical validation using Monte Carlo simulations to verify that convergence is reliably achieved within the default iteration limits across typical use cases.

Additionally, we evaluated the accuracy of the results by comparing them against reference values generated with a very high number of bootstrap iterations. This ensures that the convergence criteria not only stop at a reasonable time but also produce results with acceptable error relative to the reference values.

The results confirm that the default settings are appropriate for the vast majority of scenarios, yielding both stable and accurate outcomes.

To ensure **reproducibility** and fair comparability between the Robust Mode and the reference runs, we utilized a **fixed random seed** for all simulations. This guarantees that both the convergence checks and the reference benchmarks were performed on identical data sequences and bootstrap samples.

## Analysis Design & Constraints

To validate the robustness of the default parameters, we simulated a diverse set of conditions.

### Fixed Parameters

It is important to note that the following parameters were treated as **constants** and thus their influence was not explicitly tested in this analysis:

* **Step Size ($B_{step}$)**
* **Tolerance ($Tol$)**
* **Start Iterations ($B_{start}$)**
* **End Iterations ($B_{end}$)**

These were held constant because the three parameters we *did* test (regarding data distribution and method variability) were considered the most critical factors influencing the stability of achieving a desired tolerance.

### Evaluated Scenarios

We tested against multiple distribution types (Normal, Bimodal, Skewed, Likert) and sample sizes, as well as varying method properties.

**Data Scenarios:**
![Scenario Parameters](https://raw.githubusercontent.com/lerdmann1601/HERA-Matlab/main/tests/Robustness_Report_Example/Graphics/Param_Scenarios_Example.png)

**Method Configurations:**
![Method Parameters](https://raw.githubusercontent.com/lerdmann1601/HERA-Matlab/main/tests/Robustness_Report_Example/Graphics/Param_Methods_Example.png)

### Analysis Scope

To isolate the behavior of individual convergence algorithms, the following simplifications were applied:

* **Effect Size Metrics**: Both Cliff's Delta and Relative Difference are always calculated, and the convergence criterion is based on the averaged stability of both. Only Cliff's Delta was used as the accuracy metric for comparison against reference values. This design choice allows us to test whether the averaged convergence approach produces accurate individual estimates despite relying on a combined stability indicator.
* **Ranking Stability**: A single-metric ranking (M1 only) was used to provide a stable reference point. Multi-metric rankings (M1_M2, M1_M2_M3) may exhibit different convergence characteristics due to increased complexity.
* **Thresholds for Ranking**: Pre-calculated reference thresholds were used to isolate the ranking convergence behavior from threshold estimation variability.

## Validation Results

We evaluated the performance of the Robust Mode based on two key criteria: **Convergence Success** (did it finish?) and **Accuracy** (were the results correct?).

### 1. Convergence Frequency

In our analysis, we observed that:

* The **Robust Mode** achieved exceptionally high convergence rates.
* The **Ranking** phase was the only area where non-convergence was observed, and even then, it occurred in **< 0.3%** of ranking cases.
* This implies that for all other metrics (BCa Intervals, Thresholds), convergence was virtually 100%.
* In the rare ranking failures, the maximum number of iterations (`B_end`) was reached before the convergence criterion was satisfied. To resolve this we increased `B_end` for the default settings provided.

### 2. Accuracy Validation

Beyond simply checking for convergence, we also validated the **accuracy** of the results.

* We compared the outcomes of the Robust Mode against "Gold Standard" reference values generated with very high bootstrap iteration counts:
  * **Thresholds**: $B_{ref} = 15{,}000$
  * **BCa Confidence Intervals**: $B_{ref} = 30{,}000$
  * **Ranking Stability**: $B_{ref} = 5{,}000$
* The goal was to analyze the **distribution of errors** (deviation from reference) to ensure that the convergence method yields results that are not just stable, but also practically in line with the theoretical limits.
* The results confirmed that the Robust Mode produces estimates with acceptably small deviation from reference values across all simulated scenarios. Median absolute errors were typically well below 5% relative to the reference values. This also validates the averaged convergence approach: despite determining convergence from a combined stability indicator, the individual accuracy of Cliff's Delta remained excellent.

### Global Results Summary

**Ranking Convergence:**
![Ranking Summary](https://raw.githubusercontent.com/lerdmann1601/HERA-Matlab/main/tests/Robustness_Report_Example/Graphics/Global_Summary_Ranking_Example.png)

**BCa Confidence Interval Convergence:**
![BCa Summary](https://raw.githubusercontent.com/lerdmann1601/HERA-Matlab/main/tests/Robustness_Report_Example/Graphics/Global_Summary_BCa_Example.png)

**Threshold Calculations:**
![Thresholds Summary](https://raw.githubusercontent.com/lerdmann1601/HERA-Matlab/main/tests/Robustness_Report_Example/Graphics/Global_Summary_Thresholds_Example.png)

## Discussion: Convergence Modes

Based on these results, we can discuss the theoretical implications of the different convergence parameters one could set in HERA compared to the default settings:

* **Relaxed**: Prioritizes speed. While significantly faster, it carries a higher risk of terminating before true stability is reached, especially in complex scenarios and BCa confidence intervals.
* **Strict**: Prioritizes guaranteed stability. It enforces rigorous checks that may lead to very high iteration counts. While safer, this can be computationally expensive and potentially overkill for well-behaved datasets where the Default mode would suffice. Also it could lead to non-convergence in some cases for the desired tolerance.
* **Default**: Designed to balance efficiency and reliability, it seems to be the best option for general usage. Interestingly, in some tested scenarios it provides higher accuracy for the BCa confidence intervals compared to stricter settings when comparing to the reference values. This is likely because stricter settings with more iterations can occasionally include rare extreme bootstrap samples that increase variance without improving the central estimate.

## Conclusion

Our analysis confirms that the Default settings provided with HERA should be a safe choice for general usage, achieving convergence in 99.7% of tested Ranking cases and 100% of other cases, while maintaining high accuracy against reference standards.
  
## Computational Note

The complete analysis (7 scenarios × 50 independent simulations × 6 different candidates ranked per scenario, comprising 2,100 distinct datasets) was performed on a standard consumer laptop (Base Model Apple 16" 2021 M1 MBP, 16 GB RAM) in approximately 11.5 hours. For each simulation, high-precision reference values were independently derived (up to B=30,000 for BCa), followed by convergence tests for all three modes (Relaxed, Default, Strict) utilizing 10–40 stability trials per B-step. This demonstrates the computational efficiency of the parallelized MATLAB implementation.

## Full Report

For a comprehensive look at the data, including specific breakdowns for each distribution type, please refer to the generated PDF report:

[**Download Full Combined Report (PDF)**](https://raw.githubusercontent.com/lerdmann1601/HERA-Matlab/main/tests/Robustness_Report_Example/Full_Combined_Report_Example.pdf)

## Running Your Own Analysis

You can perform this analysis yourself to verify the robustness of your own changes or to test different parameter configurations. The analysis is integrated directly into the `start_ranking` command.

### Syntax

```matlab
% Standard executions (Interactive CLI)
HERA.start_ranking('convergence', 'true', [options])

% Advanced executions (JSON Configurations)
HERA.start_ranking('convergence', 'path/to/my_config.json')
```

### Options

* `'logPath'`: (Optional) Specifies where the results and reports should be saved.
  * `'interactive'`: Opens a folder selection dialog.
  * `'path/to/folder'`: Uses the specified path.
  * `''` (default): Automatically saves to `Documents/HERA_convergence_Log`.
* `'sims'`: (Optional) Number of simulations per scenario (Default: 15). Higher values (e.g., 50-100) provide more robust statistics but take longer to run.

### Examples

**1. Run with Default Settings:**
This runs the standard analysis (15 sims/scenario) and saves to the default Documents folder.

```matlab
HERA.start_ranking('convergence', 'true')
```

**2. Select Output Folder Interactively:**
This prompts you to choose where to save the reports.

```matlab
HERA.start_ranking('convergence', 'true', 'logPath', 'interactive')
```

**3. Run a High-Precision Study:**
Run 50 simulations per scenario for higher statistical power.

```matlab
HERA.start_ranking('convergence', 'true', 'sims', 50)
```

**4. Run with a Custom JSON Configuration:**
For full control over the analysis parameters like RAM target, random seed offsets or simulation quantities, you can build a `.json` configuration file. This skips the standard `start_ranking` options entirely.

*Note: All parameters in the JSON configuration are entirely **optional**. The minimal required JSON file is just an empty `userInput` object (`{ "userInput": {} }`), from which standard defaults will fill all missing values.*

```matlab
HERA.start_ranking('convergence', 'path/to/convergence_config.json')
```

**Example JSON Configuration file:**

```json
{
  "userInput": {
     "n_sims_per_cond": 50,
     "output_dir": "/Users/Name/Results",
     "target_memory": 200,
     "simulation_seed": 123,
     "bootstrap_seed_offset": 1000,
     "scenario_seed_offset": 10000,
     "reference_seed_offset": 5000,
     "modes": {
         "Default": {
              "thr": {
                  "n": 50,
                  "sm": 5,
                  "st": 5,
                  "tol": 0.005,
                  "step": 200
              },
              "bca": {
                  "start": 500,
                  "step": 500
              },
              "rnk": {
                  "tol": 0.001
              }
         }
     }
  }
}
```

* `target_memory`: Overrides the dynamic RAM detection algorithm with a fixed limit (in MB).
  * **What it means:** It caps the memory footprint allowed per parallel worker during bootstrap batch processing. If the estimated memory for all iterations exceeds this limit, they are automatically split into smaller chunks to prevent out-of-memory errors.
  * **Default Calculation:** Because the convergence analysis uses `parfeval` for nested parallel execution (which carries slightly higher memory overhead), it applies a more conservative baseline than regular HERA runs. If omitted, it detects the system's total physical RAM and sets a base target of **25 MB per 1 GB of RAM** (e.g., 400 MB for a 16 GB system), with an unshakeable minimum of 200 MB. This base value is then divided by the number of CPU cores to determine the final effective memory limit per batch on each worker.
* `simulation_seed`: Seeds the main study. Can be locked for strict reproduction (fallback: 123).
* `bootstrap_seed_offset`: Keeps the internal evaluations structurally independent of the outer iterations (fallback: 1000).
* `scenario_seed_offset`: Ensures different tested scenarios (like Normal vs LogNormal) don't overlap in their RNG seeds (fallback: 10000).
* `reference_seed_offset`: Establishes the gap between simulated data generated and reference calculations executed (fallback: 5000).
* `modes`: You can specify missing details (`Relaxed`, `Default`, `Strict`) specifically for each of the three algorithms (`thr` for Thresholds, `bca` for BCa CI, `rnk` for Ranking), down to single parameters (`start`, `step`, `end`, `sm`, `st`, `tol`). Missing properties will automatically fall back to the package defaults. 👉 [Bootstrap Configuration](https://lerdmann1601.github.io/HERA-Matlab/Bootstrap_Configuration)
