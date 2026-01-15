# Configuration Referenc

To run HERA in **Batch Mode**, create a `.json` file (e.g.,
`analysis_config.json`). Below is the complete list of available parameters.

## Minimal Configuration

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

## Full Configuration Example (JSON)

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

## Parameter Dictionary

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

👉 [Bootstrap Configuration (Auto-Convergence)](https://github.com/lerdmann1601/HERA-Matlab/blob/main/docs/Bootstrap_Configuration.md)
