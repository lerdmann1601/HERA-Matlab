# Configuration Reference

To run HERA with a configuration file (for a single reproducible run, automated scripts, or batch processing), create a `.json` file (e.g.,
`analysis_config.json`). Below is the complete list of available parameters.

> [!IMPORTANT]
> **Use Default Parameters:** It is strongly recommended to use the provided
> default parameters to ensure methodological consistency. If different parameters
> are used, they should be scientifically justified and clearly documented in
> your research methodology.

## Minimal Configuration

```json
{
  "userInput": {
    "folderPath": "/Users/Name/Data",
    "fileType": ".csv",
    "metric_names": ["Accuracy", "Precision", "Throughput"],
    "ranking_mode": "M1_M2_M3",
    "output_dir": "/Users/Name/Results"
  }
}
```

## Full Configuration Example (JSON)

This example shows **all** JSON configuration parameters, using `2000`, `10000`, and `500` as the recommended manual defaults for `manual_B_*` (note that `HERA.default()` initializes `manual_B_*` as `[]` to enable Auto-Convergence by default; see the Parameter Dictionary below for programmatic/in-memory options such as `custom_data`, `dataset_names`, `create_csvs`, and `quiet_mode`).
Parameters inside `system` and `bootstrap_*` must be nested correctly as shown. Note that internal runtime tuning heuristics inside `system` (`target_memory`, `jack_parfor_thr`, `jack_vec_limit`, `delta_mat_limit`, `min_batch_size`) are intentionally stripped from the exported `configuration.json` and `results.config.system` after initialization.

```json
{
  "userInput": {
    "folderPath": "/Path/To/Data",
    "fileType": ".csv",
    "metric_names": ["Metric1", "Metric2", "Metric3"],
    "output_dir": "/Path/To/Results",
    "language": "en",
    "check_updates": true,
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

    "system": {
      "target_memory": "auto",
      "jack_parfor_thr": 300,
      "jack_vec_limit": 150,
      "delta_mat_limit": 30000,
      "min_batch_size": 100
    },

    "manual_B_thr": 2000,
    "manual_B_ci": 10000,
    "manual_B_rank": 500,

    "bootstrap_seed_offset": 1000,

    "bootstrap_thresholds": {
      "B_start": 100,
      "B_step": 100,
      "B_end": 15000,
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
    },

    "selected_permutations": []
  }
}
```

## Parameter Dictionary

| Category | Parameter | Type | Default | Description |
| :--- | :--- | :--- | :--- | :--- |
| **Input/Output** | `folderPath` | string | - | Absolute path to data folder (File-Based Mode). |
| | `fileType` | string | - | `.csv` or `.xlsx` (File-Based Mode). |
| | `custom_data` | cell / list | - | **(Programmatic/In-Memory Mode)** Cell array (or Python list of arrays/DataFrames) of `[N_Subjects x N_Datasets]` matrices, bypassing file loading. |
| | `dataset_names` | array | `[]` | **(Programmatic/In-Memory Mode)** Optional column/dataset names matching columns in `custom_data` (auto-generates `D1, D2, ...` if omitted). |
| | `metric_names` | array | - | List of metric/file names in hierarchical order. *(Note: Higher values indicate better performance. Otherwise, invert metrics prior to input).* |
| | `output_dir` | string | `pwd` | Path to save results. |
| | `language` | string | `"en"` | Output language code. |
| **Logic** | `ranking_mode` | string | `"M1_M2_M3"` | Logic mode (`M1`, `M1_M2`, `M1_M3A`, `M1_M2_M3`). |
| | `run_sensitivity_analysis` | bool | `true` | Run ranking for all metric permutations. |
| | `selected_permutations` | array | `[]` | List of custom permutations to test (empty = all). |
| | `run_power_analysis` | bool | `true` | Run resampling-based empirical power analysis. |
| | `min_data_completeness` | double | `0.80` | Min fraction of valid pairs required (0.8 = 80%). |
| **Statistics** | `ci_level` | double | `0.95` | Confidence interval level (e.g., 0.95 for 95%). |
| | `alphas` | array | `[0.05, ...]` | Significance level for each metric. |
| | `power_simulations` | int | `10000` | Number of simulations for power analysis. |
| **System** | `check_updates` | bool | `true` | Check GitHub Releases for newer HERA versions on startup. |
| | `reproducible` | bool | `true` | Use fixed RNG seed. |
| | `seed` | int | `123` | RNG seed value. |
| | `bootstrap_seed_offset` | int | `1000` | Base random seed offset for parallel bootstrap operations. Prevents RNG substream collisions. |
| | `num_workers` | int/str | `"auto"` | Number of parallel workers. `"auto"` uses `parcluster('local').NumWorkers`. |
| | `system.target_memory` | int/str | `"auto"` | Target memory per chunk (MB). Computed from RAM when `"auto"`; stripped from saved `configuration.json` / `results.config.system`. |
| | `system.jack_parfor_thr` | int | `300` | Min $n$ to trigger parallel execution (stripped from saved `configuration.json` / `results.config.system`). |
| | `system.jack_vec_limit` | int | `150` | Max $n$ for vectorized Jackknife calculations (stripped from saved `configuration.json` / `results.config.system`). |
| | `system.delta_mat_limit` | int | `30000` | Max $n_x \cdot n_y$ product for matrix-based Cliff's Delta (stripped from saved `configuration.json` / `results.config.system`). |
| | `system.min_batch_size` | int | `100` | Min batch size for parallel processing (stripped from saved `configuration.json` / `results.config.system`). |
| **Output & Graphics** | `create_reports` | bool | `true` | Generate PDF reports and high-res plots (including convergence plots). If `false`, no plots are saved. |
| | `plot_theme` | string | `"light"` | `"light"`, `"dark"`, `"colourblind light"`, or `"colourblind dark"`. |
| | `create_csvs` | bool | `~quiet_mode` | Write intermediate BCa/Bootstrap CSV tables in `calculate_bca_ci` and `bootstrap_ranking` (defaults to `true` unless `quiet_mode` is `true`). |
| | `quiet_mode` | bool | `false` | Suppress verbose bootstrap console output and interactive pauses/prompts during batch or programmatic runs. |
| **Bootstrap (Manual)** | `manual_B_thr` | int/[] | `[]` (`2000` manual) | Iterations for Thresholds (`[]` = Auto-Convergence default; `2000` recommended for manual mode). |
| | `manual_B_ci` | int/[] | `[]` (`10000` manual) | Iterations for CIs (`[]` = Auto-Convergence default; `10000` recommended for manual mode). |
| | `manual_B_rank` | int/[] | `[]` (`500` manual) | Iterations for Rank Stability (`[]` = Auto-Convergence default; `500` recommended for manual mode). |
| **Bootstrap (Auto)** | `bootstrap_thresholds` | struct | (See Below) | Config for Threshold convergence. |
| | `bootstrap_ci` | struct | (See Below) | Config for CI convergence. |
| | `bootstrap_ranks` | struct | (See Below) | Config for Rank Stability convergence. |

👉 [Bootstrap Configuration (Auto-Convergence)](https://lerdmann1601.github.io/HERA-Matlab/Bootstrap_Configuration)
