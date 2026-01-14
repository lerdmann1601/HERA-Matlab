# Bootstrap Configuration (Auto-Convergence)

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

👉 [Convergence Modes and Troubleshooting](https://github.com/lerdmann1601/HERA-Matlab/blob/main/docs/Convergence_Modes.md)
