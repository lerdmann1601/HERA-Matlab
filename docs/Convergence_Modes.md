# Convergence Modes

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
> [tests/Convergence_Analysis.md](https://github.com/lerdmann1601/HERA-Matlab/blob/main/tests/Convergence_Analysis.md).

## Troubleshooting Convergence

While the automated check should work for most datasets, "difficult" data with high
variance or flat likelihood landscapes may fail to converge within
`B_end`. In this case, you can try the following:

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
