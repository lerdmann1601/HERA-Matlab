# Convergence Modes

To ensure that bootstrap estimates (Thresholds, Confidence Intervals, and Ranking) are stable and reliable without requiring the user to specify a fixed number of iterations, HERA implements an *Adaptive Convergence Algorithm*. The process monitors estimator stability and terminates automatically when a plateau is reached. Crucially, the algorithm identifies a common *global number* of iterations $B$, ensuring that within each evaluation step (Thresholds, BCa intervals, and Ranking) estimates are measured with the same level of precision.

## Iterative Procedure

The algorithm operates by incrementally expanding the computational effort until the estimates achieve the desired precision. The default parameters for all analysis phases (Thresholds, CIs, and Rankings) are defined in the [Bootstrap Configuration](https://lerdmann1601.github.io/HERA-Matlab/Bootstrap_Configuration).

1. **Stepwise Increase**: Starting from a base number (`B_start`), the count of bootstrap iterations ($B$) is increased in fixed increments (`B_step`) at each checkpoint.
2. **Independent Trials**: To quantify stability at each step, HERA performs multiple independent trials (`n_trials`) of the bootstrap procedure. This allows the algorithm to distinguish between true convergence and random sampling artifacts.
3. **Adaptive Termination**: The process continues until the stability metrics, calculated across these independent trials, satisfy the chosen convergence mode's criteria or the maximum iteration count (`B_max` / `B_end`) is reached.

## Measuring Stability

Before applying any convergence criteria, HERA quantifies the global stability of the bootstrap estimates. The calculation differs depending on the analysis phase:

* **Thresholds**: Stability is measured using the *Robust Coefficient of Quartile Variation (RCQV)* across the independent trials.
* **BCa Confidence Intervals**: HERA employs a *two-stage aggregation* for robust estimation. First, the RCQV is calculated for the CI width of each individual pairwise comparison across trials. Then, the *median of these stability values* across all pairs is taken to represent the overall metric stability.
* **Rankings**: Since ranks are discrete, the *Interquartile Range (IQR)* of the rank confidence interval widths is used directly to monitor fluctuations.
* **Global Aggregation**: Element-level stability values are aggregated into a global indicator. For Thresholds and BCa, HERA uses the *arithmetic mean* to ensure global precision. For Rankings, it uses the *maximum* value to enforce worst-case robustness.

## 1. Robust Convergence (Default)

The standard mode for scientific reporting, designed to handle potentially non-monotonic convergence behaviors and avoid premature termination.

* **Smoothing**: HERA applies a trailing *moving average* window (`smoothing_window`) to the global stability curve to filter out stochastic noise and prevent early stops due to random fluctuations.
* **Stability Plateau**: Convergence is declared only when the relative change of the smoothed stability curve remains below the `convergence_tolerance` for a sequence of `convergence_streak_needed` consecutive steps. This ensures that the algorithm has entered a true stable plateau.

## 2. Simple Convergence

A fast heuristic intended for exploratory analysis or when computational speed is prioritized.

* **Mechanism**: Monitors the *raw relative change* of the global stability indicator without smoothing.
* **Stopping Criterion**: Terminates as soon as a *single* relative change between steps falls below the `convergence_tolerance`.
* **Risk**: More susceptible to early fluctuations; a "warm-up" phase (`min_steps_for_convergence_check`) is recommended to ensure reliability.

## 3. Fallback: Elbow Method

Triggered automatically if the robust plateau criterion is not met within the maximum iteration budget ($B_{\text{max}}$).

* **Mechanism**: Identifies the *point of diminishing returns* by calculating the maximum orthogonal distance to the line connecting the first and last points of the normalized stability curve (point of maximum curvature).

> [!WARNING]
> The Ellbow Fallback is a diagnostic aid, not a guarantee!
> Results should be treated as a suggestion to be verified by visual inspection
> and for further analysis (for more details see Troubleshooting Convergence).

---

> [!NOTE]
> The different parameters used in the **Robust Mode** for
> convergence checking have been validated in a empirical Convergence Analysis.
> The default parameters were found to be sufficient for the majority of scenarios.
> However it is not a guarantee!
> Convergence failures using the elbow method were rare (0.8% of cases,
> exclusively in Ranking). Simple Convergence Mode was not assessed in this analysis!
> For detailed results, see
> [Convergence Analysis](https://lerdmann1601.github.io/HERA-Matlab/Convergence_Analysis).

## Troubleshooting Convergence

While the automated check should work for most datasets, "difficult" data with high
variance or flat likelihood landscapes may fail to converge within
`B_end`. In this case, you can try the following:

1. **Check the Elbow**: Inspect the generated stability plots. If you see a
    clear "elbow" where the curve flattens but fluctuates slightly above the
    strict `convergence_tolerance`, the convergence parameters might be too strict
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

5. **Too Early Convergence**: If the automated check finds convergence too
    early (e.g., in robust mode), you can make the parameters stricter (e.g.,
    decrease `convergence_tolerance` or increase `B_step` and/or `B_start`).

> [!IMPORTANT]
> **Reproducibility Note:** Visual inspection of convergence plots is strongly
> recommended for final reporting! All procedures should use the same fixed
> random seed for full reproducibility.
