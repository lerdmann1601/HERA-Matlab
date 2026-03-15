# Ranking Modes Explained

HERA (Hierarchical-Compensatory, Effect-Size-Driven Ranking Algorithm) implements a structured, mostly non-parametric framework to resolve trade-offs between up to three user-prioritized metrics. Unlike traditional aggregation methods that often rely on subjective weighting, HERA utilizes a sequential ranking logic. This logic combines statistical significance and effect sizes with data-driven thresholds for practical equivalence in which candidates are considered neutral if they do not show a significant and relevant difference.

For more information please refer to the [Theoretical Background](https://lerdmann1601.github.io/HERA-Matlab/Theoretical_Background) and [Methodological Guidelines](https://lerdmann1601.github.io/HERA-Matlab/Methodological_Guidelines) sections.

## The Core Philosophy: Hierarchical-Compensatory

The algorithm behaves like a multi-stage tournament. It distinguishes between *primary performance* (Metric 1), *non-negotiable concerns* such as safety or fundamental accuracy (Metric 2), and *local tie-breakers* (Metric 3).

This allows HERA to resolve trade-offs without collapsing multi-dimensional performance into a single scalar, maintaining full transparency about *why* a method was ranked in its position.

---

## 1. The "Win" Logic: Dual-Criteria Assessment

A candidate $A$ only "wins" over candidate $B$ if the difference satisfies three conjunctive criteria. If any criterion is not met, the comparison is considered *neutral*. Here Cliff's Delta ($d$) and Relative Mean Difference ($\text{RelDiff}$) are used as a complementary measure to the Statistical Significance ($p$) to provide a dual-criteria assessment of Dominance and Magnitude. This logic prevents "wins" based on outlier driven or consistent but irrelevant differences and ensures that a "win" is both significant and relevant.

1. **Statistical Significance**: $p < \alpha_{\text{Holm}}$ (Wilcoxon signed-rank test with Holm-Bonferroni correction).
2. **Stochastic Dominance**: $|d| \geq \theta_d$ (Cliff's Delta must exceed the data-driven threshold).
3. **Practical Magnitude**: $\text{RelDiff} \geq \theta_{\text{RelDiff}}$ (Relative Mean Difference must exceed the threshold).
   The Relative Mean Difference (RelDiff) is defined for group means $\bar{x}$ and $\bar{y}$ as:
   $$\text{RelDiff} = \frac{|\bar{x} - \bar{y}|}{|\frac{1}{2}(\bar{x} + \bar{y})|}$$

### Data-Driven Threshold Calculation

To eliminate subjective parameterization, HERA derives the thresholds $\theta_d$ and $\theta_{\text{RelDiff}}$ directly from the data using a *Percentile Bootstrap* procedure. Following the principle of *Empirical Benchmarks*, where the data itself provides the standard for the plausibility and magnitude of differences, the algorithm identifies a distribution-based *Minimal Important Difference (MID)*:

1. **Bootstrap Medians**: For each metric, HERA performs bootstrap iterations on the distribution of pairwise effect sizes and calculates the median for each iteration.
2. **Threshold Selection**: The threshold is defined as the lower $\alpha/2$-quantile (e.g., the 2.5th percentile) of these bootstrap medians. This identifies the lower bound of what is considered a "typical" effect in the dataset.

### SEM-Based Lower Bound for Relative Difference

For cases where the bootstrapped thresholds for the *RelDiff* becomes systematically too low, potentially allowing "wins" based on trivial differences, HERA enforces a dynamic lower bound inspired by the *Smallest Worthwhile Effect* concept. The minimum threshold is derived from the critical value of the $t$-distribution multiplied by the median Standard Error of the Mean (SEM) across all candidates, normalized by the grand mean of the metric. The final RelDiff threshold is the maximum of the empirical bootstrap estimate and this SEM-based lower bound.

> [!WARNING]
> **Safety Mechanism**: HERA implements specific numerical safeguards to handle scenarios where mathematical singularities could occur:
>
> - **RelDiff Effect Size**: The calculation is deactivated (set to 0) if the sum of the means of the compared pair is
> zero (e.g., in centered data), as the relative magnitude loses its reference point.
> - **SEM-based Threshold**: The lower bound is deactivated (set to 0) if the grand mean across all candidates is zero,
> or if the median of the data exhibits zero variance (degenerate data).
> In these cases, the algorithm conservatively falls back to the empirical bootstrap threshold ($\theta_d$) and
> significance testing to ensure a stable and deterministic result.

---

## 2. Sequential Ranking Stages

### Stage 1: Initial Tournament (Metric 1)

Datasets are primarily sorted by their *Win Count* in Metric 1.

- **Tie-Break A**: If Win Counts are equal, the pairwise Cliff's Delta ($d$) between the tied datasets decides.
- **Tie-Break B**: If $|d| < \epsilon$, the raw mean value of Metric 1 serves as the final arbiter.

### Stage 2: Global Compensatory Correction (Metric 2)

The ranking from Stage 1 is iteratively adjusted. If a lower-ranked dataset shows a significant and relevant win over a higher-ranked one in Metric 2, they are swapped.

This acts as a "veto" logic—no matter how well a method performs in Metric 1, it cannot hold its rank if it fails significantly against a competitor in a critical Metric 2 (e.g., a "Safety" concern).

### Stage 3: Local Refinement (Metric 3)

Resolves remaining "neutral" clusters where Metric 1 or 2 could not establish a robust difference:

- **Logic 3A**: Swaps adjacent datasets if Metric 2 was neutral but Metric 3 shows a clear win. In the full hierarchy mode (`M1_M2_M3`), this serves as a local secondary correction for subordinate constraints. In `M1_M3A` mode, however, it functions as a classic tie-breaker for Metric 1.
- **Logic 3B**: Swaps adjacent datasets within a cluster if *both* Metric 1 AND 2 were neutral to achieve a final ordering.

---

## 3. Operational Modes & Use Cases

| Mode | Stage Progression | Logical Focus | Recommended Use Case |
| :--- | :--- | :--- | :--- |
| **`M1`** | Stage 1 only | Single-objective | Evaluating a single primary KPI. |
| **`M1_M2`** | Stage 1 $\rightarrow$ 2 | **Correction** | Balancing performance vs. "Hard" constraints (Safety, Veto-criteria). |
| **`M1_M3A`** | Stage 1 $\rightarrow$ 3A | **Tie-Breaking** | resolving ties without overriding primary results. |
| **`M1_M2_M3`** | Stage 1 $\rightarrow$ 2 $\rightarrow$ 3 | **Full Logic** | Multi-objective sequential benchmarking where $M_1$ sorts, $M_2$ corrects and $M_3$ refines. |

---

## 4. Sensitivity & Stability Analysis

To ensure the ranking is not an artifact of a specific configuration or point estimate, HERA provides several steps to assess the sensitivity and stability of the ranking.

- **Borda Count**: Re-runs the analysis for all permutations of the metric hierarchy to help identify candidates that are robustly superior regardless of prioritization and to asses the sensitivity of a ranked candidate to the choice of metric hierarchy.
- **Cluster Bootstrap**: Resamples participants as whole clusters and performs the whole ranking process in each bootstrap sample to generate a *95% Rank Confidence Interval* with the percentile method, to assess the stability of each position and the relative frequency of each candidate for a specific rank under resampling.
- **Empirical Power Analysis**: Estimates the possibility of replicating a "win", "loss", or "neutral" result through simulation (cluster bootstrap) given the data characteristics. High empirical power (>80%) combined with a neutral result suggests that the true difference is likely negligible, whereas low power indicates that the dataset may be too small or noisy to draw robust conclusions.
