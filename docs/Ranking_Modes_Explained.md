# Ranking Modes Explained

HERA (Hierarchical-Compensatory, Effect-Size-Driven Ranking Algorithm) implements a structured, mostly non-parametric framework to resolve trade-offs between up to three user-prioritized metrics. Unlike traditional aggregation methods that often rely on subjective weighting, HERA utilizes a sequential ranking logic. This logic combines statistical significance and effect sizes with data-driven thresholds for practical equivalence in which candidates are considered neutral if they do not show a significant and relevant difference.

For more information please refer to the [Theoretical Background](https://lerdmann1601.github.io/HERA-Matlab/Methodology) and [Methodological Guidelines & Limitations](https://lerdmann1601.github.io/HERA-Matlab/Methodological_Guidelines_&_Limitations) sections.

## The Core Philosophy: Hierarchical-Compensatory

The algorithm is structured as a multi-stage tournament that distinguishes between **primary baseline performance** (Metric 1), **compensatory rank correction and strict veto constraints** (Metric 2), and **local adjustment and indifference resolution** (Metric 3).

This allows HERA to resolve trade-offs without collapsing multi-dimensional performance into a single weighted scalar, maintaining full transparency about *why* a candidate was ranked in its position.

---

## 1. The "Win" Logic: Dual-Criteria Assessment

A candidate $A$ only "wins" over candidate $B$ if the difference satisfies three conjunctive criteria. If any criterion is not met, the comparison is considered *neutral*. Here, Cliff's Delta ($d$) and the Relative Mean Difference ($\text{RelDiff}$) complement Statistical Significance ($p$) to provide a dual-criteria assessment of **Stochastic Dominance** and **Practical Magnitude**. While the Wilcoxon signed-rank test evaluates within-pair differences, Cliff's Delta uses the independent-samples formulation across all cross-pairs as a conservative, orthogonal filter. This prevents "wins" driven by isolated outliers or consistent but practically trivial differences, ensuring that every win reflects both within-pair statistical significance and global distributional superiority.

1. **Statistical Significance**: $p \leq \alpha_{\text{Holm}}$ (paired Wilcoxon signed-rank test with step-down Holm-Bonferroni correction).
2. **Stochastic Dominance**: $|d| \geq \theta_d$ (Cliff's Delta, $d = P(X > Y) - P(Y > X)$, must meet or exceed the data-driven threshold).
3. **Practical Magnitude**: $\text{RelDiff} \geq \theta_{\text{RelDiff}}$ (Relative Mean Difference must meet or exceed the threshold).
   The Relative Mean Difference ($\text{RelDiff}$) is defined for group means $\bar{x}$ and $\bar{y}$ as:
   $\text{RelDiff} = \frac{\vert\bar{x} - \bar{y}\vert}{\vert\frac{1}{2}(\bar{x} + \bar{y})\vert}$

### Data-Driven Threshold Calculation

To eliminate subjective parameterization, HERA derives the thresholds $\theta_d$ and $\theta_{\text{RelDiff}}$ directly from the data using a *Percentile Bootstrap* procedure. Following the principle of *Empirical Benchmarks*, where the data itself provides the standard for the plausibility and magnitude of differences, the algorithm identifies a distribution-based *Minimal Important Difference (MID)*:

1. **Bootstrap Medians**: For each metric, HERA performs bootstrap iterations on the distribution of pairwise effect sizes and calculates the median for each iteration.
2. **Threshold Selection**: The threshold is defined as the lower $\alpha/2$-quantile (e.g., the 2.5th percentile) of these bootstrap medians. This identifies the lower bound of what is considered a "typical" effect in the dataset.

### SEM-Based Lower Bound for Relative Difference

For cases where the bootstrapped threshold for $\text{RelDiff}$ ($\theta_{\text{RelDiff},b}$) becomes systematically too low (e.g., in homogeneous or low-variance datasets), potentially allowing "wins" based on trivial differences within measurement noise, HERA enforces a dynamic lower bound ($\theta_{\text{SEM}}$) inspired by the *Smallest Worthwhile Effect* concept. This minimum threshold is derived from the two-sided critical value of the $t$-distribution ($df = n_{\text{median}} - 1$) multiplied by the median Standard Error of the Mean (SEM) across all candidates, normalized by the absolute grand mean of the metric. The final $\text{RelDiff}$ threshold is the maximum of the empirical bootstrap estimate and this SEM-based lower bound: $\theta_{\text{RelDiff}} = \max(\theta_{\text{RelDiff},b}, \theta_{\text{SEM}})$.

> [!WARNING]
> **Safety Mechanism**: HERA implements specific numerical safeguards to handle scenarios where mathematical singularities could occur:
>
> - **RelDiff Effect Size**: To prevent division by zero when the sum of the means of a compared pair is zero ($\bar{x} + \bar{y} = 0$, e.g., in z-transformed or zero-centered data), $\text{RelDiff}$ is conservatively set to $0$ for that pair, accepting a Type II error (preventing a win unless $\theta_{\text{RelDiff}} = 0$) to preserve metric integrity.
> - **SEM-based Threshold**: The lower bound $\theta_{\text{SEM}}$ is explicitly deactivated (set to $0$) if the grand mean across all candidates is zero ($\bar{S}_{\text{total}} = 0$), and it vanishes naturally if the median of the candidate SEMs is zero (when the majority of candidates exhibit zero variance). In these cases, the algorithm conservatively falls back to the empirical bootstrap threshold ($\theta_{\text{RelDiff},b}$) to dictate the minimum relevance criterion for $\text{RelDiff}$ (and if all pairwise $\text{RelDiff}$ values are zero, the decision relies on $\theta_d$ and statistical significance).

---

## 2. Sequential Ranking Stages

### Stage 1: Primary Baseline Ranking (Metric 1)

Candidates are initially sorted in descending order by their *Win Count* ($w_i^{(1)}$) in Metric 1. Ties between adjacent candidates are resolved deterministically via a cascading hierarchy:

- **Tie-Break A (Stochastic Dominance)**: If Win Counts are equal, the pairwise Cliff's Delta ($d$) between the adjacent tied candidates decides (requires $|d| > \epsilon$, with numerical tolerance $\epsilon = 10^{-9}$).
- **Tie-Break B (Numerical Mean)**: If $|d| \leq \epsilon$, the raw mean value of Metric 1 serves as the final arbiter.

### Stage 2: Compensatory Correction (Metric 2)

The baseline ranking from Stage 1 is iteratively adjusted through global pairwise comparisons with all higher-ranked candidates. If a lower-ranked candidate demonstrates a statistically significant and practically relevant win over a higher-ranked one in Metric 2, their positions are swapped.

Stage 2 unifies compensatory and strict veto approaches: Metric 2 acts as an independent compensatory rank correction that overrides the primary ordering whenever a robust superiority in Metric 2 is established—no matter how well a candidate performed in Metric 1, it cannot retain its higher rank if it is significantly and relevantly outperformed in a critical secondary metric (e.g., safety or fundamental accuracy).

### Stage 3: Local Adjustment and Refinement (Metric 3)

Stage 3 extends beyond simple tie-breaking to resolve remaining local indifferences where Metric 1 or Metric 2 did not establish a robust difference, using two sequential sub-logics:

- **Sublogic 1 / Logic 3A (Constrained Local Correction)**: Performs a single-pass swap of adjacent candidates if the preceding comparison was neutral (Metric 2 in `M1_M2_M3` mode, or Metric 1 in `M1_M3A` mode) and Metric 3 shows a statistically significant and practically relevant win. In `M1_M2_M3` mode, this provides a bounded local correction without triggering cascading global shifts; in `M1_M3A` mode, it functions as a lexicographic tie-breaker for Metric 1.
- **Sublogic 2 / Logic 3B (Indifference Resolution)**: Performs iterative swaps of adjacent candidates within neutral clusters where *both* Metric 1 AND Metric 2 were neutral, utilizing Metric 3 wins until a stable local ordering is achieved.

### Handling Intransitivity (Cycle Detection & Fallback)

Pairwise comparisons across multiple criteria or skewed distributions can theoretically exhibit cyclic preferences ($A \succ B \succ C \succ A$). HERA addresses intransitivity through a two-layered defense strategy:

1. **Intrinsic Cycle Suppression**: Requiring every win to exceed both data-driven thresholds ($\theta_d, \theta_{\text{RelDiff}}$) creates an indifference zone for trivial effects, while $\text{RelDiff}$ operates on a transitive scale of group means. This naturally suppresses spurious, noise-induced cycles before explicit detection is required.
2. **Deterministic Cycle Detection & Rollback**: While Stage 1 and Stage 3 Sublogic 1 (3A) are inherently cycle-free by design, the iterative sorting loops in Stage 2 and Stage 3 Sublogic 2 (3B) enforce a strict iteration cap ($I_{\max} = N^2$, with a numerical safety margin). If an iterative loop does not converge within this limit, the step is flagged as intransitive (`cycle_m2` or `cycle_m3b`), the loop is aborted, and the ranking deterministically reverts (rolls back) to the validated, stable ordering of the preceding hierarchy stage.

---

## 3. Operational Modes & Use Cases

| Mode | Stage Progression | Logical Focus | Recommended Use Case |
| :--- | :--- | :--- | :--- |
| **`M1`** | Stage 1 only | **Single-Metric** | Evaluating a single primary performance metric. |
| **`M1_M2`** | Stage 1 $\rightarrow$ 2 | **Compensatory** | Balancing primary performance against a secondary compensatory or strict veto criterion. |
| **`M1_M3A`** | Stage 1 $\rightarrow$ 3A | **Lexicographic** | Resolving adjacent primary-metric indifferences (ties) without overriding established primary wins. |
| **`M1_M2_M3`** | Stage 1 $\rightarrow$ 2 $\rightarrow$ 3 (3A & 3B) | **Full Hierarchy** | Multi-criteria sequential benchmarking where $M_1$ sorts, $M_2$ compensates/corrects, and $M_3$ locally adjusts and refines. |

---

## 4. Sensitivity & Stability Analysis

To ensure the ranking is not an artifact of a specific hierarchy configuration or point estimate, HERA provides three complementary analyses to assess robustness, stability, and statistical power:

- **Borda Count (Sensitivity Analysis)**: Re-runs the ranking across all $M!$ permutations of the metric hierarchy and aggregates the outcomes into a normalized Borda consensus score (0–100%). This identifies candidates that are robustly superior regardless of metric prioritization and quantifies how sensitive a candidate's rank is to the chosen hierarchy.
- **Cluster Bootstrap (Rank Stability)**: Resamples subjects as complete clusters with replacement and repeats the full hierarchical ranking in each bootstrap sample to compute a *95% Rank Confidence Interval* via the percentile method, quantifying the stability of each rank position under sampling variability.
- **Empirical Power Analysis**: Estimates the empirical probability (relative frequency) of detecting a statistically significant and practically relevant *win* for every pairwise comparison under the observed data distribution using cluster bootstrap resampling ($B = 10{,}000$ iterations by default), contrasting with classical *post-hoc* power analyses that assume the observed effect equals the true population effect. For an observed *win*, high empirical power (e.g., $> 75\%$) confirms that superiority is robust and reproducible under resampling, whereas lower power indicates a marginally reliable result. Conversely, for an observed *neutral* result, low empirical power confirms that a win is rarely reproduced under resampling, whereas intermediate power (e.g., $> 25\%$) highlights borderline cases where data variance prevents a definitive decision.
