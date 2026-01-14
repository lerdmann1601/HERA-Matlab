# Methodological Guidelines & Limitations

The statistical rigor of HERA (e.g., Holm-Bonferroni correction, Bootstrapping)
imposes practical limits on the number of datasets (*N*) and sample size (*n*).
Therefore the following guidelines are provided as theoretical considerations
but should not be taken as strict requirements.

## Number of Datasets (*N*)

Increasing *N* quadratically increases the number of pairwise comparisons (*m* =
*N*(*N*-1)/2), which reduces statistical power due to strict corrections.

* **Minimum (*N* = 3)**: Required for a meaningful ranking. (*N* = 2 is just a
    simple comparison).
* **Optimal (*N* ≈ 8–10)**: Balances ranking depth with statistical
    power (28–45 comparisons).
* **Upper Limit (*N* ≈ 15)**: Not generally recommended. The loss of power
    from FWER corrections makes detecting true differences unlikely. However,
    it is possible to use HERA with *N* > 15 and you can just give it a try.

> **Visual Limit (*N* ≤ 20)**: While HERA should technically compute rankings
> for any *N* (exported to CSV/JSON), the generated plots (e.g. Win-Loss Matrix,
> Final Summary) visually degrade beyond *N* = 20. For *N* > 20, I recommend
> relying on the machine-readable and text-based outputs. You can disable plots
> (`create_reports: false`) to save runtime.
>
> **Recommendation:**
> If you have a large pool of candidates (*N* >> 15), it could be a good idea to
> apply a global screening method (e.g., **Friedman Test** followed by Nemenyi post-hoc)
> to identify the top tier of algorithms. Ranking the entire set with HERA may be
> overly strict; instead, select the top performing subset (e.g., the best 10-15)
> and use HERA for the final ranking of the best candidates.

## Sample Size (*n*)

A balance between statistical stability and computational feasibility is
required.

* **Minimum (*n* = 16)**: Required for the Wilcoxon test to use the Normal
    Approximation in Matlab.
* **Robust Min (*n* ≈ 25–30)**: Necessary for stable BCa confidence
    intervals and Jackknife estimates (Although it automatically switches
    to Percentil Bootstrap if Bias or Jackknife estimates become unstable).
* **Optimal (*n* ≈ 50–300)**: Best balance of power, stability, and
    runtime.
* **Upper Limit (*n* ≈ 1,000–5,000)**: Higher *n* improves statistics
    but linearly scales runtime. *n* ≫ 5,000 may be computationally
    impractical due to extensive bootstrapping.

> **Recommendation:** Perform an *a priori* power analysis to estimate the
> required *n* for your chosen *N*.

## Missing Data Handling (NaN)

HERA is robust against missing data (`NaN`) but handling it comes with trade-offs:

* **Pairwise Deletion**: HERA employs pairwise deletion to maximize data
    usage without requiring complex imputation. While this assumes data is missing
    completely at random (MCAR), it remains methodologically robust: By relying
    on discrete, independent pairwise comparisons, the algorithm avoids the
    mathematical inconsistencies (e.g., non-positive definite matrices) that
    typically compromise pairwise deletion in global multivariate statistics.
* **Robust Bootstrapping**: If `NaN`s are detected, HERA automatically switches
    to a "Robust Path". This dynamically filters invalid data for each bootstrap
    resample to ensure correctness, which **significantly increases runtime**
    especially for large sample sizes (*n*).
* **Automatic Warning**: A warning is issued if valid data drops below 80% for
    any comparison however it is not a strict requirement.

> **Recommendation**: Minimize `NaN`s to preserve statistical power and
> performance. For critical analyses with substantial data loss, use
> validated imputation methods (e.g., MICE) *before* running HERA.

## Number of Metrics (*M* ≤ 3)

HERA is designed for a maximum of 3 hierarchical metrics to maintain methodological
robustness. This limit is inherent to the hierarchical-compensatory design and is
based on the following methodological considerations:

* **Loss of Interpretability**: With every additional hierarchical level, the causal
    chain of the ranking decision becomes opaque and increasingly difficult to trace.
    Limiting the depth to 3 levels ensures that the ranking logic remains transparent
    and empirically justifiable.
* **Increased Risk of Collinearity**: Adding more metrics increases the probability
    of introducing redundant criteria (e.g., two metrics measuring valid features
    of the same underlying property). In a sequential logic, these correlates
    would be falsely treated as independent evidence, distorting the ranking.
* **Functional Saturation**: The hierarchical-compensatory logic is fully saturated
    by three levels (Sorting, Correction, Finalization). Adding a fourth metric
    yields diminishing margins of utility, as the probability of meaningful rank
    adjustments approaches zero, while the complexity of the decision model
    increases disproportionately.

> **Recommendation:**
> If you want to consider more than 3 metrics and use HERA you could first
> perform a check for collinearity (e.g., using a correlation matrix).
> Strongly correlated metrics could be aggregated into a common factor
> (e.g., via Principal Component Analysis (PCA)) before running the HERA
> analysis with up to 3 metrics.
>
> If your study design requires the **simultaneous integration** of a large number
> of metrics ($M \gg 3$) HERA is not feasible and compensatory or outranking MCDA
> methods are methodologically more appropriate. In this case, approaches like
> **TOPSIS** or **PROMETHEE** might be a better choice.
