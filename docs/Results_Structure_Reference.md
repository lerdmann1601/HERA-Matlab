# Results Structure Reference

When running `results = HERA.run_ranking(...)`, the returned structure contains:

| Field | Dimensions | Description |
| --- | --- | --- |
| **Final Results** | | |
| `final_rank` | `[N x 1]` | Final rank for each dataset (1 = Best). |
| `final_order` | `[1 x N]` | Indices of datasets sorted by rank. |
| `final_bootstrap_ranks` | `[N x B]` | Bootstrapped rank distribution for stability analysis. |
| `ci_lower_rank` | `[N x 1]` | Lower bound of Rank Confidence Interval. |
| `ci_upper_rank` | `[N x 1]` | Upper bound of Rank Confidence Interval. |
| `thresholds.d_thresh` | `[1 x M]` | Final determined thresholds for Cliff's Delta. |
| `thresholds.rel_thresh` | `[1 x M]` | Final determined thresholds for Relative Difference. |
| `thresholds.rel_thresh_b` | `[1 x M]` | Raw bootstrap thresholds for Relative Difference (before SEM-cap). |
| `thresholds.min_rel_thresh` | `[1 x M]` | Data-driven minimum thresholds for Relative Difference (SEM-cap). |
| **Statistics** | | |
| `stats.mean` | `[N x M]` | Matrix of mean values per dataset (rows) and metric (columns). |
| `stats.std` | `[N x M]` | Matrix of standard deviations per dataset (rows) and metric (columns). |
| `d_vals_all` | `[Pairs x M]` | Cliff's Delta effect sizes for all pairs/metrics. |
| `rel_vals_all` | `[Pairs x M]` | Relative Mean Differences. |
| `ci_d_all` | `[Pairs x 2 x M]` | BCa Confidence Intervals for Delta. |
| `ci_r_all` | `[Pairs x 2 x M]` | BCa Confidence Intervals for RelDiff. |
| `z0_d_all` | `[Pairs x M]` | BCa bias-correction factors (z0) for Cliff's Delta. |
| `a_d_all` | `[Pairs x M]` | BCa acceleration factors (a) for Cliff's Delta. |
| `z0_r_all` | `[Pairs x M]` | BCa bias-correction factors (z0) for RelDiff. |
| `a_r_all` | `[Pairs x M]` | BCa acceleration factors (a) for RelDiff. |
| `all_p_value_matrices` | `{1 x M}` | Cell array of raw p-values (Wilcoxon). |
| `all_alpha_matrices` | `{1 x M}` | Cell array of Holm-Bonferroni corrected alphas. |
| `all_sig_matrices` | `{1 x M}` | Logical matrices indicating significant wins. |
| **Diagnostics** | | |
| `swap_details.metric1_wins` | `[N x 1]` | Counts of how many significant wins each dataset had on metric 1. (Initial sorting basis). |
| `swap_details.pairwise_swaps_metric1` | `[Pairs x 5]` | Significant wins for Metric 1: [Winner_Idx, Loser_Idx, p_val, d_val, r_val]. |
| `swap_details.pairwise_swaps_metric2` | `[Pairs x 5]` | Significant wins for Metric 2: [Winner_Idx, Loser_Idx, p_val, d_val, r_val]. |
| `swap_details.pairwise_swaps_metric3` | `[Pairs x 5]` | Significant wins for Metric 3: [Winner_Idx, Loser_Idx, p_val, d_val, r_val]. |
| `swap_details.results_metric1` | `{Pairs x 7}` | Raw comparison data for Metric 1: [Winner_Name, Loser_Name, p_val, d_val, r_val, Winner_Idx, Loser_Idx]. |
| `swap_details.results_metric2` | `{Pairs x 7}` | Raw comparison data for Metric 2: [Winner_Name, Loser_Name, p_val, d_val, r_val, Winner_Idx, Loser_Idx]. |
| `swap_details.results_metric3` | `{Pairs x 7}` | Raw comparison data for Metric 3: [Winner_Name, Loser_Name, p_val, d_val, r_val, Winner_Idx, Loser_Idx]. |
| `swap_details.metric2_global_swaps` | `[Swaps x 2]` | Indices [Winner, Loser] of datasets swapped during the iterative sorting of Metric 2 hierarchy. |
| `swap_details.metric3_swaps_a` | `[Swaps x 2]` | Swaps performed via Logic 3A (Tie-break if Metric 1 or 2 was neutral). |
| `swap_details.metric3_swaps_b` | `[Swaps x 2]` | Swaps performed via Logic 3B (Iterative if both previous metrics neutral). |
| `intermediate_orders.after_metric1` | `[1 x N]` | Ranking order determined after processing Metric 1. |
| `intermediate_orders.after_metric2` | `[1 x N]` | Ranking order determined after processing Metric 2. |
| `intermediate_orders.after_metric3` | `[1 x N]` | Ranking order determined after processing Metric 3. |
| `borda_results.rank` | `[N x 1]` | Global Borda count consensus rank (1 = Best). |
| `borda_results.score` | `[N x 1]` | Normalized Borda consensus score (0-100%). |
| `borda_results.rank_distribution` | `{N x 1}` | Cell array containing rank frequency matrices `[Rank, Count]` per dataset. |
| `borda_results.dataset_names` | `{1 x N}` | Dataset names in the same order as provided in the analysis. |
| `power_results.power_matrices` | `{1 x M}` | Cell array containing `[Pairs x 1]` vectors of win probabilities for each metric. |
| `all_permutation_ranks` | `[N x Perms]` | Ranks for every metric permutation tested. |
| `selected_permutations` | `[Perms x M]` | Indices of metrics for each permutation. |
| **Metadata & Configuration** | | |
| `dataset_names` | `{1 x N}` | Cell array of dataset names. |
| `config` | struct | Cleaned configuration used for the analysis (all settings). |
| `meta.timestamp` | string | Global timestamp and unique run-ID. |
| `meta.version` | string | HERA algorithmic version. |
| `meta.n_datasets` | int | Number of datasets / algorithms (`N`). |
| `meta.n_subjects` | int | Total number of data points per dataset (`Clusters`/`Subjects`). |
| `meta.pair_indices` | `[Pairs x 2]` | Exact index comparisons generated matching all data tables. |
| `meta.metric_list` | `{1 x M}` | Names of the evaluated metrics in hierarchical order. |
| `meta.stability_analysis.thresholds.B_vector` | `[1 x Steps]` | Bootstrap sample sizes tested for threshold convergence. |
| `meta.stability_analysis.thresholds.global_stability` | `[1 x Steps]` | Aggregated IQR/Median stability measure per step. |
| `meta.stability_analysis.thresholds.detailed_stability` | `[2 x M x Steps]` | Per-metric stability (IQR/Median) for Cliff's Delta and RelDiff. |
| `meta.stability_analysis.thresholds.converged` | bool | True if threshold estimates safely stabilized. |
| `meta.stability_analysis.thresholds.elbow_indices` | `[1 x Curves]` | Elbow-point indices (fallback if convergence not reached). |
| `meta.stability_analysis.ci.B_vector` | `[1 x Steps]` | Bootstrap sample sizes tested for CI convergence. |
| `meta.stability_analysis.ci.global_stability` | `[1 x Steps]` | Aggregated IQR/Median stability measure per step. |
| `meta.stability_analysis.ci.detailed_stability` | `[2 x M x Steps]` | Per-metric CI width stability for each effect type. |
| `meta.stability_analysis.ci.converged` | bool | True if BCa CI widths safely stabilized. |
| `meta.stability_analysis.ci.elbow_indices` | `[1 x Curves]` | Elbow-point indices (fallback if convergence not reached). |
| `meta.stability_analysis.ranks.B_vector` | `[1 x Steps]` | Bootstrap sample sizes tested for rank convergence. |
| `meta.stability_analysis.ranks.global_stability` | `[1 x Steps]` | Max IQR of rank CI widths per step. |
| `meta.stability_analysis.ranks.detailed_stability` | `[]` | Empty (rank stability is only assessed globally). |
| `meta.stability_analysis.ranks.converged` | bool | True if rank outputs safely stabilized. |
| `meta.stability_analysis.ranks.elbow_indices` | `[1 x Curves]` | Elbow-point indices (fallback if convergence not reached). |
| `meta.bootstrap_B.thresholds` | int | Final determined B used for Thresholds computation. |
| `meta.bootstrap_B.ci` | int | Final determined B used for BCa CI computation. |
| `meta.bootstrap_B.ranks` | int | Final determined B used for Rank computation. |

## Definitions & Glossary

| Term | Full Name | Description |
| --- | --- | --- |
| `p_val` | p-value | Probability of observing the results by chance (Wilcoxon Signed-Rank Test). |
| `d_val` | Cliff's Delta | Non-parametric effect size measuring the probability of one dataset outperforming another. Range: [-1, 1]. |
| `r_val` | RelDiff | Percentage-based relative difference of the means between two datasets. |
| `sig` | Significant | Logical `true` if `p_val <= alpha` (Holm-Bonferroni corrected). |
| `rel` | Relevant | Logical `true` if `abs(d_val) >= d_thresh` AND `r_val >= rel_thresh`. |

## Ranking Logic & Hierarchy

The final ranking is achieved through a multi-stage sequential process:

1. **Initial Ranking (Metric 1)**: Datasets are sorted using a 3-step tie-break logic:
    * **Step 1: Win Count**: Datasets are primarily sorted by the number of significant and relevant wins.
    * **Step 2: Cliff's Delta (d)**: If Win Counts are equal, the pairwise stochastic dominance (`d_val`) between the tied datasets decides.
    * **Step 3: Mean Value**: If `d_val` is neutral (within epsilon), the raw mean value of Metric 1 serves as the final tie-breaker.
2. **Global Correction (Metric 2)**: The ranking from Step 1 is iteratively adjusted. If a lower-ranked dataset shows a significant and relevant win over a higher-ranked one according to Metric 2, they are swapped. This correction takes precedence over Metric 1 results treating Metric 2 as a set of non-negotiable concerns (e.g., safety, fundamental accuracy) that must be satisfied regardless of primary performance.
3. **Targeted Correction (Metric 3)**: Metric 3 (using Logic 3A and 3B) specifically adjusts pairs that were "neutral" in the preceding hierarchy:
    * **Logic 3A**: Swaps adjacent datasets if Metric 2 was neutral but Metric 3 shows a significant win.
    * **Logic 3B**: Swaps adjacent datasets if both Metric 1 AND 2 were neutral but Metric 3 shows a significant win.

For more details on the possible ranking logics, see [Ranking Modes Explained](https://lerdmann1601.github.io/HERA-Matlab/Ranking_Modes_Explained).
