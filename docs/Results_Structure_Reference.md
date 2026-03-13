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
| `swap_details.metric1_wins` | `[N x 1]` | Counts of how many significant wins each dataset had on metric 1. (Global sorting is only performed for the primary metric). |
| `swap_details.pairwise_swaps_metricX` | `[Pairs x 5]` | Details of all significant/relevant pairwise wins for metric X (`X` in 1..3). |
| `swap_details.results_metricX` | `{Pairs x 7}` | Cell array of raw comparison data for metric X (`X` in 1..3). Contains: [Winner_Name, Loser_Name, p_val, d_val, r_val, Winner_Idx, Loser_Idx]. |
| `swap_details.metric2_global_swaps` | `[Swaps x 2]` | Indices of datasets swapped during the iterative sorting of the Metric 2 hierarchy. |
| `swap_details.metric3_swaps_a` | `[Swaps x 2]` | Swaps performed via Logic 3A (Tie-break if Metric 1 or 2 was neutral). |
| `swap_details.metric3_swaps_b` | `[Swaps x 2]` | Swaps performed via Logic 3B (Iterative if both previous metrics neutral). |
| `intermediate_orders.after_metricX` | `[1 x N]` | Ranking order determined after processing metric X (`X` in 1..3). |
| `borda_results.rank` | `[N x 1]` | Global Borda count consensus rank (1 = Best). |
| `borda_results.score` | `[N x 1]` | Normalized Borda consensus score (0-100%). |
| `borda_results.rank_distribution` | `{N x 1}` | Cell array containing rank frequency matrices `[Rank, Count]` per dataset. |
| `borda_results.dataset_names` | `{1 x N}` | Dataset names in the same order as the Borda results. |
| `power_results.power_matrices` | `{1 x M}` | Cell array containing `[N x N]` matrices of win probabilities (Significant AND Relevant) for each metric. |
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
