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
| `thresholds` | struct | Thresholds for Cliff's Delta and RelDiff. |
| **Statistics** | | |
| `d_vals_all` | `[Pairs x M]` | Cliff's Delta effect sizes for all pairs/metrics. |
| `rel_vals_all` | `[Pairs x M]` | Relative Mean Differences. |
| `ci_d_all` | `[Pairs x 2 x M]` | BCa Confidence Intervals for Delta. |
| `ci_r_all` | `[Pairs x 2 x M]` | BCa Confidence Intervals for RelDiff. |
| `all_p_value_matrices` | `{1 x M}` | Cell array of raw p-values (Wilcoxon). |
| `all_alpha_matrices` | `{1 x M}` | Cell array of Holm-Bonferroni corrected alphas. |
| `all_sig_matrices` | `{1 x M}` | Logical matrices indicating significant wins. |
| **Diagnostics** | | |
| `swap_details` | struct | Log of logic-based rank swaps (M1 vs M2). |
| `intermediate_orders` | struct | Rankings after each hierarchical stage. |
| `borda_results` | struct | Consensus ranking from sensitivity analysis. |
| `power_results` | struct | Post-hoc power analysis data. |
| `all_permutation_ranks` | `[N x Perms]` | Ranks for every metric permutation tested. |
| `selected_permutations` | `[Perms x M]` | Indices of metrics for each permutation. |
