function [final_order, final_rank, all_sig_matrices, all_alpha_matrices, all_p_value_matrices, swap_details, intermediate_orders] = ...
            calculate_ranking(all_data, effect_sizes, thresholds, config, dataset_names, pair_idx_all, subset_indices)
% CALCULATE_RANKING - Calculates the ranking of datasets based on a hierarchical metric system.
%
% Syntax:
%   [final_order, final_rank, all_sig_matrices, all_alpha_matrices, all_p_value_matrices, swap_details, intermediate_orders] = ...
%       calculate_ranking(all_data, effect_sizes, thresholds, config, dataset_names, pair_idx_all, subset_indices)
%
% Description:
%   This function implements a multi-stage ranking algorithm to determine the order of datasets.
%   It dynamically adapts the logic based on the number of metrics (1, 2, or 3) and the 'ranking_mode' provided in the config.
%
% Workflow:
%   1. Pre-calculation: 
%      Performs statistical tests (Wilcoxon, Holm-Bonferroni) and effect size checks (Cliff's Delta, RelDiff) for all pairs across all loaded metrics (1, 2, or 3).
%   2. Initial Ranking (Metric 1): 
%      Creates an initial sorted order based on the number of significant and relevant "wins" according to Metric 1. This is always performed.
%   3. Rank Correction (Metric 2 - Optional): 
%      If ranking_mode is 'M1_M2' or 'M1_M2_M3', iteratively adjusts the Metric 1 ranking by swapping any pair where a lower-ranked dataset beats a higher-ranked one according to Metric 2.
%   4. Rank Correction (Metric 3 / Tie-break - Optional): 
%      - If 'M1_M3A': Performs a final tie-break using Metric 2's data based on Metric 1 neutrality.
%      - If 'M1_M2_M3': Performs a final adjustment based on Metric 3, prioritizing pairs that were "neutral" in previous metrics (Logic 3A and 3B).
%   5. Finalization: 
%      Generates the final rank vector and bundles all intermediate results and swap details for logging.
%
% Ranking Logic:
%   - 'M1':       Initial Sorting (Metric 1) only.
%   - 'M1_M2':    Initial Sorting (Metric 1) -> Iterative Correction (Metric 2).
%   - 'M1_M3A':   Initial Sorting (Metric 1) -> Single Tie-Break (Metric 2, using Logic 3A rules).
%   - 'M1_M2_M3': Initial Sorting (Metric 1) -> Iterative Correction (Metric 2) -> Full Correction (Metric 3, Logic 3A+3B).
%
% Inputs:
%   all_data             - Cell array of data matrices in hierarchical order (1, 2, or 3 cells).
%   effect_sizes         - Struct with effect sizes (d_vals_all, rel_vals_all).
%   thresholds           - Struct with thresholds (d_thresh, rel_thresh).
%   config               - Struct with general configurations (must include 'ranking_mode').
%   dataset_names        - Cell array with the names of the datasets.
%   pair_idx_all         - Indices of all pairwise comparisons.
%   subset_indices       - (Optional) Vector of row indices to use from all_data (virtual view).
%
% Outputs:
%   final_order          - Vector with the final order of dataset indices.
%   final_rank           - Vector assigning each dataset its final rank.
%   all_sig_matrices     - Cell array of significance matrices for each metric.
%   all_alpha_matrices   - Cell array of Holm-Bonferroni corrected alpha values.
%   all_p_value_matrices - Cell array of raw p-values for each comparison.
%   swap_details         - Struct with detailed log information about the swap operations.
%   intermediate_orders  - Struct with the rankings after each stage.
% 
% Author: Lukas von Erdmannsdorff

%% 1. Pre-calculation: Significance for all pairs and metrics
arguments
    all_data (1,:) cell
    effect_sizes (1,1) struct
    thresholds (1,1) struct
    config (1,1) struct
    dataset_names (1,:) cell
    pair_idx_all
    subset_indices = []
end

% Initialization of basic parameters.
num_metrics = numel(all_data);
m = size(pair_idx_all, 1);
num_datasets = size(all_data{1}, 2);
ranking_mode = config.ranking_mode;

% Initialization of cell arrays to store the results for each metric.
all_temp_results = cell(1, num_metrics);
all_sig_matrices = cell(1, num_metrics);
all_alpha_matrices = cell(1, num_metrics);
all_p_value_matrices = cell(1, num_metrics);
all_pairwise_swaps = cell(1, num_metrics);

% Loop over each metric to perform the statistical tests.
for metric_idx = 1:num_metrics
    % Extracts the relevant data and thresholds for the current metric.
    data = all_data{metric_idx};
    d_vals = effect_sizes.d_vals_all(:, metric_idx);
    rel_vals = effect_sizes.rel_vals_all(:, metric_idx);
    d_thr = thresholds.d_thresh(metric_idx);
    r_thr = thresholds.rel_thresh(metric_idx);
    alpha = config.alphas(metric_idx);
    
    temp_m = cell(m, 7);
    p_values = zeros(m, 1);
    
    % Calculates p-values (Wilcoxon signed-rank test) for all pairs.
    for k = 1:m
        i = pair_idx_all(k, 1);
        j = pair_idx_all(k, 2);
        if isempty(subset_indices)
            % Standard mode: use full columns
            col_i = data(:, i);
            col_j = data(:, j);
        else
            % Optimized mode: use virtual subset
            col_i = data(subset_indices, i);
            col_j = data(subset_indices, j);
        end
        
        % Find valid rows for this pair only (pairwise exclusion)
        valid_rows = ~isnan(col_i) & ~isnan(col_j);
        % Performs the test for paired samples if enough valid pairs are present.
        if sum(valid_rows) >= 2 % signrank needs at least 2 valid pairs
            p = signrank(col_i(valid_rows), col_j(valid_rows));
        else
            p = NaN; % Indicates that the test could not be performed
        end
        p_values(k) = p;
        if nargout >= 6
            % Stores the raw results of the comparison for later logging.
            temp_m(k, :) = {dataset_names{i}, dataset_names{j}, p, d_vals(k), rel_vals(k), i, j};
        end
    end
    all_temp_results{metric_idx} = temp_m;
    
    % Performs the Holm-Bonferroni correction to control the family-wise error rate.
    [sorted_p, sort_idx] = sort(p_values); % Sort p-values in ascending order.
    % Checks for each p-value if it meets the Holm condition using centralized logic.
    [is_sig_orig, ~] = HERA.stats.holm_bonferroni(p_values, alpha);
    % Map to sorted order for the subsequent loop structure which iterates by rank k
    holm_hypothesis = is_sig_orig(sort_idx);
    
    % Initialization of the result matrices for the current metric.
    sig_matrix = false(num_datasets);      % Stores significant "wins" (directional: sig_matrix(i,j)=true -> i beats j).
    alpha_matrix = NaN(num_datasets);      % Stores the corrected alpha values for each comparison.
    p_value_matrix = NaN(num_datasets);    % Stores the raw p-values for each comparison.
    pairwise_swaps = zeros(m, 5);          % Logs all significant comparisons for output.
    swap_count = 0;
    
    % Fills the result matrices based on the Holm-Bonferroni correction.
    for k = 1:m
        % Stores the corrected alpha value for each comparison.
        i_alpha = pair_idx_all(sort_idx(k), 1);
        j_alpha = pair_idx_all(sort_idx(k), 2);
        alpha_k = alpha / (m - k + 1);
        alpha_matrix(i_alpha, j_alpha) = alpha_k;
        alpha_matrix(j_alpha, i_alpha) = alpha_k;
        
        % Stores the raw p-value for each comparison.
        orig_idx_p = pair_idx_all(k, :);
        p_value_matrix(orig_idx_p(1), orig_idx_p(2)) = p_values(k);
        p_value_matrix(orig_idx_p(2), orig_idx_p(1)) = p_values(k);
        
        % Checks if a comparison is significant after Holm-Bonferroni.
        if holm_hypothesis(k)
            orig_idx_sorted = sort_idx(k);
            i = pair_idx_all(orig_idx_sorted, 1);
            j = pair_idx_all(orig_idx_sorted, 2);
            d_val = d_vals(orig_idx_sorted);
            r_val = rel_vals(orig_idx_sorted);
            
            % Additional check if the effect sizes exceed the thresholds (practical relevance).
            % Fix: Added epsilon (1e-9) to handle floating point issues where d_val == d_thr
            if abs(d_val) >= (d_thr - 1e-9) && r_val >= (r_thr - 1e-9)
                swap_count = swap_count + 1;
        if d_val > 0 % d_val > 0 means i is better than j (x > y in calculation).
                    sig_matrix(i, j) = true;
                    if nargout >= 6
                        pairwise_swaps(swap_count, :) = [i, j, sorted_p(k), d_val, r_val];
                    end
                else % d_val < 0 means j is better than i.
                    sig_matrix(j, i) = true;
                    if nargout >= 6
                        pairwise_swaps(swap_count, :) = [j, i, sorted_p(k), d_val, r_val];
                    end
                end
            end
        end
    end
    
    % Stores the final matrices for the current metric.
    all_sig_matrices{metric_idx} = sig_matrix;
    all_alpha_matrices{metric_idx} = alpha_matrix;
    all_p_value_matrices{metric_idx} = p_value_matrix;
    all_pairwise_swaps{metric_idx} = pairwise_swaps(1:swap_count, :);
end

%% 2. Initial Ranking based on Metric 1
% Extracts the results for the first (primary) metric.
sig_metric1_matrix = all_sig_matrices{1};

if isempty(subset_indices)
    mean_metric1 = mean(all_data{1}, 1)';
else
    mean_metric1 = mean(all_data{1}(subset_indices, :), 1)';
end

% Counts the number of "wins" (significant comparisons won) for each dataset.
metric1_wins = sum(sig_metric1_matrix, 2); 
% Starts with the initial order 1, 2, 3, ...
final_order = (1:num_datasets)';

% Iterative sorting algorithm (Bubble-Sort variant) to create the initial ranking.
swapped = true;
while swapped
    swapped = false;
    % Iterates through the current ranking and compares adjacent datasets.
    for k = 1:(num_datasets - 1)
        i = final_order(k);     % Higher-ranked dataset.
        j = final_order(k + 1); % Lower-ranked dataset.
        should_swap = false;
        
        % Primary criterion: Swap if the lower-ranked dataset has more wins.
        if metric1_wins(i) < metric1_wins(j)
            should_swap = true;
        % Secondary criterion: If the number of wins is equal.
        elseif metric1_wins(i) == metric1_wins(j)
            % Check stochastic dominance (Cliff's Delta).
            idx1 = min(i, j); idx2 = max(i, j);
            row_idx = find(pair_idx_all(:,1) == idx1 & pair_idx_all(:,2) == idx2, 1);            
            if ~isempty(row_idx)
                d_val = effect_sizes.d_vals_all(row_idx, 1);               
                % Swap if j stochastically dominates i.
                if (i == idx1 && d_val < 0) || (j == idx1 && d_val > 0)
                    should_swap = true;
                % Swap if Cliffs d is identical (d=0) but j has a higher mean.
                elseif d_val == 0 && mean_metric1(j) > mean_metric1(i)
                     should_swap = true;
                end
            end
        end
        
        if should_swap
            final_order([k, k + 1]) = final_order([k + 1, k]); % Perform swap.
            swapped = true; % Mark that another iteration of the while loop is needed.
        end
    end
end
final_order = final_order'; % Convert to a row vector.
intermediate_orders.after_metric1 = final_order; % Stores the ranking after stage 1.

%% 3. Rank Swapping based on Metric 2 
% Extracts the significance matrix for the second metric.
% This section only runs if the logic mode includes Metric 2 correction.
if num_metrics >= 2 && (strcmp(ranking_mode, 'M1_M2') || strcmp(ranking_mode, 'M1_M2_M3'))
    sig_metric2_matrix = all_sig_matrices{2};
    metric2_global_swaps = zeros(m, 2);
    swap_count_metric2 = 0;

    % Backup state for Rollback safety
    backup_order_m2 = final_order;

    % Iterative global sorting process to correct the ranking.
    changed = true;
    iteration_count = 0;
    max_iterations = num_datasets^2 + 100; % Safety limit against theoretical cycles
    
    while changed && iteration_count < max_iterations
        changed = false;
        iteration_count = iteration_count + 1;
        
        % Each dataset is compared with every higher-ranked dataset.
        for i = 1:(num_datasets - 1)
            for j = (i + 1):num_datasets
                d1_idx = final_order(i); % Higher-ranked dataset.
                d2_idx = final_order(j); % Lower-ranked dataset.
                
                % Swap if the lower-ranked dataset (d2) significantly beats the higher-ranked one (d1).
                if sig_metric2_matrix(d2_idx, d1_idx)
                    % Log the swap for output.
                    swap_count_metric2 = swap_count_metric2 + 1;
                    metric2_global_swaps(swap_count_metric2, :) = [d2_idx, d1_idx];
                    
                    % Perform swap.
                    final_order([i, j]) = final_order([j, i]);
                    changed = true;
                end
            end
        end
    end
    
    if iteration_count >= max_iterations
        warning('Ranking stopped at max iterations. Cycle detected. Reverting to initial order.');
        final_order = backup_order_m2;
        swap_count_metric2 = 0;
        metric2_global_swaps = zeros(0, 2);
    end
    
    metric2_global_swaps = unique(metric2_global_swaps(1:swap_count_metric2, :), 'rows');
    intermediate_orders.after_metric2 = final_order; % Stores the ranking after stage 2.
else
    % If logic 2 is skipped, pass the order from metric 1 and set logs to empty.
    intermediate_orders.after_metric2 = intermediate_orders.after_metric1;
    metric2_global_swaps = zeros(0, 2); % Empty array with correct shape
    
    % Define sig_metric2_matrix as empty or all false for M1_M3A logic,
    % which checks M2 neutrality (which it is, as it wasn't run)
    if num_metrics >= 2
        sig_metric2_matrix = all_sig_matrices{2}; % Still need this for M1_M3A
    else
        sig_metric2_matrix = false(num_datasets); % Placeholder if no M2
    end
end
%% 4. Rank Swapping based on Metric 3
% This section is conditional based on the ranking mode.
metric3_swaps_a = zeros(0, 2); % Initialize empty
metric3_swaps_b = zeros(0, 2); % Initialize empty
if num_metrics == 2 && strcmp(ranking_mode, 'M1_M3A')
    % Logic 3A (Tie-Break) using Metric 
    % In this mode, Metric 1 is used for neutrality check, Metric 2 is the tie-breaker.
    sig_metric1_neutrality_matrix = all_sig_matrices{1};
    sig_metric2_tiebreak_matrix = all_sig_matrices{2};
    
    metric3_swaps_a = zeros(m, 2); % Log for swap logic A.
    swap_count_metric3a = 0;
    has_moved = false(1, num_datasets); % Prevents a dataset from being moved multiple times.
    % Logic 1 (A): One-time swap if Metric 1 was neutral.
    for idx = 1:(num_datasets - 1)
        d1 = final_order(idx);
        d2 = final_order(idx + 1);
        
        % Checks the conditions for a swap.
        if ~has_moved(d1) && ~has_moved(d2) && ...
           ~sig_metric1_neutrality_matrix(d1, d2) && ~sig_metric1_neutrality_matrix(d2, d1) && ... % Metric 1 is neutral.
           sig_metric2_tiebreak_matrix(d2, d1) % Metric 2 (as tiebreaker) shows a significant win for the lower rank d2.
            
            final_order([idx, idx + 1]) = final_order([idx + 1, idx]);
            swap_count_metric3a = swap_count_metric3a + 1;
            metric3_swaps_a(swap_count_metric3a, :) = [d2, d1];
            has_moved([d1, d2]) = true;
        end
    end
    metric3_swaps_a = metric3_swaps_a(1:swap_count_metric3a, :);
    % Logic 2 (B) is explicitly skipped in this mode.
    
elseif num_metrics == 3 && strcmp(ranking_mode, 'M1_M2_M3')
    % Logic 3A and 3B using Metric 3 
    sig_metric2_neutrality_matrix = all_sig_matrices{2};
    sig_metric3_matrix = all_sig_matrices{3};
    metric3_swaps_a = zeros(m, 2); % Log for swap logic A.
    swap_count_metric3a = 0;
    has_moved = false(1, num_datasets); % Prevents a dataset from being moved multiple times.
    % Logic 1 (A): One-time swap if Metric 2 was neutral.
    for idx = 1:(num_datasets - 1)
        d1 = final_order(idx);
        d2 = final_order(idx + 1);
        
        % Checks the conditions for a swap.
        if ~has_moved(d1) && ~has_moved(d2) && ...
           ~sig_metric2_neutrality_matrix(d1, d2) && ~sig_metric2_neutrality_matrix(d2, d1) && ... % Metric 2 is neutral.
           sig_metric3_matrix(d2, d1) % Metric 3 shows a significant win for the lower rank d2.
            
            final_order([idx, idx + 1]) = final_order([idx + 1, idx]);
            swap_count_metric3a = swap_count_metric3a + 1;
            metric3_swaps_a(swap_count_metric3a, :) = [d2, d1];
            has_moved([d1, d2]) = true;
        end
    end
    metric3_swaps_a = metric3_swaps_a(1:swap_count_metric3a, :);
    
    % Logic 2 (B): Iterative swap if Metric 1 AND 2 were neutral.
    metric3_swaps_b = zeros(m, 2); % Log for swap logic B.
    swap_count_metric3b = 0;
    
    % Backup state for Rollback safety
    backup_order_3b = final_order;
    
    changed = true;
    iteration_count = 0;
    max_iterations = num_datasets^2 + 100; % Safety limit against theoretical cycles
    
    while changed && iteration_count < max_iterations
        changed = false;
        iteration_count = iteration_count + 1;
        
        for k = 1:(num_datasets - 1)
            d1 = final_order(k);
            d2 = final_order(k + 1);
            
            % Checks if the first two metrics could not make a significant decision.
            metric1_neutral = ~sig_metric1_matrix(d1, d2) && ~sig_metric1_matrix(d2, d1);
            metric2_neutral = ~sig_metric2_matrix(d1, d2) && ~sig_metric2_matrix(d2, d1);
            
            % If both are neutral and Metric 3 shows a significant difference.
            if metric1_neutral && metric2_neutral && sig_metric3_matrix(d2, d1)
                final_order([k, k + 1]) = final_order([k + 1, k]);
                swap_count_metric3b = swap_count_metric3b + 1;
                metric3_swaps_b(swap_count_metric3b, :) = [d2, d1];
                changed = true;
            end
        end
    end
    
    if iteration_count >= max_iterations
        warning('Ranking stopped at max iterations. Cycle detected. Reverting to previous order.');
        final_order = backup_order_3b;
        swap_count_metric3b = 0;
        metric3_swaps_b = zeros(0, 2);
    end
    
    metric3_swaps_b = unique(metric3_swaps_b(1:swap_count_metric3b, :), 'rows');
end
% Always set the final intermediate order, even if M3 logic was skipped.
intermediate_orders.after_metric3 = final_order;

%% 5. Finalization and Preparation for Outputs
% Creates the final rank vector (e.g., final_rank(5) = 1 means dataset 5 is rank 1).
final_rank = zeros(num_datasets, 1);
for r = 1:num_datasets
    final_rank(final_order(r)) = r;
end

% Collects all detail information for output and logging in a struct.
if nargout >= 6
    swap_details = struct();
    swap_details.metric1_wins = metric1_wins;
    swap_details.pairwise_swaps_metric1 = all_pairwise_swaps{1};
    swap_details.results_metric1 = all_temp_results{1};
    
    % Dynamically add details for M2 and M3 if they exist
    if num_metrics >= 2
        swap_details.pairwise_swaps_metric2 = all_pairwise_swaps{2};
        swap_details.results_metric2 = all_temp_results{2};
        swap_details.metric2_global_swaps = metric2_global_swaps;
    else
        % Add empty placeholders if M2 was not run
        swap_details.pairwise_swaps_metric2 = zeros(0, 5);
        swap_details.results_metric2 = cell(m, 7);
        swap_details.metric2_global_swaps = zeros(0, 2);
    end
    
    if num_metrics == 3
        swap_details.pairwise_swaps_metric3 = all_pairwise_swaps{3};
        swap_details.results_metric3 = all_temp_results{3};
    else
        % Add empty placeholders if M3 was not run
        swap_details.pairwise_swaps_metric3 = zeros(0, 5);
        swap_details.results_metric3 = cell(m, 7);
    end
    
    % These are always added, even if empty (for M1 or M1_M2 modes)
    swap_details.metric3_swaps_a = metric3_swaps_a;
    swap_details.metric3_swaps_b = metric3_swaps_b;
else
    swap_details = struct(); % Return empty if not requested
end
end