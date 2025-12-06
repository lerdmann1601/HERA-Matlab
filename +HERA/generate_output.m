function generate_output(results, thresholds, config, shared_info)
% GENERATE_OUTPUT - Visualizes analysis results and saves them to csv files.
%
% Syntax:
%   generate_output(results, thresholds, config, shared_info)
%
% Description:
%   This function processes the final results of a statistical ranking analysis.
%   It dynamically adapts all console and file outputs based on the number of metrics (1, 2, or 3) and the selected 'ranking_mode'.
%
% Workflow:
%   1. Display of statistical thresholds: 
%      Outputs the thresholds used for the analysis to the console (dynamically for 1-3 metrics).
%   2. Presentation of the hierarchical ranking: 
%      Shows the step-by-step adjustment of the ranking based on the selected logic (M1, M1_M2, M1_M3A, M1_M2_M3).
%   3. Creation of the final results table: 
%      Summarizes the final ranks, confidence intervals, and metric statistics in a table (console and CSV, dynamic columns).
%   4. Creation of the sensitivity analysis table (optional): 
%      Displays the Borda consensus rank and scores for all tested permutations (console and CSV).
%   5. Saving the results and the log: 
%      Writes the detailed log of all comparisons to a separate CSV file (dynamically includes sections for M1, M2, M3A, or M3B).
%
% Inputs:
%   results     - Struct containing all analysis results, such as ranks and swap details.
%   thresholds  - Struct with the thresholds for Cliff's Delta and relative differences.
%   config      - Struct with the global configuration settings (incl. 'ranking_mode').
%   shared_info - Struct with general information like metric and dataset names.
%
% Outputs:
%   The function does not have direct return values. Instead, it generates files:
%   Results File     - A CSV file ('*_results.csv') with the final ranking table.
%   Sensitivity File - A CSV file ('*_sensitivity.csv') with the Borda/permutation results (if applicable).
%   Log File         - A CSV file ('*_log.csv') with the detailed log of all comparisons.
%
% Author: Lukas von Erdmannsdorff

%% 1. Initialization and Data Extraction
% Unpack the passed structures into local variables for easier access.
lang = shared_info.lang;
metric_names = shared_info.metric_names;
dataset_names = shared_info.dataset_names;
num_datasets = shared_info.num_datasets;
num_probanden = shared_info.num_probanden;
pair_idx_all = shared_info.pair_idx_all;
csv_dir = shared_info.csv_dir;
log_basename = shared_info.log_basename;
[~, base_name, ~] = fileparts(log_basename);
ts = shared_info.config.timestamp;
num_metrics = numel(metric_names);
ranking_mode = config.ranking_mode;
mean_metrics = shared_info.mean_metrics;
std_metrics = shared_info.std_metrics;
selected_B_rank = shared_info.selected_B_rank;

d_thresh = thresholds.d_thresh;
rel_thresh_b = thresholds.rel_thresh_b;
rel_thresh = thresholds.rel_thresh;
min_rel_thresh = thresholds.min_rel_thresh;

final_rank = results.final_rank;
final_bootstrap_ranks = results.final_bootstrap_ranks;
all_sig_matrices = results.all_sig_matrices;
all_alpha_matrices = results.all_alpha_matrices;
all_p_value_matrices = results.all_p_value_matrices;
ci_d_all = results.ci_d_all;
ci_r_all = results.ci_r_all;
swap_details = results.swap_details;
intermediate_orders = results.intermediate_orders;
final_order = results.final_order;

% Extract the point estimate and the 95% confidence intervals of the ranks (already calculated in bootstrap_ranking)
point_estimate_ranks = final_rank; 
ci_lower_rank = results.ci_lower_rank;
ci_upper_rank = results.ci_upper_rank;

%% 2. Console Output: General Information and Thresholds
metric_name_list_str = strjoin(metric_names, ', ');
fprintf(['\n' lang.output.info.header '\n'], metric_name_list_str);
fprintf([lang.output.info.mean_sd_intro '\n']);
for i = 1:num_metrics
    fprintf(' -> %s \n', metric_names{i});
end

fprintf(['\n' lang.output.info.significance_intro '\n']);
fprintf([lang.output.info.comparison_details '\n'],...
    size(pair_idx_all, 1), num_datasets, num_probanden, metric_name_list_str);
fprintf([lang.output.info.wilcoxon_desc '\n']);
fprintf([lang.output.info.cliffs_delta_desc '\n']);
fprintf([lang.output.info.rel_diff_desc '\n']);

fprintf(['\n' lang.output.thresholds.delta_header '\n']);
for i = 1:num_metrics
    fprintf(' %s: %.3f\n', metric_names{i}, d_thresh(i));
end

fprintf(['\n' lang.output.thresholds.rel_diff_header '\n']);
for i = 1:num_metrics
    fprintf(' %s: %.3f\n', metric_names{i}, rel_thresh_b(i));
end

fprintf([lang.output.thresholds.sem_protection_header]);
fprintf(['\n-> ' lang.output.thresholds.sem_protection_desc '\n']);
for i = 1:num_metrics
    fprintf(' %s: %.3f\n', metric_names{i}, min_rel_thresh(i));
end

fprintf(lang.output.info.win_criteria_header);
fprintf(['\n 1. ' lang.output.info.win_criterion1]);
fprintf(['\n 2. ' lang.output.info.win_criterion2]);
fprintf(['\n 3. ' lang.output.info.win_criterion3]);

fprintf(['\n\n' lang.output.info.bootstrap_summary_thr '\n'], shared_info.selected_B_thresholds);
fprintf([lang.output.info.bootstrap_summary_ci '\n'], shared_info.selected_B_ci);
fprintf([lang.output.info.bootstrap_summary_rank '\n'], selected_B_rank);

%% 3. Console Output Metric 1
% Detailed info of the comparisons
% This section is always run, as Metric 1 is the basis for all logic modes.
fprintf(['\n' lang.output.metric1.header '\n'], metric_names{1});
fprintf(['\n' lang.output.metric1.logic1 '\n']);
fprintf([lang.output.metric1.criteria '\n'], ...
    d_thresh(1), rel_thresh(1));
fprintf([ lang.output.metric1.logic2 '\n'], metric_names{1});
fprintf(['\n' lang.output.metric1.swaps_header '\n\n'], metric_names{1});

% This block prepares the data for the log file by collecting the reasons for wins in Metric 1.
metric1_ranked_swaps = zeros(size(pair_idx_all, 1), 6);
swap_count_metric1 = 0;
% Logic 1: Direct win by having a higher number of "won" comparisons.
for k = 1:size(swap_details.pairwise_swaps_metric1, 1)
    i = swap_details.pairwise_swaps_metric1(k, 1);
    j = swap_details.pairwise_swaps_metric1(k, 2);
    % Checks if one dataset has more wins than the other
    if swap_details.metric1_wins(i) > swap_details.metric1_wins(j) || swap_details.metric1_wins(j) > swap_details.metric1_wins(i)
        swap_count_metric1 = swap_count_metric1 + 1;
        % Format: [winner, loser, logic_code, p, d, r]
        if swap_details.metric1_wins(i) > swap_details.metric1_wins(j)
             metric1_ranked_swaps(swap_count_metric1, :) = [i, j, 1, swap_details.pairwise_swaps_metric1(k, 3:5)];
        else
             metric1_ranked_swaps(swap_count_metric1, :) = [j, i, 1, swap_details.pairwise_swaps_metric1(k, 3:5)];
        end
    end
end
% Logic 2: Tie-breaking by Stochastic Dominance (Cliff's Delta) and if equal by mean.
for k = 1:(num_datasets - 1)
    for l = (k + 1):num_datasets
        % Checks if two datasets have the same number of wins.
        if swap_details.metric1_wins(k) == swap_details.metric1_wins(l)        
            % Determine if k or l wins based on logic 2 (Cliff's Delta -> Mean).
            is_k_better = false; is_l_better = false;            
            % Find indices in results table
            indices = cell2mat(swap_details.results_metric1(:, 6:7));
            row_idx = find(ismember(indices, [k, l; l, k], 'rows'), 1);          
            if ~isempty(row_idx)
                d_val = swap_details.results_metric1{row_idx, 4};
                idx1 = swap_details.results_metric1{row_idx, 6};              
                % Check Cliff's Delta (Logic 2a)
                if (idx1 == k && d_val > 0) || (idx1 == l && d_val < 0)
                    is_k_better = true;
                elseif (idx1 == l && d_val > 0) || (idx1 == k && d_val < 0)
                    is_l_better = true;                  
                % Fallback: Use mean (Logic 2b) if d=0
                elseif d_val == 0
                     if (mean_metrics(k, 1) > mean_metrics(l, 1)), is_k_better = true;
                     elseif (mean_metrics(l, 1) > mean_metrics(k, 1)), is_l_better = true; end
                end
                if is_k_better || is_l_better
                    winner = k; loser = l;
                    if is_l_better, winner = l; loser = k; end
                    swap_count_metric1 = swap_count_metric1 + 1;
                    metric1_ranked_swaps(swap_count_metric1, :) = [winner, loser, 2, swap_details.results_metric1{row_idx, 3:5}];
                end
            end
        end
    end
end
metric1_ranked_swaps = unique(metric1_ranked_swaps(1:swap_count_metric1, :), 'rows');

% Console output for Metric 1
if isempty(metric1_ranked_swaps)
    fprintf([lang.output.log.none_found '\n']);
else
    % Prints a formatted table of all significant comparisons for Metric 1.
    % Prepare all data and headers for the console from the language file
    header_parts = {lang.output.tables.logic_prefix, lang.output.log.headers_no_power{1}, ...
        lang.output.tables.vs, lang.output.log.headers_no_power{2}, lang.output.log.headers_no_power{4}, ...
        lang.output.tables.delta_short, lang.output.log.headers_no_power{7}, ...
        lang.output.tables.reldiff_short, lang.output.log.headers_no_power{10}};

    % Sorts the results: first by logic (1 or 2), then by the winner's index.
    [~, sort_idx] = sortrows(metric1_ranked_swaps, [3,1]);
    sorted_swaps = metric1_ranked_swaps(sort_idx,:);
    
    table_data = cell(size(sorted_swaps, 1), numel(header_parts));
    for k = 1:size(sorted_swaps, 1)
        winner = sorted_swaps(k, 1);
        loser = sorted_swaps(k, 2);
        
        pair_row_idx = find((pair_idx_all(:, 1) == winner & pair_idx_all(:, 2) == loser) | (pair_idx_all(:, 1) == loser & pair_idx_all(:, 2) == winner));
        
        table_data{k, 1} = sprintf('(%d)', sorted_swaps(k, 3));
        table_data{k, 2} = dataset_names{winner};
        table_data{k, 3} = lang.output.tables.vs;
        table_data{k, 4} = dataset_names{loser};
        table_data{k, 5} = sprintf('%.4f', sorted_swaps(k, 4));
        table_data{k, 6} = sprintf('%+.3f', sorted_swaps(k, 5));
        table_data{k, 7} = sprintf('[%+.3f, %+.3f]', ci_d_all(pair_row_idx, 1, 1), ci_d_all(pair_row_idx, 2, 1));
        table_data{k, 8} = sprintf('%.3f', sorted_swaps(k, 6));
        table_data{k, 9} = sprintf('[%.3f, %.3f]', ci_r_all(pair_row_idx, 1, 1), ci_r_all(pair_row_idx, 2, 1));
    end
    
    % Calculate column widths
    col_widths = cellfun(@strlength, header_parts);
    for r = 1:size(table_data, 1)
        for c = 1:size(table_data, 2)
            col_widths(c) = max(col_widths(c), strlength(table_data{r, c}));
        end
    end
    col_widths = col_widths + 2; % 2 spaces for padding
            
    % Define the alignment for data rows
    num_cols = numel(header_parts);
    data_alignments = repmat({'c'}, 1, num_cols);
    data_alignments{2} = 'r'; 
    data_alignments{4} = 'l';
    
    % Header fully centered
    header_line = strjoin(arrayfun(@(c) format_text(header_parts{c}, col_widths(c), 'c'), 1:num_cols, 'UniformOutput', false), '|');
    fprintf('%s\n', header_line);
    fprintf('%s\n', repmat('-', 1, strlength(header_line)));
    
    % Data rows (use data_alignments)
    for r = 1:size(table_data, 1)
        row_line = strjoin(arrayfun(@(c) format_text(table_data{r, c}, col_widths(c), data_alignments{c}), 1:num_cols, 'UniformOutput', false), '|');
        fprintf('%s\n', row_line);
    end
    fprintf('%s\n', repmat('-', 1, strlength(header_line)));
end
% Output of the initial ranking after sorting with Metric 1.
fprintf(['\n' lang.output.metric1.initial_ranking '\n\n'], metric_names{1});
for r = 1:num_datasets
    dset_idx = intermediate_orders.after_metric1(r);
    fprintf([lang.output.tables.format '\n'], r, dataset_names{dset_idx}, metric_names{1}, mean_metrics(dset_idx, 1));
end

%% 4. Console Output Metric 2
% Only show if Metric 2 logic (M1_M2 or M1_M2_M3) was applied
if num_metrics >= 2 && (strcmp(ranking_mode, 'M1_M2') || strcmp(ranking_mode, 'M1_M2_M3'))
    
    fprintf(['\n' lang.output.metric2.header '\n'], metric_names{2});
    fprintf(['\n' lang.output.metric2.logic '\n'], metric_names{2});
    fprintf([lang.output.metric2.criteria '\n'], ...
        d_thresh(2), rel_thresh(2));
    fprintf(['\n' lang.output.metric2.swaps_header '\n\n'], metric_names{2});
    pairwise_swaps_metric2 = swap_details.pairwise_swaps_metric2;
    if isempty(pairwise_swaps_metric2)
        fprintf([lang.output.log.none_found '\n']);
    else
        % Prepare all data and headers for the console from the language file
        header_parts = {lang.output.log.headers_no_power{1}, lang.output.tables.vs, ...
            lang.output.log.headers_no_power{2}, lang.output.log.headers_no_power{4}, ...
            lang.output.tables.delta_short, lang.output.log.headers_no_power{7}, ...
            lang.output.tables.reldiff_short, lang.output.log.headers_no_power{10}};
            
        table_data = cell(size(pairwise_swaps_metric2, 1), numel(header_parts));
        for k = 1:size(pairwise_swaps_metric2, 1)
            i1 = pairwise_swaps_metric2(k, 1);
            i2 = pairwise_swaps_metric2(k, 2);
            pair_row_idx = find((pair_idx_all(:, 1) == i1 & pair_idx_all(:, 2) == i2) | (pair_idx_all(:, 1) == i2 & pair_idx_all(:, 2) == i1));
            
            table_data{k, 1} = dataset_names{i1};
            table_data{k, 2} = lang.output.tables.vs;
            table_data{k, 3} = dataset_names{i2};
            table_data{k, 4} = sprintf('%.4f', pairwise_swaps_metric2(k, 3));
            table_data{k, 5} = sprintf('%+.3f', pairwise_swaps_metric2(k, 4));
            table_data{k, 6} = sprintf('[%+.3f, %+.3f]', ci_d_all(pair_row_idx, 1, 2), ci_d_all(pair_row_idx, 2, 2));
            table_data{k, 7} = sprintf('%.3f', pairwise_swaps_metric2(k, 5));
            table_data{k, 8} = sprintf('[%.3f, %.3f]', ci_r_all(pair_row_idx, 1, 2), ci_r_all(pair_row_idx, 2, 2));
        end
        % Calculate column widths
        col_widths = cellfun(@strlength, header_parts);
        for r = 1:size(table_data, 1)
            for c = 1:size(table_data, 2)
                col_widths(c) = max(col_widths(c), strlength(table_data{r, c}));
            end
        end
        col_widths = col_widths + 2; % 2 spaces for padding
                
        % Define alignment for data rows
        num_cols = numel(header_parts);
        data_alignments = repmat({'c'}, 1, num_cols);
        data_alignments{1} = 'r';
        data_alignments{3} = 'l';
        
        % Header fully centered
        header_line = strjoin(arrayfun(@(c) format_text(header_parts{c}, col_widths(c), 'c'), 1:num_cols, 'UniformOutput', false), '|');
        fprintf('%s\n', header_line);
        fprintf('%s\n', repmat('-', 1, strlength(header_line)));
        
        % Data rows (use data_alignments)
        for r = 1:size(table_data, 1)
            row_line = strjoin(arrayfun(@(c) format_text(table_data{r, c}, col_widths(c), data_alignments{c}), 1:num_cols, 'UniformOutput', false), '|');
            fprintf('%s\n', row_line);
        end
        fprintf('%s\n', repmat('-', 1, strlength(header_line)));
    end
    % Output of the ranking after correction by Metric 2.
    fprintf(['\n' lang.output.metric2.ranking_after '\n\n'], metric_names{2});
    for r = 1:num_datasets
        idx = intermediate_orders.after_metric2(r);
        fprintf([lang.output.tables.format '\n'], r, dataset_names{idx}, metric_names{2}, mean_metrics(idx, 2));
    end
end % End of conditional block for Metric 2

%% 5. Console Output Metric 3
% Conditional output based on logic mode
if num_metrics == 2 && strcmp(ranking_mode, 'M1_M3A')
    % Output for M1_M3A logic 
    fprintf(['\n' lang.output.metric3.header '\n'], metric_names{2}); % Uses Metric 2 as the "tie-breaker"
    fprintf(['\n' lang.output.metric3.logic1_3A '\n'], metric_names{1}, metric_names{2});
    fprintf([lang.output.metric3.criteria '\n'], ...
        d_thresh(2), rel_thresh(2)); % Use Metric 2 thresholds
    fprintf(['\n' lang.output.metric3.swaps_header '\n\n'], metric_names{2});
    
    % Note: In M1_M3A mode, only pairwise_swaps_metric3 (from M2 data) and metric3_swaps_a exist.
    pairwise_swaps_metric3_logicA = swap_details.pairwise_swaps_metric2; % Data comes from M2
    swaps_to_print = zeros(0,5);
    
    % Filter the swaps_metric2 list to only those that were part of Logic 3A
    if ~isempty(pairwise_swaps_metric3_logicA) && ~isempty(swap_details.metric3_swaps_a)
         % Find rows in M2-swaps that match the pairs in metric3_swaps_a
         [~, ia, ~] = intersect(pairwise_swaps_metric3_logicA(:,1:2), swap_details.metric3_swaps_a, 'rows');
         swaps_to_print = pairwise_swaps_metric3_logicA(ia,:);
    end

    if isempty(swaps_to_print)
        fprintf([lang.output.log.none_found '\n']);
    else
        % Prints a formatted table (same format as M2)
        header_parts = {lang.output.log.headers_no_power{1}, lang.output.tables.vs, ...
            lang.output.log.headers_no_power{2}, lang.output.log.headers_no_power{4}, ...
            lang.output.tables.delta_short, lang.output.log.headers_no_power{7}, ...
            lang.output.tables.reldiff_short, lang.output.log.headers_no_power{10}};
            
        table_data = cell(size(swaps_to_print, 1), numel(header_parts));
        for k = 1:size(swaps_to_print, 1)
            i1 = swaps_to_print(k, 1); i2 = swaps_to_print(k, 2);
            pair_row_idx = find((pair_idx_all(:, 1) == i1 & pair_idx_all(:, 2) == i2) | (pair_idx_all(:, 1) == i2 & pair_idx_all(:, 2) == i1));
            
            table_data{k, 1} = dataset_names{i1};
            table_data{k, 2} = lang.output.tables.vs;
            table_data{k, 3} = dataset_names{i2};
            table_data{k, 4} = sprintf('%.4f', swaps_to_print(k, 3));
            table_data{k, 5} = sprintf('%+.3f', swaps_to_print(k, 4));
            table_data{k, 6} = sprintf('[%+.3f, %+.3f]', ci_d_all(pair_row_idx, 1, 2), ci_d_all(pair_row_idx, 2, 2)); % Use M2 CI
            table_data{k, 7} = sprintf('%.3f', swaps_to_print(k, 5));
            table_data{k, 8} = sprintf('[%.3f, %.3f]', ci_r_all(pair_row_idx, 1, 2), ci_r_all(pair_row_idx, 2, 2)); % Use M2 CI
        end
        % Calculate column widths
        col_widths = cellfun(@strlength, header_parts);
        for r = 1:size(table_data, 1), for c = 1:size(table_data, 2), col_widths(c) = max(col_widths(c), strlength(table_data{r, c})); end, end
        col_widths = col_widths + 2;
        num_cols = numel(header_parts);
        data_alignments = repmat({'c'}, 1, num_cols); data_alignments{1} = 'r'; data_alignments{3} = 'l';
        % Header fully centered
        header_line = strjoin(arrayfun(@(c) format_text(header_parts{c}, col_widths(c), 'c'), 1:num_cols, 'UniformOutput', false), '|');
        fprintf('%s\n', header_line); fprintf('%s\n', repmat('-', 1, strlength(header_line)));
        % Data rows (use data_alignments)
        for r = 1:size(table_data, 1)
            row_line = strjoin(arrayfun(@(c) format_text(table_data{r, c}, col_widths(c), data_alignments{c}), 1:num_cols, 'UniformOutput', false), '|');
            fprintf('%s\n', row_line);
        end
        fprintf('%s\n', repmat('-', 1, strlength(header_line)));
    end
    % Output of the final ranking after all corrections.
    fprintf(['\n' lang.output.metric3.final_ranking '\n\n'], metric_names{2}); % Still M2 data
    for r = 1:num_datasets
        idx = final_order(r); % 'final_order' contains the final sequence of indices.
        fprintf([lang.output.tables.format '\n'], r, dataset_names{idx}, metric_names{2}, mean_metrics(idx, 2));
    end

elseif num_metrics == 3 && strcmp(ranking_mode, 'M1_M2_M3')
    %  Output for M1_M2_M3 logic 
    fprintf(['\n' lang.output.metric3.header '\n'], metric_names{3});
    fprintf(['\n' lang.output.metric3.logic1 '\n'], metric_names{2}, metric_names{3});
    fprintf([lang.output.metric3.logic2 '\n'], metric_names{1}, metric_names{2}, metric_names{3});
    fprintf([lang.output.metric3.criteria '\n'], ...
        d_thresh(3), rel_thresh(3));
    fprintf(['\n' lang.output.metric3.swaps_header '\n\n'], metric_names{3});
    pairwise_swaps_metric3 = swap_details.pairwise_swaps_metric3;
    if isempty(pairwise_swaps_metric3)
        fprintf([lang.output.log.none_found '\n']);
    else
        % Prints a formatted table of all significant comparisons for Metric 3.
        % Prepare all data and headers for the console from the language file
        header_parts = {lang.output.tables.logic_prefix, lang.output.log.headers_no_power{1}, ...
            lang.output.tables.vs, lang.output.log.headers_no_power{2}, lang.output.log.headers_no_power{4}, ...
            lang.output.tables.delta_short, lang.output.log.headers_no_power{7}, ...
            lang.output.tables.reldiff_short, lang.output.log.headers_no_power{10}};
        
        % Sorts the swaps to group them by the underlying logic (1 or 2).
        sort_key = 3 * ones(size(pairwise_swaps_metric3, 1), 1);
        is_in_a = ismember(pairwise_swaps_metric3(:, 1:2), swap_details.metric3_swaps_a, 'rows');
        is_in_b = ismember(pairwise_swaps_metric3(:, 1:2), swap_details.metric3_swaps_b, 'rows');
        sort_key(is_in_a) = 1;
        sort_key(is_in_b & sort_key > 2) = 2;
        [~, sorted_indices] = sort(sort_key);
        sorted_swaps = pairwise_swaps_metric3(sorted_indices, :);
        final_sort_key = sort_key(sorted_indices);
        
        table_data = cell(size(sorted_swaps, 1), numel(header_parts));
        for k = 1:size(sorted_swaps, 1)
            i1 = sorted_swaps(k, 1); i2 = sorted_swaps(k, 2);
            pair_row_idx = find((pair_idx_all(:, 1) == i1 & pair_idx_all(:, 2) == i2) | (pair_idx_all(:, 1) == i2 & pair_idx_all(:, 2) == i1));
            
            prefix = ' ';
            if final_sort_key(k) == 1, prefix = '(1)';
            elseif final_sort_key(k) == 2, prefix = '(2)'; end
            
            table_data{k, 1} = prefix;
            table_data{k, 2} = dataset_names{i1};
            table_data{k, 3} = lang.output.tables.vs;
            table_data{k, 4} = dataset_names{i2};
            table_data{k, 5} = sprintf('%.4f', sorted_swaps(k, 3));
            table_data{k, 6} = sprintf('%+.3f', sorted_swaps(k, 4));
            table_data{k, 7} = sprintf('[%+.3f, %+.3f]', ci_d_all(pair_row_idx, 1, 3), ci_d_all(pair_row_idx, 2, 3));
            table_data{k, 8} = sprintf('%.3f', sorted_swaps(k, 5));
            table_data{k, 9} = sprintf('[%.3f, %.3f]', ci_r_all(pair_row_idx, 1, 3), ci_r_all(pair_row_idx, 2, 3));
        end
        % Calculate column widths
        col_widths = cellfun(@strlength, header_parts);
        for r = 1:size(table_data, 1)
            for c = 1:size(table_data, 2)
                col_widths(c) = max(col_widths(c), strlength(table_data{r, c}));
            end
        end
        col_widths = col_widths + 2; % 2 spaces for padding
                
        num_cols = numel(header_parts);
        data_alignments = repmat({'c'}, 1, num_cols);
        data_alignments{2} = 'r';
        data_alignments{4} = 'l';
        
        % Header fully centered
        header_line = strjoin(arrayfun(@(c) format_text(header_parts{c}, col_widths(c), 'c'), 1:num_cols, 'UniformOutput', false), '|');
        fprintf('%s\n', header_line);
        fprintf('%s\n', repmat('-', 1, strlength(header_line)));
        
        % Data rows (use data_alignments)
        for r = 1:size(table_data, 1)
            row_line = strjoin(arrayfun(@(c) format_text(table_data{r, c}, col_widths(c), data_alignments{c}), 1:num_cols, 'UniformOutput', false), '|');
            fprintf('%s\n', row_line);
        end
        fprintf('%s\n', repmat('-', 1, strlength(header_line)));
    end
    % Output of the final ranking after all corrections.
    fprintf(['\n' lang.output.metric3.final_ranking '\n\n'], metric_names{3});
    for r = 1:num_datasets
        idx = final_order(r); % 'final_order' contains the final sequence of indices.
        fprintf([lang.output.tables.format '\n'], r, dataset_names{idx}, metric_names{3}, mean_metrics(idx, 3));
    end
end % End of conditional M3 block

%% 6. Results CSV (_results.csv)
% Final results table and saving to CSV
% Console output of the final ranking
% Dynamic header
metric_name_list_str = strjoin(metric_names, ', ');
% Build the header string manually to support 1, 2, or 3 metrics
header_string = regexprep(lang.output.final_summary.header, '%s, %s, and %s', '%s'); % Base string
header_string = regexprep(header_string, '%s, %s', '%s'); % Handle 2 metrics
fprintf(['\n' header_string '\n\n'], metric_name_list_str);


% Prepare all data and headers for the console from the language file
% Dynamically build header parts
header_parts = {lang.plots.table.rank, lang.csv.headers.dataset, lang.output.tables.rank_ci};
for i = 1:num_metrics
    header_parts{end+1} = sprintf(lang.csv.headers.mean_sd_template, metric_names{i});
end
    
% Prepare table data by sorting datasets according to their final rank
[~, output_order_idx] = sort(point_estimate_ranks);
table_data = cell(num_datasets, numel(header_parts));
for r = 1:num_datasets
    idx = output_order_idx(r);
    
    % Data for final ranking table
    table_data{r, 1} = sprintf('%d', point_estimate_ranks(idx));
    table_data{r, 2} = dataset_names{idx};
    table_data{r, 3} = sprintf('[%s, %s]', strrep(sprintf('%.2f', ci_lower_rank(idx)), ',', '.'), strrep(sprintf('%.2f', ci_upper_rank(idx)), ',', '.'));
    
    % Loop to add metric data
    for i = 1:num_metrics
        table_data{r, 3+i} = strrep(sprintf('%.3f ± %.3f', mean_metrics(idx, i), std_metrics(idx, i)), ',', '.');
    end
end

% Calculate dynamic column widths for console output
col_widths = cellfun(@strlength, header_parts);
for r = 1:size(table_data, 1)
    for c = 1:size(table_data, 2)
        col_widths(c) = max(col_widths(c), strlength(table_data{r, c}));
    end
end
col_widths = col_widths + 2; % Add 2 spaces for padding

% Define column alignments for console output (all centered)
alignments = repmat({'c'}, 1, numel(header_parts));

% Output header and data rows to console using the helper function
header_line = strjoin(arrayfun(@(c) format_text(header_parts{c}, col_widths(c), alignments{c}), 1:numel(header_parts), 'UniformOutput', false), '|');
fprintf('%s\n', header_line);
fprintf('%s\n', repmat('-', 1, strlength(header_line)));

for r = 1:size(table_data, 1)
    row_line = strjoin(arrayfun(@(c) format_text(table_data{r, c}, col_widths(c), alignments{c}), 1:size(table_data, 2), 'UniformOutput', false), '|');
    fprintf('%s\n', row_line);
end
fprintf('%s\n', repmat('-', 1, strlength(header_line)));

% Write the final results to a CSV file
[~, fName, fExt] = fileparts(lang.files.results);
fName = strrep(fName, '%s_', ''); 
csv_filename_results = fullfile(csv_dir, [fName, '_', ts, fExt]);

try
    % Attempt to open the file for writing
    fid_results = fopen(csv_filename_results, 'w');
    if fid_results == -1
        error(lang.errors.file_open_error, csv_filename_results); 
    end

    % Write header to CSV
    fprintf(fid_results, '%s\n', strjoin(header_parts, ';'));

    % Write data rows to CSV
    for r = 1:size(table_data, 1)
        % Enclose each cell content in quotes for CSV compatibility
        csv_row_cells = cellfun(@(c) ['"' c '"'], table_data(r,:), 'UniformOutput', false);
        fprintf(fid_results, '%s\n', strjoin(csv_row_cells, ';'));
    end
    
    % Close the results file successfully
    fclose(fid_results); 
    fprintf(['\n' lang.output.files.final_results_saved '\n'], csv_filename_results);

catch ME
    % Safety cleanup: Ensure file is closed if an error occurs
    if exist('fid_results', 'var') && fid_results ~= -1
        fclose(fid_results); 
    end
    % Report the error to the user
    fprintf([lang.errors.file_save_error '\n'], ME.message);
end

%% Borda Results
if isfield(results, 'borda_results') && ~isempty(results.borda_results)
    fprintf(['\n' lang.output.sensitivity.header '\n\n']);
    num_perms = size(results.all_permutation_ranks, 2);
    
    % Dynamically calculate column widths
    % Column 1: Borda Rank
    width_borda = max(strlength(lang.output.sensitivity.borda_rank), strlength(sprintf('%d', num_datasets))) + 2; % +2 for padding
    
    % Column 2: Borda Score (%)
    width_score = max(strlength(lang.output.sensitivity.borda_score), strlength('100.0%')) + 2; 
    
    % Column 3: Dataset
    width_dataset = max(strlength(lang.csv.headers.dataset), max(cellfun(@strlength, dataset_names))) + 2;
    
    % Initialize lists for header texts and their widths
    header_parts_list = {lang.output.sensitivity.borda_rank, ...
                         lang.output.sensitivity.borda_score, ...
                         lang.csv.headers.dataset};
    col_widths = [width_borda, width_score, width_dataset];
    
    % Columns 4 to end: Permutations
    perm_header_parts = cell(1, num_perms);
    for p = 1:num_perms
        perm_indices = results.selected_permutations(p, :);
        primary_hierarchy_indices = results.selected_permutations(1, :);
        % This local_indices lookup correctly handles 2 or 3 metrics
        [~, local_indices] = ismember(perm_indices, primary_hierarchy_indices);
        perm_names = metric_names(local_indices);
        col_name_full = strjoin(perm_names, ' -> ');
        perm_header_parts{p} = col_name_full;
        % Width is determined by the length of the permutation title
        width_perm = strlength(col_name_full) + 2;
        col_widths(end+1) = width_perm;
    end
    header_parts_list = [header_parts_list, perm_header_parts];

    % Prepare header for the CSV file
    header_list_csv = cellfun(@(s) sprintf('"%s"', s), header_parts_list, 'UniformOutput', false);

    % Define alignment for each column for console output
    alignments = {'c', 'c', 'l'}; % Borda Rank (c), Borda Score (c), Dataset (l)
    alignments(4:3+num_perms) = {'c'}; % All permutation columns centered
    
    % Output headers to the console using the helper function (all centered)
    header_line = strjoin(arrayfun(@(c) format_text(header_parts_list{c}, col_widths(c), 'c'), 1:numel(header_parts_list), 'UniformOutput', false), ' | ');
    fprintf('%s\n', header_line);
    fprintf('%s\n', repmat('-', 1, strlength(header_line)));

    % Prepare filename for the sensitivity report
    [~, fName, fExt] = fileparts(lang.files.sensitivity_details);
    fName = strrep(fName, '%s_', '');
    csv_filename_sensitivity = fullfile(csv_dir, [fName, '_', ts, fExt]);
    
    try
        % Attempt to open CSV file for output
        fid_sens = fopen(csv_filename_sensitivity, 'w');
        if fid_sens == -1
            error(lang.errors.file_open_error, csv_filename_sensitivity); 
        end
        
        % Write header
        fprintf(fid_sens, '%s\n', strjoin(header_list_csv, ';'));
    
        % Fill the tables (Console & CSV)
        [~, sort_idx_c] = sort(results.borda_results.rank); % Sort by consensus rank
        
        for r_idx = 1:num_datasets
            d_idx = sort_idx_c(r_idx);
            
            % Prepare data for the current row (both for console and CSV)
            table_data_row = cell(1, numel(header_parts_list));
            table_data_row{1} = sprintf('%d', r_idx); % Borda Rank
            table_data_row{2} = sprintf('%.1f%%', results.borda_results.score(d_idx)); % Borda Score
            table_data_row{3} = dataset_names{d_idx}; % Dataset
            
            for p = 1:num_perms
                table_data_row{3 + p} = sprintf('%d', results.all_permutation_ranks(d_idx, p));
            end
            
            % Format and print the console row using the helper function
            row_line_console = strjoin(arrayfun(@(c) format_text(table_data_row{c}, col_widths(c), alignments{c}), ...
                1:numel(header_parts_list), 'UniformOutput', false), ' | ');
            fprintf('%s\n', row_line_console);
            
            % Prepare and write the row to the CSV file
            csv_row_cells = {table_data_row{1}, ... % Rank (number)
                             ['"' table_data_row{2} '"'], ... % Score (string)
                             ['"' table_data_row{3} '"']};   % Dataset (string)
            csv_row_cells = [csv_row_cells, table_data_row(4:end)]; % Permutation ranks (numbers)
            fprintf(fid_sens, '%s\n', strjoin(csv_row_cells, ';'));
        end
        
        % Close the sensitivity file successfully
        fclose(fid_sens);  
        fprintf(['\n' lang.output.files.sensitivity_results_saved '\n'], csv_filename_sensitivity);

    catch ME
        % Safety cleanup: Ensure file is closed if an error occurs
        if exist('fid_sens', 'var') && fid_sens ~= -1
            fclose(fid_sens); 
        end
        fprintf([lang.errors.file_save_error '\n'], ME.message);
    end
end

%% 7. Log CSV (_log.csv)
% Generate the full filename for the log file
[~, fName, fExt] = fileparts(lang.files.log_details);
fName = strrep(fName, '%s_', '');
csv_filename_log = fullfile(csv_dir, [fName, '_', ts, fExt]);

try
    % Attempt to open the log file in write mode
    fid_log = fopen(csv_filename_log, 'w');
    if fid_log == -1
        error(lang.errors.file_open_error, csv_filename_log); 
    end
    
    % Print the header for the detailed log table to the console
    fprintf(['\n' lang.output.log.header '\n\n']);

    % Check if a power analysis was performed
    power_analysis_done = isfield(results, 'power_results') && ~isempty(results.power_results);

    % Define headers
    if power_analysis_done
        header_parts_log = lang.output.log.headers_with_power;
        csv_header = strjoin(cellfun(@(c) ['"' c '"'], header_parts_log, 'UniformOutput', false), ';');
        num_log_cols = 16; 
    else
        header_parts_log = lang.output.log.headers_no_power;
        csv_header = strjoin(cellfun(@(c) ['"' c '"'], header_parts_log, 'UniformOutput', false), ';');
        num_log_cols = 15; 
    end

    % Write the generated CSV header to the log file
    fprintf(fid_log, '%s\n', csv_header);

    % Pre-allocate a cell array to hold all log data before writing to file.
    % This allows sorting the entire log before outputting it.
    m = size(pair_idx_all, 1);
    % Pre-allocate based on the number of metrics (1, 2, or 3)
    log_data = cell(m * num_metrics, num_log_cols);
    log_idx = 0; % Initialize a counter for the rows in the log_data cell array.

    % Iterate through all pairwise comparisons to collect log data for Metric 1.
    for k = 1:m
        log_idx = log_idx + 1; % Increment row index for each new entry.
        % Extract results for the current pair from the detailed results structure.
        res = swap_details.results_metric1(k, :);
        name_i = res{1}; name_j = res{2}; p_obs = res{3}; d_obs = res{4}; rel_obs = res{5}; i_idx = res{6}; j_idx = res{7};
        % Get the specific alpha level used for this comparison.
        alpha_k = all_alpha_matrices{1}(i_idx, j_idx);
        % Determine if the difference was statistically significant in either direction.
        is_sig = all_sig_matrices{1}(i_idx, j_idx) || all_sig_matrices{1}(j_idx, i_idx);
        
        % Find the corresponding row index in the original pairwise list to fetch CI data.
        pair_row_idx = find((pair_idx_all(:, 1) == i_idx & pair_idx_all(:, 2) == j_idx) | (pair_idx_all(:, 1) == j_idx & pair_idx_all(:, 2) == i_idx), 1);
        if isempty(pair_row_idx), warning(lang.warnings.pair_index_not_found, name_i, name_j); pair_row_idx=1; end
        
        % Format the confidence intervals for Cliff's Delta and relative difference as strings.
        ci_d_str = sprintf('[%+.3f, %+.3f]', ci_d_all(pair_row_idx, 1, 1), ci_d_all(pair_row_idx, 2, 1));
        ci_r_str = sprintf('[%.3f, %.3f]', ci_r_all(pair_row_idx, 1, 1), ci_r_all(pair_row_idx, 2, 1));
        
        % Determine the reason for the ranking decision based on Metric 1's logic.
        if swap_details.metric1_wins(i_idx) ~= swap_details.metric1_wins(j_idx)
            % Case 1: One dataset won more pairwise comparisons.
            reason = lang.output.log.reason_m1_wins; sort_key = 1;
            
        elseif abs(d_obs) > 0
            % Case 2: Tie in wins, broken by Stochastic Dominance (Cliff's Delta).
            reason = sprintf(lang.output.log.reason_m1_tiebreak_d, metric_names{1}); sort_key = 2;
            
        else
            % Case 3: Tie in wins AND identical distributions (d=0), broken by Mean.
            reason = sprintf(lang.output.log.reason_m1_tiebreak_mean, metric_names{1}); sort_key = 3;
        end
        
        % Assemble the log entry as a cell array. The 'Swap' column is unused for Metric 1.
        log_entry = {name_i, name_j, metric_names{1}, p_obs, alpha_k, d_obs, ci_d_str, d_thresh(1), rel_obs, ci_r_str, rel_thresh(1), sprintf('%d', is_sig)};
        % If power analysis was done, append the power value.
        if power_analysis_done
            power_val = results.power_results.power_matrices{1}(k);
            log_entry{end+1} = sprintf('%.1f', power_val * 100);
        end
        % Append the swap flag (always empty for M1), the reason, and the sort key.
        log_entry = [log_entry, {' '}, {reason}, {sort_key}];
        log_data(log_idx, :) = log_entry;
    end

    % Iterate through all pairwise comparisons to collect log data for Metric 2 (conditional)
    if num_metrics >= 2 && (strcmp(ranking_mode, 'M1_M2') || strcmp(ranking_mode, 'M1_M2_M3'))
        for k = 1:m
            log_idx = log_idx + 1; % Increment row index.
            % Extract results for the current pair.
            res = swap_details.results_metric2(k, :);
            name_i = res{1}; name_j = res{2}; p_obs = res{3}; d_obs = res{4}; rel_obs = res{5}; i_idx = res{6}; j_idx = res{7};
            alpha_k = all_alpha_matrices{2}(i_idx, j_idx);
            is_sig = all_sig_matrices{2}(i_idx, j_idx) || all_sig_matrices{2}(j_idx, i_idx);
            % Check if this pair was involved in a rank swap based on Metric 2.
            swap_detected =... 
            any(all(swap_details.metric2_global_swaps == [j_idx, i_idx], 2)) || any(all(swap_details.metric2_global_swaps == [i_idx, j_idx], 2));
            swap_flag = sprintf('%d', swap_detected); % Create a '0' or '1' flag for the log.
    
            % Find the corresponding row index for confidence intervals.
            pair_row_idx = find((pair_idx_all(:, 1) == i_idx & pair_idx_all(:, 2) == j_idx) | (pair_idx_all(:, 1) == j_idx & pair_idx_all(:, 2) == i_idx), 1);
            if isempty(pair_row_idx), warning(lang.warnings.pair_index_not_found, name_i, name_j); pair_row_idx=1; end
    
            % Format confidence interval strings.
            ci_d_str = sprintf('[%+.3f, %+.3f]', ci_d_all(pair_row_idx, 1, 2), ci_d_all(pair_row_idx, 2, 2));
            ci_r_str = sprintf('[%.3f, %.3f]', ci_r_all(pair_row_idx, 1, 2), ci_r_all(pair_row_idx, 2, 2));
    
            % Determine the reason for the outcome based on Metric 2's logic.
            if is_sig && ~swap_detected, reason = lang.output.log.reason_m2_correct; sort_key = 2; % Significant, but no swap needed (already correct order).
            elseif is_sig, reason = ' '; sort_key = 1; % Significant and a swap was performed.
            else, reason = sprintf(lang.output.log.reason_m2_no_win, metric_names{2}); sort_key = 3; % Not significant, no action taken.
            end
            
            % Assemble the log entry.
            log_entry =... 
            {name_i, name_j, metric_names{2}, p_obs, alpha_k, d_obs, ci_d_str, d_thresh(2),rel_obs, ci_r_str, rel_thresh(2), sprintf('%d', is_sig)};
            if power_analysis_done
                power_val = results.power_results.power_matrices{2}(k);
                log_entry{end+1} = sprintf('%.1f', power_val * 100);
            end
            % Append the swap flag, the reason, and the sort key.
            log_entry = [log_entry, {swap_flag}, {reason}, {sort_key}];
            log_data(log_idx, :) = log_entry;
        end
    end % End of conditional M2 log

    % Iterate through all pairwise comparisons to collect log data for Metric 3 / Tie-break (CONDITIONAL)
    if (num_metrics == 2 && strcmp(ranking_mode, 'M1_M3A'))
        % M1_M3A Logic Log (uses M2 data) 
        for k = 1:m
            % Extract results for the current pair (using M2 data)
            res = swap_details.results_metric2(k, :);
            name_i = res{1}; name_j = res{2}; p_obs = res{3}; d_obs = res{4}; rel_obs = res{5}; i_idx = res{6}; j_idx = res{7};
            
            % Check if this pair was swapped due to Metric 3A logic
            is_swap_a = any(all(swap_details.metric3_swaps_a == [i_idx, j_idx], 2)) || any(all(swap_details.metric3_swaps_a == [j_idx, i_idx], 2));
            swap_flag = sprintf('%d', is_swap_a);
            
            % Check for neutrality in Metric 1
            metric1_neutral = ~all_sig_matrices{1}(i_idx, j_idx) && ~all_sig_matrices{1}(j_idx, i_idx);
            % Check for significance in Metric 2 (our tie-breaker)
            is_metric2_sig = all_sig_matrices{2}(i_idx, j_idx) || all_sig_matrices{2}(j_idx, i_idx);
            
            log_this_entry = false;
            if is_swap_a
                % Log if a swap was made because M1 was neutral and M2 was significant
                log_this_entry = true; sort_key = 1;
                reason = sprintf(lang.output.log.reason_m3a, metric_names{1}); % M1 was neutral
            elseif metric1_neutral && is_metric2_sig
                % Log if M1 was neutral and M2 was sig, but order was already correct
                log_this_entry = true; sort_key = 2;
                reason = lang.output.log.reason_m3_correct;
            elseif metric1_neutral && ~is_metric2_sig
                % Log if M1 was neutral and M2 was also neutral
                log_this_entry = true; sort_key = 3;
                reason = sprintf(lang.output.log.reason_m3_no_win, metric_names{2});
            end
            
            if log_this_entry
                log_idx = log_idx + 1;
                alpha_k = all_alpha_matrices{2}(i_idx, j_idx); % Alpha from M2
                signif_str = sprintf('%d', is_metric2_sig);
                pair_row_idx = find((pair_idx_all(:, 1) == i_idx & pair_idx_all(:, 2) == j_idx) | ...
                    (pair_idx_all(:, 1) == j_idx & pair_idx_all(:, 2) == i_idx), 1);
                if isempty(pair_row_idx), warning(lang.warnings.pair_index_not_found, name_i, name_j); pair_row_idx=1; end
                ci_d_str = sprintf('[%+.3f, %+.3f]', ci_d_all(pair_row_idx, 1, 2), ci_d_all(pair_row_idx, 2, 2)); % CIs from M2
                ci_r_str = sprintf('[%.3f, %.3f]', ci_r_all(pair_row_idx, 1, 2), ci_r_all(pair_row_idx, 2, 2)); % CIs from M2
                
                log_entry = {name_i, name_j, metric_names{2}, p_obs, alpha_k, d_obs, ci_d_str, d_thresh(2), rel_obs, ci_r_str, rel_thresh(2), signif_str};
                if power_analysis_done
                    power_val = results.power_results.power_matrices{2}(k);
                    log_entry{end+1} = sprintf('%.1f', power_val * 100);
                end
                log_entry = [log_entry, {swap_flag}, {reason}, {sort_key}];
                log_data(log_idx, :) = log_entry;
            end
        end
    
    elseif (num_metrics == 3 && strcmp(ranking_mode, 'M1_M2_M3'))
        % M1_M2_M3 Logic Log
        for k = 1:m
            % Extract results for the current pair.
            res = swap_details.results_metric3(k, :);
            name_i = res{1}; name_j = res{2}; p_obs = res{3}; d_obs = res{4}; rel_obs = res{5}; i_idx = res{6}; j_idx = res{7};
            % Check if this pair was swapped due to Metric 3's logic A or B.
            is_swap_a = any(all(swap_details.metric3_swaps_a == [i_idx, j_idx], 2)) || any(all(swap_details.metric3_swaps_a == [j_idx, i_idx], 2));
            is_swap_b = any(all(swap_details.metric3_swaps_b == [i_idx, j_idx], 2)) || any(all(swap_details.metric3_swaps_b == [j_idx, i_idx], 2));
            % Check for neutrality (non-significance) in the preceding metrics.
            metric1_neutral = ~all_sig_matrices{1}(i_idx, j_idx) && ~all_sig_matrices{1}(j_idx, i_idx);
            metric2_neutral = ~all_sig_matrices{2}(i_idx, j_idx) && ~all_sig_matrices{2}(j_idx, i_idx);
            is_a_peer_comparison = metric1_neutral && metric2_neutral;
            % Check for significance in Metric 3 for this pair.
            is_metric3_sig = all_sig_matrices{3}(i_idx, j_idx) || all_sig_matrices{3}(j_idx, i_idx);
            % Initialize flags to decide if an entry should be logged for this comparison.
            log_this_entry = false; swap_flag = '0';
            
            % This block determines if a comparison is relevant for the Metric 3 log and assigns a reason.
            if is_swap_a
                % Log if a swap was made because M2 was neutral and M3 was significant (Logic A).
                log_this_entry = true; sort_key = 3; swap_flag = '1';
                reason = sprintf(lang.output.log.reason_m3a, metric_names{2});
            elseif is_swap_b
                % Log if a swap was made because M1 and M2 were neutral and M3 was significant (Logic B).
                log_this_entry = true; sort_key = 4; swap_flag = '1';
                reason = sprintf(lang.output.log.reason_m3b, metric_names{1}, metric_names{2});
            elseif is_metric3_sig && ~metric2_neutral
                % Log if M3 was significant but a M2 comparison took precedence.
                log_this_entry = true; sort_key = 5;
                reason = sprintf(lang.output.log.reason_m3_no_comp1, metric_names{2});
            elseif is_metric3_sig && ~metric1_neutral
                % Log if M3 was significant but a M1 comparison took precedence.
                log_this_entry = true; sort_key = 7;
                reason = sprintf(lang.output.log.reason_m3_no_comp2, metric_names{1});
            elseif is_a_peer_comparison
                % Log if the pair was tied on M1 and M2, making this a "peer comparison".
                log_this_entry = true;
                if is_metric3_sig
                    % The tie was broken by a significant M3 result.
                    sort_key = 8; winner = i_idx;
                    if d_obs < 0, winner = j_idx; end
                    % Check if the final rank reflects the M3 outcome.
                    if point_estimate_ranks(winner) < point_estimate_ranks(res{7 - res{6} / i_idx}),...
                            reason = lang.output.log.reason_m3_correct;
                    else
                        reason = lang.output.log.reason_m3_no_direct;
                    end
                else
                    % The pair remains tied as M3 was also not significant.
                    sort_key = 11; reason = sprintf(lang.output.log.reason_m3_no_win, metric_names{3});
                end
            elseif metric2_neutral
                % Log other cases where M2 was neutral, so M3 was consulted.
                log_this_entry = true;
                if is_metric3_sig
                    % M3 was significant and confirmed the existing rank order from M1.
                    sort_key = 6; reason = lang.output.log.reason_m3_correct;
                else
                    % M3 was not significant, so no change was made.
                    sort_key = 10; reason = sprintf(lang.output.log.reason_m3_no_win, metric_names{3});
                end
            end
    
            % If the comparison was deemed relevant, create and store the log entry.
            if log_this_entry
                log_idx = log_idx + 1;
                alpha_k = all_alpha_matrices{3}(i_idx, j_idx);
                signif_str = sprintf('%d', is_metric3_sig);
                
                % Find the corresponding row index for confidence intervals.
                pair_row_idx = find((pair_idx_all(:, 1) == i_idx & pair_idx_all(:, 2) == j_idx) | ...
                    (pair_idx_all(:, 1) == j_idx & pair_idx_all(:, 2) == i_idx), 1);
                if isempty(pair_row_idx), warning(lang.warnings.pair_index_not_found, name_i, name_j); pair_row_idx=1; end
    
                % Format confidence interval strings.
                ci_d_str = sprintf('[%+.3f, %+.3f]', ci_d_all(pair_row_idx, 1, 3), ci_d_all(pair_row_idx, 2, 3));
                ci_r_str = sprintf('[%.3f, %.3f]', ci_r_all(pair_row_idx, 1, 3), ci_r_all(pair_row_idx, 2, 3));
                
                % Assemble the log entry.
                log_entry = {name_i, name_j, metric_names{3}, p_obs, alpha_k, d_obs, ci_d_str, d_thresh(3), rel_obs, ci_r_str, rel_thresh(3), signif_str};
                if power_analysis_done
                    power_val = results.power_results.power_matrices{3}(k);
                    log_entry{end+1} = sprintf('%.1f', power_val * 100);
                end
                % Append the swap flag, the reason, and the sort key.
                log_entry = [log_entry, {swap_flag}, {reason}, {sort_key}];
                log_data(log_idx, :) = log_entry;
            end
        end
    end % End of conditional M3 log
    
    % Trim the pre-allocated cell array to remove any unused rows.
    log_data = log_data(1:log_idx, :);
    
    % After collecting all data, sort it for a structured and readable output.
    if ~isempty(log_data)
        % Create a primary sort key based on the metric (1, 2, or 3).
        test_types = zeros(size(log_data, 1), 1);
        for i = 1:size(log_data, 1)
            testname = log_data{i, 3}; % Get metric name from the third column.
            % Dynamic assignment of sort key
            if strcmp(testname, metric_names{1}), test_types(i) = 1;
            elseif num_metrics >= 2 && strcmp(testname, metric_names{2}), test_types(i) = 2;
            elseif num_metrics == 3 && strcmp(testname, metric_names{3}), test_types(i) = 3;
            else, test_types(i) = 99; % Should not happen
            end
        end
        
        % Create a sorting matrix with the primary key (metric type) and secondary key (internal logic).
        sort_key_col_idx = size(log_data, 2);
        sort_matrix = [test_types, cell2mat(log_data(:, sort_key_col_idx))];
        [~, sorted_indices] = sortrows(sort_matrix);
        
        % Apply the sorted order to the log data.
        sorted_log_data = log_data(sorted_indices, :);
        
        % Restructure the header for a more intuitive console display (e.g., "Dataset A vs Dataset B").
        header_line2_parts = {lang.csv.headers.dataset1, 'vs', lang.csv.headers.dataset2, header_parts_log{3:end}};
        
        % Prepare a new cell array for console output with formatted numerical values.
        table_data_console = cell(size(sorted_log_data, 1), numel(header_line2_parts));
        for i = 1:size(sorted_log_data, 1)
            row_data_orig = sorted_log_data(i, 1:(end-1)); % Exclude the internal sort key.
            
            % Populate the first three columns: Dataset1, vs, Dataset2.
            table_data_console{i, 1} = row_data_orig{1};
            table_data_console{i, 2} = 'vs';
            table_data_console{i, 3} = row_data_orig{2};
            
            % Copy and format the remaining columns with appropriate precision.
            if power_analysis_done
                table_data_console{i, 4} = row_data_orig{3};   % Metric
                table_data_console{i, 5} = sprintf('%.4f', row_data_orig{4}); % p
                table_data_console{i, 6} = sprintf('%.4f', row_data_orig{5}); % alpha
                table_data_console{i, 7} = sprintf('%+.3f', row_data_orig{6}); % d
                table_data_console{i, 8} = row_data_orig{7};   % CI(d)
                table_data_console{i, 9} = sprintf('%.3f', row_data_orig{8}); % d_thresh
                table_data_console{i, 10} = sprintf('%.3f', row_data_orig{9}); % rel
                table_data_console{i, 11} = row_data_orig{10};  % CI(rel)
                table_data_console{i, 12} = sprintf('%.3f', row_data_orig{11}); % rel_thresh
                table_data_console{i, 13} = char(row_data_orig{12}); % Sig
                table_data_console{i, 14} = char(row_data_orig{13}); % Power
                table_data_console{i, 15} = char(row_data_orig{14}); % Swap
                table_data_console{i, 16} = char(row_data_orig{15}); % Reason
            else
                table_data_console{i, 4} = row_data_orig{3};   % Metric
                table_data_console{i, 5} = sprintf('%.4f', row_data_orig{4}); % p
                table_data_console{i, 6} = sprintf('%.4f', row_data_orig{5}); % alpha
                table_data_console{i, 7} = sprintf('%+.3f', row_data_orig{6}); % d
                table_data_console{i, 8} = row_data_orig{7};   % CI(d)
                table_data_console{i, 9} = sprintf('%.3f', row_data_orig{8}); % d_thresh
                table_data_console{i, 10} = sprintf('%.3f', row_data_orig{9}); % rel
                table_data_console{i, 11} = row_data_orig{10};  % CI(rel)
                table_data_console{i, 12} = sprintf('%.3f', row_data_orig{11}); % rel_thresh
                table_data_console{i, 13} = char(row_data_orig{12}); % Sig
                table_data_console{i, 14} = char(row_data_orig{13}); % Swap
                table_data_console{i, 15} = char(row_data_orig{14}); % Reason
            end
        end
    
        % Dynamically calculate column widths for a clean console table alignment.
        col_widths = cellfun(@strlength, header_line2_parts);
        for r = 1:size(table_data_console, 1)
            for c = 1:size(table_data_console, 2)
                col_widths(c) = max(col_widths(c), strlength(table_data_console{r, c}));
            end
        end
        col_widths = col_widths + 2; 
        
        % Define text alignment for each column in the console output.
        num_cols_console = numel(header_line2_parts);
        alignments = repmat({'c'}, 1, num_cols_console); % Default to centered.
        alignments{end} = 'l'; % Left-align the final 'Reason' column for readability.
    
        % Print the formatted header to the console.
        header_line = strjoin(arrayfun(@(c) format_text(header_line2_parts{c}, col_widths(c), alignments{c}), 1:num_cols_console, 'UniformOutput', false), '|');
        fprintf('%s\n', header_line);
        fprintf('%s\n', repmat('-', 1, strlength(header_line)));
        
        % Loop through the sorted data to output both to the CONSOLE and the CSV FILE.
        for i = 1:size(sorted_log_data, 1)
            % Print the pre-formatted row for display.
            row_line_console = strjoin(arrayfun(@(c) format_text(table_data_console{i, c}, col_widths(c), alignments{c}), 1:num_cols_console, ...
                'UniformOutput', false), '|');
            fprintf('%s\n', row_line_console);
    
            % Prepare the data row for the file, ensuring proper formatting and quoting.
            row_data = sorted_log_data(i, 1:(end-1)); % Get original data (without sort key).
            csv_cells = cell(1, numel(row_data));
            for j = 1:numel(row_data)
                cell_content = row_data{j};
                if isnumeric(cell_content)
                    % Apply specific formatting for numeric types.
                    if any(j == [4, 5]), csv_cells{j} = sprintf('%.4f', cell_content);
                    elseif any(j == [8, 9, 11]), csv_cells{j} = sprintf('%.3f', cell_content);
                    elseif j == 6, csv_cells{j} = sprintf('%+.3f', cell_content);
                    else, csv_cells{j} = num2str(cell_content);
                    end
                else
                    % Enclose string data in double quotes.
                    csv_cells{j} = ['"' strtrim(char(cell_content)) '"'];
                end
            end
            % Write the semicolon-separated row to the log file.
            fprintf(fid_log, '%s\n', strjoin(csv_cells, ';'));
        end
        % Print a final separating line to the console.
        fprintf('%s\n', repmat('-', 1, strlength(header_line)));
    end
    % Close the log file.
    fclose(fid_log);
    % Notify the user that the log file has been saved successfully.
    fprintf(['\n' lang.output.files.log_file_saved '\n'], csv_filename_log);

catch ME
    % Safety cleanup: Ensure file is closed if an error occurs
    if exist('fid_log', 'var') && fid_log ~= -1
        fclose(fid_log); 
    end
    fprintf([lang.errors.file_save_error '\n'], ME.message);
end
end

%% Helper Functions
function output = format_text(text, width, alignment)
    % Width without the buffer padding of 2, as this is only relevant for 'c'
    text_len = strlength(text);
    
    switch alignment
        case 'r'
            % Right-aligned: Uses the total width 'width' minus the text length.
            % The 2 padding characters are included here but displayed as spaces to the left of the text.
            output = [repmat(' ', 1, width - text_len), text];
        case 'l'
            % Left-aligned: Uses the total width 'width' minus the text length.
            % The 2 padding characters are included here, displayed as spaces to the right of the text.
            output = [text, repmat(' ', 1, width - text_len)];
        otherwise % 'c' (Centered)
            % Centered: Recalculates padding to center the text.
            padding = width - text_len;
            padding_left = floor(padding / 2);
            padding_right = ceil(padding / 2);
            output = [repmat(' ', 1, padding_left), text, repmat(' ', 1, padding_right)];
    end
end