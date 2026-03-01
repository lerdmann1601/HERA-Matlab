function save_log(results, thresholds, config, shared_info)
% SAVE_LOG - Generates and saves the detailed log of all comparisons.
%
% Syntax:
%   HERA.output.save_log(results, thresholds, config, shared_info)
%
% Description:
%   Creates a comprehensive log file containing every pairwise comparison, 
%   the statistical results (p-value, Cliff's Delta, etc.), and the specific 
%   reasoning used by the HERA logic for ranking decisions.
%   Saves to a CSV file and prints a formatted version to the console.
%
% Inputs:
%   results     - Struct containing all statistical results and matrices.
%   thresholds  - Struct with comparison thresholds.
%   config      - Struct with configuration (ranking_mode, etc.).
%   shared_info - Struct with general info (names, paths, etc.).
%
% Outputs:
%   None (prints to console, saves CSV).
%
% Author: Lukas von Erdmannsdorff

    import HERA.output.format_text

    % Unpack necessary variables
    lang = shared_info.lang;
    metric_names = shared_info.metric_names;
    csv_dir = shared_info.csv_dir;
    ts = shared_info.config.timestamp;
    pair_idx_all = shared_info.pair_idx_all;
    
    num_metrics = numel(metric_names);
    ranking_mode = config.ranking_mode;
    
    d_thresh = thresholds.d_thresh;
    rel_thresh = thresholds.rel_thresh;
    
    swap_details = results.swap_details;
    point_estimate_ranks = results.final_rank;
    
    all_sig_matrices = results.all_sig_matrices;
    all_alpha_matrices = results.all_alpha_matrices;
    
    ci_d_all = results.ci_d_all;
    ci_r_all = results.ci_r_all;
    
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
        fclose(fid_log); % Close it so writetable can append later

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
            end
            
            % Print a final separating line to the console.
            fprintf('%s\n', repmat('-', 1, strlength(header_line)));
            
            % Accelerated file writing using writetable
            export_cell_array = cell(size(sorted_log_data, 1), num_log_cols);
            for i = 1:size(sorted_log_data, 1)
                row_data = sorted_log_data(i, 1:(end-1)); % original data
                for j = 1:numel(row_data)
                    cell_content = row_data{j};
                    if isnumeric(cell_content)
                        % Apply specific formatting for numeric types.
                        if any(j == [4, 5]), export_cell_array{i, j} = sprintf('%.4f', cell_content);
                        elseif any(j == [8, 9, 11]), export_cell_array{i, j} = sprintf('%.3f', cell_content);
                        elseif j == 6, export_cell_array{i, j} = sprintf('%+.3f', cell_content);
                        else, export_cell_array{i, j} = num2str(cell_content);
                        end
                    else
                        export_cell_array{i, j} = string(strtrim(char(cell_content)));
                    end
                end
            end
            
            var_names = arrayfun(@(x) sprintf('Var%d', x), 1:num_log_cols, 'UniformOutput', false);
            log_table = cell2table(export_cell_array, 'VariableNames', var_names);
            
            writetable(log_table, csv_filename_log, 'Delimiter', ';', 'WriteMode', 'Append', 'WriteVariableNames', false, 'QuoteStrings', true);
        end
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
