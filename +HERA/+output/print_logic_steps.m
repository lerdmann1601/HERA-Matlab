function print_logic_steps(results, thresholds, config, shared_info)
% PRINT_LOGIC_STEPS - Visualizes the step-by-step ranking logic to the console.
%
% Syntax:
%   HERA.output.print_logic_steps(results, thresholds, config, shared_info)
%
% Description:
%   Displays the hierarchical ranking process including:
%   - Metric 1 Logic: Comparisons, tied-breaking (Cliff's Delta/Mean), and table of swaps.
%   - Metric 2 Logic (if used): Corrections based on the second metric.
%   - Metric 3 / Logic 3A/3B (if used): Further refinements and final tie-breaking.
%
% Inputs:
%   results     - Struct containing logic details (swaps, intermediate ranks).
%   thresholds  - Struct with thresholds for effect sizes and relative differences.
%   config      - Struct with analysis configuration (ranking_mode, etc.).
%   shared_info - Struct with general information like metric and dataset names.
%
% Outputs:
%   None (prints directly to console).
%
% Author: Lukas von Erdmannsdorff

    import HERA.output.format_text

    % Unpack necessary variables
    lang = shared_info.lang;
    metric_names = shared_info.metric_names;
    dataset_names = shared_info.dataset_names;
    num_datasets = shared_info.num_datasets;
    pair_idx_all = shared_info.pair_idx_all;
    
    num_metrics = numel(metric_names);
    ranking_mode = config.ranking_mode;
    
    mean_metrics = shared_info.mean_metrics;
    
    d_thresh = thresholds.d_thresh;
    rel_thresh = thresholds.rel_thresh;
    
    swap_details = results.swap_details;
    intermediate_orders = results.intermediate_orders;
    final_order = results.final_order;
    
    ci_d_all = results.ci_d_all;
    ci_r_all = results.ci_r_all;

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

end
