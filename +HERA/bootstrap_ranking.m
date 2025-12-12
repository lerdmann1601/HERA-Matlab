function [final_bootstrap_ranks, selected_B_final, stability_data_rank, h_figs_rank, h_fig_hist_rank, ci_lower_rank, ci_upper_rank] = bootstrap_ranking(all_data, ...
          thresholds, config, dataset_names, final_rank, pair_idx_all, num_probanden, graphics_dir, csv_dir, manual_B, s, styles, lang, base_name)
% BOOTSTRAP_RANKING - Performs a cluster bootstrap analysis to determine rank confidence intervals.
%
% Syntax:
%   [final_bootstrap_ranks, selected_B_final, stability_data_rank, h_figs_rank, h_fig_hist_rank, ci_lower_rank, ci_upper_rank] = bootstrap_ranking(all_data, ...
%    thresholds, config, dataset_names, final_rank, pair_idx_all, num_probanden, graphics_dir, csv_dir, manual_B, s, styles, lang, base_name)
%
% Description:
%   This function evaluates the stability of the dataset ranking using a cluster bootstrap approach.
%   It employs a memory-optimized implementation (virtual indexing) to avoid data duplication during resampling.
%   Subjects (as clusters) are drawn with replacement, and for each sample, the complete ranking process (including effect size calculation) is repeated. 
%   This process generates a distribution of possible ranks for each dataset, which is then used to assess the robustness of the primary ranking result.
%   This function dynamically handles 1, 2, or 3 metrics and respects the 'ranking_mode' from the 'config' struct when calling 'calculate_ranking'.
%
% Workflow:
%   1.  Dynamic search for the optimal number of bootstrap samples (B): 
%       Analyzes the stability of the rank confidence intervals over multiple trials and an increasing number of B-values. 
%       Stability is quantified using the ratio of the Interquartile Range (IQR) to the median of the confidence interval widths.
%   2.  Selection of the final B-value ('selected_B_final'): 
%       The optimal B is determined either when stability converges (i.e., the relative improvement falls below a tolerance) - 
%       or, if convergence is not reached, through an "elbow analysis" of the stability curve.
%   3.  Visualization (Convergence): 
%       Generates and saves the convergence plot showing the stability analysis using 'HERA.plot.rank_convergence'.
%   4.  Final bootstrap analysis: 
%       Using the optimal B-value, the function runs a final, comprehensive bootstrap analysis to generate the definitive rank distribution for each dataset.
%   5.  Result Output (Console): 
%       Calculates the frequency of each rank for each dataset and prints a formatted table to the console.
%   7.  Result Output (csv): 
%       Saves the detailed rank distribution (counts and percentages) for each dataset to a CSV file.
%   8.  Visualization (Distribution): 
%       Generates and saves histograms displaying the final rank distribution for each dataset using 'HERA.plot.rank_histograms'.
%   9.  CI Calculation:
%       Computes the user-defined confidence intervals (e.g., 95%) for the final ranks based on the bootstrap distribution.
%
% Inputs:
%   all_data        - Cell array containing the original data matrices for each metric (1, 2, or 3 cells).
%   thresholds      - Struct with the pre-calculated effect size thresholds (thresholds.d_thresh, thresholds.rel_thresh).
%   config          - Global configuration struct. Uses:
%                       .bootstrap_ranks (settings)
%                       .ranking_mode (logic)
%                       .ci_level (confidence interval width)
%   dataset_names   - Cell array of strings with the names of the datasets (required for the 'calculate_ranking' function).
%   final_rank      - Vector of the primary ranking (used for sorting the console output).
%   pair_idx_all    - Matrix of indices for all pairwise comparisons between datasets.
%   num_probanden   - Scalar, the total number of subjects (clusters) in the study.
%   graphics_dir    - Path to the directory where generated plots will be saved.
%   csv_dir         - Path to the directory where CSV outputs will be saved.
%   manual_B        - (Optional) A manually specified number of bootstrap repetitions (B). If provided, the dynamic search is skipped.
%   s               - The random number generator stream to ensure reproducibility.
%   styles          - Struct containing settings for graphical design (colors, fonts, etc.).
%   lang            - Language pack struct loaded from a JSON file for all text outputs.
%   base_name       - Base name for output files (e.g., from the log timestamp).
%
% Outputs:
%   final_bootstrap_ranks - A matrix of size [num_datasets, selected_B_final] containing the rank distribution from the final analysis. 
%                           Each column represents a full ranking result from one bootstrap iteration.
%   selected_B_final      - The optimally determined number of bootstrap samples used for the final analysis.
%   stability_data_rank   - Struct with convergence curve data for JSON export.
%   h_figs_rank           - A vector of figure handles for the convergence graphic.
%   h_fig_hist_rank       - The handle of the histogram graphic showing the final rank distributions.
%   ci_lower_rank         - [num_datasets x 1] Lower bound of the 95% Confidence Interval for the rank.
%   ci_upper_rank         - [num_datasets x 1] Upper bound of the 95% Confidence Interval for the rank.
%
% Author: Lukas von Erdmannsdorff

arguments
    all_data (1,:) cell
    thresholds (1,1) struct
    config (1,1) struct
    dataset_names (1,:) cell
    final_rank
    pair_idx_all
    num_probanden (1,1) double
    graphics_dir (1,1) string
    csv_dir (1,1) string
    manual_B
    s
    styles (1,1) struct
    lang (1,1) struct
    base_name (1,1) string
end

%% 1. Initialization and Convergence Check for Optimal B
% Initialization of local variables from input parameters.
ts = config.timestamp;
num_datasets_b = size(all_data{1}, 2);
n_subj_b = num_probanden;
num_metrics = numel(all_data); % Get dynamic number of metrics

% Create a dedicated subfolder for the ranking stability plots.
subfolder_ranking = fullfile(graphics_dir, 'Ranking');
if ~exist(subfolder_ranking, 'dir')
    mkdir(subfolder_ranking);
end
% Initialize output arrays for figure handles to prevent errors if no plots are generated.
h_figs_rank = gobjects(0);
h_fig_hist_rank = gobjects(0); 

% Check if a manual B-value was provided by the user.
if ~isempty(manual_B)
    % If so, use the provided value and skip the dynamic search.
    selected_B_final = manual_B;
    stability_data_rank = []; % No stability analysis performed
    fprintf(['\n' lang.ranking.manual_b_info '\n'], selected_B_final);
    % Skip directly to the final calculation in Section 3.
else
    % If no manual B is given, perform the dynamic search for the optimal B.
    cfg_rank = config.bootstrap_ranks;
    
    % Inform the user about the start of the search and the criteria being used.
    fprintf(['\n' lang.ranking.searching_optimal_b '\n']);
    
    % Initialize vectors for the stability analysis.
    B_vector_b = cfg_rank.B_start:cfg_rank.B_step:cfg_rank.B_end; 
    stability_vector_b = NaN(1, numel(B_vector_b)); 
    final_b_idx = 0; 
    converged = false; 
    elbow_idx_rank = [];
    
    % Loop over the different B-values to check for stability.
    for b_idx = 1:numel(B_vector_b)
        Br_b = B_vector_b(b_idx);
        fprintf([' -> ' lang.ranking.checking_stability '\n'], Br_b, cfg_rank.n_trials);
        ci_widths_b = zeros(cfg_rank.n_trials, num_datasets_b);
        
        % Perform n_trials to check the stability of the rank confidence intervals for the current B-value.
        parfor t_b = 1:cfg_rank.n_trials
            % Each parallel worker gets its own reproducible substream of the random number generator.
            s_worker = s; % Workaround for passing the stream object to parfor.
            s_worker.Substream = t_b;
            rank_tmp_b = zeros(num_datasets_b, Br_b);
            
            % 1. Perform a cluster bootstrap: Subjects (clusters) are drawn with replacement.
            % Vectorized generation of indices for all Br_b iterations in this trial [N x Br_b]
            boot_indices_block = randi(s_worker, n_subj_b, [n_subj_b, Br_b]);

            % 2. Calculate Effect Sizes (Vectorized or Loop-based)
            % We need to fill: bootstrap_d_vals_all [num_pairs, num_metrics, Br_b]
            num_pairs = size(pair_idx_all, 1);
            
            % Check for NaNs globally to decide on strategy
            % If any used data has NaNs, we must use the robust loop.
            % Optimization: Check once per dataset column used.
            has_nans = false;
            for col = 1:num_datasets_b
               for mid = 1:num_metrics
                   if any(isnan(all_data{mid}(:, col)))
                       has_nans = true; break;
                   end
               end
               if has_nans, break; end
            end

            % Pre-allocate 3D arrays for this trial block
            d_vals_3d = zeros(num_pairs, num_metrics, Br_b);
            rel_vals_3d = zeros(num_pairs, num_metrics, Br_b);
            
            if ~has_nans
                 % --- FAST PATH: Vectorized Calculation ---
                 % Process all Br_b bootstrap samples simultaneously per pair.
                 % This is significantly faster as cliffs_delta processes [N x Br_b] matrices.
                 for p_idx = 1:num_pairs
                    i = pair_idx_all(p_idx, 1);
                    j = pair_idx_all(p_idx, 2);
                    for metric_idx = 1:num_metrics
                         % Extract matrices [N x Br_b]
                         data_i = all_data{metric_idx}(:, i);
                         data_j = all_data{metric_idx}(:, j);
                         x_mat = data_i(boot_indices_block);
                         y_mat = data_j(boot_indices_block);
                         
                         % Vectorized calls
                         d_vals_3d(p_idx, metric_idx, :) = HERA.stats.cliffs_delta(x_mat, y_mat);
                         rel_vals_3d(p_idx, metric_idx, :) = HERA.stats.relative_difference(x_mat, y_mat);
                    end
                 end
            else
                 % --- ROBUST PATH: Loop Calculation ---
                 % Fallback for data with NaNs to ensure correct pairwise exclusion.
                 for bb = 1:Br_b
                     current_indices = boot_indices_block(:, bb);
                     for p_idx = 1:num_pairs
                        i = pair_idx_all(p_idx, 1);
                        j = pair_idx_all(p_idx, 2);
                        for metric_idx = 1:num_metrics
                            % Standard single-vector extraction
                            col_i = all_data{metric_idx}(current_indices, i);
                            col_j = all_data{metric_idx}(current_indices, j);
                            
                            valid = ~isnan(col_i) & ~isnan(col_j);
                            x = col_i(valid); y = col_j(valid);
                            
                            if ~isempty(x)
                                d_vals_3d(p_idx, metric_idx, bb) = HERA.stats.cliffs_delta(x, y);
                                rel_vals_3d(p_idx, metric_idx, bb) = HERA.stats.relative_difference(x, y);
                            else
                                d_vals_3d(p_idx, metric_idx, bb) = NaN;
                                rel_vals_3d(p_idx, metric_idx, bb) = NaN;
                            end
                        end
                     end
                 end
            end

            % 3. & 4. Ranking (Iterative)
            % Now run the lightweight ranking logic on the pre-calculated effect sizes
            for bb = 1:Br_b
                % Slice 2D effect sizes for this specific iteration
                es_struct = struct('d_vals_all', d_vals_3d(:, :, bb), ...
                                   'rel_vals_all', rel_vals_3d(:, :, bb));
                
                % Ranking call (using original data view + current indices for M1 means calc)
                [~, bootstrap_rank] = HERA.calculate_ranking(...
                     all_data, es_struct, thresholds, config, dataset_names, pair_idx_all, boot_indices_block(:, bb));
                
                rank_tmp_b(:, bb) = bootstrap_rank;
            end
            
            % After Br_b iterations, calculate the width of the 95% confidence interval of the ranks for this one trial.
            ci_bounds_b = quantile(rank_tmp_b, [0.025, 0.975], 2);
            ci_widths_b(t_b, :) = ci_bounds_b(:, 2) - ci_bounds_b(:, 1);
        end
        
        % Calculate the stability metric for the current B-value. This is the median of the relative IQR of CI widths across all datasets.
        med_widths = median(ci_widths_b, 1);
        iqr_widths = iqr(ci_widths_b, 1);
        
        % The final stability value is the maximum of the IQR widths (absolute stability).
        stability_vector_b(b_idx) = max(iqr_widths, [], 'all', 'omitnan');
        final_b_idx = b_idx; % Keep track of the last executed index.
        
        % Check convergence using the reusable stats function
        % stability_vector_b(1:b_idx) contains the history up to current step
        [converged, results] = HERA.stats.check_convergence(stability_vector_b(1:b_idx), cfg_rank);

        % Log details from the check
        if ~isnan(results.improvement)
             if abs(results.improvement) < cfg_rank.convergence_tolerance
                fprintf(['    ' lang.ranking.convergence_run_info '\n'], results.improvement * 100, results.streak, cfg_rank.convergence_streak_needed);
             else
                fprintf(['    ' lang.ranking.stability_change_info '\n'], results.improvement * 100);
             end
        end

        if converged
            if isfield(cfg_rank, 'convergence_streak_needed') && ~isempty(cfg_rank.convergence_streak_needed) && cfg_rank.convergence_streak_needed > 0
                fprintf([lang.ranking.convergence_reached '\n'], results.improvement * 100);
                fprintf([lang.ranking.stable_runs_info '\n'], cfg_rank.convergence_streak_needed);
            else
                fprintf([lang.ranking.convergence_reached '\n'], results.improvement * 100);
            end
            break; % Abort the stability check loop as soon as convergence is reached.
        end
        
    end % End of for-loop for stability check.
        
    % Selection of the final B-value based on the analysis outcome.
    B_tested_vector_b = B_vector_b(1:final_b_idx);
    stability_vector_b_plotted = stability_vector_b(1:final_b_idx);
        
    if converged
        % If the process converged, the last tested B-value is selected as the optimal one.
        selected_B_final = B_tested_vector_b(end);
        fprintf([lang.ranking.convergence_result '\n'], selected_B_final);
    else
        % If no convergence was reached, perform an "elbow analysis" as a fallback.
        fprintf([' ' lang.ranking.elbow_analysis_info '\n']);
        
        % Use the generic elbow finder
        % It handles the 1-based index return and flat line cases internally
        [~, elbow_idx_rank] = HERA.stats.find_elbow_point(B_tested_vector_b, stability_vector_b_plotted);
        
        selected_B_final = B_tested_vector_b(elbow_idx_rank);
        fprintf([lang.ranking.elbow_result '\n'], selected_B_final);
    end
    
    % Store stability data for JSON export
    stability_data_rank = struct();
    stability_data_rank.B_vector = B_tested_vector_b;
    stability_data_rank.global_stability = stability_vector_b_plotted;
    stability_data_rank.converged = converged;
    stability_data_rank.elbow_indices = elbow_idx_rank;
    % Rank stability only has one (global) curve, so detailed_stability is empty.
    stability_data_rank.detailed_stability = [];
    
%% 2. Create and save the convergence graphic.
    h_fig_rank = HERA.plot.rank_convergence(B_tested_vector_b, stability_vector_b_plotted, ...
                                            selected_B_final, converged, elbow_idx_rank, ...
                                            config, styles, lang);

    if isgraphics(h_fig_rank)
         % Save graphic.
        [~, fName, fExt] = fileparts(lang.files.convergence_rank_stability);
        filename = fullfile(subfolder_ranking, [fName, '_', ts, fExt]);
        exportgraphics(h_fig_rank, filename, 'Resolution', 300, 'Padding', 30);
        fprintf([lang.ranking.convergence_plot_saved '\n'], filename);

        h_figs_rank(end+1) = h_fig_rank;
    end
end

%% 3. Final Bootstrap Analysis with the Optimal Number of Repetitions
% This section performs the main bootstrap process with the final, optimally determined B-value.
% We use Batched Processing to enable Vectorized Effect Size calculation (Fast Path).
% Instead of 1:B, we iterate over batches to balance memory usage and vectorization speed.

weights = ones(1, selected_B_final); % Used just for counting
% Dynamic Batch Sizing
% We prioritize saturating the parallel workers over maximizing per-core vectorization.
% If B is small (e.g. 200) and we have 10 cores, a fixed batch of 100 leaves 8 cores idle.
% We divide the work to ensure all workers are busy.

pool = gcp('nocreate');
if isempty(pool)
    % If no pool is running, estimate based on physical cores
    % (Matlab will likely start one with this many workers)
    num_workers = feature('numcores');
else
    num_workers = pool.NumWorkers;
end

% "Sweet Spot" from benchmark was 100, but we shouldn't exceed (B / Workers).
% We also enforce a minimum batch size (e.g., 10) to avoid loop overhead dominance.
MIN_BATCH = 10;
OPTIMAL_BATCH_LIMIT = 100; % Don't go larger than this even if B is huge

% Calculate split for parallelism
calc_batch = ceil(selected_B_final / num_workers);

% Apply constraints: 
% 1. Must be at least MIN_BATCH (unless B < MIN_BATCH)
% 2. Should ideally be around calc_batch to use all cores.
% 3. Capped at OPTIMAL_BATCH_LIMIT because larger batches hurt cache.

BATCH_SIZE = max(MIN_BATCH, min(calc_batch, OPTIMAL_BATCH_LIMIT));

% Edge case: if B is very small, minimal batch is total B
if selected_B_final < BATCH_SIZE
    BATCH_SIZE = selected_B_final;
end

num_batches = ceil(selected_B_final / BATCH_SIZE);

% Pre-allocate results [NumDatasets x TotalB]
% Note: Parfor cannot write to sliced indices of a variable easily if indices are procedural.
% We collect results in a cell array of batches or linear index strategy.
% Strategy: Parfor over batches, return a matrix block, then combine.

rank_batches = cell(1, num_batches);

parfor b_idx = 1:num_batches
    s_worker = s;
    s_worker.Substream = b_idx + 1000; % Offset to avoid overlap with stability check
    
    % Determine range for this batch
    start_idx = (b_idx - 1) * BATCH_SIZE + 1;
    end_idx = min(b_idx * BATCH_SIZE, selected_B_final);
    current_batch_size = end_idx - start_idx + 1;
    
    % --- 1. Generate Indices [N x BatchSize] ---
    boot_indices_block = randi(s_worker, n_subj_b, [n_subj_b, current_batch_size]);
    
    % --- 2. Calculate Effect Sizes (Hybrid Vectorized/Loop) ---
    num_pairs = size(pair_idx_all, 1);
    
    % Check for NaNs globally (Optimization: Check first column of first metric as proxy or check all)
    % A complete check is safer.
    has_nans = false;
    for col = 1:num_datasets_b
       for mid = 1:num_metrics
           if any(isnan(all_data{mid}(:, col)))
               has_nans = true; break;
           end
       end
       if has_nans, break; end
    end
    
    d_vals_3d = zeros(num_pairs, num_metrics, current_batch_size);
    rel_vals_3d = zeros(num_pairs, num_metrics, current_batch_size);
    
    if ~has_nans
        % FAST PATH
        for p_idx = 1:num_pairs
            i = pair_idx_all(p_idx, 1);
            j = pair_idx_all(p_idx, 2);
            for metric_idx = 1:num_metrics
                 data_i = all_data{metric_idx}(:, i);
                 data_j = all_data{metric_idx}(:, j);
                 x_mat = data_i(boot_indices_block);
                 y_mat = data_j(boot_indices_block);
                 
                 d_vals_3d(p_idx, metric_idx, :) = HERA.stats.cliffs_delta(x_mat, y_mat);
                 rel_vals_3d(p_idx, metric_idx, :) = HERA.stats.relative_difference(x_mat, y_mat);
            end
        end
    else
        % ROBUST PATH
        for k = 1:current_batch_size
             current_indices = boot_indices_block(:, k);
             for p_idx = 1:num_pairs
                i = pair_idx_all(p_idx, 1);
                j = pair_idx_all(p_idx, 2);
                for metric_idx = 1:num_metrics
                    col_i = all_data{metric_idx}(current_indices, i);
                    col_j = all_data{metric_idx}(current_indices, j);
                    valid = ~isnan(col_i) & ~isnan(col_j);
                    x = col_i(valid); y = col_j(valid);
                    
                    if ~isempty(x)
                        d_vals_3d(p_idx, metric_idx, k) = HERA.stats.cliffs_delta(x, y);
                        rel_vals_3d(p_idx, metric_idx, k) = HERA.stats.relative_difference(x, y);
                    else
                        d_vals_3d(p_idx, metric_idx, k) = NaN;
                        rel_vals_3d(p_idx, metric_idx, k) = NaN;
                    end
                end
             end
        end
    end
    
    % --- 3. Ranking Loop ---
    batch_ranks = zeros(num_datasets_b, current_batch_size);
    for k = 1:current_batch_size
        es_struct = struct('d_vals_all', d_vals_3d(:, :, k), ...
                           'rel_vals_all', rel_vals_3d(:, :, k));
                       
        [~, rnk] = HERA.calculate_ranking(...
             all_data, es_struct, thresholds, config, dataset_names, pair_idx_all, boot_indices_block(:, k));
        batch_ranks(:, k) = rnk;
    end
    
    rank_batches{b_idx} = batch_ranks;
end

% Flatten results
final_bootstrap_ranks = [rank_batches{:}];

%% 4. Print Bootstrap Rank Distribution to Console
% This section calculates the rank distribution and prints a formatted table to the console, sorted by the final rank.

fprintf(lang.ranking.distribution_header, selected_B_final);

% Get the sorting order based on the final rank
[~, sort_idx] = sort(final_rank);

% Define headers and prepare for dynamic formatting
header_parts = {lang.ranking.table_rank, lang.ranking.table_dataset, lang.ranking.table_dist};
alignments = {'c', 'l', 'l'}; % Center, Left, Left
table_data = cell(num_datasets_b, 3);

% Loop over datasets in the order of their final rank
for i_sorted = 1:num_datasets_b
    i = sort_idx(i_sorted); % Get the original index
    
    % Get corresponding data
    current_final_rank = final_rank(i);
    dataset_name = dataset_names{i};
    rank_data = final_bootstrap_ranks(i, :);
    
    % Calculate the counts and percentages of unique ranks
    [unique_ranks, ~, group_idx] = unique(rank_data);
    counts = accumarray(group_idx, 1);
    percentages = (counts / selected_B_final) * 100;
    
    % Create the compact distribution string (e.g., "1: 70.0%, 2: 30.0%")
    dist_parts = cell(1, numel(unique_ranks));
    for j = 1:numel(unique_ranks)
        dist_parts{j} = sprintf('%d: %.1f%%', unique_ranks(j), percentages(j));
    end
    distribution_string = strjoin(dist_parts, ', ');
    
    % Store data for the table row
    table_data(i_sorted, :) = {sprintf('%d', current_final_rank), dataset_name, distribution_string};
end

% Calculate dynamic column widths based on content
col_widths = cellfun(@strlength, header_parts);
for r = 1:size(table_data, 1)
    for c = 1:size(table_data, 2)
        col_widths(c) = max(col_widths(c), strlength(table_data{r, c}));
    end
end
col_widths = col_widths + 2; % Add 2 spaces for padding

% Create the header row
header_line_parts = arrayfun(@(c) format_text(header_parts{c}, col_widths(c), alignments{c}), 1:numel(header_parts), 'UniformOutput', false);
header_line = ['|' strjoin(header_line_parts, '|')];

% Create the external continuous separator 
% Its length must match the total width of the other table rows.
total_width = strlength(header_line);
external_separator_line = repmat('-', 1, total_width);

% Print the formatted table
fprintf('%s\n', header_line);
fprintf('%s\n', external_separator_line); % Print continuous separator under header

% Print the table content
for r = 1:size(table_data, 1)
    row_line_parts = arrayfun(@(c) format_text(table_data{r, c}, col_widths(c), alignments{c}), 1:numel(header_parts), 'UniformOutput', false);
    % Create the data row 
    row_line = ['|' strjoin(row_line_parts, '|')];
    fprintf('%s\n', row_line);
    
    % Use the continuous external separator *between* data rows
    if r < size(table_data, 1)
        fprintf('%s\n', external_separator_line);
    end
end
fprintf('%s\n', external_separator_line); % Print continuous bottom border


%% 5. Save Bootstrap Rank Distribution to CSV
% This section calculates the frequency of each rank for each dataset from bootstrap analysis and saves it to a CSV file, sorted by the final rank.

% Inform the user about the CSV saving process.
fprintf(['\n' lang.ranking.saving_csv '\n']);
% Define the output filename.
[~, fName, fExt] = fileparts(lang.files.bootstrap_rank_csv);
fName = strrep(fName, '%s_', '');
csv_filename = fullfile(csv_dir, [fName, '_', ts, fExt]);

try
    % Open the file for writing.
    fileID = fopen(csv_filename, 'w');
    
    % Write the header row to the CSV file.
    header = {lang.ranking.table_rank, ...         
              lang.ranking.table_dataset, ...       
              lang.ranking.csv_bootstrap_rank, ...  
              lang.ranking.csv_frequency_percent, ... 
              lang.ranking.csv_frequency_count};     
    
    % Write the header using comma as delimiter (to match the data fprintf below)
    fprintf(fileID, '%s\n', strjoin(header, ','));
    
    % Get the sorting order based on the final rank
    [~, sort_idx] = sort(final_rank);

    % Loop over each dataset in the order of their final rank
    for i_sorted = 1:num_datasets_b
        i = sort_idx(i_sorted); % Get the original index
        
        % Get the corresponding final rank and dataset name.
        current_final_rank = final_rank(i);
        dataset_name = dataset_names{i};
        
        % Get all bootstrap ranks for the current dataset.
        rank_data = final_bootstrap_ranks(i, :);
        
        % Calculate the counts of unique ranks.
        [unique_ranks, ~, group_idx] = unique(rank_data);
        counts = accumarray(group_idx, 1);
        
        % Convert counts to percentages.
        percentages = (counts / selected_B_final) * 100;
        
        % Write one row in the CSV for each observed rank.
        for j = 1:numel(unique_ranks)
            rank_val = unique_ranks(j);
            count_val = counts(j);
            percent_val = percentages(j);
            
            % Write the data row.
            fprintf(fileID, '%d,"%s",%d,%.2f,%d\n', ...
                current_final_rank, dataset_name, rank_val, percent_val, count_val);
        end
    end
    
    % Close the file.
    fclose(fileID);
    fprintf([lang.ranking.csv_saved '\n'], csv_filename);
    
catch ME
    % In case of an error (e.g., file permissions), close the file and report.
    if exist('fileID', 'var') && fileID ~= -1
        fclose(fileID);
    end
    fprintf('Error saving bootstrap rank CSV: %s\n', ME.message);
end

%% 6. Create and save the histogram distribution of the final ranks
h_fig_hist_rank = HERA.plot.rank_histograms(final_bootstrap_ranks, dataset_names, ...
                                            final_rank, selected_B_final, styles, lang);

if isgraphics(h_fig_hist_rank)
    % Save the complete histogram figure to a file.
    [~, fName, fExt] = fileparts(lang.files.dist_bootstrap_ranks);
    filename = fullfile(subfolder_ranking, [fName, '_', ts, fExt]);
    exportgraphics(h_fig_hist_rank, filename, 'Resolution', 300, 'Padding', 30);
    fprintf([lang.ranking.histogram_plot_saved '\n'], filename);
end

%% 7. Calculate Rank Confidence Intervals
% Calculate Rank Confidence Intervals from the bootstrap distribution
% Respect the user-configured confidence level (e.g., 0.95)
if isfield(config, 'ci_level') && ~isempty(config.ci_level)
    target_ci = config.ci_level;
else
    target_ci = 0.95; % Default
end

alpha = 1 - target_ci;
ci_lower_bound = alpha / 2;
ci_upper_bound = 1 - (alpha / 2);

ci_final_ranks = quantile(final_bootstrap_ranks, [ci_lower_bound, ci_upper_bound], 2);
ci_lower_rank = ci_final_ranks(:, 1);
ci_upper_rank = ci_final_ranks(:, 2);

end

%% Helper Function for table formatting
function output = format_text(text, width, alignment)
    % This function formats a string to a specified width and alignment for console display.
    text_len = strlength(text);
    padding = width - text_len;
    
    switch alignment
        case 'l' % left-aligned
            output = [text, repmat(' ', 1, padding)];
        case 'r' % right-aligned
            output = [repmat(' ', 1, padding), text];
        otherwise % 'c' (centered)
            padding_left = floor(padding / 2);
            padding_right = ceil(padding / 2);
            output = [repmat(' ', 1, padding_left), text, repmat(' ', 1, padding_right)];
    end
end