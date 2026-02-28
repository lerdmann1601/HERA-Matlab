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
%   Subjects (as clusters) are drawn with replacement, and for each sample, the complete ranking process (including effect size calculation) is repeated. 
%   This process generates a distribution of possible ranks for each dataset, which is then used to assess the robustness of the primary ranking result.
%   This function dynamically handles 1, 2, or 3 metrics and respects the 'ranking_mode' from the 'config' struct when calling 'calculate_ranking'.
%
% Workflow:
%   1.  Dynamic search for the optimal number of bootstrap samples (B): 
%       Analyzes stability over multiple trials using parallelized (`parfor`) execution with 
%       memory-aware batch sizing. Stability is quantified using the maximum Interquartile Range (IQR).
%   2.  Selection of the final B-value ('selected_B_final'): 
%       The optimal B is determined either when stability converges (i.e., the relative improvement falls below a tolerance) - 
%       or, if convergence is not reached, through an "elbow analysis" of the stability curve.
%   3.  Visualization (Convergence): 
%       Generates and saves the convergence plot showing the stability analysis using 'HERA.plot.rank_convergence'.
%   4.  Final bootstrap analysis: 
%       Using the optimal B-value, the function runs a final, comprehensive bootstrap analysis to generate the definitive rank distribution for each dataset.
%       This step is parallelized and uses 'config.system.target_memory' for memory-aware batching.
%   5.  Result Output (Console): 
%       Calculates the frequency of each rank for each dataset and prints a formatted table to the console.
%   6.  Result Output (csv): 
%       Saves the detailed rank distribution (counts and percentages) for each dataset to a CSV file.
%   7.  Visualization (Distribution): 
%       Generates and saves histograms displaying the final rank distribution for each dataset using 'HERA.plot.rank_histograms'.
%   8.  CI Calculation:
%       Computes the user-defined confidence intervals (e.g., 95%) for the final ranks based on the bootstrap distribution.
%
% Inputs:
%   all_data        - Cell array containing the original data matrices for each metric (1, 2, or 3 cells).
%   thresholds      - Struct with the pre-calculated effect size thresholds (thresholds.d_thresh, thresholds.rel_thresh).
%   config          - Global configuration struct. Uses:
%                       .bootstrap_ranks (settings)
%                       .ranking_mode (logic)
%                       .ci_level (confidence interval width)
%                       .system (optional performance limits & RAM settings)
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

% Extract system limits safely (pass [] if missing/default)
delta_mat_limit = [];
if isfield(config, 'system') && isfield(config.system, 'delta_mat_limit')
    delta_mat_limit = config.system.delta_mat_limit;
end

min_batch_size = 100;
if isfield(config, 'system') && isfield(config.system, 'min_batch_size')
    min_batch_size = config.system.min_batch_size;
end

%% 1. Initialization and Convergence Check for Optimal B
% Iteratively tests increasing B-values until the rank CI widths stabilize.
%
% Workflow:
%   a) Outer loop: Iterates over B-values (B_start:B_step:B_end).
%   b) Parallel stability check: For each B, runs n_trials parallel bootstrap analyses.
%   c) Convergence detection: Uses HERA.stats.check_convergence to detect plateau.
%   d) Fallback: If no convergence, uses HERA.stats.find_elbow_point.
%
% RNG Strategy:
%   - Each trial gets a unique substream: t_b (1..n_trials)
%   - This ensures reproducibility while utilizing all available CPU cores.
%
% Initialization of local variables from input parameters.
ts = config.timestamp;
num_datasets_b = size(all_data{1}, 2);
n_subj_b = num_probanden;
num_metrics = numel(all_data); % Get dynamic number of metrics

% Create a dedicated subfolder for the bootstrap_ranking plots.
subfolder_ranking = fullfile(graphics_dir, 'Ranking');

% Initialize output arrays for figure handles to prevent errors if no plots are generated.
h_figs_rank = gobjects(0);
h_fig_hist_rank = gobjects(0); 

% Check if a manual B-value was provided by the user.
cfg_rank = config.bootstrap_ranks; % Always load config for n_trials reference

if ~isempty(manual_B)
    % If so, use the provided value and skip the dynamic search.
    selected_B_final = manual_B;
    stability_data_rank = []; % No stability analysis performed
    fprintf(['\n' lang.ranking.manual_b_info '\n'], selected_B_final);
    % Skip directly to the final calculation in Section 3.
else
    % If no manual B is given, perform the dynamic search for the optimal B.
    
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
        ci_widths_b = zeros(round(cfg_rank.n_trials), num_datasets_b);
        
        % --- Parallel Worker Limit ---
        % Allows external control of worker count for nested parallelism scenarios.
        % When config.num_workers is set, limits the parfor to that many workers.
        % Default: Uses all available pool workers.
        pool = gcp('nocreate');
        current_pool_size = Inf;
        if ~isempty(pool)
            try
                current_pool_size = pool.NumWorkers;
            catch
                % Accessing pool properties (like NumWorkers) is not allowed on workers
                % We default to Inf so that config.num_workers (if set) is respected,
                % or parfor uses all available resources.
            end
        end

        if isfield(config, 'num_workers') && isnumeric(config.num_workers) && config.num_workers > 0
            parfor_limit = min(current_pool_size, config.num_workers);
        else
            parfor_limit = current_pool_size;
        end
        
        % --- Pre-Calculation of Safety Batching Parameters ---
        % Workflow (Memory Safety):
        %   a) Memory-aware batch sizing: Splits B into chunks if RAM is tight.
        %   b) Vectorized Bootstrap: Each batch computed efficiently.
        %   c) Aggregation: Prevents OOM errors for large N.
        %
        % RNG Strategy:
        %   - Preserves bit-perfect sequences via column-wise randi generation.

        % Estimate memory per iteration (ranking involves 3D arrays [Pairs x Metrics x B])
        mem_per_iter_bytes = (n_subj_b * 4) + (size(pair_idx_all, 1) * num_metrics * 2 * 8);
        
        % Use helper to determine batch size
        [BATCH_SIZE_PAR, num_batches_par] = HERA.run.get_batch_config(config, Br_b, mem_per_iter_bytes);
        
        % Perform n_trials to check the stability of the rank confidence intervals for the current B-value.
        parfor (t_b = 1:double(int32(round(cfg_rank.n_trials))), double(int32(round(parfor_limit))))
            % Each parallel worker gets its own reproducible substream of the random number generator.
            s_worker = s;
            s_worker.Substream = t_b;
            rank_tmp_b = zeros(num_datasets_b, Br_b);
            
            for b_loc = 1:num_batches_par
                 start_idx_loc = (b_loc - 1) * BATCH_SIZE_PAR + 1;
                 end_idx_loc = min(b_loc * BATCH_SIZE_PAR, Br_b);
                 current_n_loc = end_idx_loc - start_idx_loc + 1;
                 
                 % Inner loop: Vectorized Bootstrap + Ranking (Chunked)
                 % Generate all bootstrap samples at once [N x B]
                 boot_indices_block = randi(s_worker, n_subj_b, [n_subj_b, current_n_loc]);
            
                 % Pre-allocate 3D arrays for effect sizes [Pairs x Metrics x B]
                 num_pairs = size(pair_idx_all, 1);
                 d_vals_3d = zeros(num_pairs, num_metrics, current_n_loc);
                 rel_vals_3d = zeros(num_pairs, num_metrics, current_n_loc);
                 
                 % Check NaN status for this trial's data context (checking original data once is enough)
                  has_nans = false;
                  for col = 1:num_datasets_b
                     for mid = 1:num_metrics
                         if any(isnan(all_data{mid}(:, col)))
                             has_nans = true; break;
                         end
                     end
                     if has_nans, break; end
                  end
                  
                  if ~has_nans
                     % Fast path: Vectorized calculation
                     for p_idx = 1:num_pairs
                         i = pair_idx_all(p_idx, 1);
                         j = pair_idx_all(p_idx, 2);
                         for metric_idx = 1:num_metrics
                              data_i = all_data{metric_idx}(:, i);
                              data_j = all_data{metric_idx}(:, j);
                              x_mat = data_i(boot_indices_block);
                              y_mat = data_j(boot_indices_block);
                              
                              d_vals_3d(p_idx, metric_idx, :) = HERA.stats.cliffs_delta(x_mat, y_mat, delta_mat_limit);
                              rel_vals_3d(p_idx, metric_idx, :) = HERA.stats.relative_difference(x_mat, y_mat);
                         end
                     end
                 else
                     % Robust path: Loop-based calculation for NaNs
                     for k = 1:current_n_loc
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
                                     d_vals_3d(p_idx, metric_idx, k) = HERA.stats.cliffs_delta(x, y, delta_mat_limit);
                                     rel_vals_3d(p_idx, metric_idx, k) = HERA.stats.relative_difference(x, y);
                                 else
                                     d_vals_3d(p_idx, metric_idx, k) = NaN;
                                     rel_vals_3d(p_idx, metric_idx, k) = NaN;
                                 end
                             end
                          end
                     end
                 end
                 
                 % Calculate Ranking for each B
                 for bb_b = 1:current_n_loc
                      es_struct = struct('d_vals_all', d_vals_3d(:, :, bb_b), ...
                                         'rel_vals_all', rel_vals_3d(:, :, bb_b));
                                     
                     [~, bootstrap_rank] = HERA.calculate_ranking(...
                         all_data, es_struct, thresholds, config, dataset_names, pair_idx_all, boot_indices_block(:, bb_b));
                     
                     % Flatten index for result storage
                     global_b_idx = start_idx_loc + bb_b - 1;
                     rank_tmp_b(:, global_b_idx) = bootstrap_rank;
                 end
            end
            
            % After Br_b iterations, calculate the width of the 95% confidence interval of the ranks for this one trial.
            ci_bounds_b = quantile(rank_tmp_b, [0.025, 0.975], 2);
            ci_widths_b(t_b, :) = ci_bounds_b(:, 2) - ci_bounds_b(:, 1);
        end
        
        % Calculate the stability metric for the current B-value. 
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
             if isfield(cfg_rank, 'convergence_streak_needed') && ~isempty(cfg_rank.convergence_streak_needed) && cfg_rank.convergence_streak_needed > 0
                  if abs(results.improvement) < cfg_rank.convergence_tolerance
                     fprintf(['    ' lang.ranking.convergence_run_info '\n'], results.improvement * 100, results.streak, cfg_rank.convergence_streak_needed);
                  else
                     fprintf(['    ' lang.ranking.stability_change_info '\n'], results.improvement * 100);
                  end
             else
                  % Simple mode logs
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
        
%% 2. Selection of the final B-value based on the analysis outcome.
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
    % Ensure selected_B_final is a valid integer for parfor
    selected_B_final = round(selected_B_final);
    
%% 3. Create and save the convergence graphic.
    h_fig_rank = HERA.plot.rank_convergence(B_tested_vector_b, stability_vector_b_plotted, ...
                                            selected_B_final, converged, elbow_idx_rank, ...
                                            config, styles, lang);

    if isgraphics(h_fig_rank)
         % Save graphic.
        if ~exist(subfolder_ranking, 'dir'); mkdir(subfolder_ranking); end
        [~, fName, fExt] = fileparts(lang.files.convergence_rank_stability);
        filename = fullfile(subfolder_ranking, [fName, '_', ts, fExt]);
        exportgraphics(h_fig_rank, filename, 'Resolution', 300, 'BackgroundColor', styles.colors.background, 'Padding', 30);
        fprintf([lang.ranking.convergence_plot_saved '\n'], filename);

        h_figs_rank(end+1) = h_fig_rank;
    end
end

%% 4. Final Bootstrap Analysis for Rank Distributions
% Generates the final rank distribution using the optimally determined B value.
%
% Workflow:
%   a) Memory-aware batch sizing: Splits B into chunks that fit in RAM.
%   b) Parallel bootstrap over batches: Each batch gets its own RNG substream.
%   c) Fast/Robust path: Uses vectorized calculation if no NaNs, otherwise loop-based.
%   d) Aggregation: Flattens batch results into the final [num_datasets x B] matrix.
%
% RNG Strategy:
%   - Each batch gets a unique substream: OFFSET_BOOTSTRAP + b_idx
%   - Offset ensures no overlap with Phase 1 (stability check uses substreams 1..n_trials).

% Estimate total memory needed for all iterations.
bytes_per_double = 8;
bytes_per_int = 4;
mem_per_iter_bytes = (n_subj_b * bytes_per_int) + (size(pair_idx_all, 1) * num_metrics * 2 * bytes_per_double);

% Use helper to determine batch size
[BATCH_SIZE, num_batches] = HERA.run.get_batch_config(config, selected_B_final, mem_per_iter_bytes);

rank_batches = cell(1, num_batches);

% Determine base offset from config or fallback to 1000
if isfield(config, 'bootstrap_seed_offset') && ~isempty(config.bootstrap_seed_offset)
    base_offset = config.bootstrap_seed_offset;
else
    base_offset = 1000;
end

% Use substream offset to avoid overlap with stability analysis phase.
if isfield(cfg_rank, 'n_trials')
    OFFSET_BOOTSTRAP = cfg_rank.n_trials + base_offset;
else
    OFFSET_BOOTSTRAP = base_offset;
end

parfor b_idx = 1:num_batches
    s_worker = s;
    s_worker.Substream = OFFSET_BOOTSTRAP + b_idx; % Offset to avoid overlap with stability check
    
    % Determine batch range.
    start_idx = (b_idx - 1) * BATCH_SIZE + 1;
    end_idx = min(b_idx * BATCH_SIZE, selected_B_final);
    current_batch_size = end_idx - start_idx + 1;
    
    % Generate bootstrap indices [N x BatchSize].
    boot_indices_block = randi(s_worker, n_subj_b, [n_subj_b, current_batch_size]);
    
    % Calculate effect sizes for this batch.
    num_pairs = size(pair_idx_all, 1);
    
    % Check for NaNs to decide on calculation strategy.
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
        % Fast path: Vectorized calculation.
        for p_idx = 1:num_pairs
            i = pair_idx_all(p_idx, 1);
            j = pair_idx_all(p_idx, 2);
            for metric_idx = 1:num_metrics
                 data_i = all_data{metric_idx}(:, i);
                 data_j = all_data{metric_idx}(:, j);
                 x_mat = data_i(boot_indices_block);
                 y_mat = data_j(boot_indices_block);
                 
                 d_vals_3d(p_idx, metric_idx, :) = HERA.stats.cliffs_delta(x_mat, y_mat, delta_mat_limit);
                 rel_vals_3d(p_idx, metric_idx, :) = HERA.stats.relative_difference(x_mat, y_mat);
            end
        end
    else
        % Robust path: Loop-based calculation for data with NaNs.
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
                        d_vals_3d(p_idx, metric_idx, k) = HERA.stats.cliffs_delta(x, y, delta_mat_limit);
                        rel_vals_3d(p_idx, metric_idx, k) = HERA.stats.relative_difference(x, y);
                    else
                        d_vals_3d(p_idx, metric_idx, k) = NaN;
                        rel_vals_3d(p_idx, metric_idx, k) = NaN;
                    end
                end
             end
        end
    end
    
    % Calculate ranking for each bootstrap iteration.
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

%% 4 and 5. Print & Save Bootstrap Rank Distribution
HERA.output.save_ranking_table(final_bootstrap_ranks, final_rank, dataset_names, selected_B_final, lang, csv_dir, ts);

%% 6. Create and save the histogram distribution of the final ranks
% Only create this plot if reports are enabled
if isfield(config, 'create_reports') && config.create_reports
    h_fig_hist_rank = HERA.plot.rank_histograms(final_bootstrap_ranks, dataset_names, ...
                                                final_rank, selected_B_final, styles, lang);

    if isgraphics(h_fig_hist_rank)
        % Save the complete histogram figure to a file.
        if ~exist(subfolder_ranking, 'dir'); mkdir(subfolder_ranking); end
        [~, fName, fExt] = fileparts(lang.files.dist_bootstrap_ranks);
        filename = fullfile(subfolder_ranking, [fName, '_', ts, fExt]);
        exportgraphics(h_fig_hist_rank, filename, 'Resolution', 300, 'BackgroundColor', styles.colors.background, 'Padding', 30);
        fprintf([lang.ranking.histogram_plot_saved '\n'], filename);
    end
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
