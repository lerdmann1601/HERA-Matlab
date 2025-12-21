function [d_thresh, rel_thresh, rel_thresh_b, min_rel_thresh, d_vals_all, rel_vals_all, pair_idx_all, selected_B, stability_data_thr, ...
          h_fig_thr_global, h_fig_thr_detailed, h_fig_hist_thr, h_fig_hist_raw] = ...
          calculate_thresholds(all_data, num_probanden, config, graphics_dir, manual_B, s, styles, lang)
% CALCULATE_THRESHOLDS - Determines statistical thresholds for effect sizes.
% 
% Syntax:
%   [d_thresh, rel_thresh, rel_thresh_b, min_rel_thresh, d_vals_all, rel_vals_all, pair_idx_all, ...
%    selected_B, stability_data_thr, h_fig_thr_global, h_fig_thr_detailed, h_fig_hist_thr, h_fig_hist_raw] = ...
%    calculate_thresholds(all_data, num_probanden, config, graphics_dir, manual_B, s, styles, lang)
%
% Description:
%   This function analyzes data to find robust thresholds for Cliff's Delta and the relative difference.
%
% Workflow:
%   1. Initialization & Initial Statistics: 
%      Calculates effect sizes for all pairwise comparisons (for 1, 2, or 3 metrics).
%   2. SEM-based Minimum Threshold: 
%      Determines a dynamic lower limit for the relative difference based on the SEM (for 1, 2, or 3 metrics).
%   3. Dynamic Determination of Stable Bootstrap Thresholds (for B): 
%      Iteratively checks stability using parallelized trials (`parfor`) with memory-aware 
%      batch sizing. Checks convergence via `check_convergence` or `find_elbow_point`.
%   4. Final Threshold Calculation: 
%      Computes the final thresholds with the optimal B, where the relative threshold is capped by the SEM value.
%   5. Visualization:
%      Delegates the creation of convergence and distribution plots to `HERA.plot.threshold_convergence`
%      and `HERA.plot.threshold_distributions`.
%
% Inputs:
%   all_data      - Cell array of data matrices for each metric (1, 2, or 3).
%   num_probanden - Number of subjects.
%   config        - Struct with all configuration parameters (including 'config.bootstrap_thresholds' and optional 'config.system' limits).
%   manual_B      - (Optional) Manually set number of bootstrap repetitions (B). 
%   s             - Random number stream for reproducibility.
%   styles        - Struct with plot styling information.
%   lang          - Language package loaded from a JSON file.
%
% Outputs:
%   d_thresh           - Vector of Cliff's Delta thresholds per metric.
%   rel_thresh         - Vector of relative difference thresholds per metric.
%   min_rel_thresh     - Vector of SEM-based minimum thresholds per metric.
%   d_vals_all         - Matrix of all calculated Cliff's Delta values.
%   rel_vals_all       - Matrix of all calculated relative differences.
%   pair_idx_all       - Indices of all pairwise comparisons.
%   selected_B         - Optimally determined number of bootstrap samples.
%   stability_data_thr - Struct with convergence curve data for JSON export.
%   h_fig_thr_global   - Handle of the global convergence plot.
%   h_fig_thr_detailed - Handle of the detailed convergence plot.
%   h_fig_hist_thr     - Handle of the bootstrap median distribution plot.
%   h_fig_hist_raw     - Handle of the raw effect size distribution plot.
%
% Author: Lukas von Erdmannsdorff

arguments
    all_data (1,:) cell
    num_probanden (1,1) double
    config (1,1) struct
    graphics_dir (1,1) string
    manual_B
    s
    styles (1,1) struct
    lang (1,1) struct
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

%% 1. Initialization and Calculation of Initial Statistics
% Create a dedicated subfolder for the threshold analysis plots.
subfolder_thresholds = fullfile(graphics_dir, 'Threshold_Analysis');
if ~exist(subfolder_thresholds, 'dir')
    mkdir(subfolder_thresholds);
end
% Extracts parameters from the configuration structure.
metric_names = config.metric_names;
num_metrics = numel(all_data); % Get dynamic metric count
num_datasets = size(all_data{1}, 2);
ci_level = config.ci_level;
alpha_level = 1 - ci_level;
ts = config.timestamp; 

% Calculates the means for each dataset and metric for later use.
all_means = cell(1, num_metrics);
for m = 1:num_metrics
    all_means{m} = mean(all_data{m}, 1);
end

% Generates all possible unique pair combinations of the datasets.
num_pairs = nchoosek(num_datasets, 2);
pair_idx_all = nchoosek(1:num_datasets, 2);

% Preallocation of matrices for effect sizes (Cliff's Delta and Relative Difference).
d_vals_all = zeros(num_pairs, num_metrics); 
rel_vals_all = zeros(num_pairs, num_metrics); 
% Loop over all pairs to calculate the point estimates of the effect sizes.
idx_pair = 0;
for i = 1:(num_datasets - 1)
    for j = (i + 1):num_datasets
        idx_pair = idx_pair + 1;
        % Loop over all metrics (Metric 1, 2 and 3).
        for metric_idx = 1:num_metrics % Dynamic loop
            data_metric = all_data{metric_idx};
            
            % Find valid rows for this pair only (pairwise exclusion).
            valid_rows = ~isnan(data_metric(:, i)) & ~isnan(data_metric(:, j));
            x = data_metric(valid_rows, i);
            y = data_metric(valid_rows, j);
            n_valid = size(x, 1);
            
            if n_valid > 0
                % Calculation of Cliff's Delta with the valid sample size.
                d_vals_all(idx_pair, metric_idx) = HERA.stats.cliffs_delta(x, y, delta_mat_limit);
                
                % Calculation of the relative mean difference.
                rel_vals_all(idx_pair, metric_idx) = HERA.stats.relative_difference(x, y);
            else
                % If no valid pairs exist, store NaN.
                d_vals_all(idx_pair, metric_idx) = NaN;
                rel_vals_all(idx_pair, metric_idx) = NaN;
            end
        end
    end
end

%% 2. SEM-based minimum threshold as a lower bound
% Determines the t-distribution value for the given confidence level.
degrees_of_freedom = num_probanden - 1;
multiplier = tinv(1 - alpha_level / 2, degrees_of_freedom); % t-value for the confidence interval.

% Calculates a dynamic, minimum relative threshold based on the standard error of the mean (SEM).
min_rel_dynamic = zeros(1, num_metrics); 
for metric_idx = 1:num_metrics 
    data_matrix = all_data{metric_idx};
    
    % Calculates the SEM for each dataset, ignoring NaN values.
    % Important: The sample size (n) must be determined individually for each column.
    n_per_dataset = sum(~isnan(data_matrix), 1); 
    sems = std(data_matrix, 0, 1, 'omitnan') ./ sqrt(n_per_dataset);
    
    % Uses the median of the SEMs for robustness and ignores NaN values.
    median_sem = median(sems, 'omitnan'); 
    
    % Calculates the minimum relative difference based on the median SEM.
    % Also ignores NaN values for the grand mean.
    grand_mean = mean(data_matrix(:), 'omitnan');

    % Protects against division by zero if all data are NaN.
    if grand_mean == 0 || isnan(grand_mean)
        min_rel_dynamic(metric_idx) = 0;
    else
        min_rel_dynamic(metric_idx) = (multiplier * median_sem) / abs(grand_mean);
    end
end
min_rel_thresh = min_rel_dynamic; % Assigns the result to the output variable.

%% 3. Dynamic determination of stable bootstrap thresholds (for B)
% Iteratively tests increasing B-values until the threshold estimates stabilize.
%
% Workflow:
%   a) Iterate sequentially over B-values (B_start:B_step:B_end).
%   b) Iterate sequentially over each effect type / metric.
%   c) Parallelize stability trials: Execute n_trials in parallel using `parfor`.
%      This ensures maximum CPU utilization even with few metrics.
%   d) Convergence detection: Uses HERA.stats.check_convergence.
%   e) Fallback: Uses HERA.stats.find_elbow_point.
%
% RNG Strategy:
%   - Deterministic substream assignment: (m - 1) * 1000 + trial_idx
%   - Ensures reproducibility at the trial level independent of parallel execution order.
%
% Checks if a manual B value was passed.
if ~isempty(manual_B)
    selected_B = manual_B;
    % Assign empty handles for the plot outputs that are skipped
    h_fig_thr_global = gobjects(0); 
    h_fig_thr_detailed = gobjects(0);
    stability_data_thr = []; % No stability analysis performed
    fprintf(['\n' lang.thresholds.manual_b_info '\n'], selected_B);    
    % Jumps directly to the final calculation in Section 5.
else
    % Extracts the configuration for this analysis step.
    cfg_thr = config.bootstrap_thresholds;
    % Checks if the robust convergence check (with smoothing) should be used.
    use_robust_convergence_thr = isfield(cfg_thr, 'smoothing_window') && ~isempty(cfg_thr.smoothing_window) ...
                               && isfield(cfg_thr, 'convergence_streak_needed') && ~isempty(cfg_thr.convergence_streak_needed);
    
    fprintf(['\n' lang.thresholds.searching_optimal_b '\n']);
    if use_robust_convergence_thr    
        fprintf([lang.thresholds.primary_criterion '\n'])
        fprintf([lang.thresholds.robust_convergence_info '\n'], ...
            cfg_thr.smoothing_window, ...
            cfg_thr.convergence_streak_needed, ...
            cfg_thr.convergence_tolerance * 100, ...
            cfg_thr.min_steps_for_convergence_check, ...
            cfg_thr.B_end);
    else 
        fprintf([lang.thresholds.simple_convergence_info '\n'], ...
        cfg_thr.convergence_tolerance * 100, cfg_thr.min_steps_for_convergence_check, cfg_thr.B_end);
    end
    fprintf([lang.thresholds.secondary_criterion '\n']);
    
    % Initializes vectors and matrices for the stability analysis.
    B_vector = cfg_thr.B_start:cfg_thr.B_step:cfg_thr.B_end;
    % Dynamic matrix size based on num_metrics
    stability_matrix = NaN(2, num_metrics, numel(B_vector)); % Dimensions: [d/r, metric, B_value]
    final_i = 0;
    converged = false;
    % Initialization for convergence check.
    convergence_streak_counter_thr = 0;
    overall_stability_thr = NaN(1, numel(B_vector));
    
    % Main loop: Iterates over different numbers of bootstrap samples (B).
    for i = 1:numel(B_vector)
        B_current = B_vector(i);
        fprintf([' -> ' lang.ranking.checking_stability '\n'], B_current, cfg_thr.n_trials);
        % Dynamic vector size based on num_metrics
        temp_stability_vector = zeros(1, num_metrics * 2); % Vector for [d1..dN, r1..rN]
        
        % Total tasks: all effect types × n_trials.
        n_effect_types = num_metrics * 2;
        total_tasks = n_effect_types * cfg_thr.n_trials;
        
        % Extract values once per effect type (reused for all trials).
        vals_cell = cell(1, n_effect_types);
        n_vals_vec = zeros(1, n_effect_types);
        for m = 1:n_effect_types
            is_delta = m <= num_metrics;
            actual_metric_idx = mod(m - 1, num_metrics) + 1;
            if is_delta
                all_vals = abs(d_vals_all(:, actual_metric_idx));
            else
                all_vals = rel_vals_all(:, actual_metric_idx);
            end
            vals_cell{m} = all_vals(~isnan(all_vals));
            n_vals_vec(m) = numel(vals_cell{m});
        end
        
        % --- Parallel Worker Limit ---
        % Limits the number of workers if config.num_workers is set, otherwise utilizes the full pool.
        pool = gcp('nocreate');
        current_pool_size = Inf;
        if ~isempty(pool), try, current_pool_size = pool.NumWorkers; catch, end; end

        if isfield(config, 'num_workers') && isnumeric(config.num_workers) && config.num_workers > 0
            parfor_limit = min(current_pool_size, config.num_workers);
        else
            parfor_limit = current_pool_size;
        end
        
        % Data Preparation for all effect types (minimizing broadcast overhead)
        % vals_cell and n_vals_vec are already prepared above.
        
        % Loop over effect types SEQUENTIALLY
        for m = 1:n_effect_types
            
            vals = vals_cell{m};
            n_vals = n_vals_vec(m);
            
            % Skip calculation if no valid data
            if isempty(vals) || n_vals < 2
                temp_stability_vector(m) = 0; % Or handle as needed
                continue;
            end
            
            % --- Pre-Calculation of Safety Batching Parameters ---
            % Workflow (Memory Safety):
            %   a) Memory-aware batch sizing: Splits B into chunks if RAM is tight.
            %   b) Vectorized Bootstrap: Each batch computed efficiently.
            %   c) Aggregation: Prevents OOM errors for large N.
            %
            % RNG Strategy:
            %   - Preserves bit-perfect sequences via column-wise randi generation.

            % Dynamic batch sizing based on memory configuration.
            if isfield(config, 'system') && isfield(config.system, 'target_memory')
                 TARGET_MEMORY_LOC = config.system.target_memory;
                 if numel(TARGET_MEMORY_LOC) > 1, TARGET_MEMORY_LOC = TARGET_MEMORY_LOC(1); end
            else
                 TARGET_MEMORY_LOC = 200;
            end
            
            effective_memory_loc = TARGET_MEMORY_LOC;
            bytes_per_sample = n_vals * 8; % Double precision
            total_memory_needed = (double(B_current) * double(bytes_per_sample)) / (1024^2);
            
            if total_memory_needed <= double(effective_memory_loc)
                 BATCH_SIZE_PAR = double(B_current);
            else
                 BATCH_SIZE_PAR = max(min_batch_size, min(floor((double(effective_memory_loc) * 1024^2) / double(bytes_per_sample)), 20000));
            end
            num_batches_par = double(ceil(double(B_current) ./ BATCH_SIZE_PAR));

            % --- Inner Parallel Loop (Trials) ---
            % Parallelizes the bootstrap trials to ensure maximum core utilization.
            n_trials_int = double(int32(round(cfg_thr.n_trials)));
            parfor_limit = double(int32(round(parfor_limit)));
            thr_trials = zeros(n_trials_int, 1);
            
            parfor (t = 1:n_trials_int, parfor_limit)
                % RNG: Deterministic substream based on EffectType and Trial index.
                s_worker = s;
                s_worker.Substream = (m - 1) * 1000 + t;
                
                % Initialize accumulation vector
                bootstat = zeros(1, B_current);
                
                for b_loc = 1:num_batches_par
                     start_idx_loc = (b_loc - 1) * BATCH_SIZE_PAR + 1;
                     end_idx_loc = min(b_loc * BATCH_SIZE_PAR, B_current);
                     current_n_loc = end_idx_loc - start_idx_loc + 1;
                
                     
                     % Generate indices for this batch [N x Batch]
                     boot_indices = randi(s_worker, n_vals, [n_vals, current_n_loc]);
                     
                     % Calculate stats
                     bootstat(start_idx_loc:end_idx_loc) = median(vals(boot_indices), 1);
                end
                
                % Estimate threshold (Lower CI bound)
                ci_tmp = quantile(bootstat, [alpha_level / 2, 1 - alpha_level / 2]);
                thr_trials(t) = ci_tmp(1);
            end
            
            % Calculates stability as a relative measure of dispersion (IQR / Median).
            med_thr = median(thr_trials);
            iqr_thr = iqr(thr_trials);
            
            if med_thr == 0
                if iqr_thr == 0
                    temp_stability_vector(m) = 0;
                else
                    temp_stability_vector(m) = Inf;
                end
            else
                temp_stability_vector(m) = iqr_thr / abs(med_thr);
            end
        end
        
        % Stores the stability values for each metric and effect size in the main matrix.
        % Dynamic slicing based on num_metrics
        stability_matrix(1, :, i) = temp_stability_vector(1:num_metrics); % Cliff's Delta stabilities
        stability_matrix(2, :, i) = temp_stability_vector(num_metrics+1 : end); % Relative Difference stabilities
        final_i = i; % Keep track of the last executed iteration index.
        
        % Calculate a single, overall stability value for the current B by averaging across all metrics.
        % This single value is used to check for convergence.
        overall_stability_thr(i) = mean(stability_matrix(:, :, i), 'all', 'omitnan');
       
        % Convergence check: Determines if the stability has plateaued, indicating a sufficient B value.
        [converged, conv_stats] = HERA.stats.check_convergence(overall_stability_thr(1:i), cfg_thr);
        
        if use_robust_convergence_thr
             % Robust method logging
             if i >= cfg_thr.min_steps_for_convergence_check + cfg_thr.smoothing_window
                 if abs(conv_stats.improvement) < cfg_thr.convergence_tolerance
                     fprintf(['    ' lang.thresholds.convergence_run_info '\n'], conv_stats.improvement * 100, conv_stats.streak, cfg_thr.convergence_streak_needed);
                 else
                     fprintf(['    ' lang.thresholds.stability_change_info '\n'], conv_stats.improvement * 100);
                 end
                 
                 if converged
                     fprintf([lang.thresholds.convergence_reached '\n'], conv_stats.improvement * 100);
                     fprintf([lang.thresholds.stable_runs_info '\n'], cfg_thr.convergence_streak_needed);
                 end
             end
        else
             % Simple method logging
             if i >= cfg_thr.min_steps_for_convergence_check
                 if ~isnan(conv_stats.improvement)
                     fprintf(['    ' lang.thresholds.stability_change_info '\n'], conv_stats.improvement * 100);
                 end
                 if converged
                     fprintf([lang.thresholds.convergence_reached '\n'], conv_stats.improvement * 100);
                 end
             end
        end

        if converged
            break; % Exit the main for-loop over B values if convergence has been achieved.
        end
    end % End of for-loop for stability check.
    
    % Determination of the final B value based on the analysis outcome.
    B_tested_vector = B_vector(1:final_i);
    stability_matrix = stability_matrix(:, :, 1:final_i); 
    overall_stability_thr_plotted = overall_stability_thr(1:final_i);
    if converged
        % If the process converged, the last tested B value is selected as the optimal one.
        selected_B = B_tested_vector(end);
        fprintf([lang.thresholds.convergence_result '\n'], selected_B);
        elbow_indices = [];
    else
        % Otherwise, the "elbow" of the curve is determined as the optimal point.
        fprintf([lang.thresholds.elbow_analysis_info '\n']);
        
        % Dynamically build stability_vectors cell array
        stability_vectors = cell(1, num_metrics * 2);
        for k_vec = 1:(num_metrics * 2)
            is_delta = k_vec <= num_metrics;
            actual_metric_idx = mod(k_vec - 1, num_metrics) + 1;
            if is_delta
                stability_vectors{k_vec} = squeeze(stability_matrix(1, actual_metric_idx, :));
            else
                stability_vectors{k_vec} = squeeze(stability_matrix(2, actual_metric_idx, :));
            end
        end
        
        % Use the helper function to find the elbow points
        [selected_B, elbow_indices] = HERA.stats.find_elbow_point(B_tested_vector, stability_vectors);
        fprintf([lang.thresholds.elbow_result '\n'], selected_B);
    end
    % Ensure selected_B is a valid integer for parfor
    selected_B = round(selected_B);
    % Store stability data for JSON export and plotting
    stability_data_thr = struct();
    stability_data_thr.B_vector = B_tested_vector;
    stability_data_thr.global_stability = overall_stability_thr_plotted;
    stability_data_thr.detailed_stability = stability_matrix;
    stability_data_thr.elbow_indices = elbow_indices;
    stability_data_thr.converged = converged;

%% 4. Creates and saves plots to visualize the convergence analysis
% Delegate to the plotting controller
[h_fig_thr_global, h_fig_thr_detailed] = HERA.plot.threshold_convergence(...
    B_tested_vector, stability_data_thr, selected_B, config, styles, lang, graphics_dir);

end

%% 5. Final calculation of thresholds with the optimal B
% Calculates the final thresholds for Cliff's Delta and relative difference.
%
% Workflow:
%   a) Memory-aware batch sizing: Splits B into chunks that fit in RAM.
%   b) Parallel bootstrap over batches: Each batch gets its own RNG substream.
%   c) Aggregation: Flattens batch results and calculates the lower CI bound as threshold.
%
% RNG Strategy:
%   - Unique substream per (metric, effect_type, batch): offset_d + b or offset_rel + b
%   - Offsets are calculated to prevent collisions between metrics and effect types.
%
% Initialization of the output variables.
d_thresh = zeros(1, num_metrics); 
rel_thresh = zeros(1, num_metrics); 
all_bootstat_d = cell(1, num_metrics); 
all_bootstat_rel = cell(1, num_metrics);

% Dynamic batch sizing based on memory configuration.
if isfield(config, 'system') && isfield(config.system, 'target_memory')
     TARGET_MEMORY = config.system.target_memory;
     if numel(TARGET_MEMORY) > 1
         TARGET_MEMORY = TARGET_MEMORY(1);
     end
else
     TARGET_MEMORY = 200;
end

% Get worker count for effective memory
if isfield(config, 'num_workers') && isnumeric(config.num_workers)
    num_workers = config.num_workers;
    if numel(num_workers) > 1
        num_workers = num_workers(1);
    end
else
    num_workers = feature('numcores');
end
effective_memory = TARGET_MEMORY / max(1, num_workers);

for metric_idx = 1:num_metrics
    % --- Cliff's Delta ---
    d_vals_metric_raw = d_vals_all(:, metric_idx);
    d_vals_metric_abs = abs(d_vals_metric_raw(~isnan(d_vals_metric_raw)));
    
    if ~isempty(d_vals_metric_abs)
        n_data_d = numel(d_vals_metric_abs);
        
        % Calculate batch size based on memory.
        bytes_per_sample = n_data_d * 8;
        
        % Smart batching for Final Phase with defensive scalar enforcement
        total_memory_needed = (double(selected_B) * double(bytes_per_sample)) / (1024^2);
        
        if total_memory_needed <= double(effective_memory)
            BATCH_SIZE = double(selected_B);
        else
            BATCH_SIZE = max(100, min(floor((double(effective_memory) * 1024^2) / double(bytes_per_sample)), 20000));
        end
        
        % Ensure scalar double for loop limit
        num_batches = double(ceil(double(selected_B) ./ BATCH_SIZE));
        if numel(num_batches) > 1
             num_batches = num_batches(1);
        end
        % EXTRA SAFETY for parfor
        if isempty(num_batches) || any(isnan(num_batches)) || any(isinf(num_batches))
            num_batches = 1;
        end
        num_batches = double(int32(round(num_batches)));
        
        % Substream offset for this metric.
        offset_d = 1000 + (metric_idx - 1) * 2 * (num_batches + 10);
        
        % Parallel bootstrap computation.
        results_cell_d = cell(1, num_batches);
        
        parfor b = 1:num_batches
            s_par = s;
            s_par.Substream = offset_d + b;
            
            start_idx = (b - 1) * BATCH_SIZE + 1;
            end_idx = min(b * BATCH_SIZE, selected_B);
            current_n = end_idx - start_idx + 1;
            
            % Generate indices.
            indices = randi(s_par, n_data_d, [n_data_d, current_n]);
            
            % Compute bootstrap statistic.
            results_cell_d{b} = median(d_vals_metric_abs(indices), 1);
        end
        
        bootstat_d = [results_cell_d{:}];
        d_thresh(metric_idx) = quantile(bootstat_d, alpha_level / 2);
        all_bootstat_d{metric_idx} = bootstat_d; 
    else
        d_thresh(metric_idx) = 0;
        all_bootstat_d{metric_idx} = [];
    end

    % --- Relative Difference ---
    rel_vals_metric_raw = rel_vals_all(:, metric_idx);
    rel_vals_metric = rel_vals_metric_raw(~isnan(rel_vals_metric_raw));

    if ~isempty(rel_vals_metric)
        n_data_rel = numel(rel_vals_metric);
        
        % Calculate batch size based on memory.
        bytes_per_sample = n_data_rel * 8;
        
        % Smart batching
        total_memory_needed = (double(selected_B) * double(bytes_per_sample)) / (1024^2);
        
        if total_memory_needed <= double(effective_memory)
            BATCH_SIZE = double(selected_B);
        else
            BATCH_SIZE = max(100, min(floor((double(effective_memory) * 1024^2) / double(bytes_per_sample)), 20000));
        end
        num_batches = double(ceil(double(selected_B) ./ BATCH_SIZE));
        if numel(num_batches) > 1
             num_batches = num_batches(1);
        end
        % EXTRA SAFETY for parfor
        if isempty(num_batches) || any(isnan(num_batches)) || any(isinf(num_batches))
            num_batches = 1;
        end
        num_batches = double(int32(round(num_batches)));
        
        % Substream offset (shifted from Delta).
        offset_rel = offset_d + (num_batches + 10);
        
        % Parallel bootstrap computation.
        results_cell_rel = cell(1, num_batches);
        
        parfor b = 1:num_batches
            s_par = s;
            s_par.Substream = offset_rel + b;
            
            start_idx = (b - 1) * BATCH_SIZE + 1;
            end_idx = min(b * BATCH_SIZE, selected_B);
            current_n = end_idx - start_idx + 1;
            
            % Generate indices.
            indices = randi(s_par, n_data_rel, [n_data_rel, current_n]);
            
            % Compute bootstrap statistic.
            results_cell_rel{b} = median(rel_vals_metric(indices), 1);
        end
        
        bootstat_rel = [results_cell_rel{:}];
        rel_thresh(metric_idx) = quantile(bootstat_rel, alpha_level / 2);
        all_bootstat_rel{metric_idx} = bootstat_rel;
    else
        rel_thresh(metric_idx) = 0;
        all_bootstat_rel{metric_idx} = [];
    end
end
% Stores the pure bootstrap value before it is possibly corrected.
rel_thresh_b = rel_thresh;
% Ensures that the relative threshold does not fall below the SEM-based limit.
rel_thresh = max(rel_thresh, min_rel_thresh);

%% 6. Creates and saves the histogram distribution of the bootstrap thresholds and raw effects
% Delegate to the plotting controller
[h_fig_hist_thr, h_fig_hist_raw] = HERA.plot.threshold_distributions(...
    all_bootstat_d, all_bootstat_rel, d_vals_all, rel_vals_all, d_thresh, rel_thresh, rel_thresh_b, ...
    selected_B, metric_names, styles, lang, graphics_dir, config);

end
