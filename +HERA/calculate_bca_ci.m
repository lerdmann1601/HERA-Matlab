function [B_ci, ci_d_all, ci_r_all, z0_d_all, a_d_all, z0_r_all, a_r_all, stability_data_ci, h_fig_bca_global, h_fig_bca_detailed, ...
          h_fig_hist_z0, h_fig_hist_a, h_fig_hist_widths] = ...
          calculate_bca_ci(all_data, d_vals_all, rel_vals_all, pair_idx_all, num_probanden, ...
          config, metric_names, graphics_dir, csv_dir, manual_B, s, styles, lang, base_name)
% CALCULATE_BCA_CI - Dynamically finds optimal B and calculates BCa CIs for effect sizes.
%
% Syntax:
%   [B_ci, ci_d_all, ci_r_all, z0_d_all, a_d_all, z0_r_all, a_r_all, stability_data_ci, ...
%    h_fig_bca_global, h_fig_bca_detailed, h_fig_hist_z0, h_fig_hist_a, h_fig_hist_widths] = ...
%   calculate_bca_ci(all_data, d_vals_all, rel_vals_all, pair_idx_all, num_probanden, ... 
%   config, metric_names, graphics_dir, csv_dir, manual_B, s, styles, lang, base_name)
%
% Description:
%   This function acts as the main controller for the BCa (Bias-Corrected and Accelerated) 
%   Bootstrap confidence interval calculation. It dynamically determines the optimal 
%   number of bootstrap samples (B) by analyzing the stability of the CI widths.
%
% Workflow:
%   1. Dynamic determination of the bootstrap count (B): 
%      Iterates over B-values in a parfor loop. Convergence is strictly checked using
%      `HERA.stats.check_convergence`. If no convergence is found, the optimal B is 
%      determined via elbow method using `HERA.stats.find_elbow_point`.
%   2. Creation of the convergence plot: 
%      Calls `HERA.plot.bca_convergence` to generate global and detailed stability plots.
%   3. Final calculation of BCa confidence intervals: 
%      Calculates the final CIs for Cliff's Delta and the relative difference using the 
%      determined optimal B-value. This step is parallelized and uses 'config.system.target_memory'
%      for memory-aware batching.
%   4. Calculation and output of correction factors: 
%      Summarizes the bias correction (z0) and skewness correction (a) and exports 
%      them to a CSV file.
%   5. Generation of distribution plots:
%      Calls `HERA.plot.bca_distributions` to create histograms for z0, a, and CI widths.
%
% Inputs:
%   all_data      - Cell array with the data matrices for each metric (1, 2, or 3 cells).
%   d_vals_all    - Matrix of the original Cliff's Delta effect sizes.
%   rel_vals_all  - Matrix of the original relative differences.
%   pair_idx_all  - Matrix of the indices for all pairwise comparisons.
%   num_probanden - Number of subjects.
%   config        - Struct containing configuration parameters (esp. config.bootstrap_ci 
%                   and config.system.target_memory for RAM optimization).
%   metric_names  - Cell array with the names of the metrics (1, 2, or 3 names).
%   graphics_dir  - Path to the output folder for graphics.
%   csv_dir       - Path to the output folder for CSV files.
%   manual_B      - (Optional) Manually set number of bootstrap repetitions (B).
%   s             - Random number stream for reproducibility.
%   styles        - Struct with plot styling information.
%   lang          - Language pack struct.
%   base_name     - Base name for output files (e.g., from the log timestamp).
%
% Outputs:
%   B_ci               - Optimally determined number of bootstrap samples.
%   ci_d_all           - 3D array of the BCa CI for Cliff's Delta.
%   ci_r_all           - 3D array of the BCa CI for the relative difference.
%   z0_d_all, a_d_all  - Matrices of correction factors (z0, a) for Cliff's Delta.
%   z0_r_all, a_r_all  - Matrices of correction factors (z0, a) for Rel. Diff.
%   stability_data_ci  - Struct with convergence curve data for JSON export.
%   h_fig_bca_global   - Handle of the global convergence graphic.
%   h_fig_bca_detailed - Handle of the detailed convergence graphic.
%   h_fig_hist_z0      - Handle of the z0 distribution graphic.
%   h_fig_hist_a       - Handle of the a distribution graphic.
%   h_fig_hist_widths  - Handle of the CI width distribution graphic.
%
% Author: Lukas von Erdmannsdorff

arguments
    all_data (1,:) cell
    d_vals_all
    rel_vals_all
    pair_idx_all
    num_probanden (1,1) double
    config (1,1) struct
    metric_names (1,:) cell
    graphics_dir (1,1) string
    csv_dir (1,1) string
    manual_B
    s
    styles (1,1) struct
    lang (1,1) struct
    base_name (1,1) string
end

%% 1. Dynamic determination of the optimal bootstrap count (B)
% Iteratively tests increasing B-values until the CI widths stabilize.
%
% Architecture:
%   a) Outer loop: Iterates over B-values (B_start:B_step:B_end).
%   b) Parallel stability check: For each B, parallelizes over (num_metrics * 2) effect types.
%   c) Inner loop: Each effect type runs n_trials to assess CI width variability.
%   d) Convergence detection: Uses HERA.stats.check_convergence to detect plateau.
%   e) Fallback: If no convergence, uses HERA.stats.find_elbow_point.
%
% RNG Strategy:
%   - Each effect type (metric_idx) gets a unique substream: 1..(num_metrics * 2)
%   - Trade-off: With 1 metric (2 effect types), only 2 cores are utilized.
%     This is accepted to maintain RNG consistency with the validated convergence behavior.
%
% Initialization of parameters from the configuration structure.
alpha_level = 1 - config.ci_level;
num_pairs = size(pair_idx_all, 1);
num_metrics = numel(metric_names); % Get dynamic number of metrics
ts = config.timestamp;

% Create a dedicated subfolder for the CI analysis plots.
subfolder_bca_CI = fullfile(graphics_dir, 'CI_Histograms');
if ~exist(subfolder_bca_CI, 'dir')
    mkdir(subfolder_bca_CI);
end

% Preparation of variables for the parfor loops
p_pair_idx_all = pair_idx_all; 
p_all_data = all_data; 
p_d_vals_all = d_vals_all; 
p_rel_vals_all = rel_vals_all;

if ~isempty(manual_B)
    B_ci = manual_B;
    % Assign empty handles for the plot outputs that are skipped
    h_fig_bca_global = gobjects(0); 
    h_fig_bca_detailed = gobjects(0);
    stability_data_ci = []; % No stability analysis performed
    fprintf(['\n' lang.bca.manual_b_info '\n'], B_ci);
    % Jumps directly to the final calculation in Section 3.
else
    % Extracts the configuration for this analysis step.
    cfg_ci = config.bootstrap_ci;
    % Checks if the robust convergence check (with smoothing) should be used.
    use_robust_convergence_ci = isfield(cfg_ci, 'smoothing_window') && ~isempty(cfg_ci.smoothing_window) ...
                               && isfield(cfg_ci, 'convergence_streak_needed') && ~isempty(cfg_ci.convergence_streak_needed);
    
    % Console output to inform the user about the process.
    fprintf(['\n' lang.bca.searching_optimal_b '\n']);
    if use_robust_convergence_ci
        fprintf([lang.bca.primary_criterion '\n'])
        fprintf([lang.bca.robust_convergence_info '\n'], ...
            cfg_ci.smoothing_window, ...
            cfg_ci.convergence_streak_needed, ...
            cfg_ci.convergence_tolerance * 100, ...
            cfg_ci.min_steps_for_convergence_check, ...
            cfg_ci.B_end);
    else
        fprintf([lang.bca.simple_convergence_info '\n'], ...
            cfg_ci.convergence_tolerance * 100, ...
            cfg_ci.min_steps_for_convergence_check, ...
            cfg_ci.B_end);
    end
    fprintf([lang.bca.secondary_criterion '\n']);
    
    % Preparation of variables for the parfor loop and stability analysis.
    B_vector_ci = cfg_ci.B_start:cfg_ci.B_step:cfg_ci.B_end;
    % Stability matrix is dynamically sized based on num_metrics
    stability_matrix_ci = NaN(2, num_metrics, numel(B_vector_ci));
    final_i_ci = 0;
    converged_ci = false;
    
    overall_stability_ci = NaN(1, numel(B_vector_ci));
    
    % temp vector is now sized
    temp_stability_ci_vector = zeros(1, num_metrics * 2);
    
    % Main loop: Iterates over different numbers of bootstrap samples (B).
    for i = 1:numel(B_vector_ci)
        B_ci_current = B_vector_ci(i);
        fprintf([' -> ' lang.bca.checking_stability '\n'], B_ci_current, cfg_ci.n_trials);
        % Reset temp vector for each B value
        temp_stability_ci_vector(:) = 0; 
    
        % --- Parallel Worker Limit ---
        % Allows external control of worker count for nested parallelism scenarios.
        % When config.num_workers is set, limits the parfor to that many workers.
        % Default: Uses all available pool workers.
        pool = gcp('nocreate');
        if isfield(config, 'num_workers') && isnumeric(config.num_workers) && config.num_workers > 0
            parfor_limit = min(pool.NumWorkers, config.num_workers);
        else
            parfor_limit = pool.NumWorkers;
        end

        % Parallel loop to calculate stability for all metrics and effect sizes.
        % Loop count is dynamic (num_metrics * 2)
        parfor (metric_idx = 1:(num_metrics * 2), parfor_limit)
            % Each iteration gets its own reproducible substream.
            s_worker = s; 
            s_worker.Substream = metric_idx;

            % Logic to determine effect type and metric index
            is_delta = metric_idx <= num_metrics;
            actual_metric_idx = mod(metric_idx-1, num_metrics) + 1;
            stability_all_pairs = zeros(num_pairs, 1);

            % Loop over all dataset pairs.
            for k = 1:num_pairs
                idx1=p_pair_idx_all(k,1); idx2=p_pair_idx_all(k,2);
                % Access data using dynamic index
                data_x_orig=p_all_data{actual_metric_idx}(:,idx1); 
                data_y_orig=p_all_data{actual_metric_idx}(:,idx2);
                ci_widths_trial = zeros(cfg_ci.n_trials, 1);
                
                % --- Jackknife Statistics ---
                % Pre-calculate Jackknife metrics once per pair. 
                % These are deterministic and independent of the bootstrap trials, 
                % allowing significant computational savings compared to re-evaluating them inside the loop.
                
                jack_d = []; jack_r = [];
                % We only need to compute Jackknife for the relevant metric (d or r)
                if is_delta
                   jack_d = HERA.stats.jackknife(data_x_orig, data_y_orig, 'delta');
                   mean_jack = mean(jack_d);
                   jack_vals = jack_d;
                else
                   jack_r = HERA.stats.jackknife(data_x_orig, data_y_orig, 'rel');
                   mean_jack = mean(jack_r);
                   jack_vals = jack_r;
                end
                
                a_num = sum((mean_jack - jack_vals).^3);
                a_den = 6 * (sum((mean_jack - jack_vals).^2)).^(3/2);
                if a_den == 0, a = 0; else, a = a_num / a_den; end
                if ~isfinite(a), a = 0; end
            
                % Repeats the BCa calculation 'n_trials' times to assess the stability of the CI width.
                for t = 1:cfg_ci.n_trials
                    
                    % --- Vectorized Bootstrap Sampling ---
                    % Generate all bootstrap indices and sampled data in a single operation [N x B].
                    % This avoids the overhead of repeated function calls in a loop.
                    boot_indices = randi(s_worker, num_probanden, [num_probanden, B_ci_current]);
                    
                    % Direct indexing is significantly faster than generating intermediate matrices.
                    boot_x = data_x_orig(boot_indices);
                    boot_y = data_y_orig(boot_indices);
                    
                    % --- Robust Handling of Missing Data (NaN) ---
                    % Standard vectorized operations assume consistent sample sizes. 
                    % If the original data contains NaNs, pairwise exclusion results in varying 
                    % valid sample sizes per bootstrap iteration, preventing simple 3D vectorization.
                    %
                    % Strategy:
                    % 1. Check for NaNs in the source data.
                    % 2. Fast Path: If clean, use highly efficient vectorized calculations.
                    % 3. Robust Path: If NaNs exist, fall back to a loop to correctly handle varying N.
                    
                    has_nans = any(isnan(data_x_orig)) || any(isnan(data_y_orig));
                    
                    if ~has_nans
                         % FAST PATH: No NaNs present.
                         % Calculate statistics for all B iterations simultaneously using 3D expansion.
                         if is_delta
                             boot_stats = HERA.stats.cliffs_delta(boot_x, boot_y);
                         else
                             boot_stats = HERA.stats.relative_difference(boot_x, boot_y);
                         end
                    else
                         % ROBUST PATH: NaNs present.
                         % Iterate through each bootstrap sample and perform pairwise exclusion.
                         boot_stats = zeros(B_ci_current, 1);
                         for b = 1:B_ci_current
                             bx = boot_x(:, b); 
                             by = boot_y(:, b);
                             valid = ~isnan(bx) & ~isnan(by);
                             if is_delta
                                 boot_stats(b) = HERA.stats.cliffs_delta(bx(valid), by(valid));
                             else
                                 boot_stats(b) = HERA.stats.relative_difference(bx(valid), by(valid));
                             end
                         end
                    end
                    
                    % Ensure row vector
                    boot_stats = boot_stats(:)'; 

                    % Calculate BCa correction factors (z0 for bias).
                    % Access effect size using dynamic index
                    if is_delta, theta_hat = p_d_vals_all(k, actual_metric_idx);
                    else, theta_hat = p_rel_vals_all(k, actual_metric_idx); end
                    
                    z0 = norminv(sum(boot_stats < theta_hat) / B_ci_current);
                    if ~isfinite(z0), z0 = 0; end 
    
                    % Calculate BCa interval limits as percentiles.
                    z1 = norminv(alpha_level / 2); z2 = norminv(1 - alpha_level / 2);
                    a1 = normcdf(z0 + (z0 + z1) / (1 - a * (z0 + z1)));
                    a2 = normcdf(z0 + (z0 + z2) / (1 - a * (z0 + z2)));
                    if isnan(a1), a1 = alpha_level / 2; end 
                    if isnan(a2), a2 = 1 - alpha_level / 2; end
                    
                    % Read confidence interval from the sorted bootstrap distribution and save the width.
                    sorted_boots = sort(boot_stats);
                    ci_lower = sorted_boots(max(1, floor(B_ci_current * a1)));
                    ci_upper = sorted_boots(min(B_ci_current, ceil(B_ci_current * a2)));
                    ci_widths_trial(t) = ci_upper - ci_lower;
                end
                
                % Calculate stability as a relative dispersion measure (IQR / Median) of the CI widths.
                med_w = median(ci_widths_trial); iqr_w = iqr(ci_widths_trial);
                if med_w == 0
                    if iqr_w == 0 
                        stability_all_pairs(k) = 0; 
                    else
                        stability_all_pairs(k) = Inf; 
                    end
                else
                    stability_all_pairs(k) = iqr_w / abs(med_w); 
                end
            end
            % Median of the stability values over all pairs as a measure for the current metric.
            temp_stability_ci_vector(metric_idx) = median(stability_all_pairs, 'omitnan');
        end
        % Save the stability values for each metric and effect size.
        stability_matrix_ci(1, :, i) = temp_stability_ci_vector(1:num_metrics);
        stability_matrix_ci(2, :, i) = temp_stability_ci_vector(num_metrics+1 : num_metrics*2);
        final_i_ci = i;
        
        % Save the average stability value for the overall convergence check.
        overall_stability_ci(i) = nanmean(stability_matrix_ci(:, :, i), 'all');
    
        % Convergence check via Helper Function
        [converged_ci, stats] = HERA.stats.check_convergence(overall_stability_ci(1:i), cfg_ci);
        
        rel_imp = stats.improvement;
        
        % Logging
        if use_robust_convergence_ci
            % Only print if enough steps have passed for logic to kick in
            if ~isnan(rel_imp)
                if abs(rel_imp) < cfg_ci.convergence_tolerance
                    fprintf(['    ' lang.bca.convergence_run_info '\n'],...
                    rel_imp * 100, stats.streak, cfg_ci.convergence_streak_needed);
                else
                    fprintf(['    ' lang.bca.stability_change_info '\n'], rel_imp * 100);
                end
                
                if converged_ci
                    fprintf([lang.bca.convergence_reached '\n'], rel_imp * 100);
                    fprintf([lang.bca.stable_runs_info '\n'], cfg_ci.convergence_streak_needed);
                end
            end
        else
            % Simple convergence check logs
            if ~isnan(rel_imp) && i >= cfg_ci.min_steps_for_convergence_check
                 fprintf(['    ' lang.bca.stability_change_info '\n'], rel_imp * 100);    
                 if converged_ci
                    fprintf([lang.bca.convergence_reached '\n'], rel_imp * 100);
                 end
            end
        end

        if converged_ci
            break; % Break the loop if convergence is reached.
        end
    end % End of the for-loop for stability check.
    
    % Final B-value determination
    B_tested_vector_ci = B_vector_ci(1:final_i_ci);
    overall_stability_ci_plotted = overall_stability_ci(1:final_i_ci);
    
    elbow_indices = [];
    if converged_ci
        % If converged, the last tested B-value is used.
        selected_B_ci = B_tested_vector_ci(end);
        fprintf([lang.bca.convergence_result '\n'], selected_B_ci);
    else
        % Otherwise, finding 'elbow' via Helper Function
        fprintf([lang.bca.elbow_analysis_info '\n']);
        
        % Prepare curves matrix for the helper (steps x curves)
        num_curves = num_metrics * 2;
        curves_matrix = zeros(final_i_ci, num_curves);
        for k_sv = 1:num_metrics
             curves_matrix(:, k_sv) = squeeze(stability_matrix_ci(1, k_sv, 1:final_i_ci));
             curves_matrix(:, k_sv + num_metrics) = squeeze(stability_matrix_ci(2, k_sv, 1:final_i_ci));
        end
        
        [selected_B_ci, elbow_indices] = HERA.stats.find_elbow_point(B_tested_vector_ci, curves_matrix);
        fprintf([lang.bca.elbow_result '\n'], selected_B_ci);
    end
    
    % Store stability data for JSON export
    stability_data_ci = struct();
    stability_data_ci.B_vector = B_tested_vector_ci;
    stability_data_ci.global_stability = overall_stability_ci_plotted;
    stability_data_ci.detailed_stability = stability_matrix_ci(:, :, 1:final_i_ci);
    stability_data_ci.converged = converged_ci;
    stability_data_ci.elbow_indices = elbow_indices;
    B_ci = selected_B_ci; 

    %% 2. Generate graphics to visualize the convergence analysis via Helper
    [h_fig_bca_global, h_fig_bca_detailed] = HERA.plot.bca_convergence(...
        B_tested_vector_ci, overall_stability_ci_plotted, stability_matrix_ci(:, :, 1:final_i_ci), ...
        selected_B_ci, config, metric_names, graphics_dir, styles, lang, elbow_indices);
end

%% 3. Final calculation of BCa confidence intervals with optimal B_ci
% Calculates the final BCa confidence intervals for all pairs and metrics.
%
% Architecture:
%   a) Pre-compute deterministic components: Jackknife and acceleration factor 'a' 
%      are calculated once per pair (serial, reproducible, independent of B).
%   b) Memory-aware batch sizing: Splits B into chunks that fit in RAM.
%   c) Parallel bootstrap over batches: Each (metric, batch) gets a unique RNG substream.
%   d) Aggregation: Combines batch results and calculates BCa intervals serially.
%
% RNG Strategy:
%   - Unique substream per (metric, batch): metric_offset + b
%   - metric_offset = OFFSET_BASE_CI + (metric_idx - 1) * (num_batches + 100)
%   - The +100 buffer prevents collisions even if num_batches varies between metrics.
%
% Initialization of the output matrices.
ci_d_all = NaN(num_pairs, 2, num_metrics);
ci_r_all = NaN(num_pairs, 2, num_metrics);
z0_d_all = NaN(num_pairs, num_metrics);
a_d_all = NaN(num_pairs, num_metrics);
z0_r_all = NaN(num_pairs, num_metrics);
a_r_all = NaN(num_pairs, num_metrics);


% Dynamic batch sizing based on memory configuration.
if isfield(config, 'system') && isfield(config.system, 'target_memory')
     TARGET_MEMORY = config.system.target_memory;
else
     TARGET_MEMORY = 200;
end

% Get worker count for memory estimation (parfor broadcasts data to all workers).
if isfield(config, 'num_workers') && isnumeric(config.num_workers)
    num_workers = config.num_workers;
else
    num_workers = feature('numcores');
end

% Effective memory per batch accounts for parfor broadcast overhead.
effective_memory = TARGET_MEMORY / max(1, num_workers);

% Use substream offset to avoid overlap with stability analysis phase.
OFFSET_BASE_CI = 1000;

% Loop over the metrics for the final calculation.
for metric_idx = 1:num_metrics
    fprintf(lang.bca.calculating_final_ci, metric_names{metric_idx});
    
    % --- Pre-compute data and Jackknife for all pairs (serial, deterministic) ---
    pair_data_x = cell(num_pairs, 1);
    pair_data_y = cell(num_pairs, 1);
    pair_n_valid = zeros(num_pairs, 1);
    pair_jack_d = cell(num_pairs, 1);
    pair_jack_r = cell(num_pairs, 1);
    pair_theta_hat_d = zeros(num_pairs, 1);
    pair_theta_hat_r = zeros(num_pairs, 1);
    pair_a_d = zeros(num_pairs, 1);
    pair_a_r = zeros(num_pairs, 1);
    
    for k = 1:num_pairs
        i = p_pair_idx_all(k, 1);
        j = p_pair_idx_all(k, 2);
        data_x_orig = p_all_data{metric_idx}(:, i);
        data_y_orig = p_all_data{metric_idx}(:, j);
        
        % Pairwise NaN exclusion.
        valid_mask = ~isnan(data_x_orig) & ~isnan(data_y_orig);
        pair_data_x{k} = data_x_orig(valid_mask);
        pair_data_y{k} = data_y_orig(valid_mask);
        pair_n_valid(k) = sum(valid_mask);
        
        if pair_n_valid(k) >= 2
            % Pre-compute Jackknife (deterministic, only depends on original data).
            pair_jack_d{k} = HERA.stats.jackknife(pair_data_x{k}, pair_data_y{k}, 'delta');
            pair_jack_r{k} = HERA.stats.jackknife(pair_data_x{k}, pair_data_y{k}, 'rel');
            
            % Calculate acceleration factors from Jackknife.
            mean_jack_d = mean(pair_jack_d{k});
            a_num_d = sum((mean_jack_d - pair_jack_d{k}).^3);
            a_den_d = 6 * (sum((mean_jack_d - pair_jack_d{k}).^2)).^(3/2);
            if a_den_d == 0, pair_a_d(k) = 0; else, pair_a_d(k) = a_num_d / a_den_d; end
            if ~isfinite(pair_a_d(k)), pair_a_d(k) = 0; end
            
            mean_jack_r = mean(pair_jack_r{k});
            a_num_r = sum((mean_jack_r - pair_jack_r{k}).^3);
            a_den_r = 6 * (sum((mean_jack_r - pair_jack_r{k}).^2)).^(3/2);
            if a_den_r == 0, pair_a_r(k) = 0; else, pair_a_r(k) = a_num_r / a_den_r; end
            if ~isfinite(pair_a_r(k)), pair_a_r(k) = 0; end
            
            % Store original effect sizes.
            pair_theta_hat_d(k) = p_d_vals_all(k, metric_idx);
            pair_theta_hat_r(k) = p_rel_vals_all(k, metric_idx);
        end
    end
    
    % --- Memory-aware batch sizing for bootstrap ---
    max_n_valid = max(pair_n_valid);
    if max_n_valid < 2, max_n_valid = 2; end
    bytes_per_sample = max_n_valid * 8 * 2;  % x and y bootstrap samples
    total_memory_needed = (B_ci * bytes_per_sample) / (1024^2);
    
    % Smart batching: Use one batch if total memory fits, otherwise split.
    if total_memory_needed <= effective_memory
        BATCH_SIZE = B_ci;
    else
        BATCH_SIZE = max(100, floor((effective_memory * 1024^2) / bytes_per_sample));
    end
    num_batches = ceil(B_ci / BATCH_SIZE);
    
    % Unique offset for this metric to avoid substream collisions.
    metric_offset = OFFSET_BASE_CI + (metric_idx - 1) * (num_batches + 100);
    
    % --- Parallel bootstrap computation over batches ---
    % Each batch computes a portion of the B bootstrap samples.
    boot_d_all = cell(num_pairs, 1);
    boot_r_all = cell(num_pairs, 1);
    for k = 1:num_pairs
        boot_d_all{k} = zeros(1, B_ci);
        boot_r_all{k} = zeros(1, B_ci);
    end
    
    parfor b = 1:num_batches
        % Each batch gets its own reproducible substream.
        s_worker = s;
        s_worker.Substream = metric_offset + b;
        
        % Determine batch range.
        start_idx = (b - 1) * BATCH_SIZE + 1;
        end_idx = min(b * BATCH_SIZE, B_ci);
        batch_B = end_idx - start_idx + 1;
        
        % Local storage for this batch.
        batch_d = cell(num_pairs, 1);
        batch_r = cell(num_pairs, 1);
        
        for k = 1:num_pairs
            if pair_n_valid(k) < 2
                batch_d{k} = NaN(1, batch_B);
                batch_r{k} = NaN(1, batch_B);
                continue;
            end
            
            n_valid = pair_n_valid(k);
            data_x = pair_data_x{k};
            data_y = pair_data_y{k};
            
            % Generate bootstrap indices for this batch.
            boot_indices = randi(s_worker, n_valid, [n_valid, batch_B]);
            boot_x = data_x(boot_indices);
            boot_y = data_y(boot_indices);
            
            % Calculate effect sizes.
            batch_d{k} = HERA.stats.cliffs_delta(boot_x, boot_y);
            batch_r{k} = HERA.stats.relative_difference(boot_x, boot_y);
        end
        
        % Return batch results (will be aggregated outside parfor).
        boot_d_all_batch{b} = batch_d;
        boot_r_all_batch{b} = batch_r;
    end
    
    % --- Aggregate batch results ---
    for b = 1:num_batches
        start_idx = (b - 1) * BATCH_SIZE + 1;
        end_idx = min(b * BATCH_SIZE, B_ci);
        for k = 1:num_pairs
            boot_d_all{k}(start_idx:end_idx) = boot_d_all_batch{b}{k};
            boot_r_all{k}(start_idx:end_idx) = boot_r_all_batch{b}{k};
        end
    end
    
    % --- Calculate BCa intervals from aggregated bootstrap distributions ---
    temp_ci_d = NaN(num_pairs, 2);
    temp_ci_r = NaN(num_pairs, 2);
    temp_z0_d = NaN(num_pairs, 1);
    temp_a_d = NaN(num_pairs, 1);
    temp_z0_r = NaN(num_pairs, 1);
    temp_a_r = NaN(num_pairs, 1);
    
    z1 = norminv(alpha_level / 2); 
    z2 = norminv(1 - alpha_level / 2);
    
    for k = 1:num_pairs
        if pair_n_valid(k) < 2
            continue;
        end
        
        boot_d = boot_d_all{k};
        boot_r = boot_r_all{k};
        theta_hat_d = pair_theta_hat_d(k);
        theta_hat_r = pair_theta_hat_r(k);
        a_d = pair_a_d(k);
        a_r = pair_a_r(k);
        
        % BCa for Cliff's Delta.
        z0_d = norminv(sum(boot_d < theta_hat_d) / B_ci);
        if ~isfinite(z0_d), z0_d = 0; end
        a1_d = normcdf(z0_d + (z0_d + z1) / (1 - a_d * (z0_d + z1)));
        a2_d = normcdf(z0_d + (z0_d + z2) / (1 - a_d * (z0_d + z2)));
        if isnan(a1_d), a1_d = alpha_level / 2; end
        if isnan(a2_d), a2_d = 1 - alpha_level / 2; end
        sorted_d = sort(boot_d);
        temp_ci_d(k, :) = [sorted_d(max(1, floor(B_ci * a1_d))), sorted_d(min(B_ci, ceil(B_ci * a2_d)))];
        temp_z0_d(k) = z0_d;
        temp_a_d(k) = a_d;
        
        % BCa for Relative Difference.
        z0_r = norminv(sum(boot_r < theta_hat_r) / B_ci);
        if ~isfinite(z0_r), z0_r = 0; end
        a1_r = normcdf(z0_r + (z0_r + z1) / (1 - a_r * (z0_r + z1)));
        a2_r = normcdf(z0_r + (z0_r + z2) / (1 - a_r * (z0_r + z2)));
        if isnan(a1_r), a1_r = alpha_level / 2; end
        if isnan(a2_r), a2_r = 1 - alpha_level / 2; end
        sorted_r = sort(boot_r);
        temp_ci_r(k, :) = [sorted_r(max(1, floor(B_ci * a1_r))), sorted_r(min(B_ci, ceil(B_ci * a2_r)))];
        temp_z0_r(k) = z0_r;
        temp_a_r(k) = a_r;
    end
    
    % Transfer the results from the temporary variables to the final output matrices.
    ci_d_all(:, :, metric_idx) = temp_ci_d;
    ci_r_all(:, :, metric_idx) = temp_ci_r;
    z0_d_all(:, metric_idx) = temp_z0_d;
    a_d_all(:, metric_idx) = temp_a_d;
    z0_r_all(:, metric_idx) = temp_z0_r;
    a_r_all(:, metric_idx) = temp_a_r;
end

%% 4. Output of BCa Correction Factors
% Prints a formatted table with the summarized correction factors to the console.
fprintf(lang.bca.correction_factors.header, num_pairs);

% Define headers and prepare for dynamic formatting.
header_parts = {lang.bca.correction_factors.factor, lang.csv.headers.median, lang.csv.headers.mean, lang.csv.headers.min, lang.csv.headers.max};
alignments = {'l', 'c', 'c', 'c', 'c'}; % First column left, others centered.

for metric_idx = 1:num_metrics
    fprintf('\n%s:\n', metric_names{metric_idx});

    % Collect and format all data for this metric into a cell array.
    z0_d = z0_d_all(:, metric_idx);
    a_d = a_d_all(:, metric_idx);
    z0_r = z0_r_all(:, metric_idx);
    a_r = a_r_all(:, metric_idx);

    table_data = cell(4, numel(header_parts));
    table_data(1, :) = {lang.bca.correction_factors.bias_delta, sprintf('%+.3f', median(z0_d)), sprintf('%+.3f', ...
        mean(z0_d)), sprintf('%+.3f', min(z0_d)), sprintf('%+.3f', max(z0_d))};
    table_data(2, :) = {lang.bca.correction_factors.skew_delta, sprintf('%+.3f', median(a_d)), sprintf('%+.3f', ...
        mean(a_d)), sprintf('%+.3f', min(a_d)), sprintf('%+.3f', max(a_d))};
    table_data(3, :) = {lang.bca.correction_factors.bias_rel, sprintf('%+.3f', median(z0_r)), sprintf('%+.3f', ...
        mean(z0_r)), sprintf('%+.3f', min(z0_r)), sprintf('%+.3f', max(z0_r))};
    table_data(4, :) = {lang.bca.correction_factors.skew_rel, sprintf('%+.3f', median(a_r)), sprintf('%+.3f', ...
        mean(a_r)), sprintf('%+.3f', min(a_r)), sprintf('%+.3f', max(a_r))};

    % Calculate dynamic column widths based on content.
    col_widths = cellfun(@strlength, header_parts);
    for r = 1:size(table_data, 1)
        for c = 1:size(table_data, 2)
            col_widths(c) = max(col_widths(c), strlength(table_data{r, c}));
        end
    end
    col_widths = col_widths + 2; % Add 2 spaces for padding.

    % Print the formatted table.
    header_line = strjoin(arrayfun(@(c) format_text(header_parts{c}, col_widths(c), alignments{c}), 1:numel(header_parts), 'UniformOutput', false), '|');
    fprintf('%s\n', header_line);
    fprintf('%s\n', repmat('-', 1, strlength(header_line)));

    for r = 1:size(table_data, 1)
        row_line = strjoin(arrayfun(@(c) format_text(table_data{r, c}, col_widths(c), alignments{c}), 1:numel(header_parts), 'UniformOutput', false), '|');
        fprintf('%s\n', row_line);
    end
    fprintf('%s\n', repmat('-', 1, strlength(header_line)));
end

% Save BCa factors as a CSV file
fprintf(['\n' lang.bca.saving_csv '\n']);
[~, fName, fExt] = fileparts(lang.files.bca_factors);
fName = strrep(fName, '%s_', ''); 
csv_filename = fullfile(csv_dir, [fName, '_', ts, fExt]);

try
    % Attempt to open file for writing
    fileID = fopen(csv_filename, 'w');
    if fileID == -1
        error(lang.errors.file_open_error, csv_filename); 
    end
    
    % Header for the CSV file
    header = {'Metric', lang.csv.headers.effect_size, lang.csv.headers.correction_factor, ...
              lang.csv.headers.median, lang.csv.headers.mean, ...
              lang.csv.headers.min, lang.csv.headers.max};
    fprintf(fileID, '%s,%s,%s,%s,%s,%s,%s\n', header{:});
    
    % Loop over all metrics to write data
    for metric_idx = 1:numel(metric_names)
        metric_name = metric_names{metric_idx};
        
        % Extract data for the current metric
        z0_d = z0_d_all(:, metric_idx);
        a_d  = a_d_all(:, metric_idx);
        z0_r = z0_r_all(:, metric_idx);
        a_r  = a_r_all(:, metric_idx);
        
        % Write rows to the CSV file with explicit sign formatting
        fprintf(fileID, '%s,Cliff''s Delta,%s,%+.3f,%+.3f,%+.3f,%+.3f\n', ...
            metric_name, lang.bca.correction_factors.bias, median(z0_d), mean(z0_d), min(z0_d), max(z0_d));
        fprintf(fileID, '%s,Cliff''s Delta,%s,%+.3f,%+.3f,%+.3f,%+.3f\n', ...
            metric_name, lang.bca.correction_factors.skew, median(a_d), mean(a_d), min(a_d), max(a_d));
        fprintf(fileID, '%s,Relative Difference,%s,%+.3f,%+.3f,%+.3f,%+.3f\n', ...
            metric_name, lang.bca.correction_factors.bias, median(z0_r), mean(z0_r), min(z0_r), max(z0_r));
        fprintf(fileID, '%s,Relative Difference,%s,%+.3f,%+.3f,%+.3f,%+.3f\n', ...
            metric_name, lang.bca.correction_factors.skew, median(a_r), mean(a_r), min(a_r), max(a_r));
    end
    
    % Close the file successfully
    fclose(fileID);
    fprintf([lang.bca.csv_saved '\n'], csv_filename);

catch ME
    % Safety cleanup: Ensure file is closed if an error occurs
    if exist('fileID', 'var') && fileID ~= -1
        fclose(fileID); 
    end
    fprintf([lang.errors.file_save_error '\n'], ME.message);
end

%% 5. Generate histogram distributions via Helper Function
[h_fig_hist_z0, h_fig_hist_a, h_fig_hist_widths] = HERA.plot.bca_distributions(...
    z0_d_all, a_d_all, z0_r_all, a_r_all, ci_d_all, ci_r_all, metric_names, styles, lang, graphics_dir, ts, B_ci);

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

