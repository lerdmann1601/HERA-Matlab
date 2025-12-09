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
%      Iteratively checks stability and delegates convergence checks to `HERA.stats.check_convergence` 
%      and elbow detection to `HERA.stats.find_elbow_point`.
%   4. Final Threshold Calculation: 
%      Computes the final thresholds with the optimal B, where the relative threshold is capped by the SEM value.
%   5. Visualization:
%      Delegates the creation of convergence and distribution plots to `HERA.plot.threshold_convergence`
%      and `HERA.plot.threshold_distributions`.
%
% Inputs:
%   all_data      - Cell array of data matrices for each metric (1, 2, or 3).
%   num_probanden - Number of subjects.
%   config        - Struct with all configuration parameters (especially 'config.bootstrap_thresholds').
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
                d_vals_all(idx_pair, metric_idx) = HERA.stats.cliffs_delta(x, y);
                
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
        
        % Parallel loop to calculate stability for all metrics and effect sizes.
        % Loop from 1 to num_metrics * 2
        parfor metric_loop_idx = 1:(num_metrics * 2)
            
            % Each iteration gets its own reproducible substream.
            s_worker = s;
            s_worker.Substream = metric_loop_idx;
            % Dynamic logic to determine effect type and metric index
            is_delta = metric_loop_idx <= num_metrics;
            actual_metric_idx = mod(metric_loop_idx - 1, num_metrics) + 1;
        
            % Removes NaN values from the effect size vector before bootstrapping.
            if is_delta
                all_vals = abs(d_vals_all(:, actual_metric_idx));
            else
                all_vals = rel_vals_all(:, actual_metric_idx);
            end
            vals = all_vals(~isnan(all_vals)); % Keep only the valid values.
        
            % Protects against errors if an entire vector consisted only of NaNs.
            if isempty(vals)
                temp_stability_vector(metric_loop_idx) = 0; % No variance, so stable.
                continue; % Skip to the next iteration.
            end
                    
            % Repeats the threshold calculation n_trials times to capture variability.
            thr_trials = zeros(cfg_thr.n_trials, 1);
            for t = 1:cfg_thr.n_trials
                % Performs a percentile bootstrap: resamples the effect sizes. 
                % Calculates the median and repeats this B_current times to get a distribution of medians.
                bootstat = median(vals(randi(s_worker, numel(vals), [numel(vals), B_current])), 1);
                % The lower confidence interval of this distribution is an estimate of the threshold.
                ci_tmp = quantile(bootstat, [alpha_level / 2, 1 - alpha_level / 2]);
                thr_trials(t) = ci_tmp(1);
            end
            
            % Calculates stability as a relative measure of dispersion (IQR / Median).
            med_thr = median(thr_trials);
            iqr_thr = iqr(thr_trials);
            
            % Calculate stability as a relative measure of dispersion (coefficient of quartile variation: IQR / Median).
            % This metric is robust against outliers. A lower value indicates higher stability.
            if med_thr == 0
                % Handle the special case where the median of the threshold trials is zero.
                if iqr_thr == 0
                    % If IQR is also zero, all trials resulted in the exact same threshold of 0.
                    % This represents perfect stability.
                    temp_stability_vector(metric_loop_idx) = 0;
                else
                    % If the median is zero but there is some spread (IQR > 0), the relative
                    % stability is undefined. We set it to infinity to signal extreme instability.
                    temp_stability_vector(metric_loop_idx) = Inf; 
                end
            else
                % Standard case: Calculate stability as the Interquartile Range divided by the absolute median value.
                temp_stability_vector(metric_loop_idx) = iqr_thr / abs(med_thr);
            end
        end
        
        % Stores the stability values for each metric and effect size in the main matrix.
        % Dynamic slicing based on num_metrics
        stability_matrix(1, :, i) = temp_stability_vector(1:num_metrics); % Cliff's Delta stabilities
        stability_matrix(2, :, i) = temp_stability_vector(num_metrics+1 : end); % Relative Difference stabilities
        final_i = i; % Keep track of the last executed iteration index.
        
        % Calculate a single, overall stability value for the current B by averaging across all metrics.
        % This single value is used to check for convergence.
        overall_stability_thr(i) = nanmean(stability_matrix(:, :, i), 'all');
       
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
                 fprintf(['    ' lang.thresholds.stability_change_info '\n'], conv_stats.improvement * 100);
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
% Initialization of the output variables.
d_thresh = zeros(1, num_metrics); 
rel_thresh = zeros(1, num_metrics); 
all_bootstat_d = cell(1, num_metrics); 
all_bootstat_rel = cell(1, num_metrics);

% Calculates the final thresholds with the optimally determined B value.
for metric_idx = 1:num_metrics
    % Threshold for Cliff's Delta (cleaned of NaNs).
    d_vals_metric_raw = d_vals_all(:, metric_idx);
    
    % Use the absolute values before bootstrapping to find the lower confidence bound of median cliffs d.
    d_vals_metric_abs = abs(d_vals_metric_raw(~isnan(d_vals_metric_raw)));
    
    if ~isempty(d_vals_metric_abs)
        % Resample the absolute values with replacement, calculate median, and repeat 'selected_B' times.
        bootstat_d = median(d_vals_metric_abs(randi(numel(d_vals_metric_abs), [numel(d_vals_metric_abs), selected_B])), 1);
        
        % The threshold is the lower quantile (e.g., 2.5th percentile) of this distribution of medians. 
        d_thresh(metric_idx) = quantile(bootstat_d, alpha_level / 2);
        all_bootstat_d{metric_idx} = bootstat_d; 
    else
        d_thresh(metric_idx) = 0; % Fallback if no valid values.
        all_bootstat_d{metric_idx} = [];
    end

    % Threshold for the relative difference (cleaned of NaNs).
    rel_vals_metric_raw = rel_vals_all(:, metric_idx);
    rel_vals_metric = rel_vals_metric_raw(~isnan(rel_vals_metric_raw));

    if ~isempty(rel_vals_metric)
        bootstat_rel = median(rel_vals_metric(randi(numel(rel_vals_metric), [numel(rel_vals_metric), selected_B])), 1);
        rel_thresh(metric_idx) = quantile(bootstat_rel, alpha_level / 2);
        all_bootstat_rel{metric_idx} = bootstat_rel;
    else
        rel_thresh(metric_idx) = 0; % Fallback if no valid values.
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
