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
%      Finds the optimal bootstrap count (B) for stable thresholds through convergence or elbow analysis. 
%      Dynamically analyzes stability for 1, 2, or 3 metrics).
%   4. Final Threshold Calculation: 
%      Computes the final thresholds with the optimal B, where the relative threshold is capped by the SEM value.
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
% Author:   Lukas von Erdmannsdorff
% Date:     12.10.2025
% Version:  1

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
                gt = sum(x > y', 'all'); % Number of cases where x > y.
                lt = sum(x < y', 'all'); % Number of cases where x < y.
                d_vals_all(idx_pair, metric_idx) = (gt - lt) / (n_valid^2);
                
                % Calculation of the relative mean difference.
                mean_i = mean(x);
                mean_j = mean(y);
                if (mean_i + mean_j) == 0
                    rel_vals_all(idx_pair, metric_idx) = 0;
                else
                    rel_vals_all(idx_pair, metric_idx) = abs(mean_i - mean_j) / abs(mean([mean_i, mean_j]));
                end
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
        if use_robust_convergence_thr
            % Robust method with smoothing and a "streak" check to avoid premature termination due to random fluctuations.
            warmup_period = cfg_thr.min_steps_for_convergence_check;
            smoothing_window = cfg_thr.smoothing_window;
            streak_needed = cfg_thr.convergence_streak_needed;
            
            % The check only starts after a sufficient number of iterations have passed (warmup + smoothing window).
            if i >= warmup_period + smoothing_window
                % Apply a moving average to the overall stability curve to reduce noise.
                smoothed_stability = movmean(overall_stability_thr(1:i), smoothing_window, 'omitnan');
                prev_stab_smooth = smoothed_stability(end-1);
                curr_stab_smooth = smoothed_stability(end);
                
                % Calculate the relative improvement in stability from the previous step to the current one.
                if isnan(curr_stab_smooth) || isnan(prev_stab_smooth), rel_imp = 1.0; % Assume high improvement if NaN occurs, to continue the loop.
                elseif isinf(curr_stab_smooth), rel_imp = -1.0; % A move to infinity is a negative improvement.
                elseif prev_stab_smooth == 0
                    if curr_stab_smooth == 0, rel_imp = 0; % No change from perfect stability.
                    else, rel_imp = -1.0; end % Any move away from zero is a negative improvement.
                elseif isinf(prev_stab_smooth)
                    if isinf(curr_stab_smooth), rel_imp = 0; % No change from infinite instability.
                    else, rel_imp = 1.0; end % Any move away from infinity is a 100% improvement.
                else
                    % Standard calculation: (old_value - new_value) / old_value
                    rel_imp = (prev_stab_smooth - curr_stab_smooth) / prev_stab_smooth;
                end

                % Check if the improvement is below the tolerance threshold.
                if abs(rel_imp) < cfg_thr.convergence_tolerance
                    % If it is, increment the counter for consecutive stable runs.
                    convergence_streak_counter_thr = convergence_streak_counter_thr + 1; 
                    fprintf(['    ' lang.thresholds.convergence_run_info '\n'], rel_imp * 100, convergence_streak_counter_thr, streak_needed);
                else
                    % If the change is still large, reset the streak counter.
                    convergence_streak_counter_thr = 0; 
                    fprintf(['    ' lang.thresholds.stability_change_info '\n'], rel_imp * 100);
                end
                
                % If the number of consecutive stable runs meets the required streak, convergence is reached.
                if convergence_streak_counter_thr >= streak_needed 
                    fprintf([lang.thresholds.convergence_reached '\n'], rel_imp * 100);
                    fprintf([lang.thresholds.stable_runs_info '\n'], streak_needed);
                    converged = true; 
                end
            end
        else
            % Simple method without smoothing, which is faster but more sensitive to noise.
            if i >= cfg_thr.min_steps_for_convergence_check
                prev_stab = overall_stability_thr(i - 1);
                curr_stab = overall_stability_thr(i);
            
                % Calculate relative improvement, handling edge cases as in the robust method.
                if isnan(curr_stab) || isnan(prev_stab), rel_imp = 1.0;
                elseif isinf(curr_stab), rel_imp = -1.0;
                elseif prev_stab == 0
                    if curr_stab == 0, rel_imp = 0; else, rel_imp = -1.0; end
                elseif isinf(prev_stab)
                    if isinf(curr_stab), rel_imp = 0; else, rel_imp = 1.0; end
                else
                    rel_imp = (prev_stab - curr_stab) / prev_stab;
                end
                
                fprintf(['    ' lang.thresholds.stability_change_info '\n'], rel_imp * 100);

                % If the relative improvement is below the tolerance, convergence is reached immediately.
                if abs(rel_imp) < cfg_thr.convergence_tolerance
                    fprintf([lang.thresholds.convergence_reached '\n'], rel_imp * 100);
                    converged = true;
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
    else
        % Otherwise, the "elbow" of the curve is determined as the optimal point.
        fprintf([lang.thresholds.elbow_analysis_info '\n']);
        elbow_indices = zeros(1, num_metrics * 2);
        
        % Dynamically build stability_vectors cell array
        stability_vectors = cell(1, num_metrics * 2);
        for k_vec = 1:(num_metrics * 2)
            is_delta = k_vec <= num_metrics;
            actual_metric_idx = mod(k_vec - 1, num_metrics) + 1;
            if is_delta
                stability_vectors{k_vec} = squeeze(stability_matrix(1, actual_metric_idx, :))';
            else
                stability_vectors{k_vec} = squeeze(stability_matrix(2, actual_metric_idx, :))';
            end
        end
        
        % Loop through all stability curves.
        for k_elbow = 1:(num_metrics * 2)
            x_values = B_tested_vector(:);
            y_values = stability_vectors{k_elbow}(:);
            
            % Skip if there's no data or variation to analyze.
            if any(isnan(y_values)) || numel(unique(y_values)) < 2
                % If the curve is flat (e.g., stability as IQR/Median is 0), it's stable immediately.
                % Set elbow to the first index.
                elbow_indices(k_elbow) = 1; 
                continue;
            end
            
            % Normalize the data to a [0, 1] range for consistent distance calculation.
            x_norm = (x_values - min(x_values)) / (max(x_values) - min(x_values));
            y_norm = (y_values - min(y_values)) / (max(y_values) - min(y_values));
            
            % Handle cases where normalization might result in NaN (e.g., if max equals min).
            if all(isnan(x_norm)), x_norm(:) = 0; end
            if all(isnan(y_norm)), y_norm(:) = 0; end
            
            % Find the point with the maximum distance to the line connecting the first and last points.
            line_vec = [x_norm(end) - x_norm(1); y_norm(end) - y_norm(1)];
            vec_from_first = [x_norm - x_norm(1), y_norm - y_norm(1)];
            cross_prod = vec_from_first(:, 1) * line_vec(2) - vec_from_first(:, 2) * line_vec(1);
            
            [~, elbow_idx] = max(abs(cross_prod) / norm(line_vec));
            
            if isempty(elbow_idx), elbow_idx = numel(x_values); end
            elbow_indices(k_elbow) = elbow_idx;
        end        
        % The selected B is the maximum of all determined elbow points to ensure stability for all curves.
        selected_B = max(B_tested_vector(elbow_indices));
        fprintf([lang.thresholds.elbow_result '\n'], selected_B);
    end
    % Store stability data for JSON export
    stability_data_thr = struct();
    stability_data_thr.B_vector = B_tested_vector;
    stability_data_thr.global_stability = overall_stability_thr_plotted;
    stability_data_thr.detailed_stability = stability_matrix;

%% 4. Creates and saves plots to visualize the convergence analysis
% Use the passed global styles.
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');

% Plot 1: Global Convergence Curve
% This plot shows the overall stability (mean across all metrics) that was used for the convergence check.
h_fig_thr_global = figure('Name', lang.plots.titles.threshold_convergence_global, 'Color', styles.colors.background, 'Visible', 'off');
tcl_thr_global = tiledlayout(1, 1, 'Padding', 'compact');
ax_global = nexttile;

% Plot the raw global stability curve.
p1_global = plot(ax_global, B_tested_vector, overall_stability_thr_plotted * 100, '-o', 'LineWidth', 1.5, 'MarkerSize', 5, ...
    'Color', styles.colors.blue_marker, 'MarkerFaceColor', styles.colors.blue_marker);
grid(ax_global, 'on'); 
box(ax_global, 'on'); 
hold(ax_global, 'on');
set(ax_global, 'Color', styles.colors.background, 'GridColor', styles.colors.grid_color);

% Set title for the global plot.
fig_title_str_global = sprintf([lang.plots.titles.threshold_convergence_long_n_g], config.bootstrap_thresholds.n_trials);
title(tcl_thr_global, fig_title_str_global, 'FontSize', styles.font.title, 'FontWeight', 'bold', 'Color', styles.colors.text);

% Set axis labels and limits.
xlim(ax_global, [min(B_tested_vector), max(B_tested_vector) * 1.1]);
xlabel(ax_global, lang.plots.xlabels.bootstraps, 'FontSize', styles.font.label, 'Color', styles.colors.text);
ylabel(ax_global, lang.plots.ylabels.stability, 'FontSize', styles.font.label, 'Color', styles.colors.text); 
set(ax_global, 'FontSize', styles.font.tick, 'XColor', styles.colors.text, 'YColor', styles.colors.text);

% If robust convergence was used, also plot the smoothed curve (as this was used for the check).
if use_robust_convergence_thr
    smoothing_window = config.bootstrap_thresholds.smoothing_window;
    smoothed_curve_plotted_global = movmean(overall_stability_thr_plotted, smoothing_window, 'omitnan');
    p2_global = plot(ax_global, B_tested_vector, smoothed_curve_plotted_global * 100, '-', 'LineWidth', 1.5, 'Color', [0.8500 0.3250 0.0980]);
    lgd_global = legend(ax_global, [p1_global, p2_global], {lang.plots.legend.unsmoothed, lang.plots.legend.smoothed}, 'Location', 'best', ...
        'FontSize', styles.font.small_text);
    set(lgd_global, 'Color', styles.colors.background, 'TextColor', styles.colors.text, 'EdgeColor', styles.colors.text);
    lgd_global.ItemTokenSize = [15, 18];
end

% Add a marker for the optimal B value.
selected_B_idx_for_plot_global = find(B_tested_vector == selected_B, 1);
if ~isempty(selected_B_idx_for_plot_global)
    x_pos_global = selected_B;
    y_pos_global = overall_stability_thr_plotted(selected_B_idx_for_plot_global) * 100;

    plot(ax_global, x_pos_global, y_pos_global, 'x', 'Color', styles.colors.red_marker, 'MarkerSize', 10, 'LineWidth', 2, 'HandleVisibility', 'off');
    % Checks if the marker is too close to the right edge.
    current_xlim_global = xlim(ax_global);
    if x_pos_global > (current_xlim_global(1) + (current_xlim_global(2) - current_xlim_global(1)) * 0.85)
        % Extend the x-axis to make space for the text.
        xlim(ax_global, [current_xlim_global(1), x_pos_global * 1.2]);
    end
    
    text(ax_global, x_pos_global, y_pos_global, sprintf(lang.plots.misc.optimal_b_text, selected_B), ...
         'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', styles.font.small_text, 'Color', styles.colors.red_marker);
end
hold(ax_global, 'off');

% Save global graphic.
[~, fName, fExt] = fileparts(lang.files.convergence_thresholds_global);
filename_global = fullfile(subfolder_thresholds, [fName, '_', ts, fExt]);
exportgraphics(tcl_thr_global, filename_global, 'Resolution', 300, 'Padding', 30);
fprintf([lang.thresholds.convergence_plot_saved '\n'], filename_global);


% Plot 2: Detailed Convergence Curves (per Metric) 
% This plot shows the individual stability curves for each metric (without smoothing).
h_fig_thr_detailed = figure('Name', lang.plots.titles.threshold_convergence, 'Color', styles.colors.background, 'Visible', 'off');

% Creates the Tiled Layout and applies the title to it.
num_metrics_actual = numel(metric_names);
tcl_thr_detailed = tiledlayout(2, num_metrics_actual, 'TileSpacing', 'compact', 'Padding', 'compact');
fig_title_str_detailed = sprintf([lang.plots.titles.threshold_convergence_long_n_d], config.bootstrap_thresholds.n_trials);
title(tcl_thr_detailed, fig_title_str_detailed, 'FontSize', styles.font.title, 'FontWeight', 'bold', 'Color', styles.colors.text);

% Initialization for the loop.
num_effect_types = 2;
effect_type_names = {'Cliff''s Delta', 'Rel Diff'};
% Define the color for the local elbow (same as smoothed curve in global plot)
color_local_elbow = [0.8500 0.3250 0.0980]; 

% Loop through effect types (rows) and metrics (columns).
for es_type = 1:num_effect_types
    for actual_metric_idx = 1:num_metrics_actual
        ax = nexttile;
        
        % Prepare plot handles and names for the dynamic legend 
        handles_to_legend = [];
        names_to_legend = {};
        
        current_stability_curve_plotted = squeeze(stability_matrix(es_type, actual_metric_idx, :));
        
        % Plot the raw stability curve.
        p1 = plot(ax, B_tested_vector, current_stability_curve_plotted * 100, '-o', 'LineWidth', 1.5, 'MarkerSize', 5, ...
            'Color', styles.colors.blue_marker, 'MarkerFaceColor', styles.colors.blue_marker, ...
            'DisplayName', lang.plots.legend.unsmoothed);
        handles_to_legend(end+1) = p1;
        names_to_legend{end+1} = lang.plots.legend.unsmoothed;
            
        grid(ax, 'on'); 
        box(ax, 'on'); 
        hold(ax, 'on');
        set(ax, 'Color', styles.colors.background, 'GridColor', styles.colors.grid_color);
        
        subplot_title_str = sprintf('%s - %s', metric_names{actual_metric_idx}, effect_type_names{es_type});
        title(ax, subplot_title_str, 'FontSize', styles.font.label, 'Color', styles.colors.text, 'Interpreter', 'none'); 
        
        % Sets initial axis limits.
        xlim_vals = [min(B_tested_vector), max(B_tested_vector)];
        xlim(ax, [xlim_vals(1), xlim_vals(2) * 1.1]);
        xlabel(ax, lang.plots.xlabels.bootstraps, 'FontSize', styles.font.label, 'Color', styles.colors.text);
        ylabel(ax, lang.plots.ylabels.stability, 'FontSize', styles.font.label, 'Color', styles.colors.text); 
        set(ax, 'FontSize', styles.font.tick, 'XColor', styles.colors.text, 'YColor', styles.colors.text);
        
        % Plot the LOCAL elbow (if convergence failed)
        if ~converged && ~isempty(elbow_indices)
            % Calculate the index for the 'elbow_indices' vector
            if es_type == 1
                k_elbow_idx = actual_metric_idx;
            else
                k_elbow_idx = actual_metric_idx + num_metrics_actual;
            end
            
            % Get the specific elbow index for this subplot
            local_elbow_idx = elbow_indices(k_elbow_idx);
            
            if local_elbow_idx <= numel(current_stability_curve_plotted)
                x_local = B_tested_vector(local_elbow_idx);
                y_local = current_stability_curve_plotted(local_elbow_idx) * 100;
        
                % Plot as a filled circle with a dashed line style for legend
                p_elbow = plot(ax, x_local, y_local, ':o', ... % Style includes line and marker
                    'MarkerFaceColor', color_local_elbow, ...
                    'MarkerEdgeColor', color_local_elbow, ... 
                    'MarkerSize', 6, 'HandleVisibility', 'on', ... % Slightly larger (p1 is 5)
                    'LineWidth', 1.5, ... % Match line width for legend
                    'DisplayName', lang.plots.legend.local_elbow); 
                
                % Add vertical line for the elbow 
                xline(ax, x_local, '--', 'Color', color_local_elbow, 'LineWidth', 1.0, 'HandleVisibility', 'off');

                % Add the handle to our dynamic legend
                handles_to_legend(end+1) = p_elbow;
                names_to_legend{end+1} = lang.plots.legend.local_elbow; 
            end
        end

        % Add a marker for the final, globally optimal B value (selected_B).
        selected_B_idx_for_plot = find(B_tested_vector == selected_B, 1);
        if ~isempty(selected_B_idx_for_plot)
            x_pos = selected_B;
            y_pos = current_stability_curve_plotted(selected_B_idx_for_plot) * 100;

            % Plot the 'x' marker, but hide it from the legend.
            plot(ax, x_pos, y_pos, 'x', 'Color', styles.colors.red_marker, 'MarkerSize', 10, 'LineWidth', 2, 'HandleVisibility', 'off');
            
            % Checks if the marker is too close to the right edge.
            current_xlim = xlim(ax);
            if x_pos > (current_xlim(1) + (current_xlim(2) - current_xlim(1)) * 0.85)
                % Extend the x-axis to make space for the text.
                xlim(ax, [current_xlim(1), x_pos * 1.2]);
            end
            
            % Add text label for the final B value (this explains the 'x')
            text(ax, x_pos, y_pos, sprintf(lang.plots.misc.optimal_b_text, selected_B), ...
                 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', styles.font.small_text, 'Color', styles.colors.red_marker);
        end
        
        % Create the final, dynamic legend
        lgd = legend(handles_to_legend, names_to_legend, 'Location', 'best');
        set(lgd, 'Color', styles.colors.background, 'TextColor', styles.colors.text, 'EdgeColor', styles.colors.text);
        lgd.FontSize = styles.font.small_text; % Use smaller font for subplot legends
        lgd.ItemTokenSize = [15, 18];
        
        hold(ax, 'off');
    end
end

% Save graphic.
[~, fName, fExt] = fileparts(lang.files.convergence_thresholds);
filename_detailed = fullfile(subfolder_thresholds, [fName, '_', ts, fExt]);
exportgraphics(tcl_thr_detailed, filename_detailed, 'Resolution', 300, 'Padding', 30);
fprintf([lang.thresholds.convergence_plot_saved '\n'], filename_detailed);
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

%% 6. Creates and saves the histogram distribution of the bootstrap thresholds
% Use the passed global styles.
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');
effect_type_names = {'Cliff''s Delta', 'Rel Diff'};
% Create the figure.
h_fig_hist_thr = figure('Name', lang.plots.titles.threshold_dist_name, 'Color', styles.colors.background, 'Visible', 'off');
% Dynamic tiled layout
tcl_thr = tiledlayout(2, num_metrics, 'TileSpacing', 'compact', 'Padding', 'compact');
sgtitle_str = sprintf([lang.plots.titles.threshold_dist '\n'], selected_B);
title(tcl_thr, sgtitle_str, 'FontSize', styles.font.title, 'FontWeight', 'bold', 'Color', styles.colors.text);
% Loop for the top row (Cliff's Delta).
for metric_idx = 1:num_metrics
    ax = nexttile;
    bootstat_d = all_bootstat_d{metric_idx};
    data = abs(bootstat_d);
    
    hold on;
    % Handle cases with no or little data variation for robust plotting.
    if isempty(data) || isscalar(unique(data))
        bin_center = unique(data);
        if isempty(bin_center), bin_center = 0; end
        bin_width = 0.02; 
        histogram(data, 'BinEdges', [bin_center - bin_width/2, bin_center + bin_width/2], 'Normalization', 'probability', ...
            'FaceColor', styles.colors.delta_face, 'EdgeColor', styles.colors.bar_edge);
        xlim([bin_center - bin_width*5, bin_center + bin_width*5]);
        xticks(sort([bin_center, bin_center - bin_width*2, bin_center + bin_width*2]));
    else
        % Dynamically determine "nice" bin edges and tick marks for the histogram.
        points_of_interest = [data(:); d_thresh(metric_idx)]; 
        min_val = min(points_of_interest);
        max_val = max(points_of_interest);
        data_range = max_val - min_val;
        
        if data_range > 1e-6
            num_ticks_target = 5;
            raw_step = data_range / num_ticks_target;
            power = 10^floor(log10(raw_step));
            norm_step = raw_step / power;
            if norm_step < 1.5, nice_step = 1 * power;
            elseif norm_step < 3.5, nice_step = 2 * power;
            elseif norm_step < 7, nice_step = 5 * power;
            else, nice_step = 10 * power; 
            end
            nice_min = floor(min_val / nice_step) * nice_step;
            nice_max = ceil(max_val / nice_step) * nice_step;
            ticks = nice_min:nice_step:nice_max;
        else
            ticks = unique(data);
            nice_step = 0.1 * abs(ticks(1)) + 0.01;
        end
        if isempty(ticks) || numel(ticks) < 2, ticks = linspace(min_val, max_val, 3); nice_step = ticks(2)-ticks(1); end
        
        bin_edges = (ticks(1) - nice_step/2):nice_step:(ticks(end) + nice_step/2);
        
        histogram(data, 'BinEdges', bin_edges, 'Normalization', 'probability', 'FaceColor', styles.colors.delta_face, 'EdgeColor', styles.colors.bar_edge);
        
        xlim([bin_edges(1), bin_edges(end)]);
        xticks(ticks);
    end    
    % Set final y-axis limits and ticks before adding text.
    ylim([0 1]);
    set(gca, 'YTick', 0:0.1:1);
    % Add a line and text.
    xline(d_thresh(metric_idx), '--', 'Color', styles.colors.holm_threshold, 'LineWidth', 2);
    text(d_thresh(metric_idx), 1, sprintf(lang.plots.misc.boot_thr, d_thresh(metric_idx)), 'Color', styles.colors.holm_threshold, ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'Rotation', 90, 'FontSize', styles.font.tick);
    
    grid on; box on; hold off;
    set(ax, 'Color', styles.colors.background, 'GridColor', styles.colors.grid_color);
    title(sprintf('%s - %s', metric_names{metric_idx}, effect_type_names{1}), 'FontSize', styles.font.label, ...
        'Color', styles.colors.text, 'Interpreter', 'none');
    set(gca, 'FontSize', styles.font.tick, 'XColor', styles.colors.text, 'YColor', styles.colors.text);
    ylabel(lang.plots.ylabels.rel_frequency, 'FontSize', styles.font.label, 'Color', styles.colors.text); 
    xlabel(lang.plots.xlabels.median_delta, 'FontSize', styles.font.label, 'Color', styles.colors.text);
end

% Loop for the bottom row (Relative Difference).
for metric_idx = 1:num_metrics
    ax = nexttile;
    bootstat_rel = all_bootstat_rel{metric_idx};
    data = bootstat_rel;
    
    hold on;
    if isempty(data) || isscalar(unique(data))
        bin_center = unique(data);
        if isempty(bin_center), bin_center = 0; end
        bin_width = max(abs(bin_center) * 0.1, 0.05);
        if bin_width == 0, bin_width = 0.1; end
        histogram(data, 'BinEdges', [bin_center - bin_width/2, bin_center + bin_width/2], 'Normalization', 'probability', ...
            'FaceColor', styles.colors.rel_face, 'EdgeColor', styles.colors.bar_edge);
        xticks(bin_center);
    else
        points_of_interest = [data(:); rel_thresh_b(metric_idx)]; 
        if rel_thresh(metric_idx) > rel_thresh_b(metric_idx)
            points_of_interest(end+1) = rel_thresh(metric_idx);
        end
        min_val = min(points_of_interest);
        max_val = max(points_of_interest);
        data_range = max_val - min_val;
        
        if data_range > 1e-6
            num_ticks_target = 5;
            raw_step = data_range / num_ticks_target;
            power = 10^floor(log10(raw_step));
            norm_step = raw_step / power;
            if norm_step < 1.5, nice_step = 1 * power;
            elseif norm_step < 3.5, nice_step = 2 * power;
            elseif norm_step < 7, nice_step = 5 * power;
            else, nice_step = 10 * power; 
            end
            nice_min = floor(min_val / nice_step) * nice_step;
            nice_max = ceil(max_val / nice_step) * nice_step;
            ticks = nice_min:nice_step:nice_max;
        else
            ticks = unique(data);
            nice_step = 0.1 * abs(ticks(1)) + 0.01;
        end
        if isempty(ticks) || numel(ticks) < 2, ticks = linspace(min_val, max_val, 3); nice_step = ticks(2)-ticks(1); end
        
        bin_edges = (ticks(1) - nice_step/2):nice_step:(ticks(end) + nice_step/2);
        
        histogram(data, 'BinEdges', bin_edges, 'Normalization', 'probability', 'FaceColor', styles.colors.rel_face, 'EdgeColor', styles.colors.bar_edge);
        
        xlim([bin_edges(1), bin_edges(end)]);
        xticks(ticks);
    end    
    % Set final y-axis limits and ticks before adding text.
    ylim([0 1]);
    set(gca, 'YTick', 0:0.1:1);
    % Add lines and text.
    xline(rel_thresh_b(metric_idx), '--', 'Color', styles.colors.holm_threshold, 'LineWidth', 2);
    text(rel_thresh_b(metric_idx), 1, sprintf(lang.plots.misc.boot_thr, rel_thresh_b(metric_idx)), 'Color', styles.colors.holm_threshold, ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'Rotation', 90, 'FontSize', styles.font.tick);
    if rel_thresh(metric_idx) > rel_thresh_b(metric_idx)
         xline(rel_thresh(metric_idx), '--', 'Color', styles.colors.sem_threshold, 'LineWidth', 2);
         text(rel_thresh(metric_idx), 1, sprintf(lang.plots.misc.sem_thr, rel_thresh(metric_idx)), 'Color', styles.colors.sem_threshold, ...
             'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'Rotation', 90, 'FontSize', styles.font.tick, 'FontWeight', 'bold');
    end
    
    grid on; box on; hold off;
    set(ax, 'Color', styles.colors.background, 'GridColor', styles.colors.grid_color);
    title(sprintf('%s - %s', metric_names{metric_idx}, effect_type_names{2}), ...
        'FontSize', styles.font.label, 'Color', styles.colors.text, 'Interpreter', 'none');
    set(gca, 'FontSize', styles.font.tick, 'XColor', styles.colors.text, 'YColor', styles.colors.text);
    xlabel(lang.plots.xlabels.median_rel_diff, 'FontSize', styles.font.label, 'Color', styles.colors.text);
    ylabel(lang.plots.ylabels.rel_frequency, 'FontSize', styles.font.label, 'Color', styles.colors.text); 
end
% Save graphic with padding.
[~, fName, fExt] = fileparts(lang.files.dist_bootstrap_thresholds);
filename = fullfile(subfolder_thresholds, [fName, '_', ts, fExt]);
exportgraphics(tcl_thr, filename, 'Resolution', 300, 'Padding', 30);
fprintf([lang.thresholds.histogram_plot_saved '\n'], filename);

%% 7. Creates and saves the histogram distribution of the original effect sizes
% Creates the figure.
h_fig_hist_raw = figure('Name', lang.plots.titles.raw_effect_dist_name, 'Color', styles.colors.background, 'Visible', 'off');
% Dynamic tiled layout
tcl_thr = tiledlayout(2, num_metrics, 'TileSpacing', 'compact', 'Padding', 'compact');
sgtitle_str = sprintf([lang.plots.titles.raw_effect_dist '\n']);
title(tcl_thr, sgtitle_str, 'FontSize', styles.font.title, 'FontWeight', 'bold', 'Color', styles.colors.text);

% Loop for the top row (Cliff's Delta).
for metric_idx = 1:num_metrics
    ax = nexttile;
    data = d_vals_all(:, metric_idx);
    data_clean = data(isfinite(data));
    unique_data = unique(data_clean);
    
    hold on;
    
    % Handle different data scenarios for robust plotting.
    if isempty(unique_data)
        % No data -> Empty plot.
    elseif isscalar(unique_data)
        % Exactly 1 value.
        bin_center = unique_data;
        bin_width = 0.1; 
        histogram(data, 'BinEdges', [bin_center - bin_width/2, bin_center + bin_width/2], 'Normalization', 'probability', ...
            'FaceColor', styles.colors.delta_face, 'EdgeColor', styles.colors.bar_edge);
        xlim([-1.1, 1.1]);
        xticks(sort(unique([-1, 0, 1, bin_center])));       
    elseif numel(unique_data) == 2
        % Exactly 2 values.
        val1 = unique_data(1);
        val2 = unique_data(2);
        bin_width = 0.2; 
        bin_edges = [val1 - bin_width/2, val1 + bin_width/2, val2 - bin_width/2, val2 + bin_width/2];
        histogram(data, 'BinEdges', bin_edges, 'Normalization', 'probability', ...
            'FaceColor', styles.colors.delta_face, 'EdgeColor', styles.colors.bar_edge);
        xlim([-1.1, 1.1]);
        xticks(sort([val1, val2, 0]));
    else
        % Dynamic ticks based on the raw data range.
        min_val = min(data_clean);
        max_val = max(data_clean);
        data_range = max_val - min_val;
        
        if data_range > 1e-6
            num_ticks_target = 5; raw_step = data_range / num_ticks_target;
            power = 10^floor(log10(raw_step)); norm_step = raw_step / power;
            if norm_step < 1.5, nice_step = 1 * power;
            elseif norm_step < 3.5, nice_step = 2 * power;
            elseif norm_step < 7, nice_step = 5 * power;
            else, nice_step = 10 * power; 
            end
            
            nice_min = floor(min_val / nice_step) * nice_step;
            nice_max = ceil(max_val / nice_step) * nice_step;
            ticks = nice_min:nice_step:nice_max;
        else, ticks = unique(data_clean); nice_step = 0.1; 
        end
        if isempty(ticks) || numel(ticks) < 2, ticks = linspace(min_val, max_val, 3); nice_step = ticks(2)-ticks(1); end

        bin_edges = (ticks(1) - nice_step/2):nice_step:(ticks(end) + nice_step/2);
        h_hist = histogram(data, 'BinEdges', bin_edges, 'Normalization', 'probability', ...
            'FaceColor', styles.colors.delta_face, 'EdgeColor', styles.colors.bar_edge);
        
        % Restrict the final view to the range [-1, 1].
        final_xlim_min = max(-1.1, bin_edges(1));
        final_xlim_max = min(1.1, bin_edges(end));
        
        % Adjust KDE to the final, restricted axis limits.
        evaluation_points = linspace(final_xlim_min, final_xlim_max, 200);
        [f, xi] = ksdensity(data_clean, evaluation_points, 'Bandwidth', 'normal-approx');
        % Overlay the kernel density estimate, scaled to the histogram height.
        if max(f) > 0 && ~isempty(h_hist.Values) && max(h_hist.Values) > 0
             plot(xi, f * (max(h_hist.Values)/max(f)), 'Color', styles.colors.kde_line, 'LineWidth', 1.5);
        end      
        % Set final view.
        xlim([final_xlim_min, final_xlim_max]);
        xticks(ticks);
    end
    
    % Final axis properties.
    grid on; box on; hold off;
    set(ax, 'Color', styles.colors.background, 'GridColor', styles.colors.grid_color);
    title(sprintf('%s - %s', metric_names{metric_idx}, effect_type_names{1}), 'FontSize', styles.font.label, ...
        'Color', styles.colors.text, 'Interpreter', 'none');
    set(gca, 'FontSize', styles.font.tick, 'XColor', styles.colors.text, 'YColor', styles.colors.text);
    ylabel(lang.plots.ylabels.rel_frequency, 'FontSize', styles.font.label, 'Color', styles.colors.text); 
    xlabel('Cliff''s Delta', 'FontSize', styles.font.label, 'Color', styles.colors.text);
    ylim([0 1]);
    set(gca, 'YTick', 0:0.1:1);
end

% Loop for the bottom row (Relative Difference).
for metric_idx = 1:num_metrics
    ax = nexttile;
    data = rel_vals_all(:, metric_idx);
    data_clean = data(isfinite(data));
    unique_data = unique(data_clean);
    
    hold on;

    if isempty(unique_data) || isscalar(unique_data)
        bin_center = unique(data);
        if isempty(bin_center), bin_center = 0; end
        bin_width = 0.1; 
        histogram(data, 'BinEdges', [bin_center - bin_width/2, bin_center + bin_width/2], 'Normalization', 'probability', ...
            'FaceColor', styles.colors.rel_face, 'EdgeColor', styles.colors.bar_edge);
        xlim([-0.1, 2.1]);
        xticks(sort(unique([0, 1, 2, bin_center])));

    elseif numel(unique_data) == 2
        val1 = unique_data(1); val2 = unique_data(2);
        bin_width = 0.1; 
        bin_edges = [val1 - bin_width/2, val1 + bin_width/2, val2 - bin_width/2, val2 + bin_width/2];
        histogram(data, 'BinEdges', bin_edges, 'Normalization', 'probability', 'FaceColor', styles.colors.rel_face, 'EdgeColor', styles.colors.bar_edge);
        xlim([-0.1, 2.1]);
        xticks(sort(unique([0, val1, val2, 2])));
    else
        % Dynamic ticks based on the raw data range.
        min_val = min(data_clean); max_val = max(data_clean);
        data_range = max_val - min_val;
        
        if data_range > 1e-6
            num_ticks_target = 5; raw_step = data_range / num_ticks_target;
            power = 10^floor(log10(raw_step)); norm_step = raw_step / power;
            if norm_step < 1.5, nice_step = 1 * power;
            elseif norm_step < 3.5, nice_step = 2 * power;
            elseif norm_step < 7, nice_step = 5 * power;
            else, nice_step = 10 * power; 
            end
            
            nice_min = floor(min_val / nice_step) * nice_step;
            nice_max = ceil(max_val / nice_step) * nice_step;
            ticks = nice_min:nice_step:nice_max;
        else, ticks = unique(data_clean); nice_step = 0.1; 
        end
        if isempty(ticks) || numel(ticks) < 2, ticks = linspace(min_val, max_val, 3); nice_step = ticks(2)-ticks(1); end

        bin_edges = (ticks(1) - nice_step/2):nice_step:(ticks(end) + nice_step/2);
        h_hist = histogram(data, 'BinEdges', bin_edges, 'Normalization', 'probability', ...
            'FaceColor', styles.colors.rel_face, 'EdgeColor', styles.colors.bar_edge);
        
        % Adjust KDE to the final, restricted axis limits.
        final_xlim_min = max(-0.1, bin_edges(1));
        final_xlim_max = min(2.1, bin_edges(end));
        
        evaluation_points = linspace(final_xlim_min, final_xlim_max, 200); % Uses the final limits.
        [f, xi] = ksdensity(data_clean, evaluation_points, 'Bandwidth', 'normal-approx');
        if max(f) > 0 && ~isempty(h_hist.Values) && max(h_hist.Values) > 0
            plot(xi, f * (max(h_hist.Values)/max(f)), 'Color', styles.colors.kde_line, 'LineWidth', 1.5);
        end
        
        % Set final view.
        xlim([final_xlim_min, final_xlim_max]);
        xticks(ticks);
    end
    
    % Final axis properties.
    grid on; box on; hold off;
    set(ax, 'Color', styles.colors.background, 'GridColor', styles.colors.grid_color);
    title(sprintf('%s - %s', metric_names{metric_idx}, effect_type_names{2}), ...
        'FontSize', styles.font.label, 'Color', styles.colors.text, 'Interpreter', 'none');
    set(gca, 'FontSize', styles.font.tick, 'XColor', styles.colors.text, 'YColor', styles.colors.text);
    xlabel('Relative Difference', 'FontSize', styles.font.label, 'Color', styles.colors.text);
    ylabel(lang.plots.ylabels.rel_frequency, 'FontSize', styles.font.label, 'Color', styles.colors.text); 
    ylim([0 1]);
    set(gca, 'YTick', 0:0.1:1);
end

% Save graphic with padding.
[~, fName, fExt] = fileparts(lang.files.dist_raw_effects);
filename = fullfile(subfolder_thresholds, [fName, '_', ts, fExt]);
exportgraphics(tcl_thr, filename, 'Resolution', 300, 'Padding', 30);
fprintf([lang.thresholds.effects_plot_saved '\n'], filename);
end
