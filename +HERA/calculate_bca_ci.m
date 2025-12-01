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
%   This function dynamically determines the optimal number of bootstrap samples (B) for bias-corrected and accelerated (BCa) confidence intervals.
%   Subsequently, it calculates the confidence intervals and the correction factors (z0, a) for effect sizes using the optimal B-value.
%
% Workflow:
%   1. Dynamic determination of the bootstrap count (B): 
%      Iterates over B-values and checks the stability of the CI widths to find the optimal value through convergence or an elbow analysis.
%      This works dynamically for 1, 2, or 3 metrics.
%   2. Creation of the convergence plot: 
%      Visualizes the stability analysis for each effect size type and metric and saves the graphic as a PNG file.
%   3. Final calculation of BCa confidence intervals: 
%      Calculates the final CIs for Cliff's Delta and the relative difference for all pairwise comparisons using the optimal B-value.
%   4. Calculation and output of correction factors: 
%      Determines the bias correction (z0) and skewness correction (a), outputs a summary in the console and saves the results to a CSV file.
%   5. Generation of distribution plots:
%      Creates and saves histograms for the distribution of the correction factors (z0, a) and the confidence interval widths.
%
% Inputs:
%   all_data      - Cell array with the data matrices for each metric (1, 2, or 3 cells).
%   d_vals_all    - Matrix of the original Cliff's Delta effect sizes.
%   rel_vals_all  - Matrix of the original relative differences.
%   pair_idx_all  - Matrix of the indices for all pairwise comparisons.
%   num_probanden - Number of subjects.
%   config        - Struct containing configuration parameters (esp. config.bootstrap_ci).
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
%   B_ci     - Optimally determined number of bootstrap samples.
%   ci_d_all - 3D array of the BCa confidence intervals for Cliff's Delta.
%   ci_r_all - 3D array of the BCa confidence intervals for the relative difference.
%   z0_d_all - Matrix of the bias correction factor (z0) for Cliff's Delta.
%   a_d_all  - Matrix of the skewness correction factor (a) for Cliff's Delta.
%   z0_r_all - Matrix of the bias correction factor (z0) for the relative difference.
%   a_r_all  - Matrix of the skewness correction factor (a) for the relative difference.
%   stability_data_ci - Struct with convergence curve data for JSON export.
%   h_fig_bca_global  - Handle of the global convergence graphic.
%   h_fig_bca_detailed- Handle of the detailed convergence graphic.
%   h_fig_hist_z0     - Handle of the z0 distribution graphic.
%   h_fig_hist_a      - Handle of the a distribution graphic.
%   h_fig_hist_widths - Handle of the CI width distribution graphic.
%
% Author:   Lukas von Erdmannsdorff
% Date:     12.10.2025
% Version:  1


%% 1. Dynamic determination of the optimal bootstrap count (B)
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

% Preparation of variables for the parfor loops to avoid broadcasting the entire workspace.
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
    % Initialization for convergence check.
    convergence_streak_counter_ci = 0;
    overall_stability_ci = NaN(1, numel(B_vector_ci));
    
    % temp vector is now sized
    temp_stability_ci_vector = zeros(1, num_metrics * 2);
    
    % Main loop: Iterates over different numbers of bootstrap samples (B).
    for i = 1:numel(B_vector_ci)
        B_ci_current = B_vector_ci(i);
        fprintf([' -> ' lang.bca.checking_stability '\n'], B_ci_current, cfg_ci.n_trials);
        % Reset temp vector for each B value
        temp_stability_ci_vector(:) = 0; 
    
        % Parallel loop to calculate stability for all metrics and effect sizes.
        % Loop count is dynamic (num_metrics * 2)
        parfor metric_idx = 1:(num_metrics * 2)
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
            
                % Repeats the BCa calculation 'n_trials' times to assess the stability of the CI width.
                for t = 1:cfg_ci.n_trials
                    % Complete BCa calculation for one trial
                    % Calculate bootstrap distribution of the effect size.
                    boot_stats = zeros(B_ci_current, 1);
                    for b = 1:B_ci_current
                        boot_indices = randi(s_worker, num_probanden, [num_probanden, 1]);
                        boot_x_raw = data_x_orig(boot_indices); 
                        boot_y_raw = data_y_orig(boot_indices);
                        % Pairwise exclusion within the bootstrap sample
                        valid_boot_rows = ~isnan(boot_x_raw) & ~isnan(boot_y_raw);
                        boot_x = boot_x_raw(valid_boot_rows);
                        boot_y = boot_y_raw(valid_boot_rows);
                        n_valid_boot = numel(boot_x);
                        
                        if n_valid_boot > 0
                            if is_delta
                                gt = sum(boot_x > boot_y', 'all'); lt = sum(boot_x < boot_y', 'all');
                                boot_stats(b) = (gt - lt) / (n_valid_boot^2);
                            else
                                mx = mean(boot_x); my = mean(boot_y);
                                if(mx + my) == 0, boot_stats(b) = 0;
                                else, boot_stats(b) = abs(mx - my) / abs(mean([mx, my])); end
                            end
                        else
                            boot_stats(b) = NaN;
                        end
                    end
                    boot_stats = boot_stats(~isnan(boot_stats)); % Remove NaN results
                    if isempty(boot_stats), boot_stats = 0; end % Fallback
    
                    % Calculate Jackknife statistics for the acceleration factor 'a'.
                    jack_stats = zeros(num_probanden, 1);
                    for n = 1:num_probanden
                        jack_indices = 1:num_probanden; jack_indices(n) = [];
                        jack_x_raw = data_x_orig(jack_indices); 
                        jack_y_raw = data_y_orig(jack_indices);
        
                        % Pairwise exclusion within the Jackknife sample
                        valid_jack_rows = ~isnan(jack_x_raw) & ~isnan(jack_y_raw);
                        jack_x = jack_x_raw(valid_jack_rows);
                        jack_y = jack_y_raw(valid_jack_rows);
                        n_valid_jack = numel(jack_x);
        
                        if n_valid_jack > 0
                            if is_delta
                                gt = sum(jack_x > jack_y', 'all'); lt = sum(jack_x < jack_y', 'all');
                                jack_stats(n) = (gt - lt) / (n_valid_jack^2);
                            else
                                mx = mean(jack_x); my = mean(jack_y);
                                if(mx + my) == 0, jack_stats(n) = 0;
                                else, jack_stats(n) = abs(mx - my) / abs(mean([mx, my])); end
                            end
                        else
                            jack_stats(n) = NaN;
                        end
                    end
                    jack_stats = jack_stats(~isnan(jack_stats)); % Remove NaN results
                    if isempty(jack_stats), jack_stats = 0; end % Fallback
    
                    % Calculate BCa correction factors (z0 for bias, a for acceleration/skewness).
                    % Access effect size using dynamic index
                    if is_delta, theta_hat = p_d_vals_all(k, actual_metric_idx);
                    else, theta_hat = p_rel_vals_all(k, actual_metric_idx); end
                    z0 = norminv(sum(boot_stats < theta_hat) / B_ci_current);
                    mean_jack = mean(jack_stats);
                    a_num = sum((mean_jack - jack_stats).^3);
                    a_den = 6 * (sum((mean_jack - jack_stats).^2)).^(3/2);
                    if a_den == 0, a = 0; 
                    else, a = a_num / a_den; end
                    if ~isfinite(z0), z0 = 0; end 
                    if ~isfinite(a), a = 0; end
    
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
                        stability_all_pairs(k) = 0; % Perfect stability if all widths are zero.
                    else
                        stability_all_pairs(k) = Inf; % Infinite instability if median is zero but there is spread.
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
    
        % Convergence check: Determines if the stability has plateaued.
        if use_robust_convergence_ci
            % Robust method with smoothing and a "streak" check.
            warmup_period = cfg_ci.min_steps_for_convergence_check;
            smoothing_window = cfg_ci.smoothing_window;
            streak_needed = cfg_ci.convergence_streak_needed;
            if i >= warmup_period + smoothing_window
                % Apply a moving average to reduce noise.
                smoothed_stability = movmean(overall_stability_ci(1:i), smoothing_window, 'omitnan');
                prev_stab_smooth = smoothed_stability(end-1);
                curr_stab_smooth = smoothed_stability(end);
                % Calculate relative improvement in stability.
                if isnan(curr_stab_smooth) || isnan(prev_stab_smooth), rel_imp = 1.0;
                elseif isinf(curr_stab_smooth), rel_imp = -1.0;
                elseif prev_stab_smooth == 0
                    if curr_stab_smooth == 0, rel_imp = 0; 
                    else, rel_imp = -1.0; end
                elseif isinf(prev_stab_smooth)
                    if isinf(curr_stab_smooth), rel_imp = 0; 
                    else, rel_imp = 1.0; end
                else
                    rel_imp = (prev_stab_smooth - curr_stab_smooth) / prev_stab_smooth;
                end
                % Check if the improvement is below the tolerance.
                if abs(rel_imp) < cfg_ci.convergence_tolerance
                    convergence_streak_counter_ci = convergence_streak_counter_ci + 1;
                    fprintf(['    ' lang.bca.convergence_run_info '\n'],...
                    rel_imp * 100, convergence_streak_counter_ci, streak_needed);
                else
                    % Reset streak counter if improvement is still large.
                    convergence_streak_counter_ci = 0;
                    fprintf(['    ' lang.bca.stability_change_info '\n'], rel_imp * 100);
                end
                % Check if the required number of stable runs has been achieved.
                if convergence_streak_counter_ci >= streak_needed
                    fprintf([lang.bca.convergence_reached '\n'], rel_imp * 100);
                    fprintf([lang.bca.stable_runs_info '\n'], streak_needed);
                    converged_ci = true;
                end
            end
        else
            % Simple convergence check without smoothing.
            if i >= cfg_ci.min_steps_for_convergence_check
                prev_stab = overall_stability_ci(i - 1);
                curr_stab = overall_stability_ci(i);
                
                % Calculate relative improvement, handling edge cases.
                if isnan(curr_stab) || isnan(prev_stab), rel_imp = 1.0;
                elseif isinf(curr_stab), rel_imp = -1.0;
                elseif prev_stab == 0
                    if curr_stab == 0, rel_imp = 0; 
                    else, rel_imp = -1.0; end
                elseif isinf(prev_stab)
                    if isinf(curr_stab), rel_imp = 0; 
                    else, rel_imp = 1.0; end
                else
                    rel_imp = (prev_stab - curr_stab) / prev_stab;
                end
                fprintf(['    ' lang.bca.stability_change_info '\n'], rel_imp * 100);           
                % If improvement is below tolerance, convergence is reached.
                if abs(rel_imp) < cfg_ci.convergence_tolerance
                    fprintf([lang.bca.convergence_reached '\n'], rel_imp * 100);
                    converged_ci = true;
                end
            end
        end

        if converged_ci
            break; % Break the loop if convergence is reached.
        end
    end % End of the for-loop for stability check.
    
    % Final B-value determination based on the analysis outcome.
    B_tested_vector_ci = B_vector_ci(1:final_i_ci);
    overall_stability_ci_plotted = overall_stability_ci(1:final_i_ci);
    if converged_ci
        % If converged, the last tested B-value is used.
        selected_B_ci = B_tested_vector_ci(end);
        fprintf([lang.bca.convergence_result '\n'], selected_B_ci);
    else
        % Otherwise, the "elbow" of the curve is determined as the optimal point.
        fprintf([lang.bca.elbow_analysis_info '\n']);
        % Dynamic elbow analysis
        num_curves = num_metrics * 2;
        elbow_indices = zeros(1, num_curves);
        stability_vectors = cell(1, num_curves);
        % Dynamically create stability vectors for elbow check
        for k_sv = 1:num_metrics
            stability_vectors{k_sv} = squeeze(stability_matrix_ci(1, k_sv, 1:final_i_ci)); % Delta metrics
            stability_vectors{k_sv + num_metrics} = squeeze(stability_matrix_ci(2, k_sv, 1:final_i_ci)); % RelDiff metrics
        end

        % Loop through all stability curves.
        for k_elbow = 1:num_curves
            x_values = B_tested_vector_ci(:);
            y_values = stability_vectors{k_elbow}(:); % Vectorize to column
            % Skip if there's no data or variation to analyze.
            if any(isnan(y_values)) || numel(unique(y_values)) < 2
                % If the curve is flat (e.g., stability as IQR/Median is 0), it's stable immediately.
                % Set elbow to the first index.
                elbow_indices(k_elbow) = 1; 
                continue;
            end
            % Normalize data for consistent distance calculation.
            x_norm = (x_values - min(x_values)) / (max(x_values) - min(x_values));
            y_norm = (y_values - min(y_values)) / (max(y_values) - min(y_values));
            if all(isnan(x_norm)), x_norm(:) = 0; end
            if all(isnan(y_norm)), y_norm(:) = 0; end
            % Find the point with maximum distance to the line connecting start and end points.
            line_vec = [x_norm(end) - x_norm(1); y_norm(end) - y_norm(1)];
            vec_from_first = [x_norm - x_norm(1), y_norm - y_norm(1)];
            cross_prod = vec_from_first(:, 1) * line_vec(2) - vec_from_first(:, 2) * line_vec(1);
            [~, elbow_idx] = max(abs(cross_prod) / norm(line_vec));
            if isempty(elbow_idx), elbow_idx = numel(x_values); end
            elbow_indices(k_elbow) = elbow_idx;
        end
        % The selected B is the maximum of all determined elbow points.
        selected_B_ci = max(B_tested_vector_ci(elbow_indices));
        fprintf([lang.bca.elbow_result '\n'], selected_B_ci);
    end
    % Store stability data for JSON export
    stability_data_ci = struct();
    stability_data_ci.B_vector = B_tested_vector_ci;
    stability_data_ci.global_stability = overall_stability_ci_plotted;
    stability_data_ci.detailed_stability = stability_matrix_ci(:, :, 1:final_i_ci);
    B_ci = selected_B_ci; % Assign to the output variable.

%% 2. Creates and saves graphics to visualize the convergence analysis
% Use the passed global styles.
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');

% Plot 1: Global Convergence Curve
% This plot shows the overall stability (mean across all metrics) that was used for the convergence check.
h_fig_bca_global = figure('Name', lang.plots.titles.bca_convergence_global, 'Color', styles.colors.background, 'Visible', 'off');
tcl_bca_global = tiledlayout(1, 1, 'Padding', 'compact');
ax_global = nexttile;

% Plot the raw global stability curve.
p1_global = plot(ax_global, B_tested_vector_ci, overall_stability_ci_plotted * 100, '-o', 'LineWidth', 1.5, 'MarkerSize', 5, ...
    'Color', styles.colors.blue_marker, 'MarkerFaceColor', styles.colors.blue_marker);
grid(ax_global, 'on'); 
box(ax_global, 'on'); 
hold(ax_global, 'on');
set(ax_global, 'Color', styles.colors.background, 'GridColor', styles.colors.grid_color);

% Set title for the global plot.
fig_title_str_global = sprintf([lang.plots.titles.bca_convergence_long_n_g], config.bootstrap_ci.n_trials);
title(tcl_bca_global, fig_title_str_global, 'FontSize', styles.font.title, 'FontWeight', 'bold', 'Color', styles.colors.text);

% Set axis labels and limits.
xlim(ax_global, [min(B_tested_vector_ci), max(B_tested_vector_ci) * 1.1]);
xlabel(ax_global, lang.plots.xlabels.bootstraps, 'FontSize', styles.font.label, 'Color', styles.colors.text);
ylabel(ax_global, lang.plots.ylabels.stability, 'FontSize', styles.font.label, 'Color', styles.colors.text); 
set(ax_global, 'FontSize', styles.font.tick, 'XColor', styles.colors.text, 'YColor', styles.colors.text);

% If robust convergence was used, also plot the smoothed curve (as this was used for the check).
if use_robust_convergence_ci
    smoothing_window_global = config.bootstrap_ci.smoothing_window;
    smoothed_curve_plotted_global = movmean(overall_stability_ci_plotted, smoothing_window_global, 'omitnan');
    p2_global = plot(ax_global, B_tested_vector_ci, smoothed_curve_plotted_global * 100, '-', 'LineWidth', 1.5, 'Color', [0.8500 0.3250 0.0980]);
    lgd_global = legend(ax_global, [p1_global, p2_global], {lang.plots.legend.unsmoothed, lang.plots.legend.smoothed}, 'Location', 'best', ...
        'FontSize', styles.font.small_text);
    set(lgd_global, 'Color', styles.colors.background, 'TextColor', styles.colors.text, 'EdgeColor', styles.colors.text);
    lgd_global.ItemTokenSize = [15, 18];
end

% Add a marker for the optimal B value.
selected_B_ci_idx_for_plot_global = find(B_tested_vector_ci == selected_B_ci, 1);
if ~isempty(selected_B_ci_idx_for_plot_global)
    x_pos_global = selected_B_ci;
    y_pos_global = overall_stability_ci_plotted(selected_B_ci_idx_for_plot_global) * 100;
    plot(ax_global, x_pos_global, y_pos_global, 'x', 'Color', styles.colors.red_marker, 'MarkerSize', 10, 'LineWidth', 2, 'HandleVisibility', 'off');
    current_xlim_global = xlim(ax_global);
    if x_pos_global > (current_xlim_global(1) + (current_xlim_global(2) - current_xlim_global(1)) * 0.85)
        xlim(ax_global, [current_xlim_global(1), x_pos_global * 1.2]);
    end
    text(ax_global, x_pos_global, y_pos_global, sprintf(lang.plots.misc.optimal_b_text, selected_B_ci), ...
         'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', styles.font.small_text, 'Color', styles.colors.red_marker);
end
hold(ax_global, 'off');

% Save global graphic.
[~, fName, fExt] = fileparts(lang.files.convergence_bca_global);
filename_global = fullfile(subfolder_bca_CI, [fName, '_', ts, fExt]);
exportgraphics(h_fig_bca_global, filename_global, 'Resolution', 300, 'Padding', 30);
fprintf([lang.bca.convergence_plot_saved '\n'], filename_global);


% Plot 2: Detailed Convergence Curves (per Metric) 
% This plot shows the individual stability curves for each metric (without smoothing).
h_fig_bca_detailed = figure('Name', lang.plots.titles.bca_convergence_long, 'Color', styles.colors.background, 'Visible', 'off');

% Creates the Tiled Layout and applies the title.
num_metrics_actual = numel(metric_names);
tcl_bca_detailed = tiledlayout(2, num_metrics_actual, 'TileSpacing', 'compact', 'Padding', 'compact');
fig_title_str_detailed = sprintf([lang.plots.titles.bca_convergence_long_n_d], config.bootstrap_ci.n_trials);
title(tcl_bca_detailed, fig_title_str_detailed, 'FontSize', styles.font.title, 'FontWeight', 'bold', 'Color', styles.colors.text);

% Initialization for the plotting loop.
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
        
        current_stability_curve_plotted = squeeze(stability_matrix_ci(es_type, actual_metric_idx, 1:final_i_ci));
        
        % Plot the raw stability curve.
        p1 = plot(ax, B_tested_vector_ci, current_stability_curve_plotted * 100, '-o', 'LineWidth', 1.5, 'MarkerSize', 5, ...
            'Color', styles.colors.blue_marker, 'MarkerFaceColor', styles.colors.blue_marker, ...
            'DisplayName', lang.plots.legend.unsmoothed);
        handles_to_legend(end+1) = p1;
        names_to_legend{end+1} = lang.plots.legend.unsmoothed;
            
        grid(ax, 'on'); 
        box(ax, 'on'); 
        hold(ax, 'on');
        set(ax, 'Color', styles.colors.background, 'GridColor', styles.colors.grid_color);
        
        subplot_title_str = sprintf('%s - %s\n', metric_names{actual_metric_idx}, effect_type_names{es_type});
        title(ax, subplot_title_str, 'FontSize', styles.font.label, 'Color', styles.colors.text, 'Interpreter', 'none');
        
        % Set axis labels and limits.
        xlim(ax, [min(B_tested_vector_ci), max(B_tested_vector_ci) * 1.1]);
        xlabel(ax, lang.plots.xlabels.bootstraps, 'FontSize', styles.font.label, 'Color', styles.colors.text);
        ylabel(ax, lang.plots.ylabels.stability, 'FontSize', styles.font.label, 'Color', styles.colors.text);
        set(ax, 'FontSize', styles.font.tick, 'XColor', styles.colors.text, 'YColor', styles.colors.text);
        
        % Plot the LOCAL elbow (if convergence failed)
        if ~converged_ci && ~isempty(elbow_indices)
            % Calculate the index for the 'elbow_indices' vector
            if es_type == 1
                k_elbow_idx = actual_metric_idx;
            else
                k_elbow_idx = actual_metric_idx + num_metrics_actual;
            end
            
            % Get the specific elbow index for this subplot
            local_elbow_idx = elbow_indices(k_elbow_idx);
            
            if local_elbow_idx <= numel(current_stability_curve_plotted)
                x_local = B_tested_vector_ci(local_elbow_idx);
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
        
        % Add a marker for the final, globally optimal B value (selected_B_ci).
        selected_B_ci_idx_for_plot = find(B_tested_vector_ci == selected_B_ci, 1);
        if ~isempty(selected_B_ci_idx_for_plot)
            if selected_B_ci_idx_for_plot <= numel(current_stability_curve_plotted)
                x_pos = selected_B_ci;
                y_pos = current_stability_curve_plotted(selected_B_ci_idx_for_plot) * 100;
        
                % Plot the 'x' marker, but hide it from the legend.
                plot(ax, x_pos, y_pos, 'x', 'Color', styles.colors.red_marker, 'MarkerSize', 10, 'LineWidth', 2, 'HandleVisibility', 'off');
        
                % Expand the axis margin if the marker is too close to the edge.
                current_xlim = xlim(ax);
                if x_pos > (current_xlim(1) + (current_xlim(2) - current_xlim(1)) * 0.85)
                    xlim(ax, [current_xlim(1), x_pos * 1.2]);
                end
                
                % Add text label for the final B value (this explains the 'x')
                text(ax, x_pos, y_pos, sprintf(lang.plots.misc.optimal_b_text, selected_B_ci), 'VerticalAlignment', 'bottom', ...
                    'HorizontalAlignment', 'left', 'FontSize', styles.font.small_text, 'Color', styles.colors.red_marker);
            end
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
[~, fName, fExt] = fileparts(lang.files.convergence_bca);
filename_detailed = fullfile(subfolder_bca_CI, [fName, '_', ts, fExt]);
exportgraphics(h_fig_bca_detailed, filename_detailed, 'Resolution', 300, 'Padding', 30);
fprintf([lang.bca.convergence_plot_saved '\n\n'], filename_detailed);
end

%% 3. Final calculation of BCa confidence intervals with optimal B_ci
% Initialization of the output matrices.
ci_d_all = NaN(num_pairs, 2, num_metrics);
ci_r_all = NaN(num_pairs, 2, num_metrics);
z0_d_all = NaN(num_pairs, num_metrics);
a_d_all = NaN(num_pairs, num_metrics);
z0_r_all = NaN(num_pairs, num_metrics);
a_r_all = NaN(num_pairs, num_metrics);


% Loop over the metrics for the final calculation.
for metric_idx = 1:num_metrics
    fprintf(lang.bca.calculating_final_ci, metric_names{metric_idx});
    % Temporary variables for the results of the parfor loop.
    temp_ci_d = NaN(num_pairs, 2);
    temp_ci_r = NaN(num_pairs, 2);
    temp_z0_d = NaN(num_pairs, 1);
    temp_a_d = NaN(num_pairs, 1);
    temp_z0_r = NaN(num_pairs, 1);
    temp_a_r = NaN(num_pairs, 1);
    
    % Parallel loop over all pairs for the final BCa calculation with the optimal B_ci.
    parfor k = 1:num_pairs
        % Each iteration gets its own reproducible substream.
        s_worker = s; 
        s_worker.Substream = k;

        i = p_pair_idx_all(k, 1);
        j = p_pair_idx_all(k, 2);
        data_x_orig = p_all_data{metric_idx}(:, i);
        data_y_orig = p_all_data{metric_idx}(:, j);
        
        % Generate bootstrap distribution for Delta and Relative Difference.
        boot_d = zeros(B_ci, 1);
        boot_r = zeros(B_ci, 1);
        for b = 1:B_ci
            boot_indices = randi(s_worker, num_probanden, [num_probanden, 1]);
            boot_x_raw = data_x_orig(boot_indices); 
            boot_y_raw = data_y_orig(boot_indices);
            % Pairwise exclusion for Bootstrap sample.
            valid_boot_rows = ~isnan(boot_x_raw) & ~isnan(boot_y_raw);
            boot_x = boot_x_raw(valid_boot_rows);
            boot_y = boot_y_raw(valid_boot_rows);
            n_valid_boot = numel(boot_x);
            
            if n_valid_boot > 0
                gt = sum(boot_x > boot_y', 'all'); lt = sum(boot_x < boot_y', 'all');
                boot_d(b) = (gt - lt) / (n_valid_boot^2);
                mx = mean(boot_x); my = mean(boot_y);
                if (mx + my) == 0, boot_r(b) = 0;
                else, boot_r(b) = abs(mx - my) / mean([mx, my]); end
            else
                boot_d(b) = NaN;
                boot_r(b) = NaN;
            end
        end
        boot_d = boot_d(~isnan(boot_d)); boot_r = boot_r(~isnan(boot_r));
        if isempty(boot_d), boot_d=0; end 
        if isempty(boot_r), boot_r=0; end
        
        % Generate Jackknife distribution for Delta and Relative Difference.
        jack_d = zeros(num_probanden, 1);
        jack_r = zeros(num_probanden, 1);
        for n = 1:num_probanden
            jack_indices = 1:num_probanden; jack_indices(n) = [];
            jack_x_raw = data_x_orig(jack_indices);
            jack_y_raw = data_y_orig(jack_indices);
            % Pairwise exclusion for Jackknife sample.
            valid_jack_rows = ~isnan(jack_x_raw) & ~isnan(jack_y_raw);
            jack_x = jack_x_raw(valid_jack_rows);
            jack_y = jack_y_raw(valid_jack_rows);
            n_valid_jack = numel(jack_x);
    
            if n_valid_jack > 0
                gt = sum(jack_x > jack_y', 'all'); lt = sum(jack_x < jack_y', 'all');
                jack_d(n) = (gt - lt) / (n_valid_jack^2);
                mx = mean(jack_x); my = mean(jack_y);
                if (mx + my) == 0, jack_r(n) = 0;
                else, jack_r(n) = abs(mx - my) / abs(mean([mx, my])); end
            else
                jack_d(n) = NaN;
                jack_r(n) = NaN;
            end
        end
        jack_d = jack_d(~isnan(jack_d)); jack_r = jack_r(~isnan(jack_r));
        if isempty(jack_d), jack_d=0; end 
        if isempty(jack_r), jack_r=0; end
        
        z1 = norminv(alpha_level / 2); z2 = norminv(1 - alpha_level / 2);
        
        % BCa calculation for Cliff's Delta.
        theta_hat_d = p_d_vals_all(k, metric_idx);
        z0_d = norminv(sum(boot_d < theta_hat_d) / B_ci); % Bias-correction factor.
        mean_jack_d = mean(jack_d);
        a_num_d = sum((mean_jack_d - jack_d).^3);
        a_den_d = 6 * (sum((mean_jack_d - jack_d).^2)).^(3/2);
        if a_den_d == 0, a_d = 0; 
        else, a_d = a_num_d / a_den_d; end % Acceleration factor.
        if ~isfinite(z0_d), z0_d = 0; end 
        if ~isfinite(a_d), a_d = 0; end
        % Calculate adjusted percentile indices.
        a1_d = normcdf(z0_d + (z0_d + z1) / (1 - a_d * (z0_d + z1)));
        a2_d = normcdf(z0_d + (z0_d + z2) / (1 - a_d * (z0_d + z2)));
        if isnan(a1_d), a1_d = alpha_level / 2; end 
        if isnan(a2_d), a2_d = 1 - alpha_level / 2; end
        sorted_d = sort(boot_d);
        temp_ci_d(k, :) = [sorted_d(max(1, floor(B_ci * a1_d))), sorted_d(min(B_ci, ceil(B_ci * a2_d)))];
        temp_z0_d(k) = z0_d; temp_a_d(k) = a_d;
        
        % BCa calculation for Relative Difference.
        theta_hat_r = p_rel_vals_all(k, metric_idx);
        z0_r = norminv(sum(boot_r < theta_hat_r) / B_ci); % Bias-correction factor.
        mean_jack_r = mean(jack_r);
        a_num_r = sum((mean_jack_r - jack_r).^3);
        a_den_r = 6 * (sum((mean_jack_r - jack_r).^2)).^(3/2);
        if a_den_r == 0, a_r = 0; 
        else, a_r = a_num_r / a_den_r; end % Acceleration factor.
        if ~isfinite(z0_r), z0_r = 0; end 
        if ~isfinite(a_r), a_r = 0; end
        % Calculate adjusted percentile indices.
        a1_r = normcdf(z0_r + (z0_r + z1) / (1 - a_r * (z0_r + z1)));
        a2_r = normcdf(z0_r + (z0_r + z2) / (1 - a_r * (z0_r + z2)));
        if isnan(a1_r), a1_r = alpha_level / 2; end 
        if isnan(a2_r), a2_r = 1 - alpha_level / 2; end
        sorted_r = sort(boot_r);
        temp_ci_r(k, :) = [sorted_r(max(1, floor(B_ci * a1_r))), sorted_r(min(B_ci, ceil(B_ci * a2_r)))];
        temp_z0_r(k) = z0_r; temp_a_r(k) = a_r;
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

%% 5. Creates and saves the histogram distribution of BCa correction factors
% Use the passed global styles.
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');
effect_type_names = {'Cliff''s Delta', 'Rel Diff'};

% Graphic 1: Distribution of bias correction factors (z0).
h_fig_hist_z0 = figure('Name', lang.plots.titles.bca_z0_dist_name, ...
    'Color', styles.colors.background, 'Visible', 'off');
% Layout is dynamic
tcl_z0 = tiledlayout(2, num_metrics, 'TileSpacing', 'compact', 'Padding', 'compact');
sgtitle_str_z0 = sprintf(lang.plots.titles.bca_z0_dist, num_pairs);
title(tcl_z0, sgtitle_str_z0, 'FontSize', styles.font.title, 'FontWeight', 'bold', 'Color', styles.colors.text);

for es_type = 1:2 % 1 for Delta, 2 for Rel Diff
    for metric_idx = 1:num_metrics
        ax = nexttile; 
        if es_type == 1
            data = z0_d_all(:, metric_idx);
            face_color = styles.colors.delta_face;
        else
            data = z0_r_all(:, metric_idx);
            face_color = styles.colors.rel_face;
        end
        
        hold on;
        % Handle cases with no or little data variation for robust plotting.
        if isscalar(unique(data))
            bin_center = unique(data);
            if isempty(bin_center), bin_center = 0; end % Fallback for empty data
            bin_width = 0.02; % Fixed width for consistency

            histogram(data, 'BinEdges', [bin_center - bin_width/2, bin_center + bin_width/2], 'Normalization', 'probability', 'FaceColor', ...
                face_color, 'EdgeColor', styles.colors.bar_edge);
            
            % Set wider axis limits and defined ticks.
            xlim([bin_center - bin_width*5, bin_center + bin_width*5]);
            xticks(sort([bin_center, bin_center - bin_width*2, bin_center + bin_width*2]));
        else
            [f, xi] = ksdensity(data, 'Bandwidth', 'normal-approx');
            
            % Dynamic determination of axis ticks and limits for a clean look.
            min_val = min([data(:); xi(:)]);
            max_val = max([data(:); xi(:)]);
            data_range = max_val - min_val;
            
            if data_range > 1e-6 % Only calculate if a real range exists
                num_ticks_target = 5;
                raw_step = data_range / num_ticks_target;
                power = 10^floor(log10(raw_step));
                norm_step = raw_step / power;
                if norm_step < 1.5, nice_step = 1 * power;
                elseif norm_step < 3.5, nice_step = 2 * power;
                elseif norm_step < 7, nice_step = 5 * power;
                else
                    nice_step = 10 * power; 
                end
                nice_min = floor(min_val / nice_step) * nice_step;
                nice_max = ceil(max_val / nice_step) * nice_step;
                ticks = nice_min:nice_step:nice_max;
            else
                ticks = unique(data); % Fallback for constant data
                nice_step = 0.1 * abs(ticks(1)) + 0.01;
            end
            
            if isempty(ticks), ticks = linspace(min_val, max_val, 3); nice_step = ticks(2)-ticks(1); end
            if numel(ticks) < 2, ticks = linspace(ticks(1)-nice_step, ticks(1)+nice_step, 3); end
            
            % Set the bin edges so that the bars are centered on the ticks.
            bin_edges = (ticks(1) - nice_step/2):nice_step:(ticks(end) + nice_step/2);
            
            h_hist = histogram(data, 'BinEdges', bin_edges, 'Normalization', 'probability', 'FaceColor', face_color, 'EdgeColor', styles.colors.bar_edge);
            % Overlay the kernel density estimate, scaled to the histogram height.
            if max(f) > 0, plot(xi, f * (max(h_hist.Values)/max(f)), 'Color', styles.colors.kde_line, 'LineWidth', 1.5); end
            
            xlim([bin_edges(1), bin_edges(end)]);
            xticks(ticks);
        end
        
        grid on; box on; hold off;
        set(ax, 'Color', styles.colors.background, 'GridColor', styles.colors.grid_color);
        title(sprintf('%s - %s', metric_names{metric_idx}, effect_type_names{es_type}), 'FontSize', styles.font.label, ...
            'Color', styles.colors.text, 'Interpreter', 'none');
        set(gca, 'FontSize', styles.font.tick, 'XColor', styles.colors.text, 'YColor', styles.colors.text);
        xlabel(ax, lang.plots.xlabels.bias_z0, 'FontSize', styles.font.label, 'Color', styles.colors.text);
        ylabel(lang.plots.ylabels.rel_frequency, 'FontSize', styles.font.label, 'Color', styles.colors.text);
        ylim([0 1]);
        set(gca, 'YTick', 0:0.1:1);
    end
end
% Save graphic with padding.
[~, fName, fExt] = fileparts(lang.files.dist_bca_bias_z0); 
filename_z0 = fullfile(subfolder_bca_CI, [fName, '_', ts, fExt]);
exportgraphics(h_fig_hist_z0, filename_z0, 'Resolution', 300, 'Padding', 30);
fprintf([lang.bca.z0_histogram_saved '\n'], filename_z0);

% Graphic 2: Distribution of skewness correction factors (a).
h_fig_hist_a = figure('Name', lang.plots.titles.bca_a_dist_name, ...
    'Color', styles.colors.background, 'Visible', 'off');
% Layout is dynamic
tcl_a = tiledlayout(2, num_metrics, 'TileSpacing', 'compact', 'Padding', 'compact');
sgtitle_str_a = sprintf(lang.plots.titles.bca_a_dist, num_pairs);
title(tcl_a, sgtitle_str_a, 'FontSize', styles.font.title, 'FontWeight', 'bold', 'Color', styles.colors.text);

for es_type = 1:2 % 1 for Delta, 2 for Rel Diff
    for metric_idx = 1:num_metrics
        ax = nexttile;
        
        if es_type == 1
            data = a_d_all(:, metric_idx);
            face_color = styles.colors.delta_face;
        else
            data = a_r_all(:, metric_idx);
            face_color = styles.colors.rel_face;
        end
        
        hold on;
        % Handle different data scenarios for robust plotting.
        if isscalar(unique(data))
            bin_center = unique(data);
            if isempty(bin_center), bin_center = 0; end % Fallback for empty data
            bin_width = 0.02; % Fixed width for consistency

            histogram(data, 'BinEdges', [bin_center - bin_width/2, bin_center + bin_width/2], 'Normalization', 'probability', 'FaceColor', ...
                face_color, 'EdgeColor', styles.colors.bar_edge);
            
            % Set wider axis limits and defined ticks.
            xlim([bin_center - bin_width*5, bin_center + bin_width*5]);
            xticks(sort([bin_center, bin_center - bin_width*2, bin_center + bin_width*2]));
        else
            [f, xi] = ksdensity(data, 'Bandwidth', 'normal-approx');
            
            % Dynamic determination of axis ticks and limits.
            min_val = min([data(:); xi(:)]);
            max_val = max([data(:); xi(:)]);
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
            if isempty(ticks), ticks = linspace(min_val, max_val, 3); nice_step = ticks(2)-ticks(1); end
            if numel(ticks) < 2, ticks = linspace(ticks(1)-nice_step, ticks(1)+nice_step, 3); end
            bin_edges = (ticks(1) - nice_step/2):nice_step:(ticks(end) + nice_step/2);
            
            h_hist = histogram(data, 'BinEdges', bin_edges, 'Normalization', 'probability', 'FaceColor', face_color, 'EdgeColor', styles.colors.bar_edge);
            if max(f) > 0, plot(xi, f * (max(h_hist.Values)/max(f)), 'Color', styles.colors.kde_line, 'LineWidth', 1.5); end
            
            xlim([bin_edges(1), bin_edges(end)]);
            xticks(ticks);
        end
        
        grid on; box on; hold off;
        set(ax, 'Color', styles.colors.background, 'GridColor', styles.colors.grid_color);
        title(sprintf('%s - %s', metric_names{metric_idx}, effect_type_names{es_type}), 'FontSize', styles.font.label, ...
            'Color', styles.colors.text, 'Interpreter', 'none');
        set(gca, 'FontSize', styles.font.tick, 'XColor', styles.colors.text, 'YColor', styles.colors.text);
        xlabel(lang.plots.xlabels.skewness_a, 'FontSize', styles.font.label, 'Color', styles.colors.text); 
        ylabel(lang.plots.ylabels.rel_frequency, 'FontSize', styles.font.label, 'Color', styles.colors.text);
        ylim([0 1]);
        set(gca, 'YTick', 0:0.1:1);
    end
end
% Save graphic with padding.
[~, fName, fExt] = fileparts(lang.files.dist_bca_skew_a); 
filename_a = fullfile(subfolder_bca_CI, [fName, '_', ts, fExt]);
exportgraphics(h_fig_hist_a, filename_a, 'Resolution', 300, 'Padding', 30);
fprintf([lang.bca.a_histogram_saved '\n'], filename_a);

%% 6. Creates and saves the histogram distribution of CI widths
% Calculation of CI widths from the lower and upper bounds.
ci_widths_d = ci_d_all(:, 2, :) - ci_d_all(:, 1, :);
ci_widths_r = ci_r_all(:, 2, :) - ci_r_all(:, 1, :);

% Graphic 3: Distribution of CI widths.
h_fig_hist_widths = figure('Name', lang.plots.titles.ci_width_dist_name, ...
    'Color', styles.colors.background, 'Visible', 'off');
% Layout is dynamic
tcl_widths = tiledlayout(2, num_metrics, 'TileSpacing', 'compact', 'Padding', 'compact');
sgtitle_str = sprintf(lang.plots.titles.ci_width_dist, num_pairs, B_ci);
title(tcl_widths, sgtitle_str, 'FontSize', styles.font.title, 'FontWeight', 'bold', 'Color', styles.colors.text);

for es_type = 1:2 % 1 for Delta, 2 for Rel Diff
    for metric_idx = 1:num_metrics
        ax = nexttile; 
        if es_type == 1
            data = ci_widths_d(:, 1, metric_idx);
            face_color = styles.colors.delta_face;
        else
            data = ci_widths_r(:, 1, metric_idx);
            face_color = styles.colors.rel_face;
        end
        hold on;
        % Handle different data scenarios for robust plotting.
        if isscalar(unique(data))
            bin_center = unique(data);
            if isempty(bin_center), bin_center = 0; end % Fallback for empty data
            bin_width = 0.02; % Fixed width for consistency

            histogram(data, 'BinEdges', [bin_center - bin_width/2, bin_center + bin_width/2], 'Normalization', 'probability', 'FaceColor', ...
                face_color, 'EdgeColor', styles.colors.bar_edge);
            
            % Set wider axis limits and defined ticks.
            xlim([bin_center - bin_width*5, bin_center + bin_width*5]);
            xticks(sort([bin_center, bin_center - bin_width*2, bin_center + bin_width*2]));
        else
            [f, xi] = ksdensity(data, 'Bandwidth', 'normal-approx');
            
            % Dynamic determination of axis ticks and limits.
            min_val = min([data(:); xi(:)]);
            max_val = max([data(:); xi(:)]);
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
            
            if isempty(ticks), ticks = linspace(min_val, max_val, 3); nice_step = ticks(2)-ticks(1); end
            if numel(ticks) < 2, ticks = linspace(ticks(1)-nice_step, ticks(1)+nice_step, 3); end
            
            bin_edges = (ticks(1) - nice_step/2):nice_step:(ticks(end) + nice_step/2);
            
            h_hist = histogram(data, 'BinEdges', bin_edges, 'Normalization', 'probability', 'FaceColor', face_color, 'EdgeColor', styles.colors.bar_edge);
            if max(f) > 0, plot(xi, f * (max(h_hist.Values)/max(f)), 'Color', styles.colors.kde_line, 'LineWidth', 1.5); end
            
            xlim([bin_edges(1), bin_edges(end)]);
            xticks(ticks);
        end
        
        grid on; box on; hold off;
        set(ax, 'Color', styles.colors.background, 'GridColor', styles.colors.grid_color);
        title(sprintf('%s - %s', metric_names{metric_idx}, effect_type_names{es_type}), 'FontSize', styles.font.label, ...
            'Color', styles.colors.text, 'Interpreter', 'none');
        set(gca, 'FontSize', styles.font.tick, 'XColor', styles.colors.text, 'YColor', styles.colors.text);
        xlabel(lang.plots.xlabels.ci_width, 'FontSize', styles.font.label, 'Color', styles.colors.text); 
        ylabel(lang.plots.ylabels.rel_frequency, 'FontSize', styles.font.label, 'Color', styles.colors.text); 
        ylim([0 1]);
        set(gca, 'YTick', 0:0.1:1);
    end
end

% Save graphic with padding.
[~, fName, fExt] = fileparts(lang.files.dist_ci_widths); 
filename_widths = fullfile(subfolder_bca_CI, [fName, '_', ts, fExt]);
exportgraphics(h_fig_hist_widths, filename_widths, 'Resolution', 300, 'Padding', 30);
fprintf([lang.bca.ci_width_histogram_saved '\n'], filename_widths);
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
