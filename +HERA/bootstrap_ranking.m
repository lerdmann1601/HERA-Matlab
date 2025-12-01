function [final_bootstrap_ranks, selected_B_final, stability_data_rank, h_figs_rank, h_fig_hist_rank] = bootstrap_ranking(all_data, ...
          thresholds, config, dataset_names, final_rank, pair_idx_all, num_probanden, graphics_dir, csv_dir, manual_B, s, styles, lang, base_name)
% BOOTSTRAP_RANKING - Performs a cluster bootstrap analysis to determine rank confidence intervals.
%
% Syntax:
%   [final_bootstrap_ranks, selected_B_final, stability_data_rank, h_figs_rank, h_fig_hist_rank] = bootstrap_ranking(all_data, ...
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
%       Analyzes the stability of the rank confidence intervals over multiple trials and an increasing number of B-values. 
%       Stability is quantified using the ratio of the Interquartile Range (IQR) to the median of the confidence interval widths.
%   2.  Selection of the final B-value ('selected_B_final'): 
%       The optimal B is determined either when stability converges (i.e., the relative improvement falls below a tolerance) - 
%       or, if convergence is not reached, through an "elbow analysis" of the stability curve.
%   3.  Visualization (Convergence): 
%       Generates and saves the convergence plot showing the stability analysis.
%   4.  Final bootstrap analysis: 
%       Using the optimal B-value, the function runs a final, comprehensive bootstrap analysis to generate the definitive rank distribution for each dataset.
%   5.  Result Output (Console): 
%       Calculates the frequency of each rank for each dataset and prints a formatted table to the console.
%   6.  Result Output (CSV): 
%       Saves the detailed rank distribution (counts and percentages) for each dataset to a CSV file.
%   7.  Visualization (Distribution): 
%       Generates and saves histograms displaying the final rank distribution for each dataset.
%
% Inputs:
%   all_data        - Cell array containing the original data matrices for each metric (1, 2, or 3 cells).
%   thresholds      - Struct with the pre-calculated effect size thresholds (thresholds.d_thresh, thresholds.rel_thresh).
%   config          - Global configuration struct (specifically uses 'config.bootstrap_ranks' and 'config.ranking_mode').
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
%
% Author:   Lukas von Erdmannsdorff
% Date:     12.10.2025
% Version:  1


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
    % Check if the robust convergence check (with smoothing) should be used based on the configuration.
    use_robust_convergence = isfield(cfg_rank, 'smoothing_window') && ~isempty(cfg_rank.smoothing_window) ...
                            && isfield(cfg_rank, 'convergence_streak_needed') && ~isempty(cfg_rank.convergence_streak_needed);
    
    % Inform the user about the start of the search and the criteria being used.
    fprintf(['\n' lang.ranking.searching_optimal_b '\n']);
    if use_robust_convergence
        fprintf([lang.ranking.primary_criterion '\n'])
        fprintf([lang.ranking.robust_convergence_info '\n'], ...
            cfg_rank.smoothing_window, ...
            cfg_rank.convergence_streak_needed, ...
            cfg_rank.convergence_tolerance * 100, ...
            cfg_rank.min_steps_for_convergence_check, ...
            cfg_rank.B_end);
    else
        fprintf([lang.ranking.simple_convergence_info '\n'], ...
            cfg_rank.convergence_tolerance * 100, ...
            cfg_rank.min_steps_for_convergence_check, ...
            cfg_rank.B_end);
    end
    fprintf([lang.ranking.secondary_criterion '\n']);
    
    % Initialize vectors for the stability analysis.
    B_vector_b = cfg_rank.B_start:cfg_rank.B_step:cfg_rank.B_end; 
    stability_vector_b = NaN(1, numel(B_vector_b)); 
    final_b_idx = 0; 
    converged = false; 
    % Initialization for the robust convergence check.
    convergence_streak_counter = 0;
    
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
            
            % Inner loop: Performs the bootstrap process Br_b times to generate one rank distribution for one trial.
            for bb_b = 1:Br_b
                % 1. Perform a cluster bootstrap: Subjects (clusters) are drawn with replacement.
                % The indices represent the rows (subjects) to be included in the bootstrap sample.
                boot_indices = randi(s_worker, n_subj_b, [n_subj_b, 1]);
                
                % 2. A new, bootstrapped dataset is created by sampling rows from the original data based on the drawn indices.
                % Dynamically create bootstrap_all_data cell array
                bootstrap_all_data = cell(1, num_metrics);
                for m_idx = 1:num_metrics
                    bootstrap_all_data{m_idx} = all_data{m_idx}(boot_indices, :);
                end
                
                % 3. Effect sizes (Cliff's Delta, Relative Difference) are completely recalculated for the new bootstrap sample.
                num_pairs = size(pair_idx_all, 1);
                % Dynamically size effect size matrices
                bootstrap_d_vals_all = zeros(num_pairs, num_metrics);
                bootstrap_rel_vals_all = zeros(num_pairs, num_metrics);
                
                % Loop through all pairwise comparisons to calculate effect sizes.
                % Loop from 1 to num_metrics
                for p_idx = 1:num_pairs
                    i = pair_idx_all(p_idx, 1); % Index of the first dataset in the pair.
                    j = pair_idx_all(p_idx, 2); % Index of the second dataset in the pair.
                    for metric_idx = 1:num_metrics
                        boot_data_metric = bootstrap_all_data{metric_idx};
                        % Find valid (non-NaN) rows within the bootstrap sample for this specific pair.
                        valid_boot_rows = ~isnan(boot_data_metric(:, i)) & ~isnan(boot_data_metric(:, j));
                        x = boot_data_metric(valid_boot_rows, i);
                        y = boot_data_metric(valid_boot_rows, j);
                        n_valid = size(x, 1);
        
                        if n_valid > 0
                            % Calculate Cliff's Delta using the dynamic sample size (n_valid).
                            gt = sum(x > y', 'all');
                            lt = sum(x < y', 'all');
                            bootstrap_d_vals_all(p_idx, metric_idx) = (gt - lt) / (n_valid^2);
                            
                            % Calculate the Relative Mean Difference.
                            mx = mean(x);
                            my = mean(y);
                            if (mx + my) == 0
                                rel_diff = 0; % Avoid division by zero.
                            else
                                rel_diff = abs(mx - my) / abs(mean([mx, my]));
                            end
                            bootstrap_rel_vals_all(p_idx, metric_idx) = rel_diff;
                        else
                            % If no valid data pairs exist, store NaN.
                            bootstrap_d_vals_all(p_idx, metric_idx) = NaN;
                            bootstrap_rel_vals_all(p_idx, metric_idx) = NaN;
                        end
                    end
                end
                bootstrap_effect_sizes = struct('d_vals_all', bootstrap_d_vals_all, 'rel_vals_all', bootstrap_rel_vals_all);
                
                % 4. The complete ranking algorithm is run on the bootstrapped effect sizes.
                % The 'config' struct, which contains the 'ranking_mode', is passed here.
                % 'calculate_ranking' will use this to apply the correct logic (e.g., M1, M1_M2, M1_M3A, etc.)
                [~, bootstrap_rank] = HERA.calculate_ranking(...
                    bootstrap_all_data, bootstrap_effect_sizes, thresholds, config, dataset_names, pair_idx_all);
                
                % 5. The resulting rank vector for this single bootstrap iteration is saved.
                rank_tmp_b(:, bb_b) = bootstrap_rank;
            end
            
            % After Br_b iterations, calculate the width of the 95% confidence interval of the ranks for this one trial.
            ci_bounds_b = quantile(rank_tmp_b, [0.025, 0.975], 2);
            ci_widths_b(t_b, :) = ci_bounds_b(:, 2) - ci_bounds_b(:, 1);
        end
        
        % Calculate the stability metric for the current B-value. This is the median of the relative IQR of CI widths across all datasets.
        med_widths = median(ci_widths_b, 1);
        iqr_widths = iqr(ci_widths_b, 1);
        relative_iqr = zeros(size(med_widths));
        non_zero_median_idx = med_widths > 0; % Avoid division by zero.
        relative_iqr(non_zero_median_idx) = iqr_widths(non_zero_median_idx) ./ med_widths(non_zero_median_idx);
        
        % The final stability value is the median of these relative IQRs.
        stability_vector_b(b_idx) = median(relative_iqr, 'omitnan');
        final_b_idx = b_idx; % Keep track of the last executed index.
        
        % Perform the convergence check to see if the stability has plateaued.
        if use_robust_convergence
            % This method uses a smoothed stability curve and requires multiple consecutive stable runs.
            warmup_period = cfg_rank.min_steps_for_convergence_check;
            smoothing_window = cfg_rank.smoothing_window;
            streak_needed = cfg_rank.convergence_streak_needed;
            
            % Start the check only after the warmup and smoothing window periods have passed.
            if b_idx >= warmup_period + smoothing_window
                % Apply a moving average to the stability curve to reduce noise.
                smoothed_stability = movmean(stability_vector_b(1:b_idx), smoothing_window, 'omitnan');
                prev_stab_smooth = smoothed_stability(end - 1);
                curr_stab_smooth = smoothed_stability(end);
                
                % Calculate the relative improvement based on the smoothed values, handling edge cases.
                if isnan(curr_stab_smooth) || isnan(prev_stab_smooth), relative_improvement=1.0;
                elseif isinf(curr_stab_smooth), relative_improvement=-1.0;
                elseif prev_stab_smooth==0
                    if curr_stab_smooth==0, relative_improvement=0;
                    else, relative_improvement=-1.0; end
                elseif isinf(prev_stab_smooth)
                    if isinf(curr_stab_smooth), relative_improvement=0;
                    else, relative_improvement=1.0; end
                else
                    relative_improvement=(prev_stab_smooth-curr_stab_smooth)/prev_stab_smooth;
                end
                
                % Check if the improvement is below the tolerance threshold.
                if abs(relative_improvement) < cfg_rank.convergence_tolerance
                    convergence_streak_counter = convergence_streak_counter + 1;
                    fprintf(['    ' lang.ranking.convergence_run_info '\n'], relative_improvement * 100, convergence_streak_counter, streak_needed);
                else
                    % If the change is still large, reset the streak counter.
                    convergence_streak_counter = 0;
                    fprintf(['    ' lang.ranking.stability_change_info '\n'], relative_improvement * 100);
                end
                
                % If the number of consecutive stable runs meets the required streak, convergence is reached.
                if convergence_streak_counter >= streak_needed
                    fprintf([lang.ranking.convergence_reached '\n'], relative_improvement * 100);
                    fprintf([lang.ranking.stable_runs_info '\n'], streak_needed);
                    converged = true;
                end
            end
        else
            % This is the simple convergence check, which is faster but more sensitive to noise.
            if b_idx >= cfg_rank.min_steps_for_convergence_check
                prev_stability = stability_vector_b(b_idx - 1);
                curr_stability = stability_vector_b(b_idx);
                
                % Calculate relative improvement, handling edge cases as in the robust method.
                if isnan(curr_stability) || isnan(prev_stability) 
                    relative_improvement = 1.0;
                elseif isinf(curr_stability)
                    relative_improvement = -1.0;
                elseif prev_stability == 0
                    if curr_stability == 0, relative_improvement = 0;
                    else, relative_improvement = -1.0; end
                elseif isinf(prev_stability)
                    if isinf(curr_stability), relative_improvement = 0;
                    else, relative_improvement = 1.0; end
                else
                    relative_improvement = (prev_stability - curr_stability) / prev_stability; 
                end
                fprintf(['    ' lang.ranking.stability_change_info '\n'], relative_improvement * 100);              
                % If the improvement is below the tolerance, convergence is reached immediately.
                if abs(relative_improvement) < cfg_rank.convergence_tolerance
                    fprintf([lang.ranking.convergence_reached '\n'], relative_improvement * 100);
                    converged = true;
                end
            end
        end
        
        if converged
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
        x_values_b = B_tested_vector_b(:);
        y_values_b = stability_vector_b_plotted(:);
        
        % Handle cases with no data or variation to avoid errors.
        if any(isnan(y_values_b)) || numel(unique(y_values_b)) < 2
            fprintf([' ' lang.ranking.stability_vector_warning '\n']);
            % If the curve is flat (e.g., stability as IQR/Median is 0), it's stable immediately.
            % Set elbow to the first index.
            elbow_idx_rank = 1; 
        else
            % Find the point on the curve with the greatest perpendicular distance to the line connecting the start and end points.
            % Normalize the data to a [0, 1] range for consistent distance calculation.
            x_norm_b = (x_values_b - min(x_values_b)) / (max(x_values_b) - min(x_values_b));
            y_norm_b = (y_values_b - min(y_values_b)) / (max(y_values_b) - min(y_values_b));
            % Handle cases where normalization might result in NaN (e.g., if max equals min).
            if all(isnan(x_norm_b)), x_norm_b(:) = 0; end
            if all(isnan(y_norm_b)), y_norm_b(:) = 0; end
            
            % Vector representing the line from the first to the last point.
            line_vec_b = [x_norm_b(end) - x_norm_b(1); y_norm_b(end) - y_norm_b(1)];
            % Vectors from the first point to all other points.
            vec_from_first_b = [x_norm_b - x_norm_b(1), y_norm_b - y_norm_b(1)];
            % Calculate the perpendicular distance using the cross product formula.
            cross_prod_b = vec_from_first_b(:, 1) * line_vec_b(2) - vec_from_first_b(:, 2) * line_vec_b(1);
            distances_b = abs(cross_prod_b) / norm(line_vec_b);
            [~, elbow_idx_rank] = max(distances_b); % The elbow is the point of maximum distance.
        end
        
        if isempty(elbow_idx_rank), elbow_idx_rank = numel(x_values_b); end % Fallback.
        selected_B_final = B_tested_vector_b(elbow_idx_rank);
        fprintf([lang.ranking.elbow_result '\n'], selected_B_final);
    end
    % Store stability data for JSON export
    stability_data_rank = struct();
    stability_data_rank.B_vector = B_tested_vector_b;
    stability_data_rank.global_stability = stability_vector_b_plotted;
    % Rank stability only has one (global) curve, so detailed_stability is empty.
    stability_data_rank.detailed_stability = [];
    
%% 2. Create and save the convergence graphic.
    % Set global default font properties for the plot.
    set(0, 'DefaultAxesFontName', 'Arial');
    set(0, 'DefaultTextFontName', 'Arial');
    
    % Create the figure with specified properties.
    h_fig_rank = figure('Name', lang.plots.titles.rank_stability_convergence_name, 'Color', styles.colors.background, 'Visible', 'off');
    tcl_rank = tiledlayout(1, 1, 'Padding', 'compact');
    ax = nexttile;
    
    % Prepare plot handles and names for the dynamic legend 
    handles_to_legend = [];
    names_to_legend = {};
    
    % Plot the raw stability curve.
    p1 = plot(ax, B_tested_vector_b, stability_vector_b_plotted * 100, '-o', ...
         'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', styles.colors.blue_marker, 'Color', styles.colors.blue_marker);
    handles_to_legend(end+1) = p1;
    names_to_legend{end+1} = lang.plots.legend.unsmoothed;
         
    grid(ax, 'on'); box(ax, 'on'); hold(ax, 'on'); 
    set(ax, 'Color', styles.colors.background, 'GridColor', styles.colors.grid_color);
    
    % If robust convergence was used, also plot the smoothed curve.
    if use_robust_convergence
        smoothing_window = config.bootstrap_ranks.smoothing_window;
        smoothed_stability_plotted = movmean(stability_vector_b_plotted, smoothing_window, 'omitnan');
        p2 = plot(ax, B_tested_vector_b, smoothed_stability_plotted * 100, '-', 'LineWidth', 1.5, 'Color', [0.8500 0.3250 0.0980]);
        handles_to_legend(end+1) = p2;
        names_to_legend{end+1} = lang.plots.legend.smoothed;
    end
    
    % Set titles and labels for the plot.
    fig_title_str = sprintf(lang.plots.titles.rank_stability_convergence_n, config.bootstrap_ranks.n_trials);
    title(tcl_rank, fig_title_str, 'Color', styles.colors.text, 'FontSize', styles.font.title, 'FontWeight', 'bold');
    
    xlim(ax, [min(B_tested_vector_b), max(B_tested_vector_b) * 1.2]);
    xlabel(ax, lang.plots.xlabels.bootstraps, 'Color', styles.colors.text, 'FontSize', styles.font.label);
    ylabel(ax, lang.plots.ylabels.stability, 'Color', styles.colors.text, 'FontSize', styles.font.label); 
    set(ax, 'FontSize', styles.font.tick, 'XColor', styles.colors.text, 'YColor', styles.colors.text);

    % Plot the LOCAL elbow (if convergence failed)
    if ~converged
        % Get the specific elbow index for this plot
        local_elbow_idx = elbow_idx_rank; % From the elbow analysis above
        
        if local_elbow_idx <= numel(stability_vector_b_plotted)
            x_local = B_tested_vector_b(local_elbow_idx);
            y_local = stability_vector_b_plotted(local_elbow_idx) * 100;
            
            % Define the color for the local elbow (same as smoothed curve)
            color_local_elbow = [0.8500 0.3250 0.0980]; 
    
             % Plot as a filled circle with a dashed line style for legend
            p_elbow = plot(ax, x_local, y_local, ':o', ... 
                'MarkerFaceColor', color_local_elbow, ...
                'MarkerEdgeColor', color_local_elbow, ... 
                'MarkerSize', 7, 'HandleVisibility', 'on', ... 
                'LineWidth', 1.5, ... 
                'DisplayName', lang.plots.legend.local_elbow); 
            
            % Add vertical line for the elbow 
            xline(ax, x_local, '--', 'Color', color_local_elbow, 'LineWidth', 1.0, 'HandleVisibility', 'off');
            
            % Add the handle to our dynamic legend
            handles_to_legend(end+1) = p_elbow;
            names_to_legend{end+1} = lang.plots.legend.local_elbow; 
        end
    end

    % Add a marker and text to indicate the finally selected optimal B-value.
    selected_B_final_idx = find(B_tested_vector_b == selected_B_final, 1);
    if ~isempty(selected_B_final_idx)
        x_pos = selected_B_final;
        y_pos = stability_vector_b_plotted(selected_B_final_idx) * 100;
        
        % Plot the 'x' marker, but hide it from the legend.
        plot(ax, x_pos, y_pos, 'x', 'Color', styles.colors.red_marker, ...
             'MarkerSize', 10, 'LineWidth', 2, 'HandleVisibility', 'off');
    
        % Extend the x-axis margin if the marker is too close to the edge.
        current_xlim = xlim(ax);
        if x_pos > (current_xlim(1) + (current_xlim(2) - current_xlim(1)) * 0.85)
            xlim(ax, [current_xlim(1), x_pos * 1.2]);
        end
    
        % Add text label for the final B value (this explains the 'x')
        text(ax, x_pos, y_pos, sprintf(lang.plots.misc.optimal_b, selected_B_final), 'VerticalAlignment', 'bottom', ...
             'HorizontalAlignment', 'left', 'FontSize', styles.font.tick, 'Color', styles.colors.red_marker);
    end
    
    % Create the final, dynamic legend 
    lgd = legend(handles_to_legend, names_to_legend, 'Location', 'best');
    set(lgd, 'Color', styles.colors.background, 'TextColor', styles.colors.text, 'EdgeColor', styles.colors.text);
    lgd.FontSize = styles.font.tick;
    lgd.ItemTokenSize = [15, 18];
    
    hold(ax, 'off');
    
    % Save graphic.
    [~, fName, fExt] = fileparts(lang.files.convergence_rank_stability);
    filename = fullfile(subfolder_ranking, [fName, '_', ts, fExt]);
    exportgraphics(h_fig_rank, filename, 'Resolution', 300, 'Padding', 30);
    fprintf([lang.ranking.convergence_plot_saved '\n'], filename);
   
    % Add the figure handle to the output array.
    if isgraphics(h_fig_rank)
        h_figs_rank(end+1) = h_fig_rank;
    end
end

%% 3. Final Bootstrap Analysis with the Optimal Number of Repetitions
% This section performs the main bootstrap process with the final, optimally determined B-value.
final_bootstrap_ranks = zeros(num_datasets_b, selected_B_final);
parfor bb_b = 1:selected_B_final
    % This is the same bootstrap process as in the stability check, but run 'selected_B_final' times to generate the final rank distribution.
    % Each parallel worker gets its own reproducible substream.
    s_worker = s; % Workaround for passing the stream object to parfor.
    s_worker.Substream = bb_b;
    
    % 1. Perform a cluster bootstrap: Draw subjects with replacement.
    boot_indices = randi(s_worker, n_subj_b, [n_subj_b, 1]);
    
    % 2. Create the bootstrapped dataset from the sampled subjects.
    % Dynamically create bootstrap_all_data cell array
    bootstrap_all_data = cell(1, num_metrics);
    for m_idx = 1:num_metrics
        bootstrap_all_data{m_idx} = all_data{m_idx}(boot_indices, :);
    end
    
    % 3. Recalculate the effect sizes for the new bootstrapped dataset.
    num_pairs = size(pair_idx_all, 1);
    % Dynamically size effect size matrices
    bootstrap_d_vals_all = zeros(num_pairs, num_metrics);
    bootstrap_rel_vals_all = zeros(num_pairs, num_metrics);
    
    % Loop from 1 to num_metrics
    for p_idx = 1:num_pairs
        i = pair_idx_all(p_idx, 1);
        j = pair_idx_all(p_idx, 2);
        for metric_idx = 1:num_metrics
            boot_data_metric = bootstrap_all_data{metric_idx};   
            % Find valid (non-NaN) rows within the bootstrap sample for this pair.
            valid_boot_rows = ~isnan(boot_data_metric(:, i)) & ~isnan(boot_data_metric(:, j));
            x = boot_data_metric(valid_boot_rows, i);
            y = boot_data_metric(valid_boot_rows, j);
            n_valid = size(x, 1);

            if n_valid > 0
                % Calculate Cliff's Delta with the dynamic sample size.
                gt = sum(x > y', 'all');
                lt = sum(x < y', 'all');
                bootstrap_d_vals_all(p_idx, metric_idx) = (gt - lt) / (n_valid^2);
                % Calculate the Relative Mean Difference.
                mx = mean(x);
                my = mean(y);
                if (mx + my) == 0
                    rel_diff = 0;
                else
                    rel_diff = abs(mx - my) / abs(mean([mx, my]));
                end
                bootstrap_rel_vals_all(p_idx, metric_idx) = (isnan(rel_diff)) * 0 + (~isnan(rel_diff)) * rel_diff;
            else
                bootstrap_d_vals_all(p_idx, metric_idx) = NaN;
                bootstrap_rel_vals_all(p_idx, metric_idx) = NaN;
            end
        end
    end
    bootstrap_effect_sizes = struct('d_vals_all', bootstrap_d_vals_all, 'rel_vals_all', bootstrap_rel_vals_all);
    
    % 4. Call the ranking algorithm with the bootstrapped data.
    % The 'config' struct (containing 'ranking_mode') is passed.
    [~, bootstrap_rank] = HERA.calculate_ranking(...
        bootstrap_all_data, bootstrap_effect_sizes, thresholds, config, dataset_names, pair_idx_all);
    
    % 5. Save the final rank vector of this single bootstrap iteration.
    final_bootstrap_ranks(:, bb_b) = bootstrap_rank;
end

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
fprintf(['\n' lang.bca.saving_csv '\n']);
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
% Set global default font properties for the plot.
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');
% Create the figure for the histograms.
num_datasets_plot = size(final_bootstrap_ranks, 1);
h_fig_hist_rank = figure('Name', lang.plots.titles.rank_dist_name, 'Color', styles.colors.background, 'Visible', 'off');
% Set up the title and layout for the combined plot.
tcl_hist = tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');
sgtitle_str = sprintf(lang.plots.titles.rank_dist, selected_B_final);
title(tcl_hist, sgtitle_str, 'FontSize', styles.font.title, 'FontWeight', 'bold', 'Color', styles.colors.text);

% Get the sorting order based on the final rank
[~, sort_idx] = sort(final_rank);

% Create a separate subplot with a rank histogram for each dataset.
% Loop over the sorted indices 
for plot_idx = 1:num_datasets_plot
    % Get the actual dataset index based on the sorted rank
    i = sort_idx(plot_idx); 
    
    ax = nexttile;
    ranks_data = final_bootstrap_ranks(i, :);
    
    % Handle the special case where only a single, unique rank value is present for a dataset.
    if isscalar(unique(ranks_data))
        % This prevents histogram errors and creates a single, narrow bar.
        val = ranks_data(1);
        bin_width = 0.25;
        edges = [val - bin_width/2, val + bin_width/2];
        histogram(ranks_data, 'BinEdges', edges, 'Normalization', 'probability', ...
            'FaceColor', styles.colors.bar_face, 'EdgeColor', styles.colors.bar_edge);
    else
        % Normal case for integer-valued ranks: use a bin width of 1.
        h = histogram(ranks_data, 'Normalization', 'probability', 'FaceColor', styles.colors.bar_face, 'EdgeColor', styles.colors.bar_edge);
        h.BinWidth = 1; 
        h.BinLimits = [min(ranks_data)-0.5, max(ranks_data)+0.5];
    end
    grid on; box on;
    set(ax, 'Color', styles.colors.background, 'GridColor', styles.colors.grid_color);
    
    % Set titles and labels for each subplot.
    title(dataset_names{i}, 'FontSize', styles.font.label, 'Color', styles.colors.text, 'Interpreter', 'none');
    xlabel(lang.plots.xlabels.rank, 'FontSize', styles.font.label, 'Color', styles.colors.text);
    ylabel(lang.plots.ylabels.rel_frequency, 'FontSize', styles.font.label, 'Color', styles.colors.text);

    % Set general axis properties.
    set(gca, 'FontSize', styles.font.tick, ...
             'XTick', unique(ranks_data), ...
             'XColor', styles.colors.text, ...
             'YColor', styles.colors.text);
    xlim([min(ranks_data)-1, max(ranks_data)+1]);
    % Enforce the y-axis to always span from 0 to 1.
    ylim([0 1]);    
    % Force the y-axis ticks for readable standard.
    set(gca, 'YTick', 0:0.1:1);
end
% Save the complete histogram figure to a file.
[~, fName, fExt] = fileparts(lang.files.dist_bootstrap_ranks);
filename = fullfile(subfolder_ranking, [fName, '_', ts, fExt]);
exportgraphics(h_fig_hist_rank, filename, 'Resolution', 300, 'Padding', 30);
fprintf([lang.ranking.histogram_plot_saved '\n'], filename);
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