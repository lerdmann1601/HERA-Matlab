function [h_fig_global, h_fig_detailed] = threshold_convergence(B_vector, stability_data, selected_B, config, styles, lang, graphics_dir)
% THRESHOLD_CONVERGENCE - Generates the Global and Detailed convergence plots for thresholds.
%
% Syntax:
%   [h_fig_global, h_fig_detailed] = threshold_convergence(B_vector, stability_data, selected_B, config, styles, lang, graphics_dir)
%
% Description:
%   This function creates two figures to visualize the convergence of bootstrap thresholds:
%   1. Global Convergence: Shows the overall stability (mean across all metrics) vs. B.
%      Includes smoothing if robust convergence was used.
%   2. Detailed Convergence: Shows stability curves for each metric individually.
%      Displays local "elbow" points if global convergence was not reached.
%   
%   Both plots are saved to the 'Threshold_Analysis' subdirectory within graphics_dir.
%
% Inputs:
%   B_vector       - Vector of tested bootstrap repetition counts (B).
%   stability_data - Struct containing stability analysis results:
%                    * .global_stability: Vector of overall stability values.
%                    * .detailed_stability: Matrix [2 x n_metrics x n_steps] of stability values.
%                    * .elbow_indices: (Optional) Indices of local elbows for detailed plot.
%                    * .converged: (Optional) Boolean, true if global convergence reached.
%   selected_B     - The final selected number of bootstrap samples.
%   config         - Configuration struct (specifically config.bootstrap_thresholds).
%   styles         - Struct with plot styling information (colors, fonts).
%   lang           - Language package loaded from a JSON file.
%   graphics_dir   - Directory where plots should be saved.
%
% Outputs:
%   h_fig_global   - Handle to the global convergence figure.
%   h_fig_detailed - Handle to the detailed convergence figure.
%
% Author: Lukas von Erdmannsdorff

    % Create a dedicated subfolder for the threshold analysis plots.
    subfolder_thresholds = fullfile(graphics_dir, 'Threshold_Analysis');
    if ~exist(subfolder_thresholds, 'dir')
        mkdir(subfolder_thresholds);
    end
    
    ts = config.timestamp;
    metric_names = config.metric_names;
    
    % Use the passed global styles.
    set(0, 'DefaultAxesFontName', 'Arial');
    set(0, 'DefaultTextFontName', 'Arial');

    % Extract data from stability_data
    overall_stability_thr_plotted = stability_data.global_stability;
    
    % Determine if robust convergence was used based on config
    cfg_thr = config.bootstrap_thresholds;
    use_robust_convergence_thr = isfield(cfg_thr, 'smoothing_window') && ~isempty(cfg_thr.smoothing_window) ...
                               && isfield(cfg_thr, 'convergence_streak_needed') && ~isempty(cfg_thr.convergence_streak_needed);

    %% Plot 1: Global Convergence Curve
    % This plot shows the overall stability (mean across all metrics) that was used for the convergence check.
    h_fig_global = figure('Name', lang.plots.titles.threshold_convergence_global, 'Color', styles.colors.background, 'Visible', 'off');
    tcl_thr_global = tiledlayout(1, 1, 'Padding', 'compact');
    ax_global = nexttile;

    % Plot the raw global stability curve.
    p1_global = plot(ax_global, B_vector, overall_stability_thr_plotted * 100, '-o', 'LineWidth', 1.5, 'MarkerSize', 5, ...
        'Color', styles.colors.blue_marker, 'MarkerFaceColor', styles.colors.blue_marker);
    grid(ax_global, 'on'); 
    box(ax_global, 'on'); 
    hold(ax_global, 'on');
    set(ax_global, 'Color', styles.colors.background, 'GridColor', styles.colors.grid_color);

    % Set title for the global plot.
    fig_title_str_global = sprintf([lang.plots.titles.threshold_convergence_long_n_g], config.bootstrap_thresholds.n_trials);
    title(tcl_thr_global, fig_title_str_global, 'FontSize', styles.font.title, 'FontWeight', 'bold', 'Color', styles.colors.text);

    % Set axis labels and limits.
    xlim(ax_global, [min(B_vector), max(B_vector) * 1.1]);
    xlabel(ax_global, lang.plots.xlabels.bootstraps, 'FontSize', styles.font.label, 'Color', styles.colors.text);
    ylabel(ax_global, lang.plots.ylabels.stability, 'FontSize', styles.font.label, 'Color', styles.colors.text); 
    set(ax_global, 'FontSize', styles.font.tick, 'XColor', styles.colors.text, 'YColor', styles.colors.text);

    % If robust convergence was used, also plot the smoothed curve (as this was used for the check).
    if use_robust_convergence_thr
        smoothing_window = config.bootstrap_thresholds.smoothing_window;
        smoothed_curve_plotted_global = movmean(overall_stability_thr_plotted, smoothing_window, 'omitnan');
        p2_global = plot(ax_global, B_vector, smoothed_curve_plotted_global * 100, '-', 'LineWidth', 1.5, 'Color', [0.8500 0.3250 0.0980]);
        lgd_global = legend(ax_global, [p1_global, p2_global], {lang.plots.legend.unsmoothed, lang.plots.legend.smoothed}, 'Location', 'best', ...
            'FontSize', styles.font.small_text);
        set(lgd_global, 'Color', styles.colors.background, 'TextColor', styles.colors.text, 'EdgeColor', styles.colors.text);
        lgd_global.ItemTokenSize = [15, 18];
    end

    % Add a marker for the optimal B value.
    selected_B_idx_for_plot_global = find(B_vector == selected_B, 1);
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
    exportgraphics(tcl_thr_global, filename_global, 'Resolution', 300, 'BackgroundColor', styles.colors.background, 'Padding', 30);
    fprintf([lang.thresholds.convergence_plot_saved '\n'], filename_global);

    %% Plot 2: Detailed Convergence Curves (per Metric) 
    % This plot shows the individual stability curves for each metric (without smoothing).
    h_fig_detailed = figure('Name', lang.plots.titles.threshold_convergence, 'Color', styles.colors.background, 'Visible', 'off');
    
    stability_matrix = stability_data.detailed_stability;
    elbow_indices = [];
    if isfield(stability_data, 'elbow_indices')
        elbow_indices = stability_data.elbow_indices;
    end
    converged = false;
    if isfield(stability_data, 'converged')
        converged = stability_data.converged;
    end

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
            p1 = plot(ax, B_vector, current_stability_curve_plotted * 100, '-o', 'LineWidth', 1.5, 'MarkerSize', 5, ...
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
            xlim_vals = [min(B_vector), max(B_vector)];
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
                    x_local = B_vector(local_elbow_idx);
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
            selected_B_idx_for_plot = find(B_vector == selected_B, 1);
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
    exportgraphics(tcl_thr_detailed, filename_detailed, 'Resolution', 300, 'BackgroundColor', styles.colors.background, 'Padding', 30);
    fprintf([lang.thresholds.convergence_plot_saved '\n'], filename_detailed);

end
