function [h_fig_global, h_fig_detailed] = bca_convergence(B_vector, overall_stability, stability_matrix, selected_B, config, metric_names, graphics_dir, styles, lang, elbow_indices)
% BCA_CONVERGENCE - visualizes the stability analysis of the Bootstrap process.
%
% Syntax:
%   [h_fig_global, h_fig_detailed] = bca_convergence(B_vector, overall_stability, ...
%       stability_matrix, selected_B, config, metric_names, graphics_dir, styles, lang, elbow_indices)
%
% Description:
%   This function creates two figures to diagnose the convergence of the BCa bootstrapping:
%   1. Global Convergence Plot: Shows the aggregated stability metric across all effect sizes
%      and metrics. If 'Robust Convergence' is enabled, it overlays the smoothed curve.
%   2. Detailed Convergence Plot: A tiled layout showing the stability curve for each
%      individual metric and effect type (Delta vs RelDiff). If convergence failed and
%      the elbow method was used, it visualizes the local elbow points.
%
% Inputs:
%   B_vector          - Vector of the bootstrap step counts that were tested (X-axis).
%   overall_stability - The scalar stability metric (averaged) for each step in B_vector.
%   stability_matrix  - 3D Matrix of stability values [EffectTypes x Metrics x Steps].
%   selected_B        - The final chosen B value (optimal or elbow).
%   config            - Configuration struct (must contain .bootstrap_ci fields).
%   metric_names      - Cell array of strings with names of metrics.
%   graphics_dir      - Directory path where plots will be saved.
%   styles            - Struct defining colors, fonts, and other visual properties.
%   lang              - Struct containing localized strings for titles and labels.
%   elbow_indices     - (Optional) Vector of indices indicating local elbow points for detailed plots.
%
% Outputs:
%   h_fig_global      - Handle to the Global Convergence Figure.
%   h_fig_detailed    - Handle to the Detailed Tiled Layout Figure.
%
% Author: Lukas von Erdmannsdorff

    h_fig_global = gobjects(0);
    h_fig_detailed = gobjects(0);

    if nargin < 10
        elbow_indices = [];
    end

    % Extract sub-config for convenience
    cfg_ci = config.bootstrap_ci;
    ts = config.timestamp;
    
    % Determine if we used Robust Convergence (Smoothing)
    % This affects whether we plot a second smoothed line
    use_robust = isfield(cfg_ci, 'smoothing_window') && ~isempty(cfg_ci.smoothing_window) ...
                 && isfield(cfg_ci, 'convergence_streak_needed') && ~isempty(cfg_ci.convergence_streak_needed);

    % Create output subdirectory if it doesn't exist
    subfolder_bca_CI = fullfile(graphics_dir, 'CI_Histograms');
    if ~exist(subfolder_bca_CI, 'dir')
        mkdir(subfolder_bca_CI);
    end

    % Set default fonts for this session/plots
    set(0, 'DefaultAxesFontName', 'Arial');
    set(0, 'DefaultTextFontName', 'Arial');


    %% Plot 1: Global Convergence Curve
    % This plot summarizes the entire process into one line (mean stability).
    
    h_fig_global = figure('Name', lang.plots.titles.bca_convergence_global, 'Color', styles.colors.background, 'Visible', 'off');
    tcl_bca_global = tiledlayout(1, 1, 'Padding', 'compact');
    ax_global = nexttile;

    % 1. Plot Raw Stability
    p1_global = plot(ax_global, B_vector, overall_stability * 100, '-o', 'LineWidth', 1.5, 'MarkerSize', 5, ...
        'Color', styles.colors.blue_marker, 'MarkerFaceColor', styles.colors.blue_marker);
    
    % Basic Axis Setup
    grid(ax_global, 'on'); 
    box(ax_global, 'on'); 
    hold(ax_global, 'on');
    set(ax_global, 'Color', styles.colors.background, 'GridColor', styles.colors.grid_color);

    % 2. Add Title and Labels
    fig_title_str_global = sprintf([lang.plots.titles.bca_convergence_long_n_g], cfg_ci.n_trials);
    title(tcl_bca_global, fig_title_str_global, 'FontSize', styles.font.title, 'FontWeight', 'bold', 'Color', styles.colors.text);

    xlim(ax_global, [min(B_vector), max(B_vector) * 1.1]);
    xlabel(ax_global, lang.plots.xlabels.bootstraps, 'FontSize', styles.font.label, 'Color', styles.colors.text);
    ylabel(ax_global, lang.plots.ylabels.stability, 'FontSize', styles.font.label, 'Color', styles.colors.text); 
    set(ax_global, 'FontSize', styles.font.tick, 'XColor', styles.colors.text, 'YColor', styles.colors.text);

    % 3. Optional: Plot Smoothed Curve (if Robust Mode was active)
    if use_robust
        smoothing_window_global = cfg_ci.smoothing_window;
        % Re-calculate smoothing for visualization using same window
        smoothed_curve = movmean(overall_stability, smoothing_window_global, 'omitnan');
        
        p2_global = plot(ax_global, B_vector, smoothed_curve * 100, '-', 'LineWidth', 1.5, 'Color', styles.colors.convergence);
        
        lgd_global = legend(ax_global, [p1_global, p2_global], {lang.plots.legend.unsmoothed, lang.plots.legend.smoothed}, ...
            'Location', 'best', 'FontSize', styles.font.small_text);
        set(lgd_global, 'Color', styles.colors.background, 'TextColor', styles.colors.text, 'EdgeColor', styles.colors.text);
        lgd_global.ItemTokenSize = [15, 18];
    end

    % 4. Mark the Optimally Selected B
    idx_plot = find(B_vector == selected_B, 1);
    if ~isempty(idx_plot)
        x_pos = selected_B;
        y_pos = overall_stability(idx_plot) * 100;
        
        % Plot 'X' marker
        plot(ax_global, x_pos, y_pos, 'x', 'Color', styles.colors.red_marker, 'MarkerSize', 10, 'LineWidth', 2, 'HandleVisibility', 'off');
        
        % Adjust X-limits if marker is too close to edge
        cur_xlim = xlim(ax_global);
        if x_pos > (cur_xlim(1) + (cur_xlim(2) - cur_xlim(1)) * 0.85)
            xlim(ax_global, [cur_xlim(1), x_pos * 1.2]);
        end
        % Text Annotation
        text(ax_global, x_pos, y_pos, sprintf(lang.plots.misc.optimal_b_text, selected_B), ...
             'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', styles.font.small_text, 'Color', styles.colors.red_marker);
    end
    hold(ax_global, 'off');

    % 5. Save Global Plot
    [~, fName, fExt] = fileparts(lang.files.convergence_bca_global);
    filename_global = fullfile(subfolder_bca_CI, [fName, '_', ts, fExt]);
    exportgraphics(h_fig_global, filename_global, 'Resolution', 300, 'BackgroundColor', styles.colors.background, 'Padding', 30);
    fprintf([lang.bca.convergence_plot_saved '\n'], filename_global);


    %% Plot 2: Detailed Convergence (Grid)
    h_fig_detailed = figure('Name', lang.plots.titles.bca_convergence_long, 'Color', styles.colors.background, 'Visible', 'off');
    
    num_metrics = numel(metric_names);
    tcl_detailed = tiledlayout(2, num_metrics, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    fig_title_det = sprintf([lang.plots.titles.bca_convergence_long_n_d], cfg_ci.n_trials);
    title(tcl_detailed, fig_title_det, 'FontSize', styles.font.title, 'FontWeight', 'bold', 'Color', styles.colors.text);

    num_effect_types = 2;
    effect_type_names = {'Cliff''s Delta', 'Rel Diff'};
    color_local_elbow = styles.colors.convergence;

    % Iterate: Effect Type (Rows) x Metric (Cols)
    for es_type = 1:num_effect_types
        for m_idx = 1:num_metrics
            ax = nexttile;
            
            handles_leg = [];
            names_leg = {};
            
            % Extract Curve Data
            % stability_matrix dims: [2, num_metrics, steps]
            cur_curve = squeeze(stability_matrix(es_type, m_idx, :));
            
            % Slicing safety: Ensure curve length matches B_vector length
            if numel(cur_curve) > numel(B_vector)
                cur_curve = cur_curve(1:numel(B_vector));
            elseif numel(cur_curve) < numel(B_vector)
                % Can happen if B_vector has not been sliced to final index outside
                % Handle gracefully by slicing B_vector locally
                 B_local = B_vector(1:numel(cur_curve));
            else
                 B_local = B_vector;
            end
            
            % 1. Plot Individual Function
            p1 = plot(ax, B_local, cur_curve * 100, '-o', 'LineWidth', 1.5, 'MarkerSize', 5, ...
                'Color', styles.colors.blue_marker, 'MarkerFaceColor', styles.colors.blue_marker, ...
                'DisplayName', lang.plots.legend.unsmoothed);
            
            handles_leg(end+1) = p1;
            names_leg{end+1} = lang.plots.legend.unsmoothed;
            
            % Styling
            grid(ax, 'on'); box(ax, 'on'); hold(ax, 'on');
            set(ax, 'Color', styles.colors.background, 'GridColor', styles.colors.grid_color);
            
            subvals_title = sprintf('%s - %s\n', metric_names{m_idx}, effect_type_names{es_type});
            title(ax, subvals_title, 'FontSize', styles.font.label, 'Color', styles.colors.text, 'Interpreter', 'none');
            
            xlim(ax, [min(B_local), max(B_local) * 1.1]);
            xlabel(ax, lang.plots.xlabels.bootstraps, 'FontSize', styles.font.label, 'Color', styles.colors.text);
            ylabel(ax, lang.plots.ylabels.stability, 'FontSize', styles.font.label, 'Color', styles.colors.text);
            set(ax, 'FontSize', styles.font.tick, 'XColor', styles.colors.text, 'YColor', styles.colors.text);
            
            % 2. Visualization of Local Elbow (Only if Elbow Method was used)
            if ~isempty(elbow_indices)
                % Calculate linear index in elbow vector for this subplot
                % Map: [D1 D2 D3 ... R1 R2 R3 ...]
                if es_type == 1
                    k_elbow = m_idx;
                else
                    k_elbow = m_idx + num_metrics;
                end
                
                if k_elbow <= numel(elbow_indices)
                    local_idx = elbow_indices(k_elbow);
                    if local_idx <= numel(B_local)
                        x_local = B_local(local_idx);
                        y_local = cur_curve(local_idx) * 100;
                        
                        % Dashed Marker
                        p_elbow = plot(ax, x_local, y_local, ':o', ...
                            'MarkerFaceColor', color_local_elbow, ...
                            'MarkerEdgeColor', color_local_elbow, ...
                            'MarkerSize', 6, 'HandleVisibility', 'on', ...
                            'LineWidth', 1.5, ...
                            'DisplayName', lang.plots.legend.local_elbow);
                        
                        % Vertical Guide Line
                        xline(ax, x_local, '--', 'Color', color_local_elbow, 'LineWidth', 1.0, 'HandleVisibility', 'off');
                        
                        handles_leg(end+1) = p_elbow;
                        names_leg{end+1} = lang.plots.legend.local_elbow;
                    end
                end
            end
            
            % 3. Mark Selected B (Global Decision)
            idx_sel = find(B_local == selected_B, 1);
            if ~isempty(idx_sel)
                x_pos = selected_B;
                y_pos = cur_curve(idx_sel) * 100;
                
                plot(ax, x_pos, y_pos, 'x', 'Color', styles.colors.red_marker, 'MarkerSize', 10, 'LineWidth', 2, 'HandleVisibility', 'off');
                
                cl = xlim(ax);
                if x_pos > (cl(1) + (cl(2) - cl(1)) * 0.85)
                    xlim(ax, [cl(1), x_pos * 1.2]);
                end
                text(ax, x_pos, y_pos, sprintf(lang.plots.misc.optimal_b_text, selected_B), ...
                     'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', ...
                     'FontSize', styles.font.small_text, 'Color', styles.colors.red_marker);
            end
            
            % 4. Dynamic Legend
            lgd = legend(handles_leg, names_leg, 'Location', 'best');
            set(lgd, 'Color', styles.colors.background, 'TextColor', styles.colors.text, 'EdgeColor', styles.colors.text);
            lgd.FontSize = styles.font.small_text;
            lgd.ItemTokenSize = [15, 18];
            
            hold(ax, 'off');
        end
    end

    % 5. Save Detailed Plot
    [~, fName, fExt] = fileparts(lang.files.convergence_bca);
    filename_det = fullfile(subfolder_bca_CI, [fName, '_', ts, fExt]);
    exportgraphics(h_fig_detailed, filename_det, 'Resolution', 300, 'BackgroundColor', styles.colors.background, 'Padding', 30);
    fprintf([lang.bca.convergence_plot_saved '\n\n'], filename_det);
end
