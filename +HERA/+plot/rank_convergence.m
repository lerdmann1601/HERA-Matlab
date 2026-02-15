function h_fig_rank = rank_convergence(B_vector, stability_vector, selected_B, converged, elbow_idx, config, styles, lang)
% RANK_CONVERGENCE - Generates the convergence plot for rank stability.
%
% Syntax:
%   h_fig_rank = HERA.plot.rank_convergence(B_vector, stability_vector, selected_B, ...
%                converged, elbow_idx, config, styles, lang)
%
% Description:
%   This function creates a plot showing the stability of the ranking confidence intervals
%   as a function of the number of bootstrap samples (B). It visualizes the raw stability
%   curve and, if enabled, a smoothed version. It also marks the selected optimal B-value,
%   either from convergence or the elbow method.
%
% Inputs:
%   B_vector         - Vector of tested bootstrap repetition numbers (x-axis).
%   stability_vector - Vector of calculated stability values (y-axis).
%   selected_B       - The final selected number of bootstrap samples.
%   converged        - Boolean indicating if convergence was reached.
%   elbow_idx        - (Optional) Index of the elbow point if convergence failed.
%   config           - Global configuration struct (for smoothing settings).
%   styles           - Struct containing settings for graphical design.
%   lang             - Language pack struct for text outputs.
%
% Outputs:
%   h_fig_rank       - Handle to the generated figure.
%
% Author: Lukas von Erdmannsdorff

    % Check if robust convergence logic is applicable for visualization
    cfg_rank = config.bootstrap_ranks;
    use_robust_convergence = isfield(cfg_rank, 'smoothing_window') && ~isempty(cfg_rank.smoothing_window) ...
                            && isfield(cfg_rank, 'convergence_streak_needed') && ~isempty(cfg_rank.convergence_streak_needed);

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
    p1 = plot(ax, B_vector, stability_vector, '-o', ...
         'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', styles.colors.blue_marker, 'Color', styles.colors.blue_marker);
    handles_to_legend(end+1) = p1;
    names_to_legend{end+1} = lang.plots.legend.unsmoothed;
         
    grid(ax, 'on'); box(ax, 'on'); hold(ax, 'on'); 
    set(ax, 'Color', styles.colors.background, 'GridColor', styles.colors.grid_color);
    
    % If robust convergence was used, also plot the smoothed curve.
    if use_robust_convergence
        smoothing_window = cfg_rank.smoothing_window;
        smoothed_stability = movmean(stability_vector, smoothing_window, 'omitnan');
        p2 = plot(ax, B_vector, smoothed_stability, '-', 'LineWidth', 1.5, 'Color', styles.colors.convergence);
        handles_to_legend(end+1) = p2;
        names_to_legend{end+1} = lang.plots.legend.smoothed;
    end
    
    % Set titles and labels for the plot.
    fig_title_str = sprintf(lang.plots.titles.rank_stability_convergence_n, cfg_rank.n_trials);
    title(tcl_rank, fig_title_str, 'Color', styles.colors.text, 'FontSize', styles.font.title, 'FontWeight', 'bold');
    
    xlim(ax, [min(B_vector), max(B_vector) * 1.2]);
    xlabel(ax, lang.plots.xlabels.bootstraps, 'Color', styles.colors.text, 'FontSize', styles.font.label);
    ylabel(ax, lang.plots.ylabels.stability_rank, 'Color', styles.colors.text, 'FontSize', styles.font.label); 
    set(ax, 'FontSize', styles.font.tick, 'XColor', styles.colors.text, 'YColor', styles.colors.text);

    % Plot the LOCAL elbow (if convergence failed and elbow index is valid)
    if ~converged && ~isempty(elbow_idx)
        local_elbow_idx = elbow_idx; 
        
        if local_elbow_idx <= numel(stability_vector)
            x_local = B_vector(local_elbow_idx);
            y_local = stability_vector(local_elbow_idx);
            
            % Define the color for the local elbow (same as smoothed curve)
            color_local_elbow = styles.colors.convergence; 
    
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
    selected_B_idx = find(B_vector == selected_B, 1);
    if ~isempty(selected_B_idx)
        x_pos = selected_B;
        y_pos = stability_vector(selected_B_idx);
        
        % Plot the 'x' marker, but hide it from the legend.
        plot(ax, x_pos, y_pos, 'x', 'Color', styles.colors.red_marker, ...
             'MarkerSize', 10, 'LineWidth', 2, 'HandleVisibility', 'off');
    
        % Extend the x-axis margin if the marker is too close to the edge.
        current_xlim = xlim(ax);
        if x_pos > (current_xlim(1) + (current_xlim(2) - current_xlim(1)) * 0.85)
            xlim(ax, [current_xlim(1), x_pos * 1.2]);
        end
    
        % Add text label for the final B value (this explains the 'x')
        text(ax, x_pos, y_pos, sprintf(lang.plots.misc.optimal_b, selected_B), 'VerticalAlignment', 'bottom', ...
             'HorizontalAlignment', 'left', 'FontSize', styles.font.tick, 'Color', styles.colors.red_marker);
    end
    
    % Create the final, dynamic legend 
    lgd = legend(handles_to_legend, names_to_legend, 'Location', 'best');
    set(lgd, 'Color', styles.colors.background, 'TextColor', styles.colors.text, 'EdgeColor', styles.colors.text);
    lgd.FontSize = styles.font.tick;
    lgd.ItemTokenSize = [15, 18];
    
    hold(ax, 'off');
end
