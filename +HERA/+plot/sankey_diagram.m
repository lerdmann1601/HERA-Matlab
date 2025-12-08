function h_fig_sankey = sankey_diagram(results, shared_info, styles, lang, subfolder_ranking)
% SANKEY_DIAGRAM - Visualize the changes in ranking across the ranking stages.
%
% Syntax:
%   h_fig_sankey = HERA.plot.sankey_diagram(results, shared_info, styles, lang, subfolder_ranking)
%
% Description:
%   This function creates a Sankey diagram if there is more than one ranking stage.
%   It visualizes how the position of each dataset shifts from one metric to the next.
%   The resulting figure is saved as a high-resolution PNG.
%   Returns an empty graphics object if no plot is generated (e.g., only 1 metric).
%
% Inputs:
%   results           - (struct) Contains intermediate ranking orders.
%   shared_info       - (struct) General info (metric names, ranking mode, timestamp).
%   styles            - (struct) Style definitions (colors, fonts).
%   lang              - (struct) Language pack.
%   subfolder_ranking - (string) Path to save the PNG file.
%
% Outputs:
%   h_fig_sankey      - (handle) Handle to the created figure or gobject(0).
%
% Author: Lukas von Erdmannsdorff

% Unpack necessary variables.
metric_names = shared_info.metric_names;
dataset_names = shared_info.dataset_names;
num_datasets = shared_info.num_datasets;
ts = shared_info.config.timestamp;
ranking_mode = shared_info.config.ranking_mode;
num_metrics = numel(metric_names);

% Only plot if there is more than one stage to show.
if num_metrics > 1
    
    h_fig_sankey = figure('Name', lang.plots.figure_names.rank_shifts, 'Color', styles.colors.background, 'Visible', 'off');
    tcl_sankey = tiledlayout(1, 1, 'Padding', 'compact');
    ax = nexttile;
    hold(ax, 'on');
    set(ax, 'Color', styles.colors.background);

    % Dynamically build ranks and stages based on logic
    stages = {results.intermediate_orders.after_metric1};
    x_coords = 1;
    x_labels = {sprintf(lang.plots.xlabels.sankey_m1, metric_names{1})};
    metric_name_list = {metric_names{1}};

    if num_metrics >= 2 && (strcmp(ranking_mode, 'M1_M2') || strcmp(ranking_mode, 'M1_M2_M3'))
        stages{end+1} = results.intermediate_orders.after_metric2;
        x_coords(end+1) = 2;
        x_labels{end+1} = sprintf(lang.plots.xlabels.sankey_m2_corr, metric_names{2});
        metric_name_list{end+1} = metric_names{2};
    elseif num_metrics == 2 && strcmp(ranking_mode, 'M1_M3A')
        stages{end+1} = results.intermediate_orders.after_metric3; % M3A order
        x_coords(end+1) = 2;
        x_labels{end+1} = sprintf(lang.plots.xlabels.sankey_m2_tie, metric_names{2});
        metric_name_list{end+1} = metric_names{2};
    end

    if num_metrics == 3 && strcmp(ranking_mode, 'M1_M2_M3')
        stages{end+1} = results.intermediate_orders.after_metric3;
        x_coords(end+1) = 3;
        x_labels{end+1} = sprintf(lang.plots.xlabels.sankey_m3_corr, metric_names{3});
        metric_name_list{end+1} = metric_names{3};
    end

    num_stages = numel(stages);
    all_ranks = zeros(num_datasets, num_stages);
    for s_idx = 1:num_stages
        all_ranks(stages{s_idx}, s_idx) = 1:num_datasets;
    end

    % Define color palette and band width for the Sankey diagram.
    colors = styles.sankey.colors;
    band_width = 0.4;

    % Draw the connecting bands for each dataset across the stages.
    for d = 1:num_datasets
        for s_idx = 1:(num_stages - 1)
            y1_start = all_ranks(d, s_idx) - band_width; 
            y1_end   = all_ranks(d, s_idx) + band_width;
            y2_start = all_ranks(d, s_idx+1) - band_width; 
            y2_end   = all_ranks(d, s_idx+1) + band_width;
            
            patch(ax, [x_coords(s_idx), x_coords(s_idx+1), x_coords(s_idx+1), x_coords(s_idx)], ...
                [y1_start, y2_start, y2_end, y1_end], colors(d, :), 'FaceAlpha', 0.7, 'EdgeColor', 'none');
        end
        
        % Add vertically centered labels for datasets at the start and end.
        text(ax, x_coords(1) - 0.05, all_ranks(d, 1), dataset_names{d}, 'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'middle', 'FontSize', styles.font.tick, 'FontWeight', 'bold', 'Color', styles.colors.text, 'Interpreter', 'none');
        text(ax, x_coords(end) + 0.05, all_ranks(d, end), dataset_names{d}, 'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'middle', 'FontSize', styles.font.tick, 'FontWeight', 'bold', 'Color', styles.colors.text, 'Interpreter', 'none');
    end

    % Configure the axes properties.
    set(ax, 'YDir', 'reverse', 'YTick', 1:num_datasets, 'YTickLabel', 1:num_datasets, ...
        'XTick', x_coords, 'XTickLabel', x_labels, 'XTickLabelRotation', 0, ...
        'FontSize', styles.font.tick, 'Box', 'on', ...
        'XColor', styles.colors.text, 'YColor', styles.colors.text, 'GridColor', styles.colors.grid_color, 'TickLabelInterpreter', 'none');
    xlim(ax, [min(x_coords)-0.5, max(x_coords)+0.5]); 
    ylim(ax, [0.5, num_datasets + 0.5]);

    % Set the left Y-axis label.
    ylabel(ax, lang.plots.ylabels.rank, 'FontSize', styles.font.label, 'FontWeight', 'bold', 'Color', styles.colors.text);
    % Create and configure a second Y-axis on the right for symmetry.
    yyaxis(ax, 'right');
    set(ax, 'YColor', styles.colors.text, 'YDir', 'reverse', 'YLim', [0.5, num_datasets + 0.5], ...
            'YTick', 1:num_datasets, 'YTickLabel', 1:num_datasets, 'FontWeight', 'bold');
    ylabel(ax, lang.plots.ylabels.rank, 'FontSize', styles.font.label, 'FontWeight', 'bold', 'Color', styles.colors.text);
    ax.YGrid = 'off';

    % Restore the left axis as the active one.
    yyaxis(ax, 'left');
    set(ax, 'YColor', styles.colors.text);

    % Set the final title for the Sankey diagram.
    title_str_sankey = sprintf(lang.plots.titles.sankey_diagram);
    title(tcl_sankey, sprintf([title_str_sankey '\n']), ...
        'FontSize', styles.font.title, 'FontWeight', 'bold', 'Color', styles.colors.text, 'Interpreter', 'none');
    hold(ax, 'off');

    % Save the Sankey diagram graphic.
    [~, fName, fExt] = fileparts(lang.files.sankey_rank_shifts);
    filename_sankey = fullfile(subfolder_ranking, [fName, '_', ts, fExt]);
    exportgraphics(h_fig_sankey, filename_sankey, 'Resolution', 300, 'BackgroundColor', styles.colors.background, 'Padding', 30);
    fprintf(['  ' lang.plots.log_messages.graphic_saved '\n'], filename_sankey);
else
    h_fig_sankey = gobjects(0); % Create empty handle if no plot is made
end

end
