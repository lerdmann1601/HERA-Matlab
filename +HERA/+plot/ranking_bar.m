function h_fig_rank = ranking_bar(metric_idx, order, figure_name_str, y_label_str, title_base_str, filename_base_str, shared_info, styles, lang, subfolder_path)
% RANKING_BAR - Generates and saves intermediate ranking plots.
%
% Syntax:
%   h_fig = HERA.plot.ranking_bar(metric_idx, order, figure_name_str, ...
%           y_label_str, title_base_str, filename_base_str, shared_info, styles, lang, subfolder_path)
%
% Description:
%   This function creates a horizontal bar chart visualizing the ranking of datasets
%   based on a specific metric. It displays mean values with error bars illustrating
%   standard deviations. The resulting figure is saved as a high-resolution PNG
%   and the handle is returned for later PDF compilation.
%
% Inputs:
%   metric_idx        - (int) Index of the metric to plot (1 to N).
%   order             - (array) Permutation vector defining the ranking order of datasets.
%   figure_name_str   - (string) Name/Title property for the figure window.
%   y_label_str       - (string) Label for the Y-axis (e.g., "Rank (Metric 1)").
%   title_base_str    - (string) Formatted string for the plot title (e.g., "Ranking after %s").
%   filename_base_str - (string) Base name for the output graphics file.
%   shared_info       - (struct) Contains general data (names, paths, config timestamp).
%   styles            - (struct) Contains color and style definitions.
%   lang              - (struct) Language pack for labels and messages.
%   subfolder_path    - (string) Directory path where the plot should be saved.
%
% Outputs:
%   h_fig_rank        - (handle) Handle to the created figure object.
%
% Author: Lukas von Erdmannsdorff

% Unpack necessary variables from the shared_info struct.
num_datasets = shared_info.num_datasets;
dataset_names = shared_info.dataset_names;
metric_names = shared_info.metric_names;

% Prepare the data for the plot, sorted according to the provided order.
means = shared_info.mean_metrics(order, metric_idx);
stds = shared_info.std_metrics(order, metric_idx);
names = dataset_names(order);

% Create the figure object.
h_fig_rank = figure('Name', figure_name_str, 'Color', styles.colors.background, 'Visible', 'off');

tcl = tiledlayout(1, 1, 'Padding', 'compact');
ax = nexttile;

% Generate a horizontal bar chart of the mean values.
barh(ax, 1:num_datasets, means, 'FaceColor', styles.colors.p_significant, 'BarWidth', 0.7);
hold(ax, 'on');

% Add error bars representing standard deviation.
errorbar(ax, means, 1:num_datasets, stds, styles.errorbar.intermediate{:});

% Adjust axes properties (rank 1 at the top, labels, etc.).
set(ax, 'YDir', 'reverse', 'FontSize', styles.font.tick, 'YTick', 1:num_datasets, 'YTickLabel', names, ...
    'XColor', styles.colors.text, 'YColor', styles.colors.text, 'TickLabelInterpreter', 'none');

xlabel(ax, sprintf(lang.plots.xlabels.mean_sd, metric_names{metric_idx}), 'FontSize', styles.font.label, ...
    'Color', styles.colors.text, 'Interpreter', 'none');
ylabel(ax, y_label_str, 'FontSize', styles.font.label, 'Color', styles.colors.text);

title_str = sprintf([title_base_str '\n'], metric_names{metric_idx}); 
title(tcl, title_str, 'FontSize', ...
    styles.font.title, 'FontWeight', 'bold', 'Color', styles.colors.text, 'Interpreter', 'none');

ylim(ax, [0.5, num_datasets + 0.5]);
grid(ax, 'on'); box(ax, 'on');
ax.GridLineStyle = styles.line.grid;
hold(ax, 'off');

% Save the graphic.
[~, fName, fExt] = fileparts(filename_base_str);
filename = fullfile(subfolder_path, [fName, '_', shared_info.config.timestamp, fExt]);
exportgraphics(h_fig_rank, filename, 'Resolution', 300, 'BackgroundColor', styles.colors.background, 'Padding', 30);
fprintf(['  ' lang.plots.log_messages.graphic_saved '\n'], filename);

end
