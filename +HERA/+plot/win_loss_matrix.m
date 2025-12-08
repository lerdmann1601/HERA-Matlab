function h_fig_winloss = win_loss_matrix(results, shared_info, styles, lang, subfolder_ranking)
% WIN_LOSS_MATRIX - Creates a matrix visualization summarizing all pairwise comparisons.
%
% Syntax:
%   h_fig_winloss = HERA.plot.win_loss_matrix(results, shared_info, styles, lang, subfolder_ranking)
%
% Description:
%   This function generates a win-loss matrix for each metric. The matrix visualizes
%   which dataset "won" against another based on the ranking criteria.
%   If power analysis results are available, the color intensity of the cells
%   reflects the statistical power of the comparison.
%   The resulting figure is saved as a high-resolution PNG.
%
% Inputs:
%   results           - (struct) Contains all calculation results (sig matrices, power, etc.).
%   shared_info       - (struct) General info (names, dataset count, timestamp).
%   styles            - (struct) Style definitions.
%   lang              - (struct) Language pack.
%   subfolder_ranking - (string) Path to save the PNG file.
%
% Outputs:
%   h_fig_winloss     - (handle) Handle to the created figure object.
%
% Author: Lukas von Erdmannsdorff

% Unpack necessary variables.
metric_names = shared_info.metric_names;
dataset_names = shared_info.dataset_names;
num_datasets = shared_info.num_datasets;
pair_idx_all = shared_info.pair_idx_all;
ts = shared_info.config.timestamp;
num_metrics = numel(metric_names);

% Define color levels and gray values based on style settings.
has_power_results = isfield(results, 'power_results') && ~isempty(results.power_results);
intensities = styles.win_loss.intensities;
gray_values = styles.win_loss.gray_values;
diag_color  = styles.win_loss.diag_color;
% Define edges for assigning power levels (1 to 4).
edges = [-inf, 0.25, 0.5, 0.75, inf];

% Set up the figure, tiled layout, and main title.
h_fig_winloss = figure('Name', lang.plots.figure_names.win_loss_matrix, 'Color', styles.colors.background, 'Visible', 'off');
% Tiled layout is dynamic based on num_metrics
tcl_wl = tiledlayout(1, num_metrics, 'Padding', 'loose', 'TileSpacing', 'compact');
title_str_winloss = sprintf([lang.plots.titles.win_loss_matrix '\n']);
title(tcl_wl, title_str_winloss, 'FontSize', styles.font.title, 'FontWeight', 'bold', 'Color', styles.colors.text, 'Interpreter', 'none');

% Loop through each metric to create a separate win-loss matrix.
for metric_idx = 1:num_metrics
    ax = nexttile;
    sig_matrix = results.all_sig_matrices{metric_idx};
    % Convert the power results vector into a full matrix if available.
    power_matrix_full = zeros(num_datasets);
    if has_power_results
        power_vals_vec = results.power_results.power_matrices{metric_idx};
        for p_idx = 1:size(pair_idx_all, 1)
            i_p = pair_idx_all(p_idx, 1);
            j_p = pair_idx_all(p_idx, 2);
            power_matrix_full(i_p, j_p) = power_vals_vec(p_idx);
            power_matrix_full(j_p, i_p) = power_vals_vec(p_idx);
        end
    end
    win_loss_image = ones(num_datasets, num_datasets, 3);
  
    % Iterate through each pair to determine the color of the matrix cell.
    for i = 1:num_datasets
        for j = 1:num_datasets
            if i == j
                % Set the diagonal color.
                win_loss_image(i, j, :) = diag_color;
            else
                % Determine the power level and corresponding color.
                power = power_matrix_full(i, j);
                level = find(power < edges, 1, 'first') - 1;
                if isempty(level) || level == 0, level = 1; end
                if sig_matrix(i, j) % Win
                    col = styles.colors.win * intensities(level);
                elseif sig_matrix(j, i) % Loss
                    col = styles.colors.loss * intensities(level);
                else % Neutral
                    col = ones(1,3) * gray_values(level);
                end
                win_loss_image(i, j, :) = col;
            end
        end
    end
    % Display the colored matrix as an image.
    imagesc(ax, win_loss_image);
    axis(ax, 'square');
  
    % Style the axes of the matrix and set tick labels.
    set(ax, 'XTick', 1:num_datasets, 'XTickLabel', dataset_names, 'XTickLabelRotation', 90, ...
        'YTick', 1:num_datasets, 'YTickLabel', dataset_names, ...
        'TickLabelInterpreter', 'none', ...
        'FontSize', styles.font.tick, 'XAxisLocation', 'top', ...
        'Box', 'on', ...
        'XColor', styles.colors.text, 'YColor', styles.colors.text, ...
        'GridColor', styles.colors.grid_color, 'GridAlpha', 0.5, ...
        'XGrid', 'on', 'YGrid', 'on', ...
        'Layer', 'top');
    % Set the subplot title.
    title(ax, metric_names{metric_idx}, 'FontSize', styles.font.label, 'Color', styles.colors.text, 'Interpreter', 'none');
end

% Create a legend for the power levels if power analysis was run.
if has_power_results
    % Adjust the layout to make space for the legend.
    tcl_wl.OuterPosition = [0.03, 0.30, 0.94, 0.65];
    ax_legend_y_start = 0.08;
    ax_legend_height = 0.20;
    % Create a new axes object for the legend.
    ax_legend = axes('Parent', h_fig_winloss, 'Position', [0.20, ax_legend_y_start, 0.60, ax_legend_height], ...
                      'Visible', 'off', 'XLim', [0 1], 'YLim', [0 1]);
    hold(ax_legend, 'on');
   
    % Define legend properties.
    legend_intensities = styles.win_loss.legend_intensities;
    legend_gray_values = styles.win_loss.legend_gray_values;
    n_levels = 4;
    power_bins = 0:0.25:1;
  
    % Define the layout structure for legend elements.
    layout = struct('x_cols', linspace(0.35, 0.65, 4), 'y_header', 0.95, ...
                    'y_content_start', 0.75, 'y_content_end', 0.05, 'marker_size', 200);
    layout.y_positions = linspace(layout.y_content_start, layout.y_content_end, n_levels);
    % Add legend headers.
    headers = {lang.plots.legend.power, lang.plots.legend.win, lang.plots.legend.loss, lang.plots.legend.neutral};
    for k = 1:numel(headers)
        text(ax_legend, layout.x_cols(k), layout.y_header, headers{k}, ...
            'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
            'Color', styles.colors.text, 'FontSize', styles.font.label);
    end
  
    % Add legend content (power ranges and colored markers).
    for lvl = 1:n_levels
        color_win  = styles.colors.win  * legend_intensities(lvl);
        color_loss = styles.colors.loss * legend_intensities(lvl);
        color_neu  = ones(1,3) * legend_gray_values(lvl);
        range_str = sprintf('%d - %d %%', power_bins(lvl)*100, power_bins(lvl+1)*100);
        text(ax_legend, layout.x_cols(1), layout.y_positions(lvl), range_str, 'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'middle', 'Color', styles.colors.text, 'FontSize', styles.font.tick);
        scatter(ax_legend, layout.x_cols(2), layout.y_positions(lvl), layout.marker_size, color_win,  's', 'filled', 'MarkerEdgeColor', styles.colors.text);
        scatter(ax_legend, layout.x_cols(3), layout.y_positions(lvl), layout.marker_size, color_loss, 's', 'filled', 'MarkerEdgeColor', styles.colors.text);
        scatter(ax_legend, layout.x_cols(4), layout.y_positions(lvl), layout.marker_size, color_neu,  's', 'filled', 'MarkerEdgeColor', styles.colors.text);
    end
    hold(ax_legend, 'off');
end

% Save the win-loss matrix graphic.
[~, fName, fExt] = fileparts(lang.files.win_loss_matrix);
filename_winloss = fullfile(subfolder_ranking, [fName, '_', ts, fExt]);
exportgraphics(h_fig_winloss, filename_winloss, 'Resolution', 300, 'Padding', 30);
fprintf(['  ' lang.plots.log_messages.graphic_saved '\n'], filename_winloss);

end
