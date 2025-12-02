function generate_plots(results, thresholds, shared_info, styles)
% GENERATE_PLOTS - Creates and saves all ranking graphics.
%
% Syntax:
%   generate_plots(results, thresholds, shared_info, styles)
%
% Description:
%   This function visualizes the results of the ranking algorithm by creating a series of graphics.
%   All figures are automatically saved as high-resolution PNG files.
%   Additionally, the graphics are grouped thematically and compiled into several PDF reports for easy sharing.
%
% Workflow:
%   1.  Creation of intermediate ranking plots: Visualizes the rankings after sorting by Metric 1.
%       Plots for Metric 2 and 3 are created conditionally based on the 'ranking_mode'.
%   2.  Creation of detailed comparison plots: Displays all pairwise comparisons (p-value, delta, rel. diff.) for each loaded metric.
%   3.  Creation of a win-loss matrix to visualize direct duels, dynamically showing 1, 2, or 3 matrices.
%   4.  Creation of a Sankey diagram to visualize rank shifts across the executed stages (dynamically 1, 2, or 3 stages).
%   5.  Creation of the final summary plot: Combines the final ranking (with bootstrap CI) with a table of metric values (dynamically 1-3 columns).
%   6.  Creation of the Borda plot (optional): Visualizes the result of the sensitivity analysis as a consensus ranking.
%   7.  Compilation into PDF reports: Saves the created graphics in thematic PDF files (e.g., Final_Report.pdf, Ranking_Report.pdf).
%
% Inputs:
%   results     - (struct) Contains all calculated results (ranks, CIs, p-values, etc.).
%   thresholds  - (struct) Contains the calculated statistical thresholds.
%   shared_info - (struct) Contains general data (names, paths) and configurations (incl. 'config' struct with 'ranking_mode').
%   styles      - (struct) Contains all color and style settings from the 'design' function.
%
% Outputs:
%   The function has no direct return values. 
%   Instead, it generates files:
%   Ranking plots           - PNG files of intermediate and final rankings.
%   Detail comparison plots - PNG files of pairwise comparisons for each metric.
%   Win-Loss Matrix         - PNG file of the duel overview.
%   Sankey Diagram          - PNG file of the rank shifts.
%   Borda Plot              - PNG file of the consensus ranking (if calculated).
%   PDF reports             - Several PDF files that bundle the graphics thematically.
%
% Author:   Lukas von Erdmannsdorff

%% 1. Initialization and Style Definition
% Unpack the passed structures into local variables for easier access.
metric_names = shared_info.metric_names;
dataset_names = shared_info.dataset_names;
num_datasets = shared_info.num_datasets;
pair_idx_all = shared_info.pair_idx_all;
d_thresh = thresholds.d_thresh;
rel_thresh = thresholds.rel_thresh;
output_dir = shared_info.output_dir;
graphics_dir = shared_info.graphics_dir;
pdf_dir      = shared_info.pdf_dir;
plot_theme = shared_info.plot_theme;
lang = shared_info.lang; % Load the full language pack.
ts = shared_info.config.timestamp; 

% Get dynamic metric count and ranking logic
num_metrics = numel(metric_names);
ranking_mode = shared_info.config.ranking_mode;

% Create subfolders for detailed comparison and ranking plots.
subfolder_detail_comp = fullfile(graphics_dir, 'Detail_Comparison');
if ~exist(subfolder_detail_comp, 'dir')
    mkdir(subfolder_detail_comp);
end
subfolder_ranking = fullfile(graphics_dir, 'Ranking');
if ~exist(subfolder_ranking, 'dir')
    mkdir(subfolder_ranking);
end

% Initialize arrays to collect all figure handles for PDF export.
final_report_handles = gobjects(0);
ranking_report_handles = gobjects(0);
convergence_report_handles = gobjects(0);

% Log the start of the plot creation process.
fprintf(['\n' lang.plots.log_messages.start_creation '\n'], upper(plot_theme));

% Set global background and font for consistency across all plots.
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');
set(0, 'DefaultAxesColor', styles.colors.background);

%% 2. Plot: Ranking by Metric 1
% Generate the initial ranking plot based on the primary metric.
% This is always created.
ranking_report_handles = ranking_plot(1, ...
    results.intermediate_orders.after_metric1, ...
    lang.plots.figure_names.initial_ranking, ...
    lang.plots.ylabels.rank_metric_1, ...
    lang.plots.titles.initial_ranking, ...
    lang.files.initial_ranking_metric1, ...
    shared_info, styles, lang, subfolder_detail_comp, ranking_report_handles);

%% 3. Plot: Detailed Pairwise Comparisons for Metric 1
% Generate the detailed comparison plots for the primary metric.
% This is always created.
ranking_report_handles = detail_comparison_plot(1, results, thresholds, shared_info, styles, lang, subfolder_detail_comp, ranking_report_handles);

%% 4. Plot: Ranking by Metric 2
% Conditional plotting based on logic
if num_metrics >= 2 && (strcmp(ranking_mode, 'M1_M2') || strcmp(ranking_mode, 'M1_M2_M3'))
    % Generate the ranking plot after corrections based on the second metric.
    ranking_report_handles = ranking_plot(2, ...
        results.intermediate_orders.after_metric2, ...
        lang.plots.figure_names.rank_correction_m2, ...
        lang.plots.ylabels.rank_metric_2, ...
        sprintf(lang.plots.titles.ranking_after_correction_m2, metric_names{2}), ... 
        lang.files.correction_metric2, ...
        shared_info, styles, lang, subfolder_detail_comp, ranking_report_handles);
end

%% 5. Plot: Detailed Pairwise Comparisons for Metric 2
% Conditional plotting (if metric 2 exists, but NOT for M1_M3A, which is handled in Sec 6)
if num_metrics >= 2 && ~strcmp(ranking_mode, 'M1_M3A')
    % Generate the detailed comparison plots for the second metric.
    ranking_report_handles = detail_comparison_plot(2, results, thresholds, shared_info, styles, lang, subfolder_detail_comp, ranking_report_handles);
end

%% 6. Plot: Ranking by Metric 3
% Conditional plotting based on logic
if num_metrics == 2 && strcmp(ranking_mode, 'M1_M3A')
    % Plot the final ranking after M3A logic (which used M2 data)
    ranking_report_handles = ranking_plot(2, ... % Uses M2 data
        results.intermediate_orders.after_metric3, ... % Uses final order
        lang.plots.figure_names.rank_correction_m3, ...
        lang.plots.ylabels.rank_metric_2, ...
        sprintf(lang.plots.titles.ranking_after_correction_m3a, metric_names{2}), ... 
        lang.files.correction_metric3, ... % Saves as M3 correction
        shared_info, styles, lang, subfolder_detail_comp, ranking_report_handles);
    ranking_report_handles = detail_comparison_plot(2, results, thresholds, shared_info, styles, lang, subfolder_detail_comp, ranking_report_handles);
        
elseif num_metrics == 3 && strcmp(ranking_mode, 'M1_M2_M3')
    % Plot the final ranking after M3 (A+B) logic
    % Generate the ranking plot after corrections based on the third metric.
    ranking_report_handles = ranking_plot(3, ...
        results.intermediate_orders.after_metric3, ...
        lang.plots.figure_names.rank_correction_m3, ...
        lang.plots.ylabels.rank_metric_3, ...
        sprintf(lang.plots.titles.ranking_after_correction_m3, metric_names{3}), ... 
        lang.files.correction_metric3, ...
        shared_info, styles, lang, subfolder_detail_comp, ranking_report_handles);
end

%% 7. Plot: Detailed Pairwise Comparisons for Metric 3
% Conditional plotting (if metric 3 exists, regardless of logic)
if num_metrics == 3
    % Generate the detailed comparison plots for the third metric.
    ranking_report_handles = detail_comparison_plot(3, results, thresholds, shared_info, styles, lang, subfolder_detail_comp, ranking_report_handles);
end

%% 8. Plot: Win-Loss-Matrix
% Create a matrix visualization summarizing all pairwise comparisons.
% Define color levels and gray values based on style settings.
has_power_results = isfield(results, 'power_results') && ~isempty(results.power_results);
intensities = styles.win_loss.intensities;
gray_values = styles.win_loss.gray_values;
diag_color  = styles.win_loss.diag_color;
% Define edges for assigning power levels (1 to 4).
edges = [-inf, 0.25, 0.5, 0.75, inf];

% Set up the figure, tiled layout, and main title.
h_fig_winloss = figure('Name', lang.plots.figure_names.win_loss_matrix, 'Color', styles.colors.background, 'Visible', 'off');
%Tiled layout is dynamic based on num_metrics
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

%% 9. Plot: Sankey-Diagram of Rank Shifts
% Visualize the changes in ranking across the ranking stages.
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

%% 10. Plot: Final Summary (Ranking + Metrics Table)
% Create a combined figure with the final ranking plot and a results table.
h_fig_summary = figure('Name', lang.plots.figure_names.final_ranking, 'Color', styles.colors.background, 'Visible', 'off');
tcl_summary = tiledlayout(1, 2, 'Padding', 'loose', 'TileSpacing', 'compact');
title_str_summary = sprintf([lang.plots.titles.final_ranking '\n']);
title(tcl_summary, title_str_summary, 'FontSize', styles.font.title, 'FontWeight', 'bold', 'Color', styles.colors.text, 'Interpreter', 'none');

% Tile 1: Create the results table.
ax_table = nexttile;
set(ax_table, 'Visible', 'off', 'YDir', 'reverse');
ylim(ax_table, [0.5, num_datasets + 0.5]); xlim(ax_table, [0, 1]);
[~, sort_idx] = sort(results.final_rank);
sorted_names = dataset_names(sort_idx);

% Prepare all table content as strings to calculate column widths dynamically.
num_table_cols = 2 + num_metrics;
col_headers = {lang.plots.table.rank, lang.plots.table.dataset};
for i = 1:num_metrics
    col_headers{end+1} = metric_names{i};
end
table_content = cell(num_datasets, num_table_cols);

for i = 1:num_datasets
    original_idx = sort_idx(i);
    table_content{i, 1} = sprintf('%d', i);
    table_content{i, 2} = sorted_names{i};
    for m_idx = 1:num_metrics
        table_content{i, 2+m_idx} = sprintf('%.2f ± %.2f', shared_info.mean_metrics(original_idx, m_idx), shared_info.std_metrics(original_idx, m_idx));
    end
end

% Calculate the maximum width required for each column based on its content.
max_widths = zeros(1, num_table_cols);
for j = 1:num_table_cols % Check headers first.
    t = text(ax_table, 0, 0, col_headers{j}, 'Visible', 'off', 'FontWeight', 'bold', 'FontSize', ...
        styles.font.tick, 'HorizontalAlignment', 'center', 'Interpreter', 'none');
    ext = get(t, 'Extent'); max_widths(j) = max(max_widths(j), ext(3)); delete(t);
end
for i = 1:num_datasets % Then check all content rows.
    for j = 1:num_table_cols
        fontWeight = 'normal'; if j == 1, fontWeight = 'bold'; end
        t = text(ax_table, 0, 0, table_content{i, j}, 'Visible', 'off', ...
            'FontWeight', fontWeight, 'FontSize', styles.font.text, 'HorizontalAlignment', 'center', 'Interpreter', 'none');
        ext = get(t, 'Extent'); max_widths(j) = max(max_widths(j), ext(3)); delete(t);
    end
end

% Calculate column positions based on the calculated content widths.
padding = 0.04; % Relative spacing between columns.
total_width = sum(max_widths) + padding * (length(max_widths) - 1);
start_pos = (1.0 - total_width) / 2; % Center the entire table within the tile.
col_pos_vals = zeros(1, num_table_cols);
col_pos_vals(1) = start_pos + max_widths(1)/2;
for j = 2:num_table_cols, col_pos_vals(j) = col_pos_vals(j-1) + max_widths(j-1)/2 + max_widths(j)/2 + padding; end

% Place the headers in the table.
y_header_pos = 0.5;
for j = 1:num_table_cols
    text(ax_table, col_pos_vals(j), y_header_pos, col_headers{j}, 'FontWeight', 'bold', 'FontSize', styles.font.tick, ...
        'HorizontalAlignment', 'center', 'Color', styles.colors.text, 'Interpreter', 'none');
end

% Place the content for each row in the table.
for i = 1:num_datasets
    y_pos = i; % This aligns rows perfectly with the plot ticks (1, 2, 3...).
    % Rank and Name
    text(ax_table, col_pos_vals(1), y_pos, table_content{i, 1}, 'FontSize', styles.font.text, ...
        'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Color', styles.colors.text, 'Interpreter', 'none');
    text(ax_table, col_pos_vals(2), y_pos, table_content{i, 2}, 'FontSize', styles.font.text, ...
        'HorizontalAlignment', 'center', 'Color', styles.colors.text, 'Interpreter', 'none');
    % Metric data
    for m_idx = 1:num_metrics
         text(ax_table, col_pos_vals(2+m_idx), y_pos, table_content{i, 2+m_idx}, 'FontSize', styles.font.text, ...
            'HorizontalAlignment', 'center', 'Color', styles.colors.text, 'Interpreter', 'none');
    end
end

% Tile 2: Create the final rank plot.
ax_rank = nexttile;
hold(ax_rank, 'on');
set(ax_rank, 'Color', styles.colors.background);

% Prepare data for the rank plot, sorted by the final rank.
ci_ranks = quantile(results.final_bootstrap_ranks, [0.025, 0.975], 2);
[~, sort_idx] = sort(results.final_rank);
sorted_names = dataset_names(sort_idx);
sorted_ci = ci_ranks(sort_idx, :);
sorted_ranks_final = results.final_rank(sort_idx);

% Plot the 95% confidence interval error bars as a background element.
errorbar(ax_rank, sorted_ranks_final, 1:num_datasets, ...
    sorted_ranks_final - sorted_ci(:, 1), sorted_ci(:, 2) - sorted_ranks_final, ...
    styles.errorbar.final_rank{:});

% Define parameters for the frequency bubbles.
min_bubble_size = 10;
max_bubble_size = 200;
bubble_color = styles.colors.p_significant;
bubble_edge_color = styles.colors.marker_edge;
log_offset = 0.5;
scale_factor = 2.5;

% Loop over each dataset to plot the frequency bubbles representing the bootstrap rank distribution.
for i = 1:num_datasets
    y_pos = i;
    original_idx = sort_idx(i);
    rank_distribution = results.final_bootstrap_ranks(original_idx, :);
    [unique_ranks, ~, group_idx] = unique(rank_distribution);
    counts = accumarray(group_idx, 1);
    if isempty(counts), continue; end
    
    % Scale bubble sizes logarithmically to handle skewed distributions.
    scaled_sizes = min_bubble_size + ((log(counts + log_offset) * scale_factor) / log(max(counts) + log_offset)) ...
                      * (max_bubble_size - min_bubble_size);
    scatter(ax_rank, unique_ranks, repmat(y_pos, size(unique_ranks)), scaled_sizes, 'filled', ...
        'MarkerFaceColor', bubble_color, 'MarkerEdgeColor', bubble_edge_color, 'LineWidth', 0.7, ...
        'MarkerFaceAlpha', 0.8, 'HandleVisibility', 'off');
end

% Plot the final rank marker on top of the bubbles and error bars.
plot(ax_rank, sorted_ranks_final, 1:num_datasets, styles.marker.final_rank{:});
hold(ax_rank, 'off');

% Set up axes, labels, and grid for the rank plot.
set(ax_rank, 'YDir', 'reverse', 'FontSize', styles.font.tick, 'XColor', styles.colors.text, 'YColor', styles.colors.text, ...
    'YTick', 1:num_datasets, 'YTickLabel', sorted_names, 'TickLabelInterpreter', 'none');
ylabel(ax_rank, {''}, 'FontSize', styles.font.label, 'Color', styles.colors.text); % No label, names are sufficient.
xlabel(ax_rank, lang.plots.xlabels.final_rank, 'FontSize', styles.font.label, 'Color', styles.colors.text);
set(ax_rank, 'XGrid', 'on', 'YGrid', 'on', 'GridLineStyle', styles.line.grid, 'box', 'on', 'GridColor', styles.colors.grid_color);
xticks(ax_rank, 1:num_datasets);
xlim(ax_rank, [0.5, num_datasets + 0.5]); ylim(ax_rank, [0.5, num_datasets + 0.5]);

% Save the final summary graphic.
[~, fName, fExt] = fileparts(lang.files.final_ranking);
filename_summary = fullfile(output_dir, [fName, '_', ts, fExt]);
exportgraphics(h_fig_summary, filename_summary, 'Resolution', 300, 'Padding', 30);
fprintf(['  ' lang.plots.log_messages.graphic_saved '\n'], filename_summary);

%% 11. Borda-Plot (if sensitivity analysis was performed)
if isfield(results, 'borda_results') && ~isempty(results.borda_results)
    % Create the figure for the Borda consensus ranking plot.
    h_fig_borda = figure('Name', lang.plots.figure_names.borda_plot, 'Color', styles.colors.background, 'Visible', 'off');
    tcl_borda = tiledlayout(1, 1, 'Padding', 'compact');
    ax_lollipop = nexttile;
    hold(ax_lollipop, 'on');
    
    % Prepare and sort data based on the Borda rank.
    num_perms = size(results.all_permutation_ranks, 2);
    [~, sort_idx_borda] = sort(results.borda_results.rank);
    sorted_names_borda = dataset_names(sort_idx_borda);
    sorted_scores = results.borda_results.score(sort_idx_borda);
    sorted_rank_dist = results.borda_results.rank_distribution(sort_idx_borda);
    
    % Draw lollipops and rank distribution bubbles for each dataset.
    for i = 1:num_datasets
        y_pos = i;
        % Draw the lollipop stem.
        plot(ax_lollipop, [0, sorted_scores(i)], [y_pos, y_pos], '-', 'Color', [0.8 0.8 0.8], 'LineWidth', 1.5);
        % Draw the main marker for the Borda score.
        scatter(ax_lollipop, sorted_scores(i), y_pos, 100, styles.colors.p_significant, 'filled', 'MarkerEdgeColor', styles.colors.marker_edge, 'LineWidth', 1);
        % Add a text label for the score.
        text(ax_lollipop, sorted_scores(i) + 2.5, y_pos, sprintf('%.1f%%', sorted_scores(i)), 'VerticalAlignment', 'middle', ...
            'FontSize', styles.font.tick, 'Color', styles.colors.text);
        
        % Draw bubbles representing the distribution of ranks across all permutations.
        min_bubble = 40; max_bubble = 600;
        dist_data = sorted_rank_dist{i};
        unique_ranks = dist_data(:, 1);
        counts = dist_data(:, 2);
        scaled_sizes = min_bubble + (counts / num_perms) * (max_bubble - min_bubble);
        
        % Map rank positions onto the score axis for visualization.
        rank_positions_on_score_axis = ((num_datasets - unique_ranks) / (num_datasets - 1)) * sorted_scores(i);
        scatter(ax_lollipop, rank_positions_on_score_axis, repmat(y_pos, size(unique_ranks)), scaled_sizes, ...
            'filled', 'MarkerFaceColor', styles.colors.p_significant, 'MarkerEdgeColor', styles.colors.marker_edge, 'MarkerFaceAlpha', 0.5, 'LineWidth', 0.5);
    end

    % Set up axes, labels, and grid for the Borda plot.
    set(ax_lollipop, 'Color', styles.colors.background, 'YDir', 'reverse', 'FontSize', styles.font.tick, ...
        'XColor', styles.colors.text, 'YColor', styles.colors.text, ...
        'YTick', 1:num_datasets, 'YTickLabel', sorted_names_borda, 'TickLabelInterpreter', 'none');
    xlabel(ax_lollipop, lang.plots.xlabels.borda_rank, 'FontSize', styles.font.label, 'Color', styles.colors.text);
    set(ax_lollipop, 'XGrid', 'on', 'YGrid', 'on', 'GridColor', styles.colors.grid_color, 'GridLineStyle', styles.line.grid, 'box', 'on');
    xticks(ax_lollipop, 0:10:100); xlim(ax_lollipop, [-5, 110]); ylim(ax_lollipop, [0.5, num_datasets + 0.5]);
    hold(ax_lollipop, 'off');

    % Set the main title for the plot.
    title(tcl_borda, sprintf([lang.plots.titles.sensitivity_analysis '\n']), ...
        'FontSize', styles.font.title, 'FontWeight', 'bold', 'Color', styles.colors.text, 'Interpreter', 'none');
   
    % Save the Borda plot graphic.
    [~, fName, fExt] = fileparts(lang.files.borda_plot);
    filename_borda = fullfile(subfolder_ranking, [fName, '_', ts, fExt]);
    exportgraphics(h_fig_borda, filename_borda, 'Resolution', 300, 'Padding', 30);
    fprintf(['  ' lang.plots.log_messages.graphic_saved '\n'], filename_borda);
else
    h_fig_borda = gobjects(0); % Create empty handle if no plot is made
end

%% 12. Compile all graphics into PDF files
% Collect handles for the Convergence Report PDF.
h_fig_thr_global = shared_info.h_fig_thr_global;
h_fig_thr_detailed = shared_info.h_fig_thr_detailed;
h_fig_bca_global = shared_info.h_fig_bca_global;
h_fig_bca_detailed = shared_info.h_fig_bca_detailed;
h_fig_rank_conv = shared_info.h_fig_rank_conv;
convergence_report_handles = [h_fig_thr_global, h_fig_thr_detailed, h_fig_bca_global, h_fig_bca_detailed, h_fig_rank_conv];

% Collect handles for the Bootstrap Report PDF.
h_fig_hist_thr = shared_info.h_fig_hist_thr;
h_fig_hist_raw = shared_info.h_fig_hist_raw;
h_fig_hist_z0 = shared_info.h_fig_hist_z0;
h_fig_hist_a = shared_info.h_fig_hist_a;
h_fig_hist_widths = shared_info.h_fig_hist_widths;
h_fig_hist_rank = shared_info.h_fig_hist_rank;
bootstrap_report_handles = [h_fig_hist_thr, h_fig_hist_raw, h_fig_hist_z0, h_fig_hist_a, h_fig_hist_widths];

% Collect handles for the Final Report PDF.
final_report = gobjects(0);
final_report(end+1) = h_fig_summary;
if isgraphics(h_fig_borda) % Check if borda plot was created
    final_report(end+1) = h_fig_borda;
end
final_report(end+1) = h_fig_winloss;
if isgraphics(h_fig_sankey) % Check if sankey plot was created
    final_report(end+1) = h_fig_sankey;
end
final_report(end+1) = h_fig_hist_rank;

% Collect all figure handles to close them after saving.
all_handles_to_close = unique([ranking_report_handles, final_report, convergence_report_handles, bootstrap_report_handles]);

% Define a local helper function to save an array of figures to a single PDF.
save_pdf_report = @(handles, filename) ...
    arrayfun(@(i) exportgraphics(handles(i), filename, 'ContentType', 'vector', 'Append', i > 1), 1:numel(handles));

% Generate Convergence_Report.pdf.
[~, fName, fExt] = fileparts(lang.files.convergence_report);
pdf_filename_conv = fullfile(pdf_dir, [fName, '_', ts, fExt]);
if exist(pdf_filename_conv, 'file'), delete(pdf_filename_conv); end % Delete old file if it exists.
if ~isempty(convergence_report_handles)
    save_pdf_report(convergence_report_handles, pdf_filename_conv);
    fprintf(['  ' lang.plots.log_messages.pdf_saved '\n'], pdf_filename_conv);
end

% Generate Bootstrap_Report.pdf.
[~, fName, fExt] = fileparts(lang.files.bootstrap_report);
pdf_filename_boot = fullfile(pdf_dir, [fName, '_', ts, fExt]);
if exist(pdf_filename_boot, 'file'), delete(pdf_filename_boot); end
if ~isempty(bootstrap_report_handles)
    save_pdf_report(bootstrap_report_handles, pdf_filename_boot);
    fprintf(['  ' lang.plots.log_messages.pdf_saved '\n'], pdf_filename_boot);
end

% Generate Ranking_Report.pdf.
[~, fName, fExt] = fileparts(lang.files.ranking_report);
pdf_filename_rank = fullfile(pdf_dir, [fName, '_', ts, fExt]);
if exist(pdf_filename_rank, 'file'), delete(pdf_filename_rank); end
if ~isempty(ranking_report_handles)
    save_pdf_report(ranking_report_handles, pdf_filename_rank);
    fprintf(['  ' lang.plots.log_messages.pdf_saved '\n'], pdf_filename_rank);
end

% Generate Final_Report.pdf.
[~, fName, fExt] = fileparts(lang.files.final_report);
pdf_filename_final = fullfile(output_dir, [fName, '_', ts, fExt]);
if exist(pdf_filename_final, 'file'), delete(pdf_filename_final); end
if ~isempty(final_report)
    save_pdf_report(final_report, pdf_filename_final);
    fprintf(['  ' lang.plots.log_messages.pdf_saved '\n'], pdf_filename_final);
end

% Close all generated figures to free up memory.
close(all_handles_to_close);

fprintf([lang.plots.log_messages.creation_finished '\n']);
end

%% Helper Functions

function ranking_report_handles = ranking_plot(metric_idx, order, figure_name_str, ...
    y_label_str, title_base_str, filename_base_str, shared_info, styles, lang, subfolder_path, ranking_report_handles)
% Helper function to generate and save intermediate ranking plots.

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
ranking_report_handles(end+1) = h_fig_rank; % Add handle to the collection for PDF export.
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
exportgraphics(h_fig_rank, filename, 'Resolution', 300, 'Padding', 30);
fprintf(['  ' lang.plots.log_messages.graphic_saved '\n'], filename);

end

function ranking_report_handles = detail_comparison_plot(metric_idx, results, thresholds, shared_info, styles, ...
    lang, subfolder_detail_comp, ranking_report_handles)
% Helper function to generate and save detailed pairwise comparison plots for a given metric.

% Unpack necessary variables from the shared structures for easier access.
metric_names = shared_info.metric_names;
dataset_names = shared_info.dataset_names;
pair_idx_all = shared_info.pair_idx_all;
d_thresh = thresholds.d_thresh;
rel_thresh = thresholds.rel_thresh;
num_pairs = size(pair_idx_all, 1); % Total number of pairwise comparisons.

% Define pagination: set the max number of comparisons per plot page.
items_per_page = 25;
% Calculate the total number of pages needed.
num_pages = ceil(num_pairs / items_per_page);

% Prepare and sort all plot data once into a single table for efficient handling.
plot_data = table();
% Store original indices to retrieve CI data later after sorting.
plot_data.pair_idx = (1:num_pairs)';

% Extract the p-values for the current metric from the full results matrix.
% sub2ind converts the [i, j] pairs into linear indices for the matrix.
p_values = results.all_p_value_matrices{metric_idx}(sub2ind(size(results.all_p_value_matrices{metric_idx}), ...
           pair_idx_all(:,1), pair_idx_all(:,2)));
plot_data.p_values = p_values;
% Get the sort order by p-value, needed for applying Holm-Bonferroni thresholds.
[~, sort_idx_p] = sort(p_values);

% Calculate Holm-Bonferroni corrected alpha thresholds (alpha / (N-k+1)).
holm_thresholds = (shared_info.alphas(metric_idx) ./ (num_pairs - (1:num_pairs) + 1))';
% Create a vector for the thresholds, aligned with the original data pairs.
plot_data.alpha_holm = zeros(num_pairs, 1);
% Apply the sorted thresholds back to the original pair order.
plot_data.alpha_holm(sort_idx_p) = holm_thresholds;

% Extract effect size data for the current metric.
plot_data.delta = results.d_vals_all(:, metric_idx);
plot_data.abs_delta = abs(plot_data.delta); % Used for sorting.
plot_data.reldiff = results.rel_vals_all(:, metric_idx);

% Check the final significance matrix (which includes p, delta, and rel_diff criteria).
% A pair is significant if either (i,j) or (j,i) was marked as a "win".
plot_data.is_significant_p = results.all_sig_matrices{metric_idx}(sub2ind(size(results.all_sig_matrices{metric_idx}), ...
                             pair_idx_all(:,1), pair_idx_all(:,2))) | ...
                             results.all_sig_matrices{metric_idx}(sub2ind(size(results.all_sig_matrices{metric_idx}), ...
                             pair_idx_all(:,2), pair_idx_all(:,1)));
                             
% Generate Y-axis labels for each comparison (e.g., "DatasetA vs DatasetB").
y_labels = cell(num_pairs, 1);
for i = 1:num_pairs, y_labels{i} = sprintf('%s vs %s', dataset_names{pair_idx_all(i,1)}, dataset_names{pair_idx_all(i,2)}); end
plot_data.labels = y_labels;

% Sort the data table for plotting.
% This multi-level sort ensures that:
% 1. Significant "wins" are grouped at the top ('descend').
% 2. Within groups, pairs are sorted by p-value ('ascend').
% 3. Ties are broken by effect size ('descend').
plot_data_sorted = sortrows(plot_data, {'is_significant_p', 'p_values', 'abs_delta', 'reldiff'}, {'descend', 'ascend', 'descend', 'descend'});

% Loop through each page and generate a separate figure.
for page = 1:num_pages
    % Extract the data for the current page based on pagination settings.
    start_idx = (page - 1) * items_per_page + 1;
    end_idx = min(page * items_per_page, num_pairs);
    current_indices = start_idx:end_idx;
    num_current_items = numel(current_indices);

    % Get the sorted data subset for this page.
    page_data = plot_data_sorted(current_indices, :);
    % Get the corresponding CI data using the original pair indices.
    sorted_ci_d = results.ci_d_all(page_data.pair_idx, :, metric_idx);
    sorted_ci_r = results.ci_r_all(page_data.pair_idx, :, metric_idx);
    % Get the significance flag for this page's data.
    is_fully_significant = page_data.is_significant_p;
    
    % Assign bar colors: one color for non-significant, another for significant.
    colors_combined = repmat(styles.colors.p_nonsignificant, num_current_items, 1);
    colors_combined(is_fully_significant, :) = repmat(styles.colors.p_significant, sum(is_fully_significant), 1);
    
    % Dynamically calculate figure height to fit the number of items.
    base_height_px     = 280; 
    height_per_item_px = 32;  
    
    % Standard default position [left, bottom, width, height].
    default_pos = [50, 250, 1920, 1080];
    
    % Calculate the new height based on the number of items.
    new_height = base_height_px + (num_current_items * height_per_item_px);
    
    % Ensure a minimum height.
    new_height = max(new_height, 500); 
    new_height = min(new_height, default_pos(4));
    
    % The Y-position (bottom) is adjusted to keep the top edge of the figure consistent.
    current_pos = [default_pos(1), default_pos(2) + (default_pos(4) - new_height), default_pos(3), new_height];
    
    % Create the figure with the calculated position, but keep it invisible.
    h_fig = figure('Name', sprintf(lang.plots.figure_names.comparison_by_metric, metric_names{metric_idx}, page), ...
                   'Color', styles.colors.background, 'Visible', 'off', 'Position', current_pos);
    % Add the figure handle to the list for PDF export.
    ranking_report_handles(end+1) = h_fig;
    
    % Create a 1x3 tiled layout for the three subplots.
    tcl = tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
    
    % Subplot 1: p-values (log scale)
    ax1 = nexttile;
    hold(ax1, 'on');
    % Set a floor for p-values to prevent log(0) issues and improve plot readability.
    p_min_display_floor = 1e-4;
    % Add tiny random jitter to visually separate identical p-values.
    jitter = (rand(size(page_data.p_values))) * 1e-9;
    p_values_for_plot = max(page_data.p_values + jitter, p_min_display_floor);
    % Plot the horizontal bars with significance-based colors.
    barh(ax1, 1:num_current_items, p_values_for_plot, 'FaceColor', 'flat', 'CData', colors_combined, 'BarWidth', 0.7);
    % Overlay the Holm-Bonferroni threshold markers.
    plot(ax1, page_data.alpha_holm, 1:num_current_items, styles.marker.holm{:});
    hold(ax1, 'off');

    % Subplot 2: Cliff's Delta
    ax2 = nexttile;
    hold(ax2, 'on');
    % Plot the horizontal bars for the effect size.
    barh(ax2, 1:num_current_items, page_data.delta, 'FaceColor', 'flat', 'CData', colors_combined, 'BarWidth', 0.7);
    % Overlay the 95% BCa confidence intervals as error bars.
    errorbar(ax2, page_data.delta, 1:num_current_items, ...
        page_data.delta - sorted_ci_d(:,1), sorted_ci_d(:,2) - page_data.delta, styles.errorbar.detail_ci{:});
    % Draw vertical lines for the positive and negative significance thresholds.
    plot(ax2, [d_thresh(metric_idx) d_thresh(metric_idx)], [0.5, num_current_items + 0.5], styles.line.threshold{:});
    plot(ax2, [-d_thresh(metric_idx) -d_thresh(metric_idx)], [0.5, num_current_items + 0.5], styles.line.threshold{:});
    hold(ax2, 'off');

    % Subplot 3: Relative Difference
    ax3 = nexttile;
    hold(ax3, 'on');
    % Plot the horizontal bars for the relative difference.
    barh(ax3, 1:num_current_items, page_data.reldiff, 'FaceColor', 'flat', 'CData', colors_combined, 'BarWidth', 0.7);
    % Overlay the 95% BCa confidence intervals as error bars.
    errorbar(ax3, page_data.reldiff, 1:num_current_items, ...
        page_data.reldiff - sorted_ci_r(:,1), sorted_ci_r(:,2) - page_data.reldiff, styles.errorbar.detail_ci{:});
    % Draw the vertical line for the significance threshold.
    plot(ax3, [rel_thresh(metric_idx) rel_thresh(metric_idx)], [0.5, num_current_items + 0.5], styles.line.threshold{:});
    hold(ax3, 'off');

    % Configure Subplot 1 (p-values)
    % Set Y-axis to reverse (Rank 1 at top), X-axis to log scale.
    set(ax1, 'YDir', 'reverse', 'XScale', 'log', 'YTick', 1:num_current_items, 'YTickLabel', page_data.labels, ...
        'FontSize', styles.font.tick, 'XColor', styles.colors.text, 'YColor', styles.colors.text, 'TickLabelInterpreter', 'none');
    % Set Y-axis limits to tightly fit the bars.
    ylim(ax1, [0.5, num_current_items + 0.5]);
    % Define fixed ticks for the log scale for consistency.
    p_ticks_fixed = [1e-4, 1e-3, 1e-2, 1e-1, 1];
    lower_lim_fixed = 1e-5;
    xticks(ax1, p_ticks_fixed);
    xticklabels(ax1, cellstr(num2str(p_ticks_fixed(:))));
    xlim(ax1, [lower_lim_fixed, 1.1]); % Set X-limits.
    % Set labels and title from the language pack.
    xlabel(ax1, {lang.plots.xlabels.p_value_log, lang.plots.xlabels.holm_threshold}, 'FontSize', styles.font.label, ...
        'Color', styles.colors.text, 'Interpreter', 'none');
    title(ax1, sprintf([lang.plots.titles.stat_significance '\n']), 'FontSize', styles.font.label, 'Color', styles.colors.text, 'Interpreter', 'none');
    % Apply grid and box styles.
    grid(ax1, 'on'); box(ax1, 'on');
    ax1.GridLineStyle = styles.line.grid;

    % Configure Subplot 2 (Cliff's Delta)
    % Set Y-axis to reverse and apply labels.
    set(ax2, 'YDir', 'reverse', 'YTick', 1:num_current_items, 'YTickLabel', page_data.labels, 'FontSize', styles.font.tick, ...
        'XColor', styles.colors.text, 'YColor', styles.colors.text, 'TickLabelInterpreter', 'none');
    % Set Y-axis limits.
    ylim(ax2, [0.5, num_current_items + 0.5]);
    % Set X-axis limits from -1.1 to 1.1 for full range with buffer.
    xlim(ax2, [-1.1, 1.1]);
    % Set labels and title.
    xlabel(ax2, {lang.plots.xlabels.cliffs_delta_ci, lang.plots.xlabels.red_line_threshold}, 'FontSize', styles.font.label, ...
        'Color', styles.colors.text, 'Interpreter', 'none');
    title(ax2, sprintf([lang.plots.titles.effect_size_delta '\n']), 'FontSize', styles.font.label, 'Color', styles.colors.text, 'Interpreter', 'none');
    % Apply grid and box styles.
    grid(ax2, 'on'); box(ax2, 'on');
    ax2.GridLineStyle = styles.line.grid;
    
    % Configure Subplot 3 (Relative Difference)
    % Set Y-axis to reverse and apply labels.
    set(ax3, 'YDir', 'reverse', 'YTick', 1:num_current_items, 'YTickLabel', page_data.labels, 'FontSize', styles.font.tick, ...
        'XColor', styles.colors.text, 'YColor', styles.colors.text, 'TickLabelInterpreter', 'none');
    % Set Y-axis limits.
    ylim(ax3, [0.5, num_current_items + 0.5]);
    % Set X-axis limits with buffer and define clear ticks.
    xlim(ax3, [0.0, 2.1]);
    xticks(ax3, [0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2.0]); 
    % Set labels and title.
    xlabel(ax3, {lang.plots.xlabels.rel_mean_diff_ci, lang.plots.xlabels.red_line_threshold}, 'FontSize', styles.font.label, ...
        'Color', styles.colors.text, 'Interpreter', 'none');
    title(ax3, sprintf([lang.plots.titles.rel_mean_diff '\n']), 'FontSize', styles.font.label, 'Color', styles.colors.text, 'Interpreter', 'none');
    % Apply grid and box styles.
    grid(ax3, 'on'); box(ax3, 'on');
    ax3.GridLineStyle = styles.line.grid;

    % Set the main title for the entire figure.
    % Add a page indicator if there are multiple pages.
    part_str = '';
    if num_pages > 1, part_str = sprintf(lang.plots.titles.part_indicator, page, num_pages); end
    % Format the final title string.
    title_str_detail = sprintf([lang.plots.titles.detailed_comparison '\n'], metric_names{metric_idx}, part_str);
    sgtitle(tcl, title_str_detail, 'FontSize', styles.font.title, 'FontWeight', 'bold', 'Color', styles.colors.text, 'Interpreter', 'none');

    % Force MATLAB to compute the final layout before saving.
    % This is a workaround for potential rendering issues with dynamic layouts.
    drawnow;
    pause(0.1);
    set(h_fig, 'Visible', 'off');
    
    % Save the final graphic.
    filename_part = '';
    if num_pages > 1, filename_part = ['_Part' num2str(page)]; end
    base_fname = sprintf(lang.files.detail_plot, strtrim(metric_names{metric_idx}), filename_part);
    [~, fName, fExt] = fileparts(base_fname);
    filename = fullfile(subfolder_detail_comp, [fName, '_', shared_info.config.timestamp, fExt]);
    exportgraphics(h_fig, filename, 'Resolution', 300, 'Padding', 30);
    fprintf(['  ' lang.plots.log_messages.graphic_saved '\n'], filename);
end 
end