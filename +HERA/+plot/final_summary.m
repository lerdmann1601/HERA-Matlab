function h_fig_summary = final_summary(results, shared_info, styles, lang, output_dir)
% FINAL_SUMMARY - Creates a combined figure with the final ranking plot and a results table.
%
% Syntax:
%   h_fig_summary = HERA.plot.final_summary(results, shared_info, styles, lang, output_dir)
%
% Description:
%   This function generates the final summary graphic which consists of two panels:
%   1. A data table showing the ranks, dataset names, and mean/sd for each metric.
%   2. A ranking plot summarizing the distribution of bootstrap ranks (bubbles),
%      the 95% confidence intervals (error bars), and the final rank.
%   The resulting figure is saved as a high-resolution PNG.
%
% Inputs:
%   results     - (struct) Contains final ranks, bootstrap distributions, etc.
%   shared_info - (struct) General info (names, metrics, timestamp).
%   styles      - (struct) Style definitions.
%   lang        - (struct) Language pack.
%   output_dir  - (string) Path to save the PNG file.
%
% Outputs:
%   h_fig_summary - (handle) Handle to the created figure object.
%
% Author: Lukas von Erdmannsdorff

% Unpack necessary variables.
metric_names = shared_info.metric_names;
dataset_names = shared_info.dataset_names;
num_datasets = shared_info.num_datasets;
num_metrics = numel(metric_names);
ts = shared_info.config.timestamp;

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
exportgraphics(h_fig_summary, filename_summary, 'Resolution', 300, 'BackgroundColor', styles.colors.background, 'Padding', 30);
fprintf(['  ' lang.plots.log_messages.graphic_saved '\n'], filename_summary);

end
