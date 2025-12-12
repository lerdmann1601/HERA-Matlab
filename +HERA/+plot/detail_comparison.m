function handles_list = detail_comparison(metric_idx, results, thresholds, shared_info, styles, lang, subfolder_detail_comp)
% DETAIL_COMPARISON - Generates and saves detailed pairwise comparison plots.
%
% Syntax:
%   handles_list = HERA.plot.detail_comparison(metric_idx, results, thresholds, ...
%                  shared_info, styles, lang, subfolder_detail_comp)
%
% Description:
%   This function creates detailed graphical comparisons for all dataset pairs regarding
%   a specific metric. For each pair, it visualizes:
%   1. Statistical significance (p-value with Holm-Bonferroni thresholds).
%   2. Effect size (Cliff's Delta with 95% BCa Confidence Intervals).
%   3. Relative mean difference (with 95% BCa Confidence Intervals).
%
%   Due to the potentially large number of pairs, the function supports pagination,
%   generating multiple PNG files if necessary (default: 25 comparisons per page).
%   It returns a list of figure handles for subsequent PDF report generation.
%
% Inputs:
%   metric_idx            - (int) Index of the metric to analyze (1 to N).
%   results               - (struct) Contains all calculation results (p-values, CIs, etc.).
%   thresholds            - (struct) Contains statistical thresholds (delta, rel. diff).
%   shared_info           - (struct) General info (names, global indices).
%   styles                - (struct) Style definitions for plots.
%   lang                  - (struct) Language pack.
%   subfolder_detail_comp - (string) Path to save the PNG files.
%
% Outputs:
%   handles_list          - (gobjects) Array of handles to the created figures.
%
% Author: Lukas von Erdmannsdorff

% Initialize handle array.
handles_list = gobjects(0);

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
    
    h_fig = figure('Name', sprintf(lang.plots.figure_names.comparison_by_metric, metric_names{metric_idx}, page), ...
                   'Color', styles.colors.background, 'Visible', 'off', 'Units', 'pixels', 'Position', current_pos);
    % Add the figure handle to the list for PDF export.
    handles_list(end+1) = h_fig;
    
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
    set(h_fig, 'Units', 'pixels'); % Enforce units to prevent DPI scaling collapse
    set(h_fig, 'Position', current_pos); % Re-assert position
    drawnow;
    pause(0.1);
    drawnow; % Double refresh to ensure tiledlayout expands correctly

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
