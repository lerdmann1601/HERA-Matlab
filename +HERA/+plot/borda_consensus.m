function h_fig_borda = borda_consensus(results, shared_info, styles, lang, subfolder_ranking)
% BORDA_CONSENSUS - Visualizes the result of the sensitivity analysis as a consensus ranking.
%
% Syntax:
%   h_fig_borda = HERA.plot.borda_consensus(results, shared_info, styles, lang, subfolder_ranking)
%
% Description:
%   This function creates a Borda plot if sensitivity analysis was performed.
%   It visualizes the consensus ranking using a lollipop chart and displays
%   the distribution of ranks across all permutations as bubbles.
%   The resulting figure is saved as a high-resolution PNG.
%   Returns an empty graphics object if no sensitivity analysis results are present.
%
% Inputs:
%   results           - (struct) Contains borda_results and all_permutation_ranks.
%   shared_info       - (struct) General info (names, dataset count, timestamp).
%   styles            - (struct) Style definitions.
%   lang              - (struct) Language pack.
%   subfolder_ranking - (string) Path to save the PNG file.
%
% Outputs:
%   h_fig_borda       - (handle) Handle to the created figure or gobject(0).
%
% Author: Lukas von Erdmannsdorff

% Unpack necessary variables.
dataset_names = shared_info.dataset_names;
num_datasets = shared_info.num_datasets;
ts = shared_info.config.timestamp;

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
    exportgraphics(h_fig_borda, filename_borda, 'Resolution', 300, 'BackgroundColor', styles.colors.background, 'Padding', 30);
    fprintf(['  ' lang.plots.log_messages.graphic_saved '\n'], filename_borda);
else
    h_fig_borda = gobjects(0); % Create empty handle if no plot is made
end

end
