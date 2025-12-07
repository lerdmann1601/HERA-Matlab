function h_fig_hist_rank = rank_histograms(final_bootstrap_ranks, dataset_names, final_rank, selected_B, styles, lang)
% RANK_HISTOGRAMS - Generates histograms for the final rank distributions.
%
% Syntax:
%   h_fig_hist_rank = HERA.plot.rank_histograms(final_bootstrap_ranks, dataset_names, ...
%                     final_rank, selected_B, styles, lang)
%
% Description:
%   This function visualizes the distribution of ranks obtained from the bootstrap analysis
%   for each dataset. It allows for an assessment of the certainty of each rank assignment.
%
% Inputs:
%   final_bootstrap_ranks - Matrix [num_datasets x selected_B] of bootstrap ranks.
%   dataset_names         - Cell array of dataset names.
%   final_rank            - Vector of primary ranks (used for sorting the subplots).
%   selected_B            - Number of bootstrap samples used.
%   styles                - Struct containing settings for graphical design.
%   lang                  - Language pack struct for text outputs.
%
% Outputs:
%   h_fig_hist_rank       - Handle to the generated figure.
%
% Author: Lukas von Erdmannsdorff

    % Set global default font properties for the plot.
    set(0, 'DefaultAxesFontName', 'Arial');
    set(0, 'DefaultTextFontName', 'Arial');
    
    % Create the figure for the histograms.
    num_datasets_plot = size(final_bootstrap_ranks, 1);
    h_fig_hist_rank = figure('Name', lang.plots.titles.rank_dist_name, 'Color', styles.colors.background, 'Visible', 'off');
    
    % Set up the title and layout for the combined plot.
    tcl_hist = tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');
    sgtitle_str = sprintf(lang.plots.titles.rank_dist, selected_B);
    title(tcl_hist, sgtitle_str, 'FontSize', styles.font.title, 'FontWeight', 'bold', 'Color', styles.colors.text);
    
    % Get the sorting order based on the final rank to display histograms in logical order
    [~, sort_idx] = sort(final_rank);
    
    % Create a separate subplot with a rank histogram for each dataset.
    for plot_idx = 1:num_datasets_plot
        % Get the actual dataset index based on the sorted rank
        i = sort_idx(plot_idx); 
        
        ax = nexttile;
        ranks_data = final_bootstrap_ranks(i, :);
        
        % Handle the special case where only a single, unique rank value is present for a dataset.
        if isscalar(unique(ranks_data))
            % This prevents histogram errors and creates a single, narrow bar.
            val = ranks_data(1);
            bin_width = 0.25;
            edges = [val - bin_width/2, val + bin_width/2];
            histogram(ranks_data, 'BinEdges', edges, 'Normalization', 'probability', ...
                'FaceColor', styles.colors.bar_face, 'EdgeColor', styles.colors.bar_edge);
        else
            % Normal case for integer-valued ranks: use a bin width of 1.
            h = histogram(ranks_data, 'Normalization', 'probability', 'FaceColor', styles.colors.bar_face, 'EdgeColor', styles.colors.bar_edge);
            h.BinWidth = 1; 
            h.BinLimits = [min(ranks_data)-0.5, max(ranks_data)+0.5];
        end
        grid on; box on;
        set(ax, 'Color', styles.colors.background, 'GridColor', styles.colors.grid_color);
        
        % Set titles and labels for each subplot.
        title(dataset_names{i}, 'FontSize', styles.font.label, 'Color', styles.colors.text, 'Interpreter', 'none');
        xlabel(lang.plots.xlabels.rank, 'FontSize', styles.font.label, 'Color', styles.colors.text);
        ylabel(lang.plots.ylabels.rel_frequency, 'FontSize', styles.font.label, 'Color', styles.colors.text);
    
        % Set general axis properties.
        set(gca, 'FontSize', styles.font.tick, ...
                 'XTick', unique(ranks_data), ...
                 'XColor', styles.colors.text, ...
                 'YColor', styles.colors.text);
        xlim([min(ranks_data)-1, max(ranks_data)+1]);
        % Enforce the y-axis to always span from 0 to 1.
        ylim([0 1]);    
        % Force the y-axis ticks for readable standard.
        set(gca, 'YTick', 0:0.1:1);
    end
end
