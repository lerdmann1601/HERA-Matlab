function [h_fig_hist_thr, h_fig_hist_raw] = threshold_distributions(all_bootstat_d, all_bootstat_rel, d_vals_all, rel_vals_all, d_thresh, rel_thresh, rel_thresh_b, selected_B, metric_names, styles, lang, graphics_dir, config)
% THRESHOLD_DISTRIBUTIONS - Generates histograms for Bootstrap Thresholds and Raw Effect Sizes.
%
% Syntax:
%   [h_fig_hist_thr, h_fig_hist_raw] = threshold_distributions(all_bootstat_d, all_bootstat_rel, d_vals_all, ...
%       rel_vals_all, d_thresh, rel_thresh, rel_thresh_b, selected_B, metric_names, styles, lang, graphics_dir, config)
%
% Description:
%   This function creates two figures to visualize the distribution of effect sizes and thresholds:
%   1. Bootstrap Threshold Distribution: Histograms of the bootstrapped medians for both Cliff's Delta and Rel Diff.
%      Shows the calculated threshold (CI bound).
%   2. Raw Effect Size Distribution: Histograms of the raw pairwise effect sizes.
%      Overlays Kernel Density Estimate (KDE).
%   
%   Both plots are saved to the 'Threshold_Analysis' subdirectory within graphics_dir.
%
% Inputs:
%   all_bootstat_d   - Cell array of bootstrap distributions for Cliff's Delta.
%   all_bootstat_rel - Cell array of bootstrap distributions for Relative Difference.
%   d_vals_all       - Matrix of all raw Cliff's Delta values.
%   rel_vals_all     - Matrix of all raw Relative Difference values.
%   d_thresh         - Vector of final Cliff's Delta thresholds.
%   rel_thresh       - Vector of final Relative Difference thresholds (max of bootstrap and SEM).
%   rel_thresh_b     - Vector of pure bootstrap thresholds for Relative Difference.
%   selected_B       - Number of bootstrap repetitions used.
%   metric_names     - Cell array of metric names.
%   styles           - Struct with plot styling information.
%   lang             - Language package loaded from a JSON file.
%   graphics_dir     - Directory where plots should be saved.
%   config           - Configuration struct (for timestamp).
%
% Outputs:
%   h_fig_hist_thr   - Handle to the bootstrap threshold distribution figure.
%   h_fig_hist_raw   - Handle to the raw effect size distribution figure.
%
% Author: Lukas von Erdmannsdorff

    % Create a dedicated subfolder for the threshold analysis plots.
    subfolder_thresholds = fullfile(graphics_dir, 'Threshold_Analysis');
    if ~exist(subfolder_thresholds, 'dir')
        mkdir(subfolder_thresholds);
    end
    
    ts = config.timestamp;
    
    % Use the passed global styles.
    set(0, 'DefaultAxesFontName', 'Arial');
    set(0, 'DefaultTextFontName', 'Arial');
    
    num_metrics = numel(metric_names);
    effect_type_names = {'Cliff''s Delta', 'Rel Diff'};

    %% Plot 1: Histogram Distribution of Bootstrap Thresholds
    % Create the figure.
    h_fig_hist_thr = figure('Name', lang.plots.titles.threshold_dist_name, 'Color', styles.colors.background, 'Visible', 'off');
    % Dynamic tiled layout
    tcl_thr = tiledlayout(2, num_metrics, 'TileSpacing', 'compact', 'Padding', 'compact');
    sgtitle_str = sprintf([lang.plots.titles.threshold_dist '\n'], selected_B);
    title(tcl_thr, sgtitle_str, 'FontSize', styles.font.title, 'FontWeight', 'bold', 'Color', styles.colors.text);
    
    % Loop for the top row (Cliff's Delta).
    for metric_idx = 1:num_metrics
        ax = nexttile;
        bootstat_d = all_bootstat_d{metric_idx};
        data = abs(bootstat_d);
        
        hold on;
        % Handle cases with no or little data variation for robust plotting.
        if isempty(data) || isscalar(unique(data))
            bin_center = unique(data);
            if isempty(bin_center), bin_center = 0; end
            bin_width = 0.02; 
            histogram(data, 'BinEdges', [bin_center - bin_width/2, bin_center + bin_width/2], 'Normalization', 'probability', ...
                'FaceColor', styles.colors.delta_face, 'EdgeColor', styles.colors.bar_edge);
            xlim([bin_center - bin_width*5, bin_center + bin_width*5]);
            xticks(sort([bin_center, bin_center - bin_width*2, bin_center + bin_width*2]));
        else
            % Dynamically determine "nice" bin edges and tick marks for the histogram.
            points_of_interest = [data(:); d_thresh(metric_idx)]; 
            min_val = min(points_of_interest);
            max_val = max(points_of_interest);
            data_range = max_val - min_val;
            
            if data_range > 1e-6
                num_ticks_target = 5;
                raw_step = data_range / num_ticks_target;
                power = 10^floor(log10(raw_step));
                norm_step = raw_step / power;
                if norm_step < 1.5, nice_step = 1 * power;
                elseif norm_step < 3.5, nice_step = 2 * power;
                elseif norm_step < 7, nice_step = 5 * power;
                else, nice_step = 10 * power; 
                end
                nice_min = floor(min_val / nice_step) * nice_step;
                nice_max = ceil(max_val / nice_step) * nice_step;
                ticks = nice_min:nice_step:nice_max;
            else
                ticks = unique(data);
                nice_step = 0.1 * abs(ticks(1)) + 0.01;
            end
            if isempty(ticks) || numel(ticks) < 2, ticks = linspace(min_val, max_val, 3); nice_step = ticks(2)-ticks(1); end
            
            bin_edges = (ticks(1) - nice_step/2):nice_step:(ticks(end) + nice_step/2);
            
            histogram(data, 'BinEdges', bin_edges, 'Normalization', 'probability', 'FaceColor', styles.colors.delta_face, 'EdgeColor', styles.colors.bar_edge);
            
            xlim([bin_edges(1), bin_edges(end)]);
            xticks(ticks);
        end    
        % Set final y-axis limits and ticks before adding text.
        ylim([0 1]);
        set(gca, 'YTick', 0:0.1:1);
        % Add a line and text.
        xline(d_thresh(metric_idx), '--', 'Color', styles.colors.holm_threshold, 'LineWidth', 2);
        text(d_thresh(metric_idx), 1, sprintf(lang.plots.misc.boot_thr, d_thresh(metric_idx)), 'Color', styles.colors.holm_threshold, ...
            'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'Rotation', 90, 'FontSize', styles.font.tick);
        
        grid on; box on; hold off;
        set(ax, 'Color', styles.colors.background, 'GridColor', styles.colors.grid_color);
        title(sprintf('%s - %s', metric_names{metric_idx}, effect_type_names{1}), 'FontSize', styles.font.label, ...
            'Color', styles.colors.text, 'Interpreter', 'none');
        set(gca, 'FontSize', styles.font.tick, 'XColor', styles.colors.text, 'YColor', styles.colors.text);
        ylabel(lang.plots.ylabels.rel_frequency, 'FontSize', styles.font.label, 'Color', styles.colors.text); 
        xlabel(lang.plots.xlabels.median_delta, 'FontSize', styles.font.label, 'Color', styles.colors.text);
    end
    
    % Loop for the bottom row (Relative Difference).
    for metric_idx = 1:num_metrics
        ax = nexttile;
        bootstat_rel = all_bootstat_rel{metric_idx};
        data = bootstat_rel;
        
        hold on;
        if isempty(data) || isscalar(unique(data))
            bin_center = unique(data);
            if isempty(bin_center), bin_center = 0; end
            bin_width = max(abs(bin_center) * 0.1, 0.05);
            if bin_width == 0, bin_width = 0.1; end
            histogram(data, 'BinEdges', [bin_center - bin_width/2, bin_center + bin_width/2], 'Normalization', 'probability', ...
                'FaceColor', styles.colors.rel_face, 'EdgeColor', styles.colors.bar_edge);
            xticks(bin_center);
        else
            points_of_interest = [data(:); rel_thresh_b(metric_idx)]; 
            if rel_thresh(metric_idx) > rel_thresh_b(metric_idx)
                points_of_interest(end+1) = rel_thresh(metric_idx);
            end
            min_val = min(points_of_interest);
            max_val = max(points_of_interest);
            data_range = max_val - min_val;
            
            if data_range > 1e-6
                num_ticks_target = 5;
                raw_step = data_range / num_ticks_target;
                power = 10^floor(log10(raw_step));
                norm_step = raw_step / power;
                if norm_step < 1.5, nice_step = 1 * power;
                elseif norm_step < 3.5, nice_step = 2 * power;
                elseif norm_step < 7, nice_step = 5 * power;
                else, nice_step = 10 * power; 
                end
                nice_min = floor(min_val / nice_step) * nice_step;
                nice_max = ceil(max_val / nice_step) * nice_step;
                ticks = nice_min:nice_step:nice_max;
            else
                ticks = unique(data);
                nice_step = 0.1 * abs(ticks(1)) + 0.01;
            end
            if isempty(ticks) || numel(ticks) < 2, ticks = linspace(min_val, max_val, 3); nice_step = ticks(2)-ticks(1); end
            
            bin_edges = (ticks(1) - nice_step/2):nice_step:(ticks(end) + nice_step/2);
            
            histogram(data, 'BinEdges', bin_edges, 'Normalization', 'probability', 'FaceColor', styles.colors.rel_face, 'EdgeColor', styles.colors.bar_edge);
            
            xlim([bin_edges(1), bin_edges(end)]);
            xticks(ticks);
        end    
        % Set final y-axis limits and ticks before adding text.
        ylim([0 1]);
        set(gca, 'YTick', 0:0.1:1);
        % Add lines and text.
        xline(rel_thresh_b(metric_idx), '--', 'Color', styles.colors.holm_threshold, 'LineWidth', 2);
        text(rel_thresh_b(metric_idx), 1, sprintf(lang.plots.misc.boot_thr, rel_thresh_b(metric_idx)), 'Color', styles.colors.holm_threshold, ...
            'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'Rotation', 90, 'FontSize', styles.font.tick);
        if rel_thresh(metric_idx) > rel_thresh_b(metric_idx)
             xline(rel_thresh(metric_idx), '--', 'Color', styles.colors.sem_threshold, 'LineWidth', 2);
             text(rel_thresh(metric_idx), 1, sprintf(lang.plots.misc.sem_thr, rel_thresh(metric_idx)), 'Color', styles.colors.sem_threshold, ...
                 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'Rotation', 90, 'FontSize', styles.font.tick, 'FontWeight', 'bold');
        end
        
        grid on; box on; hold off;
        set(ax, 'Color', styles.colors.background, 'GridColor', styles.colors.grid_color);
        title(sprintf('%s - %s', metric_names{metric_idx}, effect_type_names{2}), ...
            'FontSize', styles.font.label, 'Color', styles.colors.text, 'Interpreter', 'none');
        set(gca, 'FontSize', styles.font.tick, 'XColor', styles.colors.text, 'YColor', styles.colors.text);
        xlabel(lang.plots.xlabels.median_rel_diff, 'FontSize', styles.font.label, 'Color', styles.colors.text);
        ylabel(lang.plots.ylabels.rel_frequency, 'FontSize', styles.font.label, 'Color', styles.colors.text); 
    end
    % Save graphic with padding.
    [~, fName, fExt] = fileparts(lang.files.dist_bootstrap_thresholds);
    filename = fullfile(subfolder_thresholds, [fName, '_', ts, fExt]);
    exportgraphics(tcl_thr, filename, 'Resolution', 300, 'Padding', 30);
    fprintf([lang.thresholds.histogram_plot_saved '\n'], filename);

    %% Plot 2: Histogram Distribution of Raw Effect Sizes
    % Creates the figure.
    h_fig_hist_raw = figure('Name', lang.plots.titles.raw_effect_dist_name, 'Color', styles.colors.background, 'Visible', 'off');
    % Dynamic tiled layout
    tcl_thr = tiledlayout(2, num_metrics, 'TileSpacing', 'compact', 'Padding', 'compact');
    sgtitle_str = sprintf([lang.plots.titles.raw_effect_dist '\n']);
    title(tcl_thr, sgtitle_str, 'FontSize', styles.font.title, 'FontWeight', 'bold', 'Color', styles.colors.text);
    
    % Loop for the top row (Cliff's Delta).
    for metric_idx = 1:num_metrics
        ax = nexttile;
        data = d_vals_all(:, metric_idx);
        data_clean = data(isfinite(data));
        unique_data = unique(data_clean);
        
        hold on;
        
        % Handle different data scenarios for robust plotting.
        if isempty(unique_data)
            % No data -> Empty plot.
        elseif isscalar(unique_data)
            % Exactly 1 value.
            bin_center = unique_data;
            bin_width = 0.1; 
            histogram(data, 'BinEdges', [bin_center - bin_width/2, bin_center + bin_width/2], 'Normalization', 'probability', ...
                'FaceColor', styles.colors.delta_face, 'EdgeColor', styles.colors.bar_edge);
            xlim([-1.1, 1.1]);
            xticks(sort(unique([-1, 0, 1, bin_center])));       
        elseif numel(unique_data) == 2
            % Exactly 2 values.
            val1 = unique_data(1);
            val2 = unique_data(2);
            bin_width = 0.2; 
            bin_edges = [val1 - bin_width/2, val1 + bin_width/2, val2 - bin_width/2, val2 + bin_width/2];
            histogram(data, 'BinEdges', bin_edges, 'Normalization', 'probability', ...
                'FaceColor', styles.colors.delta_face, 'EdgeColor', styles.colors.bar_edge);
            xlim([-1.1, 1.1]);
            xticks(sort([val1, val2, 0]));
        else
            % Dynamic ticks based on the raw data range.
            min_val = min(data_clean);
            max_val = max(data_clean);
            data_range = max_val - min_val;
            
            if data_range > 1e-6
                num_ticks_target = 5; raw_step = data_range / num_ticks_target;
                power = 10^floor(log10(raw_step)); norm_step = raw_step / power;
                if norm_step < 1.5, nice_step = 1 * power;
                elseif norm_step < 3.5, nice_step = 2 * power;
                elseif norm_step < 7, nice_step = 5 * power;
                else, nice_step = 10 * power; 
                end
                
                nice_min = floor(min_val / nice_step) * nice_step;
                nice_max = ceil(max_val / nice_step) * nice_step;
                ticks = nice_min:nice_step:nice_max;
            else, ticks = unique(data_clean); nice_step = 0.1; 
            end
            if isempty(ticks) || numel(ticks) < 2, ticks = linspace(min_val, max_val, 3); nice_step = ticks(2)-ticks(1); end
    
            bin_edges = (ticks(1) - nice_step/2):nice_step:(ticks(end) + nice_step/2);
            h_hist = histogram(data, 'BinEdges', bin_edges, 'Normalization', 'probability', ...
                'FaceColor', styles.colors.delta_face, 'EdgeColor', styles.colors.bar_edge);
            
            % Restrict the final view to the range [-1, 1].
            final_xlim_min = max(-1.1, bin_edges(1));
            final_xlim_max = min(1.1, bin_edges(end));
            
            % Adjust KDE to the final, restricted axis limits.
            evaluation_points = linspace(final_xlim_min, final_xlim_max, 200);
            [f, xi] = ksdensity(data_clean, evaluation_points, 'Bandwidth', 'normal-approx');
            % Overlay the kernel density estimate, scaled to the histogram height.
            if max(f) > 0 && ~isempty(h_hist.Values) && max(h_hist.Values) > 0
                 plot(xi, f * (max(h_hist.Values)/max(f)), 'Color', styles.colors.kde_line, 'LineWidth', 1.5);
            end      
            % Set final view.
            xlim([final_xlim_min, final_xlim_max]);
            xticks(ticks);
        end
        
        % Final axis properties.
        grid on; box on; hold off;
        set(ax, 'Color', styles.colors.background, 'GridColor', styles.colors.grid_color);
        title(sprintf('%s - %s', metric_names{metric_idx}, effect_type_names{1}), 'FontSize', styles.font.label, ...
            'Color', styles.colors.text, 'Interpreter', 'none');
        set(gca, 'FontSize', styles.font.tick, 'XColor', styles.colors.text, 'YColor', styles.colors.text);
        ylabel(lang.plots.ylabels.rel_frequency, 'FontSize', styles.font.label, 'Color', styles.colors.text); 
        xlabel('Cliff''s Delta', 'FontSize', styles.font.label, 'Color', styles.colors.text);
        ylim([0 1]);
        set(gca, 'YTick', 0:0.1:1);
    end
    
    % Loop for the bottom row (Relative Difference).
    for metric_idx = 1:num_metrics
        ax = nexttile;
        data = rel_vals_all(:, metric_idx);
        data_clean = data(isfinite(data));
        unique_data = unique(data_clean);
        
        hold on;
    
        if isempty(unique_data) || isscalar(unique_data)
            bin_center = unique(data);
            if isempty(bin_center), bin_center = 0; end
            bin_width = 0.1; 
            histogram(data, 'BinEdges', [bin_center - bin_width/2, bin_center + bin_width/2], 'Normalization', 'probability', ...
                'FaceColor', styles.colors.rel_face, 'EdgeColor', styles.colors.bar_edge);
            xlim([-0.1, 2.1]);
            xticks(sort(unique([0, 1, 2, bin_center])));
    
        elseif numel(unique_data) == 2
            val1 = unique_data(1); val2 = unique_data(2);
            bin_width = 0.1; 
            bin_edges = [val1 - bin_width/2, val1 + bin_width/2, val2 - bin_width/2, val2 + bin_width/2];
            histogram(data, 'BinEdges', bin_edges, 'Normalization', 'probability', 'FaceColor', styles.colors.rel_face, 'EdgeColor', styles.colors.bar_edge);
            xlim([-0.1, 2.1]);
            xticks(sort(unique([0, val1, val2, 2])));
        else
            % Dynamic ticks based on the raw data range.
            min_val = min(data_clean); max_val = max(data_clean);
            data_range = max_val - min_val;
            
            if data_range > 1e-6
                num_ticks_target = 5; raw_step = data_range / num_ticks_target;
                power = 10^floor(log10(raw_step)); norm_step = raw_step / power;
                if norm_step < 1.5, nice_step = 1 * power;
                elseif norm_step < 3.5, nice_step = 2 * power;
                elseif norm_step < 7, nice_step = 5 * power;
                else, nice_step = 10 * power; 
                end
                
                nice_min = floor(min_val / nice_step) * nice_step;
                nice_max = ceil(max_val / nice_step) * nice_step;
                ticks = nice_min:nice_step:nice_max;
            else, ticks = unique(data_clean); nice_step = 0.1; 
            end
            if isempty(ticks) || numel(ticks) < 2, ticks = linspace(min_val, max_val, 3); nice_step = ticks(2)-ticks(1); end
    
            bin_edges = (ticks(1) - nice_step/2):nice_step:(ticks(end) + nice_step/2);
            h_hist = histogram(data, 'BinEdges', bin_edges, 'Normalization', 'probability', ...
                'FaceColor', styles.colors.rel_face, 'EdgeColor', styles.colors.bar_edge);
            
            % Adjust KDE to the final, restricted axis limits.
            final_xlim_min = max(-0.1, bin_edges(1));
            final_xlim_max = min(2.1, bin_edges(end));
            
            evaluation_points = linspace(final_xlim_min, final_xlim_max, 200); % Uses the final limits.
            [f, xi] = ksdensity(data_clean, evaluation_points, 'Bandwidth', 'normal-approx');
            if max(f) > 0 && ~isempty(h_hist.Values) && max(h_hist.Values) > 0
                plot(xi, f * (max(h_hist.Values)/max(f)), 'Color', styles.colors.kde_line, 'LineWidth', 1.5);
            end
            
            % Set final view.
            xlim([final_xlim_min, final_xlim_max]);
            xticks(ticks);
        end
        
        % Final axis properties.
        grid on; box on; hold off;
        set(ax, 'Color', styles.colors.background, 'GridColor', styles.colors.grid_color);
        title(sprintf('%s - %s', metric_names{metric_idx}, effect_type_names{2}), ...
            'FontSize', styles.font.label, 'Color', styles.colors.text, 'Interpreter', 'none');
        set(gca, 'FontSize', styles.font.tick, 'XColor', styles.colors.text, 'YColor', styles.colors.text);
        xlabel('Relative Difference', 'FontSize', styles.font.label, 'Color', styles.colors.text);
        ylabel(lang.plots.ylabels.rel_frequency, 'FontSize', styles.font.label, 'Color', styles.colors.text); 
        ylim([0 1]);
        set(gca, 'YTick', 0:0.1:1);
    end
    
    % Save graphic with padding.
    [~, fName, fExt] = fileparts(lang.files.dist_raw_effects);
    filename = fullfile(subfolder_thresholds, [fName, '_', ts, fExt]);
    exportgraphics(tcl_thr, filename, 'Resolution', 300, 'Padding', 30);
    fprintf([lang.thresholds.effects_plot_saved '\n'], filename);

end
