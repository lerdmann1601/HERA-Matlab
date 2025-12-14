function [h_z0, h_a, h_widths] = bca_distributions(z0_d_all, a_d_all, z0_r_all, a_r_all, ci_d_all, ci_r_all, metric_names, styles, lang, graphics_dir, timestamp, B_ci)
% BCA_DISTRIBUTIONS - Generates and saves histograms for BCa parameters and CI widths.
%
% Syntax:
%   [h_z0, h_a, h_widths] = bca_distributions(z0_d_all, a_d_all, z0_r_all, a_r_all, ...
%       ci_d_all, ci_r_all, metric_names, styles, lang, graphics_dir, timestamp, B_ci)
%
% Description:
%   This function creates three main visualization figures to assess the quality of the
%   Bootstrap Confidence Intervals:
%   1. z0 Distribution: Histograms of the bias-correction factors.
%   2. a Distribution: Histograms of the acceleration (skewness) factors.
%   3. CI Width Distribution: Histograms of the final Confidence Interval widths.
%
%   Each figure uses a dynamic tiled layout corresponding to the number of metrics and
%   effect types (Cliff's Delta vs Relative Difference).
%   It includes kernel density estimation (KDE) overlays if there is sufficient data variation,
%   or robust fallback plotting for constant/low-variance data.
%
% Inputs:
%   z0_d_all, a_d_all   - Matrices of correction factors for Cliff's Delta [Pairs x Metrics].
%   z0_r_all, a_r_all   - Matrices of correction factors for Relative Diff [Pairs x Metrics].
%   ci_d_all, ci_r_all  - 3D Arrays of Confidence Intervals [Pairs x 2 x Metrics].
%   metric_names        - Cell array of strings with names of metrics.
%   styles, lang        - Structs for styling (colors/fonts) and localization.
%   graphics_dir        - Directory path where plots will be saved.
%   timestamp           - Timestamp string for unique filenames.
%   B_ci                - Optimal B value (used in plot titles).
%
% Outputs:
%   h_z0                - Handle to the z0 Distribution Figure.
%   h_a                 - Handle to the a (Acceleration) Distribution Figure.
%   h_widths            - Handle to the CI Width Distribution Figure.
%
% Author: Lukas von Erdmannsdorff

    num_metrics = numel(metric_names);
    num_pairs = size(z0_d_all, 1);
    ts = timestamp;
    
    subfolder_bca_CI = fullfile(graphics_dir, 'CI_Histograms');
    if ~exist(subfolder_bca_CI, 'dir')
        mkdir(subfolder_bca_CI);
    end

    % Set default fonts
    set(0, 'DefaultAxesFontName', 'Arial');
    set(0, 'DefaultTextFontName', 'Arial');
    effect_type_names = {'Cliff''s Delta', 'Rel Diff'};

    %% 1. Z0 Distribution (Bias Correction)
    h_z0 = figure('Name', lang.plots.titles.bca_z0_dist_name, 'Color', styles.colors.background, 'Visible', 'off');
    tcl_z0 = tiledlayout(2, num_metrics, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    sgtitle_str_z0 = sprintf(lang.plots.titles.bca_z0_dist, num_pairs);
    title(tcl_z0, sgtitle_str_z0, 'FontSize', styles.font.title, 'FontWeight', 'bold', 'Color', styles.colors.text);

    for es_type = 1:2
        for m_idx = 1:num_metrics
            ax = nexttile(tcl_z0);
            if es_type == 1
                data = z0_d_all(:, m_idx);
                face_color = styles.colors.delta_face;
            else
                data = z0_r_all(:, m_idx);
                face_color = styles.colors.rel_face;
            end
            
            % Use helper to plot histogram + density
            plot_distribution(ax, data, face_color, styles);
            
            title(sprintf('%s - %s', metric_names{m_idx}, effect_type_names{es_type}), 'FontSize', styles.font.label, ...
                'Color', styles.colors.text, 'Interpreter', 'none');
            set(gca, 'FontSize', styles.font.tick, 'XColor', styles.colors.text, 'YColor', styles.colors.text);
            xlabel(ax, lang.plots.xlabels.bias_z0, 'FontSize', styles.font.label, 'Color', styles.colors.text);
            ylabel(lang.plots.ylabels.rel_frequency, 'FontSize', styles.font.label, 'Color', styles.colors.text);
            ylim([0 1]);
            set(gca, 'YTick', 0:0.1:1);
        end
    end
    [~, fName, fExt] = fileparts(lang.files.dist_bca_bias_z0); 
    filename_z0 = fullfile(subfolder_bca_CI, [fName, '_', ts, fExt]);
    exportgraphics(h_z0, filename_z0, 'Resolution', 300, 'BackgroundColor', styles.colors.background, 'Padding', 30);
    fprintf([lang.bca.z0_histogram_saved '\n'], filename_z0);

    %% 2. A Distribution (Acceleration/Skewness)
    h_a = figure('Name', lang.plots.titles.bca_a_dist_name, 'Color', styles.colors.background, 'Visible', 'off');
    tcl_a = tiledlayout(2, num_metrics, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    sgtitle_str_a = sprintf(lang.plots.titles.bca_a_dist, num_pairs);
    title(tcl_a, sgtitle_str_a, 'FontSize', styles.font.title, 'FontWeight', 'bold', 'Color', styles.colors.text);

    for es_type = 1:2
        for m_idx = 1:num_metrics
            ax = nexttile(tcl_a);
            if es_type == 1
                data = a_d_all(:, m_idx);
                face_color = styles.colors.delta_face;
            else
                data = a_r_all(:, m_idx);
                face_color = styles.colors.rel_face;
            end
            
            plot_distribution(ax, data, face_color, styles);
            
            title(sprintf('%s - %s', metric_names{m_idx}, effect_type_names{es_type}), 'FontSize', styles.font.label, ...
                'Color', styles.colors.text, 'Interpreter', 'none');
            set(gca, 'FontSize', styles.font.tick, 'XColor', styles.colors.text, 'YColor', styles.colors.text);
            xlabel(lang.plots.xlabels.skewness_a, 'FontSize', styles.font.label, 'Color', styles.colors.text); 
            ylabel(lang.plots.ylabels.rel_frequency, 'FontSize', styles.font.label, 'Color', styles.colors.text);
            ylim([0 1]);
            set(gca, 'YTick', 0:0.1:1);
        end
    end
    [~, fName, fExt] = fileparts(lang.files.dist_bca_skew_a); 
    filename_a = fullfile(subfolder_bca_CI, [fName, '_', ts, fExt]);
    exportgraphics(h_a, filename_a, 'Resolution', 300, 'BackgroundColor', styles.colors.background, 'Padding', 30);
    fprintf([lang.bca.a_histogram_saved '\n'], filename_a);

    %% 3. CI Widths Distribution
    ci_widths_d = ci_d_all(:, 2, :) - ci_d_all(:, 1, :);
    ci_widths_r = ci_r_all(:, 2, :) - ci_r_all(:, 1, :);

    h_widths = figure('Name', lang.plots.titles.ci_width_dist_name, 'Color', styles.colors.background, 'Visible', 'off');
    tcl_widths = tiledlayout(2, num_metrics, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    sgtitle_str = sprintf(lang.plots.titles.ci_width_dist, num_pairs, B_ci);
    title(tcl_widths, sgtitle_str, 'FontSize', styles.font.title, 'FontWeight', 'bold', 'Color', styles.colors.text);

    for es_type = 1:2
        for m_idx = 1:num_metrics
            ax = nexttile(tcl_widths);
            if es_type == 1
                data = ci_widths_d(:, 1, m_idx);
                face_color = styles.colors.delta_face;
            else
                data = ci_widths_r(:, 1, m_idx);
                face_color = styles.colors.rel_face;
            end
            
            plot_distribution(ax, data, face_color, styles);
            
            title(sprintf('%s - %s', metric_names{m_idx}, effect_type_names{es_type}), 'FontSize', styles.font.label, ...
                'Color', styles.colors.text, 'Interpreter', 'none');
            set(gca, 'FontSize', styles.font.tick, 'XColor', styles.colors.text, 'YColor', styles.colors.text);
            xlabel(lang.plots.xlabels.ci_width, 'FontSize', styles.font.label, 'Color', styles.colors.text); 
            ylabel(lang.plots.ylabels.rel_frequency, 'FontSize', styles.font.label, 'Color', styles.colors.text); 
            ylim([0 1]);
            set(gca, 'YTick', 0:0.1:1);
        end
    end
    [~, fName, fExt] = fileparts(lang.files.dist_ci_widths); 
    filename_widths = fullfile(subfolder_bca_CI, [fName, '_', ts, fExt]);
    exportgraphics(h_widths, filename_widths, 'Resolution', 300, 'BackgroundColor', styles.colors.background, 'Padding', 30);
    fprintf([lang.bca.ci_width_histogram_saved '\n'], filename_widths);

end

%% Helper Function: Robust Histogram + Density Plotting
function plot_distribution(ax, data, face_color, styles)
% PLOT_DISTRIBUTION - Automatically handles data scaling, binning, and kernel density estimation.

    hold(ax, 'on');
    grid(ax, 'on'); box(ax, 'on');
    set(ax, 'Color', styles.colors.background, 'GridColor', styles.colors.grid_color);

    % Case 1: Constant or effectively single-value data
    if isscalar(unique(data))
        bin_center = unique(data);
        if isempty(bin_center), bin_center = 0; end
        bin_width = 0.02;
        
        histogram(ax, data, 'BinEdges', [bin_center - bin_width/2, bin_center + bin_width/2], ...
            'Normalization', 'probability', 'FaceColor', face_color, 'EdgeColor', styles.colors.bar_edge);
        
        xlim(ax, [bin_center - bin_width*5, bin_center + bin_width*5]);
        xticks(ax, sort([bin_center, bin_center - bin_width*2, bin_center + bin_width*2]));
        
    % Case 2: Standard Distribution
    else
        [f, xi] = ksdensity(data, 'Bandwidth', 'normal-approx');
        
        % Nicer axis ticks (Auto-Scale)
        min_val = min([data(:); xi(:)]);
        max_val = max([data(:); xi(:)]);
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
        if isempty(ticks), ticks = linspace(min_val, max_val, 3); nice_step = ticks(2)-ticks(1); end
        if numel(ticks) < 2, ticks = linspace(ticks(1)-nice_step, ticks(1)+nice_step, 3); end
        
        % Ensure bins are centered on ticks
        bin_edges = (ticks(1) - nice_step/2):nice_step:(ticks(end) + nice_step/2);
        
        h_hist = histogram(ax, data, 'BinEdges', bin_edges, 'Normalization', 'probability', ...
            'FaceColor', face_color, 'EdgeColor', styles.colors.bar_edge);
        
        % Overlay Density Curve scaled to Histogram
        if max(f) > 0
            plot(ax, xi, f * (max(h_hist.Values)/max(f)), 'Color', styles.colors.kde_line, 'LineWidth', 1.5);
        end
        xlim(ax, [bin_edges(1), bin_edges(end)]);
        xticks(ax, ticks);
    end
    hold(ax, 'off');
end
