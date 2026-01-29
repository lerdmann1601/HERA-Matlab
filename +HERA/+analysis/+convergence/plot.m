function plot(command, varargin)
% PLOT - Plotting implementation for the convergence robustness analysis.
%
% Syntax:
%   HERA.analysis.convergence.plot('parameter_overview', params, scenarios, modes, refs, out_dir, ts, pdf)
%   HERA.analysis.convergence.plot('scientific_reports', res_struct, modes, styles, refs, limits, params, out_dir, base_name)
%
% Description:
%   This function serves as a dispatcher for the various plotting requirements of the robustness study.
%   It encapsulates the logic for generating the parameter overview (config table) and the 
%   scientific report plots (2x2 grid: Errors, Deviation, Cost, Failures).
%   This modular approach keeps the main controller clean.
%
% Commands:
%   'parameter_overview' - Generates a visual summary of the scenarios and parameter sets.
%   'scientific_reports' - Generates the main analysis plots for all scenarios and a global summary.
%
% Author: Lukas von Erdmannsdorff

    switch command
        case 'parameter_overview'
            plot_parameter_overview(varargin{:});
        case 'scientific_reports'
            generate_scientific_reports(varargin{:});
        otherwise
            error('Unknown plot command: %s', command);
    end
end

function generate_scientific_reports(res_struct, modes, styles, refs, limits, params, out_dir, base_name, ts_str_in)
    % Orchestrates the generation of all scientific report plots.
    colors = [0.8 0.3 0.3; 0.2 0.5 0.8; 0.3 0.7 0.4]; 
    
    if nargin >= 9 && ~isempty(ts_str_in)
        ts_str = ts_str_in;
    else
        ts_str = string(datetime('now'), 'yyyyMMdd_HHmmss');
    end

    timestamp_folder = [base_name, '_', char(ts_str)];
    
    % Support for incremental reporting mode or if path already contains timestamp
    % If ts_str was passed, we assume the caller might want to use out_dir directly if it matches
    is_incremental = contains(base_name, 'Incremental');
    has_timestamp_in_path = exist(out_dir, 'dir') && contains(out_dir, char(ts_str));
    
    if is_incremental || has_timestamp_in_path
        final_out_dir = out_dir; 
    else
        final_out_dir = fullfile(out_dir, timestamp_folder);
        if ~exist(final_out_dir, 'dir'), mkdir(final_out_dir); end
    end

    dir_graphics = fullfile(final_out_dir, 'Graphics'); 
    dirty_pdfs = fullfile(final_out_dir, 'PDF');
    if ~exist(dir_graphics, 'dir'), mkdir(dir_graphics); end
    if ~exist(dirty_pdfs, 'dir'), mkdir(dirty_pdfs); end

    pdf_full = fullfile(final_out_dir, ['Full_Combined_Report_', char(ts_str), '.pdf']);
    
    % NOTE: The parameter overview tables are created BEFORE simulation in convergence_analysis.m
    % and appended to pdf_full there. Do NOT call plot_parameter_overview here again,
    % as it would cause duplicate pages in the combined PDF.
    
    % Aggregate global results if we have multiple scenarios

    if length(res_struct) > 1
        glob.thr = stack_results(res_struct, 'thr');
        glob.bca = stack_results(res_struct, 'bca');
        glob.rnk = stack_results(res_struct, 'rnk');
        
        pdf_global = fullfile(dirty_pdfs, ['Global_Summary_', char(ts_str), '.pdf']);
        
        % Re-evaluate targets for global since pdf_global name is now known
        if is_incremental
             target_pdfs = {pdf_global};
        else
             target_pdfs = {pdf_global, pdf_full};
        end
        
        plot_single_report(glob, 'Global_Summary', modes, colors, refs, limits, dir_graphics, char(ts_str), target_pdfs);
    end
    
    % Create Per-Scenario Reports
    for i = 1:length(res_struct)
        sc = res_struct(i);
        safe_name = regexprep(sc.name, '[^a-zA-Z0-9]+', '_'); safe_name = regexprep(safe_name, '_$', ''); 
        pdf_sc = fullfile(dirty_pdfs, sprintf('%s_%s.pdf', safe_name, char(ts_str)));
        dat.thr = sc.thr; dat.bca = sc.bca; dat.rnk = sc.rnk;
        
        if is_incremental
             target_pdfs = {pdf_sc};
        else
             target_pdfs = {pdf_sc, pdf_full};
        end
        
        plot_single_report(dat, safe_name, modes, colors, refs, limits, dir_graphics, char(ts_str), target_pdfs);
    end
end

function stacked = stack_results(res_struct, field)
    % Aggregates results across multiple scenarios.
    stacked.err = []; stacked.cost = []; stacked.fail = [];
    for i = 1:length(res_struct)
        stacked.err  = [stacked.err;  res_struct(i).(field).err];
        stacked.cost = [stacked.cost; res_struct(i).(field).cost];
        stacked.fail = [stacked.fail; res_struct(i).(field).fail];
    end
end

function plot_single_report(data, suffix_name, modes, colors, refs, limits, out_path, ts_str, pdf_paths)
    % Generates a 2x2 summary tiled layout for a specific dataset/scenario.
    % Layout: Error (Boxplot), Deviation (Scatter), Cost (Boxplot), Failure (Bar).
    
    metrics = { ...
        'Thresholds (Delta)', data.thr, 'Error (%)',   'Thresholds', refs.thr, limits.thr; ...
        'BCa CI (Width)',     data.bca, 'Width Dev (%)', 'BCa',      refs.bca, limits.bca; ...
        'Ranking (Mean)',     data.rnk, 'Rank Dev (%)',  'Ranking',  refs.rnk, limits.rnk ...
    };
    if ~iscell(pdf_paths), pdf_paths = {pdf_paths}; end

    for i = 1:size(metrics, 1)
        name = metrics{i, 1}; d_val = metrics{i, 2}; label_str = metrics{i, 3};
        metric_tag= metrics{i, 4}; ref_B = metrics{i, 5}; limit_B = metrics{i, 6};
        
        f = figure('Name', ['Robustness: ' metric_tag ' - ' suffix_name], 'Color', 'w', ... 
                   'Position', [100, 100, 1200, 950], 'Visible', 'off');
        t = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'normal');
        
        clean_suffix = strrep(suffix_name, '_', ' ');
        if contains(suffix_name, 'Global_Summary')
            title_str = sprintf('%s Analysis: %s (Ref B = %d)', name, clean_suffix, ref_B);
        else
            title_str = sprintf('%s Analysis: %s', name, clean_suffix);
        end
        title(t, title_str, 'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Helvetica', 'Color', 'k');
          
        % 1. Error Distribution (Boxplot)
        nexttile; hold on;
        d_plot = d_val.err;
        d_plot(d_plot > 50) = 50; d_plot(d_plot < -50) = -50; 
        yline(0, '--k', 'LineWidth', 1.2, 'HandleVisibility', 'off');
        for m = 1:3
             boxchart(m * ones(size(d_plot,1),1), d_plot(:,m), 'BoxFaceColor', colors(m,:), ...
                 'BoxEdgeColor', 'k', 'WhiskerLineColor', 'k', 'MarkerStyle', '.', 'MarkerColor', 'k', 'BoxFaceAlpha', 0.6);
        end
        xticks(1:3); xticklabels(modes);
        ylabel(label_str, 'FontWeight', 'bold', 'Color', 'k');
        title('Error Distribution', 'Color', 'k', 'FontSize', 12);
        y_err_max = 10; 
        ylim([-y_err_max, y_err_max]);
        setup_axis();
        
        % 2. Deviation Per Simulation (Scatter)
        nexttile; hold on;
        n_points = size(d_plot, 1);
        if n_points < 50, sz = 30; alpha_val = 0.8;
        elseif n_points < 200, sz = 15; alpha_val = 0.6;
        else, sz = 8;  alpha_val = 0.4;
        end
        h_sc = gobjects(1,3);
        % Small horizontal offset to prevent overlapping points from hiding each other
        % Y-values stay unchanged to accurately represent the actual error
        jitter_offset = [-0.15, 0, 0.15];  % Relaxed left, Default center, Strict right
        for m = 1:3
            x_jittered = (1:n_points) + jitter_offset(m);
            h_sc(m) = scatter(x_jittered, d_plot(:,m), sz, colors(m,:), 'filled', ...
                'MarkerFaceAlpha', alpha_val, 'MarkerEdgeColor', 'none', 'DisplayName', modes{m});
        end
        yline(0, '--k', 'LineWidth', 1.5, 'HandleVisibility', 'off'); 
        if n_points <= 20, xticks(1:n_points); else, xticks('auto'); end
        ylabel(label_str, 'FontWeight', 'bold', 'Color', 'k');
        xlabel('Simulation Index', 'Color', 'k');
        title('Deviation per Simulation', 'Color', 'k', 'FontSize', 12);
        legend(h_sc, modes, 'Location', 'best', 'Box', 'off', 'TextColor', 'k');
        ylim([-y_err_max, y_err_max]);
        setup_axis(); xlim([0, n_points+1]);
        
        % 3. Computational Cost (Boxplot)
        nexttile; hold on;
        for m = 1:3
            boxchart(m * ones(size(d_val.cost,1),1), d_val.cost(:,m), 'BoxFaceColor', colors(m,:), ...
                 'BoxEdgeColor', 'k', 'WhiskerLineColor', 'k', 'MarkerStyle', '+', 'MarkerColor', 'k', 'BoxFaceAlpha', 0.6);
        end
        xticks(1:3); xticklabels(modes);
        ylabel('Bootstrap Steps (B)', 'FontWeight', 'bold', 'Color', 'k');
        title('Computational Cost', 'Color', 'k', 'FontSize', 12);
        ylim([0, max(d_val.cost(:)) * 1.15]);
        setup_axis();
        
        % 4. Failure Rate (Bar)
        nexttile;
        fail_rate = mean(d_val.fail, 1) * 100;
        b = bar(fail_rate, 'FaceColor', 'flat'); b.CData = colors;
        xticklabels(modes);
        ylabel('Failure Rate (%)', 'FontWeight', 'bold', 'Color', 'k');
        title(sprintf('Convergence Failures (Limit: %d)', limit_B), 'Color', 'k', 'FontSize', 12);
        ylim([0 115]); 
        text(b.XEndPoints, b.YEndPoints, string(round(fail_rate,1))+"%", 'HorizontalAlignment','center', 'VerticalAlignment','bottom', 'Color', 'k');
        setup_axis();
        
        % Save Graphic
        clean_metric = regexprep(metric_tag, '[^a-zA-Z0-9]', '');
        final_name = fullfile(out_path, sprintf('%s_%s_%s.png', suffix_name, clean_metric, ts_str));
        exportgraphics(f, final_name, 'Resolution', 300, 'BackgroundColor', 'w', 'Padding', 30);
        
        % Set consistent paper size for PDF export (prevents sizing inconsistencies)
        set(f, 'PaperUnits', 'inches');
        set(f, 'PaperSize', [12, 9.5]);  % Match figure aspect ratio (1200x950 px at 100 dpi)
        set(f, 'PaperPosition', [0, 0, 12, 9.5]);
        
        % Append to PDFs
        for k = 1:length(pdf_paths)
            try
                if exist(pdf_paths{k}, 'file')
                     exportgraphics(f, pdf_paths{k}, 'Append', true, 'ContentType', 'vector', 'BackgroundColor', 'w');
                else
                     exportgraphics(f, pdf_paths{k}, 'ContentType', 'vector', 'BackgroundColor', 'w');
                end
            catch ME
                 warning('Could not save to PDF %s: %s', pdf_paths{k}, ME.message);
            end
        end
        close(f);
    end
end

function setup_axis()
    set(gca, 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', 'Color', 'w', ... 
        'Box', 'on', 'LineWidth', 1.2, 'FontName', 'Helvetica', ...
        'GridColor', [0.7 0.7 0.7], 'GridAlpha', 0.6, 'GridLineStyle', '--');         
    grid on;
end

function plot_parameter_overview(params, scenarios, modes, refs, out_path, ts_str, pdf_path, pdf_dir)
    % Generates the overview table graphics (Scenarios and Method Parameters).
    % Uses onCleanup for robust figure handling to ensure RAM is freed.
    common_font = 'Helvetica'; header_box_h = 0.035; 
    pdf_config = ''; % Initialize for later use
    
    % Page 1: Scenarios
    fig_w = get_scenarios_width_px(scenarios);
    f1 = figure('Name', 'Study Config: Scenarios', 'Color', 'w', 'Position', [100, 100, fig_w, 680], 'Visible', 'off');
    cleanF1 = onCleanup(@() close_if_valid(f1)); % Ensure f1 is closed even on error
    
    axes('Position', [0 0 1 1], 'Visible', 'off'); hold on; xlim([0 1]); ylim([0 1]);
    
    y_curr = 0.75; 
    text(0.5, y_curr, 'Configuration: Data Scenarios', 'HorizontalAlignment', 'center', 'FontSize', 18, 'FontWeight', 'bold', 'FontName', common_font, 'Color', 'k');
    y_curr = y_curr - 0.08; 
    draw_data_section_clean(scenarios, y_curr);
    axis off;
    
    name_p1 = fullfile(out_path, sprintf('Param_Scenarios_%s.png', ts_str));
    exportgraphics(f1, name_p1, 'Resolution', 300, 'BackgroundColor', 'w', 'Padding', 30);
    
    % Set paper size for consistent PDF export
    set(f1, 'PaperUnits', 'inches');
    set(f1, 'PaperSize', [12, 9.5]);
    set(f1, 'PaperPosition', [0, 0, 12, 9.5]);
    
    % Export to combined PDF
    if exist(pdf_path, 'file'), exportgraphics(f1, pdf_path, 'Append', true, 'ContentType', 'vector', 'BackgroundColor', 'w');
    else, exportgraphics(f1, pdf_path, 'ContentType', 'vector', 'BackgroundColor', 'w'); end
    
    % Export separate PDF for Configuration tables (will be appended to by Methods)
    if nargin >= 8 && ~isempty(pdf_dir)
        pdf_config = fullfile(pdf_dir, sprintf('Config_Tables_%s.pdf', ts_str));
        exportgraphics(f1, pdf_config, 'ContentType', 'vector', 'BackgroundColor', 'w');
    end
    
    % Explicitly close f1 and clear cleanup
    close(f1); clear cleanF1;
    
    % Page 2: Methods
    f2 = figure('Name', 'Study Config: Methods', 'Color', 'w', 'Position', [150, 150, 850, 1100], 'Visible', 'off');
    cleanF2 = onCleanup(@() close_if_valid(f2)); % Ensure f2 is closed even on error
    
    axes('Position', [0 0 1 1], 'Visible', 'off'); hold on; xlim([0 1]); ylim([0 1]);
    
    y_curr = 0.75;
    text(0.5, y_curr, 'Configuration: Method Parameters', 'HorizontalAlignment', 'center', 'FontSize', 18, 'FontWeight', 'bold', 'FontName', common_font, 'Color', 'k');
    y_curr = y_curr - 0.04; 
    ref_str = sprintf('References (B):   Thresholds = %d   |   BCa = %d   |   Ranking = %d', refs.thr, refs.bca, refs.rnk);
    text(0.5, y_curr, ref_str, 'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold', 'Color', 'k', 'FontName', common_font);
    
    % Dynamic Width Calculation for Methods Table
    char_w = 0.012; pad = 0.02; min_lbl = 0.20; min_col = 0.20;
    
    % 1. Max Label Width
    lbls = {'Parameter Set', 'Trials (n)', 'Smoothing (sm)', 'Streak (st)', 'Tolerance (tol)', 'Step Size (B)', 'B Range'};
    max_lbl = 0;
    for i = 1:length(lbls), max_lbl = max(max_lbl, length(lbls{i})); end
    col_lbl_w = max(min_lbl, max_lbl * char_w + pad);
    
    % 2. Max Value Width (Modes + Params)
    max_val = 0;
    for i = 1:length(modes), max_val = max(max_val, length(modes{i})); end
    
    % Helper to check param values
    check_p = @(p) check_param_struct_len(p);
    max_val = max(max_val, check_p(params.thr));
    max_val = max(max_val, check_p(params.bca));
    max_val = max(max_val, check_p(params.rnk));
    
    col_w = max(min_col, max_val * char_w + pad);
    num_modes = length(modes);
    total_w = col_lbl_w + num_modes * col_w;
    
    % Center the table
    if total_w < 0.9
        x_start = (1 - total_w) / 2;
    else
        x_start = 0.05;
        % Squeeze if too wide
        factor = 0.9 / total_w;
        col_lbl_w = col_lbl_w * factor;
        col_w = col_w * factor;
    end
    
    y_curr = y_curr - 0.04; 
    
    rectangle('Position', [x_start, y_curr - 0.012, total_w, header_box_h], 'FaceColor', [0.95 0.95 0.95], 'EdgeColor', 'none');
    text(x_start + col_lbl_w/2, y_curr + 0.005, 'Parameter Set', 'FontWeight', 'bold', 'FontSize', 11, 'FontName', common_font, 'Color', 'k', 'HorizontalAlignment', 'center');

    for m = 1:num_modes
        text(x_start + col_lbl_w + (m-1)*col_w + col_w/2, y_curr + 0.005, modes{m}, 'FontWeight', 'bold', 'FontSize', 11, 'HorizontalAlignment', 'center', 'FontName', common_font, 'Color', 'k');
    end
    plot([x_start, x_start + total_w], [y_curr - 0.015, y_curr - 0.015], 'k-', 'LineWidth', 1.5);
    y_curr = y_curr - 0.035;  
    y_curr = draw_param_section_clean('1. Bootstrap Thresholds', params.thr, y_curr, x_start, col_lbl_w, col_w);
    y_curr = y_curr - 0.005; 
    y_curr = draw_param_section_clean('2. BCa Confidence Intervals', params.bca, y_curr, x_start, col_lbl_w, col_w);
    y_curr = y_curr - 0.005;
    y_curr = draw_param_section_clean('3. Ranking Stability', params.rnk, y_curr, x_start, col_lbl_w, col_w);
    axis off;
    
    name_p2 = fullfile(out_path, sprintf('Param_Methods_%s.png', ts_str));
    exportgraphics(f2, name_p2, 'Resolution', 300, 'BackgroundColor', 'w', 'Padding', 30);
    
    % Set paper size for consistent PDF export
    set(f2, 'PaperUnits', 'inches');
    set(f2, 'PaperSize', [12, 9.5]);
    set(f2, 'PaperPosition', [0, 0, 12, 9.5]);
    
    if exist(pdf_path, 'file'), exportgraphics(f2, pdf_path, 'Append', true, 'ContentType', 'vector', 'BackgroundColor', 'w'); end
    
    % Append Methods table to the same Config PDF
    if nargin >= 8 && ~isempty(pdf_dir) && ~isempty(pdf_config)
        exportgraphics(f2, pdf_config, 'Append', true, 'ContentType', 'vector', 'BackgroundColor', 'w');
    end
    
    % Explicitly close f2 and clear cleanup
    close(f2); clear cleanF2;
end

    function draw_data_section_clean(scenarios, y_top)
    font_name = 'Helvetica'; row_h = 0.065;
    
    % Dynamically calculate column widths based on content
    pos = get(gcf, 'Position');
    fig_w = pos(3);
    char_w_norm = 7.0 / fig_w; % 7px per char normalized
    pad_norm = 20 / fig_w;     % 20px padding (standardized)
    
    % 1. Calculate MAX widths for each column
    l_idx = length('0'); % Approx width for index
    l_name = length('Scenario Name');
    l_n    = length('N');
    l_dist = length('Distribution');
    l_sum  = length('Data Summary');
    
    for i = 1:length(scenarios)
        l_name = max(l_name, length(scenarios(i).name));
        l_n    = max(l_n,    length(sprintf('%d', scenarios(i).N)));
        l_dist = max(l_dist, length(scenarios(i).Dist));
        l_sum  = max(l_sum,  length(scenarios(i).DataSummary));
    end
    
    % Define explicit column widths (Content + Padding)
    w_idx  = 0.04; % Fixed small width for index #
    w_name = (l_name * char_w_norm) + pad_norm;
    w_n    = (l_n    * char_w_norm) + pad_norm;
    w_dist = (l_dist * char_w_norm) + pad_norm;
    w_sum  = (l_sum  * char_w_norm) + 0.005; % Minimal padding at end
    
    % 2. Calculate Total Width
    total_w = w_idx + w_name + w_n + w_dist + w_sum;
    
    % 3. Center Table
    if total_w < 0.94
        x_start = (1 - total_w) / 2;
    else
        x_start = 0.03;
        % Squeeze if necessary (should happen rarely with get_scenarios_width_px)
        factor = 0.94 / total_w;
        w_idx = w_idx * factor; w_name = w_name * factor; w_n = w_n * factor;
        w_dist = w_dist * factor; w_sum = w_sum * factor;
        total_w = 0.94;
    end
    x_end = x_start + total_w;
    
    % 4. Define Column Positions (Start of each column)
    c1 = x_start + 0.01;           % Index (with slight offset inside box)
    c2 = x_start + w_idx;          % Name start
    c3 = c2 + w_name;              % N start
    c4 = c3 + w_n;                 % Dist start
    c5 = c4 + w_dist;              % Summary start
    
    cols = {c1, c2, c3, c4, c5}; 
    
    % Draw header background
    rectangle('Position', [x_start, y_top - 0.025, total_w, 0.05], 'FaceColor', [0.95 0.95 0.95], 'EdgeColor', 'none');

    headers = {'', 'Scenario Name', 'N', 'Distribution', 'Data Summary'};
    for i = 1:length(headers)
        text(cols{i}, y_top, headers{i}, 'FontWeight', 'bold', 'FontSize', 11, 'FontName', font_name, 'Color', 'k');
    end
    y = y_top - 0.04; 
    plot([x_start, x_end], [y, y], 'k-', 'LineWidth', 1.2); 
    y = y - 0.035; 
    for i = 1:length(scenarios)
        sc = scenarios(i); y_text = y + 0.005; 
        text(cols{1}, y_text, sprintf('%d', i), 'FontSize', 10, 'FontName', font_name, 'Color', 'k');
        text(cols{2}, y_text, sc.name, 'FontSize', 10, 'FontWeight', 'bold', 'FontName', font_name, 'Interpreter', 'none', 'Color', 'k');
        text(cols{3}, y_text, sprintf('%d', sc.N), 'FontSize', 10, 'FontName', font_name, 'Color', 'k');
        text(cols{4}, y_text, sc.Dist, 'FontSize', 10, 'FontName', font_name, 'Color', 'k');
        text(cols{5}, y_text, sc.DataSummary, 'FontSize', 10, 'Color', 'k', 'FontName', font_name);
        plot([x_start, x_end], [y-0.025, y-0.025], 'Color', [0.85 0.85 0.85], 'LineWidth', 0.5);
        y = y - row_h;
    end
end

function y_next = draw_param_section_clean(title_str, p, y, x, lbl_w, val_w)
    row_h = 0.026; font_name = 'Helvetica';
    total_w = lbl_w + length(p) * val_w;
    text(x, y, title_str, 'FontWeight', 'bold', 'FontSize', 12, 'Color', 'k', 'FontName', font_name);
    y = y - 0.015; plot([x, x + total_w], [y, y], 'k-', 'LineWidth', 0.8); y = y - 0.020; 
    lbls = {'Trials (n)', 'Smoothing (sm)', 'Streak (st)', 'Tolerance (tol)', 'Step Size (B)', 'B Range'};
    fields = {'n', 'sm', 'st', 'tol', 'step', 'range'};
    for r = 1:length(lbls)
        y_text = y + 0.003; 
        text(x + 0.02, y_text, lbls{r}, 'FontSize', 10, 'Color', 'k', 'FontName', font_name);
        num_p = numel(p);
        for m = 1:num_p
            val = p{m};
            if strcmp(fields{r}, 'range')
                str = sprintf('[%d-%d]', val.start, val.end);
            else
                v = val.(fields{r});
                if abs(v) < 1 && v > 0
                     str = sprintf('%.3f', v);
                else
                     str = sprintf('%d', v);
                end
            end
            text(x + lbl_w + (m-1)*val_w + val_w/2, y_text, str, 'HorizontalAlignment', 'center', 'FontSize', 10, 'Color', 'k', 'FontName', font_name);
        end
        plot([x, x + total_w], [y-0.012, y-0.012], 'Color', [0.85 0.85 0.85], 'LineWidth', 0.5);
        y = y - row_h;
    end
    y_next = y - 0.015; 
end

function l = check_param_struct_len(p_cell)
    l = 0;
    fields = {'n', 'sm', 'st', 'tol', 'step', 'range'};
    for m = 1:length(p_cell)
        val = p_cell{m};
        for f = 1:length(fields)
             if strcmp(fields{f}, 'range')
                str = sprintf('[%d-%d]', val.start, val.end);
            else
                v = val.(fields{f});
                if abs(v) < 1 && v > 0
                     str = sprintf('%.3f', v);
                else
                     str = sprintf('%d', v);
                end
            end
            l = max(l, length(str));
        end
    end
end

function w = get_scenarios_width_px(scenarios)
    l_name = length('Scenario Name');
    l_n    = length('N');
    l_dist = length('Distribution');
    l_sum  = length('Data Summary');
    
    for i = 1:length(scenarios)
        l_name = max(l_name, length(scenarios(i).name));
        l_n    = max(l_n,    length(sprintf('%d', scenarios(i).N)));
        l_dist = max(l_dist, length(scenarios(i).Dist));
        l_sum  = max(l_sum, length(scenarios(i).DataSummary));
    end
    
    px_per_char = 7.0; % Match draw_data_section_clean
    pad_px = 30;
    margin_left = 60; % Space for row numbers and left margin
    
    % Match the draw_data_section_clean calculation:
    % c2=0.08, c3=c2+name+pad, c4=c3+n+pad, c5=c4+dist+pad, x_end=c5+sum
    % Total content width = margin + name + pad + n + pad + dist + pad + sum
    req_px = margin_left + (l_name * px_per_char) + pad_px + (l_n * px_per_char) + pad_px + (l_dist * px_per_char) + pad_px + (l_sum * px_per_char);
    
    w = max(720, ceil(req_px));
end

function close_if_valid(fig)
    % Helper function for onCleanup - safely closes a figure if it still exists
    if isvalid(fig)
        close(fig);
    end
end
