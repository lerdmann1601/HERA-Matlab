function results = convergence_robustness(n_sims_per_cond)
% ANALYZE_BOOTSTRAP_ROBUSTNESS - Scientific validation of bootstrap parameters.
%
% Syntax:
%   results = HERA.analysis.convergence_robustness(n_sims_per_cond)
%
% Description:
%   Performs a comprehensive Monte Carlo study across multiple data conditions
%   (Sample Sizes, Distributions) to validate HERA's default bootstrap parameters.
%
%   Study Design:
%   - Scenarios: 
%       1. N=15 (Small), Normal Dist
%       2. N=50 (Medium), Normal Dist
%       3. N=100 (Large), Normal Dist
%       4. N=30, LogNormal (Skewed) Dist
%   - For each scenario, runs 'n_sims_per_cond' simulations.
%   - Calculates "Gold Standard" (High B) reference values for every dataset.
%   - Tests "Loose", "Default", and "Strict" configurations against this reference.
%
%   Key Feature:
%   - Ranking Stability is tested using the "Gold Standard" Thresholds. 
%     This isolates the stability of the ranking bootstrap from threshold noise.
%
% Outputs:
%   results - Struct containing aggregated metrics across all scenarios.
%
% Plots:
%   Generates a summary figure showing Error Distributions and Failure Rates
%   aggregated across the entire study.
%
% Author: Lukas von Erdmannsdorff

    if nargin < 1 || isempty(n_sims_per_cond), n_sims_per_cond = 2; end

    fprintf('\n=======================================================================\n');
    fprintf('   Scientific Bootstrap Robustness Study (Sims/Cond=%d)\n', n_sims_per_cond);
    fprintf('=======================================================================\n');

    %% 1. Study Configuration
    modes = {'Weak Stability', 'Default', 'High Stability'};
    
    % Define Scenarios
    scenarios = struct();
    scenarios(1).name = 'N=15 (Normal)';   scenarios(1).N = 15; scenarios(1).Dist = 'Normal';
    scenarios(2).name = 'N=50 (Normal)';   scenarios(2).N = 50; scenarios(2).Dist = 'Normal';
    scenarios(3).name = 'N=100 (Normal)';   scenarios(3).N = 100; scenarios(2).Dist = 'Normal';
    scenarios(4).name = 'N=30 (Skewed)';   scenarios(4).N = 30; scenarios(4).Dist = 'LogNormal';
    
    % Parameter Sets (Scientific: Fixed Tolerance, Varying Stability)
    % Weak:    sm=1, st=1 (Fast, less robust)
    % Default: sm=3, st=3 (Balanced)
    % High:    sm=5, st=5 (Very robust, slower)
    
    % 1. Thresholds (Fixed Tol = 0.01)
    p_thr{1} = struct('n', 25, 'sm', 1, 'st', 1, 'tol', 0.01, 'start', 100, 'step', 150, 'end', 10000); 
    p_thr{2} = struct('n', 25, 'sm', 3, 'st', 3, 'tol', 0.01, 'start', 100, 'step', 150, 'end', 10000); % Default
    p_thr{3} = struct('n', 25, 'sm', 5, 'st', 5, 'tol', 0.01, 'start', 100, 'step', 150, 'end', 10000);
    
    % 2. BCa (Fixed Tol = 0.05)
    p_bca{1} = struct('n', 30, 'sm', 1, 'st', 1, 'tol', 0.05, 'start', 100, 'step', 200, 'end', 20000); 
    p_bca{2} = struct('n', 30, 'sm', 3, 'st', 3, 'tol', 0.05, 'start', 100, 'step', 200, 'end', 20000); % Default
    p_bca{3} = struct('n', 30, 'sm', 5, 'st', 5, 'tol', 0.05, 'start', 100, 'step', 200, 'end', 20000); 
    
    % 3. Ranking (Fixed Tol = 0.005)
    p_rank{1} = struct('n', 20, 'sm', 1, 'st', 1, 'tol', 0.005, 'start', 50, 'step', 25, 'end', 1500); 
    p_rank{2} = struct('n', 20, 'sm', 3, 'st', 3, 'tol', 0.005, 'start', 50, 'step', 25, 'end', 1500); % Default
    p_rank{3} = struct('n', 20, 'sm', 5, 'st', 5, 'tol', 0.005, 'start', 50, 'step', 25, 'end', 1500); 

    % Reference Settings (The "Truth")
    ref_B_thr = 15000;
    ref_B_bca = 30000;
    ref_B_rnk = 5000;

    %% 2. Setup Environment
    temp_dir = tempname; 
    if ~exist(temp_dir, 'dir'), mkdir(temp_dir); end
    cleanup = onCleanup(@() rmdir(temp_dir, 's'));
    
    try
        lang = HERA.get_language();
        styles = HERA.design('light', 2, true); 
    catch
        lang = struct('calculating', 'Calc', 'done', 'Done'); 
        styles = struct(); 
    end
    
    cfg_base = HERA.default();
    cfg_base.create_reports = false; % Disable reports for efficiency
    cfg_base.metric_names = {'SimMetric'};
    cfg_base.timestamp = 'RobustnessStudy';

    %% 3. Execution Loop (Scenarios -> Simulations)
    
    % Storage for aggregate results (stacking all scenarios)
    total_sims = length(scenarios) * n_sims_per_cond;
    agg_res.thr = init_storage(total_sims);
    agg_res.bca = init_storage(total_sims);
    agg_res.rnk = init_storage(total_sims);
    
    global_idx = 1;
    
    hWait = waitbar(0, 'Running Robustness Study...');
    cleanupWait = onCleanup(@() delete(hWait));
    
    for sc_idx = 1:length(scenarios)
        sc = scenarios(sc_idx);
        fprintf('   Scenario %d/%d: %s...', sc_idx, length(scenarios), sc.name);
        
        for s = 1:n_sims_per_cond
            % Check for Cancel
             if ~isempty(hWait) && isvalid(hWait)
                if getappdata(hWait, 'canceling')
                    error('User Cancelled');
                end
                waitbar(global_idx / total_sims, hWait, ...
                    sprintf('Scenario %d/%d (Sim %d/%d)', sc_idx, length(scenarios), s, n_sims_per_cond));
            end

            % Unique Seed per sim/scenario
            simStream = RandStream('mlfg6331_64', 'Seed', 10000*sc_idx + s);
            RandStream.setGlobalStream(simStream);
            
            % A) Generate Data (3 Datasets for better pair coverage)
            if strcmp(sc.Dist, 'Normal')
                d1 = 11 + 2.0 * randn(sc.N, 1);
                d2 = 12 + 2.0 * randn(sc.N, 1);
                d3 = 11.5 + 2.0 * randn(sc.N, 1); % Overlapping middle ground
            else % LogNormal (Skewed)
                d1 = exp(2 + 0.4 * randn(sc.N, 1)); % Mean ~8
                d2 = exp(2.2 + 0.4 * randn(sc.N, 1)); % Mean ~10
                d3 = exp(2.1 + 0.4 * randn(sc.N, 1)); % Intermediate
            end
            all_data = {[d1, d2, d3]};
            ds_names = {'C1', 'C2', 'C3'};
            p_idx = [1 2 3];
            eff = HERA.test.TestHelper.calculate_real_effects(all_data, 1);
            
            % B) Compute Gold Standard (Reference) for THIS data
            % Thresholds Ref
            ref_vals = struct();
            cmd_ref_t = ['[d_t_ref, r_t_ref, ~, ~, ~, ~, ~, ~] = HERA.calculate_thresholds(' ...
                 'all_data, sc.N, cfg_base, temp_dir, ref_B_thr, simStream, styles, lang);'];
            evalc(cmd_ref_t); 
            ref_vals.thr.d = d_t_ref(1);
            ref_vals.thr.r = r_t_ref(1);
            
            % BCa Ref
            cmd_ref_bca = ['[~, ci_d_ref, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = HERA.calculate_bca_ci(' ...
                   'all_data, eff.d_vals_all, eff.rel_vals_all, p_idx, sc.N, ' ...
                   'cfg_base, cfg_base.metric_names, temp_dir, temp_dir, ref_B_bca, simStream, styles, lang, ''Ref'');'];
            evalc(cmd_ref_bca);
            ref_vals.bca.width = ci_d_ref(1,2) - ci_d_ref(1,1);
            
            % Ranking Ref (Using Ref Thresholds!)
            % IMPORTANT: We use the calculated Ref thresholds here
            ref_thr_struct = struct('d_thresh', ref_vals.thr.d, 'rel_thresh', ref_vals.thr.r);
            [~, base_rank] = HERA.calculate_ranking(all_data, eff, ref_thr_struct, cfg_base, ds_names, p_idx);
            
            cmd_ref_rnk = ['[boot_r_ref, ~, ~, ~, ~] = HERA.bootstrap_ranking(' ...
                   'all_data, ref_thr_struct, cfg_base, ds_names, base_rank, p_idx, sc.N, ' ...
                   'temp_dir, temp_dir, ref_B_rnk, simStream, styles, lang, ''Ref'');'];
            evalc(cmd_ref_rnk);
            ref_vals.rnk.mean = mean(boot_r_ref(2,:));
            
            % C) Test Methods (Loose/Default/Strict)
            
            % 1. Thresholds
            for m = 1:3
                c = cfg_base; c.bootstrap_thresholds = map_params(p_thr{m});
                cmd = ['[d_t, ~, ~, ~, ~, ~, ~, conv_B] = HERA.calculate_thresholds(' ...
                     'all_data, sc.N, c, temp_dir, [], simStream, styles, lang);'];
                evalc(cmd);
                
                err_d = (d_t(1) - ref_vals.thr.d) / ref_vals.thr.d * 100;
                agg_res.thr.err(global_idx, m) = err_d;
                agg_res.thr.cost(global_idx, m) = conv_B;
                agg_res.thr.fail(global_idx, m) = (conv_B >= c.bootstrap_thresholds.B_end);
            end
            
            % 2. BCa
            for m = 1:3
                c = cfg_base; c.bootstrap_ci = map_params(p_bca{m});
                 cmd = ['[conv_B, ci_d, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = HERA.calculate_bca_ci(' ...
                       'all_data, eff.d_vals_all, eff.rel_vals_all, p_idx, sc.N, ' ...
                       'c, c.metric_names, temp_dir, temp_dir, [], simStream, styles, lang, ''Sim'');'];
                evalc(cmd);
                width_test = ci_d(1,2) - ci_d(1,1);
                err_w = (width_test - ref_vals.bca.width) / ref_vals.bca.width * 100;
                
                agg_res.bca.err(global_idx, m) = err_w;
                agg_res.bca.cost(global_idx, m) = conv_B;
                agg_res.bca.fail(global_idx, m) = (conv_B >= c.bootstrap_ci.B_end);
            end
            
            % 3. Ranking (Using REF Thresholds to isolate ranking stability)
            for m = 1:3
                c = cfg_base; c.bootstrap_ranks = map_params(p_rank{m});
                cmd = ['[boot_r, conv_B, ~, ~, ~] = HERA.bootstrap_ranking(' ...
                       'all_data, ref_thr_struct, c, ds_names, base_rank, p_idx, sc.N, ' ...
                       'temp_dir, temp_dir, [], simStream, styles, lang, ''Sim'');'];
                evalc(cmd);
                mean_rank = mean(boot_r(2, :));
                err_r = (mean_rank - ref_vals.rnk.mean) / ref_vals.rnk.mean * 100;
                
                agg_res.rnk.err(global_idx, m) = err_r;
                agg_res.rnk.cost(global_idx, m) = conv_B;
                agg_res.rnk.fail(global_idx, m) = (conv_B >= c.bootstrap_ranks.B_end);
            end
             
            global_idx = global_idx + 1;
            if mod(global_idx, 20) == 0, fprintf('.'); end
        end
        fprintf(' Done.\n');
    end
    
    %% 4. Visualization
    fprintf('All simulations completed.\nGenerating scientific report and saving plots...\n');
    create_scientific_plot(agg_res, modes, styles, total_sims, length(scenarios));
    results = agg_res;
end

% --- Helpers ---

function s = init_storage(n)
    s.err  = zeros(n, 3);
    s.cost = zeros(n, 3);
    s.fail = false(n, 3); 
end

function so = map_params(pi)
    so.n_trials = pi.n;
    so.smoothing_window = pi.sm;
    so.convergence_streak_needed = pi.st;
    so.convergence_tolerance = pi.tol;
    so.B_start = pi.start;
    so.B_step = pi.step;
    so.B_end = pi.end;
    so.min_steps_for_convergence_check = 2;
end

function create_scientific_plot(res, modes, styles, total_N, n_scenarios)
    f = figure('Name', 'Bootstrap Robustness Study', 'Color', 'w', 'Position', [100, 100, 1400, 800], 'Visible', 'on');
    cols = [0.8 0.3 0.3; 0.2 0.5 0.8; 0.3 0.7 0.4];
    
    t = tiledlayout(3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % Row 1: Thresholds
    plot_row(t, res.thr, 'Thresholds (Delta)', 'Error (%)', cols, modes);
    % Row 2: BCa
    plot_row(t, res.bca, 'BCa CI (Width)', 'Width Dev (%)', cols, modes);
    % Row 3: Ranking
    plot_row(t, res.rnk, 'Ranking (Mean)', 'Rank Dev (%)', cols, modes);
    
    title_str = sprintf('Robustness Study: Aggregated across %d Scenarios', n_scenarios);
    title(t, title_str, 'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Arial', 'Interpreter', 'none');

    % Save Result Interactively
    [file, path] = uiputfile({'*.png', 'PNG Image (*.png)'; '*.pdf', 'PDF File (*.pdf)'}, ...
                             'Save Robustness Report', 'HERA_Robustness_Report.png');
    
    if isequal(file, 0)
        fprintf('Report saving cancelled by user.\n');
    else
        out_file = fullfile(path, file);
        fprintf('Saving plot to: %s\n', out_file);
        try
            exportgraphics(f, out_file, 'Resolution', 300);
        catch
            saveas(f, out_file);
        end
    end
    fprintf('Analysis Complete.\n');
end

function plot_row(t, data, name, y_lab, cols, modes)
    % 1. Accuracy
    nexttile; hold on;
    d_plot = data.err;
    % Clip extreme outliers for plotting range
    d_plot(d_plot > 25) = 25; d_plot(d_plot < -25) = -25;
    
    % Use boxchart for better stability and aesthetic
    for i = 1:size(d_plot, 2)
       boxchart(i * ones(size(d_plot,1),1), d_plot(:,i), 'BoxFaceColor', cols(i,:), ...
           'MarkerStyle', '.', 'MarkerColor', 'k', 'BoxFaceAlpha', 0.6);
    end
    
    yline(0, '--k', 'Alpha', 0.5, 'LineWidth', 1.5);
    ylabel(y_lab, 'FontSize', 10, 'FontWeight', 'bold');
    ylim([-26 26]); 
    xticks(1:3); xticklabels(modes);
    title([name ' - Accuracy'], 'FontSize', 11);
    set(gca, 'FontSize', 10, 'LineWidth', 1.2, 'Box', 'on', 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', 'k', 'GridAlpha', 0.15);
    grid on;
    
    % 2. Efficiency
    nexttile; hold on;
    for i = 1:size(data.cost, 2)
       boxchart(i * ones(size(data.cost,1),1), data.cost(:,i), 'BoxFaceColor', cols(i,:), ...
           'MarkerStyle', '+', 'MarkerColor', 'k', 'BoxFaceAlpha', 0.6);
    end
    ylabel('Steps (B)', 'FontSize', 10, 'FontWeight', 'bold');
    xticks(1:3); xticklabels(modes);
    title([name ' - Cost'], 'FontSize', 11);
    set(gca, 'FontSize', 10, 'LineWidth', 1.2, 'Box', 'on', 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', 'k', 'GridAlpha', 0.15);
    grid on;
    
    % 3. Reliability
    nexttile;
    fail_rate = mean(data.fail, 1) * 100;
    b = bar(fail_rate, 'FaceColor', 'flat');
    b.CData = cols;
    xticklabels(modes);
    ylabel('Limit Hit (%)', 'FontSize', 10, 'FontWeight', 'bold');
    ylim([0 105]);
    title([name ' - Elbow Usage'], 'FontSize', 11);
    set(gca, 'FontSize', 10, 'LineWidth', 1.2, 'Box', 'on', 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', 'k', 'GridAlpha', 0.15);
    grid on;
    
    labels = string(round(fail_rate, 1)) + "%";
    text(b.XEndPoints, b.YEndPoints, labels, 'HorizontalAlignment','center', 'VerticalAlignment','bottom', 'FontSize', 9);
end

function ppm = ParforProgressMonitor(title, n, ~)
    ppm = struct();
end
