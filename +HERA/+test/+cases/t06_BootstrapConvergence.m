function passed = t06_BootstrapConvergence(default_config, ~, ~, styles, lang)
% T06_BOOTSTRAPCONVERGENCE - Test Case 6: Bootstrap Convergence Algorithms
%
% Description:
%   Validates the mathematical stopping criteria for the 3 bootstrap processes:
%   1. Thresholds (Simple Convergence)
%   2. Confidence Intervals (BCa, Robust Convergence)
%   3. Ranking (Elbow Method/Simple)
%
%   Verifies that the implementations respect tolerances and max iterations.
%
% Inputs:
%   default_config - (struct) Base configuration
%   styles         - (struct) Plot styles
%   lang           - (struct) Language struct
%
% Returns:
%   passed - (logical) True if test passes

    import HERA.test.TestHelper
    import HERA.*
    
    passed = false;
    tests_passed = 0;
    
    title_str = 'Test 6: Bootstrap Convergence Algorithms (Global Verification)';
    fprintf('\n%s\n', title_str);
    fprintf('%s\n', repmat('-', 1, strlength(title_str)));
    
    % Hypothesis: The bootstrap process uses different convergence criteria. 
    % We need to verify that all 3 implementations (Thresholds, BCa, Ranking) 
    % correctly respect their configuration (Simple vs Robust vs Elbow).
    
    fprintf('[Test] Verification of stopping criteria for bootstrap processes.\n');
    fprintf('[Setup] Testing 3 algorithms (Simple, Robust, Elbow) across 3 analysis types with synthetic Data.\n');
    
    % 1. Setup Fixed Data for Convergence Test
    temp_dir = tempname; if ~exist(temp_dir, 'dir'), mkdir(temp_dir); end
    cleanUpTemp = onCleanup(@() rmdir(temp_dir, 's'));
    s_data = RandStream('mlfg6331_64', 'Seed', 666);
    RandStream.setGlobalStream(s_data);
    
    n_conv = 30;
    d1_conv = TestHelper.generate_exact_data(n_conv, [10, 15], 3); 
    all_data_conv = {d1_conv};
    ds_names_conv = {'C1', 'C2'};
    p_idx_conv = [1 2];
    
    % Input Table 
    fprintf('[Input]\n');
    h_in = {'Parameter', 'Value'}; 
    d_align = {'l', 'l'}; 
    h_align = {'c', 'c'}; 
    table_data = {
        'Sample Size (n)', num2str(n_conv);
        'Means', '[10, 15]';
        'SD', '3'
    };
    TestHelper.print_auto_table(h_in, table_data, d_align, h_align);
    
    % Pre-calculation needed for Bootstrap Ranking call
    eff_conv = TestHelper.calculate_real_effects(all_data_conv, 1);
    thr_dummy = struct('d_thresh', [0.15], 'rel_thresh', [0.05]); 
    config_conv = default_config;
    config_conv.metric_names = {'TestMetric'};
    [~, rank_base] = calculate_ranking(all_data_conv, eff_conv, thr_dummy, config_conv, ds_names_conv, p_idx_conv);

    % Define Base Configurations for the 3 algorithms
    bs_thr = struct('B_start', 100, 'B_step', 100, 'B_end', 5000, 'n_trials', 25, 'min_steps_for_convergence_check', 3, 'smoothing_window', 3, ...
        'convergence_streak_needed', 3, 'convergence_tolerance', 0.01);    
    bs_bca = struct('B_start', 100, 'B_step', 200, 'B_end', 15000, 'n_trials', 30, 'min_steps_for_convergence_check', 3, 'smoothing_window', 3, ...
        'convergence_streak_needed', 3, 'convergence_tolerance', 0.05);     
    bs_rank = struct('B_start', 50, 'B_step', 25, 'B_end', 1000, 'n_trials',15, 'min_steps_for_convergence_check', 3, 'smoothing_window', 3, ...
        'convergence_streak_needed', 3, 'convergence_tolerance', 0.005);    

    % Struct array to loop through methods
    test_methods = struct();
    test_methods(1).name = 'Thresholds'; test_methods(1).cfg = bs_thr; test_methods(1).lang_sec = 'thresholds';
    test_methods(2).name = 'BCa';        test_methods(2).cfg = bs_bca; test_methods(2).lang_sec = 'bca';
    test_methods(3).name = 'Ranking';    test_methods(3).cfg = bs_rank; test_methods(3).lang_sec = 'ranking';

    fprintf('\n[Bootstrap Configurations]\n');
    h_cfg = {'Algorithm', 'B_{start}', 'B_{step}', 'B_{end}', 'Trials', 'Min. Steps', 'Window', 'Streak', 'Tolerance'}; 
    d_cfg = {'l', 'c', 'c', 'c', 'c', 'c', 'c', 'c', 'c'}; 
    h_cfg_a = {'c', 'c', 'c', 'c', 'c', 'c', 'c', 'c', 'c'}; 
    t_cfg = {
        'Thresholds', num2str(bs_thr.B_start), num2str(bs_thr.B_step), num2str(bs_thr.B_end), num2str(bs_thr.n_trials), num2str(bs_thr.min_steps_for_convergence_check), num2str(bs_thr.smoothing_window), num2str(bs_thr.convergence_streak_needed), sprintf('%.1f%%', bs_thr.convergence_tolerance*100);
        'BCa CI', num2str(bs_bca.B_start), num2str(bs_bca.B_step), num2str(bs_bca.B_end), num2str(bs_bca.n_trials), num2str(bs_bca.min_steps_for_convergence_check), num2str(bs_bca.smoothing_window), num2str(bs_bca.convergence_streak_needed), sprintf('%.1f%%', bs_bca.convergence_tolerance*100);
        'Ranking', num2str(bs_rank.B_start), num2str(bs_rank.B_step), num2str(bs_rank.B_end), num2str(bs_rank.n_trials), num2str(bs_rank.min_steps_for_convergence_check), num2str(bs_rank.smoothing_window), num2str(bs_rank.convergence_streak_needed), sprintf('%.1f%%', bs_rank.convergence_tolerance*100)
    };
    TestHelper.print_auto_table(h_cfg, t_cfg, d_cfg, h_cfg_a);

    % Result Table Header
    fprintf('\n[Result]\n');
    header_parts = {'Method', 'Mode', 'Result', 'Status', 'Logic Check'};
    d_align = {'l', 'l', 'l', 'c', 'l'}; 
    h_align = {'c', 'c', 'c', 'c', 'c'};
    
    table_data = cell(9, 5); % 3 methods * 3 modes
    row_idx = 0;

    seed_boot = 999;
    test6_all_passed = true;

    % Loop through Method -> Mode
    for m = 1:3
        curr_method = test_methods(m).name;
        base_cfg    = test_methods(m).cfg;
        lang_sec    = lang.(test_methods(m).lang_sec);

        for mode_idx = 1:3
            row_idx = row_idx + 1;
            run_cfg = base_cfg;
            
            % Configure the specific Convergence Mode
            if mode_idx == 1 % Simple
                mode_name = 'Simple';
                % Disable smoothing/streak to force simple check
                run_cfg.smoothing_window = []; run_cfg.convergence_streak_needed = [];  
                config_str = sprintf('Tol=%.1f%%, (No Smooth/Streak)', run_cfg.convergence_tolerance*100);
            elseif mode_idx == 2 % Robust
                mode_name = 'Robust';
                config_str = sprintf('Tol=%.1f%%, Sm=%d, St=%d', ...
                    run_cfg.convergence_tolerance*100, run_cfg.smoothing_window, run_cfg.convergence_streak_needed);
            elseif mode_idx == 3 % Elbow
                mode_name = 'Elbow';
                % Tolerance -1 forces the algorithm to run until B_end for Elbow analysis
                run_cfg.convergence_tolerance = -1.0; 
                config_str = sprintf('Tol=-1.0, B_max=%d (Force Fail)', run_cfg.B_end);
            end
            
            % Reset RNG for fairness
            myStream = RandStream('mlfg6331_64', 'Seed', seed_boot);
            cfg_run = config_conv;
            cmd = '';
            
            % Execute Command dynamically and capture stability_data
            % Note: We use evalc simply to silence output, but we extract variables via the command logic
            try
                if strcmp(curr_method, 'Thresholds')
                    cfg_run.bootstrap_thresholds = run_cfg;
                    % Output 9 is stability_data_thr, Output 8 is selected_B
                    cmd = ['[~, ~, ~, ~, ~, ~, ~, B_res, stability_data, ~, ~, ~, ~] = ' ...
                        'HERA.calculate_thresholds(all_data_conv, n_conv, cfg_run, temp_dir, [], myStream, styles, lang);'];
                elseif strcmp(curr_method, 'BCa')
                    cfg_run.bootstrap_ci = run_cfg;
                    % Output 8 is stability_data_ci, Output 1 is selected_B
                    cmd = ['[B_res, ~, ~, ~, ~, ~, ~, stability_data, ~, ~, ~, ~, ~] = ' ...
                        'HERA.calculate_bca_ci(all_data_conv, eff_conv.d_vals_all, eff_conv.rel_vals_all, p_idx_conv, n_conv, cfg_run, cfg_run.metric_names,' ...
                        'temp_dir, temp_dir, [], myStream, styles, lang, ''Test6'');'];
                elseif strcmp(curr_method, 'Ranking')
                    cfg_run.bootstrap_ranks = run_cfg;
                    % Output 3 is stability_data_rank, Output 2 is selected_B
                    cmd = ['[~, B_res, stability_data, ~, ~] = ' ...
                        'HERA.bootstrap_ranking(all_data_conv, thr_dummy, cfg_run, ds_names_conv, rank_base, p_idx_conv, n_conv,' ...
                        'temp_dir, temp_dir, [], myStream, styles, lang, ''Test6'');'];
                end
                
                % Execute and capture output (silencing console)
                T = evalc(cmd); 
                
                is_valid = false;
                logic_msg = '';
                
                if ~exist('B_res', 'var') || ~exist('stability_data', 'var')
                    res_txt = 'Data missing';
                    B_res = -1; stability_data = struct();
                else
                    res_txt = sprintf('B = %d', B_res);
                end
                
                % Validate Result against Logic using Structured Data
                
                if strcmp(mode_name, 'Elbow')
                    % Elbow Logic: 
                    % 1. Must NOT have converged (hit max or tol=-1) -> stability_data.converged == false
                    % 2. Must have found elbow points -> ~isempty(stability_data.elbow_indices)
                    b_ok = (B_res >= run_cfg.B_start); % Should be some valid B
                    conv_state = stability_data.converged;
                    elbow_ok = isfield(stability_data, 'elbow_indices') && ~isempty(stability_data.elbow_indices);
                    
                    if b_ok && ~conv_state && elbow_ok
                        is_valid = true;
                        logic_msg = 'Elbow detected';
                    else
                        logic_msg = sprintf('B_ok=%d, Conv=%d, Elbow=%d', b_ok, conv_state, elbow_ok);
                    end
                    
                elseif strcmp(mode_name, 'Simple')
                    % Simple Logic:
                    % 1. Convergence flag must be true
                    % 2. B should be less than max (optimization worked)
                    b_ok = (B_res < run_cfg.B_end);
                    conv_state = stability_data.converged;
                    
                    if b_ok && conv_state
                        is_valid = true;
                        logic_msg = 'Converged (Simple)';
                    else
                        logic_msg = sprintf('B_ok=%d, Conv=%d', b_ok, conv_state);
                    end
                    
                elseif strcmp(mode_name, 'Robust')
                    % Robust Logic:
                    % 1. Convergence flag must be true
                    % 2. Must have used robust calculation.
                    b_ok = (B_res < run_cfg.B_end);
                    conv_state = stability_data.converged;
                    
                    if b_ok && conv_state
                        is_valid = true;
                        logic_msg = 'Converged (Robust)';
                    else
                        logic_msg = sprintf('B_ok=%d, Conv=%d', b_ok, conv_state);
                    end
                end
                
                if is_valid
                    stat = 'PASS'; 
                else
                    stat = 'FAIL'; 
                    test6_all_passed = false; 
                end
                
                % Store Data Row
                table_data(row_idx, :) = {curr_method, mode_name, res_txt, stat, logic_msg};
                
            catch ME
                table_data(row_idx, :) = {curr_method, mode_name, 'ERR', 'FAIL', ME.message};
                test6_all_passed = false;
            end
        end
    end
    % Print Full Table
    TestHelper.print_auto_table(header_parts, table_data, d_align, h_align);
    
    if test6_all_passed
        tests_passed = tests_passed + 1;
        fprintf('\n[Status] PASS: All convergence algorithms behave correct.\n');
    else
        fprintf('\n[Status] FAIL: Anomalies detected in convergence logic.\n');
    end
    
    if tests_passed == 1
        passed = true;
    end
end
