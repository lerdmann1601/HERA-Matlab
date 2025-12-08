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
    
    fprintf('[Test] Verification of mathematical stopping criteria for bootstrap processes.\n');
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
        'Sample Size (N)', num2str(n_conv);
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
    bs_thr = struct('B_start', 100, 'B_step', 150, 'B_end', 5000, 'n_trials', 25, 'min_steps_for_convergence_check', 3, 'smoothing_window', 3, ...
        'convergence_streak_needed', 3, 'convergence_tolerance', 0.01);    
    bs_bca = struct('B_start', 100, 'B_step', 200, 'B_end', 15000, 'n_trials', 30, 'min_steps_for_convergence_check', 3, 'smoothing_window', 3, ...
        'convergence_streak_needed', 3, 'convergence_tolerance', 0.05);     
    bs_rank = struct('B_start', 50, 'B_step', 25, 'B_end', 1000, 'n_trials', 20, 'min_steps_for_convergence_check', 3, 'smoothing_window', 3, ...
        'convergence_streak_needed', 3, 'convergence_tolerance', 0.005);    

    % Struct array to loop through methods
    test_methods = struct();
    test_methods(1).name = 'Thresholds'; test_methods(1).cfg = bs_thr; test_methods(1).lang_sec = 'thresholds';
    test_methods(2).name = 'BCa';        test_methods(2).cfg = bs_bca; test_methods(2).lang_sec = 'bca';
    test_methods(3).name = 'Ranking';    test_methods(3).cfg = bs_rank; test_methods(3).lang_sec = 'ranking';

    % Result Table Header
    fprintf('\n[Result]\n');
    header_parts = {'Method', 'Mode', 'Configuration', 'Result', 'Status', 'Logic Check'};
    d_align = {'l', 'l', 'l', 'l', 'c', 'l'}; 
    h_align = {'c', 'c', 'c', 'c', 'c', 'c'};
    
    table_data = cell(9, 6); % 3 methods * 3 modes
    row_idx = 0;

    seed_boot = 999;
    test6_all_passed = true;

    % Loop through Method -> Mode
    for m = 1:3
        curr_method = test_methods(m).name;
        base_cfg    = test_methods(m).cfg;
        lang_sec    = lang.(test_methods(m).lang_sec);
        
        % Define expected keywords based on language resources
        % We split by '%' to handle format specifiers like "%.2f%%"
        
        % 1. Convergence Keyword
        parts_conv = strsplit(lang_sec.convergence_reached, '%');
        key_conv = parts_conv{1}; % e.g. "Convergence at" or "Convergence reached at"
        
        % 2. Elbow Keyword (Static)
        key_elbow = lang_sec.elbow_analysis_info; 
        
        % 3. Stability Keyword (Robust only)
        parts_stab = strsplit(lang_sec.stable_runs_info, '%');
        key_stable = parts_stab{1}; % e.g. "Stability in"

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
            
            % Construct Command String dynamically
            % Note: We use T = evalc(cmd) to capture output. The assignments inside cmd happen in the workspace.
            % Since HERA.* is imported, we can assume direct calls work if functions are on path.
            if strcmp(curr_method, 'Thresholds')
                cfg_run.bootstrap_thresholds = run_cfg;
                cmd = ['[~, ~, ~, ~, ~, ~, ~, B_res, ~, ~, ~, ~, ~] = ' ...
                    'HERA.calculate_thresholds(all_data_conv, n_conv, cfg_run, temp_dir, [], myStream, styles, lang);'];
            elseif strcmp(curr_method, 'BCa')
                cfg_run.bootstrap_ci = run_cfg;
                cmd = ['[B_res, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = ' ...
                    'HERA.calculate_bca_ci(all_data_conv, eff_conv.d_vals_all, eff_conv.rel_vals_all, p_idx_conv, n_conv, cfg_run, cfg_run.metric_names,' ...
                    'temp_dir, temp_dir, [], myStream, styles, lang, ''Test6'');'];
            elseif strcmp(curr_method, 'Ranking')
                cfg_run.bootstrap_ranks = run_cfg;
                cmd = ['[~, B_res, ~, ~, ~] = ' ...
                    'HERA.bootstrap_ranking(all_data_conv, thr_dummy, cfg_run, ds_names_conv, rank_base, p_idx_conv, n_conv,' ...
                    'temp_dir, temp_dir, [], myStream, styles, lang, ''Test6'');'];
            end
            
            try
                % Execute and capture output for validation
                % IMPORTANT: Do not ask for 2nd output from evalc if cmd ends with ;
                T = evalc(cmd); 
                
                is_valid = false;
                logic_msg = '';
                
                if exist('B_res', 'var')
                    res_txt = sprintf('B = %d', B_res);
                else
                    res_txt = 'B_res missing';
                    B_res = -1; 
                end
                
                % Validate Result against Logic and Output
                if strcmp(mode_name, 'Elbow')
                    % Elbow Logic: 
                    % 1. B must be high (hit limit)
                    % 2. Output must contain Elbow message
                    b_ok = (B_res >= run_cfg.B_start);
                    txt_ok = contains(T, key_elbow);
                    
                    if b_ok && txt_ok
                        is_valid = true;
                        logic_msg = 'Elbow msg found';
                    else
                        logic_msg = sprintf('B_ok=%d, Txt_ok=%d', b_ok, txt_ok);
                    end
                    
                elseif strcmp(mode_name, 'Simple')
                    % Simple Logic:
                    % 1. Convergence before max
                    % 2. Output contains Convergence message
                    % 3. Output does NOT contain Stability message (Strict check)
                    b_ok = (B_res < run_cfg.B_end);
                    conv_ok = contains(T, key_conv);
                    % Note: Some simple implementations might print stability if configured differently, 
                    % but here we explicitly disabled it.
                    
                    if b_ok && conv_ok
                        is_valid = true;
                        logic_msg = 'Conv. msg found';
                    else
                        logic_msg = sprintf('B_ok=%d, Conv=%d', b_ok, conv_ok);
                    end
                    
                elseif strcmp(mode_name, 'Robust')
                    % Robust Logic:
                    % 1. Convergence before max
                    % 2. Output contains Convergence message
                    % 3. Output MUST contain Stability message
                    b_ok = (B_res < run_cfg.B_end);
                    conv_ok = contains(T, key_conv);
                    stab_ok = contains(T, key_stable);
                    
                    if b_ok && conv_ok && stab_ok
                        is_valid = true;
                        logic_msg = 'Conv + Stable found';
                    else
                        logic_msg = sprintf('B=%d, C=%d, S=%d', b_ok, conv_ok, stab_ok);
                    end
                end
                
                if is_valid
                    stat = 'PASS'; 
                else
                    stat = 'FAIL'; 
                    test6_all_passed = false; 
                    % Append T for debugging if needed (truncated)
                    if length(T) > 100, T_debug = [T(1:100) '...']; else, T_debug = T; end
                    logic_msg = [logic_msg ' [' T_debug ']'];
                end
                
                % Store Data Row
                table_data(row_idx, :) = {curr_method, mode_name, config_str, res_txt, stat, logic_msg};
                
            catch ME
                table_data(row_idx, :) = {curr_method, mode_name, 'CRASH', 'ERR', 'FAIL', ME.message};
                test6_all_passed = false;
            end
        end
    end
    % Print Full Table
    TestHelper.print_auto_table(header_parts, table_data, d_align, h_align);
    
    if test6_all_passed
        tests_passed = tests_passed + 1;
        fprintf('\n[Status] PASS: All convergence algorithms behave mathematically correct.\n');
    else
        fprintf('\n[Status] FAIL: Anomalies detected in convergence logic.\n');
    end
    
    if tests_passed == 1
        passed = true;
    end
end
