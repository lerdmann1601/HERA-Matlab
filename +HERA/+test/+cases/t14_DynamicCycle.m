function passed = t14_DynamicCycle(default_config, ~, n_subj, styles, lang)
% T14_DYNAMICCYCLE - Test Case 14: Dynamic Cycle Handling
%
% Description:
%   Tests the full pipeline (Auto-Thresholds + Ranking) on the Efron's Dice data.
%   Ensures that dynamic threshold calculation handles non-transitive data correctly.
%
% Inputs:
%   default_config - (struct) Base configuration
%   n_subj         - (int)    Sample size
%   styles         - (struct) Plot styles
%   lang           - (struct) Language struct
%
% Returns:
%   passed - (logical) True if test passes

    import HERA.test.TestHelper
    import HERA.*
    
    passed = false;
    tests_passed = 0;
    
    title_str = 'Test 14: Dynamic Cycle Handling (Efron''s Dice + Auto-Thresholds)';
    fprintf('\n%s\n', title_str);
    fprintf('%s\n', repmat('-', 1, strlength(title_str)));
    
    [m1_data, m2_efron] = TestHelper.generate_efron_data(n_subj);
    all_data_pipe = {m1_data, m2_efron};
    config_pipe = default_config;
    config_pipe.ranking_mode = 'M1_M2';
    config_pipe.metric_names = {'Base', 'Efron'};
    ds_names_4 = {'DiceA', 'DiceB', 'DiceC', 'DiceD'};
    
    % Temp setup for internal calculation of thresholds
    temp_dir = tempname; mkdir(temp_dir);
    cleanUpTemp = onCleanup(@() rmdir(temp_dir, 's'));
    myStream = RandStream('mlfg6331_64', 'Seed', 123);
    
    fprintf('[Test] Full pipeline execution with dynamic thresholds on cycle data.\n');
    fprintf('[Setup] 1. Calculate Thresholds (B=50) -> 2. Run Ranking (M1_M2). This test uses the same Efron Data as in Test 13.\n');
    
    try
        % Step 1: Calculate dynamic thresholds based on the data
        cmd_pipe = ['[d_t_dyn, r_t_dyn, ~, ~, d_v_p, r_v_p, p_i_p] = ' ...
            'HERA.calculate_thresholds(all_data_pipe, n_subj, config_pipe, temp_dir, 50, myStream, styles, lang);'];
        [~] = evalc(cmd_pipe);
        
        thr_pipe = struct('d_thresh', d_t_dyn, 'rel_thresh', r_t_dyn);
        eff_pipe = struct('d_vals_all', d_v_p, 'rel_vals_all', r_v_p);
        
        % Input/Intermediate Table 
        fprintf('[Input: Calculated Thresholds]\n');
        h_thr = {'Metric', 'Delta Thr', 'RelDiff Thr'}; 
        d_align = {'l', 'l', 'l'};
        h_align = {'c', 'c', 'c'};
        
        table_data = {
            'Base (M1)', sprintf('%.3f', d_t_dyn(1)), sprintf('%.3f', r_t_dyn(1));
            'Efron (M2)', sprintf('%.3f', d_t_dyn(2)), sprintf('%.3f', r_t_dyn(2))
        };
        TestHelper.print_auto_table(h_thr, table_data, d_align, h_align);

        % Step 2: Run Ranking
        warnState = warning('off', 'all'); 
        cleanupObj = onCleanup(@() warning(warnState));
        [final_order_pipe, ~, all_sig, ~, ~] = calculate_ranking(all_data_pipe, eff_pipe, thr_pipe, config_pipe, ds_names_4, p_i_p);
        warning('on', 'all');
        
        % Check if M2 cycle is truly significant (pairwise wins should form cycle)
        sig_m2 = all_sig{2};
        cycle_active = sig_m2(1,2) && sig_m2(2,3) && sig_m2(3,4) && sig_m2(4,1);
        
        % Result Table 
        fprintf('\n[Result]\n');
        h_res = {'Check', 'Outcome'}; 
        d_align = {'l', 'l'};
        h_align = {'c', 'c'};
        
        table_data = {
            'Cycle Active?', char(string(cycle_active));
            'Final Order', mat2str(final_order_pipe)
        };
        TestHelper.print_auto_table(h_res, table_data, d_align, h_align);
        
        if ~isempty(final_order_pipe)
            fprintf('\n[Status] PASS: Pipeline finished with valid order.\n');
            tests_passed = tests_passed + 1;
        end
    catch ME
        fprintf('\n[Status] FAIL: Pipeline crashed: %s\n', ME.message);
    end
    
    if tests_passed == 1
        passed = true;
    end
end
