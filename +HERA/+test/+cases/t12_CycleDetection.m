function passed = t12_CycleDetection(default_config, thresholds, n_subj, ~, ~)
% T12_CYCLEDETECTION - Test Case 12: Cycle Detection Mechanism
%
% Description:
%   Tests the ability of the algorithm to detect and break infinite loops 
%   (Rock-Paper-Scissors scenarios) in the ranking logic.
%   Forces an artificial cycle via effect size injection.
%
% Inputs:
%   default_config - (struct) Base configuration
%   thresholds     - (struct) Base thresholds
%   n_subj         - (int)    Sample size
%
% Returns:
%   passed - (logical) True if test passes

    import HERA.test.TestHelper
    import HERA.*
    
    passed = false;
    tests_passed = 0;
    
    title_str = 'Test 12: Cycle Detection (Artificial Injection)';
    fprintf('\n%s\n', title_str);
    fprintf('%s\n', repmat('-', 1, strlength(title_str)));
    config = default_config;
    ds_names = {'D1', 'D2', 'D3'};
    
    m1_data = TestHelper.generate_exact_data(n_subj, [30, 20, 10], 1.0); % Base Order [1, 2, 3]
    m2_data = TestHelper.generate_exact_data(n_subj, [0, 0, 0], 1.0);
    
    % Inject a Logical Cycle in M2 by forcing effect sizes:
    % 2 beats 1, 3 beats 2, 1 beats 3 (Rock-Paper-Scissors)
    eff = TestHelper.calculate_real_effects({m1_data, m2_data}, 2);
    eff.d_vals_all(:, 2) = [-0.9; 0.9; -0.9]; 
    
    fprintf('[Test] Cycle Detection and Fallback Logic.\n');
    fprintf('[Setup] M1: Transitive. M2: Artificial Cycle injected via Cliff''s Delta.\n');
    
    % Input Table 
    fprintf('[Input: M2 Injection]\n');
    h_cyc = {'Pair', 'Delta', 'Direction'}; 
    d_align = {'l', 'c', 'l'}; 
    h_align = {'c', 'c', 'c'};
    table_data = {
        '2 vs 1', '-0.9', '2 > 1';
        '3 vs 2', '-0.9', '3 > 2';
        '1 vs 3', '+0.9', '1 > 3'
    };
    TestHelper.print_auto_table(h_cyc, table_data, d_align, h_align);
    
    config.ranking_mode = 'M1_M2';
    
    % Suppress warning for expected cycle detection
    warnState = warning('off', 'all'); 
    cleanupObj = onCleanup(@() warning(warnState));
    [final_order, ~, ~, ~, ~] = calculate_ranking({m1_data, m2_data}, eff, thresholds, config, ds_names, nchoosek(1:3, 2));
    warning('on', 'all');
    
    % Result
    fprintf('\n[Result]\n');
    h_res = {'Check', 'Outcome'}; 
    d_align = {'l', 'c'}; 
    h_align = {'c', 'c'}; 
    
    % The algorithm should detect the infinite loop in swaps and revert to M1 order
    cycle_handled = isequal(final_order(:)', [1, 2, 3]);

    table_data = {
        'Cycle Detected?', char(string(cycle_handled));
        'Final Order', mat2str(final_order(:)')
    };
    TestHelper.print_auto_table(h_res, table_data, d_align, h_align);

    % Assertions
    fprintf('\n')
    if TestHelper.check_result(final_order, [1, 2, 3], 'Cycle Handled (Revert to M1)')
        tests_passed = tests_passed + 1;
        fprintf('[Status] PASS: Infinite loop detected, M2 swaps aborted.\n');
    else
        fprintf('[Status] FAIL: Cycle detection failed.\n');
    end
    
    if tests_passed == 1
        passed = true;
    end
end
