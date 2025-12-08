function passed = t09_ModeM1M3A(default_config, thresholds, n_subj, ~, ~)
% T09_MODEM1M3A - Test Case 9: Mode M1_M3A Tie-Break Logic
%
% Description:
%   Tests the Tie-Break Logic A (Single Swap limit).
%   In M1_M3A mode, if M1 is neutral (tied), M3 is used to break ties. 
%   However, the movement is limited to ONE step per dataset to prevent 
%   bubbling unstable results too far.
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
    
    title_str = 'Test 9: Mode M1_M3A (Tie-Break Logic A)';
    fprintf('\n%s\n', title_str);
    fprintf('%s\n', repmat('-', 1, strlength(title_str)));
    
    config = default_config;
    ds_names_3 = {'A', 'B', 'C'};
    
    fprintf('[Test] Verifying M3A Tie-Break Logic (using M2 data).\n');
    fprintf('[Setup] M1: Neutral (A>B>C). M2: Strong C > B > A.\n');
    
    % M1: Very close means (1.0 vs 0.98), large SD -> Neutral / No Signif. diff
    m1_data = TestHelper.generate_exact_data(n_subj, [1.0 0.98, 0.0], 0.0);
    
    % M2: Significant differences (C wins)
    m2_data = TestHelper.generate_exact_data(n_subj, [0, 10, 20], 1.0);
    
    config.ranking_mode = 'M1_M3A'; 
    
    % Input Table 
    fprintf('[Input]\n');
    h_in = {'Metric', 'Data Structure', 'Role'}; 
    d_align = {'l', 'l', 'l'};
    h_align = {'c', 'c', 'c'};
    table_data = {
        'M1', 'Means ~1.0', 'Primary (Neutral)';
        'M2', 'C(20) > B(10) > A(0)', 'Tie-Breaker (Strong)'
    };
    TestHelper.print_auto_table(h_in, table_data, d_align, h_align);
    
    [~, final_order] = TestHelper.run_single_test({m1_data, m2_data}, thresholds, config, ds_names_3);
    
    % Result Table
    fprintf('\n[Result]\n');
    h_res = {'Check', 'Final Order'}; 
    d_align = {'l', 'l'}; 
    h_align = {'c', 'c'};
    table_data = {
        'M3A Logic Applied', mat2str(final_order)
    };
    TestHelper.print_auto_table(h_res, table_data, d_align, h_align);

    % Expected: [2, 1, 3] (Indices for B, A, C) because M3A logic limits movement 
    % to ONE step per dataset to prevent bubbling unstable results.
    fprintf('\n');
    if TestHelper.check_result(final_order, [2, 1, 3], 'M1_M3A Single Tie-Break (B, A, C)')
        tests_passed = tests_passed + 1;
        fprintf('   (Validated that M3A strictly limits moves to 1 per dataset).\n');
    else
        fprintf('\n[Status] FAIL: Tie-Break logic failed.\n');
    end
    
    if tests_passed == 1
        passed = true;
    end
end
