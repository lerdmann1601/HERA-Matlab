function passed = t10_ModeM1M2M3(default_config, thresholds, n_subj, ~, ~)
% T10_MODEM1M2M3 - Test Case 10: Mode M1_M2_M3 Fallback Logic
%
% Description:
%   Tests Tie-Break Logic 3B.
%   Hypothesis: If M1 AND M2 are neutral (tied), sort strictly by M3 (Gradient).
%   This is the fallback "tie-breaker of last resort".
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
    
    title_str = 'Test 10: Mode M1_M2_M3 (Logic 3B Fallback)';
    fprintf('\n%s\n', title_str);
    fprintf('%s\n', repmat('-', 1, strlength(title_str)));
    
    ds_names = {'D1', 'D2', 'D3'};
    config = default_config;
    
    % Logic 3B: If M1 AND M2 are neutral (tied), sort by M3.
    % Ensure Clean RNG for this sensitive test to prevent noise flipping D1/D2
    s_local = RandStream('mlfg6331_64', 'Seed', 999);
    s_global = RandStream.setGlobalStream(s_local);

    fprintf('[Test] Verifying Logic 3B (M3 Sorting when M1/M2 tied).\n');
    fprintf('[Setup] M1/M2 Neutral (~10). M3 Strong D3(20) > D1/D2(10).\n');
    
    % M1/M2: Small differences, statistically neutral
    m1_data = TestHelper.generate_exact_data(n_subj, [10.4, 10.2, 10.0], 1.0);
    m2_data = TestHelper.generate_exact_data(n_subj, [10.4, 10.2, 10.0], 1.0);
    
    % M3: D3 (20) is significantly better than D1/D2
    m3_data = TestHelper.generate_exact_data(n_subj, [10.0, 10.0, 20.0], 1.0);
    
    % Input Table for Means 
    fprintf('[Input]\n');
    h_m = {'Dataset', 'M1 Mean', 'M2 Mean', 'M3 Mean'}; 
    d_align = {'c', 'c', 'c', 'c'}; 
    h_align = {'c', 'c', 'c', 'c'};
    table_data = {
        'D1', sprintf('%.2f', mean(m1_data(:,1))), sprintf('%.2f', mean(m2_data(:,1))), sprintf('%.2f', mean(m3_data(:,1)));
        'D2', sprintf('%.2f', mean(m1_data(:,2))), sprintf('%.2f', mean(m2_data(:,2))), sprintf('%.2f', mean(m3_data(:,2)));
        'D3', sprintf('%.2f', mean(m1_data(:,3))), sprintf('%.2f', mean(m2_data(:,3))), sprintf('%.2f', mean(m3_data(:,3)))
    };
    TestHelper.print_auto_table(h_m, table_data, d_align, h_align);
    
    config.ranking_mode = 'M1_M2_M3';
    [~, final_order] = TestHelper.run_single_test({m1_data, m2_data, m3_data}, thresholds, config, ds_names);
    
    % Result Table
    fprintf('\n[Result]\n');
    h_res = {'Check', 'Final Order'}; 
    d_align = {'l', 'l'}; 
    h_align = {'c', 'c'};
    table_data = {
        'Logic 3B Applied', mat2str(final_order)
    };
    TestHelper.print_auto_table(h_res, table_data, d_align, h_align);

    % Expectation: D3 (Index 3) moves to top because M3 decides the tie.
    fprintf('\n');
    if TestHelper.check_result(final_order, [3, 1, 2], 'M1_M2_M3 Logic B Swap (D3 moves to top)')
        tests_passed = tests_passed + 1;
    else
        fprintf('\n[Status] FAIL: Logic 3B fallback failed.\n');
    end
    
    % Restore Global RNG
    RandStream.setGlobalStream(s_global);
    
    if tests_passed == 1
        passed = true;
    end
end
