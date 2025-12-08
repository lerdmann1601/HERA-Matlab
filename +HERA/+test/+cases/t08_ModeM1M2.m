function passed = t08_ModeM1M2(default_config, thresholds, n_subj, ~, ~)
% T08_MODEM1M2 - Test Case 8: Mode M1_M2 Logic
%
% Description:
%   Tests the Rank Correction Swap logic.
%   In M1_M2 mode, Metric 1 establishes the base order. Metric 2 acts as a 
%   correction layer: if a lower-ranked dataset significantly beats a 
%   higher-ranked one in M2, they swap places.
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
    
    title_str = 'Test 8: Mode M1_M2 (Rank Correction Swap)';
    fprintf('\n%s\n', title_str);
    fprintf('%s\n', repmat('-', 1, strlength(title_str)));
    
    ds_names = {'D1', 'D2', 'D3'};
    config = default_config;
    
    % Hypothesis: In M1_M2 mode, Metric 1 establishes the base order. 
    % Metric 2 acts as a correction layer: if a lower-ranked dataset significantly 
    % beats a higher-ranked one in M2, they swap places.
    
    fprintf('[Test] Verifying M2 Correction Swap logic.\n');
    fprintf('[Setup] M1 Order: D1 > D2 > D3. M2 Signal: D2 > D1 (Contradiction).\n');
    
    % M1: D1 wins (10 > 5 > 0)
    m1_data = TestHelper.generate_exact_data(n_subj, [10, 5, 0], 1.0);
    % M2: D2 wins against D1 (10 > 0) -> Should trigger swap
    m2_data = TestHelper.generate_exact_data(n_subj, [0, 10, 0], 1.0);
    
    config.ranking_mode = 'M1_M2';

    % Input Table 
    fprintf('[Input]\n');
    h_in = {'Metric', 'Signal Structure', 'Implication'}; 
    d_align = {'l', 'l', 'l'};
    h_align = {'c', 'c', 'c'};
    table_data = {
        'M1', 'D1(10) > D2(5)', 'D1 wins (Initial)';
        'M2', 'D2(10) > D1(0)', 'D2 wins (Strong Correction)'
    };
    TestHelper.print_auto_table(h_in, table_data, d_align, h_align);
    
    [final_order, ~, eff, p_val] = TestHelper.run_single_test_full({m1_data, m2_data}, thresholds, config, ds_names);
    
    % Display M2 values specifically to prove correction validity
    m2_d = eff.d_vals_all(:, 2);
    m2_r = eff.rel_vals_all(:, 2);
    m2_p = p_val{2}; 
    p_val_d1_d2 = m2_p(1, 2);

    % Result Table for M2 Stats 
    fprintf('\n[Result: M2 Stats D1 vs D2]\n');
    h_m2 = {'Metric', 'Value', 'Threshold', 'Significant?'}; 
    d_align = {'l', 'l', 'l', 'c'}; 
    h_align = {'c', 'c', 'c', 'c'};
    
    table_data = {
        'Delta', sprintf('%.2f', abs(m2_d(1))), sprintf('%.2f', thresholds.d_thresh(2)), char(string(abs(m2_d(1))>thresholds.d_thresh(2)));
        'RelDiff', sprintf('%.2f', m2_r(1)), sprintf('%.2f', thresholds.rel_thresh(2)), char(string(m2_r(1)>=thresholds.rel_thresh(2)));
        'P-Value', sprintf('%.4e', p_val_d1_d2), '0.05', char(string(p_val_d1_d2 < 0.05))
    };
    TestHelper.print_auto_table(h_m2, table_data, d_align, h_align);
    
    % Expectation: D2 swaps with D1 -> Order: D2, D1, D3 (2, 1, 3)
    fprintf('\n');
    if TestHelper.check_result(final_order, [2, 1, 3], 'M1_M2 Swap Applied (D2 overtakes D1)')
        tests_passed = tests_passed + 1;
    end
    
    if tests_passed == 1
        passed = true;
    end
end
