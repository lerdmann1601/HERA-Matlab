function passed = t13_PathologicalData(default_config, thresholds, n_subj, ~, ~)
% T13_PATHOLOGICALDATA - Test Case 13: Pathological Data (Efron)
%
% Description:
%   Tests the algorithm using Efron's Nontransitive Dice data.
%   Real-world cycle scenario: A>B, B>C, C>D, but D>A.
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
    
    title_str = 'Test 13: Pathological Data (Efron''s Dice)';
    fprintf('\n%s\n', title_str);
    fprintf('%s\n', repmat('-', 1, strlength(title_str)));
    config = default_config;
    
    % Efron's Dice is a famous mathematical set of dice that are non-transitive.
    % A > B, B > C, C > D, but D > A.
    [m1_data, m2_data] = TestHelper.generate_efron_data(n_subj);
    ds_names_4 = {'DiceA', 'DiceB', 'DiceC', 'DiceD'};
    
    fprintf('[Test] Handling of Nontransitive Dice Data (Real-world Cycle).\n');
    fprintf('[Setup] M1: Transitive [40, 30, 20, 10]. M2: Efron''s Dice.\n');
    
    % Input Table 
    fprintf('[Input]\n');
    h_ef = {'Metric', 'Data Structure'}; 
    d_align = {'l', 'l'};
    h_align = {'c', 'c'};
    table_data = {
        'M1', '40 > 30 > 20 > 10';
        'M2', 'A>B (2/3), B>C (2/3), C>D (2/3), D>A (2/3)'
    };
    TestHelper.print_auto_table(h_ef, table_data, d_align, h_align);
    
    config.ranking_mode = 'M1_M2';
    % Suppress warning for expected cycle detection
    warnState = warning('off', 'all'); 
    cleanupObj = onCleanup(@() warning(warnState));
    [~, final_order] = TestHelper.run_single_test({m1_data, m2_data}, thresholds, config, ds_names_4);
    warning('on', 'all');
    
    % We accept any valid permutation of 4 datasets as long as it doesn't crash or hang
    if ~isempty(final_order) && length(unique(final_order)) == 4
        fprintf('\n[Result] Final Order: %s\n', mat2str(final_order));
        fprintf('[Status] PASS: Algorithm terminated with valid rank (Stable Fallback).\n');
        tests_passed = tests_passed + 1;
    else
        fprintf('[Status] FAIL: Cycle Handling failed (Result empty or invalid).\n');
    end
    
    if tests_passed == 1
        passed = true;
    end
end
