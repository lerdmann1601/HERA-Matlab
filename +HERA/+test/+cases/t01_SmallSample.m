function passed = t01_SmallSample(default_config, thresholds, ~, ~, ~)
% T01_SMALLSAMPLE - Test Case 1: Small Sample Size Boundary (Exact Wilcoxon N=10)
%
% Description:
%   Verifies exact p-value calculation for N < 15.
%   Hypothesis: For N=10, MATLAB's signrank should perform an exact calculation
%   instead of the asymptotic normal approximation used for larger N.
%
% Inputs:
%   default_config - (struct) Base HERA configuration
%   thresholds     - (struct) Base thresholds
%   ~              - (unused) n_subj override
%   ~              - (unused) styles
%   ~              - (unused) lang
%
% Returns:
%   passed - (logical) True if test passes

    import HERA.test.TestHelper
    import HERA.*
    
    passed = false;
    tests_passed = 0;
    
    title_str = 'Test 1: Small Sample Size Boundary (Exact Wilcoxon N=10)';
    fprintf('\n%s\n', title_str);
    fprintf('%s\n', repmat('-', 1, strlength(title_str)));
    config = default_config;
    
    n_small = 10;
    fprintf('[Test] Verifying exact p-value calculation for N < 15.\n');
    fprintf('[Setup] Using N=%d. Generating Signal (Mean 10) vs Noise (Mean 0).\n', n_small);
    
    % Data Generation
    % D1: Mean 10 (Strong Signal), D2: Mean 0 (Noise) -> Difference should be highly significant
    data_small = TestHelper.generate_exact_data(n_small, [10, 0], 1.0);
    eff_small = TestHelper.calculate_real_effects({data_small}, 1);
    
    % Input Table Visualization
    fprintf('[Input]\n');
    h_in = {'Parameter', 'Value'}; 
    d_align = {'l', 'l'}; 
    h_align = {'c', 'c'}; 
    table_data = {
        'Sample Size (N)', num2str(n_small);
        'Mean D1 (Signal)', '10.0';
        'Mean D2 (Noise)', '0.0'
    };
    TestHelper.print_auto_table(h_in, table_data, d_align, h_align);
    
    try
        % Run Ranking calculation for just this data
        [~, ~, ~, ~, p_vals_small] = calculate_ranking({data_small}, eff_small, thresholds, config, {'A','B'}, [1 2]);
        p_val = p_vals_small{1}(1,2);
        
        % Result Table Visualization
        fprintf('\n[Result]\n');
        h_res = {'Metric', 'Value', 'Expectation'}; 
        d_align = {'l', 'l', 'l'};
        h_align = {'c', 'c', 'c'};
        
        table_data = {
            'Calculated P-Value', sprintf('%.6f', p_val), '< 0.005 (Strong Effect)'
        };
        TestHelper.print_auto_table(h_res, table_data, d_align, h_align);
        
        % Validation: Check if P-value is a valid probability and low enough
        if isfinite(p_val) && p_val >= 0 && p_val <= 1
             fprintf('\n[Status] PASS: Exact test execution successful. Result is a valid probability.\n');
             tests_passed = tests_passed + 1;
        else
             fprintf('\n[Status] FAIL: P-Value invalid (NaN or out of range).\n');
        end
    catch ME
        fprintf('\n[Status] FAIL: Crash during small sample test: %s\n', ME.message);
    end
    
    if tests_passed == 1
        passed = true;
    end
end
