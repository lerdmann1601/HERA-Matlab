function passed = t05_Stability(default_config, thresholds, n_subj, ~, ~)
% T05_STABILITY - Test Case 5: Stability Check (Identical Data)
%
% Description:
%   Tests algorithm stability when absolute zero-variance inputs are provided.
%   Hypothesis: The algorithm must not crash or produce NaNs when comparing 
%   two completely identical, zero-variance datasets (Singularity Check).
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
    
    title_str = 'Test 5: Stability Check (Identical Data / Zero Variance)';
    fprintf('\n%s\n', title_str);
    fprintf('%s\n', repmat('-', 1, strlength(title_str)));
    config = default_config;
    
    fprintf('[Test] Algorithm stability with absolute zero-variance inputs.\n');
    
    data_zero = zeros(n_subj, 3); 
    eff_zero = TestHelper.calculate_real_effects({data_zero}, 1);
    
    % Input Table 
    fprintf('[Input]\n');
    h_in = {'Data Type', 'Size', 'Content'}; 
    d_align = {'l', 'l', 'l'};
    h_align = {'c', 'c', 'c'};
    
    table_data = {
        'Zero Matrix', sprintf('%dx3', n_subj), 'All 0'
    };
    TestHelper.print_auto_table(h_in, table_data, d_align, h_align);
    
    ds_names = {'D1', 'D2', 'D3'}; 
    
    try
        % Run ranking
        [final_order_z, ~] = calculate_ranking({data_zero}, eff_zero, thresholds, config, ds_names, nchoosek(1:3, 2));
        
        % Result Table
        fprintf('\n[Result]\n');
        h_res = {'Outcome', 'Fallback Order'}; 
        d_align = {'l', 'l'};
        h_align = {'c', 'c'};
        
        table_data = {
            'No Crash', mat2str(final_order_z)
        };
        TestHelper.print_auto_table(h_res, table_data, d_align, h_align);
        
        fprintf('\n[Status] PASS: Stability Check successful.\n');
        tests_passed = tests_passed + 1;
    catch ME
        fprintf('\n[Status] FAIL: Crash detected: %s\n', ME.message);
    end
    
    if tests_passed == 1
        passed = true;
    end
end
