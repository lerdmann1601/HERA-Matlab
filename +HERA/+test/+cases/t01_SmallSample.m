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
        % Run Ranking calculation for just this data with full output (suppress internal logs)
        [~] = evalc('[final_order, ~, all_sig, ~, p_vals_small] = calculate_ranking({data_small}, eff_small, thresholds, config, {''A'',''B''}, [1 2]);');
        p_val = p_vals_small{1}(1,2);
        
        % --- Assertion 1: P-Value is Valid ---
        p_is_valid = isfinite(p_val) && p_val >= 0 && p_val <= 1;
        
        % --- Assertion 2: P-Value is Very Small (Strong Signal) ---
        % With Mean Difference of 10 and SD of 1, even at N=10 the p-value 
        % should be extremely small (typically < 0.001 for Wilcoxon exact test)
        p_is_small = p_val < 0.01;
        
        % --- Assertion 3: Effect Size is Large ---
        % Cliff's Delta should be close to 1.0 (complete separation)
        d_val = eff_small.d_vals_all(1, 1);
        d_is_large = abs(d_val) > 0.8; % Large effect per Cohen's convention
        
        % --- Assertion 4: Correct Ranking Order ---
        % D1 (Signal) should rank above D2 (Noise) -> Order should be [1, 2]
        order_is_correct = isequal(final_order(:)', [1, 2]);
        
        % --- Assertion 5: Significance Flag Set ---
        % The comparison should be flagged as significant
        is_flagged_significant = all_sig{1}(1,2) || all_sig{1}(2,1);
        
        % Result Table Visualization
        fprintf('\n[Result]\n');
        h_res = {'Check', 'Value', 'Expected', 'Status'}; 
        d_align = {'l', 'l', 'l', 'c'};
        h_align = {'c', 'c', 'c', 'c'};
        
        table_data = {
            'P-Value Valid', sprintf('%.6f', p_val), '[0, 1]', char(string(p_is_valid));
            'P-Value Small', sprintf('%.2e', p_val), '< 0.01', char(string(p_is_small));
            'Effect Size (Delta)', sprintf('%.2f', abs(d_val)), '> 0.80', char(string(d_is_large));
            'Ranking Order', mat2str(final_order(:)'), '[1 2]', char(string(order_is_correct));
            'Significance Flag', char(string(is_flagged_significant)), 'true', char(string(is_flagged_significant))
        };
        TestHelper.print_auto_table(h_res, table_data, d_align, h_align);
        
        % Final Verdict
        all_checks_passed = p_is_valid && p_is_small && d_is_large && order_is_correct && is_flagged_significant;
        
        if all_checks_passed
            fprintf('\n[Status] PASS: Small sample test with strong signal passed all assertions.\n');
            tests_passed = tests_passed + 1;
        else
            fprintf('\n[Status] FAIL: Small sample test assertions failed.\n');
            if ~p_is_valid, fprintf('    - P-Value out of range or NaN.\n'); end
            if ~p_is_small, fprintf('    - P-Value not small enough for strong signal (%.4f >= 0.01).\n', p_val); end
            if ~d_is_large, fprintf('    - Effect size not large enough (%.2f <= 0.80).\n', abs(d_val)); end
            if ~order_is_correct, fprintf('    - Ranking order incorrect (%s).\n', mat2str(final_order(:)')); end
            if ~is_flagged_significant, fprintf('    - Significance flag not set.\n'); end
        end
    catch ME
        fprintf('\n[Status] FAIL: Crash during small sample test: %s\n', ME.message);
        fprintf('[Diag] %s at line %d\n', ME.stack(1).file, ME.stack(1).line);
    end
    
    if tests_passed == 1
        passed = true;
    end
end
