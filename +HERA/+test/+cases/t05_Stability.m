function passed = t05_Stability(default_config, thresholds, n_subj, ~, ~)
% T05_STABILITY - Test Case 5: Stability Check (Identical Data)
%
% Description:
%   Tests algorithm stability when absolute zero-variance inputs are provided.
%   Tests algorithm stability when absolute zero-variance inputs are provided.
%   Hypothesis: The algorithm must not crash and must produce a consistent,
%   deterministic output. NaNs in p-values are expected (singularity) but
%   should not trigger significance flags.
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
        % Run ranking with full output capture (suppress internal logs)
        [~] = evalc('[final_order_z, ~, all_sig_z, ~, p_vals_z] = calculate_ranking({data_zero}, eff_zero, thresholds, config, ds_names, nchoosek(1:3, 2));');
        
        % --- Assertion 1: Effect Sizes must be Zero/Negligible ---
        % For identical zero-variance data, Cliff's Delta and RelDiff must be 0
        d_vals_zero = eff_zero.d_vals_all;
        r_vals_zero = eff_zero.rel_vals_all;
        
        effects_are_zero = all(abs(d_vals_zero(:)) < 1e-12) && all(abs(r_vals_zero(:)) < 1e-12);
        
        % --- Assertion 2: No Significance Flags ---
        % With zero variance, no comparison can be significant
        no_significant_wins = ~any(all_sig_z{1}(:));
        
        % --- Assertion 3: P-Values Safe ---
        % For 0 variance, p-values might be NaN (undefined z-score). 
        % This is acceptable as long as it doesn't trigger significance (checked in Assertion 2).
        p_mat = p_vals_z{1};
        p_values_safe = all(isnan(p_mat(:)) | (p_mat(:) >= 0 & p_mat(:) <= 1));
        
        % --- Assertion 4: Valid Permutation Output ---
        % Result must be a valid permutation of [1,2,3]
        is_valid_permutation = isequal(sort(final_order_z(:)'), [1 2 3]);
        
        % --- Assertion 5: Deterministic Result ---
        % Re-run should produce identical output (suppress internal logs)
        [~] = evalc('[final_order_z2, ~] = calculate_ranking({data_zero}, eff_zero, thresholds, config, ds_names, nchoosek(1:3, 2));');
        is_deterministic = isequal(final_order_z, final_order_z2);
        
        % Result Table
        fprintf('\n[Result]\n');
        h_res = {'Check', 'Value', 'Expected', 'Status'}; 
        d_align = {'l', 'l', 'l', 'c'};
        h_align = {'c', 'c', 'c', 'c'};
        
        table_data = {
            'Effect Sizes', sprintf('d=%.2e, r=%.2e', max(abs(d_vals_zero(:))), max(abs(r_vals_zero(:)))), '~0', char(string(effects_are_zero));
            'Significance Flags', sprintf('%d wins', sum(all_sig_z{1}(:))), '0', char(string(no_significant_wins));
            'P-Values Safe', sprintf('%d NaNs', sum(isnan(p_mat(:)))), 'NaN or [0,1]', char(string(p_values_safe));
            'Valid Permutation', mat2str(final_order_z(:)'), 'Any permutation', char(string(is_valid_permutation));
            'Deterministic', char(string(is_deterministic)), 'true', char(string(is_deterministic))
        };
        TestHelper.print_auto_table(h_res, table_data, d_align, h_align);
        
        % Final Verdict
        all_checks_passed = effects_are_zero && no_significant_wins && p_values_safe && is_valid_permutation && is_deterministic;
        
        if all_checks_passed
            fprintf('\n[Status] PASS: Stability verified: Deterministic output and safe handling of singularities (NaNs allowed but no significance).\n');
            tests_passed = tests_passed + 1;
        else
            fprintf('\n[Status] FAIL: Stability assertions failed.\n');
            if ~effects_are_zero, fprintf('    - Effect sizes not zero.\n'); end
            if ~no_significant_wins, fprintf('    - Unexpected significance flags.\n'); end
            if ~p_values_safe, fprintf('    - P-Values contain invalid numbers (outside 0-1 and not NaN).\n'); end
            if ~is_valid_permutation, fprintf('    - Invalid output permutation.\n'); end
            if ~is_deterministic, fprintf('    - Non-deterministic result.\n'); end
        end
    catch ME
        fprintf('\n[Status] FAIL: Crash detected: %s\n', ME.message);
        fprintf('[Diag] %s at line %d\n', ME.stack(1).file, ME.stack(1).line);
    end
    
    if tests_passed == 1
        passed = true;
    end
end
