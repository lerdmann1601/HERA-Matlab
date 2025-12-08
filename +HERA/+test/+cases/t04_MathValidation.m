function passed = t04_MathValidation(default_config, ~, ~, styles, lang)
% T04_MATHVALIDATION - Test Case 4: Mathematical Validation
%
% Description:
%   Validates mathematical corner cases including:
%   - Zero Variance inputs (Safety check)
%   - BCa Logic (Finite bias/skewness)
%   - SEM Logic (Threshold lower bounds)
%
% Inputs:
%   default_config - (struct) Base configuration
%   styles         - (struct) Plot styles (mock)
%   lang           - (struct) Language struct (mock)
%
% Returns:
%   passed - (logical) True if test passes

    import HERA.test.TestHelper
    import HERA.*
    
    passed = false;
    tests_passed = 0;
    
    title_str = 'Test 4: Mathematical Validation (Thresholds/SEM/BCa)';
    fprintf('\n%s\n', title_str);
    fprintf('%s\n', repmat('-', 1, strlength(title_str)));
    
    fprintf('[Test] Validation of mathematical corner cases (Zero Variance, SEM Logic, BCa Factors).\n');
    
    % Setup: Simple data for math check
    % D1: Zero variance (all 0) -> Should trigger safety mechanisms (NaN protection)
    d1 = zeros(20, 3); 
    % D2: Small variance signal -> Should have normal thresholds
    d2 = repmat([10, 10.1, 10.2], 20, 1) + randn(20, 3) * 0.5; 
    
    config_math = default_config;
    config_math.metric_names = {'ZeroVar', 'SmallVar'};
    
    % Setup Table 
    fprintf('[Setup]\n');
    h_set = {'Metric Name', 'Data Characteristics', 'Expected Behavior'}; 
    d_align = {'l', 'l', 'l'};
    h_align = {'c', 'c', 'c'};
    
    table_data = {
        'ZeroVar', 'All zeros (0 variance)', 'Thresholds=0, BCa Factors=0';
        'SmallVar', 'Mean~10, SD~0.5', 'Valid Thresholds, Finite BCa'
    };
    TestHelper.print_auto_table(h_set, table_data, d_align, h_align);

    % Register automatic cleanup for temp files
    temp_dir = tempname; 
    if ~exist(temp_dir, 'dir'), mkdir(temp_dir); end
    cleanUpTemp = onCleanup(@() rmdir(temp_dir, 's'));
    myStream = RandStream('mlfg6331_64', 'Seed', 1);
    
    try
        % Step 1: Calculate Thresholds (Percentile & SEM)
        % We use evalc to capture/suppress the console output of the function
        cmd_thr = ['[d_thr, r_thr, r_thr_b, min_r_thr, d_v, r_v, p_idx] = ' ...
            'HERA.calculate_thresholds({d1, d2}, 20, config_math, temp_dir, 50, myStream, styles, lang);'];
        [~] = evalc(cmd_thr);
        
        % Step 2: BCa Logic Validation (Bootstrap Confidence Intervals)
        cmd_bca = ['[~, ~, ~, z0_d, a_d, ~, ~, ~, ~, ~, ~, ~, ~] = ' ...
             'HERA.calculate_bca_ci({d1, d2}, d_v, r_v, p_idx, 20, config_math, config_math.metric_names,' ...
             '         temp_dir, temp_dir, 50, myStream, styles, lang, ''Test10'');'];
        [~] = evalc(cmd_bca);
        
        % Checks & Assertions
        % 1. SEM Logic: Final threshold should be max of (Bootstrap, SEM)
        sem_logic_ok = abs(r_thr(2) - max(r_thr_b(2), min_r_thr(2))) < 1e-9;
        
        % 2. BCa Normal: Bias (z0) and Skewness (a) should be finite numbers
        z0_val = z0_d(1, 2); a_val = a_d(1, 2);
        bca_valid_normal = isfinite(z0_val) && isfinite(a_val);
        
        % 3. Zero Variance Safety: z0 and a should be 0 (no bias/skew if no variance)
        z0_zero = z0_d(1, 1); a_zero  = a_d(1, 1);
        protection_ok = (z0_zero == 0 && a_zero == 0);

        % Result Table 
        fprintf('\n[Result]\n');
        h_res = {'Check', 'Value', 'Reference', 'Status'}; 
        d_align = {'l', 'l', 'l', 'l'};
        h_align = {'c', 'c', 'c', 'c'};
        
        table_data = {
            'ZeroVar: Delta Thresh', sprintf('%.4e', d_thr(1)), '< 1e-12', 'OK';
            'ZeroVar: RelDiff Boot', sprintf('%.4e', r_thr_b(1)), '0 (No Var)', '-';
            'ZeroVar: SEM Thr',      sprintf('%.4e', min_r_thr(1)), '0 (No Var)', '-';
            'ZeroVar: Final RelDiff',sprintf('%.4e', r_thr(1)), 'max(Boot, SEM)', '-';
            'SmallVar: Delta Thresh', sprintf('%.4f', d_thr(2)), '> 0', 'OK';
            'SmallVar: RelDiff Boot', sprintf('%.4f', r_thr_b(2)), 'Bootstrap', '-';
            'SmallVar: SEM Thr',      sprintf('%.4f', min_r_thr(2)), 'SEM Lower Bound', '-';
            'SmallVar: Final RelDiff',sprintf('%.4f', r_thr(2)), 'max(Boot, SEM)', char(string(sem_logic_ok));
            'BCa Normal (SmallVar)', sprintf('%.2f / %.2f', z0_val, a_val), 'Finite', char(string(bca_valid_normal));
            'BCa ZeroVar (Safety)', sprintf('%.2f / %.2f', z0_zero, a_zero), '0.00 / 0.00', char(string(protection_ok))
        };
        TestHelper.print_auto_table(h_res, table_data, d_align, h_align);

        % Final Verification
        if d_thr(1) < 1e-12 && d_thr(2) > 0 && sem_logic_ok && bca_valid_normal && protection_ok
             fprintf('\n[Status] PASS: All Mathematical Validations successful.\n');
             tests_passed = tests_passed + 1;
        else
             fprintf('\n[Status] FAIL: Mathematical validation failed.\n');
        end
    catch ME
        fprintf('\n[Status] FAIL: Error: %s\n', ME.message);
    end
    
    if tests_passed == 1
        passed = true;
    end
end
