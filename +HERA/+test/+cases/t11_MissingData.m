function passed = t11_MissingData(default_config, thresholds, ~, ~, ~)
% T11_MISSINGDATA - Test Case 11: Systematic Missing Data Handling
%
% Description:
%   Tests the robustness of the ranking algorithm against varying levels 
%   of missing data (NaN injection).
%
% Inputs:
%   default_config - (struct) Base configuration
%   thresholds     - (struct) Base thresholds
%
% Returns:
%   passed - (logical) True if test passes

    import HERA.test.TestHelper
    import HERA.*
    
    passed = false;
    tests_passed = 0;
    global_seed = 123;
    
    title_str = 'Test 11: Systematic Missing Data (NaN) Handling';
    fprintf('\n%s\n', title_str);
    fprintf('%s\n', repmat('-', 1, strlength(title_str)));
    
    ds_names_5 = {'D1','D2','D3','D4','D5'};
    n_test11 = 200;
    % Base data: clear ladder 10,8,6,4,2
    m1_full = TestHelper.generate_exact_data(n_test11, [10, 8, 6, 4, 2], 2.0);
    config = default_config;
    config.ranking_mode = 'M1';
    pairs_5 = nchoosek(1:5, 2);
    
    fprintf('[Test] Algorithm can handle systematic missing data (NaN).\n');
    
    % Reference Run (Full Data)
    eff_full = TestHelper.calculate_real_effects({m1_full}, 1);
    [order_ref, ~] = calculate_ranking({m1_full}, eff_full, thresholds, config, ds_names_5, pairs_5);
    
    % Setup Table 
    fprintf('[Setup]\n');
    h_set = {'Parameter', 'Description'}; 
    d_align = {'l', 'l'}; 
    h_align = {'c', 'c'};
    table_data = {
        'Data Structure', sprintf('5 Datasets (D1-D5), N=%d per group', n_test11);
        'Missing Data Levels', '20%, 40%, 60%, 80%'
    };
    TestHelper.print_auto_table(h_set, table_data, d_align, h_align);

    % Input Table 
    fprintf('\n[Input]\n');
    h_in = {'Parameter', 'Description'}; 
    d_align = {'l', 'l'}; 
    h_align = {'c', 'c'};
    table_data = {
        'Signal Profile', 'Means: [10, 8, 6, 4, 2], SD: 2.0';
        'Total Data Points', sprintf('%d (Matrix: %dx5)', numel(m1_full), n_test11);
        'Reference Order', mat2str(order_ref)
    };
    TestHelper.print_auto_table(h_in, table_data, d_align, h_align);

    missing_levels = [0.20, 0.40, 0.60, 0.80];
    robustness_passed = true;
    
    % Result Table 
    fprintf('\n[Result]\n');
    h_nan = {'Missing %', 'Order Match', 'Status'}; 
    d_align = {'c', 'l', 'l'}; 
    h_align = {'c', 'c', 'c'};
    
    table_data = cell(length(missing_levels), 3);
    idx_row = 0;

    for lvl = missing_levels
        idx_row = idx_row + 1;
        m1_nan = m1_full;
        
        % Inject NaNs randomly
        num_missing = round(numel(m1_nan) * lvl);
        rng(global_seed + round(lvl*100)); 
        m1_nan(randperm(numel(m1_nan), num_missing)) = NaN;
        
        % Recalculate with missing data
        eff_nan = TestHelper.calculate_real_effects({m1_nan}, 1);
        [order_nan, ~] = calculate_ranking({m1_nan}, eff_nan, thresholds, config, ds_names_5, pairs_5);
        
        match = isequal(order_nan, order_ref);
        status_str = 'PASS';
        if ~match
            % At 80% missing data, degradation is expected/allowed
            if lvl >= 0.80
                status_str = 'INFO: Degradation Exp.';
            else
                status_str = 'FAIL: Unexpected';
                robustness_passed = false;
            end
        end
        
        table_data(idx_row, :) = {sprintf('%d%%', round(lvl*100)), mat2str(order_nan), status_str};
    end
    TestHelper.print_auto_table(h_nan, table_data, d_align, h_align);
    
    if robustness_passed
        fprintf('\n[Status] PASS: Algorithm works robustly with missing data.\n');
        tests_passed = tests_passed + 1;
    else
        fprintf('\n[Status] FAIL: Robustness check failed (unexpected degradation).\n');
    end
    
    if tests_passed == 1
        passed = true;
    end
end
