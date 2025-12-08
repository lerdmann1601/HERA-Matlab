function passed = t07_InitialSorting(default_config, thresholds, n_subj, ~, ~)
% T07_INITIALSORTING - Test Case 7: Initial Sorting Logic
%
% Description:
%   Tests the hierarchy of sorting criteria in Mode 1 (M1):
%   1. Significant Win Count
%   2. Stochastic Dominance (Sum of Cliff's Deltas)
%   3. Arithmetic Mean (Tie-Breaker)
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
    
    title_str = 'Test 7: Initial Sorting Logic (Dominance, Delta, Mean)';
    fprintf('\n%s\n', title_str);
    fprintf('%s\n', repmat('-', 1, strlength(title_str)));
    
    ds_names = {'D1', 'D2', 'D3'};
    
    % The initial sorting uses 3 criteria in order:
    % 1. Win Count (Number of significant wins against others)
    % 2. Stochastic Dominance (Sum of Cliff's Deltas in direct comparison)
    % 3. Arithmetic Mean (Tie-breaker if d=0)
    
    % Step 1: Standard Dominance (Clear Winner)
    fprintf('[Test] Step 1: Testing Standard Dominance (M1).\n');
    means_1 = [10, 5, 0]; sd_1 = 2.0;
    fprintf('[Setup] Datasets: D1=%.1f, D2=%.1f, D3=%.1f (SD=%.1f)\n', means_1(1), means_1(2), means_1(3), sd_1);
    m1_data = TestHelper.generate_exact_data(n_subj, means_1, sd_1); 
    config = default_config;
    config.ranking_mode = 'M1';
    
    [final_order_1, ~, eff_1, p_val_1] = TestHelper.run_single_test_full({m1_data}, thresholds, config, ds_names);
    
    % Combined Diagnostic Table 
    fprintf('\n[Result]\n');
    h_diag = {'Pair', 'P-Value', 'Delta', 'RelDiff'};
    d_align = {'l', 'l', 'l', 'l'}; 
    h_align = {'c', 'c', 'c', 'c'};
    
    pairs_idx = nchoosek(1:3, 2); 
    table_data = cell(size(pairs_idx, 1), 4);
    for k = 1:size(pairs_idx, 1)
        i = pairs_idx(k, 1); j = pairs_idx(k, 2);
        pair_label = sprintf('%s vs %s', ds_names{i}, ds_names{j});
        
        table_data(k, :) = {
            pair_label, ...
            sprintf('%.4e', p_val_1{1}(i, j)), ...
            sprintf('%.4f', eff_1.d_vals_all(k, 1)), ...
            sprintf('%.4f', eff_1.rel_vals_all(k, 1))
        };
    end
    TestHelper.print_auto_table(h_diag, table_data, d_align, h_align);
    fprintf('\n');
    
    pass_1 = TestHelper.check_result(final_order_1, [1, 2, 3], 'M1 Clear Winner (D1>D2>D3)');
    
    % Step  2: Tie-Break via Delta Injection 
    fprintf('\n[Test] Step 2: Testing Tie-Break via Delta Injection (M1).\n');
    fprintf('[Setup] Generating Noise Data (Means=0). Forcing Ties on Win-Counts.\n');
    % Data: Identical noise (Neutral p-values -> Win Counts all equal)
    m1_data_noise = TestHelper.generate_exact_data(n_subj, [0, 0, 0], 1.0);
    effect_sizes = TestHelper.calculate_real_effects({m1_data_noise}, 1);
    
    % Injection: Manually override Cliff's Delta to simulate subtle dominance
    % D3 > D2 (-0.4), D3 > D1 (-0.8), D2 > D1 (-0.4).
    % Note: Negative Delta typically means Group 2 > Group 1. 
    effect_sizes.d_vals_all(:, 1) = [-0.4; -0.8; -0.4]; 
    
    % Input Table for Injection 
    fprintf('[Input]\n');
    h_inj = {'Pair', 'Injected Delta', 'Implication'}; 
    d_align = {'l', 'c', 'l'}; 
    h_align = {'c', 'c', 'c'};
    table_data = {
        'D1 vs D2', '-0.4', 'D2 > D1';
        'D1 vs D3', '-0.8', 'D3 > D1';
        'D2 vs D3', '-0.4', 'D3 > D2'
    };
    TestHelper.print_auto_table(h_inj, table_data, d_align, h_align);
    
    import HERA.calculate_ranking
    [final_order_2, ~, ~, ~, ~] = calculate_ranking({m1_data_noise}, effect_sizes, thresholds, config, ds_names, nchoosek(1:3, 2));
    
    % Result Table
    fprintf('\n[Result]\n');
    h_res = {'Metric', 'Final Order'}; 
    d_align = {'l', 'l'}; 
    h_align = {'c', 'c'};
    
    table_data = {
        'Tie-Break (Delta)', mat2str(final_order_2)
    };
    TestHelper.print_auto_table(h_res, table_data, d_align, h_align);

    % Expect: 3, 2, 1 (Due to injected Deltas)
    fprintf('\n');
    pass_2 = TestHelper.check_result(final_order_2, [3, 2, 1], 'M1 Tie-Break via Delta');
    
    % Step 3: Tie-Break via Mean 
    fprintf('\n[Test] Step 3: Testing Tie-Break via Mean (M1).\n');
    % Means are very close, but distinct.
    means_3 = [0.01, 0.02, 0.00];
    
    % Input Table
    fprintf('[Input]\n');
    h_in = {'Dataset', 'Mean', 'Condition'}; 
    d_align = {'l', 'l', 'l'}; 
    h_align = {'c', 'c', 'c'};
    
    table_data = {
        'D1', sprintf('%.2f', means_3(1)), 'Delta Forced to 0';
        'D2', sprintf('%.2f', means_3(2)), 'Delta Forced to 0';
        'D3', sprintf('%.2f', means_3(3)), 'Delta Forced to 0'
    };
    TestHelper.print_auto_table(h_in, table_data, d_align, h_align);

    m1_data_mean = TestHelper.generate_exact_data(n_subj, means_3, 2.0);
    eff_mean = TestHelper.calculate_real_effects({m1_data_mean}, 1);
    eff_mean.d_vals_all(:) = 0; % Force Delta to 0 to simulate perfect stochastic equality
    
    [final_order_3, ~, ~, ~, ~] = calculate_ranking({m1_data_mean}, eff_mean, thresholds, config, ds_names, nchoosek(1:3, 2));
    
    % Result Table
    fprintf('\n[Result]\n');
    h_res = {'Metric', 'Final Order'}; 
    d_align = {'l', 'l'}; 
    h_align = {'c', 'c'};
    
    table_data = {
        'Tie-Break (Mean)', mat2str(final_order_3)
    };
    TestHelper.print_auto_table(h_res, table_data, d_align, h_align);
    
    % Expect: D2 (0.02) > D1 (0.01) > D3 (0.00) -> [2, 1, 3]
    fprintf('\n');
    pass_3 = TestHelper.check_result(final_order_3, [2, 1, 3], 'M1 Tie-Break via Mean (D2>D1>D3)');
    
    if pass_1 && pass_2 && pass_3
        tests_passed = tests_passed + 1;
    end
    
    if tests_passed == 1
        passed = true;
    end
end
