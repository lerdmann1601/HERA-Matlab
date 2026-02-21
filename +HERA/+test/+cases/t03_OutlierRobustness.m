function passed = t03_OutlierRobustness(default_config, thresholds, n_subj, ~, ~)
% T03_OUTLIERROBUSTNESS - Test Case 3: Outlier Robustness & Threshold Logic
%
% Description:
%   Tests the "AND" condition in ranking logic:
%   Win = (p < alpha) AND (|d| > d_thresh) AND (r > r_thresh).
%   Verifies that missing just one criterion results in Neutral rank.
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
    
    global_seed = 123; % From main script context
    
    title_str = 'Test 3: Outlier Robustness & Threshold Logic (Divergence Check)';
    fprintf('\n%s\n', title_str);
    fprintf('%s\n', repmat('-', 1, strlength(title_str)));
    config = default_config;
    
    fprintf('[Test] Win condition = (p < alpha) AND (|d| > d_thr) AND (r > r_thr).\n');
    fprintf('[Setup] N=4 Candidates, n=%d Subjects, Fixed Thresholds (B=N/A).\n', n_subj);
    fprintf('        Generating synthetic edge cases to test partial criteria matches.\n');
    
    % Case 1: High Delta, but trivial RelDiff (Microscopic difference with 0 variance)
    d1_micro = ones(n_subj, 1) * 1.001;
    d2_micro = ones(n_subj, 1) * 1.000;
    
    % Case 2: High RelDiff, but low Delta (Caused by a single massive outlier)
    rng(global_seed);
    d3_noise = randn(n_subj, 1);
    d4_outlier = d3_noise; 
    d4_outlier(1) = 1000; % Massive Mean shift due to one value, Median unchanged
    
    data_robust = [d1_micro, d2_micro, d3_noise, d4_outlier];
    names_robust = {'Micro1', 'Micro2', 'Noise', 'Outlier'};
    
    % Set explicit thresholds for verification
    thr_robust = thresholds;
    thr_robust.d_thresh = [0.2, 0.2, 0.2]; % Delta must be > 0.2
    thr_robust.rel_thresh = [0.05, 0.05, 0.05]; % RelDiff must be > 5%
    
    % Input Table 
    fprintf('[Input]\n');
    h_in = {'Comparison', 'Scenario', 'Target Failure'}; 
    d_align = {'l', 'l', 'l'};
    h_align = {'c', 'c', 'c'};
    
    table_data = {
        'Micro1 vs Micro2', 'High Delta (determin.), Trivial RelDiff', 'Fail RelDiff';
        'Noise vs Outlier', 'High Mean Diff, Low Median Diff (d)', 'Fail Delta';
        'Thresholds', 'Delta > 0.20, RelDiff > 0.05', '-'
    };
    TestHelper.print_auto_table(h_in, table_data, d_align, h_align);
    
    eff_robust = TestHelper.calculate_real_effects({data_robust}, 1);
    
    % Analyze Pair 1: Micro1 vs Micro2
    d_p1 = eff_robust.d_vals_all(1, 1); 
    r_p1 = eff_robust.rel_vals_all(1, 1);
    
    % Analyze Pair 2: Noise vs Outlier (Index 3 vs 4)
    d_p2 = eff_robust.d_vals_all(6, 1); 
    r_p2 = eff_robust.rel_vals_all(6, 1);
    
    % Run Ranking
    [~, ~, all_sig_rob, ~, ~] = calculate_ranking({data_robust}, eff_robust, thr_robust, config, names_robust, nchoosek(1:4, 2));
    
    % Check Results: Did any pair register a "Win"?
    win_micro = all_sig_rob{1}(1, 2) || all_sig_rob{1}(2, 1);
    win_outlier = all_sig_rob{1}(3, 4) || all_sig_rob{1}(4, 3);
    
    % Result Table
    fprintf('\n[Result]\n');
    h_res = {'Pair', 'Delta', 'RelDiff', 'Win (Exp: 0)'}; 
    d_align = {'l', 'c', 'c', 'c'}; 
    h_align = {'c', 'c', 'c', 'c'};
        
    table_data = {
        'Micro1 vs Micro2', sprintf('%.2f', abs(d_p1)), sprintf('%.4f', r_p1), sprintf('%d', win_micro);
        'Noise vs Outlier', sprintf('%.2f', abs(d_p2)), sprintf('%.2f', r_p2), sprintf('%d', win_outlier)
    };
    TestHelper.print_auto_table(h_res, table_data, d_align, h_align);
    
    % Validation
    if ~win_micro && ~win_outlier
        fprintf('\n[Status] PASS: Partial criteria matches correctly rejected.\n');
        tests_passed = tests_passed + 1;
    else
        fprintf('\n[Status] FAIL: Robustness logic failed. MicroWin: %d, OutlierWin: %d\n', win_micro, win_outlier);
    end
    
    if tests_passed == 1
        passed = true;
    end
end
