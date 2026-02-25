function passed = t15_SystemTest(default_config, thresholds, n_subj, ~, ~)
% T15_SYSTEMTEST - Test Case 15: Global System Test (Waterfall)
%
% Description:
%   Integration test using a "Waterfall" data profile to test all logic tiers:
%   - M2 Correction (Top Tier)
%   - M1 Base Ranking (Mid Tier)
%   - M3 Sorting (Bottom Tier)
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
    
    title_str = 'Test 15: Global System Test (Static Thresholds)';
    fprintf('\n%s\n', title_str);
    fprintf('%s\n', repmat('-', 1, strlength(title_str)));
    
    % Integration test for the entire hierarchy logic using "Waterfall" data.
    fprintf('[Test] Multi-Tier "Waterfall" sorting (M2 Correction, M1 Base, M3 Tail).\n');
    
    % Generate 15 datasets split into Tiers (see helper function)
    [m1_data, m2_data, m3_data, ds_names_15] = TestHelper.generate_waterfall_data(n_subj);
    
    % Calculate actual means and SDs for transparency in the log
    mu1 = mean(m1_data); mu2 = mean(m2_data); mu3 = mean(m3_data);
    sd1 = std(m1_data);  sd2 = std(m2_data);  sd3 = std(m3_data);
    
    % Input Table (Detailed Data Profile)
    % We show representatives for each Tier to verify the generated signal values.
    fprintf('[Input: Data Profile (Waterfall Structure)]\n');
    h_in = {'Tier', 'Rep.', 'M1 (Mean)', 'M1 (SD)', 'M2 (Mean)', 'M2 (SD)', 'M3 (Mean)', 'M3 (SD)', 'Intended Logic'}; 
    d_align = {'l', 'c', 'c', 'c', 'c', 'c', 'c', 'c', 'l'}; 
    h_align = {'c', 'c', 'c', 'c', 'c', 'c', 'c', 'c', 'c'};
    
    table_data = {
        'Top (Tier 2)', 'D6',  sprintf('%.1f', mu1(6)),  sprintf('%.2f', sd1(6)),  sprintf('%.1f', mu2(6)),  sprintf('%.2f', sd2(6)),  sprintf('%.1f', mu3(6)),  sprintf('%.2f', sd3(6)),  'M2 High -> Correction';
        'Mid (Tier 1)', 'D1',  sprintf('%.1f', mu1(1)),  sprintf('%.2f', sd1(1)),  sprintf('%.1f', mu2(1)),  sprintf('%.2f', sd2(1)),  sprintf('%.1f', mu3(1)),  sprintf('%.2f', sd3(1)),  'M1 High -> Base Rank';
        'Low (Tier 3)', 'D11', sprintf('%.1f', mu1(11)), sprintf('%.2f', sd1(11)), sprintf('%.1f', mu2(11)), sprintf('%.2f', sd2(11)), sprintf('%.1f', mu3(11)), sprintf('%.2f', sd3(11)), 'M3A Swap Logic';
        'Bot (Tier 4)', 'D13', sprintf('%.1f', mu1(13)), sprintf('%.2f', sd1(13)), sprintf('%.1f', mu2(13)), sprintf('%.2f', sd2(13)), sprintf('%.1f', mu3(13)), sprintf('%.2f', sd3(13)), 'M3B Sorting Logic'
    };
    TestHelper.print_auto_table(h_in, table_data, d_align, h_align);
    
    % Run with Static Thresholds (Standard Config)
    config = default_config;
    eff = TestHelper.calculate_real_effects({m1_data, m2_data, m3_data}, 3);
    pairs_15 = nchoosek(1:15, 2);
    
    [final_order, ~, ~, ~, ~] = calculate_ranking({m1_data, m2_data, m3_data}, eff, thresholds, config, ds_names_15, pairs_15);
    
    % Validation by Tier logic
    % Calculate average rank of each tier to verify separation
    tier2_ranks = []; % Should be top (1-5)
    for i = 6:10, tier2_ranks = [tier2_ranks, find(final_order == i)]; end
    
    tier1_ranks = []; % Should be mid (6-10)
    for i = 1:5, tier1_ranks = [tier1_ranks, find(final_order == i)]; end
    
    tier3_ranks = []; % Should be bottom (11-15)
    for i = 11:15, tier3_ranks = [tier3_ranks, find(final_order == i)]; end
    
    mean_rank_t2 = mean(tier2_ranks);
    mean_rank_t1 = mean(tier1_ranks);
    mean_rank_t3 = mean(tier3_ranks);
    
    % Result Table 
    fprintf('\n[Result]\n');
    h_res = {'Group', 'Mean Rank', 'Expectation', 'Status'}; 
    d_align = {'l', 'c', 'l', 'c'}; 
    h_align = {'c', 'c', 'c', 'c'};
    
    % Check logic
    t2_ok = mean_rank_t2 < mean_rank_t1; % Tier 2 (M2) beat Tier 1 (M1)
    t1_ok = mean_rank_t1 < mean_rank_t3; % Tier 1 (M1) beat Tier 3 (M3)
    
    table_data = {
        'Tier 2 (M2)', sprintf('%.1f', mean_rank_t2), 'Top (~3.0)', char(string(t2_ok));
        'Tier 1 (M1)', sprintf('%.1f', mean_rank_t1), 'Mid (~8.0)', char(string(t1_ok));
        'Tier 3 (M3)', sprintf('%.1f', mean_rank_t3), 'Bot (~13.0)', '-'
    };
    TestHelper.print_auto_table(h_res, table_data, d_align, h_align);
    
    if t2_ok && t1_ok
        fprintf('\n[Status] PASS: Full Hierarchy respected (M2 > M1 > M3).\n');
        tests_passed = tests_passed + 1;
    else
         fprintf('\n[Status] FAIL: Hierarchy violation.\n');
         if ~t2_ok, fprintf('    - M2 Correction failed (Tier 2 not above Tier 1).\n'); end
         if ~t1_ok, fprintf('    - M1 Base Logic failed (Tier 1 not above Tier 3).\n'); end
    end
    
    if tests_passed == 1
        passed = true;
    end
end
