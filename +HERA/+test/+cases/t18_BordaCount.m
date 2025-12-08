function passed = t18_BordaCount(default_config, thresholds, n_subj, ~, ~)
% T18_BORDACOUNT - Test Case 18: Borda Count Logic
%
% Description:
%   Tests simple Borda Count aggregation, which is used to combine 
%   rankings from Sensitivity Analysis simulations.
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
    
    title_str = 'Test 18: Borda Count Aggregation Logic Check';
    fprintf('\n%s\n', title_str);
    fprintf('%s\n', repmat('-', 1, strlength(title_str)));
    
    % Borda Count aggregates rankings from Sensitivity Analysis.
    
    % Simulation: 3 Datasets, 2 Permutations
    sim_ranks = [1, 2; 2, 1; 3, 3]; % Rows=Datasets, Cols=Permutations
    sim_names = {'D1', 'D2', 'D3'};
    
    fprintf('[Test] Correctness of Borda score calculation.\n');
    fprintf('[Setup] Simulating specific rank distribution for 3 datasets over 2 permutations.\n');
    
    % Input Table 
    fprintf('[Input]\n');
    h_in = {'Dataset', 'Ranks (P1, P2)', 'Points Calc'}; 
    d_align = {'c', 'c', 'l'}; 
    h_align = {'c', 'c', 'c'};
    table_data = {
        'D1', '1, 2', '2 + 1 = 3';
        'D2', '2, 1', '1 + 2 = 3';
        'D3', '3, 3', '0 + 0 = 0'
    };
    TestHelper.print_auto_table(h_in, table_data, d_align, h_align);
    
    borda_res = HERA.borda_ranking(sim_ranks, sim_names);
    
    score_d1 = borda_res.score(1);
    score_d3 = borda_res.score(3);
    
    % Result Table 
    fprintf('\n[Result]\n');
    h_b = {'Dataset', 'Score', 'Expected'}; 
    d_align = {'c', 'c', 'c'}; 
    h_align = {'c', 'c', 'c'};
    
    table_data = {
        'D1', sprintf('%.1f%%', score_d1), '75.0%';
        'D3', sprintf('%.1f%%', score_d3), '0.0%'
    };
    TestHelper.print_auto_table(h_b, table_data, d_align, h_align);
    
    if abs(score_d1 - 75.0) < 1e-6 && abs(score_d3 - 0.0) < 1e-6
        fprintf('\n[Status] PASS: Borda calculation mathematically correct.\n');
        tests_passed = tests_passed + 1;
    else
        fprintf('\n[Status] FAIL: Borda calculation mismatch.\n');
    end
    
    if tests_passed == 1
        passed = true;
    end
end
