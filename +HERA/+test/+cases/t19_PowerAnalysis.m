function passed = t19_PowerAnalysis(default_config, thresholds, ~, ~, lang)
% T19_POWERANALYSIS - Test Case 19: Post-Hoc Power Analysis
%
% Description:
%   Sanity check to verify that Power Analysis returns plausible values 
%   (0 <= P <= 1) and yields >0.8 for a purposely Strong Signal case.
%
% Inputs:
%   default_config - (struct) Base configuration
%   thresholds     - (struct) Base thresholds
%   lang           - (struct) Language struct
%
% Returns:
%   passed - (logical) True if test passes

    import HERA.test.TestHelper
    import HERA.*
    
    passed = false;
    tests_passed = 0;
    global_seed = 123;
    s_global = RandStream('mlfg6331_64', 'Seed', global_seed);
    
    title_str = 'Test 19: Post-Hoc Power Analysis (Sanity Check)';
    fprintf('\n%s\n', title_str);
    fprintf('%s\n', repmat('-', 1, strlength(title_str)));
    
    fprintf('[Test] Verifying that Power Analysis returns plausible values (0 <= Power <= 1).\n');
    
    % Setup: Strong Signal (Effect Size large) -> High Power expected
    % D1 Mean=10, D2 Mean=0. SD=2. N=20. Cliff's Delta will be approx 1.0.
    n_pwr = 20;
    data_pwr = TestHelper.generate_exact_data(n_pwr, [10, 0], 2.0);
    all_data_pwr = {data_pwr};
    ds_names_pwr = {'HighSignal', 'LowSignal'};
    pairs_pwr = [1 2];
    
    % Setup Table
    fprintf('[Setup]\n');
    h_set = {'Dataset', 'Mean', 'SD', 'Description'};
    d_align = {'l', 'c', 'c', 'l'};
    h_align = {'c', 'c', 'c', 'c'};
    
    table_data = {
        'HighSignal', '10.0', '2.0', 'Strong Effect Source';
        'LowSignal',  '0.0',  '2.0', 'Noise/Reference Group'
    };
    TestHelper.print_auto_table(h_set, table_data, d_align, h_align);

    % Config & Prerequisites
    config_pwr = default_config;
    config_pwr.ranking_mode = 'M1';
    eff_pwr = TestHelper.calculate_real_effects(all_data_pwr, 1);
    
    % Calculate Alpha Matrices first (needed as input for power_analysis)
    [~, ~, ~, all_alphas_pwr, ~] = calculate_ranking(all_data_pwr, eff_pwr, thresholds, config_pwr, ds_names_pwr, pairs_pwr);
    
    % Simulation Settings
    n_sims = 100; % Low number for unit test speed
    
    % Input Table (Simulation Parameters & Win Criteria)
    fprintf('\n[Input]\n');
    h_in = {'Parameter', 'Value', 'Condition'}; 
    d_align = {'l', 'l', 'l'}; 
    h_align = {'c', 'c', 'c'};
    
    % Extract values for display
    d_thr_disp = thresholds.d_thresh(1);
    r_thr_disp = thresholds.rel_thresh(1);
    alpha_disp = all_alphas_pwr{1}(1,2);

    table_data = {
        'Subjects (n)', num2str(n_pwr), 'Resampling Size';
        'Simulations', num2str(n_sims), 'Bootstrap Iterations';
        'Alpha Limit', sprintf('%.3f', alpha_disp), 'Significant if p < Alpha';
        'Effect Limit', sprintf('d > %.2f, r > %.2f', d_thr_disp, r_thr_disp), 'Relevant if Eff > Limit'
    };
    TestHelper.print_auto_table(h_in, table_data, d_align, h_align);
    
    try
        % Call the function
        cmd_pwr = ['power_results = HERA.power_analysis(' ...
                   'all_data_pwr, config_pwr, thresholds, n_pwr, n_sims, pairs_pwr, s_global, lang);'];
        [~] = evalc(cmd_pwr);
        
        % Extract Result
        pwr_val = power_results.power_matrices{1}(1);
        
        % Result Table
        fprintf('\n[Result]\n');
        h_res = {'Metric', 'Value', 'Plausibility'}; 
        d_align = {'l', 'l', 'l'}; 
        h_align = {'c', 'c', 'c'};
        
        % Expectation: With such a strong signal (Mean 10 vs 0), power should be perfect or near perfect.
        is_plausible = (pwr_val >= 0.8) && (pwr_val <= 1.0);
        
        table_data = {
            'Calculated Power', sprintf('%.2f', pwr_val), '> 0.80 (Expected)'
        };
        TestHelper.print_auto_table(h_res, table_data, d_align, h_align);
        
        if is_plausible
            fprintf('\n[Status] PASS: Power Analysis produced a logical high value for strong signal.\n');
            tests_passed = tests_passed + 1;
        else
            fprintf('\n[Status] FAIL: Power Value implausible (Expected > 0.8 for strong signal).\n');
        end
        
    catch ME
        fprintf('\n[Status] FAIL: Power Analysis crashed: %s\n', ME.message);
        fprintf('\n[Diag] %s\n', ME.stack(1).name);
    end
    
    if tests_passed == 1
        passed = true;
    end
end
