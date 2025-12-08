function passed = t16_FullIntegration(default_config, ~, n_subj, styles, lang)
% T16_FULLINTEGRATION - Test Case 16: Full System Integration (Dynamic)
%
% Description:
%   Tests the full pipeline (Auto-Thresholds + Ranking) on the complex 
%   "Waterfall" dataset (15 Datasets, 4 Tiers).
%   Ensures that dynamic thresholding works correctly in conjunction with 
%   the tiered ranking logic.
%
% Inputs:
%   default_config - (struct) Base configuration
%   n_subj         - (int)    Sample size
%   styles         - (struct) Plot styles
%   lang           - (struct) Language struct
%
% Returns:
%   passed - (logical) True if test passes

    import HERA.test.TestHelper
    import HERA.*
    
    passed = false;
    tests_passed = 0;
    
    title_str = 'Test 16: Full System Integration (15 Datasets + Auto-Thresholds)';
    fprintf('\n%s\n', title_str);
    fprintf('%s\n', repmat('-', 1, strlength(title_str)));
    
    fprintf('[Test] Validating the full logic chain (M2 correction, M3A tie-break, M3B sort).\n');
    fprintf('[Setup] Re-generating "Waterfall" data (same profile as Test 15) to test Dynamic Thresholds.\n');
    
    % Scenario for full logic chain validation (same as T15 but with calculated thresholds)
    [m1_data, m2_data, m3_data, ds_names_full] = TestHelper.generate_waterfall_data(n_subj);
    all_data_full = {m1_data, m2_data, m3_data};
    config_full = default_config;
    temp_dir = tempname; mkdir(temp_dir);
    cleanUpTemp = onCleanup(@() rmdir(temp_dir, 's'));
    myStream = RandStream('mlfg6331_64', 'Seed', 123);
    
    % Input Table 
    fprintf('[Input: Tiers & Logic Targets]\n');
    h_in = {'Tier', 'Datasets', 'Target Logic'}; 
    d_align = {'l', 'l', 'l'};
    h_align = {'c', 'c', 'c'};
    table_data = {
        'Tier 1', 'D1-D5', 'Base Rank (Mid) -> Displaced by Tier 2';
        'Tier 2', 'D6-D10', 'Correction (Top) -> Wins via M2';
        'Tier 3', 'D11-D12', 'Tie-Break -> Logic 3A (M3 Swap)';
        'Tier 4', 'D13-D15', 'Sorting -> Logic 3B (M3 Gradient)'
    };
    TestHelper.print_auto_table(h_in, table_data, d_align, h_align);

    try
        % Step 1: Auto-Calculate Thresholds
        cmd_full = ['[d_t_full, r_t_full, ~, ~, d_v_full, r_v_full, p_i_full] = ' ...
            'HERA.calculate_thresholds(all_data_full, n_subj, config_full, temp_dir, 50, myStream, styles, lang);'];
        [~] = evalc(cmd_full);
        
        % Step 2: Run Full Ranking 
        thr_full = struct('d_thresh', d_t_full, 'rel_thresh', r_t_full);
        eff_full = struct('d_vals_all', d_v_full, 'rel_vals_all', r_v_full);
        
        [final_order_full, ~, ~, ~, ~] = calculate_ranking(all_data_full, eff_full, thr_full, config_full, ds_names_full, p_i_full);
        
        % Step 3: Display Results
        % Re-calculate means for display consistency
        mean_m1 = mean(m1_data); mean_m2 = mean(m2_data); mean_m3 = mean(m3_data);
        
        % Result Table 
        fprintf('\n[Result]\n');
        header_parts = {'Rank', 'Dataset', 'M1 Mean', 'M2 Mean', 'M3 Mean'};
        d_align = {'c', 'c', 'c', 'c', 'c'}; 
        h_align = {'c', 'c', 'c', 'c', 'c'};
        
        table_data = cell(15, 5);
        for r = 1:15
            idx = final_order_full(r);
            table_data(r, :) = {
                sprintf('%d', r), ds_names_full{idx}, ...
                sprintf('%.1f', mean_m1(idx)), ...
                sprintf('%.1f', mean_m2(idx)), ...
                sprintf('%.1f', mean_m3(idx))
            };
        end
        TestHelper.print_auto_table(header_parts, table_data, d_align, h_align);
        fprintf('\n');

        % Validation checks:
        top_5 = final_order_full(1:5);
        tier2_success = all(ismember(top_5, 6:10)); % Check if Tier 2 (M2 strong) is on top
        
        mid_5 = final_order_full(6:10);
        tier1_success = all(ismember(mid_5, 1:5)); % Check if Tier 1 (M1 strong) is 2nd
        
        rank_d11 = find(final_order_full == 11);
        rank_d12 = find(final_order_full == 12);
        logic3a_success = rank_d12 < rank_d11; % Check if M3A flipped D11/D12
        
        rank_d13 = find(final_order_full == 13);
        rank_d14 = find(final_order_full == 14);
        rank_d15 = find(final_order_full == 15);
        logic3b_success = (rank_d15 < rank_d14) && (rank_d14 < rank_d13); % Check M3B sort
        
        if tier2_success && tier1_success && logic3a_success && logic3b_success
             fprintf('[Status] PASS: Full System Validation Successful.\n');
             tests_passed = tests_passed + 1;
        else
             fprintf('[Status] FAIL: Logic Validation Failed.\n');
             if ~tier2_success, fprintf('    - M2 Correction failed (Tier 2 not at top).\n'); end
             if ~tier1_success, fprintf('    - Tier 1 not in mid.\n'); end
             if ~logic3a_success, fprintf('    - Logic 3A failed (D12 did not beat D11).\n'); end
             if ~logic3b_success, fprintf('    - Logic 3B failed (D13-D15 sorting incorrect).\n'); end
        end
        
    catch ME
        fprintf('[Status] FAIL: Pipeline crashed: %s\n', ME.message);
        fprintf('[Diag] %s\n', ME.stack(1).name);
    end
    
    if tests_passed == 1
        passed = true;
    end
end
