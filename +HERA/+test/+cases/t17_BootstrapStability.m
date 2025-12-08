function passed = t17_BootstrapStability(default_config, ~, n_subj, styles, lang)
% T17_BOOTSTRAPSTABILITY - Test Case 17: Bootstrap Uncertainty & Stability
%
% Description:
%   Runs a Cluster Bootstrap to analyze the stability of ranks.
%   Ideally, distinct tiers (Top vs Bottom) should have non-overlapping 
%   Confidence Intervals, while close ranks should overlap.
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
    bs_config.B_start = 50; bs_config.B_step = 50; bs_config.B_end = 150;
    bs_config.n_trials = 3; bs_config.convergence_tolerance = 0.05;
    bs_config.min_steps_for_convergence_check = 2; bs_config.smoothing_window = 3; bs_config.convergence_streak_needed = 2;
    
    title_str = 'Test 17: Bootstrap Rank Stability Analysis';
    fprintf('\n%s\n', title_str);
    fprintf('%s\n', repmat('-', 1, strlength(title_str)));
    
    fprintf('[Test] Verification of rank stability and confidence interval separation.\n');
    fprintf('[Setup] Re-generating "Waterfall" data (same profile as Test 15) for Stability Analysis.\n');
    
    % 1. Setup Data
    [m1_data, m2_data, m3_data, ds_names_full] = TestHelper.generate_waterfall_data(n_subj);
    all_data_full = {m1_data, m2_data, m3_data};
    
    config_full = default_config;
    config_full.bootstrap_ranks = bs_config;
    
    % Setup Table 
    fprintf('[Input]\n');
    h_in = {'Parameter', 'Value'}; 
    d_align = {'l', 'l'}; 
    h_align = {'c', 'c'};
    table_data = {
        'Data Profile', 'Waterfall (see Test 15)';
        'Bootstrap Method', 'Cluster Bootstrap';
        'Iterations (B)', '100'
    };
    TestHelper.print_auto_table(h_in, table_data, d_align, h_align);
    
    % Use a temporary directory for intermediate files (automatically deleted later)
    temp_dir = tempname; 
    mkdir(temp_dir);
    cleanUpTemp = onCleanup(@() rmdir(temp_dir, 's'));
    
    try
        % 2. Calculate Thresholds & Base Ranking (Required prerequisites)
        s_boot = RandStream('mlfg6331_64', 'Seed', 999); 
        
        cmd_thr = ['[d_t, r_t, ~, ~, d_v, r_v, p_idx] = HERA.calculate_thresholds(' ...
                   'all_data_full, n_subj, config_full, temp_dir, 50, s_boot, styles, lang);'];
        [~] = evalc(cmd_thr);
            
        thr_full = struct('d_thresh', d_t, 'rel_thresh', r_t);
        eff_full = struct('d_vals_all', d_v, 'rel_vals_all', r_v);
        
        [~, final_rank, ~, ~, ~] = calculate_ranking(...
            all_data_full, eff_full, thr_full, config_full, ds_names_full, p_idx);
        
        % 3. Run Cluster Bootstrap
        % Measures variability of ranks
        cmd_boot = ['[boot_ranks, ~, ~, ~, ~] = HERA.bootstrap_ranking(' ...
                    'all_data_full, thr_full, config_full, ds_names_full, final_rank, p_idx, n_subj, ' ...
                    'temp_dir, temp_dir, 100, s_boot, styles, lang, ''UnitTest_Boot'');'];
        [~] = evalc(cmd_boot);
            
        % 4. Print Percentage Rank Distributions 
        fprintf('\n[Result: Rank Distributions]\n');
        header_parts = {'Rank', 'Dataset', 'Distribution (B=100)'};
        d_align = {'c', 'c', 'l'}; 
        h_align = {'c', 'c', 'c'};
        
        % Inline Helper to format distribution string
        get_dist_str = @(r) strjoin(arrayfun(@(u) sprintf('%d: %.0f%%', u, sum(r==u)/numel(r)*100), ...
                                    unique(r), 'UniformOutput', false), ', ');
        
        [~, sort_idx] = sort(final_rank);
        
        table_data = cell(numel(ds_names_full), 3);
        
        for i = 1:numel(ds_names_full)
            idx = sort_idx(i); 
            name = ds_names_full{idx};
            ranks = boot_ranks(idx, :);
            dist_str = get_dist_str(ranks);
            
            table_data(i, :) = {sprintf('%d', i), name, dist_str};
        end
        TestHelper.print_auto_table(header_parts, table_data, d_align, h_align);
        fprintf('\n');

        % 5. Check Inter-Tier Separation (Representatives)
        idx_d6 = 6; idx_d1 = 1; idx_d11 = 11;
        ci_d6  = quantile(boot_ranks(idx_d6, :), [0.025, 0.975]);
        ci_d1  = quantile(boot_ranks(idx_d1, :), [0.025, 0.975]);
        ci_d11 = quantile(boot_ranks(idx_d11, :), [0.025, 0.975]);
        
        % CIs should NOT overlap between tiers
        sep_top_mid = ci_d6(2) < ci_d1(1); 
        sep_mid_bot = ci_d1(2) < ci_d11(1);
        
        if sep_top_mid && sep_mid_bot
            fprintf('[Status] PASS: Separation confirmed: Robust distinction between Tiers.\n');
        else
            fprintf('[Status] INFO: Overlap detected between Tiers (acceptable for low B/noise).\n');
        end
        % 6. Check Intra-Tier Overlap
        idx_d2 = 2; 
        ranks_d2 = boot_ranks(idx_d2, :);
        ci_d2 = quantile(ranks_d2, [0.025, 0.975]);
        
        % Check for Overlap: Max(LowerBound) < Min(UpperBound)
        overlap_ok = max(ci_d1(1), ci_d2(1)) < min(ci_d1(2), ci_d2(2));
        
        if overlap_ok
            fprintf('[Status] PASS: Overlap confirmed: Correctly identifies D1 and D2 have overlapping CIs.\n');
        else
            fprintf('[Status] FAIL: Unexpected separation of identical datasets.\n');
        end
        
        if sep_top_mid && sep_mid_bot && overlap_ok
            tests_passed = tests_passed + 1;
        end
        
    catch ME
        fprintf('\n[Status] FAIL: Bootstrap execution failed: %s\n', ME.message);
        fprintf('[Diag] %s at line %d\n', ME.stack(1).file, ME.stack(1).line);
    end
    
    if tests_passed == 1
        passed = true;
    end
end
