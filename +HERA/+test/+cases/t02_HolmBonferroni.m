function passed = t02_HolmBonferroni(default_config, thresholds, n_subj, ~, ~)
% T02_HOLMBONFERRONI - Test Case 2: Statistical Rigor (Holm-Bonferroni Correction)
%
% Description:
%   Verifies dynamic alpha adaptation for multiple comparisons.
%   Hypothesis: The alpha threshold must adapt dynamically to the number of 
%   pairwise comparisons (m).
%   Formula: Alpha_corrected = Alpha_base / (m - rank + 1)
%
% Inputs:
%   default_config - (struct) Base HERA configuration
%   thresholds     - (struct) Base thresholds
%   n_subj         - (int)    Sample size
%
% Returns:
%   passed - (logical) True if test passes

    import HERA.test.TestHelper
    import HERA.*
    
    passed = false;
    tests_passed = 0;
    
    title_str = 'Test 2: Statistical Rigor (Holm-Bonferroni Correction)';
    fprintf('\n%s\n', title_str);
    fprintf('%s\n', repmat('-', 1, strlength(title_str)));
    
    % PART A: Basic Functionality (N=5 Datasets)
    fprintf('[Test] Checking dynamic alpha adaptation for N=5 datasets (10 pairs).\n');
    n_ds = 5;
    ds_names_hb = {'D1', 'D2', 'D3', 'D4', 'D5'};
    
    fprintf('[Setup] Generating random noise data. Alpha Base = 0.05.\n');
    % We use noise data (means all 0) so P-values are random/uniform under Null Hypothesis.
    data_hb = TestHelper.generate_exact_data(n_subj, zeros(1, n_ds), 1.0);
    
    eff_hb = TestHelper.calculate_real_effects({data_hb}, 1);
    config_hb = default_config;
    config_hb.ranking_mode = 'M1';
    config_hb.alphas = [0.05];
    
    % Create all unique pairs (n choose 2)
    pairs_hb = nchoosek(1:n_ds, 2);
    n_pairs_hb = size(pairs_hb, 1); % 10 pairs
    
    % Input Table 
    fprintf('[Input]\n');
    h_in = {'Parameter', 'Value'}; 
    d_align = {'l', 'l'};
    h_align = {'c', 'c'};
    table_data = {
        'Datasets', num2str(n_ds);
        'Comparisons (Pairs)', num2str(n_pairs_hb);
        'Base Alpha', sprintf('%.2f', config_hb.alphas(1))
    };
    TestHelper.print_auto_table(h_in, table_data, d_align, h_align);

    % Run ranking on noise data
    % We are interested in 'all_alphas' output
    [~, ~, ~, all_alphas, all_p_vals] = calculate_ranking({data_hb}, eff_hb, thresholds, config_hb, ds_names_hb, pairs_hb);
    
    % Validation 1: The strictest alpha should be Base / Num_Pairs (Classical Check)
    min_alpha_observed = min(all_alphas{1}(:));
    expected_min_alpha = 0.05 / n_pairs_hb; % 0.005
    
    % Validation 2: Verify Exact Correspondence with Centralized Function
    % This ensures that calculate_ranking is consistently using the HERA.stats.holm_bonferroni logic.
    % We must extract strictly the unique pairs (upper triangle) to match the M hypotheses hypothesis.
    p_vals_vec = zeros(n_pairs_hb, 1);
    observed_alphas_vec = zeros(n_pairs_hb, 1);
    for k = 1:n_pairs_hb
        p_vals_vec(k) = all_p_vals{1}(pairs_hb(k,1), pairs_hb(k,2));
        observed_alphas_vec(k) = all_alphas{1}(pairs_hb(k,1), pairs_hb(k,2));
    end
    
    [~, expected_alphas_vec] = HERA.stats.holm_bonferroni(p_vals_vec, config_hb.alphas(1));
    
    % Calculate sorting match
    % The helper returns alphas mapped to the input vector.
    alphas_match = max(abs(observed_alphas_vec - expected_alphas_vec)) < 1e-10;

    % Result Table 
    fprintf('\n[Result]\n');
    h_res = {'Metric', 'Value'}; 
    d_align = {'l', 'l'};
    h_align = {'c', 'c'};
    
    table_data = {
        'Expected Min Alpha', sprintf('%.6f', expected_min_alpha);
        'Observed Min Alpha', sprintf('%.6f', min_alpha_observed);
        'Difference', sprintf('%.1e', abs(min_alpha_observed - expected_min_alpha));
        'Full Vector Match', string(alphas_match)
    };
    TestHelper.print_auto_table(h_res, table_data, d_align, h_align);
    
    if abs(min_alpha_observed - expected_min_alpha) < 1e-6 && alphas_match
        fprintf('\n[Status] PASS: Correction correctly applied (Verified against HERA.stats.holm_bonferroni).\n');
        tests_passed = tests_passed + 1;
    else
        fprintf('\n[Status] FAIL: Alpha correction mismatch or helper function divergence.\n');
    end
    
    if tests_passed == 1
        passed = true;
    end
end
