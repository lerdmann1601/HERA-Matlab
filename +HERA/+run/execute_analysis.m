function analysis_results = execute_analysis(all_data, num_probanden, dataset_names, pair_idx_all, userInput, setupData, styles)
% EXECUTE_ANALYSIS - Perform statistical analysis and ranking.
%
% Syntax:
%   analysis_results = execute_analysis(all_data, num_probanden, dataset_names, pair_idx_all, userInput, setupData, styles)
%
% Description:
%   This function handles Steps 6-10 of the ranking analysis:
%   6. Statistical Analysis (Thresholds, CIs)
%   7. Ranking Calculation (+ Sensitivity)
%   8. Borda Ranking
%   9. Bootstrap Analysis of Ranks
%   10. Power Analysis
%
% Inputs:
%   all_data       - Cell array of data matrices.
%   num_probanden  - Number of subjects.
%   dataset_names  - Cell array of dataset names.
%   pair_idx_all   - Matrix of pair indices.
%   userInput      - User input structure (for ranking_mode, sensitivity flags).
%   setupData      - Environment setup structure (config, paths, lang, rng).
%   styles         - Plotting styles structure.
%
% Outputs:
%   analysis_results - Structure containing all calculation results.
%
% Author: Lukas von Erdmannsdorff

    import HERA.*
    
    % Unpack necessary variables from setupData
    config = setupData.config;
    lang = setupData.lang;
    graphics_dir = setupData.graphics_dir;
    csv_dir = setupData.csv_dir;
    s = setupData.s;
    base_name = setupData.base_name;
    
    metric_names = config.metric_names;

    %% 6. Statistical Analysis and CI Calculation
    % Determine the statistical thresholds for effect sizes (Cliff's Delta, RelDiff).
    fprintf(['\n' lang.run_ranking.calc_thresholds '\n']);
    [d_thresh, rel_thresh, rel_thresh_b, min_rel_thresh, d_vals_all, rel_vals_all, pair_idx_all, selected_B, stability_data_thr, ...
        h_fig_thr_global, h_fig_thr_detailed, h_fig_hist_thr, h_fig_hist_raw] = ...
         calculate_thresholds(all_data, num_probanden, config, graphics_dir, config.manual_B_thr, s, styles, lang);

    % Calculate the Bias-Corrected and accelerated (BCa) confidence intervals.
    fprintf(['\n' lang.run_ranking.calc_bca '\n']);
    [selected_B_ci, ci_d_all, ci_r_all, z0_d_all, a_d_all, z0_r_all, a_r_all, stability_data_ci, h_fig_bca_global, h_fig_bca_detailed, ...
        h_fig_hist_z0, h_fig_hist_a, h_fig_hist_widths] = ...
        calculate_bca_ci(all_data, d_vals_all, rel_vals_all, pair_idx_all, ...
            num_probanden, config, metric_names, graphics_dir, csv_dir, config.manual_B_ci, s, styles, lang, base_name);

    % Bundle results into structs for clean passing to other functions.
    effect_sizes = struct('d_vals_all', d_vals_all, 'rel_vals_all', rel_vals_all);
    thresholds   = struct('d_thresh', d_thresh, 'rel_thresh', rel_thresh, 'rel_thresh_b', rel_thresh_b, 'min_rel_thresh', min_rel_thresh);

    %% 7. Ranking Calculation (+ Sensitivity Analysis)
    % The ranking is calculated for the primary hierarchy and all other selected permutations.
    num_permutations = size(userInput.selected_permutations, 1);
    num_datasets = length(dataset_names); % Derived from names
    all_permutation_ranks = zeros(num_datasets, num_permutations);
    
    fprintf(['\n' lang.run_ranking.calc_ranking_permutations '\n'], num_permutations);

    % Pre-declare variables to avoid warnings (though not strictly necessary in MATLAB if assigned in loop p=1)
    final_order = []; final_rank = []; all_sig_matrices = []; all_alpha_matrices = []; all_p_value_matrices = []; 
    swap_details = []; intermediate_orders = [];

    for p = 1:num_permutations
        % Define the primary hierarchy from the first permutation entry.
        primary_hierarchy_indices = userInput.selected_permutations(1, :);
        % Get the current metric order for this iteration.
        current_permutation = userInput.selected_permutations(p, :);

        % Find the mapping from the current permutation back to the primary order.
        [~, local_indices] = ismember(current_permutation, primary_hierarchy_indices);
        
        % Reorder all necessary data and configurations for the current ranking run.
        perm_metric_names = metric_names(local_indices);
        perm_all_data = all_data(local_indices);
        perm_d_vals = d_vals_all(:, local_indices);
        perm_rel_vals = rel_vals_all(:, local_indices);
        perm_d_thresh = d_thresh(local_indices);
        perm_rel_thresh = rel_thresh(local_indices);
        perm_config = config;
        perm_config.alphas = config.alphas(local_indices);
        perm_config.ranking_mode = userInput.ranking_mode; 
        
        fprintf([' -> ' lang.run_ranking.permutation_progress '\n'], p, num_permutations, strjoin(perm_metric_names, ' -> '));
        
        % Prepare the reordered data structs for the ranking function.
        perm_effect_sizes = struct('d_vals_all', perm_d_vals, 'rel_vals_all', perm_rel_vals);
        perm_thresholds = struct('d_thresh', perm_d_thresh, 'rel_thresh', perm_rel_thresh);
        
        % Calculate the ranking for the current permutation and store the result.
        [~, single_run_rank] = calculate_ranking(perm_all_data, perm_effect_sizes, perm_thresholds, perm_config, dataset_names, pair_idx_all);
        all_permutation_ranks(:, p) = single_run_rank;

        % The first permutation (p=1) is the primary hierarchy; store its detailed results.
        if p == 1
             % Call with the original (primary) data and config (which now includes ranking_mode)
             [final_order, final_rank, all_sig_matrices, all_alpha_matrices, all_p_value_matrices, swap_details, intermediate_orders] = ...
                 calculate_ranking(all_data, effect_sizes, thresholds, config, dataset_names, pair_idx_all);
        end
    end

    %% 8. Borda Ranking Calculation
    % If sensitivity analysis is enabled, calculate the consensus rank using the Borda method.
    if userInput.run_sensitivity_analysis
        fprintf(['\n' lang.run_ranking.calc_borda '\n']);
        borda_results = borda_ranking(all_permutation_ranks, dataset_names);
    else
        borda_results = []; % Set to empty if not performed.
    end

    %% 9. Bootstrap Analysis of Ranks (for Primary Hierarchy)
    % Assess the stability of the primary ranking using a cluster bootstrap.
    fprintf(['\n' lang.run_ranking.bootstrap_ranks '\n']);
    [final_bootstrap_ranks, selected_B_rank, stability_data_rank, h_figs_rank, h_fig_hist_rank, ci_lower_rank, ci_upper_rank] = ...
        bootstrap_ranking(all_data, thresholds, config, dataset_names, final_rank, pair_idx_all, num_probanden, graphics_dir, csv_dir,...
        config.manual_B_rank, s, styles, lang, base_name);

    %% 10. Power Analysis
    % If enabled, perform a post-hoc, non-parametric power analysis.
    if userInput.run_power_analysis
        fprintf(['\n' lang.run_ranking.calc_power '\n']);
        power_results = power_analysis(all_data, config, thresholds, num_probanden, userInput.power_simulations, pair_idx_all, s, lang);
    else
        power_results = []; % Set to empty if not performed.
    end

    %% Pack Results
    analysis_results.d_vals_all = d_vals_all;
    analysis_results.rel_vals_all = rel_vals_all;
    analysis_results.thresholds = thresholds;
    analysis_results.pair_idx_all = pair_idx_all;
    analysis_results.selected_B = selected_B;
    analysis_results.stability_data_thr = stability_data_thr;
    
    analysis_results.selected_B_ci = selected_B_ci;
    analysis_results.ci_d_all = ci_d_all;
    analysis_results.ci_r_all = ci_r_all;
    analysis_results.z0_d_all = z0_d_all;
    analysis_results.a_d_all = a_d_all;
    analysis_results.z0_r_all = z0_r_all;
    analysis_results.a_r_all = a_r_all;
    analysis_results.stability_data_ci = stability_data_ci;
    
    analysis_results.final_order = final_order;
    analysis_results.final_rank = final_rank;
    analysis_results.all_sig_matrices = all_sig_matrices;
    analysis_results.all_alpha_matrices = all_alpha_matrices;
    analysis_results.all_p_value_matrices = all_p_value_matrices;
    analysis_results.swap_details = swap_details;
    analysis_results.intermediate_orders = intermediate_orders;
    
    analysis_results.all_permutation_ranks = all_permutation_ranks;
    analysis_results.borda_results = borda_results;
    
    analysis_results.final_bootstrap_ranks = final_bootstrap_ranks;
    analysis_results.selected_B_rank = selected_B_rank;
    analysis_results.stability_data_rank = stability_data_rank;
    analysis_results.ci_lower_rank = ci_lower_rank;
    analysis_results.ci_upper_rank = ci_upper_rank;
    
    analysis_results.power_results = power_results;
    
    % Pack Handles
    analysis_results.handles.h_fig_thr_global = h_fig_thr_global;
    analysis_results.handles.h_fig_thr_detailed = h_fig_thr_detailed;
    analysis_results.handles.h_fig_hist_thr = h_fig_hist_thr;
    analysis_results.handles.h_fig_hist_raw = h_fig_hist_raw;
    analysis_results.handles.h_fig_bca_global = h_fig_bca_global;
    analysis_results.handles.h_fig_bca_detailed = h_fig_bca_detailed;
    analysis_results.handles.h_fig_hist_z0 = h_fig_hist_z0;
    analysis_results.handles.h_fig_hist_a = h_fig_hist_a;
    analysis_results.handles.h_fig_hist_widths = h_fig_hist_widths;
    analysis_results.handles.h_figs_rank = h_figs_rank;
    analysis_results.handles.h_fig_hist_rank = h_fig_hist_rank;

end
