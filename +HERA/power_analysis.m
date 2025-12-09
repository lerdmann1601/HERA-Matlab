function power_results = power_analysis(all_data, config, thresholds, num_probanden, num_simulations, pair_idx_all, s, lang)
% POWER_ANALYSIS - Performs a non-parametric post-hoc power analysis.
%
% Syntax:
%   power_results = power_analysis(all_data, config, thresholds, num_probanden, num_simulations, pair_idx_all, s, lang)
%
% Description:
%   This function estimates the statistical power for each pairwise comparison to detect a "win" (significant p-value AND exceeding the effect size thresholds).
%
% Workflow (per metric):
%   1.  Bootstrap Simulation Loop (Parallelized):
%       - Resample subjects with replacement.
%       - Calculate p-values for ALL pairs in this sample.
%       - Apply Holm-Bonferroni correction (Step-Down) to this sample's p-values.
%       - Check Effect Sizes (Cliff's Delta, RelDiff).
%       - Determine "Wins" for this specific sample.
%   2.  Aggregation:
%       - Count how often each pair "won" across all simulations.
%   3.  Power Calculation:
%       - Power = (Win Count / Num Simulations).
%
% Inputs:
%   all_data           - Original data matrices.
%   config             - Configuration struct (contains alphas).
%   thresholds         - Effect size thresholds.
%   num_probanden      - Number of subjects.
%   num_simulations    - Number of bootstrap repetitions for the estimation.
%   pair_idx_all       - Indices of all pairwise comparisons.
%   s                  - RNG stream for reproducibility.
%   lang               - Language pack structure.
%
% Outputs:
%   power_results    - Structure containing the power matrix for each metric.
%
% Author: Lukas von Erdmannsdorff

arguments
    all_data (1,:) cell
    config (1,1) struct
    thresholds (1,1) struct
    num_probanden (1,1) double
    num_simulations (1,1) double
    pair_idx_all (:,2) double
    s
    lang (1,1) struct
end

    num_metrics = numel(all_data);
    num_pairs = size(pair_idx_all, 1);
    power_matrices = cell(1, num_metrics);
    
    fprintf(['  ' lang.power.bootstrap_info '\n'], num_simulations);

    for m = 1:num_metrics
        % Initialize Win Counters for all pairs
        % We count wins for (i,j) and (j,i) separately
        % Map pair index k to (i,j)
        win_counts = zeros(num_pairs, 1); 
        
        data_metric = all_data{m};
        alpha_base = config.alphas(m);
        d_thr = thresholds.d_thresh(m);
        r_thr = thresholds.rel_thresh(m);
        
        % Parallelize the Simulations.
        % This is necessary because Holm-Bonferroni requires all pairs to be calculated together.
        parfor sim = 1:num_simulations
            s_worker = s; % Worker-specific stream for reproducibility
            s_worker.Substream = sim;
            
            % 1. Draw Bootstrap Sample
            boot_indices = randi(s_worker, num_probanden, [num_probanden, 1]);
            
            % 2. Calculate P-Values for all pairs in this sample
            current_p_values = zeros(num_pairs, 1);
            current_d_vals = zeros(num_pairs, 1);
            current_r_vals = zeros(num_pairs, 1);
            valid_comparison = false(num_pairs, 1);
            
            for k = 1:num_pairs
                idx1 = pair_idx_all(k, 1);
                idx2 = pair_idx_all(k, 2);
                
                boot_x_raw = data_metric(boot_indices, idx1);
                boot_y_raw = data_metric(boot_indices, idx2);
                
                % Pairwise deletion
                valid_boot_rows = ~isnan(boot_x_raw) & ~isnan(boot_y_raw);
                boot_x = boot_x_raw(valid_boot_rows);
                boot_y = boot_y_raw(valid_boot_rows);
                n_valid = numel(boot_x);
                
                if n_valid >= 2
                    valid_comparison(k) = true;
                    % P-Value
                    current_p_values(k) = signrank(boot_x, boot_y);
                    
                    % Effect Sizes
                    current_d_vals(k) = HERA.stats.cliffs_delta(boot_x, boot_y);
                    current_r_vals(k) = HERA.stats.relative_difference(boot_x, boot_y);
                else
                    current_p_values(k) = 1.0; % No significance if no data
                    valid_comparison(k) = false;
                end
            end
            
            % 3. Apply Holm-Bonferroni Correction (Step-Down)
            [sorted_p, sort_idx] = sort(current_p_values);
            is_significant = false(num_pairs, 1);
            
            for k = 1:num_pairs
                alpha_k = alpha_base / (num_pairs - k + 1);
                if sorted_p(k) <= alpha_k
                    % Mark the original index as significant
                    is_significant(sort_idx(k)) = true;
                else
                    % Stop Rule: All subsequent p-values are also not significant
                    break;
                end
            end
            
            % 4. Check Win Conditions (Significance + Effect Size)
            % We accumulate wins for the specific pair index k
            current_wins = zeros(num_pairs, 1);
            
            for k = 1:num_pairs
                if valid_comparison(k) && is_significant(k)
                    d_val = current_d_vals(k);
                    r_val = current_r_vals(k);
                    
                    % Check Effect Size Thresholds
                    if abs(d_val) >= (d_thr - 1e-9) && r_val >= (r_thr - 1e-9)
                         current_wins(k) = 1;
                    end
                end
            end
            
            % Accumulate results from this simulation
            win_counts = win_counts + current_wins;
        end
        
        % Calculate Power
        power_matrix = win_counts / num_simulations;
        power_matrices{m} = power_matrix;
    end    
    power_results = struct('power_matrices', {power_matrices});
end
