function power_results = power_analysis(all_data, all_alpha_matrices, thresholds, num_probanden, num_simulations, pair_idx_all, s, lang)
% POWER_ANALYSIS - Performs a non-parametric post-hoc power analysis.
%
% Syntax:
%   power_results = power_analysis(all_data, all_alpha_matrices, thresholds, num_probanden, num_simulations, pair_idx_all, s, lang)
%
% Description:
%   This function estimates the statistical power for each pairwise comparison to detect a "win" (significant p-value AND exceeding the effect size thresholds).
%
% Workflow (per comparison):
%   1.  Bootstrap Simulation: `num_simulations` experiments are simulated by drawing
%       subject pairs from the original dataset with replacement.
%   2.  Win Check: For each simulation, the complete criteria check
%       (p-value, Cliff's Delta, RelDiff) is performed.
%   3.  Power Calculation: Power is the percentage of simulations in which
%       a "win" was successfully detected.
%
% Inputs:
%   all_data           - Original data matrices.
%   all_alpha_matrices - Holm-corrected alpha values for each comparison.
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
% Author:   Lukas von Erdmannsdorff
% Date:     12.10.2025
% Version:  1

    num_metrics = numel(all_data);
    num_pairs = size(pair_idx_all, 1);
    power_matrices = cell(1, num_metrics);
    
    fprintf(['  ' lang.power.bootstrap_info '\n'], num_simulations);

    for m = 1:num_metrics
        power_matrix = NaN(num_pairs, 1);
        data_metric = all_data{m};
        
        % parfor can be used here to process pairs in parallel
        parfor p = 1:num_pairs
            s_worker = s; % Worker-specific stream for reproducibility
            s_worker.Substream = p;
            
            idx1 = pair_idx_all(p, 1);
            idx2 = pair_idx_all(p, 2);
            
            win_counter = 0;
            for j = 1:num_simulations
                % Draw a cluster bootstrap sample (may contain NaNs)
                boot_indices = randi(s_worker, num_probanden, [num_probanden, 1]);
                boot_x_raw = data_metric(boot_indices, idx1);
                boot_y_raw = data_metric(boot_indices, idx2);
           
                % Perform pairwise deletion within the simulation
                valid_boot_rows = ~isnan(boot_x_raw) & ~isnan(boot_y_raw);
                boot_x = boot_x_raw(valid_boot_rows);
                boot_y = boot_y_raw(valid_boot_rows);
                n_valid = numel(boot_x);
                % Only perform tests if valid pairs exist
                if n_valid >= 2 % signrank requires at least 2 pairs
                    % Calculate p-value for this sample
                    p_val = signrank(boot_x, boot_y);
        
                    % Calculate effect sizes for this sample
                    gt = sum(boot_x > boot_y', 'all');
                    lt = sum(boot_x < boot_y', 'all');
                    d_val = (gt - lt) / (n_valid^2); % Use dynamic n_valid
                    
                    mx = mean(boot_x);
                    my = mean(boot_y);
                    r_val = abs(mx - my) / abs(mean([mx, my]));
                    if isnan(r_val), r_val = 0; end
        
                    % Check all win criteria
                    alpha_holm = all_alpha_matrices{m}(idx1, idx2);
                    if p_val <= alpha_holm && abs(d_val) >= thresholds.d_thresh(m) && r_val >= thresholds.rel_thresh(m)
                        win_counter = win_counter + 1;
                    end
                end
            end
            
            % Power is the proportion of "won" simulations
            power_matrix(p) = win_counter / num_simulations;
        end
                power_matrices{m} = power_matrix;
    end    
    power_results = struct('power_matrices', {power_matrices});
end
