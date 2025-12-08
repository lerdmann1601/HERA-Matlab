function [is_significant, adjusted_alphas] = holm_bonferroni(p_values, alpha_global)
% HOLM_BONFERRONI - Performs Holm-Bonferroni step-down correction.
%
% Syntax:
%   [is_significant, adjusted_alphas] = HERA.stats.holm_bonferroni(p_values, alpha_global)
%
% Description:
%   Applies the Holm-Bonferroni method to control the family-wise error rate (FWER).
%   Calculates adjusted thresholds for ALL hypotheses, even non-significant ones.
%
% Inputs:
%   p_values     - Vector of p-values.
%   alpha_global - Scalar global significance level (e.g., 0.05).
%
% Outputs:
%   is_significant  - Logical vector indicating rejected hypotheses.
%   adjusted_alphas - Vector of critical alpha thresholds for each test.
%
% Author: Lukas von Erdmannsdorff

    % Ensure column vector for stable processing
    p_values = p_values(:); 
    m = numel(p_values);
    
    % Sort p-values in ascending order
    [sorted_p, sort_idx] = sort(p_values);
    
    % 1. Calculate thresholds for ALL steps (Vectorized)
    % Formula: alpha / (m - rank + 1)
    k_vec = (1:m)';
    adjusted_alphas_sorted = alpha_global ./ (m - k_vec + 1);
    
    % 2. Step-Down Procedure (Decision Logic)
    is_significant_sorted = false(m, 1);
    
    for k = 1:m
        if sorted_p(k) <= adjusted_alphas_sorted(k)
            is_significant_sorted(k) = true;
        else
            % Stop Rule: If a hypothesis k is not rejected, all subsequent 
            % hypotheses (k+1...m) must also be retained.
            break; 
        end
    end
    
    % 3. Map results back to original order
    is_significant = false(size(p_values));
    adjusted_alphas = zeros(size(p_values));
    
    is_significant(sort_idx) = is_significant_sorted;
    adjusted_alphas(sort_idx) = adjusted_alphas_sorted;
end