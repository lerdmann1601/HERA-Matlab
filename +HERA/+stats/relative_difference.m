function rel_diff = relative_difference(x, y)
% RELATIVE_DIFFERENCE - Calculates the relative difference between independent means.
%
% Syntax:
%   rel_diff = HERA.stats.relative_difference(x, y)
%
% Description:
%   Computes the relative difference between the means of two samples.
%   Assumes clean data (NaN handling via pairwise exclusion is caller's responsibility).
%
%   Formula: val = abs(mx - my) / abs(mean([mx, my]))
%
% Inputs:
%   x - Column vector of the first sample.
%   y - Column vector of the second sample.
%
% Outputs:
%   rel_diff - Scalar relative difference value.
%
% Author: Lukas von Erdmannsdorff

    % Force column vectors if input is a simple row vector
    if isrow(x) && isvector(x), x = x(:); end
    if isrow(y) && isvector(y), y = y(:); end

    % Calculate means (assumes data is clean, NaN handling is caller's responsibility)
    mx = mean(x, 1);
    my = mean(y, 1);
    
    % Check for zero-sum to avoid division by zero.
    sum_means = mx + my;
    
    % Initialize result
    rel_diff = zeros(size(mx));
    
    % Mask for valid (non-zero) denominators
    valid_mask = sum_means ~= 0;
    
    if any(valid_mask)
         % Calculate relative difference for valid columns.
         % Formula: |mean(x) - mean(y)| / |mean(mean(x) + mean(y))|
         mean_pair = abs((mx(valid_mask) + my(valid_mask)) / 2); 
         diff_val = abs(mx(valid_mask) - my(valid_mask));
         
         % Safety check: Ensure non-zero denominator.
         safe_mask = mean_pair ~= 0;
         
         % Map back to original indices
         valid_indices = find(valid_mask);
         final_indices = valid_indices(safe_mask);
         
         rel_diff(final_indices) = diff_val(safe_mask) ./ mean_pair(safe_mask);
    end
end
