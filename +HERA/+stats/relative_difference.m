function rel_diff = relative_difference(x, y)
% RELATIVE_DIFFERENCE - Calculates the relative difference between independent means.
%
% Syntax:
%   rel_diff = HERA.stats.relative_difference(x, y)
%
% Description:
%   Computes the relative difference between the means of two samples.
%   Includes a safety check for division by zero if the sum of means is zero.
%
%   Features:
%   - Supports Vectorized NaN Handling (using omitnan and pairwise exclusion).
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

    % Check for NaNs and enforce pairwise exclusion.
    % If a subject is missing in x, it is also ignored in y (and vice versa).
    if any(isnan(x), 'all') || any(isnan(y), 'all')
        nan_mask = isnan(x) | isnan(y);
        x(nan_mask) = NaN;
        y(nan_mask) = NaN;
        
        % Use omitnan to ignore the masked values
        mx = mean(x, 1, 'omitnan');
        my = mean(y, 1, 'omitnan');
    else
        % Calculate means (Standard path)
        mx = mean(x, 1);
        my = mean(y, 1);
    end
    
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
