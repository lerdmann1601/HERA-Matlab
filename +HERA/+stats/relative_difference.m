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

    % Calculate means (column-wise if matrices)
    mx = mean(x, 1);
    my = mean(y, 1);
    
    % Robustly handle the case where the denominator might be zero
    sum_means = mx + my;
    
    % Initialize result
    rel_diff = zeros(size(mx));
    
    % Mask for non-zero denominators
    valid_mask = sum_means ~= 0;
    
    if any(valid_mask)
         % Calculate relative difference for valid columns
         % Uses the absolute difference relative to the absolute mean of the pair
         mean_pair = abs((mx(valid_mask) + my(valid_mask)) / 2); %Equivalent to mean([mx, my])
         diff_val = abs(mx(valid_mask) - my(valid_mask));
         
         % Check for 0 in mean_pair (double safety, though coverage by sum_means!=0 usually sufficient for pos numbers)
         safe_mask = mean_pair ~= 0;
         
         % We need to map back to the original full vector indices
         valid_indices = find(valid_mask);
         final_indices = valid_indices(safe_mask);
         
         rel_diff(final_indices) = diff_val(safe_mask) ./ mean_pair(safe_mask);
    end
end
