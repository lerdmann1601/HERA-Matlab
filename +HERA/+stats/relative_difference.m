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

    % Calculate means
    mx = mean(x);
    my = mean(y);
    
    % Robustly handle the case where the denominator might be zero
    if (mx + my) == 0
        rel_diff = 0;
    else
        % Calculate relative difference
        % Uses the absolute difference relative to the absolute mean of the pair
        rel_diff = abs(mx - my) / abs(mean([mx, my]));
    end
end
