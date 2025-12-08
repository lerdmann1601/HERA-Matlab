function d = cliffs_delta(x, y)
% CLIFFS_DELTA - Calculates Cliff's Delta effect size.
%
% Syntax:
%   d = HERA.stats.cliffs_delta(x, y)
%
% Description:
%   Computes the non-parametric effect size Cliff's Delta.
%   It measures the degree of overlap ("stochastic dominance") between two distributions.
%
%   Note on Independence:
%   Standard Cliff's Delta assumes independent samples. In HERA, data is typically
%   paired (same subjects). We use this standard formulation intentionally as a
%   conservative estimate of dominance to assess practical relevance.
%   Statistical significance for the paired structure is separately handled
%   by the Wilcoxon signed-rank test.
%
%   Formula: d = (#(x > y) - #(x < y)) / (nx * ny)
%
% Inputs:
%   x - Column vector of the first sample.
%   y - Column vector of the second sample.
%
% Outputs:
%   d - Scalar value of Cliff's Delta (-1 to +1).
%
% Author: Lukas von Erdmannsdorff

    % Ensure column vectors to guarantee correct matrix expansion
    % This prevents errors if row vectors are passed accidentally
    x = x(:);
    y = y(:);

    % Calculate sample sizes
    nx = size(x, 1);
    ny = size(y, 1);
    
    % Vectorized calculation of dominance statistics
    % x (col) > y' (row) creates an (nx x ny) logical matrix of all-pairs comparisons
    % sum(..., 'all') counts the total number of true comparisons
    gt = sum(x > y', 'all'); % Number of cases where x_i > y_j
    lt = sum(x < y', 'all'); % Number of cases where x_i < y_j
    
    % Final calculation
    % Note: If nx or ny is 0, this returns NaN
    d = (gt - lt) / (nx * ny);
end