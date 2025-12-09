function d = cliffs_delta(x, y)
% CLIFFS_DELTA - Calculates Cliff's Delta effect size using efficient rank logic.
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
%   Implementation Note:
%   This implementation uses an efficient O(N log N) sorting logic
%   to avoid the quadratic complexity of naive pairwise comparisons.
%   The calculation uses (2*U - N) / N, which is mathematically and numerically
%   exactly the definition (GT - LT) / N. This avoids floating point
%   discrepancies that would occur with (2*U)/N - 1.
%
% Inputs:
%   x - Column vector of the first sample.
%   y - Column vector of the second sample.
%
% Outputs:
%   d - Scalar value of Cliff's Delta (-1 to +1).
%
% Author: Lukas von Erdmannsdorff

    % Ensure column vectors
    x = x(:);
    y = y(:);
    
    nx = numel(x);
    ny = numel(y);
    
    % Handle empty inputs
    if nx == 0 || ny == 0
        d = NaN;
        return;
    end
    
    %% Efficient Rank Calculation 
    all_data = [x; y];
    
    % 1. Find unique values and their indices
    % 'idx_map' maps each value to its index in the sorted list of unique values
    [~, ~, idx_map] = unique(all_data);
    
    % 2. Count frequency of each value
    counts = accumarray(idx_map, 1);
    
    % 3. Compute ranks (Mid-Ranks for ties)
    % The average rank of a group of identical values is calculated from the
    % cumulative count minus half the group size.
    cum_counts = cumsum(counts);
    avg_ranks_unique = cum_counts - (counts - 1) / 2;
    
    % 4. Map ranks back to original data
    ranks = avg_ranks_unique(idx_map);
    
    %% U-Statistic and Cliff's Delta
    
    % 1. Sum of ranks for the first group (x)
    sum_ranks_x = sum(ranks(1:nx));
    
    % 2. Calculate U-statistic: U = R_x - n_x*(n_x + 1)/2
    % This is always an integer (or x.5, which multiplied by 2 becomes an integer)
    U = sum_ranks_x - (nx * (nx + 1) / 2);
    
    % 3. Product of sample sizes
    N = nx * ny;
    
    % 4. Bit-perfect calculation:
    % We move the subtraction into the numerator to avoid floating-point rounding
    % errors that occur with the Division-then-Subtract order.
    % (2*U - N) is mathematically identical to (GT - LT).
    d = (2 * U - N) / N;
end