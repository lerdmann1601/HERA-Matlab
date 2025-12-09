function d = cliffs_delta(x, y)
% CLIFFS_DELTA - Calculates Cliff's Delta effect size using hybrid efficient logic.
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
%   This is a hybrid implementation that automatically selects the fastest method:
%   1. For small samples (N*M < 30,000): Uses vectorized matrix comparison (O(N^2)).
%      This is faster for small N due to lower overhead.
%   2. For large samples: Uses efficient rank-based logic (O(N log N)).
%      This avoids the quadratic complexity scaling.
%
%   Both methods compute (GT - LT) / N. The rank method uses (2*U - N) / N
%   which is mathematically equivalent and numerically stable.
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
    
    % Heuristic threshold for switching algorithms
    % Benchmark on M1 Mac (2025) showed crossover around N=175 (nx*ny ~ 30000)
    % For smaller samples, the vectorization overhead of sort/unique/accumarray 
    % outweighs the O(N^2) cost of simple matrix comparison.
    IS_SMALL_SAMPLE = (nx * ny) < 30000;
    
    if IS_SMALL_SAMPLE
        %% Method A: Naive Matrix Comparison (O(Nx * Ny))
        % Faster for small sample sizes due to low overhead.
        
        % Vectorized calculation of dominance statistics
        % x (col) > y' (row) creates an (nx x ny) logical matrix
        gt = sum(x > y', 'all'); % Number of cases where x_i > y_j
        lt = sum(x < y', 'all'); % Number of cases where x_i < y_j
        
        d = (gt - lt) / (nx * ny);
        
    else
        %% Method B: Efficient Rank-Based Calculation (O(N log N))
        % efficient for large sample sizes.
        
        all_data = [x; y];
        
        % 1. Find unique values and their indices
        [~, ~, idx_map] = unique(all_data);
        
        % 2. Count frequency of each value
        counts = accumarray(idx_map, 1);
        
        % 3. Compute ranks (Mid-Ranks for ties)
        cum_counts = cumsum(counts);
        avg_ranks_unique = cum_counts - (counts - 1) / 2;
        
        % 4. Map ranks back to original data
        ranks = avg_ranks_unique(idx_map);
        
        % 5. U-Statistic calculation
        sum_ranks_x = sum(ranks(1:nx));
        
        % U = R_x - n_x*(n_x + 1)/2
        U = sum_ranks_x - (nx * (nx + 1) / 2);
        
        % 6. Cliff's Delta
        % (2*U - N) / N is algebraically equivalent to (GT - LT) / N
        N_prod = nx * ny;
        d = (2 * U - N_prod) / N_prod;
    end
end