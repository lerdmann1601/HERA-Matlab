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
%   Features:
%   - Supports Vectorized Matrix Inputs (N x B).
%   - Assumes clean data (NaN handling via pairwise exclusion is caller's responsibility).
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

    % Check for Matrix (Vectorized) Mode
    % If both inputs are matrices with same number of columns > 1, treat as B independent pairs
    [rx, cx] = size(x);
    [ry, cy] = size(y);
    
    IS_VECTORIZED = (cx > 1) || (cy > 1);
    
    if ~IS_VECTORIZED
        %% STANDARD MODE (Vectors)
        x = x(:);
        y = y(:);
        nx = numel(x);
        ny = numel(y);
        
        if nx == 0 || ny == 0
            d = NaN;
            return;
        end
        B = 1;
    else
        %% VECTORIZED MODE (Matrices)
        % Expecting x and y to support column-wise operations
        % Case 1: x is Nx1, y is NxB (One sample vs Many) -> Expand x
        % Case 2: x is NxB, y is NxB (Many pairwise) -> Column-by-Column
        
        if cx == 1 && cy > 1
            x = repmat(x, 1, cy);
            cx = cy;
        elseif cy == 1 && cx > 1
            y = repmat(y, 1, cx);
            cy = cx;
        end
        
        if cx ~= cy
            error('Dimension mismatch: x and y must have same number of columns for vectorized mode.');
        end
        
        nx = rx; % Rows in x
        ny = ry; % Rows in y
        B = cx;  % Number of batches/bootstraps
        
        if nx == 0 || ny == 0
            d = NaN(1, B);
            return;
        end
    end
    
    % Heuristic threshold for switching algorithms.
    % Empirical testing shows crossover point around N*M = 30,000.
    % Bejond this point, rank-based logic is generally faster.
    IS_SMALL_SAMPLE = (nx * ny) < 30000;
    
    if IS_SMALL_SAMPLE
        %% Method A: Matrix Comparison (O(N*M))
        % efficient for small sample sizes due to low overhead.
        
        if ~IS_VECTORIZED
            % Standard Vector inputs
            % x (col) > y' (row) creates an (nx x ny) logical matrix
            gt = sum(x > y', 'all'); % Number of cases where x_i > y_j
            lt = sum(x < y', 'all'); % Number of cases where x_i < y_j
            
            d = (gt - lt) / (nx * ny);
        else
            % Vectorized Matrix inputs (Column-wise)
            % Assumes data is clean (NaN handling is caller's responsibility).
            % Use 3D expansion: 
            % X: nx x 1 x B
            % Y: 1 x ny x B
            
            X_3D = reshape(x, nx, 1, B);
            Y_3D = reshape(y, 1, ny, B);
            
            % Comparisons (nx x ny x B)
            gt_3D = X_3D > Y_3D;
            lt_3D = X_3D < Y_3D;
            
            % Sum over pages (dim 1 and 2)
            gt = reshape(sum(sum(gt_3D, 1), 2), 1, B);
            lt = reshape(sum(sum(lt_3D, 1), 2), 1, B);
            
            d = (gt - lt) / (nx * ny);
        end
        
    else
        %% Method B: Efficient Rank-Based Calculation (O(N log N))
        % efficient for large sample sizes.
        
        if IS_VECTORIZED
             % Fallback for vectorized large samples (rare in this context)
             % Because 'unique' and 'accumarray' are hard to vectorize across columns without loop.
             % Given B is large, a simple loop over B is acceptable here if N is large.
             d = zeros(1, B);
             for b = 1:B
                 d(b) = HERA.stats.cliffs_delta(x(:,b), y(:,b));
             end
             return;
        end

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