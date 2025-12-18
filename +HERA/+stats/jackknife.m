function [vals, a] = jackknife(x, y, metric_type, vec_limit)
% JACKKNIFE - Calculates Jackknife statistics for Cliff's Delta or Relative Difference.
%
% Syntax:
%   [vals, a] = HERA.stats.jackknife(x, y, metric_type, [vec_limit])
%
% Description:
%   Computes the Jackknife (leave-one-out) statistics for a given metric.
%   It iteratively removes one subject (pair) and calculates the effect size
%   on the remaining N-1 subjects.
%
%   Additionally, it calculates the acceleration factor (a) for BCa confidence
%   intervals, utilizing the 3rd moment of the jackknife distribution.
%
%   This function strictly handles missing data (NaNs) by performing pairwise exclusion
%   within each Jackknife subsample. It is optimized for performance by avoiding
%   redundant memory allocations where possible.
%
%   Implementation Note:
%   This is a hybrid implementation that automatically selects the fastest method.
%   Empirical benchmarks (M1 MBP, Dec 2025) determined N ~ 150 as the crossover point:
%   1. For standard samples (N <= vec_limit) without NaNs: Uses vectorized matrix indexing.
%      This is significantly faster for small N due to vectorized BLAS operations.
%   2. For large samples (N > vec_limit) or with NaNs: Uses the loop-based approach.
%      This avoids the O(N^3) memory/compute scaling of the matrix expansion,
%      which becomes a severe bottleneck at large N.
%
%   Both methods produce bit-identical results.
%
% Inputs:
%   x           - Column vector of the first sample.
%   y           - Column vector of the second sample.
%   metric_type - String/Char 'delta' for Cliff's Delta, 'rel' for Rel. Diff.
%                 Or boolean true for Delta, false for Rel Diff (internal optimization).
%   vec_limit   - (Optional) Threshold for N to switch algorithms. 
%                 Pass [] to use internal default (150).
%
% Outputs:
%   vals - Column vector of Jackknife statistics (length depends on valid input pairs).
%   a    - Acceleration factor (scalar) for BCa calculation.
%
% Author: Lukas von Erdmannsdorff

    % Input validation
    x = x(:); 
    y = y(:);
    N = numel(x);
    
    if numel(y) ~= N
        error('Input vectors must have the same length.');
    end
    
    vals = zeros(N, 1);
    
    % Determine metric type safely
    is_delta = false;
    if (ischar(metric_type) || isstring(metric_type))
        if strcmpi(metric_type, 'delta')
            is_delta = true;
        elseif strcmpi(metric_type, 'rel')
            is_delta = false;
        else
            error('Unknown metric type. Use ''delta'' or ''rel''.');
        end
    elseif islogical(metric_type) || isnumeric(metric_type)
         is_delta = logical(metric_type);
    end
    
    % --- Hybrid Algorithm Selection ---
    if nargin < 4 || isempty(vec_limit)
        vec_limit = 150;
    end
    MAX_VECTORIZED_SAMPLE_SIZE = vec_limit;
    
    % Check for NaN presence in source data
    has_nans = any(isnan(x)) || any(isnan(y));
    
    if N > MAX_VECTORIZED_SAMPLE_SIZE || has_nans
        %% Method A: Original Loop-Based Approach (O(N^2) - O(N^3))
        % Preferred for:
        % 1. Large samples (N > 150): Avoids O(N^3) memory allocation of matrix method.
        % 2. Data with NaNs: Pairwise exclusion requires variable sample sizes per iteration.
        
        for n = 1:N
            % Leave-one-out indexing
            % We create a logical mask to exclude the n-th element
            idx = true(N, 1); 
            idx(n) = false;
            
            jx = x(idx); 
            jy = y(idx);
            
            % Robust handling of NaNs (Pairwise Exclusion)
            % We only use pairs where both x and y are valid numbers
            valid = ~isnan(jx) & ~isnan(jy);
            jx_valid = jx(valid); 
            jy_valid = jy(valid);
            
            % Check if sufficient data remains
            if ~isempty(jx_valid)
                if is_delta
                    vals(n) = HERA.stats.cliffs_delta(jx_valid, jy_valid);
                else
                    vals(n) = HERA.stats.relative_difference(jx_valid, jy_valid);
                end
            else
                vals(n) = NaN; % Determine behavior for empty set
            end
        end
        
    else
        %% Method B: Vectorized Matrix Approach (Small/Medium Samples, No NaNs)
        % Efficient for N <= 150.
        % Strategy:
        % 1. Create an index matrix where each column represents a leave-one-out sample.
        % 2. Use matrix indexing to extract all jackknife samples simultaneously.
        % 3. Call the effect size function once with matrix inputs.
        
        % Build leave-one-out index matrix: (N-1) x N
        all_indices = repmat((1:N)', 1, N);
        diag_mask = ~eye(N, 'logical');
        loo_indices = reshape(all_indices(diag_mask), N-1, N);
        
        % Generate all leave-one-out samples: (N-1) x N matrices
        x_loo = x(loo_indices);
        y_loo = y(loo_indices);
        
        % Call effect size function with matrix inputs (vectorized mode)
        % The function handles N x B matrices and returns 1 x N output.
        if is_delta
            vals = HERA.stats.cliffs_delta(x_loo, y_loo);
        else
            vals = HERA.stats.relative_difference(x_loo, y_loo);
        end
        
        % Ensure column vector output (functions return row for matrix input)
        vals = vals(:);
    end
    
    % Post-processing: remove NaNs from the result to return a clean vector of statistics
    % This is standard for BCa acceleration calculation (a) where we need the moments of valid stats.
    vals = vals(~isnan(vals));
    
    if isempty(vals)
        vals = 0; % Safety fallback to avoid empty returns causing issues in mean/sum
    end

    % --- Optional: Calculate BCa Acceleration Factor (a) ---
    if nargout > 1
        mean_jack = mean(vals);
        a_num = sum((mean_jack - vals).^3);
        a_den = 6 * (sum((mean_jack - vals).^2)).^(3/2);
        
        if a_den == 0
            a = 0;
        else
            a = a_num / a_den;
        end
        
        if ~isfinite(a)
            a = 0;
        end
    end
end
