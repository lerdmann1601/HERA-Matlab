function vals = jackknife(x, y, metric_type)
% JACKKNIFE - Calculates Jackknife statistics for Cliff's Delta or Relative Difference.
%
% Syntax:
%   vals = HERA.stats.jackknife(x, y, metric_type)
%
% Description:
%   Computes the Jackknife (leave-one-out) statistics for a given metric.
%   It iteratively removes one subject (pair) and calculates the effect size
%   on the remaining N-1 subjects.
%
%   This function strictly handles missing data (NaNs) by performing pairwise exclusion
%   within each Jackknife subsample. It is optimized for performance by avoiding
%   redundant memory allocations where possible.
%
% Inputs:
%   x           - Column vector of the first sample.
%   y           - Column vector of the second sample.
%   metric_type - String/Char 'delta' for Cliff's Delta, 'rel' for Rel. Diff.
%                 Or boolean true for Delta, false for Rel Diff (internal optimization).
%
% Outputs:
%   vals - Column vector of Jackknife statistics (length depends on valid input pairs).
%          If input has N elements, output has N elements (some might be NaN if N-1 sample is invalid).
%          Note: The function returns clean values (NaNs removed from result) to match
%          BCa requirements, or zeros if empty.
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
    
    % Post-processing: remove NaNs from the result to return a clean vector of statistics
    % This is standard for BCa acceleration calculation (a) where we need the moments of valid stats.
    vals = vals(~isnan(vals));
    
    if isempty(vals)
        vals = 0; % Safety fallback to avoid empty returns causing issues in mean/sum
    end
end
