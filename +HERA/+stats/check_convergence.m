function [is_converged, stats] = check_convergence(stability_history, config)
% CHECK_CONVERGENCE - Checks if the stability curve has converged using robust or simple criteria.
%
% Syntax:
%   [is_converged, stats] = check_convergence(stability_history, config)
%
% Description:
%   This function analyzes the history of stability values (e.g., dispersion of CI widths)
%   to determine if the bootstrap process has stabilized. It supports two modes:
%   1. Robust Mode (Smoothing): Applies a moving average to smooth out noise and requires
%      a consistent "streak" of small relative changes to declare convergence.
%   2. Simple Mode: Checks the immediate relative change between the last two points.
%
%   The function handles edge cases such as NaN, Infinite, or Zero stability values
%   to prevent runtime errors and ensure reliable convergence detection.
%
% Inputs:
%   stability_history - Vector of stability values calculated up to the current bootstrap step.
%                       (Format: [1 x N] double)
%   config            - Configuration struct containing convergence parameters:
%                       * config.min_steps_for_convergence_check: Minimum steps before checking.
%                       * config.convergence_tolerance: Threshold for relative improvement (e.g., 0.01 for 1%).
%                       * (Optional for Robust) config.smoothing_window: Window size for `movmean`.
%                       * (Optional for Robust) config.convergence_streak_needed: Number of consecutive stable steps required.
%
% Outputs:
%   is_converged      - Boolean flag. True if convergence criteria are met, False otherwise.
%   stats             - Struct containing debugging statistics:
%                       * stats.improvement: The calculated relative improvement (signed).
%                       * stats.streak: Current number of consecutive stable steps (Robust mode only).
%                       * stats.smoothed_value: The smoothed (or raw) stability value used for the check.
%
% logic:
%   - For Robust Mode:
%     1. Applies a moving average (sliding window) to the full history.
%     2. Calculates relative improvement between the current and previous smoothed values.
%     3. checks backwards through history to count how many consecutive steps have met the tolerance.
%   - For Simple Mode:
%     1. Calculates relative improvement between the last two raw values.
%     2. Checks if improvement is within tolerance.
%
% Author: Lukas von Erdmannsdorff

    % Default output initialization
    is_converged = false;
    stats = struct('improvement', NaN, 'streak', 0, 'smoothed_value', NaN);

    n = length(stability_history);
    cfg = config;
    
    % Determine if robust convergence logic (w/ smoothing) is requested based on config fields
    use_robust = isfield(cfg, 'smoothing_window') && ~isempty(cfg.smoothing_window) && ...
                 isfield(cfg, 'convergence_streak_needed') && ~isempty(cfg.convergence_streak_needed);

    if use_robust
        warmup = cfg.min_steps_for_convergence_check;
        window = cfg.smoothing_window;
        streak_needed = cfg.convergence_streak_needed;
        
        % Only proceed if enough data points exist (Warmup + Window)
        if n >= warmup + window
            
        
            %% Robust Convergence Logic (Vectorized Trailing Window)
            % Strategy:
            % We use a specific Trailing Moving Average [window-1, 0].
            % This ensures that the smoothed value at index 'i' acts purely on
            % past data (i-(w-1) to i). This is crucial for:
            % -> Kausal integrity (no influence from future data).
            % -> Invariance (historical stability does not change with new data).
            
            % 1. Apply efficient vectorized smoothing to the entire history.
            % 'omitnan' handles initial steps gracefully.
            smoothed_full = movmean(stability_history, [window-1, 0], 'omitnan');
            
            % 2. Extract relevant smoothed values for improvement calculation
            %    We need to check the streak ending at the current step 'n'.
            %    We look back 'streak_needed' steps.
            
            % Indices to check: From n down to (n - streak + 1)
            % For each check at index k, we compare smooth(k) vs smooth(k-1).
            
            % Create index vectors for the 'current' and 'previous' values in the streak
            % We need at least 'streak_needed' comparisons.
            % e.g. if streak=3, we need pairs: (n, n-1), (n-1, n-2), (n-2, n-3).
            
            check_indices = n : -1 : max(warmup + window, n - streak_needed + 1);
            
            if isempty(check_indices)
                stats.streak = 0;
                return;
            end
            
            % Vectorized extraction of smoothed values
            curr_vals = smoothed_full(check_indices);
            prev_vals = smoothed_full(check_indices - 1);
            
            % 3. Calculate Relative Improvements (Vectorized)
            %    Handle zeros and infinities element-wise.
            
            % Pre-allocate
            rel_imps = zeros(size(curr_vals));
            
            % Logic matches calculate_relative_improvement but vectorized
            % Note: Since we are inside a heavily called function, inline vectorization is faster.
            
            mask_valid = ~isnan(curr_vals) & ~isnan(prev_vals) & ~isinf(curr_vals) & ~isinf(prev_vals) & (prev_vals ~= 0);
            
            % Standard case
            rel_imps(mask_valid) = (prev_vals(mask_valid) - curr_vals(mask_valid)) ./ prev_vals(mask_valid);
            
            % Edge cases (simplified for speed, matching helper logic broadly)
            % If prev=0 and curr=0 -> 0 change. If prev=0 and curr!=0 -> Infinite change (-1 or +1)
            rel_imps(prev_vals == 0 & curr_vals == 0) = 0;
            rel_imps(prev_vals == 0 & curr_vals ~= 0) = -1.0; % High change
            rel_imps(isinf(curr_vals) | isinf(prev_vals)) = 1.0; % Infinite -> unstable
            
            % 4. Check Tolerance
            is_stable = abs(rel_imps) < cfg.convergence_tolerance;
            
            % 5. Determine Streak
            %    Count consecutive true values starting from the first element (which corresponds to step 'n')
            %    find the first index where is_stable is false.
            first_fail = find(~is_stable, 1, 'first');
            
            if isempty(first_fail)
                % All checked values are stable
                current_streak = numel(is_stable);
            else
                % Streak is up to the failure ( indices 1 to first_fail-1 )
                current_streak = first_fail - 1;
            end
            
            % 6. Update Output Stats
            stats.streak = current_streak;
            stats.smoothed_value = smoothed_full(n);
            stats.improvement = rel_imps(1); % Improvement of the most recent step (n vs n-1)
            
            % 7. Final Decision
            if current_streak >= streak_needed
                is_converged = true;
            end
        end
    else
        %% Simple Convergence Logic 
        % Simple convergence check (No smoothing, immediate check)
        % We strictly need at least 2 data points to calculate a change.
        if n >= max(2, cfg.min_steps_for_convergence_check)
            % Ensure scalar values for stability
            curr_stab = stability_history(end);
            if numel(curr_stab) > 1, curr_stab = curr_stab(1); end
            
            prev_stab = stability_history(end-1);
            if numel(prev_stab) > 1, prev_stab = prev_stab(1); end
            
            % Calculate relative improvement on raw values
            rel_imp = calculate_relative_improvement(prev_stab, curr_stab);
            
            stats.improvement = rel_imp;
            stats.smoothed_value = curr_stab; % No smoothing applied
            stats.streak = 0; % Not applicable
            
            % If absolute improvement is below tolerance, we converge immediately
            if abs(rel_imp) < cfg.convergence_tolerance
                is_converged = true;
            end
        end
    end
end

function rel_imp = calculate_relative_improvement(prev, curr)
% CALCULATE_RELATIVE_IMPROVEMENT - Helper to calculate improvement safely.
% Handles NaN, Inf, and Zero cases to avoid runtime errors or misleading results.

    if any(isnan(curr), 'all') || any(isnan(prev), 'all')
        % If data is missing, assume unstable/high change (100%) to force continuation
        rel_imp = 1.0;
    elseif isinf(curr)
        % If current became infinite (instability), return -1.0 (worsened)
        rel_imp = -1.0;
    elseif prev == 0
        % Prevent division by zero
        if curr == 0
            rel_imp = 0; % No change
        else
            rel_imp = -1.0; % From zero to non-zero is technically infinite worsening
        end
    elseif isinf(prev)
        % From infinite to something else?
        if isinf(curr)
            rel_imp = 0; % Still infinite, no change
        else
            rel_imp = 1.0; % From infinite to finite is 100% improvement
        end
    else
        % Standard relative improvement calculation: (Old - New) / Old
        % Positive value means stability decreased (value got smaller, which is good for IQR/Median)
        % Note: Stability metric is IQR/Median or IQR, so lower is usually "more stable" or "tighter".
        % However, "convergence" here means the value stops CHANGING, not necessarily that it goes to zero.
        % The generic formula (prev - curr) / prev checks for decay.
        % But actually we take ABS(rel_imp) < tol in the check, so direction doesn't strictly matter for the 'stop' trigger,
        % only magnitude of change.
        rel_imp = (prev - curr) / prev;
    end
end
