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
            % Apply moving average to the entire history to reduce noise
            % 'omitnan' ensures single NaNs don't break the window
            smoothed_full = movmean(stability_history, window, 'omitnan');
            
            curr_smooth = smoothed_full(end);
            prev_smooth = smoothed_full(end-1);
            
            % Calculate relative improvement with safety checks
            rel_imp = calculate_relative_improvement(prev_smooth, curr_smooth);
            
            % Update stats output
            stats.improvement = rel_imp;
            stats.smoothed_value = curr_smooth;

            % Calculate the current "streak" of stability.
            % Since this function is stateless, we look backwards from the current step 'n'
            % to count how many consecutive steps (back down to the warmup period) 
            % have satisfied the tolerance condition.
            current_streak = 0;
            
            % Loop backwards from current step
            for k = n:-1:(warmup + window)
                s_curr = smoothed_full(k);
                s_prev = smoothed_full(k-1);
                
                % Re-calculate improvement for past steps
                imp = calculate_relative_improvement(s_prev, s_curr);
                
                % Check tolerance condition
                if abs(imp) < cfg.convergence_tolerance
                    current_streak = current_streak + 1;
                else
                    % Streak is broken
                    break;
                end
            end
            stats.streak = current_streak;
            
            % Check if the streak requirement is satisfied
            if current_streak >= streak_needed
                is_converged = true;
            end
        end
    else
        % Simple convergence check (No smoothing, immediate check)
        if n >= cfg.min_steps_for_convergence_check
            curr_stab = stability_history(end);
            prev_stab = stability_history(end-1);
            
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

    if isnan(curr) || isnan(prev)
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
