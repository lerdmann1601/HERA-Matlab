function calc_pooled_csv(results, modes, out_dir, ts_str)
% CALC_POOLED_CSV - Calculates and saves pooled convergence stats across all scenarios.
%
% Syntax:
%   HERA.analysis.convergence.calc_pooled_csv(results, modes, out_dir, ts_str)
%
% Description:
%   Aggregates raw simulation data (Errors, Costs, Fails) from all tested 
%   scenarios into a single global summary CSV file for final reporting.
%
% Author: Lukas von Erdmannsdorff


    if isempty(results)
        return;
    end
    
    pooled_filename = fullfile(out_dir, sprintf('Pooled_Results_%s.csv', ts_str));
    pooled_fid = fopen(pooled_filename, 'w');
    
    if pooled_fid == -1
        warning('Could not create/open Pooled Results CSV: %s', pooled_filename);
        return;
    end
    
    fprintf(pooled_fid, 'Metric,Mode,Median_Error_Percent,IQR_Error,Error_Q1,Error_Q3,CI95_Lower,CI95_Upper,Median_Cost_B,IQR_Cost_B,Cost_Q1,Cost_Q3,Cost_CI95_Lower,Cost_CI95_Upper,Failure_Rate_Percent\n');
    fclose(pooled_fid);
    
    metric_names = {'Thresholds', 'BCa', 'Ranking'};
    metric_fields = {'thr', 'bca', 'rnk'};
    
    for m_id = 1:length(metric_fields)
        field = metric_fields{m_id};
        m_name = metric_names{m_id};
        
        [~, n_modes] = size(results(1).(field).err);
        
        % Pre-allocate arrays for RAM efficiency
        n_scenarios = length(results);
        n_sims_per_scenario = size(results(1).(field).err, 1);
        total_sims = n_scenarios * n_sims_per_scenario;
        
        for mode_idx = 1:n_modes
            if mode_idx <= length(modes)
                mode_name = modes{mode_idx};
            else
                mode_name = sprintf('Mode_%d', mode_idx);
            end
            
            pooled_err = NaN(total_sims, 1);
            pooled_cost = NaN(total_sims, 1);
            pooled_fail = false(total_sims, 1);
            
            % Pool all data across all scenarios
            idx_start = 1;
            for s = 1:n_scenarios
                n_current_sims = size(results(s).(field).err, 1);
                idx_end = idx_start + n_current_sims - 1;
                
                pooled_err(idx_start:idx_end) = results(s).(field).err(:, mode_idx);
                pooled_cost(idx_start:idx_end) = results(s).(field).cost(:, mode_idx);
                pooled_fail(idx_start:idx_end) = logical(results(s).(field).fail(:, mode_idx));
                
                idx_start = idx_end + 1;
            end
            
            % Calculate Stats
            err_vals = pooled_err(~isnan(pooled_err));
            cost_vals = pooled_cost(~isnan(pooled_cost));
            
            iqr_err = iqr(err_vals);
            if ~isempty(err_vals)
                 median_err= median(err_vals);
                 err_q1   = prctile(err_vals, 25);
                 err_q3   = prctile(err_vals, 75);
                 ci_lower = prctile(err_vals, 2.5);
                 ci_upper = prctile(err_vals, 97.5);
            else
                 median_err= NaN;
                 err_q1   = NaN;
                 err_q3   = NaN;
                 ci_lower = NaN;
                 ci_upper = NaN;
            end
            
            iqr_cost = iqr(cost_vals);
            if ~isempty(cost_vals)
                 median_cost   = median(cost_vals);
                 cost_q1       = prctile(cost_vals, 25);
                 cost_q3       = prctile(cost_vals, 75);
                 cost_ci_lower = prctile(cost_vals, 2.5);
                 cost_ci_upper = prctile(cost_vals, 97.5);
            else
                 median_cost   = NaN;
                 cost_q1       = NaN;
                 cost_q3       = NaN;
                 cost_ci_lower = NaN;
                 cost_ci_upper = NaN;
            end
            
            valid_sims = length(pooled_fail);
            if valid_sims > 0
                fail_rate = (sum(pooled_fail) / valid_sims) * 100;
            else
                fail_rate = NaN;
            end
            
            % Form the single row as a cell array and append using writetable
            row_data = {m_name, mode_name, round(median_err, 4), round(iqr_err, 4), round(err_q1, 4), ...
                round(err_q3, 4), round(ci_lower, 4), round(ci_upper, 4), round(median_cost, 1), ...
                round(iqr_cost, 1), round(cost_q1, 1), round(cost_q3, 1), round(cost_ci_lower, 1), ...
                round(cost_ci_upper, 1), round(fail_rate, 1)};
            
            T = cell2table(row_data);
            writetable(T, pooled_filename, 'Delimiter', ',', 'WriteMode', 'Append', 'WriteVariableNames', false, 'QuoteStrings', true);
        end
    end
end
