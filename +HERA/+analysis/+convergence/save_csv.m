function save_csv(command, varargin)
% SAVE_CSV - Saves the raw simulation results of the convergence analysis to CSV files.
%
% Syntax:
%   HERA.analysis.convergence.save_csv('config', modes, scenarios, params, out_dir, ts_str)
%   HERA.analysis.convergence.save_csv('results', results, modes, out_dir, ts_str)
%
% Description:
%   Handles writing configuration setups and raw simulation results to structured CSV files.
%
% Author: Lukas von Erdmannsdorff

    if isstruct(command)
        % Backward compatibility for calls without the 'results' string prefix
        save_results_csv(command, varargin{:});
        return;
    end

    switch command
        case 'config'
            save_config_csv(varargin{:});
        case 'results'
            save_results_csv(varargin{:});
        otherwise
            error('Unknown save_csv command: %s', command);
    end
end

function save_config_csv(modes, scenarios, params, out_dir, ts_str)
    % Creates Configuration CSVs (Scenarios & Methods)
    scenarios_filename = fullfile(out_dir, sprintf('Scenarios_%s.csv', ts_str));
    if ~exist(scenarios_filename, 'file')
        try
            s_fid = fopen(scenarios_filename, 'w');
            if s_fid ~= -1
                fprintf(s_fid, 'Index,Scenario_Name,N,Distribution,Data_Summary\n');
                fclose(s_fid);
                
                s_T = cell(length(scenarios), 5);
                for i = 1:length(scenarios)
                    sc = scenarios(i);
                    clean_name = strrep(sc.name, ',', ';');
                    clean_dist = strrep(sc.Dist, ',', ';');
                    clean_summary = strrep(sc.DataSummary, ',', ';');
                    s_T(i, :) = {i, clean_name, sc.N, clean_dist, clean_summary};
                end
                writetable(cell2table(s_T), scenarios_filename, 'Delimiter', ',', 'WriteMode', 'Append', 'WriteVariableNames', false, 'QuoteStrings', true);
            end
            
            p_filename = fullfile(out_dir, sprintf('Method_Parameters_%s.csv', ts_str));
            p_fid = fopen(p_filename, 'w');
            if p_fid ~= -1
                mode_str = strjoin(modes, ',');
                fprintf(p_fid, 'Metric,Parameter,%s\n', mode_str);
                fclose(p_fid);
                
                param_metrics = {'Bootstrap Thresholds', 'BCa Confidence Intervals', 'Ranking Stability'};
                param_structs = {params.thr, params.bca, params.rnk};
                lbls = {'Trials (n)', 'Smoothing (sm)', 'Streak (st)', 'Tolerance (tol)', 'Step Size (B)', 'B Range'};
                fields = {'n', 'sm', 'st', 'tol', 'step', 'range'};
                
                num_rows = length(param_metrics) * length(lbls);
                num_cols = 2 + length(modes);
                p_T = cell(num_rows, num_cols);
                row_idx = 1;
                
                for m_idx = 1:length(param_metrics)
                    m_name = param_metrics{m_idx};
                    p_cell = param_structs{m_idx};
                    
                    for r = 1:length(lbls)
                        p_T{row_idx, 1} = m_name;
                        p_T{row_idx, 2} = lbls{r};
                        for c = 1:length(p_cell)
                            val = p_cell{c};
                            if strcmp(fields{r}, 'range')
                                str = sprintf('[%d-%d]', val.start, val.end);
                            else
                                v = val.(fields{r});
                                if abs(v) < 1 && v > 0
                                     str = sprintf('%.3f', v);
                                else
                                     str = sprintf('%d', v);
                                end
                            end
                            p_T{row_idx, 2+c} = str;
                        end
                        row_idx = row_idx + 1;
                    end
                end
                writetable(cell2table(p_T), p_filename, 'Delimiter', ',', 'WriteMode', 'Append', 'WriteVariableNames', false, 'QuoteStrings', true);
            end
        catch ME_cfg
            warning('Could not write configuration CSVs: %s', ME_cfg.message);
        end
    end
end

function save_results_csv(results, modes, out_dir, ts_str)
    % Writes detailed and aggregated results for each scenario
    arguments
        results (:,1) struct
        modes (1,:) cell
        out_dir (1,1) string
        ts_str (1,1) string
    end

    % Ensure base output directory exists
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end
    
    % 1. Setup Folders
    % Create a single subfolder for all detailed results
    results_subfolder = fullfile(out_dir, sprintf('Results_%s', ts_str));
    if ~exist(results_subfolder, 'dir'), mkdir(results_subfolder); end
    
    % 2. Create Global Summary CSV (Aggregated)
    % This file will stay in the main CSV folder for easy access
    global_filename = fullfile(out_dir, sprintf('Global_Summary_%s.csv', ts_str));
    
    % Check if file exists to determine if we need a header
    write_header = ~exist(global_filename, 'file');
    
    global_fid = fopen(global_filename, 'a'); % Open in Append mode
    if global_fid ~= -1
        if write_header
            % Header for Aggregated Data
            fprintf(global_fid, 'Scenario,Metric,Mode,Median_Error_Percent,IQR_Error,Error_Q1,Error_Q3,CI95_Lower,CI95_Upper,Median_Cost_B,IQR_Cost_B,Cost_Q1,Cost_Q3,Cost_CI95_Lower,Cost_CI95_Upper,Failure_Rate_Percent\n');
        end
        fclose(global_fid);
    else
        warning('Could not create/open Global Summary CSV: %s', global_filename);
    end

    metric_names = {'Thresholds', 'BCa', 'Ranking'};
    metric_fields = {'thr', 'bca', 'rnk'};
    
    try
        % Iterate over each scenario
        for s = 1:length(results)
            sc = results(s);
            
            % Generate Safe Scenario Name
            safe_name = regexprep(sc.name, '[^a-zA-Z0-9]+', '_'); 
            safe_name = regexprep(safe_name, '_$', ''); 
            
            % Scenario Filename (Saved in the Results_timestamp subfolder)
            filename = fullfile(results_subfolder, sprintf('%s.csv', safe_name));
            
            fileID = fopen(filename, 'w');
            if fileID == -1
                warning('Could not open file for writing: %s', filename);
                continue;
            end

            % Write Header for Detailed Scenario File (Raw Data)
            fprintf(fileID, 'SimID,Metric,Mode,Error,Cost,Fail\n');
            fclose(fileID);
            
            % Iterate over each metric type
            for m_id = 1:length(metric_fields)
                field = metric_fields{m_id};
                m_name = metric_names{m_id};
                
                % Get data matrices: [Sims x Modes]
                err_mat = sc.(field).err;
                cost_mat = sc.(field).cost;
                fail_mat = sc.(field).fail;
                
                [n_sims, n_modes] = size(err_mat);
                
                % Process each mode
                for mode_idx = 1:n_modes
                    % Map mode index to name
                    if mode_idx <= length(modes)
                        mode_name = modes{mode_idx};
                    else
                        mode_name = sprintf('Mode_%d', mode_idx);
                    end
                    
                    % --- A. Write Raw Data to Scenario File ---
                    t_T_raw = cell(n_sims, 6);
                    for sim_idx = 1:n_sims
                        val_err = err_mat(sim_idx, mode_idx);
                        val_cost = cost_mat(sim_idx, mode_idx);
                        val_fail = fail_mat(sim_idx, mode_idx);
                        
                        t_T_raw(sim_idx, :) = {sim_idx, m_name, mode_name, round(val_err, 4), val_cost, val_fail};
                    end
                    writetable(cell2table(t_T_raw), filename, 'Delimiter', ',', 'WriteMode', 'Append', 'WriteVariableNames', false, 'QuoteStrings', true);
                    
                    % --- B. Calculate & Write Aggregated Stats to Global Summary ---
                    if global_fid ~= -1
                        % Calculate Median, IQR, and 95% CI for Errors (approximate normal)
                        err_vals = err_mat(:, mode_idx);
                        err_vals = err_vals(~isnan(err_vals)); % Remove NaNs for stat calculations
                        
                        cost_vals = cost_mat(:, mode_idx);
                        cost_vals = cost_vals(~isnan(cost_vals));
                        
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
                        
                        fail_rate = (sum(fail_mat(:, mode_idx)) / n_sims) * 100;
                        
                        row_data = {safe_name, m_name, mode_name, round(median_err, 4), round(iqr_err, 4), round(err_q1, 4), ...
                            round(err_q3, 4), round(ci_lower, 4), round(ci_upper, 4), round(median_cost, 1), round(iqr_cost, 1), ...
                            round(cost_q1, 1), round(cost_q3, 1), round(cost_ci_lower, 1), round(cost_ci_upper, 1), round(fail_rate, 1)};
                        writetable(cell2table(row_data), global_filename, 'Delimiter', ',', 'WriteMode', 'Append', 'WriteVariableNames', false, 'QuoteStrings', true);
                    end
                end
            end
        end
        
    catch ME
        fprintf('Error saving CSV results: %s\n', ME.message);
    end
end
