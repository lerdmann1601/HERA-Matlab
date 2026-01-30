function save_csv(results, modes, out_dir, ts_str)
% SAVE_CSV - Saves the raw simulation results of the convergence analysis to CSV files.
%
% Syntax:
%   HERA.analysis.convergence.save_csv(results, modes, out_dir, ts_str)
%
% Description:
%   Iterates through all simulated scenarios and saves the detailed results
%   (Error, Cost, Failure) for each metric (Threshold, BCa, Ranking) and 
%   parameter mode (Relaxed, Default, Strict) into structured CSV files.
%
% Inputs:
%   results  - Struct array containing the simulation data for all scenarios.
%   modes    - Cell array of mode names (e.g., {'Relaxed', 'Default', 'Strict'}).
%   out_dir  - Directory where the CSV files will be saved.
%   ts_str   - Timestamp string to ensure unique filenames.
%
% Author: Lukas von Erdmannsdorff

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
            fprintf(global_fid, 'Scenario,Metric,Mode,Mean_Error_Percent,Std_Error,Mean_Cost_B,Failure_Rate_Percent\n');
        end
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
                    for sim_idx = 1:n_sims
                        val_err = err_mat(sim_idx, mode_idx);
                        val_cost = cost_mat(sim_idx, mode_idx);
                        val_fail = fail_mat(sim_idx, mode_idx);
                        
                        fprintf(fileID, '%d,%s,%s,%.4f,%d,%d\n', ...
                            sim_idx, m_name, mode_name, val_err, val_cost, val_fail);
                    end
                    
                    % --- B. Calculate & Write Aggregated Stats to Global Summary ---
                    if global_fid ~= -1
                        mean_err = mean(err_mat(:, mode_idx), 'omitnan');
                        std_err = std(err_mat(:, mode_idx), 'omitnan');
                        mean_cost = mean(cost_mat(:, mode_idx), 'omitnan');
                        fail_rate = (sum(fail_mat(:, mode_idx)) / n_sims) * 100;
                        
                        fprintf(global_fid, '%s,%s,%s,%.4f,%.4f,%.1f,%.1f\n', ...
                            safe_name, m_name, mode_name, mean_err, std_err, mean_cost, fail_rate);
                    end
                end
            end
            
            fclose(fileID);
        end
        
        % Close Global CSV
        if global_fid ~= -1
            fclose(global_fid);
        end
        
    catch ME
        if exist('fileID', 'var') && fileID ~= -1, fclose(fileID); end
        if exist('global_fid', 'var') && global_fid ~= -1, fclose(global_fid); end
        fprintf('Error saving CSV results: %s\n', ME.message);
    end
end
