function results = convergence_analysis(n_sims_per_cond, log_path_or_mode)
% CONVERGENCE_ANALYSIS - Scientific validation of bootstrap and ranking parameters.
%
% Syntax:
%   results = HERA.analysis.convergence_analysis(n_sims_per_cond)
%   results = HERA.analysis.convergence_analysis(n_sims_per_cond, log_path_or_mode)
%
% Description:
%   This function performs a Monte Carlo study to validate the robustness and accuracy
%   of the HERA bootstrap algorithms (Thresholds, BCa CI, Ranking). It simulates various 
%   data scenarios (Normal, LogNormal, Likert, Bimodal) and compares three different 
%   parameter configurations (Relaxed, Default, Strict) against a high-precision 
%   Reference ("Truth").
%
% Workflow:
%   1. User Input & Initialization:
%      - Determines output location via Auto-Log logic or User Dialog.
%   2. Configuration:
%      - Defines data scenarios (N=25/50/100, different distributions).
%      - Defines parameter sets for Thresholds, BCa, and Ranking.
%   3. Setup Environment:
%      - Initializes the parallel pool and temporary directories.
%   4. Simulation Loop:
%      - Generates synthetic data for each simulation/scenario.
%      - Calculates the "Truth" (Reference) using very high bootstrap counts (B).
%      - Iteratively tests the three parameter configurations.
%      - Resets random seeds before every step to ensure strict reproducibility.
%   5. Reporting:
%      - Generates scientific plots for every scenario.
%      - Exports detailed simulation results to structured CSV files.
%      - Compiles a global summary PDF report.
%
% Inputs:
%   n_sims_per_cond  - (Optional) Number of simulations per scenario (Default: 15).
%                      Higher values (e.g., 50-100) yield more robust statistics.
%                      Alternatively, this can be a path to a .json configuration file.
%   log_path_or_mode - (Optional) Output configuration (Default: "" -> Auto-Log).
%                      - ""            : Auto-Log to 'Documents/HERA_convergence_Log'.
%                      - "interactive" : Opens folder selection dialog.
%                      - [Path String] : Uses the specified custom path.
%                      Note: Overridden in JSON mode by `output_dir`.
%
% Outputs:
%   results          - Struct containing the raw simulation data (Errors, Costs, Failures)
%                      for all scenarios and methods.
%
% Author: Lukas von Erdmannsdorff

    arguments
        n_sims_per_cond = 15
        log_path_or_mode (1,1) string = ""
    end

    clc;
    
    %% Parse JSON Config if provided
    customConfig = struct();
    if ischar(n_sims_per_cond) || isstring(n_sims_per_cond)
        configFile = string(n_sims_per_cond);
        if endsWith(configFile, ".json", "IgnoreCase", true)
            if exist(configFile, 'file')
                try
                    json_text = fileread(configFile);
                    loadedData = jsondecode(json_text);
                    if isfield(loadedData, 'userInput')
                        customConfig = loadedData.userInput;
                    else
                        customConfig = loadedData;
                    end
                    
                    % Extract n_sims_per_cond 
                    if isfield(customConfig, 'n_sims_per_cond') && isnumeric(customConfig.n_sims_per_cond)
                        n_sims_per_cond = customConfig.n_sims_per_cond;
                    else
                        n_sims_per_cond = 15;
                    end
                    
                    % Extract output_dir
                    if isfield(customConfig, 'output_dir') && (ischar(customConfig.output_dir) || isstring(customConfig.output_dir))
                        log_path_or_mode = string(customConfig.output_dir);
                    end
                catch ME
                    error('Error reading configuration file: %s', ME.message);
                end
            else
                error('Configuration file not found: %s', configFile);
            end
        else
            error('First argument must be numeric (n_sims_per_cond) or a path to a .json config file.');
        end
    else
        % Validate n_sims_per_cond as numeric scalar
        if ~isnumeric(n_sims_per_cond) || ~isscalar(n_sims_per_cond) || n_sims_per_cond <= 0 || mod(n_sims_per_cond, 1) ~= 0
            error('n_sims_per_cond must be a positive integer.');
        end
    end
    
    %% 1. User Input & Initialization
    % Handle arguments for custom path or mode
    if log_path_or_mode == ""
        % Default: Auto-Log Mode
        [out_dir, path_source] = get_writable_convergence_log_path();
        fprintf('Auto-Log Mode enabled. Saving to: %s (%s)\n', out_dir, path_source);
    elseif strcmpi(log_path_or_mode, 'interactive')
        % Interactive Mode
        out_dir = uigetdir(pwd, 'Select Folder to Save Robustness Reports');
        if isequal(out_dir, 0)
            fprintf('Cancelled by user.\n');
            results = [];
            return;
        end
    else
        % Custom Path Mode
        if ~exist(log_path_or_mode, 'dir')
            mkdir(log_path_or_mode);
        end
        out_dir = log_path_or_mode;
        fprintf('Custom Log Path enabled. Saving to: %s\n', out_dir);
    end
    
    base_name = 'Robustness_Report';
    t_start = tic;

    import HERA.analysis.convergence.*

    %% 2. Configuration
    % Retrieve all config parameters
    [n_datasets, modes, scenarios, params, refs, limits, cfg_base, colors, ram_gb] = config(n_sims_per_cond, customConfig);

    %% 3. Setup Environment
    temp_dir = tempname; 
    if ~exist(temp_dir, 'dir'), mkdir(temp_dir); end
    cleanTemp = onCleanup(@() rmdir(temp_dir, 's'));
    
    try
        lang = HERA.get_language();
        styles = HERA.design('light', 2, true); 
    catch
        lang = struct('calculating', 'Calc', 'done', 'Done'); 
        styles = struct(); 
    end

    % Output Directories setup
    ts_str = string(datetime('now'), 'yyyyMMdd_HHmmss');
    timestamp_folder = [base_name, '_', char(ts_str)];
    final_out_dir = fullfile(out_dir, timestamp_folder);
    if ~exist(final_out_dir, 'dir'), mkdir(final_out_dir); end
    
    dir_graphics = fullfile(final_out_dir, 'Graphics'); 
    dirty_pdfs = fullfile(final_out_dir, 'PDF');
    if ~exist(dir_graphics, 'dir'), mkdir(dir_graphics); end
    if ~exist(dirty_pdfs, 'dir'), mkdir(dirty_pdfs); end
    
    dir_csv = fullfile(final_out_dir, 'CSV');
    if ~exist(dir_csv, 'dir'), mkdir(dir_csv); end
    
    % Start logging terminal output to file
    log_filename = fullfile(final_out_dir, ['conv_log_', char(ts_str), '.txt']);
    diary(log_filename);
    cleanupDiary = onCleanup(@() diary('off'));  % Ensures diary closes on ANY exit (incl. Ctrl+C)

    pdf_full = fullfile(final_out_dir, ['Full_Combined_Report_', char(ts_str), '.pdf']);
    
    % Print Header and Config Summary to Console & Log
    fprintf('\n==========================================================\n');
    fprintf('   Scientific Bootstrap Robustness Study (Sims/Cond=%d)\n', n_sims_per_cond);
    fprintf('==========================================================\n');
    if isfield(customConfig, 'target_memory') && isnumeric(customConfig.target_memory)
        fprintf(' Memory Limit:          %d MB (User Override)\n', cfg_base.system.target_memory);
    else
        fprintf(' Detected RAM:          %.1f GB\n', ram_gb);
        fprintf(' Memory Limit:          %d MB\n', cfg_base.system.target_memory);
    end
    if isfield(cfg_base, 'simulation_seed'), sseed = cfg_base.simulation_seed; else, sseed = 123; end
    if isfield(cfg_base, 'scenario_seed_offset'), sce_o = cfg_base.scenario_seed_offset; else, sce_o = 10000; end
    if isfield(cfg_base, 'reference_seed_offset'), ref_o = cfg_base.reference_seed_offset; else, ref_o = 1; end
    fprintf(' Simulation Base Seed:  %d\n', sseed);
    fprintf(' Bootstrap Offset:      %d\n', cfg_base.bootstrap_seed_offset);
    fprintf(' Scenario Offset:       %d\n', sce_o);
    fprintf(' Reference Offset:      %d\n', ref_o);
    fprintf('==========================================================\n\n');
    
    % Overview Plot (Created BEFORE simulation to free RAM for batch processing)
    plot('parameter_overview', params, scenarios, modes, refs, dir_graphics, char(ts_str), pdf_full, dirty_pdfs);
    
    % Write Configuration CSVs (Created BEFORE simulation)
    save_csv('config', modes, scenarios, params, dir_csv, ts_str);
    
    % Save JSON Config
    try
        json_file = fullfile(final_out_dir, ['configuration_', char(ts_str), '.json']);
        
        cfg_out = struct();
        cfg_out.n_sims_per_cond = n_sims_per_cond;
        cfg_out.output_dir = log_path_or_mode;
        
        % Ensure we save the effective configuration applied
        cfg_out.target_memory = cfg_base.system.target_memory;
        if isfield(cfg_base, 'simulation_seed')
            cfg_out.simulation_seed = cfg_base.simulation_seed;
        else
            cfg_out.simulation_seed = 123; 
        end
        cfg_out.bootstrap_seed_offset = cfg_base.bootstrap_seed_offset;
        
        if isfield(cfg_base, 'scenario_seed_offset')
            cfg_out.scenario_seed_offset = cfg_base.scenario_seed_offset;
        else
            cfg_out.scenario_seed_offset = 10000;
        end
        
        if isfield(cfg_base, 'reference_seed_offset')
            cfg_out.reference_seed_offset = cfg_base.reference_seed_offset;
        else
            cfg_out.reference_seed_offset = 1;
        end
        
        if exist('customConfig', 'var') && isstruct(customConfig)
            % Ensure we don't duplicate existing top-level fields
            fnames = fieldnames(customConfig);
            for i=1:length(fnames)
                cfg_out.(fnames{i}) = customConfig.(fnames{i});
            end
        end
        
        json_text = jsonencode(cfg_out, 'PrettyPrint', true);                
        fid = fopen(json_file, 'w');
        if fid ~= -1
            fprintf(fid, '%s', json_text);
            fclose(fid);
        end
    catch ME
        fprintf('Warning: Could not save configuration JSON: %s\n', ME.message);
    end
    
    % Clear variables used only for preliminary setup
    clear json_text;
    if exist('customConfig', 'var'), clear customConfig; end
    if exist('fnames', 'var'), clear fnames; end
    if exist('cfg_out', 'var'), clear cfg_out; end
    
    try
        %% 4. Simulation
        hWait = [];
        if usejava('desktop')
            hWait = waitbar(0, 'Running Robustness Study...');
            
            % Resize waitbar: Sligthly wider for 2-line text, and slightly taller
            pos = hWait.Position;
            width_increase = 5; 
            height_increase = 15;
            hWait.Position = [pos(1)-(width_increase/2), pos(2), pos(3)+width_increase, pos(4)+height_increase];
        end
        
        cleanWait = onCleanup(@() delete_valid(hWait));
        
        % Parallel Pool
        pool = gcp('nocreate');
        if isempty(pool)
            fprintf('Starting new parallel pool...\n');
            pool = parpool('SpmdEnabled', false);
        else
            fprintf('Reusing existing parallel pool (%d workers).\n', pool.NumWorkers);
        end
        % num_workers = pool.NumWorkers; % Unused
        cfg_base.num_workers = pool.NumWorkers;
        fprintf('\n');

        % Run Simulation
        results = simulate(scenarios, params, n_sims_per_cond, refs, cfg_base, temp_dir, styles, lang, hWait, out_dir, ts_str, final_out_dir, colors, modes, limits, n_datasets);

        %% 5. Reporting
        fprintf('All simulations completed. Generating global summary...\n');
        if ~isempty(hWait) && isvalid(hWait)
             waitbar(1, hWait, 'Generating Global Report...'); drawnow;
        end
        
        plot('scientific_reports', results, modes, styles, refs, limits, params, final_out_dir, base_name, char(ts_str));
        
        % Final duration logging
        t_duration = toc(t_start);
        hrs = floor(t_duration / 3600); mins = floor(mod(t_duration, 3600) / 60); secs = round(mod(t_duration, 60));
        
        if hrs > 0, time_str = sprintf('%d hours, %d minutes, %d seconds', hrs, mins, secs);
        elseif mins > 0, time_str = sprintf('%d minutes, %d seconds', mins, secs);
        else, time_str = sprintf('%d seconds', secs);
        end
        fprintf('\nTotal Study Duration: %s\n', time_str);
        
        % Optional: Clean up parallel pool to free memory
        pool = gcp('nocreate');
        if ~isempty(pool)
            fprintf('Parallel pool active (%d workers). Run \"delete(gcp)\" to free memory.\n', pool.NumWorkers);
        end
        
        % --- Final Summary of Saved Files ---
        fprintf('\n');
        fprintf('==================================\n');
        fprintf('   Study Completed Successfully\n');
        fprintf('==================================\n');
        fprintf('Summary Report (PDF):  %s\n', pdf_full);
        fprintf('Configuration (JSON):  %s\n', json_file);
        fprintf('Global Data (CSV):     %s\n', fullfile(out_dir, ['Global_Summary_' char(ts_str) '.csv']));
        fprintf('Detailed Results:      %s\n', fullfile(final_out_dir, 'CSV', ['Results_' char(ts_str)]));
        fprintf('Graphics Folder:       %s\n', dir_graphics);
        fprintf('PDF Folder:            %s\n', dirty_pdfs);
        fprintf('Log File:              %s\n', log_filename);

    catch ME
        % Note: diary is automatically closed by cleanupDiary (onCleanup)
        fprintf('\nError occurred during study: %s\n', ME.message);
        rethrow(ME);
    end
end

function [log_folder, source] = get_writable_convergence_log_path()
    % Determine a valid location for log files.
    % Priority: 1. User Documents (Persistent), 2. TempDir (Fallback)
    
    app_folder_name = 'HERA_convergence_Log';
    
    % Try standard Documents folder first
    if ispc
        docs_path = fullfile(getenv('USERPROFILE'), 'Documents');
    else
        docs_path = fullfile(getenv('HOME'), 'Documents');
    end
    
    target_path = fullfile(docs_path, app_folder_name);
    
    % Attempt to create/access Documents folder
    try
        if ~exist(target_path, 'dir')
            mkdir(target_path);
        end
        % Simple write test to confirm permissions
        testFile = fullfile(target_path, 'write_test.tmp');
        fid = fopen(testFile, 'w');
        if fid == -1
            error('No write access');
        end
        fclose(fid);
        delete(testFile);
        
        log_folder = target_path;
        source = 'Documents (Persistent)';
        return;
    catch
        % Fallback to tempdir if Documents fails (e.g. Restricted User)
        log_folder = fullfile(tempdir, app_folder_name);
        if ~exist(log_folder, 'dir')
            mkdir(log_folder);
        end
        source = 'TempDir (Fallback)';
    end
end

function delete_valid(h)
    if ~isempty(h) && isvalid(h)
        delete(h);
    end
end
