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
%   log_path_or_mode - (Optional) Output configuration (Default: "" -> Auto-Log).
%                      - ""            : Auto-Log to 'Documents/HERA_convergence_Log'.
%                      - "interactive" : Opens folder selection dialog.
%                      - [Path String] : Uses the specified custom path.
%
% Outputs:
%   results          - Struct containing the raw simulation data (Errors, Costs, Failures)
%                      for all scenarios and methods.
%
% Author: Lukas von Erdmannsdorff

    arguments
        n_sims_per_cond (1,1) double {mustBePositive, mustBeInteger} = 15
        log_path_or_mode (1,1) string = ""
    end

    clc;
    
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
    [n_datasets, modes, scenarios, params, refs, limits, cfg_base, colors] = config(n_sims_per_cond);
    
    fprintf('\n==========================================================\n');
    fprintf('   Scientific Bootstrap Robustness Study (Sims/Cond=%d)\n', n_sims_per_cond);
    fprintf('==========================================================\n');

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
    
    % Overview Plot (Created BEFORE simulation to free RAM for batch processing)
    plot('parameter_overview', params, scenarios, modes, refs, dir_graphics, char(ts_str), pdf_full, dirty_pdfs);
    
    % Write Configuration CSVs (Created BEFORE simulation)
    save_csv('config', modes, scenarios, params, dir_csv, ts_str);

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
        
        % Note: diary is automatically closed by cleanupDiary (onCleanup)
        % --- Final Summary of Saved Files ---
        fprintf('\n');
        fprintf('==================================\n');
        fprintf('   Study Completed Successfully\n');
        fprintf('==================================\n');
        fprintf('Summary Report (PDF):  %s\n', pdf_full);
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
