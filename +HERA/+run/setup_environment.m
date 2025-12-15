function [userInput, setupData] = setup_environment(userInput)
% SETUP_ENVIRONMENT - Initialize the analysis environment.
%
% Syntax:
%   [userInput, setupData] = setup_environment(userInput)
%
% Description:
%   This function handles the initialization steps (1-4) of the ranking analysis.
%   It sets up defaults, creates output directories, configures the random number generator,
%   and initializes the parallel pool.
%
% Inputs:
%   userInput - (struct) Initial configuration structure provided by the user.
%
% Outputs:
%   userInput - (struct) Updated configuration structure with defaults applied.
%   setupData - (struct) Structure containing environment variables:
%       .output_dir   - Path to the main output directory.
%       .graphics_dir - Path to the graphics subdirectory.
%       .pdf_dir      - Path to the PDF subdirectory (if reports enabled).
%       .csv_dir      - Path to the CSV output subdirectory.
%       .log_path     - Full path to the log file.
%       .base_name    - Base name for files (derived from log filename).
%       .lang         - Language structure.
%       .config       - Validated configuration struct.
%       .s            - Random number stream.
%       .timestamp    - Timestamp string used for folder naming.
%
% Author: Lukas von Erdmannsdorff

    %% 1. Configuration Setup
    
    % Import the HERA namespace to find internal functions
    import HERA.*

    % Check if at least one data source is provided (File or Memory).
    if ~isfield(userInput, 'folderPath') && ~isfield(userInput, 'custom_data')
        error('HERA:InvalidInput', 'Input structure must contain either "folderPath" (for file mode) or "custom_data" (for developer mode).');
    end

    % Load Central Defaults 
    defaults = HERA.default(); 

    % Ensure language is set immediately as it is needed for logging/error messages
    if ~isfield(userInput, 'language')
        userInput.language = defaults.language;
    end

    % Load the specified language pack to localize all console and plot outputs.
    lang = language_code(userInput.language); 

    % Validate Metric Names (Critical: Cannot be defaulted)
    if ~isfield(userInput, 'metric_names')
        error(lang.errors.metric_names_missing);
    end
    num_metrics = numel(userInput.metric_names);

    % Only fill fields that are missing
    fields = fieldnames(defaults);
    for i = 1:numel(fields)
        fieldName = fields{i};
        % Only apply default if the user/wrapper didn't provide it yet
        if ~isfield(userInput, fieldName)
            userInput.(fieldName) = defaults.(fieldName);
        end
    end

    % Automatic Target Memory Calculation
    % If target_memory is empty (default), calculate it based on system RAM
    if isfield(userInput, 'system') && isstruct(userInput.system)
        % Check if user manually set it (non-empty)
        if isfield(userInput.system, 'target_memory') && ~isempty(userInput.system.target_memory)
             fprintf([lang.run_ranking.ram_manual '\n'], userInput.system.target_memory);
        else
            % Auto calculation
            [calc_mem, ram_gb, ram_status] = HERA.run.get_target_memory();
            userInput.system.target_memory = calc_mem;
            
            if isnan(ram_gb)
                 % Fallback case
                 fprintf([lang.run_ranking.ram_fallback '\n'], ram_status, calc_mem);
            else
                 % Success case
                 fprintf([lang.run_ranking.ram_auto '\n'], ram_gb, calc_mem);
            end
        end
    else
         % Should not happen if defaults are loaded correctly, but safe fallback
         [calc_mem, ram_gb, ram_status] = HERA.run.get_target_memory();
         userInput.system.target_memory = calc_mem;
         if isnan(ram_gb)
             fprintf([lang.run_ranking.ram_fallback '\n'], ram_status, calc_mem);
         else
             fprintf([lang.run_ranking.ram_auto '\n'], ram_gb, calc_mem);
         end
    end

    % Logic Auto-Correction (Dependencies that depend on User Input)
    % Ensure Alphas match metric count (if generic default was loaded)
    if numel(userInput.alphas) ~= num_metrics
        % Expand scalar default to vector
        userInput.alphas = repmat(userInput.alphas(1), 1, num_metrics);
    end

    % Ensure Permutations exist if Sensitivity is on
    if userInput.run_sensitivity_analysis && isempty(userInput.selected_permutations)
        userInput.selected_permutations = perms(1:num_metrics);
    end

    % Ensure compatibility for direct calls. If input has no .config, wrap it to match expected structure.
    if ~isfield(userInput, 'config')
        userInput.config = userInput;
    end
    
    % Ensure metric_names are synced to config (critical for interactive mode)
    userInput.config.metric_names = userInput.metric_names;
    
    % Ensure ranking_mode is synced to config (critical for calculate_ranking logic)
    if isfield(userInput, 'ranking_mode')
        userInput.config.ranking_mode = userInput.ranking_mode;
    end

    %% 2. Environment Initialization
    % Create a unique output folder name using the current timestamp.
    timestamp_folder = string(datetime('now'), 'yyyyMMdd_HHmmss');
    output_dir_name = "Ranking_" + timestamp_folder;
    output_dir = fullfile(userInput.output_dir, output_dir_name);
    mkdir(output_dir);

    % Create subdirectories for organized output.
    % The 'Graphics' folder is always created to store diagnostic convergence plots.
    graphics_dir = fullfile(output_dir, "Graphics");
    mkdir(graphics_dir);
    % The 'PDF' folder is only created if full reports are requested.
    if userInput.create_reports
        pdf_dir = fullfile(output_dir, "PDF");
        mkdir(pdf_dir);
    else
        pdf_dir = ''; % Leave empty if not needed
    end
    csv_dir = fullfile(output_dir, "Output");
    mkdir(csv_dir);

    % Set up the diary to log all console output to a text file.
    log_filename = "Ranking_" + timestamp_folder + ".txt";
    log_path = fullfile(output_dir, log_filename);
    [~, base_name, ~] = fileparts(log_filename);
    diary(char(log_path));
    
    fprintf('=======================\n');
    fprintf('%s\n', lang.run_ranking.header);
    fprintf('=======================\n');
    fprintf('Version: %s\n', HERA.get_version());
    pause(0.5);

    % Save the complete user configuration to a JSON file for reproducibility.
    config_path = fullfile(output_dir, "configuration.json"); 
    try
        userInputToSave = userInput;
        if isfield(userInputToSave, 'system') && isfield(userInputToSave.system, 'target_memory')
            userInputToSave.system = rmfield(userInputToSave.system, 'target_memory');
        end
        % Also check config.system if it exists separately
        if isfield(userInputToSave, 'config') && isfield(userInputToSave.config, 'system') && isfield(userInputToSave.config.system, 'target_memory')
            userInputToSave.config.system = rmfield(userInputToSave.config.system, 'target_memory');
        end
        
        data_to_save = struct('userInput', userInputToSave);
        json_text = jsonencode(data_to_save, 'PrettyPrint', true);     
        
        fid = fopen(config_path, 'w'); 
        if fid == -1
            error(lang.errors.file_open_error, config_path); 
        end
        
        fprintf(fid, '%s', json_text); 
        fclose(fid); 
    catch ME
        if exist('fid', 'var') && fid ~= -1
            fclose(fid); 
        end
        fprintf([lang.errors.file_save_error '\n'], ME.message);
    end

    % Extract the configuration struct for easier access.
    config = userInput.config;
    fprintf('\n%s\n', lang.run_ranking.config_loaded);

    % Display a summary of the chosen bootstrap configurations in the console.
    fprintf('--------------------------------------------------\n');
    fprintf([lang.run_ranking.cpu_cores '\n'], num2str(userInput.num_workers));
    fprintf([lang.run_ranking.reproducibility '\n'], mat2str(userInput.reproducible), userInput.seed); pause(0.5);

    % Display settings for threshold calculation.
    cfg = config; % local alias
    fprintf(['\n-> ' lang.run_ranking.percentile_bootstrap '\n']); pause(0.5);
    if ~isempty(cfg.manual_B_thr)
        fprintf(['   ' lang.run_ranking.method_manual '\n'], cfg.manual_B_thr); pause(0.5);
    else
        cfg_thr = cfg.bootstrap_thresholds;
        fprintf(['   ' lang.run_ranking.method_auto '\n']); pause(0.5);
        if isempty(cfg_thr.smoothing_window)
            fprintf(['   ' lang.run_ranking.convergence_simple '\n']); pause(0.5);
        else
            fprintf(['   ' lang.run_ranking.convergence_robust '\n']); pause(0.5);
            fprintf(['   ' lang.run_ranking.smoothing_window '\n'], cfg_thr.smoothing_window); pause(0.5);
            fprintf(['   ' lang.run_ranking.stable_runs '\n'], cfg_thr.convergence_streak_needed); pause(0.5);
        end
        fprintf(['   ' lang.run_ranking.b_range '\n'], cfg_thr.B_start, cfg_thr.B_step, cfg_thr.B_end); pause(0.5);
        fprintf(['   ' lang.run_ranking.trials_per_b '\n'], cfg_thr.n_trials); pause(0.5);
        fprintf(['   ' lang.run_ranking.tolerance '\n'], cfg_thr.convergence_tolerance * 100); pause(0.5);
    end

    % Display settings for BCa confidence interval calculation.
    fprintf(['\n-> ' lang.run_ranking.bca_bootstrap '\n']); pause(0.5);
    if ~isempty(cfg.manual_B_ci)
        fprintf(['   ' lang.run_ranking.method_manual '\n'], cfg.manual_B_ci); pause(0.5);
    else
        cfg_ci = cfg.bootstrap_ci;
        fprintf(['   ' lang.run_ranking.method_auto '\n']); pause(0.5);
        if isempty(cfg_ci.smoothing_window)
            fprintf(['   ' lang.run_ranking.convergence_simple '\n']); pause(0.5);
        else
            fprintf(['   ' lang.run_ranking.convergence_robust '\n']); pause(0.5);
            fprintf(['   ' lang.run_ranking.smoothing_window '\n'], cfg_ci.smoothing_window); pause(0.5);
            fprintf(['   ' lang.run_ranking.stable_runs '\n'], cfg_ci.convergence_streak_needed); pause(0.5);
        end
        fprintf(['   ' lang.run_ranking.b_range '\n'], cfg_ci.B_start, cfg_ci.B_step, cfg_ci.B_end); pause(0.5);
        fprintf(['   ' lang.run_ranking.trials_per_b '\n'], cfg_ci.n_trials); pause(0.5);
        fprintf(['   ' lang.run_ranking.tolerance '\n'], cfg_ci.convergence_tolerance * 100); pause(0.5);
    end

    % Display settings for rank stability analysis.
    fprintf(['\n-> ' lang.run_ranking.cluster_bootstrap '\n']); pause(0.5);
    if ~isempty(cfg.manual_B_rank)
        fprintf(['   ' lang.run_ranking.method_manual '\n'], cfg.manual_B_rank); pause(0.5);
    else
        cfg_rank = cfg.bootstrap_ranks;
        fprintf(['   ' lang.run_ranking.method_auto '\n']); pause(0.5);
        if isempty(cfg_rank.smoothing_window)
            fprintf(['   ' lang.run_ranking.convergence_simple '\n']); pause(0.5);
        else
            fprintf(['   ' lang.run_ranking.convergence_robust '\n']); pause(0.5);
            fprintf(['   ' lang.run_ranking.smoothing_window '\n'], cfg_rank.smoothing_window); pause(0.5);
            fprintf(['   ' lang.run_ranking.stable_runs '\n'], cfg_rank.convergence_streak_needed); pause(0.5);
        end
        fprintf(['   ' lang.run_ranking.b_range '\n'], cfg_rank.B_start, cfg_rank.B_step, cfg_rank.B_end); pause(0.5);
        fprintf(['   ' lang.run_ranking.trials_per_b '\n'], cfg_rank.n_trials); pause(0.5);
        fprintf(['   ' lang.run_ranking.tolerance '\n'], cfg_rank.convergence_tolerance * 100); pause(0.5);
    end

    fprintf('--------------------------------------------------\n');
    
    % Update config with timestamp
    config.timestamp = char(timestamp_folder);
    
    %% 3. Reproducibility Setup
    s = []; % Initialize the random stream variable.
    if userInput.reproducible
        % Use a fixed seed for reproducible results.
        s = RandStream('mlfg6331_64', 'Seed', userInput.seed);
        fprintf(['\n' lang.run_ranking.rng_reproducible '\n'], userInput.seed);
    else
        % Use a shuffled seed for non-deterministic results.
        s = RandStream('mlfg6331_64', 'Seed', 'shuffle');
        fprintf(['\n' lang.run_ranking.rng_shuffle '\n']);
    end
    RandStream.setGlobalStream(s); % Set the configured stream as the global default.

    %% 4. Parallel Processing Initialization
    % Determine intended number of workers
    if ischar(userInput.num_workers) && strcmp(userInput.num_workers, 'auto')
        target_workers = max(1, parcluster('local').NumWorkers);
        is_auto = true;
    else
        target_workers = userInput.num_workers;
        is_auto = false;
    end

    % Check if a pool is already running 
    currentPool = gcp('nocreate');

    if ~isempty(currentPool)
        if currentPool.NumWorkers ~= target_workers
            % Mismatch -> Restart
            fprintf([lang.run_ranking.pool_restart_mismatch '\n'], currentPool.NumWorkers, target_workers);
            delete(currentPool);
            currentPool = []; % Ensure it enters the creation block below
        else
            % Match -> Reuse
            fprintf([lang.run_ranking.pool_active_skip '\n'], currentPool.NumWorkers);
        end
    end

    if isempty(currentPool)
        % No pool exists (or was just deleted) -> Start one
        if is_auto
            fprintf([lang.run_ranking.parallel_start_auto '\n'], target_workers);
        else
            fprintf([lang.run_ranking.parallel_start_manual '\n'], target_workers);
        end
        % SpmdEnabled=false reduces overhead since HERA only uses parfor, not SPMD features.
        % This optimization is recommended for pure parfor workloads.
        currentPool = parpool(target_workers, 'SpmdEnabled', false);
    end
    
    % Store actual worker count in config for use by all bootstrap functions.
    config.num_workers = currentPool.NumWorkers;
    
    % Explicitly tell the parallel workers where the +HERA package is located.
    if ~isempty(currentPool)
        try
            % Get the directory containing this file (which is inside +HERA/+run)
            [functionPath, ~, ~] = fileparts(mfilename('fullpath'));
            % We need the Parent of Parent of functionPath to get the root of package
            % Setup file is in .../+HERA/+run/
            % functionPath is .../+HERA/+run
            % parent is .../+HERA
            % grandParent is .../
            packageParentDir = fileparts(fileparts(functionPath)); 
            
            % Ensure it is on the Client path
            addpath(packageParentDir);
            
            % Force all workers to add this path AND WAIT until they are done
            F = parfevalOnAll(currentPool, @addpath, 0, packageParentDir);
            wait(F); 
            
        catch ME
            fprintf([lang.run_ranking.path_warning '\n'], ME.message);
        end
    end
    
    % Pack outputs
    setupData.output_dir = output_dir;
    setupData.graphics_dir = graphics_dir;
    setupData.pdf_dir = pdf_dir;
    setupData.csv_dir = csv_dir;
    setupData.log_path = log_path;
    setupData.base_name = base_name;
    setupData.lang = lang;
    setupData.config = config;
    setupData.s = s;
    setupData.timestamp = timestamp_folder;
    setupData.num_workers = currentPool.NumWorkers; % Store actual worker count for bootstrap functions

end

function lang = language_code(language_code)
    % A simplified helper function to load the language JSON file.
    % This is used in non-interactive contexts like the main run script.
    % Get path relative to this function.
    % Note: since this file is in +HERA/+run, we need to go up one level to +HERA, then to language.
    % Original: was in +HERA/run_ranking.m -> base_path was +HERA.
    % Now: in +HERA/+run/setup_environment.m -> base_path is +HERA/+run.
    base_path = fileparts(mfilename('fullpath'));
    % language folder is in +HERA/language.
    % base_path is .../+HERA/+run
    % fileparts(base_path) is .../+HERA
    hera_path = fileparts(base_path);
    file_path = fullfile(hera_path, 'language', [language_code, '.json']);
    if ~exist(file_path, 'file')
        error('Language file for code "%s" not found at %s.', language_code, file_path);
    end
    json_text = fileread(file_path);
    lang = jsondecode(json_text);
end
