function start_ranking(varargin)
% START_RANKING - Interactive script to configure and start the ranking.
%
% Syntax:
%   HERA.start_ranking()                                    1. Interactive Mode
%   HERA.start_ranking('configFile', 'path/to/config.json') 2. Batch/Server Mode
%   HERA.start_ranking('runtest', 'true')                   3. Unit Test Mode (Auto Log)
%   HERA.start_ranking('runtest', 'true', 'logPath', '...') 4. Unit Test Mode (Custom Log)

% Description:
%   This function serves as the primary interface for configuring and starting the ranking process.
%   In interactive mode, it guides the user through all necessary settings. 
%   Alternatively, a configuration file can be provided to start the analysis non-interactively (e.g., on a server).
%   All settings are collected in a structure and passed to the main function 'run_ranking.m'.
%
%   It supports three operation modes:
%
%   1. Interactive: Guides the user through settings via CLI prompts (default).
%   2. Config File: Loads settings from a JSON file for automated runs.
%   3. Unit Tests:  Triggers the HERA test suite.
%
% Workflow:
%   1.  Interactive Mode Initialization: Defines default configuration parameters in a central helper function.
%
%   2.  Mode Selection: Offers a quick start (uses defaults), manual configuration, or loading an existing configuration file.
%
%   3.  Manual Configuration (optional): Guides the user through all detailed settings.
%       - Reproducibility: Option to use a fixed seed for the random number generator (RNG) for repeatable results.
%       - Parallel Processing: Sets the number of CPU cores to use for calculations ('auto' or a specific number).
%       - Bootstrap Method: For three separate analysis (thresholds, confidence intervals, rank stability), the user can choose
%                           Manual: A fixed number of bootstrap samples (B-value).
%                           Automatic: A convergence search (robust or simple) to find a stable B-value automatically.
%
%   4.  Data Selection:
%       - Selects the file type of the metric data (.csv or .xlsx).
%       - Loads all metric files from a user-selected folder.
%       - Selects the number of metrics (1, 2, or 3) for the hierarchy.
%       - Defines the primary hierarchical order of the selected metrics.
%       - Selects the ranking logic (e.g., M1+M2 or M1+M3A) if 2 metrics are chosen.
%
%   5.  Additional Analyses Configuration:
%       - Sensitivity Analysis:
%         If enabled (only for >= 2 metrics), the user can select which alternative metric hierarchies should be analyzed.
%       - Power Analysis:
%         If enabled, the user specifies the number of simulations for the post-hoc power analysis.
%
%   6.  Statistical Parameters:
%       - Adjusts the significance level (alpha values) individually for each metric in the hierarchy.
%       - Sets the confidence level for the confidence interval calculations.
%
%   7.  Graphics & Output:
%       - Selects the graphical theme for all plots ('dark' or 'light').
%       - Report Generation: Option to enable/disable the creation of PDF reports 
%         and complex summary plots ('create_reports'). 
%         Disabling this puts the tool in "Diagnostic Batch Mode" (faster execution, 
%         only JSON/CSV data and essential convergence plots are saved).
%       - Defines the main output folder where all results will be saved.
%       - Offers to save the current settings as a reusable .json configuration file.
%
%   8.  Final Review & Start:
%       - Displays a complete summary of all chosen settings for a final check.
%       - Prompts the user to start the analysis, restart the configuration, or abort the process.
%       - Calls the main function 'run_ranking.m' with the final configuration structure.
%
% Inputs:
%   The function can receive optional input arguments to control it non-interactively.
%   This is useful for automated calls, e.g., from the terminal.
%
%   'configFile' - (char) Optional input argument.
%                  Provide the path to a .json file containing the 'userInput' structure.
%                  This will skip the interactive mode and start the analysis directly.
%                  Example call: start_ranking('configFile', 'C:\configs\my_analysis.json')
%   'logPath'    - (char) Optional input argument. Specifies the output folder for unit test logs.
%
% Outputs:
%   The function has no direct return values. Its purpose is to configure and start the ranking process.
%
% Example for creating a configuration file for non-interactive mode:
%   The script can be run without interactive prompts by providing a .json configuration file.
%   This is ideal for automated or server-based analysis.
%
%   Method 1 for Non-Matlab Users (Creating the JSON file manually):
%   Create a plain text file with a .json extension (e.g., 'my_analysis.json') using any text editor.
%   The file must contain a 'userInput' object with the following minimal required fields.
%   The script will automatically fill in all other missing options with default values.
%
%   Minimal example for a .json file:
%   {
%     "userInput": {
%       "folderPath": "C:\\Data\\My_Study",
%       "fileType": ".csv",
%       "metric_names": [
%         "Metric 1",
%         "Metric 2",
%         "Metric 3"
%       ],
%       "ranking_mode": "M1_M2_M3",
%       "output_dir": "C:\\Results\\My_Study"
%     }
%   }
%
%   Important:
%   - In JSON, Windows paths must use double backslashes (\\). macOS/Linux paths use single forward slashes (/).
%   - 'metric_names' defines the hierarchical order of the analysis.
%   - 'ranking_mode' must be specified ('M1', 'M1_M2', 'M1_M3A', 'M1_M2_M3') and must match the number of metrics.
%
%   You can optionally add more settings to override the defaults. For example:
%   {
%     "userInput": {
%       "folderPath": "C:\\Data\\My_Study",
%       "fileType": ".csv",
%       "metric_names": [ "Metric 1", "Metric 2" ],
%       "ranking_mode": "M1_M2",
%       "output_dir": "C:\\Results\\My_Study",
%       "reproducible": true,
%       "seed": 42,
%       "plot_theme": "dark"
%     }
%   }
%
%   Method 2 for Matlab Users (Creating the JSON file programmatically):
%
%   % 1. Create the main structure
%   userInput = struct();
%
%   % 2. Fill the necessary fields
%   userInput.folderPath   = 'C:\Data\My_Study'; % Use your path
%   userInput.fileType     = '.csv';
%   userInput.metric_names = {'Metric 1', 'Metric 2'}; % Note the order!
%   userInput.ranking_mode = 'M1_M2'; % Must match metric_names
%   userInput.output_dir   = 'C:\Results\My_Study'; % Use your path
%
%   % 3. Optionally add more settings to override defaults
%   userInput.reproducible = true;
%   userInput.seed = 42;
%   userInput.run_power_analysis = true;
%   userInput.plot_theme = 'dark';
%
%   % 4. Save the structure in a .json file
%   data_to_save = struct('userInput', userInput); % Wrap in a struct for compatibility
%   json_text = jsonencode(data_to_save, 'PrettyPrint', true);
%   fid = fopen('C:\Configs\my_analysis.json', 'w');
%   fprintf(fid, '%s', json_text);
%   fclose(fid);
%
% Author: Lukas von Erdmannsdorff

% Import the HERA namespace to find internal functions
import HERA.*
% Load the language file containing all user-facing strings.
lang = language_code('en');

% HERA Dependency Check
% (Only check in MATLAB environment; Runtime has them bundled)
if ~isdeployed
    required_toolboxes = {'Statistics and Machine Learning Toolbox', 'Parallel Computing Toolbox'};
    installed_ver = ver;
    installed_toolboxes = {installed_ver.Name};
    missing_toolboxes = {};

    for i = 1:length(required_toolboxes)
        if ~any(strcmp(installed_toolboxes, required_toolboxes{i}))
            missing_toolboxes{end+1} = required_toolboxes{i};
        end
    end

    if ~isempty(missing_toolboxes)
        error_msg = sprintf('HERA requires the following missing toolboxes:\n');
        for i = 1:length(missing_toolboxes)
            error_msg = sprintf('%s - %s\n', error_msg, missing_toolboxes{i});
        end
        error('%s\nPlease install them via the Add-On Explorer.', error_msg);
    end
end

% Check if the function is called with the 'userInput' struct.
if nargin == 1 && isstruct(varargin{1})
    userInput = varargin{1};
else
    % Set up an argument parser to handle optional command-line arguments.
    % This is primarily used for passing a configuration file path.
    p = inputParser;
    % Defines 'configFile' as an optional string argument.
    addParameter(p, 'configFile', '', @ischar); % Defines 'configFile' as an optional string argument.
    % Define 'runtest' as Trigger for Unit Tests.
    addParameter(p, 'runtest', '', @ischar);
    % Define 'logPath' to override default log location for tests
    addParameter(p, 'logPath', '', @ischar);
    parse(p, varargin{:}); % Parses the input arguments.

    %% Batch/Server Mode
    % If a configuration file is provided, load it and start the analysis non-interactively.
    if ~isempty(p.Results.configFile)
        % Check if the specified configuration file exists.
        if exist(p.Results.configFile, 'file')
            fprintf('Loading configuration from: %s\n', p.Results.configFile);
            try
                % Read the entire content of the JSON file into a string.
                json_text = fileread(p.Results.configFile); % Corrected: fileread
                % Decode the JSON string into a MATLAB structure.
                loadedData = jsondecode(json_text);            
                % Check if the loaded data contains the required 'userInput' structure.
                if isfield(loadedData, 'userInput')
                    % Load defaults from central library
                    userInput_defaults = HERA.default();
                    % Merge the loaded settings with the defaults, filling in any missing fields.
                    userInput = fill_defaults(loadedData.userInput, userInput_defaults);
                    
                    % Verify that essential paths and metric names are not empty.
                    % These cannot be determined non-interactively if missing.
                    if isempty(userInput.folderPath) || isempty(userInput.metric_names) || isempty(userInput.output_dir) || isempty(userInput.ranking_mode)
                        fprintf(['\n' lang.errors.missing_critical_batch '\n']);
                        error(lang.errors.batch_config_incomplete);
                    end
                    try
                        files_check = dir(fullfile(userInput.folderPath, ['*', userInput.fileType]));
                        available_metrics_check = cellfun(@(c) regexprep(c, userInput.fileType, ''), {files_check.name}, 'UniformOutput', false);
                        userInput.available_metrics = available_metrics_check;
                    catch ME_term
                        fprintf([lang.errors.metrics_read_error '\n'], userInput.folderPath, ME_term.message);
                        userInput.available_metrics = {}; 
                    end
                    fprintf([lang.start_ranking.batch_defaults_filled '\n']);
                    % Call the main analysis function with the complete configuration.
                    run_ranking(userInput);
                else
                    % Error if the JSON file is not structured correctly.
                    fprintf([lang.errors.config_file_not_found '\n'], p.Results.configFile);
                end
            catch ME
                % Catch and display any errors during file reading or parsing.
                fprintf([lang.errors.load_error '\n'], ME.message);
            end
        else
            % Error if the specified file path is invalid.
            fprintf([lang.errors.config_file_not_found '\n'], p.Results.configFile);
        end
        return; % Exit the script after non-interactive execution.
    end

    %% Unit Test Mode
    % Check if the 'runtest' flag is set 
    if ~isempty(p.Results.runtest)
        fprintf('=======================\n');
        fprintf('Starting HERA Unit Test\n');
        fprintf('=======================\n');
        
        try
            % Call the test suite. 
            % If logPath argument is provided, pass it through.
            % If empty, call without args to trigger auto-detection in run_unit_test.
            if ~isempty(p.Results.logPath)
                run_unit_test(p.Results.logPath);
            else
                run_unit_test(); 
            end
            
            fprintf('\nAll unit tests completed.\n');
            exit_code = 0; % Success
        catch ME
            % Handle and report any errors during testing
            fprintf('\nCritical Error during Unit Test:\n%s\n', ME.message);
            fprintf('In file: %s (Line %d)\n', ME.stack(1).file, ME.stack(1).line);
            exit_code = 1; % Error
        end
        
        % Determine environment and exit accordingly
        if isdeployed
            % In compiled runtime: Quit application with status code
            quit(exit_code);
        else
            % In MATLAB Editor: Just return to prevent running the main app
            return;
        end
    end
    
    % Safety catch: If we are here, neither struct, nor configFile, nor runtest was valid.
    if nargin > 0
        error('Invalid arguments. Please provide a userInput struct, ''configFile'', or ''runtest''.');
    end
end

% Start of the interactive main loop for user configuration.
% The loop allows the user to restart the configuration process.
while true
    %% Interactive Mode Initialization and Default Configuration
    clc; 
    
    % Display the main title of the script.
    fprintf('=====================\n');
    fprintf('%s\n', lang.start_ranking.main_title);
    fprintf('=====================\n');
    
    % Initialize structures to hold the user's configuration.
    defaults = HERA.default(); % Get all default settings from default function.
    userInput = defaults; % Start with a full set of default values.
    config = defaults; % A temporary structure for settings that might be nested later.
    configLoadedFromFile = false; % Flag to track if the configuration was loaded from a file.
    
    %% Standard vs Manual vs Config Configuration
    % Display a summary of the standard (default) settings to the user.
    fprintf('\n%s\n', lang.start_ranking.default_intro); pause(0.5);
    fprintf([' - ' lang.start_ranking.default_reproducibility '\n'], mat2str(defaults.reproducible), defaults.seed); pause(0.5);
    fprintf([' - ' lang.start_ranking.default_ci_level '\n'], defaults.ci_level); pause(0.5);
    fprintf([' - ' lang.start_ranking.default_alpha '\n']); pause(0.5);
    fprintf([' - ' lang.start_ranking.default_bootstrap '\n']); pause(0.5);
    fprintf([' - ' lang.start_ranking.default_analyses '\n']); pause(0.5);
    fprintf([' - ' lang.start_ranking.default_report '\n']); pause(0.5);
    fprintf([' - ' lang.start_ranking.default_theme '\n'], defaults.plot_theme); pause(0.5);
    
    % Loop to get a valid main action choice from the user (standard, manual, load).
    while true
        fprintf('\n%s\n', lang.start_ranking.prompt_hint); pause(0.5); 
        % Format and display the main prompt using strings from the language file.
        main_prompt_str = sprintf(lang.start_ranking.main_prompt_formatted, ...
            lang.start_ranking.main_prompt, ...
            lang.start_ranking.action_standard, ...
            lang.start_ranking.action_manual, ...
            lang.start_ranking.action_load, ...
            lang.general.standard_char);
        main_choice = input(main_prompt_str, 's'); % Get user input as a string.
        if isempty(main_choice), main_choice = lang.general.standard_char; end % Default to 'standard' if input is empty.
    
        % Process the user's choice.
        switch lower(main_choice)
            % Case 1: Use standard settings.
            case lang.general.standard_char
                % For standard mode, ensure manual bootstrap values are cleared to trigger automatic mode.
                config.manual_B_thr = [];
                config.manual_B_ci = [];
                config.manual_B_rank = [];
                
                % Explicitly set key parameters to their default values.
                userInput.reproducible = defaults.reproducible;
                userInput.seed = defaults.seed;
                userInput.num_workers = defaults.num_workers;
                userInput.run_sensitivity_analysis = defaults.run_sensitivity_analysis;
                userInput.run_power_analysis = defaults.run_power_analysis;
                userInput.power_simulations = defaults.power_simulations;
                userInput.plot_theme = defaults.plot_theme;
                userInput.ranking_mode = defaults.ranking_mode; 
                break; % Exit the choice loop.
            % Case 2: Manual configuration.
            case lang.general.manual_char
                fprintf('\n--- %s ---\n', lang.start_ranking.manual_header);
                fprintf('\n%s\n', lang.start_ranking.prompt_hint);
                pause(0.5);
                
                % Guide the user through detailed settings using helper functions.
                % Reproducibility and parallel processing settings.
                userInput.reproducible = get_yes_no_input(lang.prompts.reproducible, defaults.reproducible, lang);
                if userInput.reproducible
                    userInput.seed = get_numeric_input(lang.prompts.seed, defaults.seed, true, lang);
                else
                    userInput.seed = []; % Seed is not needed if not reproducible.
                end
                userInput.num_workers = get_worker_input(lang.prompts.workers, defaults.num_workers, lang);
                
                % Individual bootstrap settings for each analysis step.
                fprintf('\n%s:\n', lang.prompts.bootstrap_thresholds);
                [config.manual_B_thr, config.bootstrap_thresholds] = configure_bootstrap_step(defaults.bootstrap_thresholds, ...
                    defaults.manual_B_thr, lang.prompts.bootstrap_thresholds, lang);
                
                fprintf('\n%s:\n', lang.prompts.bootstrap_bca);
                [config.manual_B_ci, config.bootstrap_ci] = configure_bootstrap_step(defaults.bootstrap_ci, ...
                    defaults.manual_B_ci, lang.prompts.bootstrap_bca, lang);
                
                fprintf('\n%s:\n', lang.prompts.bootstrap_ranking);
                [config.manual_B_rank, config.bootstrap_ranks] = configure_bootstrap_step(defaults.bootstrap_ranks, ...
                    defaults.manual_B_rank, lang.prompts.bootstrap_ranking, lang);
                break; % Exit the choice loop.
            % Case 3: Load configuration from a file.
            case lang.general.load_char
                % Open a file dialog to select a .json configuration file.
                [file, path] = uigetfile('*.json', lang.prompts.load_config_title);
                % If the user cancels the dialog, restart the choice loop.
                if isequal(file, 0)
                    fprintf('%s\n', lang.errors.selection_cancelled);
                    continue; 
                end                
                fullFilePath = fullfile(path, file);
                fprintf(['\n' lang.prompts.loading_config '\n'], fullFilePath);
                try
                    % Read and decode the selected JSON file.
                    json_text = fileread(fullFilePath); % Corrected: fileread
                    loadedData = jsondecode(json_text);
                    % Check for the required 'userInput' field.
                    if isfield(loadedData, 'userInput')
                        % Load defaults and merge the loaded configuration.
                        userInput_defaults = HERA.default();
                        userInput = fill_defaults(loadedData.userInput, userInput_defaults);
                        % Ensure compatibility for direct calls. If input has no .config, wrap it to match expected structure.
                        if ~isfield(userInput, 'config')
                            userInput.config = userInput;
                        end
                        % Check if critical, environment-specific fields are present.
                        if isempty(userInput.folderPath) || isempty(userInput.metric_names) || isempty(userInput.output_dir) ...
                           || isempty(userInput.ranking_mode) % Check ranking_mode
                            fprintf('%s\n', lang.warnings.incomplete_config);
                            configLoadedFromFile = false; % Treat as incomplete, forcing user to select paths.
                            % Reset dynamic paths to ensure they are re-queried.
                            userInput.folderPath = '';
                            userInput.metric_names = {};
                            userInput.output_dir = '';
                            userInput.ranking_mode = ''; % Reset ranking_mode
                        else
                            fprintf('%s\n', lang.prompts.load_success);
                            configLoadedFromFile = true; % Flag that config is complete.
                        end                    
                        % Restore the nested config structure from the loaded file.
                        config = userInput.config;
                        break; % Exit the choice loop.
                    else
                        % If file is invalid, inform user and restart the choice loop.
                        fprintf('%s\n', lang.errors.invalid_config_file);
                        continue; 
                    end
                catch ME
                    % Handle file reading/parsing errors.
                    fprintf([lang.errors.load_error '\n'], ME.message);
                    continue; 
                end
            otherwise
                % Handle invalid input for the main action choice.
                fprintf('%s\n', sprintf(lang.errors.invalid_input_sml, lang.general.standard_char, lang.general.manual_char, lang.general.load_char));
        end
    end
    
    %% Load data and set order (skipped if a complete config was loaded)
    if ~configLoadedFromFile
        % This block prompts for data paths and metric order if not provided by a config file.
        while true
            fprintf('\n%s\n', lang.prompts.select_data);
            % Prompt the user to select the data file type (.csv or .xlsx).
            prompt = sprintf(lang.prompts.file_type_formatted, ...
                lang.prompts.file_type, lang.general.csv_char, lang.general.excel_char, lang.general.csv_char, lang.general.excel_char, lang.general.csv_char);
            choice = input(prompt, 's');
            if isempty(choice), choice = lang.general.csv_char; 
            end % Default to .csv.   
            % Set file type and mask for file search based on user choice.
            if strcmpi(choice, lang.general.csv_char)
                userInput.fileType = '.csv';
                file_mask = '*.csv';
                break;
            elseif strcmpi(choice, lang.general.excel_char)
                userInput.fileType = '.xlsx';
                file_mask = '*.xlsx';
                break;
            else
                fprintf('%s\n', sprintf(lang.errors.invalid_input_ce, lang.general.csv_char, lang.general.excel_char));
            end
        end
        % Loop until a valid folder containing the specified file type is selected.
        while true
            folder_prompt = sprintf(lang.prompts.select_folder, userInput.fileType);
            folderPath = uigetdir(pwd, folder_prompt); % Open folder selection dialog.
            if isequal(folderPath, 0), fprintf('%s\n', lang.errors.selection_cancelled); return; end % Exit if cancelled.
            
            files = dir(fullfile(folderPath, file_mask)); % Find files matching the mask.
            if ~isempty(files)
                userInput.folderPath = folderPath; % Store the valid path.
                break; % Exit the loop.
            else
                fprintf(lang.errors.no_files_found, userInput.fileType); % Inform user if no files are found.
            end
        end
        % Extract metric names from the filenames.
        available_metrics = cellfun(@(c) regexprep(c, userInput.fileType, ''), {files.name}, 'UniformOutput', false);
        userInput.available_metrics = available_metrics;
        fprintf('\n%s\n', lang.prompts.found_metrics);
        % Display the available metrics to the user.
        for i = 1:numel(available_metrics), fprintf('  %d: %s\n', i, available_metrics{i}); 
            pause(0.5);
        end
        
        % Get the number of metrics for the hierarchy. 
        while true
            num_metrics_str = input(sprintf('%s [3]: ', lang.start_ranking.num_metrics_prompt), 's');
            if isempty(num_metrics_str), num_metrics = 3; break; end
            num_metrics = str2double(num_metrics_str);
            if ~isnan(num_metrics) && isscalar(num_metrics) && any(num_metrics == [1, 2, 3])
                break;
            else
                fprintf('%s\n\n', lang.errors.invalid_metric_count);
            end
        end
        
        % Set default order string based on number of metrics.
        default_order_examples = {lang.prompts.metric_order_example_1, lang.prompts.metric_order_example_2, lang.prompts.metric_order_example_3};
        default_order_str = default_order_examples{num_metrics};
        % Loop to get a valid hierarchical order for the metrics.
        while true
            prompt_text = sprintf(lang.prompts.metric_order_dynamic, num_metrics, default_order_str);
            order_choice_str = input(prompt_text, 's');
            order_choice = str2num(order_choice_str);   % Convert string input to a numeric array.        
            % Validate the input: must be N unique numbers within the available range.
            if isempty(order_choice) || numel(order_choice) ~= num_metrics || max(order_choice) > numel(available_metrics) ...
                || min(order_choice) < 1 || numel(unique(order_choice)) ~= num_metrics
                fprintf('%s\n\n', sprintf(lang.errors.invalid_metric_order_dynamic, num_metrics));
            else
                break; % Exit if valid.
            end
        end       
        % Store the selected metric hierarchy.
        userInput.metric_names = available_metrics(order_choice);
        fprintf([lang.prompts.ranking_order_confirm '\n'], strjoin(userInput.metric_names, ' -> '));
        pause(0.5);
        
        % Select ranking logic based on number of metrics. 
        if num_metrics == 1
            userInput.ranking_mode = 'M1';
            userInput.run_sensitivity_analysis = false; % Sensitivity analysis not possible with 1 metric.
        elseif num_metrics == 2
            fprintf('%s\n', lang.start_ranking.metric_logic_2);
            fprintf('  %s\n', lang.start_ranking.metric_logic_2_opt1);
            fprintf('  %s\n', lang.start_ranking.metric_logic_2_opt2);
            while true
                logic_choice_str = input(sprintf('%s [1]: ', lang.start_ranking.metric_logic_2_prompt), 's');
                if isempty(logic_choice_str), logic_choice_str = '1'; end
                if strcmp(logic_choice_str, '1')
                    userInput.ranking_mode = 'M1_M2';
                    break;
                elseif strcmp(logic_choice_str, '2')
                    userInput.ranking_mode = 'M1_M3A';
                    break;
                else
                    fprintf('%s\n', lang.errors.invalid_metric_logic_2);
                end
            end
        else % num_metrics == 3
             % If standard mode, use default. If manual, it's also M1_M2_M3
            userInput.ranking_mode = 'M1_M2_M3';
        end
        
        %% Configuration of Sensitivity Analysis
        % If in manual mode, ask whether to run sensitivity analysis (only if num_metrics > 1).
        if num_metrics > 1
            if strcmpi(main_choice, lang.general.manual_char)
                userInput.run_sensitivity_analysis = get_yes_no_input(lang.prompts.run_sensitivity, defaults.run_sensitivity_analysis, lang);
            end
            % Note: In 'standard' mode, run_sensitivity_analysis is already true by default.
        else
            userInput.run_sensitivity_analysis = false; % Override default if num_metrics is 1
        end
        
        if userInput.run_sensitivity_analysis % This is only true if num_metrics > 1
            % In standard mode, automatically include all possible permutations.
            if strcmpi(main_choice, lang.general.standard_char) 
                all_perms = perms(order_choice);
                userInput.selected_permutations = all_perms;
            else % In manual mode, let the user select which permutations to test.
                all_perms = perms(order_choice);
                % Exclude the primary hierarchy from the list of alternatives.
                all_perms(ismember(all_perms, order_choice, 'rows'), :) = [];
                fprintf('%s\n', lang.prompts.select_permutations);
                % Display the alternative permutations.
                for p_idx = 1:size(all_perms, 1)
                    fprintf('  [%d] %s\n', p_idx, strjoin(available_metrics(all_perms(p_idx, :)), ' -> '));
                    pause(0.5);
                end          
                % Loop to get a valid selection of permutations.
                while true
                    prompt = sprintf(lang.prompts.permutation_choice, lang.general.all);
                    choice_str = input(prompt, 's');
                    if isempty(choice_str), choice_str = lang.general.all; 
                    end % Default to 'all'.     
                    if strcmpi(choice_str, lang.general.all)
                        selected_indices = 1:size(all_perms, 1); % Select all permutations.
                        break;
                    else
                        selected_indices = str2num(choice_str); % Convert string to numeric indices.
                        % Validate the selected indices.
                        if ~isempty(selected_indices) && all(selected_indices >= 1) && all(selected_indices <= size(all_perms, 1))
                            break;
                        else
                            fprintf('%s\n', lang.errors.invalid_input);
                        end
                    end
                end
                % Store the primary hierarchy plus the selected alternative permutations.
                userInput.selected_permutations = [order_choice; all_perms(selected_indices, :)];
            end
        else
            % If sensitivity analysis is disabled, only the primary hierarchy is stored.
            userInput.selected_permutations = order_choice;
        end
        
        %% Statistical Parameters and Power Analysis
        % This section is skipped in 'standard' mode.
        if ~strcmpi(main_choice, lang.general.standard_char)
            fprintf('%s\n', lang.prompts.stat_params_header);
            % Get the confidence level for CIs.
            ci_val = input(sprintf('%s [%.2f]: ', lang.prompts.ci_level, defaults.ci_level));
            if ~isempty(ci_val), config.ci_level = ci_val; end
        
            % Get the alpha level for each metric in the hierarchy.
            alphas = [];
            for i = 1:num_metrics % Use dynamic num_metrics
                alpha_prompt = sprintf(lang.prompts.alpha_level, i, userInput.metric_names{i});
                alpha_val = input(sprintf('%s [%.2f]: ', alpha_prompt, defaults.alphas(1)));
                if isempty(alpha_val), alpha_val = defaults.alphas(1); end % Use default if empty.
                alphas(i) = alpha_val;
            end
            config.alphas = alphas; 
            
            % Configuration of Power Analysis.
            userInput.run_power_analysis = get_yes_no_input(lang.prompts.run_power, defaults.run_power_analysis, lang);
            if userInput.run_power_analysis
                fprintf('%s\n', lang.prompts.power_description);
                pause(0.5);
                % Get the number of simulations for the power analysis.
                userInput.power_simulations = get_numeric_input(lang.prompts.power_sims, defaults.power_simulations, true, lang);
            else
                userInput.power_simulations = []; % Set to empty if disabled.
            end

            %% Report Generation and Graphics Style
            % For Batch usage detail plots and pdf reports can be skipped.  
            % Yes: Generates PDF reports and summary plots (time-consuming).
            % No: "Diagnostic Mode" - Skips PDFs and heavy plots, but keeps convergence plots (in Graphics/) and all CSV/JSON data.
            userInput.create_reports = get_yes_no_input(lang.prompts.create_reports, defaults.create_reports, lang);
            
            % Logic: Only ask for Theme if we are creating reports, OR if we are doing automatic convergence (which creates plots)
            % If create_reports is FALSE, we check if any Manual B value is EMPTY (meaning Auto mode is on).
            ask_for_theme = false;
            if userInput.create_reports
                ask_for_theme = true;
            elseif isempty(config.manual_B_thr) || isempty(config.manual_B_ci) || isempty(config.manual_B_rank)
                ask_for_theme = true; % Need theme for convergence plots
            end

            userInput.plot_theme = 'light'; % Default theme.
            if ask_for_theme
                fprintf('%s\n', lang.prompts.graphics_theme);
                while true
                    % Prompt for theme choice.
                    theme_prompt = sprintf('%s (%s/%s) [%s]: ', ...
                        lang.prompts.theme_choice, lang.general.dark_char, lang.general.light_char, lang.general.light_char);
                    theme_choice = input(theme_prompt, 's');
                    if isempty(theme_choice), theme_choice = lang.general.light_char; end           
                    % Set theme based on valid input.
                    if strcmpi(theme_choice, lang.general.dark_char)
                        userInput.plot_theme = 'dark';
                        break;
                    elseif strcmpi(theme_choice, lang.general.light_char)
                        userInput.plot_theme = 'light';
                        break;
                    else
                        fprintf('%s\n', sprintf(lang.errors.invalid_input_dl, lang.general.dark_char, lang.general.light_char));
                    end
                end
            end
        end
        
        %% Set Path for Saving
        fprintf('%s\n', lang.prompts.select_output_folder);
        
        % Loop until a valid output directory with write permissions is selected.
        while true
            output_dir = uigetdir(pwd, lang.prompts.output_folder_title); % Open folder selection dialog.
            if isequal(output_dir, 0), fprintf('%s\n', lang.errors.selection_cancelled); return; 
            end % Exit if cancelled.
            % Test for write permissions by creating and deleting a temporary file.
            test_filename = fullfile(output_dir, 'permission_test.tmp');
            try
                fid = fopen(test_filename, 'w');
                if fid == -1, error(lang.errors.no_write_permission); end % Check if file could be opened.
                fclose(fid);
                delete(test_filename); % Clean up the temporary file.
                fprintf([lang.prompts.output_folder_success '\n'], output_dir);
                userInput.output_dir = output_dir; % Store the valid path.
                break; % Exit the loop.
            catch ME
                % Inform user if write access is denied.
                fprintf(['\n' lang.errors.write_permission_denied '\n'], output_dir);
            end
        end
    end % End of if ~configLoadedFromFile
    
    %% Finalize Configuration, Optionally Save, and Start Analysis
    % Nest the detailed 'config' structure inside the main 'userInput' structure.
    userInput.config = config;
    % If the configuration was not loaded, ask the user if they want to save it.
    if ~configLoadedFromFile
        while true
            save_prompt = sprintf('%s (%s/%s) [%s]: ', lang.prompts.save_config, lang.general.yes_char, lang.general.no_char, lang.general.no_char);
            save_choice = input(save_prompt, 's');
            if isempty(save_choice), save_choice = lang.general.no_char; end % Default to 'no'.
            if any(strcmpi(save_choice, {lang.general.yes_char, lang.general.no_char})), break;
            else, fprintf('%s\n', lang.errors.invalid_input); end
        end
        if strcmpi(save_choice, lang.general.yes_char)
            % Open a file dialog to specify the save location for the .json file.
            [file, path] = uiputfile('ranking_config.json', lang.prompts.save_config_title);
            if ischar(file) % Check if a valid file name was provided (not cancelled).
                try
                    % Wrap the userInput struct in another struct to match the loading format.
                    data_to_save = struct('userInput', userInput);
                    % Encode the MATLAB struct into a nicely formatted JSON string.
                    json_text = jsonencode(data_to_save, 'PrettyPrint', true);                
                    % Write the JSON string to the selected file.
                    fid = fopen(fullfile(path, file), 'w');
                    fprintf(fid, '%s', json_text);
                    fclose(fid);
                    fprintf([lang.prompts.save_success '\n'], fullfile(path, file));
                catch ME
                    % Handle any errors during file writing.
                    if exist('fid', 'var') && fid ~= -1, fclose(fid); end 
                    fprintf([lang.errors.save_error '\n'], ME.message);
                end
            else
                % Inform the user if the save operation was cancelled.
                fprintf('%s\n', lang.prompts.save_cancelled);
            end
        end
    end
    
    %% Configuration Summary
    % Display a final summary of all settings before starting the analysis.    
    % Get the logic mode string from language file
    logic_mode_str = '';
    if isfield(userInput, 'ranking_mode') % Check if field exists (for robustness)
        switch userInput.ranking_mode
            case 'M1'
                logic_mode_str = lang.summary.logic_M1;
            case 'M1_M2'
                logic_mode_str = lang.summary.logic_M1_M2;
            case 'M1_M3A'
                logic_mode_str = lang.summary.logic_M1_M3A;
            case 'M1_M2_M3'
                logic_mode_str = lang.summary.logic_M1_M2_M3;
            otherwise
                logic_mode_str = 'Unknown';
        end
    end
    % Determine if the graphics theme is relevant for the summary.
    % It is shown if reports are enabled OR if automatic convergence is used (which generates plots).
    % If everything is manual and reports are off, no graphics are created, so the theme is irrelevant.
    show_theme_summary = userInput.create_reports || ...
                         (isempty(userInput.config.manual_B_thr) || ...
                          isempty(userInput.config.manual_B_ci) || ...
                          isempty(userInput.config.manual_B_rank));            
    if ~configLoadedFromFile
        fprintf('\n%s\n', lang.summary.title);
        pause(0.5);
        fprintf(' -> %s: %s\n', lang.summary.data_folder, userInput.folderPath);
        pause(0.5);
        fprintf(' -> %s:\n%s\n', lang.summary.bootstrap_config, get_bootstrap_string(userInput.config, lang)); 
        pause(0.5);
        fprintf(' -> %s:\n', sprintf(lang.summary.metric_hierarchy, logic_mode_str)); 
        fprintf('       %s\n', strjoin(userInput.metric_names, ' -> '));
        pause(0.5);
        % Display sensitivity analysis details if it's enabled.
        if userInput.run_sensitivity_analysis && size(userInput.selected_permutations, 1) > 1
            fprintf([' -> ' lang.summary.sensitivity_analysis '\n'], size(userInput.selected_permutations, 1));
            % Loop through all selected permutations and display them.
            for p = 1:size(userInput.selected_permutations, 1)
                perm_indices = userInput.selected_permutations(p, :);
                perm_names = available_metrics(perm_indices);
                fprintf('       %s\n', strjoin(perm_names, ' -> '));
                pause(0.5);
            end
        end
        % Display power analysis details if it's enabled.
        if userInput.run_power_analysis
            fprintf([' -> ' lang.summary.power_analysis '\n'], userInput.power_simulations);
            pause(0.5);
        end        
        % Display the chosen graphics theme only if graphics are actually generated.
        if show_theme_summary
            if strcmpi(userInput.plot_theme, 'dark')
                theme_str = 'Dark';
            else
                theme_str = 'Light';
            end
            fprintf(' -> %s: %s\n', lang.summary.graphics_theme, theme_str);
            pause(0.5);
        end
        fprintf(' -> %s: %s\n', lang.summary.output_dir, userInput.output_dir);
        pause(0.5);
        % Display if with report or not.
        if userInput.create_reports
            report_str = lang.summary.report_yes;
        else
            report_str = lang.summary.report_no;
        end
        fprintf(' -> %s: %s\n', lang.summary.create_reports, report_str);
        pause(0.5);        
    else % Display summary for a configuration loaded from a file.
        fprintf(['\n' lang.start_ranking.config_loaded_from_file '\n']);
        pause(0.5);
        fprintf(' -> %s: %s\n', lang.summary.data_folder, userInput.folderPath);
        pause(0.5);
        fprintf(' -> %s:\n%s\n', lang.summary.bootstrap_config, get_bootstrap_string(userInput.config, lang)); 
        pause(0.5);
        fprintf(' -> %s:\n', sprintf(lang.summary.metric_hierarchy, logic_mode_str)); 
        fprintf('       %s\n', strjoin(userInput.metric_names, ' -> '));
        pause(0.5);
        % Display sensitivity analysis details if enabled.
        if userInput.run_sensitivity_analysis && size(userInput.selected_permutations, 1) > 1
            fprintf([' -> ' lang.summary.sensitivity_analysis '\n'], size(userInput.selected_permutations, 1));
            % Re-read metric names from the folder to display permutations correctly.
            try
                files_check = dir(fullfile(userInput.folderPath, ['*', userInput.fileType]));
                available_metrics_check = cellfun(@(c) regexprep(c, userInput.fileType, ''), {files_check.name}, 'UniformOutput', false);
                userInput.available_metrics = available_metrics_check;
                for p = 1:size(userInput.selected_permutations, 1)
                    perm_indices = userInput.selected_permutations(p, :);
                    perm_names = available_metrics_check(perm_indices);
                    fprintf('       %s\n', strjoin(perm_names, ' -> '));
                    pause(0.5);
                end
            catch
                % Handle case where files might have changed since config was saved.
                userInput.available_metrics = {};
                fprintf([lang.start_ranking.files_changed_warning '\n']);
                pause(0.5);
            end
        end
        % Display power analysis details if enabled.
        if userInput.run_power_analysis
            fprintf([' -> ' lang.summary.power_analysis '\n'], userInput.power_simulations);
            pause(0.5);
        end        
        % Display theme only if graphics are actually generated.
        if show_theme_summary
            if strcmpi(userInput.plot_theme, 'dark')
                theme_str = 'Dark';
            else
                theme_str = 'Light';
            end
            fprintf(' -> %s: %s\n', lang.summary.graphics_theme, theme_str);
            pause(0.5);
        end

        fprintf(' -> %s: %s\n', lang.summary.output_dir, userInput.output_dir);
        pause(0.5);
        % Display report.
        if userInput.create_reports
            report_str = lang.summary.report_yes;
        else
            report_str = lang.summary.report_no;
        end
        fprintf(' -> %s: %s\n', lang.summary.create_reports, report_str);
        pause(0.5);       
    end

    %% Final Action
    % Format the final action prompt strings.
    start_text = sprintf('[%s]%s', lang.general.start_char, lang.prompts.action_start(2:end));
    new_text = sprintf('[%s]%s', lang.general.new_char, lang.prompts.action_new(2:end));
    abort_text = sprintf('[%s]%s', lang.general.abort_char, lang.prompts.action_abort(2:end));
    % Prompt the user to start, restart, or abort.
    final_prompt = sprintf('\n%s: %s, %s, %s [%s]: ', ...
        lang.prompts.final_action, start_text, new_text, abort_text, lang.general.start_char);
    final_choice = input(final_prompt, 's');
    if isempty(final_choice), final_choice = lang.general.start_char; end % Default to 'start'.  
    switch lower(final_choice)
        case lang.general.start_char
            fprintf([lang.prompts.starting_analysis '\n\n']);
            try
                % Call the main analysis function.
                run_ranking(userInput);
            catch ME
                % Catch and display any errors that occur during the analysis.
                fprintf('%s\n', lang.errors.analysis_aborted);
                fprintf('%s: %s\n', lang.errors.error_message, ME.message);
                fprintf('%s: %s, %s: %d\n', lang.errors.in_file, ME.stack(1).file, lang.errors.at_line, ME.stack(1).line);
            end
            break; % Exit the main configuration loop.
        case lang.general.new_char
            % Restart the main while loop for a new configuration.
            fprintf('%s\n', lang.prompts.restarting_config);
        case lang.general.abort_char
            % Abort the process and exit the script.
            fprintf('%s\n', lang.prompts.ranking_aborted);
            return; 
    end
end % End of "while true" loop

%% Helper functions for user input
% Function to configure a single bootstrap step (manual or automatic).
function [manual_B, auto_config] = configure_bootstrap_step(default_auto_cfg, default_manual_B, name, lang)
    manual_B = []; % Initialize manual B value.
    auto_config = default_auto_cfg; % Start with default automatic settings.
    while true
        % Format the prompt for choosing the bootstrap method.
        manual_text = sprintf('[%s]%s', lang.general.manual_char, lang.prompts.choice_manual(2:end));
        robust_text = sprintf('[%s]%s', lang.general.robust_char, lang.prompts.choice_robust(2:end));
        simple_text = sprintf('[%s]%s', lang.general.simple_char, lang.prompts.choice_simple(2:end));

        prompt = sprintf('%s %s, %s, %s? [%s]: ', ...
            lang.prompts.choose_for, manual_text, robust_text, simple_text, lang.general.robust_char);
        choice = input(prompt, 's');
        if isempty(choice), choice = lang.general.robust_char; end % Default to robust auto.
        
        if strcmpi(choice, lang.general.manual_char)
            % Get a fixed number of bootstrap samples from the user.
            manual_B = get_numeric_input(lang.prompts.manual_b_value, default_manual_B, true, lang);
            auto_config = struct(); % Clear auto settings.
            break;
        elseif strcmpi(choice, lang.general.robust_char)
            % Select robust automatic search.
            fprintf([lang.prompts.robust_search_selected '\n'], name);
            % Ask if user wants to adjust the advanced convergence parameters.
            if get_yes_no_input(lang.prompts.adjust_convergence, false, lang)
                auto_config = get_convergence_params(default_auto_cfg, lang);
            end
            break;
        elseif strcmpi(choice, lang.general.simple_char)
            % Select simple automatic search (no smoothing).
            fprintf([lang.prompts.simple_search_selected '\n'], name);
            auto_config.smoothing_window = [];
            auto_config.convergence_streak_needed = [];
            % Ask if user wants to adjust the advanced convergence parameters.
            if get_yes_no_input(lang.prompts.adjust_convergence, false, lang)
                auto_config = get_convergence_params(auto_config, lang);
            end
            break;
        else
            fprintf('%s\n', lang.errors.invalid_input);
        end
    end
end

% Function to get a validated Yes/No input from the user.
function choice = get_yes_no_input(prompt, default_val, lang)
    if default_val, default_str = lang.general.yes_char; else, default_str = lang.general.no_char; end
    while true
        user_input = input(sprintf('%s (%s/%s) [%s]: ', prompt, lang.general.yes_char, lang.general.no_char, default_str), 's');
        if isempty(user_input), user_input = default_str; end % Apply default if empty.
        % Check if the input is one of the valid choices.
        if any(strcmpi(user_input, {lang.general.yes_char, lang.general.no_char}))
            choice = strcmpi(user_input, lang.general.yes_char); % Return boolean true for 'yes'.
            break;
        else
            fprintf('%s\n', sprintf(lang.errors.invalid_input_yn, lang.general.yes_char, lang.general.no_char));
        end
    end
end

% Function to get a validated numeric (integer or float) input.
function val = get_numeric_input(prompt, default_val, is_integer, lang)
    while true
        val_str = input(sprintf('%s [%g]: ', prompt, default_val), 's');
        if isempty(val_str), val = default_val; break; end % Apply default if empty.    
        val = str2double(val_str); % Convert string to double.
        % Check for validity (is a number, is scalar, is non-negative).
        is_valid = ~isnan(val) && isscalar(val) && val >= 0;
        if is_integer
            is_valid = is_valid && (val == floor(val)); % Also check if it's a whole number.
        end
        
        if is_valid, break;
        else
            if is_integer, fprintf('%s\n', lang.errors.positive_integer);
            else, fprintf('%s\n', lang.errors.positive_number); end
        end
    end
end

% Function to get the number of parallel workers ('auto' or an integer).
function num_workers = get_worker_input(prompt, default_val, lang)
    while true
        user_input = input(sprintf('%s [%s]: ', prompt, default_val), 's');
        if isempty(user_input), user_input = default_val; end % Apply default if empty.        
        % Check for the 'auto' keyword.
        if strcmpi(user_input, lang.general.auto)
            num_workers = lang.general.auto;
            break;
        end       
        % Check for a valid positive integer.
        num_val = str2double(user_input);
        if ~isnan(num_val) && isscalar(num_val) && num_val > 0 && (num_val == floor(num_val))
            num_workers = num_val;
            break;
        else
            fprintf('%s\n', sprintf(lang.errors.integer_or_auto, lang.general.auto));
        end
    end
end

% Function to get detailed convergence parameters from the user.
function cfg_custom = get_convergence_params(default_cfg, lang)
    cfg_custom = default_cfg; % Start with default values.
    % Prompt for each parameter individually.
    cfg_custom.B_start = get_numeric_input(['    ' lang.params.b_start], default_cfg.B_start, true, lang);
    while true
        b_end = get_numeric_input(['    ' lang.params.b_end], default_cfg.B_end, true, lang);
        if b_end > cfg_custom.B_start % B_end must be greater than B_start.
            cfg_custom.B_end = b_end;
            break;
        else
            fprintf([lang.errors.b_end_error '\n'], cfg_custom.B_start);
        end
    end
    cfg_custom.B_step = get_numeric_input(['    ' lang.params.b_step], default_cfg.B_step, true, lang);
    cfg_custom.n_trials = get_numeric_input(['    ' lang.params.n_trials], default_cfg.n_trials, true, lang);
    cfg_custom.convergence_tolerance = get_numeric_input(['    ' lang.params.tolerance], default_cfg.convergence_tolerance, false, lang);
    % Only ask for robust search parameters if they are applicable.
    if ~isempty(cfg_custom.smoothing_window)
        cfg_custom.smoothing_window = get_numeric_input(['    ' lang.params.smoothing], default_cfg.smoothing_window, true, lang);
        cfg_custom.convergence_streak_needed = get_numeric_input(['    ' lang.params.streak], default_cfg.convergence_streak_needed, true, lang);
    end
end

% Function to generate a formatted string summarizing the bootstrap settings.
function summary_str = get_bootstrap_string(config, lang)
    parts = {}; % Initialize a cell array to hold parts of the string.    
    % Add a line for Percentile-Bootstrap (Thresholds).
    if isfield(config, 'manual_B_thr') && ~isempty(config.manual_B_thr)
        parts{end+1} = sprintf(['    ' lang.summary.percentile_manual], config.manual_B_thr);
    elseif isfield(config, 'bootstrap_thresholds') && isstruct(config.bootstrap_thresholds)
        if isempty(config.bootstrap_thresholds.smoothing_window)
            parts{end+1} = ['    ' lang.summary.percentile_simple];
        else
            parts{end+1} = ['    ' lang.summary.percentile_robust];
        end
    end
    % Add a line for BCa-Bootstrap (CI).
    if isfield(config, 'manual_B_ci') && ~isempty(config.manual_B_ci)
        parts{end+1} = sprintf(['    ' lang.summary.bca_manual], config.manual_B_ci);
    elseif isfield(config, 'bootstrap_ci') && isstruct(config.bootstrap_ci)
        if isempty(config.bootstrap_ci.smoothing_window)
            parts{end+1} = ['    ' lang.summary.bca_simple];
        else
            parts{end+1} = ['    ' lang.summary.bca_robust];
        end
    end   
    % Add a line for Cluster-Bootstrap (Ranks).
    if isfield(config, 'manual_B_rank') && ~isempty(config.manual_B_rank)
        parts{end+1} = sprintf(['    ' lang.summary.cluster_manual], config.manual_B_rank);
    elseif isfield(config, 'bootstrap_ranks') && isstruct(config.bootstrap_ranks)
        if isempty(config.bootstrap_ranks.smoothing_window)
            parts{end+1} = ['    ' lang.summary.cluster_simple];
        else
            parts{end+1} = ['    ' lang.summary.cluster_robust];
        end
    end
    % Join the parts with newlines to create the final summary string.
    summary_str = strjoin(parts, '\n');
end

%% Function to merge a loaded configuration with the default settings.
function loadedInput = fill_defaults(loadedInput, default)
    % Ensures loadedInput is a struct even if the file is faulty.
    if ~isstruct(loadedInput) || isempty(loadedInput), loadedInput = struct(); end 
    
    fields = fieldnames(default); % Get all field names from the defaults struct.
    for j = 1:length(fields)
        fieldName = fields{j};
        defaultVal = default.(fieldName);
        
        % Case 1: If a field is missing in the loaded config, add it from defaults.
        if ~isfield(loadedInput, fieldName)
            loadedInput.(fieldName) = defaultVal;
        
        % Case 2: If a field is a struct in both, recurse to fill nested fields.
        elseif isstruct(loadedInput.(fieldName)) && isstruct(defaultVal)
            loadedInput.(fieldName) = fill_defaults(loadedInput.(fieldName), defaultVal);
        
        % Case 3: If a field exists but is empty, fill it with the default value.
        elseif isempty(loadedInput.(fieldName))
             % Exception: Do not overwrite critical, dynamic path fields if they are intentionally empty.
             is_dynamic_path = any(strcmpi(fieldName, {'output_dir', 'folderPath', 'metric_names', 'fileType', 'selected_permutations', 'ranking_mode'})); 
             
             if ~is_dynamic_path
                loadedInput.(fieldName) = defaultVal;
             end
        end
    end
end

%% Function to load the language JSON file.
function lang = language_code(language_code)
    % Define potential paths for the language file
    possible_paths = {};
    
    if isdeployed
        % 1. Root of the CTF (common for AdditionalFiles)
        possible_paths{end+1} = fullfile(ctfroot, 'language');
        % 2. Inside +HERA package in CTF
        possible_paths{end+1} = fullfile(ctfroot, '+HERA', 'language');
        % 3. Relative to the compiled function location
        possible_paths{end+1} = fullfile(fileparts(mfilename('fullpath')), 'language');
    else
        % Standard development path
        possible_paths{end+1} = fullfile(fileparts(mfilename('fullpath')), 'language');
    end
    
    file_path = '';
    found = false;
    
    for i = 1:length(possible_paths)
        candidate = fullfile(possible_paths{i}, [language_code, '.json']);
        if exist(candidate, 'file')
            file_path = candidate;
            found = true;
            break;
        end
    end
    
    if ~found
        % Debugging output for deployed mode
        if isdeployed
            fprintf('DEBUG: ctfroot = %s\n', ctfroot);
            fprintf('DEBUG: mfilename path = %s\n', fileparts(mfilename('fullpath')));
            fprintf('DEBUG: Searched locations:\n');
            for i = 1:length(possible_paths)
                fprintf('  - %s\n', possible_paths{i});
            end
            % List contents of ctfroot to help debugging
            fprintf('DEBUG: Contents of ctfroot:\n');
            d = dir(ctfroot);
            for k=1:length(d), fprintf('  %s\n', d(k).name); end
        end
        error('Language file for code "%s" not found.', language_code);
    end
    
    % Read and decode the language file.
    json_text = fileread(file_path); 
    lang = jsondecode(json_text);
end
end
