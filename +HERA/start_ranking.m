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
import HERA.start.MainConfig
import HERA.start.DataSelection
import HERA.start.Statistics
import HERA.start.Summary
import HERA.start.ConfigValidator
import HERA.start.Utils
% Load the language file containing all user-facing strings.
lang = Utils.language_code('en');

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
    % Parse optional command-line arguments using a local helper with arguments block
    args = parse_args(varargin{:});

    %% Batch/Server Mode
    % If a configuration file is provided, load it and start the analysis non-interactively.
    if args.configFile ~= ""
        % Check if the specified configuration file exists.
        if exist(args.configFile, 'file')
            fprintf('Loading configuration from: %s\n', args.configFile);
            try
                % Read the entire content of the JSON file into a string.
                json_text = fileread(args.configFile); 
                % Decode the JSON string into a MATLAB structure.
                loadedData = jsondecode(json_text);            
                % Check if the loaded data contains the required 'userInput' structure.
                if isfield(loadedData, 'userInput')
                    % Load defaults from central library
                    userInput_defaults = HERA.default();
                    % Merge the loaded settings with the defaults, filling in any missing fields.
                    userInput = Utils.fill_defaults(loadedData.userInput, userInput_defaults);
                    
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
                    fprintf([lang.errors.config_file_not_found '\n'], args.configFile);
                end
            catch ME
                % Catch and display any errors during file reading or parsing.
                fprintf([lang.errors.load_error '\n'], ME.message);
            end
        else
            % Error if the specified file path is invalid.
            fprintf([lang.errors.config_file_not_found '\n'], args.configFile);
        end
        return; % Exit the script after non-interactive execution.
    end

    %% Unit Test Mode
    % Check if the 'runtest' flag is set 
    if args.runtest ~= ""
        fprintf('=======================\n');
        fprintf('Starting HERA Unit Test\n');
        fprintf('=======================\n');
        
        try
            % Call the test suite. 
            % If logPath argument is provided, pass it through.
            % If empty, call without args to trigger auto-detection in run_unit_test.
            if args.logPath ~= ""
                run_unit_test(args.logPath);
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
    fprintf('Version: %s\n', HERA.get_version());
    pause(0.5);
    
    % Initialize structures to hold the user's configuration.
    defaults = HERA.default(); % Get all default settings from default function.
    userInput = defaults; % Start with a full set of default values.
    config = defaults; % A temporary structure for settings that might be nested later.
    configLoadedFromFile = false; % Flag to track if the configuration was loaded from a file.
    
    try
        %% 1. Main Configuration (Standard vs Manual vs Config)
        [userInput, config, configLoadedFromFile, main_choice] = MainConfig(defaults, lang);
        
        %% 2. Data Selection (if not loaded from file)
        userInput = DataSelection(userInput, configLoadedFromFile, main_choice, defaults, lang);
        % Check if selection was successful
        if ~configLoadedFromFile && isempty(userInput.folderPath)
            return; 
        end
    
        %% 3. Statistics & Output Configuration
        [userInput, config] = Statistics(userInput, config, configLoadedFromFile, main_choice, defaults, lang);
        if ~configLoadedFromFile && isempty(userInput.output_dir)
            return; 
        end
        
        %% Finalize Configuration, Optionally Save, and Start Analysis
        % Nest the detailed 'config' structure inside the main 'userInput' structure.
        userInput.config = config;
        % If the configuration was not loaded, ask the user if they want to save it.
        if ~configLoadedFromFile
            while true
                save_prompt = sprintf('%s (%s/%s) [%s]: ', lang.prompts.save_config, lang.general.yes_char, lang.general.no_char, lang.general.no_char);
                user_input = input(save_prompt, 's');
                HERA.start.UserInterface.check_exit_command(user_input, lang);
                
                [isValid, error_msg, val] = ConfigValidator.validate_save_choice(user_input, lang);
                
                if isValid
                    save_choice = val;
                    break;
                else
                    fprintf('%s\n', error_msg);
                end
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
        
        %% 4. Configuration Summary
        Summary(userInput, configLoadedFromFile, lang);
    
        %% Final Action
        % Format the final action prompt strings.
        start_text = sprintf('[%s]%s', lang.general.start_char, lang.prompts.action_start(2:end));
        new_text = sprintf('[%s]%s', lang.general.new_char, lang.prompts.action_new(2:end));
        abort_text = sprintf('[%s]%s', lang.general.abort_char, lang.prompts.action_abort(2:end));
        % Prompt the user to start, restart, or abort.
        while true
            final_prompt = sprintf('\n%s: %s, %s, %s [%s]: ', ...
                lang.prompts.final_action, start_text, new_text, abort_text, lang.general.start_char);
            user_input = input(final_prompt, 's');
            HERA.start.UserInterface.check_exit_command(user_input, lang);
            
            [isValid, error_msg, val] = ConfigValidator.validate_final_action(user_input, lang);
            
            if isValid
                final_choice = val;
                break;
            else
                fprintf('%s\n', error_msg);
            end
        end
        
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
    catch ME
        % Graceful exit catch block
        if strcmp(ME.identifier, 'HERA:UserExit')
            fprintf('\n%s\n', ME.message);
            return;
        else
            rethrow(ME);
        end
    end
end % End of "while true" loop
end % End of start_ranking function

% Helper Function to parse command-line arguments
function args = parse_args(options)
    arguments
        options.configFile (1,1) string = ""
        options.runtest (1,1) string = ""
        options.logPath (1,1) string = ""
    end
    args = options;
end
