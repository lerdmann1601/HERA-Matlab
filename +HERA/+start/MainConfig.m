function [userInput, config, configLoadedFromFile, main_choice] = MainConfig(defaults, lang)
% MAINCONFIG - Handles the main menu selection (Standard, Manual, Load).
%
% Syntax:
%   [userInput, config, configLoadedFromFile, main_choice] = MainConfig(defaults, lang)
%
% Description:
%   Display the main menu to the user, offering "Standard", "Manual", or "Load Config" modes.
%   Initializes the `userInput` and `config` structures based on the selection.
%
% Inputs:
%   defaults - (struct) Default configuration values.
%   lang     - (struct) Language pack for prompts and messages.
%
% Outputs:
%   userInput            - (struct) Initialized user input structure.
%   config               - (struct) Initialized configuration structure.
%   configLoadedFromFile - (logical) True if a configuration was loaded from simple file.
%   main_choice          - (char) The user's main menu choice (e.g., 's', 'm', 'l').
%
% Author: Lukas von Erdmannsdorff

    import HERA.start.ConfigValidator
    import HERA.start.UserInterface
    import HERA.start.Utils

    userInput = defaults;
    config = defaults;
    configLoadedFromFile = false;

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
        % Process the user's choice.
        user_input = input(main_prompt_str, 's'); 
        UserInterface.check_exit_command(user_input, lang); 
        [isValid, error_msg, main_choice] = ConfigValidator.validate_main_action(user_input, lang);
        
        if ~isValid
            fprintf('%s\n', error_msg);
            continue;
        end
        
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
                userInput.reproducible = UserInterface.get_yes_no_input(lang.prompts.reproducible, defaults.reproducible, lang);
                if userInput.reproducible
                    userInput.seed = UserInterface.get_numeric_input(lang.prompts.seed, defaults.seed, true, lang);
                else
                    userInput.seed = []; % Seed is not needed if not reproducible.
                end
                userInput.num_workers = UserInterface.get_worker_input(lang.prompts.workers, defaults.num_workers, lang);
                
                % Individual bootstrap settings for each analysis step.
                fprintf('\n%s:\n', lang.prompts.bootstrap_thresholds);
                [config.manual_B_thr, config.bootstrap_thresholds] = UserInterface.configure_bootstrap_step(defaults.bootstrap_thresholds, ...
                    defaults.manual_B_thr, lang.prompts.bootstrap_thresholds, lang);
                
                fprintf('\n%s:\n', lang.prompts.bootstrap_bca);
                [config.manual_B_ci, config.bootstrap_ci] = UserInterface.configure_bootstrap_step(defaults.bootstrap_ci, ...
                    defaults.manual_B_ci, lang.prompts.bootstrap_bca, lang);
                
                fprintf('\n%s:\n', lang.prompts.bootstrap_ranking);
                [config.manual_B_rank, config.bootstrap_ranks] = UserInterface.configure_bootstrap_step(defaults.bootstrap_ranks, ...
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
                    json_text = fileread(fullFilePath); 
                    loadedData = jsondecode(json_text);
                    % Check for the required 'userInput' field.
                    if isfield(loadedData, 'userInput')
                        % Load defaults and merge the loaded configuration.
                        userInput_defaults = HERA.default();
                        userInput = Utils.fill_defaults(loadedData.userInput, userInput_defaults);
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
                            % userInput.ranking_mode = ''; % Keep loaded ranking_mode if possible, or reset? logic says reset.
                            % user code reset it.
                            userInput.ranking_mode = '';
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
end
