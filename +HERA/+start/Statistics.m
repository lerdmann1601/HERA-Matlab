function [userInput, config] = Statistics(userInput, config, configLoadedFromFile, main_choice, defaults, lang)
% STATISTICS - Handles statistics, power analysis, reports, and output folder.
%
% Syntax:
%   [userInput, config] = Statistics(userInput, config, configLoadedFromFile, main_choice, defaults, lang)
%
% Description:
%   Configures statistical parameters (alpha, CI), power analysis settings,
%   reporting options, and prompts the user to select the output directory.
%
% Inputs:
%   userInput            - (struct) Current user input structure.
%   config               - (struct) Current detailed configuration.
%   configLoadedFromFile - (logical) Flag if config was loaded.
%   main_choice          - (char) Main menu choice.
%   defaults             - (struct) Default settings.
%   lang                 - (struct) Language pack.
%
% Outputs:
%   userInput            - (struct) Updated user input with stats and output dir.
%   config               - (struct) Updated configuration with stats.
%
% Author: Lukas von Erdmannsdorff

    import HERA.start.ConfigValidator
    import HERA.start.UserInterface

    %% Statistical Parameters and Power Analysis
    % This section is skipped in 'standard' mode OR if a config was loaded.
    if ~strcmpi(main_choice, lang.general.standard_char) && ~configLoadedFromFile
        fprintf('%s\n', lang.prompts.stat_params_header);
        % Get the confidence level for CIs.
        config.ci_level = UserInterface.get_numeric_input(lang.prompts.ci_level, defaults.ci_level, false, lang);
    
        % Get the alpha level for each metric in the hierarchy.
        num_metrics = numel(userInput.metric_names);
        alphas = [];
        for i = 1:num_metrics
            alpha_prompt = sprintf(lang.prompts.alpha_level, i, userInput.metric_names{i});
            alphas(i) = UserInterface.get_numeric_input(alpha_prompt, defaults.alphas(1), false, lang);
        end
        config.alphas = alphas; 
        
        % Configuration of Power Analysis.
        userInput.run_power_analysis = UserInterface.get_yes_no_input(lang.prompts.run_power, defaults.run_power_analysis, lang);
        if userInput.run_power_analysis
            fprintf('%s\n', lang.prompts.power_description);
            pause(0.5);
            % Get the number of simulations for the power analysis.
            userInput.power_simulations = UserInterface.get_numeric_input(lang.prompts.power_sims, defaults.power_simulations, true, lang);
        else
            userInput.power_simulations = []; % Set to empty if disabled.
        end

        %% Report Generation and Graphics Style
        userInput.create_reports = UserInterface.get_yes_no_input(lang.prompts.create_reports, defaults.create_reports, lang);
        
        % Logic: Only ask for Theme if we are creating reports, OR if we are doing automatic convergence
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
                theme_prompt = sprintf('%s (%s/%s/%s/%s) [%s]: ', ...
                    lang.prompts.theme_choice, lang.general.dark_char, lang.general.light_char, lang.general.cb_light_char, lang.general.cb_dark_char, lang.general.light_char);
                user_input = input(theme_prompt, 's');
                HERA.start.UserInterface.check_exit_command(user_input, lang);
                
                [isValid, error_msg, val] = ConfigValidator.validate_theme_choice(user_input, lang);
                
                if isValid
                    userInput.plot_theme = val;
                    break;
                else
                    fprintf('%s\n', error_msg);
                end
            end
        end
    end
    
    %% Set Path for Saving
    if ~configLoadedFromFile
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
    end
end
