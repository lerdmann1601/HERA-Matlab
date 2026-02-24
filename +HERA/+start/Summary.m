function Summary(userInput, configLoadedFromFile, lang)
% SUMMARY - Displays a final summary of the configuration before execution.
%
% Syntax:
%   Summary(userInput, configLoadedFromFile, lang)
%
% Description:
%   Prints a comprehensive summary of all settings (data path, metrics, statistics, output)
%   to the console, allowing the user to review their choices.
%
% Inputs:
%   userInput            - (struct) Fully configured user input.
%   configLoadedFromFile - (logical) Flag if config was loaded.
%   lang                 - (struct) Language pack.
%
% Outputs:
%   None. Display only.
%
% Author: Lukas von Erdmannsdorff

    import HERA.start.UserInterface

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
        fprintf(' -> %s:\n%s\n', lang.summary.bootstrap_config, UserInterface.get_bootstrap_string(userInput.config, lang)); 
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
                perm_names = userInput.available_metrics(perm_indices);
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
            switch lower(userInput.plot_theme)
                case 'dark'
                    theme_str = 'Dark';
                case 'colourblind dark'
                    theme_str = 'Colourblind Dark';
                case 'colourblind light'
                    theme_str = 'Colourblind Light';
                otherwise
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
        fprintf(' -> %s:\n%s\n', lang.summary.bootstrap_config, UserInterface.get_bootstrap_string(userInput.config, lang)); 
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
                % userInput.available_metrics = available_metrics_check; % Do not overwrite, just use local
                for p = 1:size(userInput.selected_permutations, 1)
                    perm_indices = userInput.selected_permutations(p, :);
                    perm_names = available_metrics_check(perm_indices);
                    fprintf('       %s\n', strjoin(perm_names, ' -> '));
                    pause(0.5);
                end
            catch
                % Handle case where files might have changed since config was saved.
                % userInput.available_metrics = {};
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
            switch lower(userInput.plot_theme)
                case 'dark'
                    theme_str = 'Dark';
                case 'colourblind dark'
                    theme_str = 'Colourblind Dark';
                case 'colourblind light'
                    theme_str = 'Colourblind Light';
                otherwise
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
end
