function userInput = DataSelection(userInput, configLoadedFromFile, main_choice, defaults, lang)
% DATASELECTION - Handles data file, folder, and metric order selection.
%
% Syntax:
%   userInput = DataSelection(userInput, configLoadedFromFile, main_choice, defaults, lang)
%
% Description:
%   Prompts the user to select the data format (CSV/Excel), choose the source folder,
%   identifies available metrics, and allows the user to define the hierarchical order.
%   Also handles sensitivity analysis configuration (permutations) if applicable.
%
% Inputs:
%   userInput            - (struct) Current user input structure.
%   configLoadedFromFile - (logical) Flag if config was loaded.
%   main_choice          - (char) Main menu choice ('s', 'm', 'l').
%   defaults             - (struct) Default settings.
%   lang                 - (struct) Language pack.
%
% Outputs:
%   userInput            - (struct) Updated user input structure with data paths and metrics.
%
% Author: Lukas von Erdmannsdorff

    import HERA.start.ConfigValidator
    import HERA.start.UserInterface

    %% Load data and set order (skipped if a complete config was loaded)
    if ~configLoadedFromFile
        % This block prompts for data paths and metric order if not provided by a config file.
        while true
            fprintf('\n%s\n', lang.prompts.select_data);
            % Prompt the user to select the data file type (.csv or .xlsx).
            prompt = sprintf(lang.prompts.file_type_formatted, ...
                lang.prompts.file_type, lang.general.csv_char, lang.general.excel_char, lang.general.csv_char, lang.general.excel_char, lang.general.csv_char);
            
            user_input = input(prompt, 's');
            [isValid, error_msg, val] = ConfigValidator.validate_file_type(user_input, lang);
            
            if isValid
                userInput.fileType = val;
                file_mask = ['*' val];
                break;
            else
                fprintf('%s\n', error_msg);
            end
        end
        % Loop until a valid folder containing the specified file type is selected.
        while true
            folder_prompt = sprintf(lang.prompts.select_folder, userInput.fileType);
            folderPath = uigetdir(pwd, folder_prompt); % Open folder selection dialog.
            
            [isValid, error_msg, files] = ConfigValidator.validate_folder_content(folderPath, userInput.fileType, lang);
            
            if isValid
                userInput.folderPath = folderPath; % Store the valid path.
                break; % Exit the loop.
            else
                fprintf('%s\n', error_msg);
                if isequal(folderPath, 0), return; end % Exit if cancelled (handled by caller? check return)
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
            user_input = input(sprintf('%s [3]: ', lang.start_ranking.num_metrics_prompt), 's');
            [isValid, error_msg, val] = ConfigValidator.validate_metric_count(user_input, lang);
            
            if isValid
                num_metrics = val;
                break;
            else
                fprintf('%s\n\n', error_msg);
            end
        end
        
        % Set default order string based on number of metrics.
        default_order_examples = {lang.prompts.metric_order_example_1, lang.prompts.metric_order_example_2, lang.prompts.metric_order_example_3};
        default_order_str = default_order_examples{num_metrics};
        % Loop to get a valid hierarchical order for the metrics.
        while true
            prompt_text = sprintf(lang.prompts.metric_order_dynamic, num_metrics, default_order_str);
            user_input = input(prompt_text, 's');
            
            [isValid, error_msg, val] = ConfigValidator.validate_metric_order(user_input, num_metrics, numel(available_metrics), lang);
            
            if isValid
                order_choice = val;
                break;
            else
                fprintf('%s\n\n', error_msg);
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
                user_input = input(sprintf('%s [1]: ', lang.start_ranking.metric_logic_2_prompt), 's');
                [isValid, error_msg, val] = ConfigValidator.validate_ranking_mode_2_metrics(user_input, lang);
                
                if isValid
                    userInput.ranking_mode = val;
                    break;
                else
                    fprintf('%s\n', error_msg);
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
                userInput.run_sensitivity_analysis = UserInterface.get_yes_no_input(lang.prompts.run_sensitivity, defaults.run_sensitivity_analysis, lang);
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
                    user_input = input(prompt, 's');
                    
                    [isValid, error_msg, val] = ConfigValidator.validate_permutation_choice(user_input, size(all_perms, 1), lang);
                    
                    if isValid
                        if strcmp(val, 'all')
                            selected_indices = 1:size(all_perms, 1);
                        else
                            selected_indices = val;
                        end
                        break;
                    else
                        fprintf('%s\n', error_msg);
                    end
                end
                % Store the primary hierarchy plus the selected alternative permutations.
                userInput.selected_permutations = [order_choice; all_perms(selected_indices, :)];
            end
        else
            % If sensitivity analysis is disabled, only the primary hierarchy is stored.
            userInput.selected_permutations = order_choice;
        end
    end
end
