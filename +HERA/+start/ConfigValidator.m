classdef ConfigValidator
    % CONFIGVALIDATOR - Static class for validating user inputs in HERA.start_ranking.
    %
    % Description:
    %   This class serves as a centralized validation logic provider for the interactive configuration
    %   prompts in 'start_ranking.m'. It strictly separates the logic of verifying inputs from the
    %   user interface code.
    %
    %   Each method is designed to be static and stateless, taking a raw user input (string) and
    %   returning a standardized valid/invalid status, an error message (if applicable), and the
    %   processed value (converted to the correct data type).
    %
    % Syntax:
    %   [isValid, error_msg, value] = HERA.ConfigValidator.method_name(input, lang_struct, [options]);
    %
    % Outputs:
    %   isValid   - (logical) True if the input is valid, False otherwise.
    %   error_msg - (char)    A localized error message from the 'lang' structure if invalid.
    %   value     - (varies)  The processed and converted value (e.g., double, char, boolean) ready for use.
    %
    % Usage:
    %   This class is intended to be used within 'while true' loops in the main script:
    %
    %   while true
    %       user_input = input(prompt, 's');
    %       [isValid, error_msg, value] = HERA.ConfigValidator.validate_something(user_input, lang);
    %       if isValid
    %           config.value = value;
    %           break;
    %       else
    %           fprintf('%s\n', error_msg);
    %       end
    %   end
    %
    % Author: Lukas von Erdmannsdorff
    
    methods (Static)
        
        %% Main Menu / General Actions
        
        function [isValid, error_msg, value] = validate_main_action(user_input, lang)
            % validate_main_action - Validates the main menu choice (Standard, Manual, Load).
            %
            % Inputs:
            %   user_input - (char) User input string (s/m/l or similar).
            %   lang       - (struct) Language structure for error messages.
            
            isValid = false;
            error_msg = '';
            value = '';
            
            % Default behavior: Empty input defaults to 'Standard'
            if isempty(user_input)
                value = lang.general.standard_char;
                isValid = true;
                return;
            end
            
            % Check against allowed characters (Standard, Manual, Load)
            if any(strcmpi(user_input, {lang.general.standard_char, lang.general.manual_char, lang.general.load_char}))
                value = lower(user_input);
                isValid = true;
            else
                % Construct error message for invalid selection
                error_msg = sprintf(lang.errors.invalid_input_sml, lang.general.standard_char, lang.general.manual_char, lang.general.load_char);
            end
        end

        function [isValid, error_msg, value] = validate_final_action(user_input, lang)
            % validate_final_action - Validates the final action choice (Start, New, Abort).
            
            isValid = false;
            error_msg = '';
            value = '';
            
            % Default behavior: Empty input defaults to 'Start'
            if isempty(user_input)
                value = lang.general.start_char;
                isValid = true;
                return;
            end
            
            % Check against allowed characters (Start, New, Abort)
            if any(strcmpi(user_input, {lang.general.start_char, lang.general.new_char, lang.general.abort_char}))
                value = lower(user_input);
                isValid = true;
            else
                error_msg = lang.errors.invalid_input;
            end
        end

        function [isValid, error_msg, value] = validate_save_choice(user_input, lang)
            % validate_save_choice - Validates the Yes/No input for saving the configuration.
            
            isValid = false;
            error_msg = '';
            value = '';
            
            % Default behavior: Empty input defaults to 'No' (safety first)
            if isempty(user_input)
                value = lang.general.no_char;
                isValid = true;
                return;
            end
            
            % Check against Yes/No characters
            if any(strcmpi(user_input, {lang.general.yes_char, lang.general.no_char}))
                value = user_input;
                isValid = true;
            else
                error_msg = lang.errors.invalid_input;
            end
        end

        %% Generic Helpers (Numeric & Boolean)
        
        function [isValid, error_msg, value] = validate_yes_no(user_input, default_val, lang)
            % validate_yes_no - General validator for Yes/No questions.
            %
            % Returns:
            %   value - (logical) true for 'yes', false for 'no'.
            
            isValid = false;
            error_msg = '';
            value = default_val;
            
            if default_val, default_str = lang.general.yes_char; else, default_str = lang.general.no_char; end
            
            % Use default if input is empty
            if isempty(user_input)
                value = default_val;
                isValid = true;
                return;
            end
            
            % Check against valid characters
            if any(strcmpi(user_input, {lang.general.yes_char, lang.general.no_char}))
                value = strcmpi(user_input, lang.general.yes_char);
                isValid = true;
            else
                error_msg = sprintf(lang.errors.invalid_input_yn, lang.general.yes_char, lang.general.no_char);
            end
        end
        
        function [isValid, error_msg, value] = validate_numeric(user_input, default_val, is_integer, lang)
            % validate_numeric - check if input is a valid number (int or float).
            %
            % Inputs:
            %   user_input - (char) Raw string input.
            %   default_val - (double) Value to use if input is empty.
            %   is_integer - (logical) If true, enforces integer constraint.
            
            isValid = false;
            error_msg = '';
            value = default_val;
            
            % Use default
            if isempty(user_input)
                value = default_val;
                isValid = true;
                return;
            end
            
            % Convert to number
            val = str2double(user_input);
            % Basic validity check: Not NaN, Scalar, Non-negative (assumed requirement)
            check_valid = ~isnan(val) && isscalar(val) && val >= 0;
            
            % Strict integer check if required
            if is_integer && check_valid
                check_valid = check_valid && (val == floor(val));
            end
            
            if check_valid
                value = val;
                isValid = true;
            else
                if is_integer
                    error_msg = lang.errors.positive_integer;
                else
                    error_msg = lang.errors.positive_number;
                end
            end
        end
        
        function [isValid, error_msg, value] = validate_worker_count(user_input, default_val, lang)
            % validate_worker_count - Validates parallel worker count or 'auto'.
            
            isValid = false;
            error_msg = '';
            value = default_val;
            
            if isempty(user_input)
                value = default_val;
                isValid = true;
                return;
            end
            
            % Special case: 'auto'
            if strcmpi(user_input, lang.general.auto)
                value = lang.general.auto;
                isValid = true;
                return;
            end
            
            % Otherwise, must be a positive integer
            num_val = str2double(user_input);
            if ~isnan(num_val) && isscalar(num_val) && num_val > 0 && (num_val == floor(num_val))
                value = num_val;
                isValid = true;
            else
                error_msg = sprintf(lang.errors.integer_or_auto, lang.general.auto);
            end
        end
        
        %% Data & Folder Selection
        
        function [isValid, error_msg, value] = validate_file_type(user_input, lang)
            % validate_file_type - Validates user selection for .csv or .xlsx.
            
            isValid = false;
            error_msg = '';
            value = '';
            
            % Default to .csv
            if isempty(user_input)
                value = '.csv';
                isValid = true;
                return;
            end
            
            if strcmpi(user_input, lang.general.csv_char)
                value = '.csv';
                isValid = true;
            elseif strcmpi(user_input, lang.general.excel_char)
                value = '.xlsx';
                isValid = true;
            else
                error_msg = sprintf(lang.errors.invalid_input_ce, lang.general.csv_char, lang.general.excel_char);
            end
        end
        
        function [isValid, error_msg, files] = validate_folder_content(folderPath, fileType, lang)
            % validate_folder_content - Checks if the selected folder contains valid files.
            %
            % Inputs:
            %   folderPath - (char) Path returned by uigetdir.
            %   fileType   - (char) Extension to look for (e.g., '.csv').
            
            isValid = false;
            error_msg = '';
            files = [];
            
            % Check if user cancelled the dialog (uigetdir returns 0)
            if isequal(folderPath, 0)
                error_msg = sprintf('%s\n', lang.errors.selection_cancelled);
                return;
            end
            
            % Check for files matching the mask
            mask = ['*' fileType];
            files = dir(fullfile(folderPath, mask));
            
            if ~isempty(files)
                isValid = true;
            else
                error_msg = sprintf(lang.errors.no_files_found, fileType);
            end
        end

        function [isValid, error_msg, value] = validate_metric_count(user_input, lang)
            % validate_metric_count - Validates the number of metrics (1, 2, or 3).
            
            isValid = false;
            error_msg = '';
            value = 3;
            
            if isempty(user_input)
                value = 3;
                isValid = true;
                return;
            end
            
            num_metrics = str2double(user_input);
            % Must be one of the allowed specific integers
            if ~isnan(num_metrics) && isscalar(num_metrics) && any(num_metrics == [1, 2, 3])
                value = num_metrics;
                isValid = true;
            else
                error_msg = lang.errors.invalid_metric_count;
            end
        end
        
        function [isValid, error_msg, value] = validate_metric_order(user_input, num_metrics, max_index, lang)
            % validate_metric_order - Validates the hierarchical order input.
            %
            % Inputs:
            %   num_metrics - (int) how many metrics are being selected (target length).
            %   max_index   - (int) maximum available metric index (source size).
            
            isValid = false;
            error_msg = '';
            value = [];
            
            val = str2num(user_input); %#ok<ST2NM>
            
            % Validation checklist:
            % 1. Not empty
            % 2. Correct number of elements
            % 3. No index out of bounds
            % 4. No index less than 1
            % 5. No duplicates (all unique)
            if isempty(val) || numel(val) ~= num_metrics || max(val) > max_index ...
                || min(val) < 1 || numel(unique(val)) ~= num_metrics
                error_msg = sprintf(lang.errors.invalid_metric_order_dynamic, num_metrics);
            else
                value = val;
                isValid = true;
            end
        end
        
        function [isValid, error_msg, value] = validate_ranking_mode_2_metrics(user_input, lang)
            % validate_ranking_mode_2_metrics - Selects ranking logic for 2-metric case.
            
            isValid = false;
            error_msg = '';
            value = '';
            
            % Default option 1
            if isempty(user_input)
                value = 'M1_M2';
                isValid = true;
                return;
            end
            
            if strcmp(user_input, '1')
                value = 'M1_M2';
                isValid = true;
            elseif strcmp(user_input, '2')
                value = 'M1_M3A';
                isValid = true;
            else
                error_msg = lang.errors.invalid_metric_logic_2;
            end
        end
        
        %% Sensitivity & Statistics
        
        function [isValid, error_msg, value] = validate_permutation_choice(user_input, num_started_perms, lang)
            % validate_permutation_choice - Validates selection of alternative hierarchies.
            
            isValid = false;
            error_msg = '';
            value = [];
            
            % Default or 'all' selection
            if isempty(user_input) || strcmpi(user_input, lang.general.all)
                value = 'all'; 
                isValid = true;
                return;
            end
            
            % Parse numeric indices
            selected_indices = str2num(user_input); %#ok<ST2NM>
            
            % Check bounds
            if ~isempty(selected_indices) && all(selected_indices >= 1) && all(selected_indices <= num_started_perms)
                value = selected_indices;
                isValid = true;
            else
                error_msg = lang.errors.invalid_input;
            end
        end
        
        function [isValid, error_msg, value] = validate_bootstrap_method(user_input, lang)
            % validate_bootstrap_method - Validates choice: manual, robust, or simple.
            
            isValid = false;
            error_msg = '';
            value = '';
            
            if isempty(user_input)
                value = lang.general.robust_char;
                isValid = true;
                return;
            end
            
            if any(strcmpi(user_input, {lang.general.manual_char, lang.general.robust_char, lang.general.simple_char}))
                value = lower(user_input);
                isValid = true;
            else
                error_msg = lang.errors.invalid_input;
            end
        end
        
        function [isValid, error_msg, value] = validate_convergence_b_end(b_end, b_start, lang)
            % validate_convergence_b_end - Ensures B_end is strictly greater than B_start.
            
            isValid = false;
            error_msg = '';
            value = b_end;
            
            if b_end > b_start
                isValid = true;
            else
                error_msg = sprintf([lang.errors.b_end_error], b_start);
            end
        end
        
        %% Graphics
        
        function [isValid, error_msg, value] = validate_theme_choice(user_input, lang)
            % validate_theme_choice - Validates Dark vs Light theme.
            
            isValid = false;
            error_msg = '';
            value = '';
            
            % Default to light
            if isempty(user_input)
                value = 'light';
                isValid = true;
                return;
            end
            
            clean_input = lower(strtrim(user_input));
            if any(strcmp(clean_input, {lang.general.dark_char, 'dark'}))
                value = 'dark';
                isValid = true;
            elseif any(strcmp(clean_input, {lang.general.light_char, 'light'}))
                value = 'light';
                isValid = true;
            elseif any(strcmp(clean_input, {lang.general.cb_light_char, 'colourblind light', 'colourblind_light'}))
                value = 'colourblind light';
                isValid = true;
            elseif any(strcmp(clean_input, {lang.general.cb_dark_char, 'colourblind dark', 'colourblind_dark'}))
                value = 'colourblind dark';
                isValid = true;
            else
                error_msg = sprintf(lang.errors.invalid_input_theme, lang.general.light_char, lang.general.dark_char, lang.general.cb_light_char, lang.general.cb_dark_char);
            end
        end
        
    end
end
