classdef UserInterface
    % USERINTERFACE - Static class for handling user interactions in HERA.start_ranking.
    %
    % Description:
    %   This class serves as the interface layer between the user and the configuration logic.
    %   It contains the loops and prompts necessary to collect valid input from the command line.
    %   It delegates all validation to the 'HERA.start.ConfigValidator' class.
    %
    % Author: Lukas von Erdmannsdorff
    
    methods (Static)
        
        function choice = get_yes_no_input(prompt, default_val, lang)
            % get_yes_no_input - Prompts the user for a Yes/No answer until valid input is received.
            %
            % Syntax:
            %   choice = HERA.start.UserInterface.get_yes_no_input(prompt, default_val, lang);
            %
            % Inputs:
            %   prompt      - (char) The question to ask the user.
            %   default_val - (logical) Default return value if input is empty.
            %   lang        - (struct) Language structure for localized strings.
            %
            % Outputs:
            %   choice      - (logical) true for 'yes', false for 'no'.
            
            % Determine the default string for display (e.g., 'Y' or 'N').
            if default_val, default_str = lang.general.yes_char; else, default_str = lang.general.no_char; end
            
            % Loop until a valid choice is made.
            while true
                % Display prompt with current default.
                user_input = input(sprintf('%s (%s/%s) [%s]: ', prompt, lang.general.yes_char, lang.general.no_char, default_str), 's');
                HERA.start.UserInterface.check_exit_command(user_input, lang);
                
                % Delegate validation to ConfigValidator.
                [isValid, error_msg, val] = HERA.start.ConfigValidator.validate_yes_no(user_input, default_val, lang);
                
                if isValid
                    choice = val;
                    break;
                else
                    fprintf('%s\n', error_msg);
                end
            end
        end
        
        function val = get_numeric_input(prompt, default_val, is_integer, lang)
            % get_numeric_input - Prompts the user for a numeric value.
            %
            % Syntax:
            %   val = HERA.start.UserInterface.get_numeric_input(prompt, default_val, is_integer, lang);
            %
            % Inputs:
            %   prompt      - (char) The input prompt description.
            %   default_val - (double) Default value if input is empty.
            %   is_integer  - (logical) Flag to enforce integer input.
            %   lang        - (struct) Language structure.
            %
            % Outputs:
            %   val         - (double) The validated numeric value.
            
            while true
                val_str = input(sprintf('%s [%g]: ', prompt, default_val), 's');
                HERA.start.UserInterface.check_exit_command(val_str, lang);
                
                % Delegate validation.
                [isValid, error_msg, value] = HERA.start.ConfigValidator.validate_numeric(val_str, default_val, is_integer, lang);
                
                if isValid
                    val = value;
                    break;
                else
                    fprintf('%s\n', error_msg);
                end
            end
        end
        
        function num_workers = get_worker_input(prompt, default_val, lang)
            % get_worker_input - Prompts for number of parallel workers or 'auto'.
            %
            % Syntax:
            %   num_workers = HERA.start.UserInterface.get_worker_input(prompt, default_val, lang);
            
            while true
                user_input = input(sprintf('%s [%s]: ', prompt, default_val), 's');
                HERA.start.UserInterface.check_exit_command(user_input, lang);
                
                [isValid, error_msg, value] = HERA.start.ConfigValidator.validate_worker_count(user_input, default_val, lang);
                
                if isValid
                    num_workers = value;
                    break;
                else
                    fprintf('%s\n', error_msg);
                end
            end
        end
        
        function [manual_B, auto_config] = configure_bootstrap_step(default_auto_cfg, default_manual_B, name, lang)
            % configure_bootstrap_step - Interactive configuration for a specific bootstrap analysis step.
            %
            % Description:
            %   Allows the user to choose between 'Manual' (fixed samples), 'Robust' (auto w/ smoothing),
            %   or 'Simple' (auto w/o smoothing) bootstrap convergence methods.
            %
            % Inputs:
            %   default_auto_cfg - (struct) Default settings for automatic convergence.
            %   default_manual_B - (int) Default number of samples for manual mode.
            %   name             - (char) Name of the analysis step (e.g., 'CI', 'Thresholds').
            %   lang             - (struct) Language structure.
            %
            % Outputs:
            %   manual_B    - (int/empty) User-selected fixed sample size (empty if auto selected).
            %   auto_config - (struct) Configuration for automatic convergence (defaults/overridden).
            
            manual_B = []; 
            auto_config = default_auto_cfg; 
            
            while true
                % Format the options for the prompt.
                manual_text = sprintf('[%s]%s', lang.general.manual_char, lang.prompts.choice_manual(2:end));
                robust_text = sprintf('[%s]%s', lang.general.robust_char, lang.prompts.choice_robust(2:end));
                simple_text = sprintf('[%s]%s', lang.general.simple_char, lang.prompts.choice_simple(2:end));

                prompt = sprintf('%s %s, %s, %s? [%s]: ', ...
                    lang.prompts.choose_for, manual_text, robust_text, simple_text, lang.general.robust_char);
                user_input = input(prompt, 's');
                HERA.start.UserInterface.check_exit_command(user_input, lang);
                
                % Validate the method choice first.
                [isValid, error_msg, choice] = HERA.start.ConfigValidator.validate_bootstrap_method(user_input, lang);
                
                if isValid
                    if strcmpi(choice, lang.general.manual_char)
                        % If Manual: Ask for the fixed B value.
                        manual_B = HERA.start.UserInterface.get_numeric_input(lang.prompts.manual_b_value, default_manual_B, true, lang);
                        auto_config = struct(); % Clear auto settings.
                        break;
                    elseif strcmpi(choice, lang.general.robust_char)
                        % If Robust: Confirm and optionally adjust advanced parameters.
                        fprintf([lang.prompts.robust_search_selected '\n'], name);
                        if HERA.start.UserInterface.get_yes_no_input(lang.prompts.adjust_convergence, false, lang)
                            auto_config = HERA.start.UserInterface.get_convergence_params(default_auto_cfg, lang);
                        end
                        break;
                    elseif strcmpi(choice, lang.general.simple_char)
                        % If Simple: Disable smoothing/streak and optionally adjust parameters.
                        fprintf([lang.prompts.simple_search_selected '\n'], name);
                        auto_config.smoothing_window = [];
                        auto_config.convergence_streak_needed = [];
                        if HERA.start.UserInterface.get_yes_no_input(lang.prompts.adjust_convergence, false, lang)
                            auto_config = HERA.start.UserInterface.get_convergence_params(auto_config, lang);
                        end
                        break;
                    end
                else
                    fprintf('%s\n', error_msg);
                end
            end
        end
        
        function cfg_custom = get_convergence_params(default_cfg, lang)
            % get_convergence_params - detailed configuration for convergence parameters.
            %
            % Description:
            %   Prompts the user for advanced settings like B_start, B_end, steps, tolerances, etc.
            
            cfg_custom = default_cfg; 
            
            % 1. Start Value
            cfg_custom.B_start = HERA.start.UserInterface.get_numeric_input(['    ' lang.params.b_start], default_cfg.B_start, true, lang);
            
            % 2. End Value (Must be > Start Value)
            while true
                b_end_val = HERA.start.UserInterface.get_numeric_input(['    ' lang.params.b_end], default_cfg.B_end, true, lang);
                [isValid, error_msg, val] = HERA.start.ConfigValidator.validate_convergence_b_end(b_end_val, cfg_custom.B_start, lang);
                if isValid
                    cfg_custom.B_end = val;
                    break;
                else
                    fprintf('%s\n', error_msg);
                end
            end
            
            % 3. Step size and trials
            cfg_custom.B_step = HERA.start.UserInterface.get_numeric_input(['    ' lang.params.b_step], default_cfg.B_step, true, lang);
            cfg_custom.n_trials = HERA.start.UserInterface.get_numeric_input(['    ' lang.params.n_trials], default_cfg.n_trials, true, lang);
            cfg_custom.convergence_tolerance = HERA.start.UserInterface.get_numeric_input(['    ' lang.params.tolerance], default_cfg.convergence_tolerance, false, lang);
            
            % 4. Robust/Smoothing parameters (only if applicable)
            if ~isempty(cfg_custom.smoothing_window)
                cfg_custom.smoothing_window = HERA.start.UserInterface.get_numeric_input(['    ' lang.params.smoothing], default_cfg.smoothing_window, true, lang);
                cfg_custom.convergence_streak_needed = HERA.start.UserInterface.get_numeric_input(['    ' lang.params.streak], default_cfg.convergence_streak_needed, true, lang);
            end
        end
        
        function summary_str = get_bootstrap_string(config, lang)
            % get_bootstrap_string - Formats the bootstrap settings into a readable summary string.
            %
            % Inputs:
            %   config - (struct) The current configuration subsection.
            %   lang   - (struct) Language structure.
            %
            % Outputs:
            %   summary_str - (char) Multi-line string summarizing the chosen methods.
            
            parts = {}; 
            
            % --- Thresholds Section ---
            if isfield(config, 'manual_B_thr') && ~isempty(config.manual_B_thr)
                parts{end+1} = sprintf(['    ' lang.summary.percentile_manual], config.manual_B_thr);
            elseif isfield(config, 'bootstrap_thresholds') && isstruct(config.bootstrap_thresholds)
                if isempty(config.bootstrap_thresholds.smoothing_window)
                    parts{end+1} = ['    ' lang.summary.percentile_simple];
                else
                    parts{end+1} = ['    ' lang.summary.percentile_robust];
                end
            end
            
            % --- Confidence Intervals Section ---
            if isfield(config, 'manual_B_ci') && ~isempty(config.manual_B_ci)
                parts{end+1} = sprintf(['    ' lang.summary.bca_manual], config.manual_B_ci);
            elseif isfield(config, 'bootstrap_ci') && isstruct(config.bootstrap_ci)
                if isempty(config.bootstrap_ci.smoothing_window)
                    parts{end+1} = ['    ' lang.summary.bca_simple];
                else
                    parts{end+1} = ['    ' lang.summary.bca_robust];
                end
            end   
            
            % --- Ranks Section ---
            if isfield(config, 'manual_B_rank') && ~isempty(config.manual_B_rank)
                parts{end+1} = sprintf(['    ' lang.summary.cluster_manual], config.manual_B_rank);
            elseif isfield(config, 'bootstrap_ranks') && isstruct(config.bootstrap_ranks)
                if isempty(config.bootstrap_ranks.smoothing_window)
                    parts{end+1} = ['    ' lang.summary.cluster_simple];
                else
                    parts{end+1} = ['    ' lang.summary.cluster_robust];
                end
            end
            
            % Combine all parts.
            summary_str = strjoin(parts, '\n');
        end

        % --- Exit Command Section ---
        function check_exit_command(user_input, lang)
            % check_exit_command - Checks if the user wants to abort the process.
            % Throws a specific error 'HERA:UserExit' if detected.
            
            if any(strcmpi(user_input, {'exit', 'quit', 'q'}))
                 % Use a custom error ID to catch it cleanly in the main script.
                 error('HERA:UserExit', 'User aborted the configuration.');
            end
        end

    end
end
