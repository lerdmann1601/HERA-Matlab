classdef Utils
    % UTILS - Static utility functions for HERA.start_ranking.
    %
    % Description:
    %   This class contains general utility functions used during the configuration phase
    %   of the ranking process. It handles tasks such as merging user configurations with
    %   default values and loading localization (language) files.
    %
    %   These methods are stateless and support both standard MATLAB execution and
    %   deployed environments (MATLAB Compiler).
    %
    % Author: Lukas von Erdmannsdorff
    
    methods (Static)
        
        function loadedInput = fill_defaults(loadedInput, default)
            % fill_defaults - Recursively merges a loaded configuration with default settings.
            %
            % Description:
            %   This function ensures that a configuration structure (e.g., loaded from a JSON file)
            %   contains all necessary fields required by the HERA application. It takes a loaded
            %   structure and a default structure. Any fields present in 'default' but missing
            %   in 'loadedInput' are added. 
            %
            %   It handles nested structures recursively.
            %   It also intelligently handles empty fields: standard fields are filled with defaults
            %   if empty, but specific "dynamic" fields (like file paths) are left empty to allow
            %   re-prompting or error handling later.
            %
            % Syntax:
            %   loadedInput = HERA.start.Utils.fill_defaults(loadedInput, default);
            %
            % Inputs:
            %   loadedInput - (struct) The configuration structure loaded from a file (or part of it).
            %   default     - (struct) The reference structure containing default values.
            %
            % Outputs:
            %   loadedInput - (struct) The merged structure with defaults filled in.
            
            % Ensure loadedInput is a struct even if the input was empty or invalid.
            if ~isstruct(loadedInput) || isempty(loadedInput), loadedInput = struct(); end 
            
            % Iterate through all fields in the default structure.
            fields = fieldnames(default); 
            for j = 1:length(fields)
                fieldName = fields{j};
                defaultVal = default.(fieldName);
                
                % Case 1: If a field is completely missing in the loaded config, add it from defaults.
                if ~isfield(loadedInput, fieldName)
                    loadedInput.(fieldName) = defaultVal;
                
                % Case 2: If a field is a struct in both, recurse to fill nested fields.
                elseif isstruct(loadedInput.(fieldName)) && isstruct(defaultVal)
                    loadedInput.(fieldName) = HERA.start.Utils.fill_defaults(loadedInput.(fieldName), defaultVal);
                
                % Case 3: If a field exists but is empty, we decide whether to fill it.
                elseif isempty(loadedInput.(fieldName))
                     % Define a list of "dynamic" fields that should NOT be auto-filled if they are empty.
                     % These fields usually require specific user input (paths, file types, specific choices).
                     is_dynamic_path = any(strcmpi(fieldName, {'output_dir', 'folderPath', 'metric_names', 'fileType', 'selected_permutations', 'ranking_mode'})); 
                     
                     % Only fill with default if it is NOT a dynamic path field.
                     if ~is_dynamic_path
                        loadedInput.(fieldName) = defaultVal;
                     end
                end
            end
        end
        
        function lang = language_code(language_code)
            % language_code - Loads the JSON language file for user-facing strings.
            %
            % Description:
            %   Loads the specified language definition file (e.g., 'en.json') and returns
            %   it as a MATLAB structure. This function handles the logic of finding the
            %   file in different environments (standard MATLAB vs. Deployed/Compiled).
            %
            % Syntax:
            %   lang = HERA.start.Utils.language_code('en');
            %
            % Inputs:
            %   language_code - (char) The language code to load (e.g., 'en').
            %
            % Outputs:
            %   lang - (struct) A structure containing all text strings.
            
            possible_paths = {};
            
            % Determine base path relative to THIS file's location (+HERA/+start/Utils.m)
            % We want to reach the +HERA/language folder.
            % 1. Get the directory of this file: .../HERA-Matlab/+HERA/+start
            current_path = fileparts(mfilename('fullpath')); 
            % 2. Get the parent directory: .../HERA-Matlab/+HERA
            hera_path = fileparts(current_path); 
            
            % Define search paths based on execution mode
            if isdeployed
                % In deployed mode (CTF), paths can vary.
                % 1. Root of the CTF (common for AdditionalFiles inclusion)
                possible_paths{end+1} = fullfile(ctfroot, 'language');
                % 2. Inside +HERA package in CTF
                possible_paths{end+1} = fullfile(ctfroot, '+HERA', 'language');
                % 3. Relative to the function location (most robust internal check)
                possible_paths{end+1} = fullfile(hera_path, 'language');
            else
                % Standard development path: Inside +HERA/language relative to this file.
                possible_paths{end+1} = fullfile(hera_path, 'language');
            end
            
            file_path = '';
            found = false;
            
            % Iterate through all possible paths to find the file.
            for i = 1:length(possible_paths)
                candidate = fullfile(possible_paths{i}, [language_code, '.json']);
                if exist(candidate, 'file')
                    file_path = candidate;
                    found = true;
                    break;
                end
            end
            
            % Error handling if file is not found.
            if ~found
                if isdeployed
                    % Provide debug info in deployed mode to help diagnose path issues.
                    fprintf('DEBUG: searched %s\n', strjoin(possible_paths, ', '));
                end
                error('Language file for code "%s" not found.', language_code);
            end
            
            % Read and decode the JSON content.
            json_text = fileread(file_path); 
            lang = jsondecode(json_text);
        end

        function s = clean_struct(s)
            % clean_struct - Recursively sanitizes struct fields to pure doubles.
            %
            % Description:
            %   Helper to recursively sanitize struct fields.
            %   Handles structs and cell arrays to ensure deep sanitization of the configuration.
            %   Crucially, it converts numeric-strings (e.g. "3") to real doubles, preventing
            %   type mismatch errors in parallel computing or logical comparisons.
            %
            % Syntax:
            %   cleaned = HERA.start.Utils.clean_struct(inputStruct);
            %
            % Inputs:
            %   s - (struct) The input structure to clean.
            %
            % Outputs:
            %   s - (struct) The sanitized structure.
            
            fields = fieldnames(s);
            for i = 1:numel(fields)
                val = s.(fields{i});
              
                if isstruct(val)
                    % Recursively clean struct content (handle struct arrays correctly)
                    for j = 1:numel(val)
                        val(j) = HERA.start.Utils.clean_struct(val(j));
                    end
                    s.(fields{i}) = val;
                    
                elseif iscell(val)
                    % Recursively clean cell array content
                    for k = 1:numel(val)
                        c_val = val{k};
                        if isstruct(c_val)
                            val{k} = HERA.start.Utils.clean_struct(c_val);
                        elseif isnumeric(c_val)
                            val{k} = double(c_val);
                        elseif isstring(c_val) || ischar(c_val)
                             % Convert numeric strings (e.g. "3") to doubles
                             num_val = str2double(c_val);
                             if ~isnan(num_val) && ~isinf(num_val) && isreal(num_val)
                                 val{k} = num_val;
                             end
                        end
                    end
                    s.(fields{i}) = val;
                    
                elseif isnumeric(val)
                    s.(fields{i}) = double(val);
                
                elseif isstring(val) || ischar(val)
                     % Convert numeric strings (e.g. "3") to doubles
                     num_val = str2double(val);
                     if ~isnan(num_val) && ~isinf(num_val) && isreal(num_val)
                         s.(fields{i}) = num_val;
                     end
                end
            end
        end
        
    end
end
