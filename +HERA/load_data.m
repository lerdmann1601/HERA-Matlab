function [is_valid, all_data, dataset_names, num_probanden, num_datasets, pair_idx_all] = load_data(userInput, lang)
% LOAD_DATA - Loads, validates, and checks the quality of the input data.
%
% Syntax:
%   [is_valid, all_data, dataset_names, num_probanden, num_datasets, pair_idx_all] = load_data(userInput, lang)
%
% Description:
%   This function handles the entire process of loading and validating the data.
%   It supports two modes of data ingestion:
%   1. Direct Injection: Using pre-loaded matrices provided in 'userInput.custom_data'.
%   2. File Loading: Reading CSV/Excel files from a specified folder if no custom data is present.
%
%   It checks for consistency (same number of subjects and columns across metrics) and missing values.
%   In case of issues, it interacts with the user.
%
% Workflow:
%   1. In-Memory Check: Checks if 'userInput.custom_data' is provided. If so, validates dimensions and 
%      consistency of the provided matrices directly, skipping file operations.
%   2. File Loading (Fallback): If no custom data is found, loops through each metric file specified 
%      in 'userInput', reads the data, and validates file consistency.
%   3. Validation (Missing Values): Checks every pairwise comparison (e.g., A vs B) for the number 
%      of valid (non-NaN) subject pairs. It warns the user interactively if pairs have low data coverage 
%      (based on 'min_data_completeness' in config) and aborts if any pair has zero valid subjects.
%
% Inputs:
%   userInput     - (struct) Contains settings (folderPath, fileType, metric_names) or direct data (custom_data).
%   lang          - (struct) The language pack containing all strings for displaying messages.
%   dataset_names - (cell, optional) Custom names for datasets when using custom_data.
%
% Outputs:
%   is_valid      - (boolean) A flag indicating whether the data is valid and the analysis should proceed.
%   all_data      - (cell) A cell array where each cell contains a data matrix for one metric.
%   dataset_names - (cell) A cell array of strings containing the names of the datasets.
%   num_probanden - (integer) The number of subjects (rows) found in the data.
%   num_datasets  - (integer) The number of datasets (columns) being compared.
%   pair_idx_all  - (matrix) An Nx2 matrix containing the indices for all pairwise comparisons.
%
% Author: Lukas von Erdmannsdorff

arguments
    userInput (1,1) struct
    lang (1,1) struct
end

    % Initialize return values for the error case
    is_valid = false;
    all_data = {};
    dataset_names = {};
    num_probanden = 0;
    num_datasets = 0;
    pair_idx_all = [];

    %% Check for In-Memory Data 
    % Allows calling the toolbox directly with matrices instead of files.
    if isfield(userInput, 'custom_data') && ~isempty(userInput.custom_data)
        fprintf('\nDetected in-memory data input. Skipping file load.\n');
        
        try
            all_data = userInput.custom_data;
            
            % Validate structure (must be cell array of matrices)
            if ~iscell(all_data) || ~all(cellfun(@isnumeric, all_data))
                error('Input "custom_data" must be a cell array of numeric matrices.');
            end
            
            % Validate metric count matches config
            if numel(all_data) ~= numel(userInput.metric_names)
                error('Mismatch: %d metrics provided in data, but %d names configured.', ...
                      numel(all_data), numel(userInput.metric_names));
            end
            
            % Set Dimensions
            [num_probanden, num_datasets] = size(all_data{1});
            
            % Handle Dataset Names (Optional in struct, generate if missing)
            if isfield(userInput, 'dataset_names') && ~isempty(userInput.dataset_names)
                dataset_names = userInput.dataset_names;
                if numel(dataset_names) ~= num_datasets
                    error('Length of "dataset_names" does not match number of columns.');
                end
            else
                % Generate generic names (D1, D2...)
                dataset_names = arrayfun(@(x) sprintf('D%d', x), 1:num_datasets, 'UniformOutput', false);
                fprintf('Note: No dataset names provided. Generated: %s\n', strjoin(dataset_names, ', '));
            end
            
            % Validate consistency across metrics (same dimensions)
            for m = 2:numel(all_data)
                [r, c] = size(all_data{m});
                if r ~= num_probanden || c ~= num_datasets
                    error('Inconsistent dimensions in metric %d. Expected [%d x %d].', ...
                          m, num_probanden, num_datasets);
                end
            end
            
            % Finalize Setup for In-Memory
            pair_idx_all = nchoosek(1:num_datasets, 2);
            % Note: is_valid is not set to true here yet, it will be set at the very end of the function
            
        catch ME
            fprintf('\nError processing in-memory data: %s\n', ME.message);
            is_valid = false;
            return;
        end
        
    end
    
    %% Data loading and validation (File-based)
    % Only runs if 'all_data' is empty (i.e., no custom data was provided)
    if isempty(all_data)
        fprintf(['\n' lang.load_data.import_header '\n'], userInput.folderPath);
        file_extension = userInput.fileType;
        try
            % Initialize variables to store properties of the first loaded file for consistency checks
            first_dataset_names = {};
            first_num_probanden = 0;
    
            % Loop through each specified metric file
            for i = 1:numel(userInput.metric_names)
                metric_name = userInput.metric_names{i};
                % Construct the full path to the current metric file
                filepath = fullfile(userInput.folderPath, [metric_name, file_extension]);
                fprintf([' -> ' lang.load_data.reading_metric '\n'], i, filepath);
                pause(0.5); % Brief pause for better readability in the console
                
                % Automatically detect import options for robust file reading
                opts = detectImportOptions(filepath);
                % Ensure original column headers are preserved
                opts.VariableNamingRule = 'preserve'; 
                % Read the data from the file into a table
                T = readtable(filepath, opts);
                
                % Extract dataset names (column headers), assuming the first column contains subject IDs
                current_dataset_names = T.Properties.VariableNames(2:end);
                % Get the number of subjects (rows) from the current table
                [current_num_probanden, ~] = size(T{:, 2:end});
                
                % Extract the raw data matrix
                metric_data = T{:, 2:end};
                % Handle cases where data might be loaded as a cell array (e.g., with empty cells)
                if iscell(metric_data)
                     % Replace empty cells with NaN for consistent numerical processing
                     metric_data(cellfun(@(c) isempty(c) || (ischar(c) && isempty(strtrim(c))), metric_data)) = {NaN};
                     % Convert the cell array to a standard numeric matrix
                     metric_data = cell2mat(metric_data);
                end
                
                % Verify that the data matrix contains only numeric or NaN values
                is_all_numeric_or_nan = all(isnumeric(metric_data) | isnan(metric_data), 'all');
                if ~is_all_numeric_or_nan
                    % Throw an error if non-numeric data is found
                    error(lang.load_data.error_non_numeric, [metric_name, file_extension]);
                end
    
                % For the first file, store its properties as the reference standard
                if i == 1
                    first_dataset_names = current_dataset_names;
                    first_num_probanden = current_num_probanden;
                    dataset_names = first_dataset_names;
                    num_probanden = first_num_probanden;
                else
                    % For all subsequent files, validate them against the first file
                    % Check if the number of subjects is consistent across all files
                    if current_num_probanden ~= first_num_probanden
                        error(lang.load_data.error_inconsistent_subjects, ...
                              [metric_name, file_extension], current_num_probanden, first_num_probanden);
                    end
                    % Check if the dataset names (and their order) are identical across all files
                    if ~isequal(current_dataset_names, first_dataset_names)
                        error(lang.load_data.error_inconsistent_columns, [metric_name, file_extension]);
                    end
                end
                % Append the validated data matrix for the current metric to the main cell array
                all_data{end+1} = metric_data;
            end
            % Get the number of datasets from the dimensions of the first loaded matrix
            [~, num_datasets] = size(all_data{1});
            
        catch ME
            % Catch any error that occurs during the loading process
            fprintf(['\n' lang.load_data.error_loading_failed '\n']);
            fprintf([lang.load_data.error_reason '\n'], ME.message);
            return; % Exit function on error; is_valid remains false
        end
        
        % Display a summary of the successfully loaded data
        fprintf(['\n' lang.load_data.load_success '\n']);
        fprintf('-------------------------------------------\n');
        fprintf([lang.load_data.summary_subjects '\n'], num_probanden);
        pause(0.5);
        fprintf([lang.load_data.summary_datasets '\n'], num_datasets);
        pause(0.5);
        fprintf(['\n' lang.load_data.summary_datasets_analyzed '\n']);
        for i = 1:numel(dataset_names), fprintf('  -> %s\n', dataset_names{i}); pause(0.5); end
        fprintf(['\n' lang.load_data.summary_primary_metrics '\n']);
        fprintf('%s\n', strjoin(userInput.metric_names, ' -> '));
        pause(0.5);
    end

    %% Validate data for missing values
    fprintf(['\n' lang.load_data.validation_header '\n']);
    fprintf([lang.load_data.checking_completeness '\n']);
    
    % Define threshold based on user config (removes Magic Number)
    % Default is 0.80 (80% valid pairs required) if not specified.
    if isfield(userInput, 'min_data_completeness')
        min_n_threshold = userInput.min_data_completeness;
    else
        min_n_threshold = 0.80; % Fallback for legacy configs or direct calls
    end
    
    % A cell array to store warning messages for the final summary
    low_n_pairs_warnings = {};
    % A flag to track if a critical warning (a pair with zero valid subjects) occurs
    has_critical_warning = false;
    
    % Generate all unique pairwise combinations of dataset indices
    pair_idx_all = nchoosek(1:num_datasets, 2);
    
    % Iterate through each metric to check for missing values
    for m_idx = 1:numel(userInput.metric_names) 
        fprintf(['  ' lang.load_data.metric_label '\n'], userInput.metric_names{m_idx});
        pause(0.5);
        data_metric = all_data{m_idx};
        
        % Iterate through each pair of datasets
        for p_idx = 1:size(pair_idx_all, 1)
            % Get the column indices for the current pair
            i = pair_idx_all(p_idx, 1);
            j = pair_idx_all(p_idx, 2);
            % Find rows where data is present for BOTH datasets in the pair (this is pairwise deletion)
            valid_rows = ~isnan(data_metric(:, i)) & ~isnan(data_metric(:, j));
            % Count the number of subjects with complete data for this specific pair
            effective_n = sum(valid_rows);
            
            % Critical check: if no subjects have valid data for this pair, the analysis is impossible
            if effective_n == 0
                has_critical_warning = true;
                % Generate and store a critical warning message
                warning_msg = sprintf([' -> ' lang.load_data.warning_critical], dataset_names{i}, dataset_names{j});
                fprintf('%s\n', warning_msg);
                low_n_pairs_warnings{end+1} = warning_msg;
            
            % Warning check: if the number of valid subjects is below the defined threshold
            elseif (effective_n / num_probanden) < min_n_threshold
                % Generate and store a low-N warning message
                warning_msg = sprintf([' -> ' lang.load_data.warning_low_n], ...
                    dataset_names{i}, dataset_names{j}, effective_n, num_probanden, (effective_n / num_probanden)*100);
                fprintf('%s\n', warning_msg);
                low_n_pairs_warnings{end+1} = warning_msg;
            end
        end
    end

    % If any warnings were generated, display a summary and prompt the user for action
    if ~isempty(low_n_pairs_warnings)
        fprintf(['\n' lang.load_data.summary_header '\n']);
        pause(0.5);
        % Print each collected warning message
        for k = 1:numel(low_n_pairs_warnings)
            fprintf('%s\n', low_n_pairs_warnings{k}); 
            pause(0.5);
        end
        
        % If a critical warning occurred, the analysis cannot continue. Abort the process.
        if has_critical_warning
            fprintf(['\n' lang.load_data.error_critical_warning '\n']);
            return; % Exit function on error; is_valid remains false
        end
        
        % Inform the user about the pairwise deletion approach and its implications
        fprintf(['\n' lang.load_data.note_pairwise_deletion '\n'])
        pause(0.5);
        
        % Start a loop to get user confirmation to proceed despite the warnings
        while true
            fprintf(['\n' lang.load_data.note_enter_key '\n']);
            % Prompt the user to confirm if they want to continue
            prompt = sprintf(lang.load_data.prompt_continue, lang.general.yes_char, lang.general.no_char, lang.general.yes_char);
            choice = input(prompt, 's');
            % Set the default choice to 'yes' if the user just presses Enter
            if isempty(choice), choice = lang.general.yes_char; end
            
            % If the user agrees, break the loop and continue the analysis
            if strcmpi(choice, lang.general.yes_char)
                fprintf([lang.load_data.analysis_continuing '\n']);
                break;
            % If the user chooses to abort, exit the function
            elseif strcmpi(choice, lang.general.no_char)
                fprintf([lang.load_data.analysis_aborted '\n']);
                return; % Exit function on abort; is_valid remains false
            % Handle invalid user input
            else
                fprintf([lang.load_data.error_invalid_input '\n'], lang.general.yes_char, lang.general.no_char);
            end
        end
    end
    fprintf('-------------------------------------------\n');

    % If all checks are passed and the user agrees to continue, set the validity flag to true
    is_valid = true;
    fprintf([lang.load_data.validation_success '\n']);
end