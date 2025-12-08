function save_sensitivity(results, shared_info, metric_names)
% SAVE_SENSITIVITY - Displays and saves the Borda count sensitivity analysis.
%
% Syntax:
%   HERA.output.save_sensitivity(results, shared_info, metric_names)
%
% Description:
%   If the Borda analysis was performed (permutation of metric orders), this 
%   function outputs the consensus Borda rank and score, alongside the ranks 
%   for each permutation. Saves to a CSV file.
%
% Inputs:
%   results      - Struct containing borda_results and permutations.
%   shared_info  - Struct with file paths and general info.
%   metric_names - Cell array of metric names.
%
% Outputs:
%   None (prints to console, saves CSV).
%
% Author: Lukas von Erdmannsdorff

    import HERA.output.format_text

    % Unpack necessary variables
    lang = shared_info.lang;
    num_datasets = shared_info.num_datasets;
    dataset_names = shared_info.dataset_names;
    csv_dir = shared_info.csv_dir;
    ts = shared_info.config.timestamp;
    
    %% Borda Results
    if isfield(results, 'borda_results') && ~isempty(results.borda_results)
        fprintf(['\n' lang.output.sensitivity.header '\n\n']);
        num_perms = size(results.all_permutation_ranks, 2);
        
        % Dynamically calculate column widths
        % Column 1: Borda Rank
        width_borda = max(strlength(lang.output.sensitivity.borda_rank), strlength(sprintf('%d', num_datasets))) + 2; % +2 for padding
        
        % Column 2: Borda Score (%)
        width_score = max(strlength(lang.output.sensitivity.borda_score), strlength('100.0%')) + 2; 
        
        % Column 3: Dataset
        width_dataset = max(strlength(lang.csv.headers.dataset), max(cellfun(@strlength, dataset_names))) + 2;
        
        % Initialize lists for header texts and their widths
        header_parts_list = {lang.output.sensitivity.borda_rank, ...
                             lang.output.sensitivity.borda_score, ...
                             lang.csv.headers.dataset};
        col_widths = [width_borda, width_score, width_dataset];
        
        % Columns 4 to end: Permutations
        perm_header_parts = cell(1, num_perms);
        for p = 1:num_perms
            perm_indices = results.selected_permutations(p, :);
            primary_hierarchy_indices = results.selected_permutations(1, :);
            % This local_indices lookup correctly handles 2 or 3 metrics
            [~, local_indices] = ismember(perm_indices, primary_hierarchy_indices);
            perm_names = metric_names(local_indices);
            col_name_full = strjoin(perm_names, ' -> ');
            perm_header_parts{p} = col_name_full;
            % Width is determined by the length of the permutation title
            width_perm = strlength(col_name_full) + 2;
            col_widths(end+1) = width_perm;
        end
        header_parts_list = [header_parts_list, perm_header_parts];
    
        % Prepare header for the CSV file
        header_list_csv = cellfun(@(s) sprintf('"%s"', s), header_parts_list, 'UniformOutput', false);
    
        % Define alignment for each column for console output
        alignments = {'c', 'c', 'l'}; % Borda Rank (c), Borda Score (c), Dataset (l)
        alignments(4:3+num_perms) = {'c'}; % All permutation columns centered
        
        % Output headers to the console using the helper function (all centered)
        header_line = strjoin(arrayfun(@(c) format_text(header_parts_list{c}, col_widths(c), 'c'), 1:numel(header_parts_list), 'UniformOutput', false), ' | ');
        fprintf('%s\n', header_line);
        fprintf('%s\n', repmat('-', 1, strlength(header_line)));
    
        % Prepare filename for the sensitivity report
        [~, fName, fExt] = fileparts(lang.files.sensitivity_details);
        fName = strrep(fName, '%s_', '');
        csv_filename_sensitivity = fullfile(csv_dir, [fName, '_', ts, fExt]);
        
        try
            % Attempt to open CSV file for output
            fid_sens = fopen(csv_filename_sensitivity, 'w');
            if fid_sens == -1
                error(lang.errors.file_open_error, csv_filename_sensitivity); 
            end
            
            % Write header
            fprintf(fid_sens, '%s\n', strjoin(header_list_csv, ';'));
        
            % Fill the tables (Console & CSV)
            [~, sort_idx_c] = sort(results.borda_results.rank); % Sort by consensus rank
            
            for r_idx = 1:num_datasets
                d_idx = sort_idx_c(r_idx);
                
                % Prepare data for the current row (both for console and CSV)
                table_data_row = cell(1, numel(header_parts_list));
                table_data_row{1} = sprintf('%d', r_idx); % Borda Rank
                table_data_row{2} = sprintf('%.1f%%', results.borda_results.score(d_idx)); % Borda Score
                table_data_row{3} = dataset_names{d_idx}; % Dataset
                
                for p = 1:num_perms
                    table_data_row{3 + p} = sprintf('%d', results.all_permutation_ranks(d_idx, p));
                end
                
                % Format and print the console row using the helper function
                row_line_console = strjoin(arrayfun(@(c) format_text(table_data_row{c}, col_widths(c), alignments{c}), ...
                    1:numel(header_parts_list), 'UniformOutput', false), ' | ');
                fprintf('%s\n', row_line_console);
                
                % Prepare and write the row to the CSV file
                csv_row_cells = {table_data_row{1}, ... % Rank (number)
                                 ['"' table_data_row{2} '"'], ... % Score (string)
                                 ['"' table_data_row{3} '"']};   % Dataset (string)
                csv_row_cells = [csv_row_cells, table_data_row(4:end)]; % Permutation ranks (numbers)
                fprintf(fid_sens, '%s\n', strjoin(csv_row_cells, ';'));
            end
            
            % Close the sensitivity file successfully
            fclose(fid_sens);  
            fprintf(['\n' lang.output.files.sensitivity_results_saved '\n'], csv_filename_sensitivity);
    
        catch ME
            % Safety cleanup: Ensure file is closed if an error occurs
            if exist('fid_sens', 'var') && fid_sens ~= -1
                fclose(fid_sens); 
            end
            fprintf([lang.errors.file_save_error '\n'], ME.message);
        end
    end
end
