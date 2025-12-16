function save_ranking_table(final_bootstrap_ranks, final_rank, dataset_names, selected_B_final, lang, csv_dir, ts)
% SAVE_RANKING_TABLE - Displays and saves the bootstrap ranking distribution.
%
% Syntax:
%   HERA.output.save_ranking_table(final_bootstrap_ranks, final_rank, dataset_names, selected_B_final, lang, csv_dir, ts)
%
% Description:
%   Prints the distribution of bootstrap ranks for each dataset to the console
%   and saves the detailed stats to a CSV file.
%
% Inputs:
%   final_bootstrap_ranks - [num_datasets x B] matrix of rank distributions.
%   final_rank            - Vector of primary ranking for sorting output.
%   dataset_names         - Cell array of dataset names.
%   selected_B_final      - Number of bootstrap samples.
%   lang                  - Language structure.
%   csv_dir               - Directory to save the CSV file.
%   ts                    - Timestamp string.
%
% Author: Lukas von Erdmannsdorff

    import HERA.output.format_text

    num_datasets_b = length(dataset_names);
    
    %% Console Output
    % This section calculates the rank distribution and prints a formatted table to the console, sorted by the final rank.

    fprintf(lang.ranking.distribution_header, selected_B_final);

    % Get the sorting order based on the final rank
    [~, sort_idx] = sort(final_rank);

    % Define headers and prepare for dynamic formatting
    header_parts = {lang.ranking.table_rank, lang.ranking.table_dataset, lang.ranking.table_dist};
    alignments = {'c', 'l', 'l'}; % Center, Left, Left
    table_data = cell(num_datasets_b, 3);

    % Loop over datasets in the order of their final rank
    for i_sorted = 1:num_datasets_b
        i = sort_idx(i_sorted); % Get the original index
        
        % Get corresponding data
        current_final_rank = final_rank(i);
        dataset_name = dataset_names{i};
        rank_data = final_bootstrap_ranks(i, :);
        
        % Calculate the counts and percentages of unique ranks
        [unique_ranks, ~, group_idx] = unique(rank_data);
        counts = accumarray(group_idx, 1);
        percentages = (counts / selected_B_final) * 100;
        
        % Create the compact distribution string (e.g., "1: 70.0%, 2: 30.0%")
        dist_parts = cell(1, numel(unique_ranks));
        for j = 1:numel(unique_ranks)
            dist_parts{j} = sprintf('%d: %.1f%%', unique_ranks(j), percentages(j));
        end
        distribution_string = strjoin(dist_parts, ', ');
        
        % Store data for the table row
        table_data(i_sorted, :) = {sprintf('%d', current_final_rank), dataset_name, distribution_string};
    end

    % Calculate dynamic column widths based on content
    col_widths = cellfun(@strlength, header_parts);
    for r = 1:size(table_data, 1)
        for c = 1:size(table_data, 2)
            col_widths(c) = max(col_widths(c), strlength(table_data{r, c}));
        end
    end
    col_widths = col_widths + 2; % Add 2 spaces for padding

    % Create the header row
    header_line_parts = arrayfun(@(c) format_text(header_parts{c}, col_widths(c), alignments{c}), 1:numel(header_parts), 'UniformOutput', false);
    header_line = ['|' strjoin(header_line_parts, '|')];

    % Create the external continuous separator 
    % Its length must match the total width of the other table rows.
    total_width = strlength(header_line);
    external_separator_line = repmat('-', 1, total_width);

    % Print the formatted table
    fprintf('%s\n', header_line);
    fprintf('%s\n', external_separator_line); % Print continuous separator under header

    % Print the table content
    for r = 1:size(table_data, 1)
        row_line_parts = arrayfun(@(c) format_text(table_data{r, c}, col_widths(c), alignments{c}), 1:numel(header_parts), 'UniformOutput', false);
        % Create the data row 
        row_line = ['|' strjoin(row_line_parts, '|')];
        fprintf('%s\n', row_line);
        
        % Use the continuous external separator *between* data rows
        if r < size(table_data, 1)
            fprintf('%s\n', external_separator_line);
        end
    end
    fprintf('%s\n', external_separator_line); % Print continuous bottom border

    %% CSV Output
    % This section calculates the frequency of each rank for each dataset from bootstrap analysis and saves it to a CSV file.

    % Inform the user about the CSV saving process.
    fprintf(['\n' lang.ranking.saving_csv '\n']);
    % Define the output filename.
    [~, fName, fExt] = fileparts(lang.files.bootstrap_rank_csv);
    fName = strrep(fName, '%s_', '');
    csv_filename = fullfile(csv_dir, [fName, '_', ts, fExt]);

    try
        % Open the file for writing.
        fileID = fopen(csv_filename, 'w');
        
        % Write the header row to the CSV file.
        header = {lang.ranking.table_rank, ...         
                  lang.ranking.table_dataset, ...       
                  lang.ranking.csv_bootstrap_rank, ...  
                  lang.ranking.csv_frequency_percent, ... 
                  lang.ranking.csv_frequency_count};     
        
        % Write the header using comma as delimiter (to match the data fprintf below)
        fprintf(fileID, '%s\n', strjoin(header, ','));
        
        % Get the sorting order based on the final rank
        [~, sort_idx] = sort(final_rank);

        % Loop over each dataset in the order of their final rank
        for i_sorted = 1:num_datasets_b
            i = sort_idx(i_sorted); % Get the original index
            
            % Get the corresponding final rank and dataset name.
            current_final_rank = final_rank(i);
            dataset_name = dataset_names{i};
            
            % Get all bootstrap ranks for the current dataset.
            rank_data = final_bootstrap_ranks(i, :);
            
            % Calculate the counts of unique ranks.
            [unique_ranks, ~, group_idx] = unique(rank_data);
            counts = accumarray(group_idx, 1);
            
            % Convert counts to percentages.
            percentages = (counts / selected_B_final) * 100;
            
            % Write one row in the CSV for each observed rank.
            for j = 1:numel(unique_ranks)
                rank_val = unique_ranks(j);
                count_val = counts(j);
                percent_val = percentages(j);
                
                % Write the data row.
                fprintf(fileID, '%d,"%s",%d,%.2f,%d\n', ...
                    current_final_rank, dataset_name, rank_val, percent_val, count_val);
            end
        end
        
        % Close the file.
        fclose(fileID);
        fprintf([lang.ranking.csv_saved '\n'], csv_filename);
        
    catch ME
        % In case of an error (e.g., file permissions), close the file and report.
        if exist('fileID', 'var') && fileID ~= -1
            fclose(fileID);
        end
        fprintf('Error saving bootstrap rank CSV: %s\n', ME.message);
    end
end
