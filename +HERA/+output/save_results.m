function save_results(results, shared_info, metric_names)
% SAVE_RESULTS - Generates the final results table and saves to CSV.
%
% Syntax:
%   HERA.output.save_results(results, shared_info, metric_names)
%
% Description:
%   Formats the final ranking with confidence intervals and statistics (Mean/SD)
%   for all metrics. Displays this table in the console and writes it to a 
%   CSV file defined in shared_info.
%
% Inputs:
%   results      - Struct containing final ranks and CIs.
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
    
    mean_metrics = shared_info.mean_metrics;
    std_metrics = shared_info.std_metrics;
    num_metrics = numel(metric_names);
    
    point_estimate_ranks = results.final_rank; 
    ci_lower_rank = results.ci_lower_rank;
    ci_upper_rank = results.ci_upper_rank;
    
    %% Final results table and saving to CSV
    % Console output of the final ranking
    % Dynamic header
    metric_name_list_str = strjoin(metric_names, ', ');
    % Build the header string manually to support 1, 2, or 3 metrics
    header_string = regexprep(lang.output.final_summary.header, '%s, %s, and %s', '%s'); % Base string
    header_string = regexprep(header_string, '%s, %s', '%s'); % Handle 2 metrics
    fprintf(['\n' header_string '\n\n'], metric_name_list_str);

    % Prepare all data and headers for the console from the language file
    % Dynamically build header parts
    header_parts = {lang.plots.table.rank, lang.csv.headers.dataset, lang.output.tables.rank_ci};
    for i = 1:num_metrics
        header_parts{end+1} = sprintf(lang.csv.headers.mean_sd_template, metric_names{i});
    end
        
    % Prepare table data by sorting datasets according to their final rank
    [~, output_order_idx] = sort(point_estimate_ranks);
    table_data = cell(num_datasets, numel(header_parts));
    for r = 1:num_datasets
        idx = output_order_idx(r);
        
        % Data for final ranking table
        table_data{r, 1} = sprintf('%d', point_estimate_ranks(idx));
        table_data{r, 2} = dataset_names{idx};
        table_data{r, 3} = sprintf('[%s, %s]', strrep(sprintf('%.2f', ci_lower_rank(idx)), ',', '.'), strrep(sprintf('%.2f', ci_upper_rank(idx)), ',', '.'));
        
        % Loop to add metric data
        for i = 1:num_metrics
            table_data{r, 3+i} = strrep(sprintf('%.3f ± %.3f', mean_metrics(idx, i), std_metrics(idx, i)), ',', '.');
        end
    end

    % Calculate dynamic column widths for console output
    col_widths = cellfun(@strlength, header_parts);
    for r = 1:size(table_data, 1)
        for c = 1:size(table_data, 2)
            col_widths(c) = max(col_widths(c), strlength(table_data{r, c}));
        end
    end
    col_widths = col_widths + 2; % Add 2 spaces for padding

    % Define column alignments for console output (all centered)
    alignments = repmat({'c'}, 1, numel(header_parts));

    % Output header and data rows to console using the helper function
    header_line = strjoin(arrayfun(@(c) format_text(header_parts{c}, col_widths(c), alignments{c}), 1:numel(header_parts), 'UniformOutput', false), '|');
    fprintf('%s\n', header_line);
    fprintf('%s\n', repmat('-', 1, strlength(header_line)));

    for r = 1:size(table_data, 1)
        row_line = strjoin(arrayfun(@(c) format_text(table_data{r, c}, col_widths(c), alignments{c}), 1:size(table_data, 2), 'UniformOutput', false), '|');
        fprintf('%s\n', row_line);
    end
    fprintf('%s\n', repmat('-', 1, strlength(header_line)));

    % Write the final results to a CSV file
    [~, fName, fExt] = fileparts(lang.files.results);
    fName = strrep(fName, '%s_', ''); 
    csv_filename_results = fullfile(csv_dir, [fName, '_', ts, fExt]);

    try
        % Attempt to open the file for writing
        fid_results = fopen(csv_filename_results, 'w');
        if fid_results == -1
            error(lang.errors.file_open_error, csv_filename_results); 
        end

        % Write header to CSV
        fprintf(fid_results, '%s\n', strjoin(header_parts, ';'));
        fclose(fid_results);
        
        % Write data rows to CSV using writetable
        % table_data is already a cell array of strings. 
        % We convert it to a table and append.
        var_names = arrayfun(@(x) sprintf('Var%d', x), 1:numel(header_parts), 'UniformOutput', false);
        T = cell2table(table_data, 'VariableNames', var_names);
        
        writetable(T, csv_filename_results, 'Delimiter', ';', 'WriteMode', 'Append', 'WriteVariableNames', false, 'QuoteStrings', true);
        
        fprintf(['\n' lang.output.files.final_results_saved '\n'], csv_filename_results);

    catch ME
        % Safety cleanup: Ensure file is closed if an error occurs
        if exist('fid_results', 'var') && fid_results ~= -1
            fclose(fid_results); 
        end
        % Report the error to the user
        fprintf([lang.errors.file_save_error '\n'], ME.message);
    end
end
