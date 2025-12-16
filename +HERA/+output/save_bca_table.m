function save_bca_table(z0_d_all, a_d_all, z0_r_all, a_r_all, metric_names, lang, num_pairs, csv_dir, ts)
% SAVE_BCA_TABLE - Displays and saves the BCa correction factors table.
%
% Syntax:
%   HERA.output.save_bca_table(z0_d_all, a_d_all, z0_r_all, a_r_all, metric_names, lang, num_pairs, csv_dir, ts)
%
% Description:
%   Prints a structured table of BCa correction factors (Bias z0, Skewness a)
%   to the console and saves the data to a CSV file.
%
% Inputs:
%   z0_d_all     - Bias correction factors for Cliff's Delta.
%   a_d_all      - Skewness correction factors for Cliff's Delta.
%   z0_r_all     - Bias correction factors for Relative Difference.
%   a_r_all      - Skewness correction factors for Relative Difference.
%   metric_names - Cell array of metric names.
%   lang         - Language structure for text output.
%   num_pairs    - Number of pairs in the analysis.
%   csv_dir      - Directory to save the CSV file.
%   ts           - Timestamp string for the filename.
%
% Author: Lukas von Erdmannsdorff

    import HERA.output.format_text

    num_metrics = numel(metric_names);

    %% Console Output
    % Prints a formatted table with the summarized correction factors to the console.
    fprintf(lang.bca.correction_factors.header, num_pairs);

    % Define headers and prepare for dynamic formatting.
    header_parts = {lang.bca.correction_factors.factor, lang.csv.headers.median, lang.csv.headers.mean, lang.csv.headers.min, lang.csv.headers.max};
    alignments = {'l', 'c', 'c', 'c', 'c'}; % First column left, others centered.

    for metric_idx = 1:num_metrics
        fprintf('\n%s:\n', metric_names{metric_idx});

        % Collect and format all data for this metric into a cell array.
        z0_d = z0_d_all(:, metric_idx);
        a_d = a_d_all(:, metric_idx);
        z0_r = z0_r_all(:, metric_idx);
        a_r = a_r_all(:, metric_idx);

        table_data = cell(4, numel(header_parts));
        table_data(1, :) = {lang.bca.correction_factors.bias_delta, sprintf('%+.3f', median(z0_d)), sprintf('%+.3f', ...
            mean(z0_d)), sprintf('%+.3f', min(z0_d)), sprintf('%+.3f', max(z0_d))};
        table_data(2, :) = {lang.bca.correction_factors.skew_delta, sprintf('%+.3f', median(a_d)), sprintf('%+.3f', ...
            mean(a_d)), sprintf('%+.3f', min(a_d)), sprintf('%+.3f', max(a_d))};
        table_data(3, :) = {lang.bca.correction_factors.bias_rel, sprintf('%+.3f', median(z0_r)), sprintf('%+.3f', ...
            mean(z0_r)), sprintf('%+.3f', min(z0_r)), sprintf('%+.3f', max(z0_r))};
        table_data(4, :) = {lang.bca.correction_factors.skew_rel, sprintf('%+.3f', median(a_r)), sprintf('%+.3f', ...
            mean(a_r)), sprintf('%+.3f', min(a_r)), sprintf('%+.3f', max(a_r))};

        % Calculate dynamic column widths based on content.
        col_widths = cellfun(@strlength, header_parts);
        for r = 1:size(table_data, 1)
            for c = 1:size(table_data, 2)
                col_widths(c) = max(col_widths(c), strlength(table_data{r, c}));
            end
        end
        col_widths = col_widths + 2; % Add 2 spaces for padding.

        % Print the formatted table.
        header_line = strjoin(arrayfun(@(c) format_text(header_parts{c}, col_widths(c), alignments{c}), 1:numel(header_parts), 'UniformOutput', false), '|');
        fprintf('%s\n', header_line);
        fprintf('%s\n', repmat('-', 1, strlength(header_line)));

        for r = 1:size(table_data, 1)
            row_line = strjoin(arrayfun(@(c) format_text(table_data{r, c}, col_widths(c), alignments{c}), 1:numel(header_parts), 'UniformOutput', false), '|');
            fprintf('%s\n', row_line);
        end
        fprintf('%s\n', repmat('-', 1, strlength(header_line)));
    end

    %% CSV Output
    % Save BCa factors as a CSV file
    fprintf(['\n' lang.bca.saving_csv '\n']);
    [~, fName, fExt] = fileparts(lang.files.bca_factors);
    fName = strrep(fName, '%s_', ''); 
    csv_filename = fullfile(csv_dir, [fName, '_', ts, fExt]);

    try
        % Attempt to open file for writing
        fileID = fopen(csv_filename, 'w');
        if fileID == -1
            error(lang.errors.file_open_error, csv_filename); 
        end
        
        % Header for the CSV file
        header = {'Metric', lang.csv.headers.effect_size, lang.csv.headers.correction_factor, ...
                  lang.csv.headers.median, lang.csv.headers.mean, ...
                  lang.csv.headers.min, lang.csv.headers.max};
        fprintf(fileID, '%s,%s,%s,%s,%s,%s,%s\n', header{:});
        
        % Loop over all metrics to write data
        for metric_idx = 1:numel(metric_names)
            metric_name = metric_names{metric_idx};
            
            % Extract data for the current metric
            z0_d = z0_d_all(:, metric_idx);
            a_d  = a_d_all(:, metric_idx);
            z0_r = z0_r_all(:, metric_idx);
            a_r  = a_r_all(:, metric_idx);
            
            % Write rows to the CSV file with explicit sign formatting
            fprintf(fileID, '%s,Cliff''s Delta,%s,%+.3f,%+.3f,%+.3f,%+.3f\n', ...
                metric_name, lang.bca.correction_factors.bias, median(z0_d), mean(z0_d), min(z0_d), max(z0_d));
            fprintf(fileID, '%s,Cliff''s Delta,%s,%+.3f,%+.3f,%+.3f,%+.3f\n', ...
                metric_name, lang.bca.correction_factors.skew, median(a_d), mean(a_d), min(a_d), max(a_d));
            fprintf(fileID, '%s,Relative Difference,%s,%+.3f,%+.3f,%+.3f,%+.3f\n', ...
                metric_name, lang.bca.correction_factors.bias, median(z0_r), mean(z0_r), min(z0_r), max(z0_r));
            fprintf(fileID, '%s,Relative Difference,%s,%+.3f,%+.3f,%+.3f,%+.3f\n', ...
                metric_name, lang.bca.correction_factors.skew, median(a_r), mean(a_r), min(a_r), max(a_r));
        end
        
        % Close the file successfully
        fclose(fileID);
        fprintf([lang.bca.csv_saved '\n'], csv_filename);

    catch ME
        % Safety cleanup: Ensure file is closed if an error occurs
        if exist('fileID', 'var') && fileID ~= -1
            fclose(fileID); 
        end
        fprintf([lang.errors.file_save_error '\n'], ME.message);
    end
end
