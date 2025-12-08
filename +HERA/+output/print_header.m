function print_header(thresholds, shared_info)
% PRINT_HEADER - Prints the general information and thresholds to the console.
%
% Syntax:
%   HERA.output.print_header(thresholds, shared_info)
%
% Description:
%   Displays the initial header information for the HERA analysis output, including:
%   - Metric names and mean/SD introduction.
%   - Comparison details (number of datasets, probanden, etc.).
%   - Descriptions of statistical tests (Wilcoxon, Cliff's Delta, etc.).
%   - Threshold values for Cliff's Delta, relative differences, and SEM protection.
%   - Win criteria explanations.
%   - Bootstrap parameters summary.
%
% Inputs:
%   thresholds  - Struct with the thresholds for Cliff's Delta and relative differences.
%   shared_info - Struct with general information like metric and dataset names.
%
% Outputs:
%   None (prints directly to console).
%
% Author: Lukas von Erdmannsdorff

    % Unpack necessary variables
    lang = shared_info.lang;
    metric_names = shared_info.metric_names;
    num_metrics = numel(metric_names);
    pair_idx_all = shared_info.pair_idx_all;
    num_datasets = shared_info.num_datasets;
    num_probanden = shared_info.num_probanden;
    
    d_thresh = thresholds.d_thresh;
    rel_thresh_b = thresholds.rel_thresh_b;
    min_rel_thresh = thresholds.min_rel_thresh;
    
    %% Console Output: General Information and Thresholds
    metric_name_list_str = strjoin(metric_names, ', ');
    fprintf(['\n' lang.output.info.header '\n'], metric_name_list_str);
    fprintf([lang.output.info.mean_sd_intro '\n']);
    for i = 1:num_metrics
        fprintf(' -> %s \n', metric_names{i});
    end

    fprintf(['\n' lang.output.info.significance_intro '\n']);
    fprintf([lang.output.info.comparison_details '\n'],...
        size(pair_idx_all, 1), num_datasets, num_probanden, metric_name_list_str);
    fprintf([lang.output.info.wilcoxon_desc '\n']);
    fprintf([lang.output.info.cliffs_delta_desc '\n']);
    fprintf([lang.output.info.rel_diff_desc '\n']);

    fprintf(['\n' lang.output.thresholds.delta_header '\n']);
    for i = 1:num_metrics
        fprintf(' %s: %.3f\n', metric_names{i}, d_thresh(i));
    end

    fprintf(['\n' lang.output.thresholds.rel_diff_header '\n']);
    for i = 1:num_metrics
        fprintf(' %s: %.3f\n', metric_names{i}, rel_thresh_b(i));
    end

    fprintf([lang.output.thresholds.sem_protection_header]);
    fprintf(['\n-> ' lang.output.thresholds.sem_protection_desc '\n']);
    for i = 1:num_metrics
        fprintf(' %s: %.3f\n', metric_names{i}, min_rel_thresh(i));
    end

    fprintf(lang.output.info.win_criteria_header);
    fprintf(['\n 1. ' lang.output.info.win_criterion1]);
    fprintf(['\n 2. ' lang.output.info.win_criterion2]);
    fprintf(['\n 3. ' lang.output.info.win_criterion3]);

    fprintf(['\n\n' lang.output.info.bootstrap_summary_thr '\n'], shared_info.selected_B_thresholds);
    fprintf([lang.output.info.bootstrap_summary_ci '\n'], shared_info.selected_B_ci);
    fprintf([lang.output.info.bootstrap_summary_rank '\n'], shared_info.selected_B_rank);

end
