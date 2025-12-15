function results = export_results(analysis_results, all_data, dataset_names, num_probanden, num_datasets, userInput, setupData, styles)
% EXPORT_RESULTS - Aggregate results and generate outputs.
%
% Syntax:
%   results = export_results(analysis_results, all_data, dataset_names, num_probanden, num_datasets, userInput, setupData, styles)
%
% Description:
%   This function handles Steps 11-12 of the ranking analysis:
%   11. Result Aggregation & Output (CSV, Plots, PDF)
%   12. JSON Export
%
% Inputs:
%   analysis_results - Struct from execute_analysis containing calculation results.
%   all_data         - Raw metric data.
%   dataset_names    - Names of datasets.
%   num_probanden    - Number of subjects.
%   num_datasets     - Number of datasets.
%   userInput        - Original user input.
%   setupData        - Environment setup data.
%   styles           - Plot styles.
%
% Outputs:
%   results          - Final results structure.
%
% Author: Lukas von Erdmannsdorff

    import HERA.*

    % Unpack setupData
    lang = setupData.lang;
    config = setupData.config;
    output_dir = setupData.output_dir;
    graphics_dir = setupData.graphics_dir;
    if isfield(setupData, 'pdf_dir')
        pdf_dir = setupData.pdf_dir;
    else
        % Fallback if not set (e.g. create_reports false)
        pdf_dir = ''; 
    end
    csv_dir = setupData.csv_dir;
    base_name = setupData.base_name;
    % setupData.config should be used, but userInput.available_metrics might be needed.
    
    metric_names = config.metric_names; 
    num_metrics = numel(metric_names);

    %% 11. Result Aggregation and Output
    fprintf(['\n' lang.run_ranking.creating_output '\n']);
    
    % Calculate mean and standard deviation of the raw metrics for the final table.
    mean_metrics = zeros(num_datasets, num_metrics);
    std_metrics = zeros(num_datasets, num_metrics);
    for i = 1:num_metrics
        mean_metrics(:, i) = mean(all_data{i}, 1)';
        std_metrics(:, i) = std(all_data{i}, 0, 1)';
    end
    
    % Aggregate all analysis results into a single 'results' struct.
    % We pull fields from analysis_results
    results = struct();
    results.final_rank = analysis_results.final_rank;
    results.ci_lower_rank = analysis_results.ci_lower_rank;
    results.ci_upper_rank = analysis_results.ci_upper_rank;
    results.final_bootstrap_ranks = analysis_results.final_bootstrap_ranks;
    results.all_sig_matrices = analysis_results.all_sig_matrices;
    results.all_alpha_matrices = analysis_results.all_alpha_matrices;
    results.all_p_value_matrices = analysis_results.all_p_value_matrices;
    results.d_vals_all = analysis_results.d_vals_all;
    results.rel_vals_all = analysis_results.rel_vals_all;
    results.ci_d_all = analysis_results.ci_d_all;
    results.ci_r_all = analysis_results.ci_r_all;
    results.swap_details = analysis_results.swap_details;
    results.intermediate_orders = analysis_results.intermediate_orders;
    results.final_order = analysis_results.final_order;
    results.borda_results = analysis_results.borda_results;
    results.power_results = analysis_results.power_results;
    results.all_permutation_ranks = analysis_results.all_permutation_ranks;
    results.selected_permutations = userInput.selected_permutations;
    
    % Aggregate all shared information (paths, names, etc.) into a 'shared_info' struct.
    shared_info = struct();
    shared_info.output_dir = output_dir;
    shared_info.graphics_dir = graphics_dir;
    shared_info.pdf_dir = pdf_dir;
    shared_info.csv_dir = csv_dir;
    shared_info.log_basename = base_name;
    shared_info.dataset_names = dataset_names;
    
    % Handle 'available_metrics' for metadata export.
    % In File-Based Mode (via start_ranking), this contains all metrics found in the folder.
    % In Direct-Developer Mode, this field might be missing. Fallback to the metrics actually analyzed.
    if isfield(userInput, 'available_metrics')
        shared_info.available_metrics = userInput.available_metrics;
    else
        shared_info.available_metrics = metric_names;
    end
    
    shared_info.metric_names = metric_names;
    shared_info.num_datasets = num_datasets;
    shared_info.num_probanden = num_probanden;
    shared_info.pair_idx_all = analysis_results.pair_idx_all;
    shared_info.all_data = all_data;
    shared_info.mean_metrics = mean_metrics;
    shared_info.std_metrics = std_metrics;
    shared_info.z0_d_all = analysis_results.z0_d_all;
    shared_info.a_d_all = analysis_results.a_d_all;
    shared_info.z0_r_all = analysis_results.z0_r_all;
    shared_info.a_r_all = analysis_results.a_r_all;
    shared_info.alphas = config.alphas;
    
    shared_info.selected_B_thresholds = analysis_results.selected_B;
    shared_info.selected_B_ci = analysis_results.selected_B_ci;
    shared_info.selected_B_rank = analysis_results.selected_B_rank;
    
    % Unpack handles
    handles = analysis_results.handles;
    shared_info.h_fig_thr_global = handles.h_fig_thr_global;
    shared_info.h_fig_thr_detailed = handles.h_fig_thr_detailed;
    shared_info.h_fig_bca_global = handles.h_fig_bca_global;
    shared_info.h_fig_bca_detailed = handles.h_fig_bca_detailed;
    shared_info.h_fig_rank_conv = handles.h_figs_rank;
    shared_info.h_fig_hist_thr = handles.h_fig_hist_thr;
    shared_info.h_fig_hist_raw = handles.h_fig_hist_raw;
    shared_info.h_fig_hist_z0 = handles.h_fig_hist_z0;
    shared_info.h_fig_hist_a = handles.h_fig_hist_a;
    shared_info.h_fig_hist_widths = handles.h_fig_hist_widths;
    shared_info.h_fig_hist_rank = handles.h_fig_hist_rank;
    
    shared_info.plot_theme = userInput.plot_theme;
    shared_info.lang = lang; 
    shared_info.config = config; 
    
    % Call functions to generate all final outputs.
    generate_output(results, analysis_results.thresholds, config, shared_info);
    
    % Generate comprehensive graphical reports only if requested.
    if userInput.create_reports
        generate_plots(results, analysis_results.thresholds, shared_info, styles);
    else
        %If reports are disabled, generate_plots is not called (which closes figures).
        % We must close them manually here to prevent RAM leaks in batch mode.
        handles = analysis_results.handles;
        all_handles = [handles.h_fig_thr_global, handles.h_fig_thr_detailed, ...
                       handles.h_fig_bca_global, handles.h_fig_bca_detailed, ...
                       handles.h_fig_hist_thr, handles.h_fig_hist_raw, ...
                       handles.h_fig_hist_z0, handles.h_fig_hist_a, ...
                       handles.h_fig_hist_widths, handles.h_figs_rank, ...
                       handles.h_fig_hist_rank];
                   
        % Filter specifically for valid graphics objects to avoid errors with empty/structs
        valid_handles = all_handles(isgraphics(all_handles));
        close(valid_handles);
    end
    
    %% 12. JSON Export for Machine Readability and AI Training
    fprintf(['\n' lang.run_ranking.json_saving '\n']);
    try
        % Bundle all relevant data structures into a single export struct.
    
        % dataset_names: Contains index of dataset names
        json_export_data.dataset_names = shared_info.dataset_names;
        
        % config: Contains all initial settings
        json_export_data.config = config;

        % Remove target_memory from saved config
        if isfield(json_export_data.config, 'system') && isfield(json_export_data.config.system, 'target_memory')
            json_export_data.config.system = rmfield(json_export_data.config.system, 'target_memory');
        end
        
        % meta: Contains key execution parameters
        json_export_data.meta = struct();
        json_export_data.meta.n_subjects = shared_info.num_probanden;
        json_export_data.meta.n_datasets = shared_info.num_datasets;
        json_export_data.meta.version = HERA.get_version(); 
        json_export_data.meta.timestamp = shared_info.log_basename;
        
        json_export_data.meta.pair_indices = shared_info.pair_idx_all;
        json_export_data.meta.metric_list = shared_info.available_metrics;
        
        % bootstrap_B
        json_export_data.meta.stability_analysis = struct();
        json_export_data.meta.stability_analysis.thresholds = analysis_results.stability_data_thr;
        json_export_data.meta.stability_analysis.ci = analysis_results.stability_data_ci;
        json_export_data.meta.stability_analysis.ranks = analysis_results.stability_data_rank;
        json_export_data.meta.bootstrap_B = struct();
        json_export_data.meta.bootstrap_B.thresholds = shared_info.selected_B_thresholds;
        json_export_data.meta.bootstrap_B.ci = shared_info.selected_B_ci;
        json_export_data.meta.bootstrap_B.ranks = shared_info.selected_B_rank;
    
        % stats
        json_export_data.stats = struct();
        json_export_data.stats.mean = shared_info.mean_metrics;
        json_export_data.stats.std = shared_info.std_metrics;
    
        % thresholds
        json_export_data.thresholds = analysis_results.thresholds;
        
        % stats_calcs
        stats_calcs_export = struct();
        stats_calcs_export.d_vals_all = results.d_vals_all;
        stats_calcs_export.rel_vals_all = results.rel_vals_all;
        stats_calcs_export.ci_d_all = results.ci_d_all;
        stats_calcs_export.ci_r_all = results.ci_r_all;
        stats_calcs_export.z0_d_all = shared_info.z0_d_all;
        stats_calcs_export.a_d_all = shared_info.a_d_all;
        stats_calcs_export.z0_r_all = shared_info.z0_r_all;
        stats_calcs_export.a_r_all = shared_info.a_r_all;
        stats_calcs_export.all_sig_matrices = results.all_sig_matrices;
        stats_calcs_export.all_alpha_matrices = results.all_alpha_matrices;
        stats_calcs_export.all_p_value_matrices = results.all_p_value_matrices;
        json_export_data.stats_calcs = stats_calcs_export;
    
        % results export
        results_export = struct();
        results_export.final_rank = results.final_rank;
        results_export.ci_lower_rank = results.ci_lower_rank;
        results_export.ci_upper_rank = results.ci_upper_rank;
        results_export.final_order = results.final_order;
        results_export.intermediate_orders = results.intermediate_orders;
        results_export.swap_details = results.swap_details;
        results_export.final_bootstrap_ranks = results.final_bootstrap_ranks;
        results_export.borda_results = results.borda_results;
        results_export.power_results = results.power_results;
        results_export.all_permutation_ranks = results.all_permutation_ranks;
        results_export.selected_permutations = results.selected_permutations;
        json_export_data.results = results_export;
        
        % Generate JSON filename.
        [~, fName, fExt] = fileparts(lang.files.json_analysis_data);
        fName = strrep(fName, '%s_', '');
        json_filename = fullfile(csv_dir, [fName, '_', char(setupData.timestamp), fExt]);
        
        % Encode the main struct into a readable JSON string.
        json_text = jsonencode(json_export_data, 'PrettyPrint', true);
        
        % Write the JSON string to the file.
        fid = fopen(json_filename, 'w');
        if fid == -1
            error(lang.errors.file_open_error, json_filename); 
        end
        fprintf(fid, '%s', json_text);
        fclose(fid);
        
        fprintf([lang.run_ranking.json_save_success '\n'], json_filename);
    
    catch ME
        % Catch and display any errors during the JSON export.
        if exist('fid', 'var') && fid ~= -1
            fclose(fid);
        end
        fprintf([lang.errors.file_save_error_details '\n'], ME.message); 
    end
end
