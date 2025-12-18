function results = run_ranking(userInput)
% RUN_RANKING - Main controller for the entire ranking analysis.
%
% Syntax:
%   HERA.run_ranking(userInput)           Run analysis, save files only
%   results = HERA.run_ranking(userInput) Run analysis, return results struct
%
% Description:
%   This function controls the algorithm (HERA) for evaluating and comparing data across quality metrics.
%   It is called by the 'start_ranking.m' and receives all necessary parameters in the 'userInput' structure.
%   It loads all external data, performs a complete statistical analysis, and saves all results in a dedicated output folder.
%   The entire process is designed to produce reproducible results.
%
%   It can be used in two usage scenarios:
%   1. Standalone/Batch: Called by 'HERA.start_ranking' or manually with file paths. 
%      Loads data from disk and saves reports.
%   2. Integration/Developer: Called directly with data matrices in memory. 
%      Returns a structured results object for further processing in MATLAB.
%
%
% Workflow:
%   1.  Configuration Setup: 
%       Loads defaults, validates metric names, and adjusts system settings.
%   2.  Environment Initialization: 
%       Creates output directories, sets up the log file, and saves the configuration.
%   3.  Reproducibility Setup: 
%       Configures the random number generator.
%   4.  Parallel Processing: 
%       Initializes the parallel pool for computation.
%   5.  Data Loading: 
%       Imports and validates metric data (1, 2, or 3 metrics) from specified files.
%   6.  Statistical Analysis: 
%       Calculates effect sizes, determines statistical thresholds, and computes BCa confidence intervals for all loaded metrics.
%   7.  Ranking Calculation:
%       Computes the final ranking based on the primary metric hierarchy and the selected ranking_mode (e.g., M1, M1_M2, M1_M3A, or M1_M2_M3). 
%       A sensitivity analysis is performed by repeating the ranking for all selected metric permutations.
%   8.  Consensus Ranking: 
%       If sensitivity analysis is active, a Borda count is used to determine a consensus rank.
%   9.  Stability Analysis: 
%       Assesses the rank stability of the primary hierarchy using a cluster bootstrap (which respects the selected ranking_mode).
%   10. Power Analysis: 
%       If enabled, conducts a post-hoc power analysis.
%   11. Result Aggregation & Output: 
%       Bundles all data and calls functions to generate logs, CSV tables, graphics, and PDF reports, all dynamically adapting to the number of metrics.
%   12. JSON Export: 
%       Saves all configurations, statistical calculations, and final results into a JSON file for machine readability and meta-analysis.
%   13. Finalization: 
%       Stops the log (Parallel pool remains active for subsequent runs).
%
% Inputs:
%   userInput - (struct) Configuration structure. Must contain ONE of the following data sources:
%
%       A) File-Based Mode (Standard):
%          .folderPath   - (char) Path to folder containing CSV/XLSX files.
%          .fileType     - (char) File extension (e.g., '.csv').
%          .metric_names - (cell) Names of metrics to analyze (must match filenames).
%
%       B) Direct Injection Mode (Developer):
%          .custom_data  - (cell) {M x 1} Cell array of matrices. Each cell contains an [N_Subjects x N_Datasets] matrix.
%          .metric_names - (cell) Names of metrics (must match number of matrices).
%          .system       - (struct) Optional performance settings (e.g., .target_memory, .jack_vec_limit).
%          .dataset_names- (cell, optional) Names of the columns/datasets.
%
%       Common Settings (Optional overrides):
%          .ranking_mode - (char) 'M1', 'M1_M2', 'M1_M3A', or 'M1_M2_M3' (Default: M1_M2_M3).
%          .output_dir   - (char) Path to save results (Default: current dir).
%          .num_workers  - (int) Number of parallel workers or 'auto'.
%          .ci_level     - (double) Confidence level for Intervals (Default: 0.95).
%
% Outputs:
%   results - (Optional) Struct containing detailed analysis results:
%       .final_rank            - [N x 1] Final rank for each dataset (1 = best).
%       .final_order           - [1 x N] Indices of datasets in ranked order.
%       .final_bootstrap_ranks - [N x B] Bootstrapped rank distribution for stability analysis.
%       .ci_lower_rank         - [N x 1] Lower bound of Rank CI.
%       .ci_upper_rank         - [N x 1] Upper bound of Rank CI.
%
%       .all_sig_matrices      - {1 x M} Cell array of logical matrices (significant wins).
%       .all_p_value_matrices  - {1 x M} Cell array of raw p-values (Wilcoxon).
%       .all_alpha_matrices    - {1 x M} Cell array of Holm-Bonferroni corrected alphas.
%
%       .d_vals_all            - [Pairs x M] Effect sizes (Cliff's Delta).
%       .rel_vals_all          - [Pairs x M] Relative differences.
%       .ci_d_all              - [Pairs x 2 x M] BCa Confidence Intervals for Delta.
%       .ci_r_all              - [Pairs x 2 x M] BCa Confidence Intervals for RelDiff.
%
%       .swap_details          - Struct containing logs of logic-based rank swaps.
%       .intermediate_orders   - Struct containing rankings after each stage (M1, M2...).
%
%       .borda_results         - (If sensitivity active) Struct with consensus ranking scores.
%       .power_results         - (If power active) Struct with post-hoc power estimates.
%       .all_permutation_ranks - [N x Perms] Ranks for all sensitivity permutations.
%       .selected_permutations - [Perms x Metrics] Indices of metric permutations.
%
%   Additionally, the function generates the following files in 'userInput.output_dir':
%       - Log file (.txt) with console output.
%       - CSV files with detailed results tables.
%       - JSON file with machine-readable export.
%       - (Optional) PDF Report and high-res figures.
%
% Author: Lukas von Erdmannsdorff

    %% 1-4. Initialization and Setup
    % Call submodule for Steps 1-4
    import HERA.*
    
    tic; % Start a timer to measure the total execution time.
    % Note: setup_environment calls HERA.run.setup_environment if we didn't import HERA.run.
    % We will assume HERA.run package is accessible via HERA.run.* or import.
    % Since run is a package under HERA, we can use HERA.run.<func>.
    
    [userInput, setupData] = HERA.run.setup_environment(userInput);
    
    % Ensure resources (diary) is cleaned up automatically when this function exits.
    % setup_environment enables the diary. We must ensure it closes.
    cleanupObj = onCleanup(@() diary('off'));

    %% 5. Load and Validate Data
    % Call the data loading function to import and validate the metric data.
    % This function is dynamic and handles 1, 2, or 3 metrics.
    [is_data_valid, all_data, dataset_names, num_probanden, num_datasets, pair_idx_all] = load_data(userInput, setupData.lang);
    
    % Abort if data validation fails.
    if ~is_data_valid
        diary off;
        delete(gcp('nocreate')); % Close parallel pool on error.
        return;
    end
    
    %% 5b. Define Styles based on User Input
    % Load the color and style settings for all plots based on the selected theme.
    styles = design(userInput.plot_theme, num_datasets, userInput.run_power_analysis);

    %% 6-10. Statistical Analysis and Ranking
    % Call submodule for Steps 6-10
    analysis_results = HERA.run.execute_analysis(all_data, num_probanden, dataset_names, pair_idx_all, userInput, setupData, styles);

    %% 11-12. Results and Export
    % Call submodule for Steps 11-12
    results = HERA.run.export_results(analysis_results, all_data, dataset_names, num_probanden, num_datasets, userInput, setupData, styles);

    %% 13. Finalize
    % Stop the timer and display the total execution time.
    fprintf(['\n' setupData.lang.run_ranking.ranking_finished '\n'], setupData.output_dir);
    total_duration_seconds = toc;
    minutes = floor(total_duration_seconds / 60);
    seconds = rem(total_duration_seconds, 60);
    fprintf('\nTotal execution time: %d minutes and %.2f seconds.\n', minutes, seconds);
end