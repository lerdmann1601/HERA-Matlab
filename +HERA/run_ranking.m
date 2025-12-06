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
%   1.  Initialization: 
%       Creates output directories and sets up log file.
%   2.  Configuration Loading: 
%       Loads all settings (incl. dynamic metric count and ranking_mode) and graphic styles.
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
%       Closes the parallel pool and stops the log.
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

%% 1. Initialization and Setup

% Input Validation & Default Injection (Developer Support)
% This ensures run_ranking works for both:
% 1. start_ranking calls (defaults already set -> ignored here)
% 2. Direct developer calls (defaults missing -> injected here)

% Import the HERA namespace to find internal functions
import HERA.*

% Ensure the input is a structure as expected.
validateattributes(userInput, {'struct'}, {'nonempty', 'scalar'}, 'run_ranking', 'userInput');

% Check if at least one data source is provided (File or Memory).
if ~isfield(userInput, 'folderPath') && ~isfield(userInput, 'custom_data')
    error('HERA:InvalidInput', 'Input structure must contain either "folderPath" (for file mode) or "custom_data" (for developer mode).');
end

% Load the specified language pack to localize all console and plot outputs.
lang = language_code(userInput.language);

% Load Central Defaults 
defaults = HERA.default(); 

% Validate Metric Names (Critical: Cannot be defaulted)
if ~isfield(userInput, 'metric_names')
    error(lang.errors.metric_names_missing);
end
num_metrics = numel(userInput.metric_names);

% Only fill fields that are missing
fields = fieldnames(defaults);
for i = 1:numel(fields)
    fieldName = fields{i};
    % Only apply default if the user/wrapper didn't provide it yet
    if ~isfield(userInput, fieldName)
        userInput.(fieldName) = defaults.(fieldName);
    end
end

% Logic Auto-Correction (Dependencies that depend on User Input)
% Ensure Alphas match metric count (if generic default was loaded)
if numel(userInput.alphas) ~= num_metrics
    % Expand scalar default to vector
    userInput.alphas = repmat(userInput.alphas(1), 1, num_metrics);
end

% Ensure Permutations exist if Sensitivity is on
if userInput.run_sensitivity_analysis && isempty(userInput.selected_permutations)
    userInput.selected_permutations = perms(1:num_metrics);
end

% Ensure compatibility for direct calls. If input has no .config, wrap it to match expected structure.
if ~isfield(userInput, 'config')
    userInput.config = userInput;
end

% Create a unique output folder name using the current timestamp.
timestamp_folder = string(datetime('now'), 'yyyyMMdd_HHmmss');
output_dir_name = "Ranking_" + timestamp_folder;
output_dir = fullfile(userInput.output_dir, output_dir_name);
mkdir(output_dir);

% Create subdirectories for organized output.
% The 'Graphics' folder is always created to store diagnostic convergence plots.
graphics_dir = fullfile(output_dir, "Graphics");
mkdir(graphics_dir);
% The 'PDF' folder is only created if full reports are requested.
if userInput.create_reports
    pdf_dir = fullfile(output_dir, "PDF");
    mkdir(pdf_dir);
else
    pdf_dir = ''; % Leave empty if not needed
end
csv_dir = fullfile(output_dir, "Output");
mkdir(csv_dir);

% Ensure resources (diary) is cleaned up automatically when this function exits.
cleanupObj = onCleanup(@() cleanup_resources());

    function cleanup_resources()
        % Helper function to close logs and pool
        diary off;
    end

% Set up the diary to log all console output to a text file.
log_filename = "Ranking_" + timestamp_folder + ".txt";
log_path = fullfile(output_dir, log_filename);
[~, base_name, ~] = fileparts(log_filename);
diary(char(log_path));
tic; % Start a timer to measure the total execution time.
fprintf('=======================\n');
fprintf('%s\n', lang.run_ranking.header);
fprintf('=======================\n');

% Save the complete user configuration to a JSON file for reproducibility.
config_path = fullfile(output_dir, "configuration.json"); 
try
    data_to_save = struct('userInput', userInput);
    json_text = jsonencode(data_to_save, 'PrettyPrint', true);     
    
    fid = fopen(config_path, 'w'); 
    if fid == -1
        error(lang.errors.file_open_error, config_path); 
    end
    
    fprintf(fid, '%s', json_text); 
    fclose(fid); 
catch ME
    if exist('fid', 'var') && fid ~= -1
        fclose(fid); 
    end
    fprintf([lang.errors.file_save_error '\n'], ME.message);
end

% Extract the configuration struct for easier access.
cfg = userInput.config;
fprintf('\n%s\n', lang.run_ranking.config_loaded);

% Display a summary of the chosen bootstrap configurations in the console.
fprintf('--------------------------------------------------\n');
fprintf([lang.run_ranking.cpu_cores '\n'], num2str(userInput.num_workers));
fprintf([lang.run_ranking.reproducibility '\n'], mat2str(userInput.reproducible), userInput.seed); pause(0.5);

% Display settings for threshold calculation.
fprintf(['\n-> ' lang.run_ranking.percentile_bootstrap '\n']); pause(0.5);
if ~isempty(cfg.manual_B_thr)
    fprintf(['   ' lang.run_ranking.method_manual '\n'], cfg.manual_B_thr); pause(0.5);
else
    cfg_thr = cfg.bootstrap_thresholds;
    fprintf(['   ' lang.run_ranking.method_auto '\n']); pause(0.5);
    if isempty(cfg_thr.smoothing_window)
        fprintf(['   ' lang.run_ranking.convergence_simple '\n']); pause(0.5);
    else
        fprintf(['   ' lang.run_ranking.convergence_robust '\n']); pause(0.5);
        fprintf(['   ' lang.run_ranking.smoothing_window '\n'], cfg_thr.smoothing_window); pause(0.5);
        fprintf(['   ' lang.run_ranking.stable_runs '\n'], cfg_thr.convergence_streak_needed); pause(0.5);
    end
    fprintf(['   ' lang.run_ranking.b_range '\n'], cfg_thr.B_start, cfg_thr.B_step, cfg_thr.B_end); pause(0.5);
    fprintf(['   ' lang.run_ranking.trials_per_b '\n'], cfg_thr.n_trials); pause(0.5);
    fprintf(['   ' lang.run_ranking.tolerance '\n'], cfg_thr.convergence_tolerance * 100); pause(0.5);
end

% Display settings for BCa confidence interval calculation.
fprintf(['\n-> ' lang.run_ranking.bca_bootstrap '\n']); pause(0.5);
if ~isempty(cfg.manual_B_ci)
    fprintf(['   ' lang.run_ranking.method_manual '\n'], cfg.manual_B_ci); pause(0.5);
else
    cfg_ci = cfg.bootstrap_ci;
    fprintf(['   ' lang.run_ranking.method_auto '\n']); pause(0.5);
    if isempty(cfg_ci.smoothing_window)
        fprintf(['   ' lang.run_ranking.convergence_simple '\n']); pause(0.5);
    else
        fprintf(['   ' lang.run_ranking.convergence_robust '\n']); pause(0.5);
        fprintf(['   ' lang.run_ranking.smoothing_window '\n'], cfg_ci.smoothing_window); pause(0.5);
        fprintf(['   ' lang.run_ranking.stable_runs '\n'], cfg_ci.convergence_streak_needed); pause(0.5);
    end
    fprintf(['   ' lang.run_ranking.b_range '\n'], cfg_ci.B_start, cfg_ci.B_step, cfg_ci.B_end); pause(0.5);
    fprintf(['   ' lang.run_ranking.trials_per_b '\n'], cfg_ci.n_trials); pause(0.5);
    fprintf(['   ' lang.run_ranking.tolerance '\n'], cfg_ci.convergence_tolerance * 100); pause(0.5);
end

% Display settings for rank stability analysis.
fprintf(['\n-> ' lang.run_ranking.cluster_bootstrap '\n']); pause(0.5);
if ~isempty(cfg.manual_B_rank)
    fprintf(['   ' lang.run_ranking.method_manual '\n'], cfg.manual_B_rank); pause(0.5);
else
    cfg_rank = cfg.bootstrap_ranks;
    fprintf(['   ' lang.run_ranking.method_auto '\n']); pause(0.5);
    if isempty(cfg_rank.smoothing_window)
        fprintf(['   ' lang.run_ranking.convergence_simple '\n']); pause(0.5);
    else
        fprintf(['   ' lang.run_ranking.convergence_robust '\n']); pause(0.5);
        fprintf(['   ' lang.run_ranking.smoothing_window '\n'], cfg_rank.smoothing_window); pause(0.5);
        fprintf(['   ' lang.run_ranking.stable_runs '\n'], cfg_rank.convergence_streak_needed); pause(0.5);
    end
    fprintf(['   ' lang.run_ranking.b_range '\n'], cfg_rank.B_start, cfg_rank.B_step, cfg_rank.B_end); pause(0.5);
    fprintf(['   ' lang.run_ranking.trials_per_b '\n'], cfg_rank.n_trials); pause(0.5);
    fprintf(['   ' lang.run_ranking.tolerance '\n'], cfg_rank.convergence_tolerance * 100); pause(0.5);
end

fprintf('--------------------------------------------------\n');

% Transfer user-defined settings into local variables.
config = userInput.config;
config.timestamp = char(timestamp_folder); 
metric_names = userInput.metric_names;
config.metric_names = metric_names; 
config.ranking_mode = userInput.ranking_mode; 
num_metrics = numel(metric_names); 

%% 2. Reproducibility of the Random Number Generator
s = []; % Initialize the random stream variable.
if userInput.reproducible
    % Use a fixed seed for reproducible results.
    s = RandStream('mlfg6331_64', 'Seed', userInput.seed);
    fprintf(['\n' lang.run_ranking.rng_reproducible '\n'], userInput.seed);
else
    % Use a shuffled seed for non-deterministic results.
    s = RandStream('mlfg6331_64', 'Seed', 'shuffle');
    fprintf(['\n' lang.run_ranking.rng_shuffle '\n']);
end
RandStream.setGlobalStream(s); % Set the configured stream as the global default.

%% 3. Initialize Parallel Processing
% Check if a pool is already running to avoid expensive restarts in batch mode.
currentPool = gcp('nocreate');

if isempty(currentPool)
    % No pool exists -> Start one according to user config.
    if ischar(userInput.num_workers) && strcmp(userInput.num_workers, 'auto')
        % Automatically determine optimal core count.
        n_cores = max(1, parcluster('local').NumWorkers);
        fprintf([lang.run_ranking.parallel_start_auto '\n'], n_cores);
        currentPool = parpool(n_cores); % Capture the pool object!
    else
        % Use user-specified number.
        n_cores = userInput.num_workers;
        fprintf([lang.run_ranking.parallel_start_manual '\n'], n_cores);
        currentPool = parpool(n_cores); % Capture the pool object!
    end
else
    % Pool exists -> Reuse it for performance.
    fprintf([lang.run_ranking.pool_active_skip '\n'], currentPool.NumWorkers);
    % Optional: Check if worker count matches request (warn only, don't kill).
    if isnumeric(userInput.num_workers) && currentPool.NumWorkers ~= userInput.num_workers
        fprintf([lang.run_ranking.pool_size_mismatch '\n'], currentPool.NumWorkers, userInput.num_workers);
    end
end

% Explicitly tell the parallel workers where the +HERA package is located.
if ~isempty(currentPool)
    try
        % Get the directory containing this file (which is inside +HERA)
        [functionPath, ~, ~] = fileparts(mfilename('fullpath'));
        % We need the Parent directory of +HERA to attach the package correctly
        packageParentDir = fileparts(functionPath); 
        
        % Ensure it is on the Client path
        addpath(packageParentDir);
        
        % Force all workers to add this path AND WAIT until they are done
        F = parfevalOnAll(currentPool, @addpath, 0, packageParentDir);
        wait(F); 
        
    catch ME
        fprintf([lang.run_ranking.path_warning '\n'], ME.message);
    end
end

%% 4. Load and Validate Data
% Call the data loading function to import and validate the metric data.
% This function is already dynamic and handles 1, 2, or 3 metrics.
[is_data_valid, all_data, dataset_names, num_probanden, num_datasets, pair_idx_all] = load_data(userInput, lang);

% Abort if data validation fails.
if ~is_data_valid
    diary off;
    delete(gcp('nocreate')); % Close parallel pool on error.
    return;
end

% config.metric_names = userInput.metric_names; % Redundant, already set in Sec 1

%% 5. Define Styles based on User Input
% Load the color and style settings for all plots based on the selected theme.
styles = design(userInput.plot_theme, num_datasets, userInput.run_power_analysis);

%% 6. Statistical Analysis and CI Calculation
% Determine the statistical thresholds for effect sizes (Cliff's Delta, RelDiff).
% This function loops 1:numel(all_data), so it's already dynamic.
fprintf(['\n' lang.run_ranking.calc_thresholds '\n']);
[d_thresh, rel_thresh, rel_thresh_b, min_rel_thresh, d_vals_all, rel_vals_all, pair_idx_all, selected_B, stability_data_thr, ...
    h_fig_thr_global, h_fig_thr_detailed, h_fig_hist_thr, h_fig_hist_raw] = ...
     calculate_thresholds(all_data, num_probanden, config, graphics_dir, config.manual_B_thr, s, styles, lang);

% Calculate the Bias-Corrected and accelerated (BCa) confidence intervals.
% This function is also dynamic based on numel(all_data).
fprintf(['\n' lang.run_ranking.calc_bca '\n']);
[selected_B_ci, ci_d_all, ci_r_all, z0_d_all, a_d_all, z0_r_all, a_r_all, stability_data_ci, h_fig_bca_global, h_fig_bca_detailed, ...
    h_fig_hist_z0, h_fig_hist_a, h_fig_hist_widths] = ...
    calculate_bca_ci(all_data, d_vals_all, rel_vals_all, pair_idx_all, ...
        num_probanden, config, metric_names, graphics_dir, csv_dir, config.manual_B_ci, s, styles, lang, base_name);

% Bundle results into structs for cleaner passing to other functions.
effect_sizes = struct('d_vals_all', d_vals_all, 'rel_vals_all', rel_vals_all);
thresholds   = struct('d_thresh', d_thresh, 'rel_thresh', rel_thresh, 'rel_thresh_b', rel_thresh_b, 'min_rel_thresh', min_rel_thresh);

%% 7. Ranking Calculation (+ Sensitivity Analysis)
% The ranking is calculated for the primary hierarchy and all other selected permutations.
num_permutations = size(userInput.selected_permutations, 1);
all_permutation_ranks = zeros(num_datasets, num_permutations);
fprintf(['\n' lang.run_ranking.calc_ranking_permutations '\n'], num_permutations);

for p = 1:num_permutations
    % Define the primary hierarchy from the first permutation entry.
    primary_hierarchy_indices = userInput.selected_permutations(1, :);
    % Get the current metric order for this iteration.
    current_permutation = userInput.selected_permutations(p, :);

    % Find the mapping from the current permutation back to the primary order.
    % This correctly handles permutations of 1, 2, or 3 metrics.
    [~, local_indices] = ismember(current_permutation, primary_hierarchy_indices);
    
    % Reorder all necessary data and configurations for the current ranking run.
    perm_metric_names = metric_names(local_indices);
    perm_all_data = all_data(local_indices);
    perm_d_vals = d_vals_all(:, local_indices);
    perm_rel_vals = rel_vals_all(:, local_indices);
    perm_d_thresh = d_thresh(local_indices);
    perm_rel_thresh = rel_thresh(local_indices);
    perm_config = config;
    perm_config.alphas = config.alphas(local_indices);
    perm_config.ranking_mode = userInput.ranking_mode; 
    
    fprintf([' -> ' lang.run_ranking.permutation_progress '\n'], p, num_permutations, strjoin(perm_metric_names, ' -> '));
    
    % Prepare the reordered data structs for the ranking function.
    perm_effect_sizes = struct('d_vals_all', perm_d_vals, 'rel_vals_all', perm_rel_vals);
    perm_thresholds = struct('d_thresh', perm_d_thresh, 'rel_thresh', perm_rel_thresh);
    
    % Calculate the ranking for the current permutation and store the result.
    [~, single_run_rank] = calculate_ranking(perm_all_data, perm_effect_sizes, perm_thresholds, perm_config, dataset_names, pair_idx_all);
    all_permutation_ranks(:, p) = single_run_rank;

    % The first permutation (p=1) is the primary hierarchy; store its detailed results.
    if p == 1
         % Call with the original (primary) data and config (which now includes ranking_mode)
         [final_order, final_rank, all_sig_matrices, all_alpha_matrices, all_p_value_matrices, swap_details, intermediate_orders] = ...
             calculate_ranking(all_data, effect_sizes, thresholds, config, dataset_names, pair_idx_all);
    end
end

%% 8. Borda Ranking Calculation
% If sensitivity analysis is enabled, calculate the consensus rank using the Borda method.
if userInput.run_sensitivity_analysis
    fprintf(['\n' lang.run_ranking.calc_borda '\n']);
    borda_results = borda_ranking(all_permutation_ranks, dataset_names);
else
    borda_results = []; % Set to empty if not performed.
end

%% 9. Bootstrap Analysis of Ranks (for Primary Hierarchy)
% Assess the stability of the primary ranking using a cluster bootstrap.
% This function now receives the 'config' struct which contains the 'ranking_mode'.
fprintf(['\n' lang.run_ranking.bootstrap_ranks '\n']);
% This function now receives the 'config' struct which contains the 'ranking_mode'.
fprintf(['\n' lang.run_ranking.bootstrap_ranks '\n']);
[final_bootstrap_ranks, selected_B_rank, stability_data_rank, h_figs_rank, h_fig_hist_rank, ci_lower_rank, ci_upper_rank] = ...
    bootstrap_ranking(all_data, thresholds, config, dataset_names, final_rank, pair_idx_all, num_probanden, graphics_dir, csv_dir,...
    config.manual_B_rank, s, styles, lang, base_name);

%% 10. Power Analysis
% If enabled, perform a post-hoc, non-parametric power analysis.
if userInput.run_power_analysis
    fprintf(['\n' lang.run_ranking.calc_power '\n']);
    power_results = power_analysis(all_data, config, thresholds, num_probanden, userInput.power_simulations, pair_idx_all, s, lang);
else
    power_results = []; % Set to empty if not performed.
end

%% 11. Result Aggregation and Output
fprintf(['\n' lang.run_ranking.creating_output '\n']);
% Calculate mean and standard deviation of the raw metrics for the final table.
% This loop is already dynamic based on num_metrics.
mean_metrics = zeros(num_datasets, num_metrics);
std_metrics = zeros(num_datasets, num_metrics);
for i = 1:num_metrics
    mean_metrics(:, i) = mean(all_data{i}, 1)';
    std_metrics(:, i) = std(all_data{i}, 0, 1)';
end

% Aggregate all analysis results into a single 'results' struct.
results = struct();
results.final_rank = final_rank;
results.ci_lower_rank = ci_lower_rank;
results.ci_upper_rank = ci_upper_rank;
results.final_bootstrap_ranks = final_bootstrap_ranks;
results.all_sig_matrices = all_sig_matrices;
results.all_alpha_matrices = all_alpha_matrices;
results.all_p_value_matrices = all_p_value_matrices;
results.d_vals_all = d_vals_all;
results.rel_vals_all = rel_vals_all;
results.ci_d_all = ci_d_all;
results.ci_r_all = ci_r_all;
results.swap_details = swap_details;
results.intermediate_orders = intermediate_orders;
results.final_order = final_order;
results.borda_results = borda_results;
results.power_results = power_results;
results.all_permutation_ranks = all_permutation_ranks;
results.selected_permutations = userInput.selected_permutations;

% Aggregate all shared information (paths, names, etc.) into a 'shared_info' struct.
shared_info = struct();
shared_info.output_dir = output_dir;
shared_info.graphics_dir = graphics_dir;
shared_info.pdf_dir = pdf_dir;
shared_info.csv_dir = csv_dir;
shared_info.log_basename = base_name;
shared_info.dataset_names = dataset_names;
shared_info.available_metrics = userInput.available_metrics;
shared_info.metric_names = metric_names;
shared_info.num_datasets = num_datasets;
shared_info.num_probanden = num_probanden;
shared_info.pair_idx_all = pair_idx_all;
shared_info.all_data = all_data;
shared_info.mean_metrics = mean_metrics;
shared_info.std_metrics = std_metrics;
shared_info.z0_d_all = z0_d_all;
shared_info.a_d_all = a_d_all;
shared_info.z0_r_all = z0_r_all;
shared_info.a_r_all = a_r_all;
shared_info.alphas = config.alphas;
shared_info.selected_B_thresholds = selected_B;
shared_info.selected_B_ci = selected_B_ci;
shared_info.selected_B_rank = selected_B_rank;
shared_info.h_fig_thr_global = h_fig_thr_global;
shared_info.h_fig_thr_detailed = h_fig_thr_detailed;
shared_info.h_fig_bca_global = h_fig_bca_global;
shared_info.h_fig_bca_detailed = h_fig_bca_detailed;
shared_info.h_fig_rank_conv = h_figs_rank;
shared_info.h_fig_hist_thr = h_fig_hist_thr;
shared_info.h_fig_hist_raw = h_fig_hist_raw;
shared_info.h_fig_hist_z0 = h_fig_hist_z0;
shared_info.h_fig_hist_a = h_fig_hist_a;
shared_info.h_fig_hist_widths = h_fig_hist_widths;
shared_info.h_fig_hist_rank = h_fig_hist_rank;
shared_info.plot_theme = userInput.plot_theme;
shared_info.lang = lang; 
shared_info.config = config; 


% Call functions to generate all final outputs.
% These functions receive 'config' or 'shared_info.config' and will adapt.
generate_output(results, thresholds, config, shared_info);

% Generate comprehensive graphical reports only if requested.
% Diagnostic plots (convergence) are already saved by the calculation functions.
if userInput.create_reports
    generate_plots(results, thresholds, shared_info, styles);
end

%% 12. JSON Export for Machine Readability and AI Training
% This block saves all critical data into a single JSON file. 
% This creates a record of the entire analysis, ideal for downstream processing, meta-analysis, or as a structured dataset for training AI models.
fprintf(['\n' lang.run_ranking.json_saving '\n']);
try
    % Bundle all relevant data structures into a single export struct.

    % dataset_names: Contains index of dataset names
    % This represents the "original" loaded index of datasets
    json_export_data.dataset_names = shared_info.dataset_names;
    
    % config: Contains all initial settings (alphas, CI level, bootstrap params, ranking_mode)
    % This represents the "input configuration" of the analysis.
    json_export_data.config = config;
    
    % meta: Contains key execution parameters (metadata) for full reproducibilityand context for interpreting the results data.
    json_export_data.meta = struct();
    json_export_data.meta.n_subjects = shared_info.num_probanden;
    json_export_data.meta.n_datasets = shared_info.num_datasets;
    json_export_data.meta.version = HERA.get_version(); % Automatic version detection
    json_export_data.meta.timestamp = shared_info.log_basename;
    
    % pair_indices: Maps the N-row pairwise arrays (e.g., in stats_calcs) to their respective [Dataset A, Dataset B] indices.
    json_export_data.meta.pair_indices = shared_info.pair_idx_all;
    
    % metric_list: The master list of all metrics found in the data folder, in their original load order. 
    % This list is the target for the indices used in 'results.selected_permutations'.
    json_export_data.meta.metric_list = shared_info.available_metrics;
    
    % bootstrap_B: The final number of iterations (B) selected by the convergence analysis (or manually set) for each of the three steps.
    json_export_data.meta.stability_analysis = struct();
    json_export_data.meta.stability_analysis.thresholds = stability_data_thr;
    json_export_data.meta.stability_analysis.ci = stability_data_ci;
    json_export_data.meta.stability_analysis.ranks = stability_data_rank;
    json_export_data.meta.bootstrap_B = struct();
    json_export_data.meta.bootstrap_B.thresholds = shared_info.selected_B_thresholds;
    json_export_data.meta.bootstrap_B.ci = shared_info.selected_B_ci;
    json_export_data.meta.bootstrap_B.ranks = shared_info.selected_B_rank;

    % stats: Contains basic descriptive statistics (Mean, StdDev) for each dataset, corresponding to the original metric order.
    json_export_data.stats = struct();
    json_export_data.stats.mean = shared_info.mean_metrics;
    json_export_data.stats.std = shared_info.std_metrics;

    % thresholds: Contains all calculated statistical thresholds (d_thresh, rel_thresh, etc.).
    % This represents one part of the "statistical basis" for the decisions.
    json_export_data.thresholds = thresholds;
    
    % Restructuring for optimal JSON output:
    % The 'results' struct (from Section 11) is split into two parts for a clean separation of statistics and ranking targets.
    
    % stats_calcs: Contains all pairwise statistical calculations (Effect Sizes, CIs, p-values).
    % These are the "features" for downstream analysis.
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

    % results: Contains all final ranking logic and high-level results.
    % These are the "targets" for downstream analysis.
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
    
    % Note: shared_info.all_data (the raw subject data) is intentionally excluded to keep the JSON file size manageable. 
    % The relevant features (effect sizes, p-values) are already contained in 'stats_calcs'.

    % Generate JSON filename.
    [~, fName, fExt] = fileparts(lang.files.json_analysis_data);
    fName = strrep(fName, '%s_', '');
    json_filename = fullfile(csv_dir, [fName, '_', char(timestamp_folder), fExt]);
    
    % Encode the main struct into a readable JSON string.
    json_text = jsonencode(json_export_data, 'PrettyPrint', true);
    
    % Write the JSON string to the file.
    fid = fopen(json_filename, 'w');
    if fid == -1
        % Throw an error if the file cannot be opened.
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

%% 13. Finalize
% Stop the timer and display the total execution time.
fprintf(['\n' lang.run_ranking.ranking_finished '\n'], output_dir);
total_duration_seconds = toc;
minutes = floor(total_duration_seconds / 60);
seconds = rem(total_duration_seconds, 60);
fprintf('\nTotal execution time: %d minutes and %.2f seconds.\n', minutes, seconds);

% Close the parallel pool and stop logging.
fprintf(['\n' lang.run_ranking.closing_pool '\n']);
delete(gcp('nocreate'));
diary off;
end

function lang = language_code(language_code)
    % A simplified helper function to load the language JSON file.
    % This is used in non-interactive contexts like the main run script.
    % Get path relative to this function.
    % This works in both MATLAB and deployed mode as the file is extracted.
    base_path = fileparts(mfilename('fullpath'));
    file_path = fullfile(base_path, 'language', [language_code, '.json']);
    if ~exist(file_path, 'file')
        error('Language file for code "%s" not found.', language_code);
    end
    json_text = fileread(file_path);
    lang = jsondecode(json_text);
end