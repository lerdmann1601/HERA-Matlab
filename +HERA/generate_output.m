function generate_output(results, thresholds, config, shared_info)
% GENERATE_OUTPUT - Visualizes analysis results and saves them to csv files.
%
% Syntax:
%   generate_output(results, thresholds, config, shared_info)
%
% Description:
%   This function processes the final results of a statistical ranking analysis.
%   It dynamically adapts all console and file outputs based on the number of metrics (1, 2, or 3) and the selected 'ranking_mode'.
%
% Workflow:
%   1. Display of statistical thresholds: 
%      Outputs the thresholds used for the analysis to the console (dynamically for 1-3 metrics).
%   2. Presentation of the hierarchical ranking: 
%      Shows the step-by-step adjustment of the ranking based on the selected logic (M1, M1_M2, M1_M3A, M1_M2_M3).
%   3. Creation of the final results table: 
%      Summarizes the final ranks, confidence intervals, and metric statistics in a table (console and CSV, dynamic columns).
%   4. Creation of the sensitivity analysis table (optional): 
%      Displays the Borda consensus rank and scores for all tested permutations (console and CSV).
%   5. Saving the results and the log: 
%      Writes the detailed log of all comparisons to a separate CSV file (dynamically includes sections for M1, M2, M3A, or M3B).
%
% Inputs:
%   results     - Struct containing all analysis results, such as ranks and swap details.
%   thresholds  - Struct with the thresholds for Cliff's Delta and relative differences.
%   config      - Struct with the global configuration settings (incl. 'ranking_mode').
%   shared_info - Struct with general information like metric and dataset names.
%
% Outputs:
%   The function does not have direct return values. Instead, it generates files:
%   Results File     - A CSV file ('*_results.csv') with the final ranking table.
%   Sensitivity File - A CSV file ('*_sensitivity.csv') with the Borda/permutation results (if applicable).
%   Log File         - A CSV file ('*_log.csv') with the detailed log of all comparisons.
%
% Author: Lukas von Erdmannsdorff

arguments
    results (1,1) struct
    thresholds (1,1) struct
    config (1,1) struct
    shared_info (1,1) struct
end

%% 1. Initialization and Data Extraction
% Unpack the passed structures into local variables for easier access.
lang = shared_info.lang;
metric_names = shared_info.metric_names;

% Note: We no longer need to fully unpack everything here as most logic is moved to sub-functions.
% However, we keep the imports for clarity on what module is being used.

%% 2. Console Output: General Information and Thresholds
HERA.output.print_header(thresholds, shared_info);

%% 3-5. Console Output: Logic Steps for Metric 1, 2, and 3
HERA.output.print_logic_steps(results, thresholds, config, shared_info);

%% 6. Results CSV (_results.csv)
HERA.output.save_results(results, shared_info, metric_names);

%% Borda Results
HERA.output.save_sensitivity(results, shared_info, metric_names);

%% 7. Log CSV (_log.csv)
HERA.output.save_log(results, thresholds, config, shared_info);

end