function defaults = default()
% DEFAULT - Returns the global default configuration for HERA.
%
% Syntax:
%   defaults = HERA.default()
%
% Description:
%   This function serves as the "Single Source of Truth" for all configuration
%   parameters. It is used by:
%   1. HERA.start_ranking (to populate the interactive GUI defaults)
%   2. HERA.run_ranking   (to fill missing fields in developer mode)
%
%   Centralizing these values ensures consistency across batch, interactive,
%   and developer usage scenarios.
%
% Outputs:
%   defaults - (struct) A structure containing all necessary configuration
%              fields with their standard values.
%
% Author:   Lukas von Erdmannsdorff
% Date:     12.10.2025
% Version:  1

    defaults = struct();
    
    %% General Settings
    defaults.reproducible = true;
    defaults.seed = 123;
    defaults.num_workers = 'auto';
    defaults.create_reports = true; 
    defaults.plot_theme = 'light';
    defaults.language = 'en';    
    
    %% Statistical Parameters
    defaults.ci_level = 0.95;
    % Fallback for alphas, will be dynamically expanded to match metric count later
    defaults.alphas = [0.05, 0.05, 0.05]; 
    
    %% Logic & Analysis Settings
    % Default ranking logic (Full hierarchy)
    defaults.ranking_mode = 'M1_M2_M3';
    defaults.run_sensitivity_analysis = true;
    defaults.run_power_analysis = true;
    defaults.power_simulations = 10000;
    
    % Data quality settings
    % Warning threshold for missing data (0.80 = 80% valid pairs required)
    defaults.min_data_completeness = 0.80;

    %% Bootstrap Configuration
    % Manual Bootstrap Defaults (Empty = Automatic Mode)
    defaults.manual_B_thr = 2000;
    defaults.manual_B_ci = 5000;
    defaults.manual_B_rank = 500;

    
    % Thresholds (Percentile Bootstrap)
    cfg_thr = struct();
    cfg_thr.B_start = 100;
    cfg_thr.B_step = 100;
    cfg_thr.B_end = 10000;
    cfg_thr.n_trials = 25;
    cfg_thr.min_steps_for_convergence_check = 1;
    cfg_thr.convergence_tolerance = 0.005;
    cfg_thr.smoothing_window = 3;
    cfg_thr.convergence_streak_needed = 3;
    defaults.bootstrap_thresholds = cfg_thr;
    
    % Confidence Intervals (BCa Bootstrap) - Needs more iterations
    cfg_ci = struct();
    cfg_ci.B_start = 100;
    cfg_ci.B_step = 200;
    cfg_ci.B_end = 20000;
    cfg_ci.n_trials = 25;
    cfg_ci.min_steps_for_convergence_check = 1;
    cfg_ci.convergence_tolerance = 0.01;
    cfg_ci.smoothing_window = 4;
    cfg_ci.convergence_streak_needed = 3;
    defaults.bootstrap_ci = cfg_ci;
    
    % Rank Stability (Cluster Bootstrap) - Discrete distribution converges faster
    cfg_rank = struct();
    cfg_rank.B_start = 50;
    cfg_rank.B_step = 10;
    cfg_rank.B_end = 1500;
    cfg_rank.n_trials = 15;
    cfg_rank.min_steps_for_convergence_check = 1;
    cfg_rank.convergence_tolerance = 0.005;
    cfg_rank.smoothing_window = 3;
    cfg_rank.convergence_streak_needed = 3;
    defaults.bootstrap_ranks = cfg_rank;

    %% Placeholders for Input Data
    % These fields are initialized empty and must be filled by the user/loader
    defaults.folderPath = '';
    defaults.fileType = '';
    defaults.available_metrics = {};
    defaults.metric_names = {};
    defaults.output_dir = '';
    defaults.selected_permutations = [];    
end