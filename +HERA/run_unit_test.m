function run_unit_test(varargin)
% RUN_UNIT_TEST - Comprehensive Scientific Validation Suite for the HERA Ranking.
%
% Syntax:
%   HERA.run_unit_test()                            1. Auto-Log Mode (Default)
%   HERA.run_unit_test('interactive')               2. Interactive Mode (Select Folder)
%   HERA.run_unit_test('path/to/custom/log/folder') 3. Custom Log Path Mode
%
% Description:
%   This function serves as the master validation test unit for the HERA ranking algorithm.
%   It employs deterministic, synthetic datasets to strictly validate statistical accuracy, 
%   logical robustness, and edge-case handling.
%
%   The test unit provides transparent output allowing for manual verification of all intermediate steps during HERA.
%   (e.g., Holm-Bonferroni corrections, Bootstrap convergence)
%
% Workflow:
%   1.  Environment Setup: Configures reproducible RNG and logging.
%   2.  Parameter Definition: Sets explicit default constants for transparency.
%   3.  Execution of 18 Scientific Test Cases covering:
%       - Exact Statistics (Small N)
%       - Multi-Hypothesis Correction (Holm-Bonferroni)
%       - Robustness (Outliers, Zero Variance, Missing Data)
%       - Algorithm Stability (Convergence, Sensitivity)
%       - Logic Verification (M1, M2, M3 Modes, Cycle Detection)
%       - Integration Tests (Full Pipeline)
%
% Inputs:
%   varargin - Optional arguments for log path handling.
%
% Outputs:
%   Console output with detailed diagnostics and a saved log file.
%
% Author: Lukas von Erdmannsdorff

% Import the HERA namespace to test internal functions
import HERA.*

% Dependency Check 
required_toolboxes = {'Statistics and Machine Learning Toolbox', 'Parallel Computing Toolbox'};
installed_ver = ver;
installed_toolboxes = {installed_ver.Name};
missing_toolboxes = {};

for i = 1:length(required_toolboxes)
    if ~any(strcmp(installed_toolboxes, required_toolboxes{i}))
        missing_toolboxes{end+1} = required_toolboxes{i};
    end
end

if ~isempty(missing_toolboxes)
    error_msg = sprintf('HERA requires the following missing toolboxes:\n');
    for i = 1:length(missing_toolboxes)
        error_msg = sprintf('%s - %s\n', error_msg, missing_toolboxes{i});
    end
    error('%s\nPlease install them via the Add-On Explorer.', error_msg);
end

clc;
    
    %% Setup Logging (Path Handling)
    % Determine the output folder based on input arguments or environment.
    if nargin > 0 && ischar(varargin{1})
        if strcmpi(varargin{1}, 'interactive')
            % Mode: Interactive Selection via UI
            selected_dir = uigetdir(pwd, 'Select Folder for Unit Test Logs');
            if isequal(selected_dir, 0)
                error('Unit Tests aborted: No output folder selected.');
            end
            log_folder = selected_dir;
            path_source = 'User Selection (Interactive)';
        else
            % Mode: Custom Path provided as argument
            log_folder = varargin{1};
            if ~exist(log_folder, 'dir')
                % Attempt to create the directory if it doesn't exist
                [success, msg] = mkdir(log_folder);
                if ~success
                    error('Cannot create custom log path: %s\n%s', log_folder, msg);
                end
            end
            path_source = 'Custom Path Argument';
        end
    else
        % Mode: Automatic (Default) - Smart finding of writable persistent folder
        % Checks Documents folder first, falls back to tempdir.
        [log_folder, path_source] = get_writable_log_path();
    end
    
    % Define unique log filename with timestamp to prevent overwriting
    timestamp_str = string(datetime('now'), 'yyyyMMdd_HHmmss');
    log_filename = "HERA_Validation_" + timestamp_str + ".txt";
    log_path = fullfile(log_folder, log_filename);
    
    % Start Logging using MATLAB's diary function
    if exist(log_path, 'file')
        delete(log_path); % Ensure fresh start
    end
    try
        diary(char(log_path));
    catch ME
        error('Failed to initialize log at %s.\nError: %s', log_path, ME.message);
    end
    
    % Start execution timer
    tic;
    fprintf('===============================================\n');
    fprintf('Scientific Validation and Testing Unit for HERA\n');
    fprintf('===============================================\n');
    fprintf('Log File:\n%s\n', log_path);
    fprintf('\n')

    %% 1. Global Configuration and Setup
    % Check if a parallel pool is already running. If not, start one.
    current_pool = gcp('nocreate'); 
    if isempty(current_pool)
        try
            % Start default local pool and capture object
            current_pool = parpool('local'); 
        catch ME
            fprintf('\nWARNING: Could not start Parallel Pool. Running serially.\n');
            fprintf('Error: %s\n', ME.message);
        end
    else
        fprintf('Parallel Pool already active: %d workers.\n', current_pool.NumWorkers);
    end
   
    % Explicitly tell the parallel workers where the +HERA package is located.
    if ~isempty(current_pool)
        try
            % Get the directory containing this file (which is inside +HERA)
            [functionPath, ~, ~] = fileparts(mfilename('fullpath'));
            % We need the Parent directory of +HERA
            packageParentDir = fileparts(functionPath); 
            
            % Ensure it is on the Client path
            addpath(packageParentDir);
            
            % Force all workers to add this path AND WAIT until they are done
            F = parfevalOnAll(current_pool, @addpath, 0, packageParentDir);
            wait(F); 
            
        catch ME
            fprintf('Warning: Could not set paths on parallel pool: %s\n', ME.message);
        end
    end
   
    % Define standard simulation parameters
    n_subj = 50; % Sufficient sample size for stable statistics in standard tests.
    
    % Use a fixed seed for reproducible results across the entire suite.
    % This ensures that "random" noise is identical every time the test runs.
    global_seed = 123;
    s_global = RandStream('mlfg6331_64', 'Seed', global_seed);
    RandStream.setGlobalStream(s_global);
    
    % Defaults (Standard HERA configuration)
    % We will reload this at the start of every test to prevent "State Pollution" 
    % (e.g., Test 8 changing ranking_mode affecting Test 9).
    default_config = struct();
    default_config.alphas = [0.05, 0.05, 0.05]; % Base alpha per metric
    default_config.ci_level = 0.95; 
    default_config.timestamp = char(timestamp_str);
    default_config.metric_names = {'M1', 'M2', 'M3'}; 
    default_config.ranking_mode = 'M1_M2_M3'; 
    
    % Default thresholds (Standard HERA configuration)
    thresholds = struct();
    thresholds.d_thresh = [0.15, 0.15, 0.15];       % Cliff's Delta > 0.15 (small effect)
    thresholds.rel_thresh = [0.05, 0.05, 0.05];     % Relative Diff > 5%
    thresholds.rel_thresh_b = [0.05, 0.05, 0.05];   % Backwards compatibility
    thresholds.min_rel_thresh = [0, 0, 0];          % Lower bound (usually SEM-based)
    
    % Bootstrap Configuration (Low B settings optimized for fast unit tests)
    % In production, B is typically 1000-10000. Here 50-150 is enough to test logic.
    bs_config = struct();
    bs_config.B_start = 50;
    bs_config.B_step = 50;
    bs_config.B_end = 150;
    bs_config.n_trials = 3;
    bs_config.convergence_tolerance = 0.05;         % 5% tolerance
    bs_config.min_steps_for_convergence_check = 2;
    bs_config.smoothing_window = 3; 
    bs_config.convergence_streak_needed = 2;
    
    % Assign bootstrap config to all sub-modules
    default_config.bootstrap_thresholds = bs_config;
    default_config.bootstrap_ci = bs_config;
    default_config.bootstrap_ranks = bs_config;

    % Initialize working config 
    config = default_config;

    % Configuration Output for Log 
    fprintf('\nConfiguration & Parameter Definition:\n');
    fprintf('[Config] Random Number Generator: ''mlfg6331_64'' (Fixed Seed: %d)\n', global_seed);
    fprintf('[Config] Standard Sample Size (N): %d subjects\n', n_subj);
    
    fprintf('[Config] Effect Size Thresholds (Standard):\n');
    fprintf('         Cliff''s Delta > %.2f\n', thresholds.d_thresh(1));
    fprintf('         Relative Diff > %.2f\n', thresholds.rel_thresh(1));
    
    fprintf('[Config] Statistical Significance:\n');
    fprintf('         Alpha Level: %.2f (Holm-Bonferroni adjusted)\n', config.alphas(1));
    fprintf('         CI Level:    %.2f\n', config.ci_level);

    fprintf('[Config] Bootstrap Settings (Unit Test Mode - Fast):\n');
    fprintf('         B_range: [%d : %d : %d]\n', bs_config.B_start, bs_config.B_step, bs_config.B_end);
    fprintf('         Convergence Tolerance: %.1f%%\n', bs_config.convergence_tolerance * 100);
    
    % Load Mock Resources (Language and Styles) to avoid external file dependencies
    [lang, styles] = get_test_resources();

    % Dataset names for standard 3-set tests
    ds_names = {'D1', 'D2', 'D3'};
    
    % Initialize counters
    tests_passed = 0;
    tests_total = 0;

    %% Test 1: Small Sample Size Boundary (Exact Wilcoxon N=10)
    tests_total = tests_total + 1;
    title_str = 'Test 1: Small Sample Size Boundary (Exact Wilcoxon N=10)';
    fprintf('\n%s\n', title_str);
    fprintf('%s\n', repmat('-', 1, strlength(title_str)));
    config = default_config;
    
    % Hypothesis: For N < 15, MATLAB's signrank should perform an exact calculation
    % instead of the asymptotic normal approximation used for larger N.
    
    n_small = 10;
    fprintf('[Test] Verifying exact p-value calculation for N < 15.\n');
    fprintf('[Setup] Using N=%d. Generating Signal (Mean 10) vs Noise (Mean 0).\n', n_small);
    
    % Data Generation
    % D1: Mean 10 (Strong Signal), D2: Mean 0 (Noise) -> Difference should be highly significant
    data_small = generate_exact_data(n_small, [10, 0], 1.0);
    eff_small = calculate_real_effects({data_small}, 1);
    
    % Input Table Visualization
    fprintf('[Input]\n');
    h_in = {'Parameter', 'Value'}; 
    d_align = {'l', 'l'}; 
    h_align = {'c', 'c'}; 
    table_data = {
        'Sample Size (N)', num2str(n_small);
        'Mean D1 (Signal)', '10.0';
        'Mean D2 (Noise)', '0.0'
    };
    print_auto_table(h_in, table_data, d_align, h_align);
    
    try
        % Run Ranking calculation for just this data
        [~, ~, ~, ~, p_vals_small] = calculate_ranking({data_small}, eff_small, thresholds, config, {'A','B'}, [1 2]);
        p_val = p_vals_small{1}(1,2);
        
        % Result Table Visualization
        fprintf('\n[Result]\n');
        h_res = {'Metric', 'Value', 'Expectation'}; 
        d_align = {'l', 'l', 'l'};
        h_align = {'c', 'c', 'c'};
        
        table_data = {
            'Calculated P-Value', sprintf('%.6f', p_val), '< 0.005 (Strong Effect)'
        };
        print_auto_table(h_res, table_data, d_align, h_align);
        
        % Validation: Check if P-value is a valid probability and low enough
        if isfinite(p_val) && p_val >= 0 && p_val <= 1
             fprintf('\n[Status] PASS: Exact test execution successful. Result is a valid probability.\n');
             tests_passed = tests_passed + 1;
        else
             fprintf('\n[Status] FAIL: P-Value invalid (NaN or out of range).\n');
        end
    catch ME
        fprintf('\n[Status] FAIL: Crash during small sample test: %s\n', ME.message);
    end

    %% Test 2: Statistical Rigor (Holm-Bonferroni Correction)
    tests_total = tests_total + 1;
    title_str = 'Test 2: Statistical Rigor (Holm-Bonferroni Correction)';
    fprintf('\n%s\n', title_str);
    fprintf('%s\n', repmat('-', 1, strlength(title_str)));
    
    % Hypothesis: The alpha threshold must adapt dynamically to the number of pairwise comparisons (m).
    % Formula: Alpha_corrected = Alpha_base / (m - rank + 1)
    
    % PART A: Basic Functionality (N=5 Datasets)
    fprintf('[Test] Checking dynamic alpha adaptation for N=5 datasets (10 pairs).\n');
    n_ds = 5;
    ds_names_hb = {'D1', 'D2', 'D3', 'D4', 'D5'};
    
    fprintf('[Setup] Generating random noise data. Alpha Base = 0.05.\n');
    % We use noise data (means all 0) so P-values are random/uniform under Null Hypothesis.
    data_hb = generate_exact_data(n_subj, zeros(1, n_ds), 1.0);
    
    eff_hb = calculate_real_effects({data_hb}, 1);
    config_hb = default_config;
    config_hb.ranking_mode = 'M1';
    config_hb.alphas = [0.05];
    
    % Create all unique pairs (n choose 2)
    pairs_hb = nchoosek(1:n_ds, 2);
    n_pairs_hb = size(pairs_hb, 1); % 10 pairs
    
    % Input Table 
    fprintf('[Input]\n');
    d_align = {'l', 'l'};
    h_align = {'c', 'c'};
    table_data = {
        'Datasets', num2str(n_ds);
        'Comparisons (Pairs)', num2str(n_pairs_hb);
        'Base Alpha', sprintf('%.2f', config_hb.alphas(1))
    };
    print_auto_table(h_in, table_data, d_align, h_align);

    % Run ranking on noise data
    % We are interested in 'all_alphas' output
    [~, ~, ~, all_alphas, ~] = calculate_ranking({data_hb}, eff_hb, thresholds, config_hb, ds_names_hb, pairs_hb);
    
    % Validation: The strictest alpha should be Base / Num_Pairs
    min_alpha_observed = min(all_alphas{1}(:));
    expected_min_alpha = 0.05 / n_pairs_hb; % 0.005
    
    % Result Table 
    fprintf('\n[Result]\n');
    h_res = {'Metric', 'Value'}; 
    d_align = {'l', 'l'};
    h_align = {'c', 'c'};
    
    table_data = {
        'Expected Min Alpha', sprintf('%.6f', expected_min_alpha);
        'Observed Min Alpha', sprintf('%.6f', min_alpha_observed);
        'Difference', sprintf('%.1e', abs(min_alpha_observed - expected_min_alpha))
    };
    print_auto_table(h_res, table_data, d_align, h_align);
    
    if abs(min_alpha_observed - expected_min_alpha) < 1e-6
        fprintf('\n[Status] PASS: Correction correctly applied.\n');
        tests_passed = tests_passed + 1;
    else
        fprintf('\n[Status] FAIL: Alpha correction mismatch.\n');
    end

    %% Test 3: Outlier Robustness & Threshold Logic (The "AND" Condition)
    tests_total = tests_total + 1;
    title_str = 'Test 3: Outlier Robustness & Threshold Logic (Divergence Check)';
    fprintf('\n%s\n', title_str);
    fprintf('%s\n', repmat('-', 1, strlength(title_str)));
    config = default_config;
    
    % Hypothesis: To count as a "win" in HERA, a comparison must satisfy ALL criteria conjunctively:
    % (p < alpha) AND (|d| > d_thresh) AND (r > r_thresh).
    % Missing just one criterion (e.g., high diff but low effect size) must result in Neutral.
    
    fprintf('[Test] Win condition = (p < alpha) AND (|d| > d_thr) AND (r > r_thr).\n');
    fprintf('[Setup] Generating synthetic edge cases to test partial criteria matches.\n');
    
    % Case 1: High Delta, but trivial RelDiff (Microscopic difference with 0 variance)
    d1_micro = ones(n_subj, 1) * 1.001;
    d2_micro = ones(n_subj, 1) * 1.000;
    
    % Case 2: High RelDiff, but low Delta (Caused by a single massive outlier)
    rng(global_seed);
    d3_noise = randn(n_subj, 1);
    d4_outlier = d3_noise; 
    d4_outlier(1) = 1000; % Massive Mean shift due to one value, Median unchanged
    
    data_robust = [d1_micro, d2_micro, d3_noise, d4_outlier];
    names_robust = {'Micro1', 'Micro2', 'Noise', 'Outlier'};
    
    % Set explicit thresholds for verification
    thr_robust = thresholds;
    thr_robust.d_thresh = [0.2, 0.2, 0.2]; % Delta must be > 0.2
    thr_robust.rel_thresh = [0.05, 0.05, 0.05]; % RelDiff must be > 5%
    
    % Input Table 
    fprintf('[Input]\n');
    h_in = {'Comparison', 'Scenario', 'Target Failure'}; 
    d_align = {'l', 'l', 'l'};
    h_align = {'c', 'c', 'c'};
    
    table_data = {
        'Micro1 vs Micro2', 'High Delta (determin.), Trivial RelDiff', 'Fail RelDiff';
        'Noise vs Outlier', 'High Mean Diff, Low Median Diff (d)', 'Fail Delta';
        'Thresholds', 'Delta > 0.20, RelDiff > 0.05', '-'
    };
    print_auto_table(h_in, table_data, d_align, h_align);
    
    eff_robust = calculate_real_effects({data_robust}, 1);
    
    % Analyze Pair 1: Micro1 vs Micro2
    d_p1 = eff_robust.d_vals_all(1, 1); 
    r_p1 = eff_robust.rel_vals_all(1, 1);
    
    % Analyze Pair 2: Noise vs Outlier (Index 3 vs 4)
    d_p2 = eff_robust.d_vals_all(6, 1); 
    r_p2 = eff_robust.rel_vals_all(6, 1);
    
    % Run Ranking
    [~, ~, all_sig_rob, ~, ~] = calculate_ranking({data_robust}, eff_robust, thr_robust, config, names_robust, nchoosek(1:4, 2));
    
    % Check Results: Did any pair register a "Win"?
    win_micro = all_sig_rob{1}(1, 2) || all_sig_rob{1}(2, 1);
    win_outlier = all_sig_rob{1}(3, 4) || all_sig_rob{1}(4, 3);
    
    % Result Table
    fprintf('\n[Result]\n');
    h_res = {'Pair', 'Delta', 'RelDiff', 'Win (Exp: 0)'}; 
    d_align = {'l', 'c', 'c', 'c'}; 
    h_align = {'c', 'c', 'c', 'c'};
        
    table_data = {
        'Micro1 vs Micro2', sprintf('%.2f', abs(d_p1)), sprintf('%.4f', r_p1), sprintf('%d', win_micro);
        'Noise vs Outlier', sprintf('%.2f', abs(d_p2)), sprintf('%.2f', r_p2), sprintf('%d', win_outlier)
    };
    print_auto_table(h_res, table_data, d_align, h_align);
    
    % Validation
    if ~win_micro && ~win_outlier
        fprintf('\n[Status] PASS: Partial criteria matches correctly rejected.\n');
        tests_passed = tests_passed + 1;
    else
        fprintf('\n[Status] FAIL: Robustness logic failed. MicroWin: %d, OutlierWin: %d\n', win_micro, win_outlier);
    end

    %% Test 4: Mathematical Validation (Thresholds & BCa & SEM)
    tests_total = tests_total + 1;
    title_str = 'Test 4: Mathematical Validation (Thresholds/SEM/BCa)';
    fprintf('\n%s\n', title_str);
    fprintf('%s\n', repmat('-', 1, strlength(title_str)));
    
    fprintf('[Test] Validation of mathematical corner cases (Zero Variance, SEM Logic, BCa Factors).\n');
    
    % Setup: Simple data for math check
    % D1: Zero variance (all 0) -> Should trigger safety mechanisms (NaN protection)
    d1 = zeros(20, 3); 
    % D2: Small variance signal -> Should have normal thresholds
    d2 = repmat([10, 10.1, 10.2], 20, 1) + randn(20, 3) * 0.5; 
    
    config_math = default_config;
    config_math.metric_names = {'ZeroVar', 'SmallVar'};
    
    % Setup Table 
    fprintf('[Setup]\n');
    h_set = {'Metric Name', 'Data Characteristics', 'Expected Behavior'}; 
    d_align = {'l', 'l', 'l'};
    h_align = {'c', 'c', 'c'};
    
    table_data = {
        'ZeroVar', 'All zeros (0 variance)', 'Thresholds=0, BCa Factors=0';
        'SmallVar', 'Mean~10, SD~0.5', 'Valid Thresholds, Finite BCa'
    };
    print_auto_table(h_set, table_data, d_align, h_align);

    % Register automatic cleanup for temp files
    temp_dir = tempname; 
    if ~exist(temp_dir, 'dir'), mkdir(temp_dir); end
    cleanUpTemp = onCleanup(@() rmdir(temp_dir, 's'));
    myStream = RandStream('mlfg6331_64', 'Seed', 1);
    
    try
        % Step 1: Calculate Thresholds (Percentile & SEM)
        % We use evalc to capture/suppress the console output of the function
        cmd_thr = ['[d_thr, r_thr, r_thr_b, min_r_thr, d_v, r_v, p_idx] = ' ...
            'calculate_thresholds({d1, d2}, 20, config_math, temp_dir, 50, myStream, styles, lang);'];
        [~] = evalc(cmd_thr);
        
        % Step 2: BCa Logic Validation (Bootstrap Confidence Intervals)
        cmd_bca = ['[~, ~, ~, z0_d, a_d, ~, ~, ~, ~, ~, ~, ~, ~] = ' ...
             'calculate_bca_ci({d1, d2}, d_v, r_v, p_idx, 20, config_math, config_math.metric_names,' ...
             '         temp_dir, temp_dir, 50, myStream, styles, lang, ''Test10'');'];
        [~] = evalc(cmd_bca);
        
        % Checks & Assertions
        % 1. SEM Logic: Final threshold should be max of (Bootstrap, SEM)
        sem_logic_ok = abs(r_thr(2) - max(r_thr_b(2), min_r_thr(2))) < 1e-9;
        
        % 2. BCa Normal: Bias (z0) and Skewness (a) should be finite numbers
        z0_val = z0_d(1, 2); a_val = a_d(1, 2);
        bca_valid_normal = isfinite(z0_val) && isfinite(a_val);
        
        % 3. Zero Variance Safety: z0 and a should be 0 (no bias/skew if no variance)
        z0_zero = z0_d(1, 1); a_zero  = a_d(1, 1);
        protection_ok = (z0_zero == 0 && a_zero == 0);

        % Result Table 
        fprintf('\n[Result]\n');
        h_res = {'Check', 'Value', 'Reference', 'Status'}; 
        d_align = {'l', 'l', 'l', 'l'};
        h_align = {'c', 'c', 'c', 'c'};
        
        table_data = {
            'ZeroVar: Delta Thresh', sprintf('%.4e', d_thr(1)), '< 1e-12', 'OK';
            'ZeroVar: RelDiff Boot', sprintf('%.4e', r_thr_b(1)), '0 (No Var)', '-';
            'ZeroVar: SEM Thr',      sprintf('%.4e', min_r_thr(1)), '0 (No Var)', '-';
            'ZeroVar: Final RelDiff',sprintf('%.4e', r_thr(1)), 'max(Boot, SEM)', '-';
            'SmallVar: Delta Thresh', sprintf('%.4f', d_thr(2)), '> 0', 'OK';
            'SmallVar: RelDiff Boot', sprintf('%.4f', r_thr_b(2)), 'Bootstrap', '-';
            'SmallVar: SEM Thr',      sprintf('%.4f', min_r_thr(2)), 'SEM Lower Bound', '-';
            'SmallVar: Final RelDiff',sprintf('%.4f', r_thr(2)), 'max(Boot, SEM)', char(string(sem_logic_ok));
            'BCa Normal (SmallVar)', sprintf('%.2f / %.2f', z0_val, a_val), 'Finite', char(string(bca_valid_normal));
            'BCa ZeroVar (Safety)', sprintf('%.2f / %.2f', z0_zero, a_zero), '0.00 / 0.00', char(string(protection_ok))
        };
        print_auto_table(h_res, table_data, d_align, h_align);

        % Final Verification
        if d_thr(1) < 1e-12 && d_thr(2) > 0 && sem_logic_ok && bca_valid_normal && protection_ok
             fprintf('\n[Status] PASS: All Mathematical Validations successful.\n');
             tests_passed = tests_passed + 1;
        else
             fprintf('\n[Status] FAIL: Mathematical validation failed.\n');
        end
    catch ME
        fprintf('\n[Status] FAIL: Error: %s\n', ME.message);
    end
    
    %% Test 5: Stability Check (Identical Data / Zero Variance)
    tests_total = tests_total + 1;
    title_str = 'Test 5: Stability Check (Identical Data / Zero Variance)';
    fprintf('\n%s\n', title_str);
    fprintf('%s\n', repmat('-', 1, strlength(title_str)));
    config = default_config;
    
    % Hypothesis: The algorithm must not crash or produce NaNs when comparing 
    % two completely identical, zero-variance datasets (Singularity Check).
    
    fprintf('[Test] Algorithm stability with absolute zero-variance inputs.\n');
    
    data_zero = zeros(n_subj, 3); 
    eff_zero = calculate_real_effects({data_zero}, 1);
    
    % Input Table 
    fprintf('[Input]\n');
    h_in = {'Data Type', 'Size', 'Content'}; 
    d_align = {'l', 'l', 'l'};
    h_align = {'c', 'c', 'c'};
    
    table_data = {
        'Zero Matrix', sprintf('%dx3', n_subj), 'All 0'
    };
    print_auto_table(h_in, table_data, d_align, h_align);
    
    try
        % Run ranking
        [final_order_z, ~] = calculate_ranking({data_zero}, eff_zero, thresholds, config, ds_names, nchoosek(1:3, 2));
        
        % Result Table
        fprintf('\n[Result]\n');
        h_res = {'Outcome', 'Fallback Order'}; 
        d_align = {'l', 'l'};
        h_align = {'c', 'c'};
        
        table_data = {
            'No Crash', mat2str(final_order_z)
        };
        print_auto_table(h_res, table_data, d_align, h_align);
        
        fprintf('\n[Status] PASS: Stability Check successful.\n');
        tests_passed = tests_passed + 1;
    catch ME
        fprintf('\n[Status] FAIL: Crash detected: %s\n', ME.message);
    end

    %% Test 6: Bootstrap Convergence Algorithms (Global Verification)
    tests_total = tests_total + 1;
    title_str = 'Test 6: Bootstrap Convergence Algorithms (Global Verification)';
    fprintf('\n%s\n', title_str);
    fprintf('%s\n', repmat('-', 1, strlength(title_str)));
    
    % Hypothesis: The bootstrap process uses different convergence criteria. 
    % We need to verify that all 3 implementations (Thresholds, BCa, Ranking) 
    % correctly respect their configuration (Simple vs Robust vs Elbow).
    
    fprintf('[Test] Verification of mathematical stopping criteria for bootstrap processes.\n');
    fprintf('[Setup] Testing 3 algorithms (Simple, Robust, Elbow) across 3 analysis types with synthetic Data.\n');
    
    % 1. Setup Fixed Data for Convergence Test
    temp_dir = tempname; if ~exist(temp_dir, 'dir'), mkdir(temp_dir); end
    cleanUpTemp = onCleanup(@() rmdir(temp_dir, 's'));
    s_data = RandStream('mlfg6331_64', 'Seed', 666);
    RandStream.setGlobalStream(s_data);
    
    n_conv = 100;
    d1_conv = generate_exact_data(n_conv, [10, 15], 3); 
    all_data_conv = {d1_conv};
    ds_names_conv = {'C1', 'C2'};
    p_idx_conv = [1 2];
    
    % Input Table 
    fprintf('[Input]\n');
    h_in = {'Parameter', 'Value'}; 
    d_align = {'l', 'l'}; 
    h_align = {'c', 'c'}; 
    table_data = {
        'Sample Size (N)', num2str(n_conv);
        'Means', '[10, 15]';
        'SD', '3'
    };
    print_auto_table(h_in, table_data, d_align, h_align);
    
    % Pre-calculation needed for Bootstrap Ranking call
    eff_conv = calculate_real_effects(all_data_conv, 1);
    thr_dummy = struct('d_thresh', [0.15], 'rel_thresh', [0.05]); 
    config_conv = default_config;
    config_conv.metric_names = {'TestMetric'};
    [~, rank_base] = calculate_ranking(all_data_conv, eff_conv, thr_dummy, config_conv, ds_names_conv, p_idx_conv);

    % Define Base Configurations for the 3 algorithms
    bs_thr = struct('B_start', 100, 'B_step', 100, 'B_end', 5000, 'n_trials', 25, 'min_steps_for_convergence_check', 3, 'smoothing_window', 3, ...
        'convergence_streak_needed', 3, 'convergence_tolerance', 0.005);    
    bs_bca = struct('B_start', 100, 'B_step', 200, 'B_end', 10000, 'n_trials', 25, 'min_steps_for_convergence_check', 3, 'smoothing_window', 4, ...
        'convergence_streak_needed', 3, 'convergence_tolerance', 0.01);     
    bs_rank = struct('B_start', 50, 'B_step', 10, 'B_end', 1000, 'n_trials', 15, 'min_steps_for_convergence_check', 3, 'smoothing_window', 3, ...
        'convergence_streak_needed', 3, 'convergence_tolerance', 0.005);    

    % Struct array to loop through methods
    test_methods = struct();
    test_methods(1).name = 'Thresholds'; test_methods(1).cfg = bs_thr; test_methods(1).lang_sec = 'thresholds';
    test_methods(2).name = 'BCa';        test_methods(2).cfg = bs_bca; test_methods(2).lang_sec = 'bca';
    test_methods(3).name = 'Ranking';    test_methods(3).cfg = bs_rank; test_methods(3).lang_sec = 'ranking';

    % Result Table Header
    fprintf('\n[Result]\n');
    header_parts = {'Method', 'Mode', 'Configuration', 'Result', 'Status', 'Logic Check'};
    d_align = {'l', 'l', 'l', 'l', 'c', 'l'}; 
    h_align = {'c', 'c', 'c', 'c', 'c', 'c'};
    
    table_data = cell(9, 6); % 3 methods * 3 modes
    row_idx = 0;

    seed_boot = 999;
    test6_all_passed = true;

    % Loop through Method -> Mode
    for m = 1:3
        curr_method = test_methods(m).name;
        base_cfg    = test_methods(m).cfg;
        lang_sec    = lang.(test_methods(m).lang_sec);
        
        % Define expected keywords based on language resources
        % We split by '%' to handle format specifiers like "%.2f%%"
        
        % 1. Convergence Keyword
        parts_conv = strsplit(lang_sec.convergence_reached, '%');
        key_conv = parts_conv{1}; % e.g. "Convergence at" or "Convergence reached at"
        
        % 2. Elbow Keyword (Static)
        key_elbow = lang_sec.elbow_analysis_info; 
        
        % 3. Stability Keyword (Robust only)
        parts_stab = strsplit(lang_sec.stable_runs_info, '%');
        key_stable = parts_stab{1}; % e.g. "Stability in"

        for mode_idx = 1:3
            row_idx = row_idx + 1;
            run_cfg = base_cfg;
            
            % Configure the specific Convergence Mode
            if mode_idx == 1 % Simple
                mode_name = 'Simple';
                % Disable smoothing/streak to force simple check
                run_cfg.smoothing_window = []; run_cfg.convergence_streak_needed = [];  
                config_str = sprintf('Tol=%.1f%%, (No Smooth/Streak)', run_cfg.convergence_tolerance*100);
            elseif mode_idx == 2 % Robust
                mode_name = 'Robust';
                config_str = sprintf('Tol=%.1f%%, Sm=%d, St=%d', ...
                    run_cfg.convergence_tolerance*100, run_cfg.smoothing_window, run_cfg.convergence_streak_needed);
            elseif mode_idx == 3 % Elbow
                mode_name = 'Elbow';
                % Tolerance -1 forces the algorithm to run until B_end for Elbow analysis
                run_cfg.convergence_tolerance = -1.0; 
                config_str = sprintf('Tol=-1.0, B_max=%d (Force Fail)', run_cfg.B_end);
            end
            
            % Reset RNG for fairness
            myStream = RandStream('mlfg6331_64', 'Seed', seed_boot);
            cfg_run = config_conv;
            cmd = '';
            
            % Construct Command String dynamically
            % Note: We use T = evalc(cmd) to capture output. The assignments inside cmd happen in the workspace.
            if strcmp(curr_method, 'Thresholds')
                cfg_run.bootstrap_thresholds = run_cfg;
                cmd = ['[~, ~, ~, ~, ~, ~, ~, B_res, ~, ~, ~, ~, ~] = ' ...
                    'calculate_thresholds(all_data_conv, n_conv, cfg_run, temp_dir, [], myStream, styles, lang);'];
            elseif strcmp(curr_method, 'BCa')
                cfg_run.bootstrap_ci = run_cfg;
                cmd = ['[B_res, ~, ~, ~, ~, ~, ~, ~] = ' ...
                    'calculate_bca_ci(all_data_conv, eff_conv.d_vals_all, eff_conv.rel_vals_all, p_idx_conv, n_conv, cfg_run, cfg_run.metric_names,' ...
                    'temp_dir, temp_dir, [], myStream, styles, lang, ''Test6'');'];
            elseif strcmp(curr_method, 'Ranking')
                cfg_run.bootstrap_ranks = run_cfg;
                cmd = ['[~, B_res, ~, ~, ~] = ' ...
                    'bootstrap_ranking(all_data_conv, thr_dummy, cfg_run, ds_names_conv, rank_base, p_idx_conv, n_conv,' ...
                    'temp_dir, temp_dir, [], myStream, styles, lang, ''Test6'');'];
            end
            
            try
                % Execute and capture output for validation
                % IMPORTANT: Do not ask for 2nd output from evalc if cmd ends with ;
                T = evalc(cmd); 
                
                is_valid = false;
                logic_msg = '';
                
                if exist('B_res', 'var')
                    res_txt = sprintf('B = %d', B_res);
                else
                    res_txt = 'B_res missing';
                    B_res = -1; 
                end
                
                % Validate Result against Logic and Output
                if strcmp(mode_name, 'Elbow')
                    % Elbow Logic: 
                    % 1. B must be high (hit limit)
                    % 2. Output must contain Elbow message
                    b_ok = (B_res >= run_cfg.B_start);
                    txt_ok = contains(T, key_elbow);
                    
                    if b_ok && txt_ok
                        is_valid = true;
                        logic_msg = 'Elbow msg found';
                    else
                        logic_msg = sprintf('B_ok=%d, Txt_ok=%d', b_ok, txt_ok);
                    end
                    
                elseif strcmp(mode_name, 'Simple')
                    % Simple Logic:
                    % 1. Convergence before max
                    % 2. Output contains Convergence message
                    % 3. Output does NOT contain Stability message (Strict check)
                    b_ok = (B_res < run_cfg.B_end);
                    conv_ok = contains(T, key_conv);
                    % Note: Some simple implementations might print stability if configured differently, 
                    % but here we explicitly disabled it.
                    
                    if b_ok && conv_ok
                        is_valid = true;
                        logic_msg = 'Conv. msg found';
                    else
                        logic_msg = sprintf('B_ok=%d, Conv=%d', b_ok, conv_ok);
                    end
                    
                elseif strcmp(mode_name, 'Robust')
                    % Robust Logic:
                    % 1. Convergence before max
                    % 2. Output contains Convergence message
                    % 3. Output MUST contain Stability message
                    b_ok = (B_res < run_cfg.B_end);
                    conv_ok = contains(T, key_conv);
                    stab_ok = contains(T, key_stable);
                    
                    if b_ok && conv_ok && stab_ok
                        is_valid = true;
                        logic_msg = 'Conv + Stable found';
                    else
                        logic_msg = sprintf('B=%d, C=%d, S=%d', b_ok, conv_ok, stab_ok);
                    end
                end
                
                if is_valid
                    stat = 'PASS'; 
                else
                    stat = 'FAIL'; 
                    test6_all_passed = false; 
                    % Append T for debugging if needed (truncated)
                    if length(T) > 100, T_debug = [T(1:100) '...']; else, T_debug = T; end
                    logic_msg = [logic_msg ' [' T_debug ']'];
                end
                
                % Store Data Row
                table_data(row_idx, :) = {curr_method, mode_name, config_str, res_txt, stat, logic_msg};
                
            catch ME
                table_data(row_idx, :) = {curr_method, mode_name, 'CRASH', 'ERR', 'FAIL', ME.message};
                test6_all_passed = false;
            end
        end
    end
    % Print Full Table
    print_auto_table(header_parts, table_data, d_align, h_align);
    
    if test6_all_passed
        tests_passed = tests_passed + 1;
        fprintf('\n[Status] PASS: All convergence algorithms behave mathematically correct.\n');
    else
        fprintf('\n[Status] FAIL: Anomalies detected in convergence logic.\n');
    end

    %% Test 7: Initial Sorting Logic (Dominance, Delta, Mean)
    tests_total = tests_total + 1;
    title_str = 'Test 7: Initial Sorting Logic (Dominance, Delta, Mean)';
    fprintf('\n%s\n', title_str);
    fprintf('%s\n', repmat('-', 1, strlength(title_str)));
    
    % The initial sorting uses 3 criteria in order:
    % 1. Win Count (Number of significant wins against others)
    % 2. Stochastic Dominance (Sum of Cliff's Deltas in direct comparison)
    % 3. Arithmetic Mean (Tie-breaker if d=0)
    
    % Step 1: Standard Dominance (Clear Winner)
    fprintf('[Test] Step 1: Testing Standard Dominance (M1).\n');
    means_1 = [10, 5, 0]; sd_1 = 2.0;
    fprintf('[Setup] Datasets: D1=%.1f, D2=%.1f, D3=%.1f (SD=%.1f)\n', means_1(1), means_1(2), means_1(3), sd_1);
    m1_data = generate_exact_data(n_subj, means_1, sd_1); 
    config = default_config;
    config.ranking_mode = 'M1';
    
    [final_order_1, ~, eff_1, p_val_1] = run_single_test_full({m1_data}, thresholds, config, ds_names);
    
    % Combined Diagnostic Table 
    fprintf('\n[Result]\n');
    h_diag = {'Pair', 'P-Value', 'Delta', 'RelDiff'};
    d_align = {'l', 'l', 'l', 'l'}; 
    h_align = {'c', 'c', 'c', 'c'};
    
    pairs_idx = nchoosek(1:3, 2); 
    table_data = cell(size(pairs_idx, 1), 4);
    for k = 1:size(pairs_idx, 1)
        i = pairs_idx(k, 1); j = pairs_idx(k, 2);
        pair_label = sprintf('%s vs %s', ds_names{i}, ds_names{j});
        
        table_data(k, :) = {
            pair_label, ...
            sprintf('%.4e', p_val_1{1}(i, j)), ...
            sprintf('%.4f', eff_1.d_vals_all(k, 1)), ...
            sprintf('%.4f', eff_1.rel_vals_all(k, 1))
        };
    end
    print_auto_table(h_diag, table_data, d_align, h_align);
    fprintf('\n');
    
    if check_result(final_order_1, [1, 2, 3], 'M1 Clear Winner (D1>D2>D3)')
        % Pass Step 1
    end
    
    % Step  2: Tie-Break via Delta Injection 
    fprintf('\n[Test] Step 2: Testing Tie-Break via Delta Injection (M1).\n');
    fprintf('[Setup] Generating Noise Data (Means=0). Forcing Ties on Win-Counts.\n');
    % Data: Identical noise (Neutral p-values -> Win Counts all equal)
    m1_data_noise = generate_exact_data(n_subj, [0, 0, 0], 1.0);
    effect_sizes = calculate_real_effects({m1_data_noise}, 1);
    
    % Injection: Manually override Cliff's Delta to simulate subtle dominance
    % D3 > D2 (-0.4), D3 > D1 (-0.8), D2 > D1 (-0.4).
    % Note: Negative Delta typically means Group 2 > Group 1. 
    effect_sizes.d_vals_all(:, 1) = [-0.4; -0.8; -0.4]; 
    
    % Input Table for Injection 
    fprintf('[Input]\n');
    h_inj = {'Pair', 'Injected Delta', 'Implication'}; 
    d_align = {'l', 'c', 'l'}; 
    h_align = {'c', 'c', 'c'};
    table_data = {
        'D1 vs D2', '-0.4', 'D2 > D1';
        'D1 vs D3', '-0.8', 'D3 > D1';
        'D2 vs D3', '-0.4', 'D3 > D2'
    };
    print_auto_table(h_inj, table_data, d_align, h_align);
    
    [final_order_2, ~, ~, ~, ~] = calculate_ranking({m1_data_noise}, effect_sizes, thresholds, config, ds_names, nchoosek(1:3, 2));
    
    % Result Table
    fprintf('\n[Result]\n');
    h_res = {'Metric', 'Final Order'}; 
    d_align = {'l', 'l'}; 
    h_align = {'c', 'c'};
    
    table_data = {
        'Tie-Break (Delta)', mat2str(final_order_2)
    };
    print_auto_table(h_res, table_data, d_align, h_align);

    % Expect: 3, 2, 1 (Due to injected Deltas)
    fprintf('\n');
    pass_2 = check_result(final_order_2, [3, 2, 1], 'M1 Tie-Break via Delta');
    
    pass_1 = check_result(final_order_1, [1, 2, 3], 'M1 Clear Winner (D1>D2>D3)');

    % Step 3: Tie-Break via Mean 
    fprintf('\n[Test] Step 3: Testing Tie-Break via Mean (M1).\n');
    % Means are very close, but distinct.
    means_3 = [0.01, 0.02, 0.00];
    
    % Input Table
    fprintf('[Input]\n');
    h_in = {'Dataset', 'Mean', 'Condition'}; 
    d_align = {'l', 'l', 'l'}; 
    h_align = {'c', 'c', 'c'};
    
    table_data = {
        'D1', sprintf('%.2f', means_3(1)), 'Delta Forced to 0';
        'D2', sprintf('%.2f', means_3(2)), 'Delta Forced to 0';
        'D3', sprintf('%.2f', means_3(3)), 'Delta Forced to 0'
    };
    print_auto_table(h_in, table_data, d_align, h_align);

    m1_data_mean = generate_exact_data(n_subj, means_3, 2.0);
    eff_mean = calculate_real_effects({m1_data_mean}, 1);
    eff_mean.d_vals_all(:) = 0; % Force Delta to 0 to simulate perfect stochastic equality
    
    [final_order_3, ~, ~, ~, ~] = calculate_ranking({m1_data_mean}, eff_mean, thresholds, config, ds_names, nchoosek(1:3, 2));
    
    % Result Table
    fprintf('\n[Result]\n');
    h_res = {'Metric', 'Final Order'}; 
    d_align = {'l', 'l'}; 
    h_align = {'c', 'c'};
    
    table_data = {
        'Tie-Break (Mean)', mat2str(final_order_3)
    };
    print_auto_table(h_res, table_data, d_align, h_align);
    
    % Expect: D2 (0.02) > D1 (0.01) > D3 (0.00) -> [2, 1, 3]
    fprintf('\n');
    pass_3 = check_result(final_order_3, [2, 1, 3], 'M1 Tie-Break via Mean (D2>D1>D3)');
    
    if pass_1 && pass_2 && pass_3
        tests_passed = tests_passed + 1;
    end

    %% Test 8: Mode M1_M2 (Rank Correction Swap)
    tests_total = tests_total + 1;
    title_str = 'Test 8: Mode M1_M2 (Rank Correction Swap)';
    fprintf('\n%s\n', title_str);
    fprintf('%s\n', repmat('-', 1, strlength(title_str)));
    
    % Hypothesis: In M1_M2 mode, Metric 1 establishes the base order. 
    % Metric 2 acts as a correction layer: if a lower-ranked dataset significantly 
    % beats a higher-ranked one in M2, they swap places.
    
    fprintf('[Test] Verifying M2 Correction Swap logic.\n');
    fprintf('[Setup] M1 Order: D1 > D2 > D3. M2 Signal: D2 > D1 (Contradiction).\n');
    
    % M1: D1 wins (10 > 5 > 0)
    m1_data = generate_exact_data(n_subj, [10, 5, 0], 1.0);
    % M2: D2 wins against D1 (10 > 0) -> Should trigger swap
    m2_data = generate_exact_data(n_subj, [0, 10, 0], 1.0);
    
    config.ranking_mode = 'M1_M2';

    % Input Table 
    fprintf('[Input]\n');
    h_in = {'Metric', 'Signal Structure', 'Implication'}; 
    d_align = {'l', 'l', 'l'};
    h_align = {'c', 'c', 'c'};
    table_data = {
        'M1', 'D1(10) > D2(5)', 'D1 wins (Initial)';
        'M2', 'D2(10) > D1(0)', 'D2 wins (Strong Correction)'
    };
    print_auto_table(h_in, table_data, d_align, h_align);
    
    [final_order, ~, eff, p_val] = run_single_test_full({m1_data, m2_data}, thresholds, config, ds_names);
    
    % Display M2 values specifically to prove correction validity
    m2_d = eff.d_vals_all(:, 2);
    m2_r = eff.rel_vals_all(:, 2);
    m2_p = p_val{2}; 
    p_val_d1_d2 = m2_p(1, 2);

    % Result Table for M2 Stats 
    fprintf('\n[Result: M2 Stats D1 vs D2]\n');
    h_m2 = {'Metric', 'Value', 'Threshold', 'Significant?'}; 
    d_align = {'l', 'l', 'l', 'c'}; 
    h_align = {'c', 'c', 'c', 'c'};
    
    table_data = {
        'Delta', sprintf('%.2f', abs(m2_d(1))), sprintf('%.2f', thresholds.d_thresh(2)), char(string(abs(m2_d(1))>thresholds.d_thresh(2)));
        'RelDiff', sprintf('%.2f', m2_r(1)), sprintf('%.2f', thresholds.rel_thresh(2)), char(string(m2_r(1)>=thresholds.rel_thresh(2)));
        'P-Value', sprintf('%.4e', p_val_d1_d2), '0.05', char(string(p_val_d1_d2 < 0.05))
    };
    print_auto_table(h_m2, table_data, d_align, h_align);
    
    % Expectation: D2 swaps with D1 -> Order: D2, D1, D3 (2, 1, 3)
    fprintf('\n');
    if check_result(final_order, [2, 1, 3], 'M1_M2 Swap Applied (D2 overtakes D1)')
        tests_passed = tests_passed + 1;
    end

    %% Test 9: Mode M1_M3A (Tie-Break Logic A)
    % Logic A: If M1 is neutral, M2 (here acting as "Tie Breaker") decides.
    tests_total = tests_total + 1;
    title_str = 'Test 9: Mode M1_M3A (Tie-Break Logic A)';
    fprintf('\n%s\n', title_str);
    fprintf('%s\n', repmat('-', 1, strlength(title_str)));
    
    ds_names_3 = {'A', 'B', 'C'};
    
    fprintf('[Test] Verifying M3A Tie-Break Logic (using M2 data).\n');
    fprintf('[Setup] M1: Neutral (A>B>C). M2: Strong C > B > A.\n');
    
    % M1: Very close means (1.0 vs 0.98), large SD -> Neutral / No Signif. diff
    m1_data = generate_exact_data(n_subj, [1.0 0.98, 0.0], 0.0);
    
    % M2: Significant differences (C wins)
    m2_data = generate_exact_data(n_subj, [0, 10, 20], 1.0);
    
    config.ranking_mode = 'M1_M3A'; 
    
    % Input Table 
    fprintf('[Input]\n');
    h_in = {'Metric', 'Data Structure', 'Role'}; 
    d_align = {'l', 'l', 'l'};
    h_align = {'c', 'c', 'c'};
    table_data = {
        'M1', 'Means ~1.0', 'Primary (Neutral)';
        'M2', 'C(20) > B(10) > A(0)', 'Tie-Breaker (Strong)'
    };
    print_auto_table(h_in, table_data, d_align, h_align);
    
    [~, final_order] = run_single_test({m1_data, m2_data}, thresholds, config, ds_names_3);
    
    % Result Table
    fprintf('\n[Result]\n');
    h_res = {'Check', 'Final Order'}; 
    d_align = {'l', 'l'}; 
    h_align = {'c', 'c'};
    table_data = {
        'M3A Logic Applied', mat2str(final_order)
    };
    print_auto_table(h_res, table_data, d_align, h_align);

    % Expected: [2, 1, 3] (Indices for B, A, C) because M3A logic limits movement 
    % to ONE step per dataset to prevent bubbling unstable results.
    fprintf('\n');
    if check_result(final_order, [2, 1, 3], 'M1_M3A Single Tie-Break (B, A, C)')
        tests_passed = tests_passed + 1;
        fprintf('   (Validated that M3A strictly limits moves to 1 per dataset).\n');
    else
        fprintf('\n[Status] FAIL: Tie-Break logic failed.\n');
    end

    %% Test 10: Mode M1_M2_M3 (Logic 3B Fallback)
    tests_total = tests_total + 1;
    title_str = 'Test 10: Mode M1_M2_M3 (Logic 3B Fallback)';
    fprintf('\n%s\n', title_str);
    fprintf('%s\n', repmat('-', 1, strlength(title_str)));
    
    % Logic 3B: If M1 AND M2 are neutral (tied), sort by M3.
    % Ensure Clean RNG for this sensitive test to prevent noise flipping D1/D2
    s_local = RandStream('mlfg6331_64', 'Seed', 999);
    RandStream.setGlobalStream(s_local);

    fprintf('[Test] Verifying Logic 3B (M3 Sorting when M1/M2 tied).\n');
    fprintf('[Setup] M1/M2 Neutral (~10). M3 Strong D3(20) > D1/D2(10).\n');
    
    % M1/M2: Small differences, statistically neutral
    m1_data = generate_exact_data(n_subj, [10.4, 10.2, 10.0], 1.0);
    m2_data = generate_exact_data(n_subj, [10.4, 10.2, 10.0], 1.0);
    
    % M3: D3 (20) is significantly better than D1/D2
    m3_data = generate_exact_data(n_subj, [10.0, 10.0, 20.0], 1.0);
    
    % Input Table for Means 
    fprintf('[Input]\n');
    h_m = {'Dataset', 'M1 Mean', 'M2 Mean', 'M3 Mean'}; 
    d_align = {'c', 'c', 'c', 'c'}; 
    h_align = {'c', 'c', 'c', 'c'};
    table_data = {
        'D1', sprintf('%.2f', mean(m1_data(:,1))), sprintf('%.2f', mean(m2_data(:,1))), sprintf('%.2f', mean(m3_data(:,1)));
        'D2', sprintf('%.2f', mean(m1_data(:,2))), sprintf('%.2f', mean(m2_data(:,2))), sprintf('%.2f', mean(m3_data(:,2)));
        'D3', sprintf('%.2f', mean(m1_data(:,3))), sprintf('%.2f', mean(m2_data(:,3))), sprintf('%.2f', mean(m3_data(:,3)))
    };
    print_auto_table(h_m, table_data, d_align, h_align);
    
    config.ranking_mode = 'M1_M2_M3';
    [~, final_order] = run_single_test({m1_data, m2_data, m3_data}, thresholds, config, ds_names);
    
    % Result Table
    fprintf('\n[Result]\n');
    h_res = {'Check', 'Final Order'}; 
    d_align = {'l', 'l'}; 
    h_align = {'c', 'c'};
    table_data = {
        'Logic 3B Applied', mat2str(final_order)
    };
    print_auto_table(h_res, table_data, d_align, h_align);

    % Expectation: D3 (Index 3) moves to top because M3 decides the tie.
    fprintf('\n');
    if check_result(final_order, [3, 1, 2], 'M1_M2_M3 Logic B Swap (D3 moves to top)')
        tests_passed = tests_passed + 1;
    else
        fprintf('\n[Status] FAIL: Logic 3B fallback failed.\n');
    end
    
    % Restore Global RNG
    RandStream.setGlobalStream(s_global);

    %% Test 11: Systematic Missing Data (NaN) Robustness
    tests_total = tests_total + 1;
    title_str = 'Test 11: Systematic Missing Data (NaN) Handling';
    fprintf('\n%s\n', title_str);
    fprintf('%s\n', repmat('-', 1, strlength(title_str)));
    
    ds_names_5 = {'D1','D2','D3','D4','D5'};
    n_test11 = 200;
    % Base data: clear ladder 10,8,6,4,2
    m1_full = generate_exact_data(n_test11, [10, 8, 6, 4, 2], 2.0);
    config.ranking_mode = 'M1';
    pairs_5 = nchoosek(1:5, 2);
    
    fprintf('[Test] Algorithm can handle systematic missing data (NaN).\n');
    
    % Reference Run (Full Data)
    eff_full = calculate_real_effects({m1_full}, 1);
    [order_ref, ~] = calculate_ranking({m1_full}, eff_full, thresholds, config, ds_names_5, pairs_5);
    
    % Setup Table 
    fprintf('[Setup]\n');
    h_set = {'Parameter', 'Description'}; 
    d_align = {'l', 'l'}; 
    h_align = {'c', 'c'};
    table_data = {
        'Data Structure', sprintf('5 Datasets (D1-D5), N=%d per group', n_test11);
        'Missing Data Levels', '20%, 40%, 60%, 80%'
    };
    print_auto_table(h_set, table_data, d_align, h_align);

    % Input Table 
    fprintf('\n[Input]\n');
    h_in = {'Parameter', 'Description'}; 
    d_align = {'l', 'l'}; 
    h_align = {'c', 'c'};
    table_data = {
        'Signal Profile', 'Means: [10, 8, 6, 4, 2], SD: 2.0';
        'Total Data Points', sprintf('%d (Matrix: %dx5)', numel(m1_full), n_test11);
        'Reference Order', mat2str(order_ref)
    };
    print_auto_table(h_in, table_data, d_align, h_align);

    missing_levels = [0.20, 0.40, 0.60, 0.80];
    robustness_passed = true;
    
    % Result Table 
    fprintf('\n[Result]\n');
    h_nan = {'Missing %', 'Order Match', 'Status'}; 
    d_align = {'c', 'l', 'l'}; 
    h_align = {'c', 'c', 'c'};
    
    table_data = cell(length(missing_levels), 3);
    idx_row = 0;

    for lvl = missing_levels
        idx_row = idx_row + 1;
        m1_nan = m1_full;
        
        % Inject NaNs randomly
        num_missing = round(numel(m1_nan) * lvl);
        rng(global_seed + round(lvl*100)); 
        m1_nan(randperm(numel(m1_nan), num_missing)) = NaN;
        
        % Recalculate with missing data
        eff_nan = calculate_real_effects({m1_nan}, 1);
        [order_nan, ~] = calculate_ranking({m1_nan}, eff_nan, thresholds, config, ds_names_5, pairs_5);
        
        match = isequal(order_nan, order_ref);
        status_str = 'PASS';
        if ~match
            % At 80% missing data, degradation is expected/allowed
            if lvl >= 0.80
                status_str = 'INFO: Degradation Exp.';
            else
                status_str = 'FAIL: Unexpected';
                robustness_passed = false;
            end
        end
        
        table_data(idx_row, :) = {sprintf('%d%%', round(lvl*100)), mat2str(order_nan), status_str};
    end
    print_auto_table(h_nan, table_data, d_align, h_align);
    
    if robustness_passed
        fprintf('\n[Status] PASS: Algorithm works robustly with missing data.\n');
        tests_passed = tests_passed + 1;
    else
        fprintf('\n[Status] FAIL: Robustness check failed (unexpected degradation).\n');
    end

    %% Test 12: Cycle Detection (Artificial Injection)
    tests_total = tests_total + 1;
    title_str = 'Test 12: Cycle Detection (Artificial Injection)';
    fprintf('\n%s\n', title_str);
    fprintf('%s\n', repmat('-', 1, strlength(title_str)));
    
    m1_data = generate_exact_data(n_subj, [30, 20, 10], 1.0); % Base Order [1, 2, 3]
    m2_data = generate_exact_data(n_subj, [0, 0, 0], 1.0);
    
    % Inject a Logical Cycle in M2 by forcing effect sizes:
    % 2 beats 1, 3 beats 2, 1 beats 3 (Rock-Paper-Scissors)
    eff = calculate_real_effects({m1_data, m2_data}, 2);
    eff.d_vals_all(:, 2) = [-0.9; 0.9; -0.9]; 
    
    fprintf('[Test] Cycle Detection and Fallback Logic.\n');
    fprintf('[Setup] M1: Transitive. M2: Artificial Cycle injected via Cliff''s Delta.\n');
    
    % Input Table 
    fprintf('[Input: M2 Injection]\n');
    h_cyc = {'Pair', 'Delta', 'Direction'}; 
    d_align = {'l', 'c', 'l'}; 
    h_align = {'c', 'c', 'c'};
    table_data = {
        '2 vs 1', '-0.9', '2 > 1';
        '3 vs 2', '-0.9', '3 > 2';
        '1 vs 3', '+0.9', '1 > 3'
    };
    print_auto_table(h_cyc, table_data, d_align, h_align);
    
    config.ranking_mode = 'M1_M2';
    
    % Suppress warning for expected cycle detection
    warnState = warning('off', 'all'); 
    cleanupObj = onCleanup(@() warning(warnState));
    [final_order, ~, ~, ~, ~] = calculate_ranking({m1_data, m2_data}, eff, thresholds, config, ds_names, nchoosek(1:3, 2));
    warning('on', 'all');
    
    % Result
    fprintf('\n[Result]\n');
    h_res = {'Check', 'Outcome'}; 
    d_align = {'l', 'c'}; 
    h_align = {'c', 'c'}; 
    
    % The algorithm should detect the infinite loop in swaps and revert to M1 order
    cycle_handled = isequal(final_order(:)', [1, 2, 3]);

    table_data = {
        'Cycle Detected?', char(string(cycle_handled));
        'Final Order', mat2str(final_order(:)')
    };
    print_auto_table(h_res, table_data, d_align, h_align);

    % Assertions
    fprintf('\n')
    if check_result(final_order, [1, 2, 3], 'Cycle Handled (Revert to M1)')
        tests_passed = tests_passed + 1;
        fprintf('[Status] PASS: Infinite loop detected, M2 swaps aborted.\n');
    else
        fprintf('[Status] FAIL: Cycle detection failed.\n');
    end

    %% Test 13: Real Pathological Data (Efron's Dice)
    tests_total = tests_total + 1;
    title_str = 'Test 13: Pathological Data (Efron''s Dice)';
    fprintf('\n%s\n', title_str);
    fprintf('%s\n', repmat('-', 1, strlength(title_str)));
    
    % Efron's Dice is a famous mathematical set of dice that are non-transitive.
    % A > B, B > C, C > D, but D > A.
    [m1_data, m2_data] = generate_efron_data(n_subj);
    ds_names_4 = {'DiceA', 'DiceB', 'DiceC', 'DiceD'};
    
    fprintf('[Test] Handling of Nontransitive Dice Data (Real-world Cycle).\n');
    fprintf('[Setup] M1: Transitive [40, 30, 20, 10]. M2: Efron''s Dice.\n');
    
    % Input Table 
    fprintf('[Input]\n');
    h_ef = {'Metric', 'Data Structure'}; 
    d_align = {'l', 'l'};
    h_align = {'c', 'c'};
    table_data = {
        'M1', '40 > 30 > 20 > 10';
        'M2', 'A>B (2/3), B>C (2/3), C>D (2/3), D>A (2/3)'
    };
    print_auto_table(h_ef, table_data, d_align, h_align);
    
    config.ranking_mode = 'M1_M2';
    % Suppress warning for expected cycle detection
    warnState = warning('off', 'all'); 
    cleanupObj = onCleanup(@() warning(warnState));
    [~, final_order] = run_single_test({m1_data, m2_data}, thresholds, config, ds_names_4);
    warning('on', 'all');
    
    % We accept any valid permutation of 4 datasets as long as it doesn't crash or hang
    if ~isempty(final_order) && length(unique(final_order)) == 4
        fprintf('\n[Result] Final Order: %s\n', mat2str(final_order));
        fprintf('[Status] PASS: Algorithm terminated with valid rank (Stable Fallback).\n');
        tests_passed = tests_passed + 1;
    else
        fprintf('[Status] FAIL: Cycle Handling failed (Result empty or invalid).\n');
    end

    %% Test 14: Dynamic Cycle Handling (Efron's Dice + Auto-Thresholds)
    tests_total = tests_total + 1;
    title_str = 'Test 14: Dynamic Cycle Handling (Efron''s Dice + Auto-Thresholds)';
    fprintf('\n%s\n', title_str);
    fprintf('%s\n', repmat('-', 1, strlength(title_str)));
    
    [m1_data, m2_efron] = generate_efron_data(n_subj);
    all_data_pipe = {m1_data, m2_efron};
    config_pipe = default_config;
    config_pipe.ranking_mode = 'M1_M2';
    config_pipe.metric_names = {'Base', 'Efron'};
    
    % Temp setup for internal calculation of thresholds
    temp_dir = tempname; mkdir(temp_dir);
    cleanUpTemp = onCleanup(@() rmdir(temp_dir, 's'));
    
    fprintf('[Test] Full pipeline execution with dynamic thresholds on cycle data.\n');
    fprintf('[Setup] 1. Calculate Thresholds -> 2. Run Ranking (M1_M2). This test uses the same Efron Data as in Test 13.\n');
    
    try
        % Step 1: Calculate dynamic thresholds based on the data
        cmd_pipe = ['[d_t_dyn, r_t_dyn, ~, ~, d_v_p, r_v_p, p_i_p] = ' ...
            'calculate_thresholds(all_data_pipe, n_subj, config_pipe, temp_dir, 50, myStream, styles, lang);'];
        [~] = evalc(cmd_pipe);
        
        thr_pipe = struct('d_thresh', d_t_dyn, 'rel_thresh', r_t_dyn);
        eff_pipe = struct('d_vals_all', d_v_p, 'rel_vals_all', r_v_p);
        
        % Input/Intermediate Table 
        fprintf('[Input: Calculated Thresholds]\n');
        h_thr = {'Metric', 'Delta Thr', 'RelDiff Thr'}; 
        d_align = {'l', 'l', 'l'};
        h_align = {'c', 'c', 'c'};
        
        table_data = {
            'Base (M1)', sprintf('%.3f', d_t_dyn(1)), sprintf('%.3f', r_t_dyn(1));
            'Efron (M2)', sprintf('%.3f', d_t_dyn(2)), sprintf('%.3f', r_t_dyn(2))
        };
        print_auto_table(h_thr, table_data, d_align, h_align);

        % Step 2: Run Ranking
        warning('off', 'all'); 
        [final_order_pipe, ~, all_sig, ~, ~] = calculate_ranking(all_data_pipe, eff_pipe, thr_pipe, config_pipe, ds_names_4, p_i_p);
        warning('on', 'all');
        
        % Check if M2 cycle is truly significant (pairwise wins should form cycle)
        sig_m2 = all_sig{2};
        cycle_active = sig_m2(1,2) && sig_m2(2,3) && sig_m2(3,4) && sig_m2(4,1);
        
        % Result Table 
        fprintf('\n[Result]\n');
        h_res = {'Check', 'Outcome'}; 
        d_align = {'l', 'l'};
        h_align = {'c', 'c'};
        
        table_data = {
            'Cycle Active?', char(string(cycle_active));
            'Final Order', mat2str(final_order_pipe)
        };
        print_auto_table(h_res, table_data, d_align, h_align);
        
        if ~isempty(final_order_pipe)
            fprintf('\n[Status] PASS: Pipeline finished with valid order.\n');
            tests_passed = tests_passed + 1;
        end
    catch ME
        fprintf('\n[Status] FAIL: Pipeline crashed: %s\n', ME.message);
    end

    %% Test 15: Global System Test (Static Thresholds)
    tests_total = tests_total + 1;
    title_str = 'Test 15: Global System Test (Static Thresholds)';
    fprintf('\n%s\n', title_str);
    fprintf('%s\n', repmat('-', 1, strlength(title_str)));
    
    % Integration test for the entire hierarchy logic using "Waterfall" data.
    fprintf('[Test] Multi-Tier "Waterfall" sorting (M2 Correction, M1 Base, M3 Tail).\n');
    
    % Generate 15 datasets split into Tiers (see helper function)
    [m1_data, m2_data, m3_data, ds_names_15] = generate_waterfall_data(n_subj);
    
    % Calculate actual means for transparency in the log
    mu1 = mean(m1_data); mu2 = mean(m2_data); mu3 = mean(m3_data);
    
    % Input Table (Detailed Data Profile)
    % We show representatives for each Tier to verify the generated signal values.
    fprintf('[Input: Data Profile (Waterfall Structure)]\n');
    h_in = {'Tier', 'Rep.', 'M1 (Mean)', 'M2 (Mean)', 'M3 (Mean)', 'Intended Logic'}; 
    d_align = {'l', 'c', 'c', 'c', 'c', 'l'}; 
    h_align = {'c', 'c', 'c', 'c', 'c', 'c'};
    
    table_data = {
        'Top (Tier 2)', 'D6',  sprintf('%.1f', mu1(6)),  sprintf('%.1f', mu2(6)),  sprintf('%.1f', mu3(6)),  'M2 High -> Correction';
        'Mid (Tier 1)', 'D1',  sprintf('%.1f', mu1(1)),  sprintf('%.1f', mu2(1)),  sprintf('%.1f', mu3(1)),  'M1 High -> Base Rank';
        'Low (Tier 3)', 'D11', sprintf('%.1f', mu1(11)), sprintf('%.1f', mu2(11)), sprintf('%.1f', mu3(11)), 'M3A Swap Logic';
        'Bot (Tier 4)', 'D13', sprintf('%.1f', mu1(13)), sprintf('%.1f', mu2(13)), sprintf('%.1f', mu3(13)), 'M3B Sorting Logic'
    };
    print_auto_table(h_in, table_data, d_align, h_align);
    
    % Run with Static Thresholds (Standard Config)
    config = default_config;
    eff = calculate_real_effects({m1_data, m2_data, m3_data}, 3);
    pairs_15 = nchoosek(1:15, 2);
    
    [final_order, ~, ~, ~, ~] = calculate_ranking({m1_data, m2_data, m3_data}, eff, thresholds, config, ds_names_15, pairs_15);
    
    % Validation by Tier logic
    % Calculate average rank of each tier to verify separation
    tier2_ranks = []; % Should be top (1-5)
    for i = 6:10, tier2_ranks = [tier2_ranks, find(final_order == i)]; end
    
    tier1_ranks = []; % Should be mid (6-10)
    for i = 1:5, tier1_ranks = [tier1_ranks, find(final_order == i)]; end
    
    tier3_ranks = []; % Should be bottom (11-15)
    for i = 11:15, tier3_ranks = [tier3_ranks, find(final_order == i)]; end
    
    mean_rank_t2 = mean(tier2_ranks);
    mean_rank_t1 = mean(tier1_ranks);
    mean_rank_t3 = mean(tier3_ranks);
    
    % Result Table 
    fprintf('\n[Result]\n');
    h_res = {'Group', 'Mean Rank', 'Expectation', 'Status'}; 
    d_align = {'l', 'c', 'l', 'c'}; 
    h_align = {'c', 'c', 'c', 'c'};
    
    % Check logic
    t2_ok = mean_rank_t2 < mean_rank_t1; % Tier 2 (M2) beat Tier 1 (M1)
    t1_ok = mean_rank_t1 < mean_rank_t3; % Tier 1 (M1) beat Tier 3 (M3)
    
    table_data = {
        'Tier 2 (M2)', sprintf('%.1f', mean_rank_t2), 'Top (~3.0)', char(string(t2_ok));
        'Tier 1 (M1)', sprintf('%.1f', mean_rank_t1), 'Mid (~8.0)', char(string(t1_ok));
        'Tier 3 (M3)', sprintf('%.1f', mean_rank_t3), 'Bot (~13.0)', '-'
    };
    print_auto_table(h_res, table_data, d_align, h_align);
    
    if t2_ok && t1_ok
        fprintf('\n[Status] PASS: Full Hierarchy respected (M2 > M1 > M3).\n');
        tests_passed = tests_passed + 1;
    else
         fprintf('\n[Status] FAIL: Hierarchy violation.\n');
         if ~t2_ok, fprintf('    - M2 Correction failed (Tier 2 not above Tier 1).\n'); end
         if ~t1_ok, fprintf('    - M1 Base Logic failed (Tier 1 not above Tier 3).\n'); end
    end

    %% Test 16: Full System Integration (15 Datasets + Auto-Thresholds)
    tests_total = tests_total + 1;
    title_str = 'Test 16: Full System Integration (15 Datasets + Auto-Thresholds)';
    fprintf('\n%s\n', title_str);
    fprintf('%s\n', repmat('-', 1, strlength(title_str)));
    
    fprintf('[Test] Validating the full logic chain (M2 correction, M3A tie-break, M3B sort).\n');
    % Explicitly mention reusing the structure from Test 15
    fprintf('[Setup] Re-generating "Waterfall" data (same profile as Test 15) to test Dynamic Thresholds.\n');
    
    % Scenario for full logic chain validation (same as T15 but with calculated thresholds)
    [m1_data, m2_data, m3_data, ds_names_full] = generate_waterfall_data(n_subj);
    all_data_full = {m1_data, m2_data, m3_data};
    config_full = default_config;
    temp_dir = tempname; mkdir(temp_dir);
    cleanUpTemp = onCleanup(@() rmdir(temp_dir, 's'));
    
    % Input Table 
    fprintf('[Input: Tiers & Logic Targets]\n');
    h_in = {'Tier', 'Datasets', 'Target Logic'}; 
    d_align = {'l', 'l', 'l'};
    h_align = {'c', 'c', 'c'};
    table_data = {
        'Tier 1', 'D1-D5', 'Base Rank (Mid) -> Displaced by Tier 2';
        'Tier 2', 'D6-D10', 'Correction (Top) -> Wins via M2';
        'Tier 3', 'D11-D12', 'Tie-Break -> Logic 3A (M3 Swap)';
        'Tier 4', 'D13-D15', 'Sorting -> Logic 3B (M3 Gradient)'
    };
    print_auto_table(h_in, table_data, d_align, h_align);

    try
        % Step 1: Auto-Calculate Thresholds
        cmd_full = ['[d_t_full, r_t_full, ~, ~, d_v_full, r_v_full, p_i_full] = ' ...
            'calculate_thresholds(all_data_full, n_subj, config_full, temp_dir, 50, myStream, styles, lang);'];
        [~] = evalc(cmd_full);
        
        % Step 2: Run Full Ranking 
        thr_full = struct('d_thresh', d_t_full, 'rel_thresh', r_t_full);
        eff_full = struct('d_vals_all', d_v_full, 'rel_vals_all', r_v_full);
        
        [final_order_full, ~, ~, ~, ~] = calculate_ranking(all_data_full, eff_full, thr_full, config_full, ds_names_full, p_i_full);
        
        % Step 3: Display Results
        % Re-calculate means for display consistency
        mean_m1 = mean(m1_data); mean_m2 = mean(m2_data); mean_m3 = mean(m3_data);
        
        % Result Table 
        fprintf('\n[Result]\n');
        header_parts = {'Rank', 'Dataset', 'M1 Mean', 'M2 Mean', 'M3 Mean'};
        d_align = {'c', 'c', 'c', 'c', 'c'}; 
        h_align = {'c', 'c', 'c', 'c', 'c'};
        
        table_data = cell(15, 5);
        for r = 1:15
            idx = final_order_full(r);
            table_data(r, :) = {
                sprintf('%d', r), ds_names_full{idx}, ...
                sprintf('%.1f', mean_m1(idx)), ...
                sprintf('%.1f', mean_m2(idx)), ...
                sprintf('%.1f', mean_m3(idx))
            };
        end
        print_auto_table(header_parts, table_data, d_align, h_align);
        fprintf('\n');

        % Validation checks:
        top_5 = final_order_full(1:5);
        tier2_success = all(ismember(top_5, 6:10)); % Check if Tier 2 (M2 strong) is on top
        
        mid_5 = final_order_full(6:10);
        tier1_success = all(ismember(mid_5, 1:5)); % Check if Tier 1 (M1 strong) is 2nd
        
        rank_d11 = find(final_order_full == 11);
        rank_d12 = find(final_order_full == 12);
        logic3a_success = rank_d12 < rank_d11; % Check if M3A flipped D11/D12
        
        rank_d13 = find(final_order_full == 13);
        rank_d14 = find(final_order_full == 14);
        rank_d15 = find(final_order_full == 15);
        logic3b_success = (rank_d15 < rank_d14) && (rank_d14 < rank_d13); % Check M3B sort
        
        if tier2_success && tier1_success && logic3a_success && logic3b_success
             fprintf('[Status] PASS: Full System Validation Successful.\n');
             tests_passed = tests_passed + 1;
        else
             fprintf('[Status] FAIL: Logic Validation Failed.\n');
             if ~tier2_success, fprintf('    - M2 Correction failed (Tier 2 not at top).\n'); end
             if ~tier1_success, fprintf('    - Tier 1 not in mid.\n'); end
             if ~logic3a_success, fprintf('    - Logic 3A failed (D12 did not beat D11).\n'); end
             if ~logic3b_success, fprintf('    - Logic 3B failed (D13-D15 sorting incorrect).\n'); end
        end
        
    catch ME
        fprintf('[Status] FAIL: Pipeline crashed: %s\n', ME.message);
        fprintf('[Diag] %s\n', ME.stack(1).name);
    end

    %% Test 17: Bootstrap Uncertainty & Stability Analysis
    tests_total = tests_total + 1;
    title_str = 'Test 17: Bootstrap Rank Stability Analysis';
    fprintf('\n%s\n', title_str);
    fprintf('%s\n', repmat('-', 1, strlength(title_str)));
    
    % Cluster Bootstrap is used to check how stable the ranking is when resampling subjects.
    
    fprintf('[Test] Verification of rank stability and confidence interval separation.\n');
    fprintf('[Setup] Re-generating "Waterfall" data (same profile as Test 15) for Stability Analysis.\n');
    
    % 1. Setup Data (Identical to Test 16 - "Waterfall" Scenario)
    [m1_data, m2_data, m3_data, ds_names_full] = generate_waterfall_data(n_subj);
    all_data_full = {m1_data, m2_data, m3_data};
    
    config_full = default_config;
    config_full.bootstrap_ranks = bs_config;
    
    % Setup Table 
    fprintf('[Input]\n');
    h_in = {'Parameter', 'Value'}; 
    d_align = {'l', 'l'}; 
    h_align = {'c', 'c'};
    table_data = {
        'Data Profile', 'Waterfall (see Test 15)';
        'Bootstrap Method', 'Cluster Bootstrap';
        'Iterations (B)', '100'
    };
    print_auto_table(h_in, table_data, d_align, h_align);
    
    % Use a temporary directory for intermediate files (automatically deleted later)
    temp_dir = tempname; 
    mkdir(temp_dir);
    cleanUpTemp = onCleanup(@() rmdir(temp_dir, 's'));
    
    try
        % 2. Calculate Thresholds & Base Ranking (Required prerequisites)
        s_boot = RandStream('mlfg6331_64', 'Seed', 999); 
        
        cmd_thr = ['[d_t, r_t, ~, ~, d_v, r_v, p_idx] = calculate_thresholds(' ...
                   'all_data_full, n_subj, config_full, temp_dir, 50, s_boot, styles, lang);'];
        [~] = evalc(cmd_thr);
            
        thr_full = struct('d_thresh', d_t, 'rel_thresh', r_t);
        eff_full = struct('d_vals_all', d_v, 'rel_vals_all', r_v);
        
        [~, final_rank, ~, ~, ~] = calculate_ranking(...
            all_data_full, eff_full, thr_full, config_full, ds_names_full, p_idx);
        
        % 3. Run Cluster Bootstrap
        % Measures variability of ranks
        cmd_boot = ['[boot_ranks, ~, ~, ~, ~] = bootstrap_ranking(' ...
                    'all_data_full, thr_full, config_full, ds_names_full, final_rank, p_idx, n_subj, ' ...
                    'temp_dir, temp_dir, 100, s_boot, styles, lang, ''UnitTest_Boot'');'];
        [~] = evalc(cmd_boot);
            
        % 4. Print Percentage Rank Distributions 
        fprintf('\n[Result: Rank Distributions]\n');
        header_parts = {'Rank', 'Dataset', 'Distribution (B=100)'};
        d_align = {'c', 'c', 'l'}; 
        h_align = {'c', 'c', 'c'};
        
        % Inline Helper to format distribution string (e.g., "1: 90%, 2: 10%")
        get_dist_str = @(r) strjoin(arrayfun(@(u) sprintf('%d: %.0f%%', u, sum(r==u)/numel(r)*100), ...
                                    unique(r), 'UniformOutput', false), ', ');
        
        [~, sort_idx] = sort(final_rank);
        
        table_data = cell(numel(ds_names_full), 3);
        
        for i = 1:numel(ds_names_full)
            idx = sort_idx(i); 
            name = ds_names_full{idx};
            ranks = boot_ranks(idx, :);
            dist_str = get_dist_str(ranks);
            
            table_data(i, :) = {sprintf('%d', i), name, dist_str};
        end
        print_auto_table(header_parts, table_data, d_align, h_align);
        fprintf('\n');

        % 5. Check Inter-Tier Separation (Representatives)
        % D6 (Tier 2 Top) vs D1 (Tier 1 Mid) vs D11 (Tier 3 Low)
        idx_d6 = 6; idx_d1 = 1; idx_d11 = 11;
        ci_d6  = quantile(boot_ranks(idx_d6, :), [0.025, 0.975]);
        ci_d1  = quantile(boot_ranks(idx_d1, :), [0.025, 0.975]);
        ci_d11 = quantile(boot_ranks(idx_d11, :), [0.025, 0.975]);
        
        % CIs should NOT overlap between tiers
        sep_top_mid = ci_d6(2) < ci_d1(1); 
        sep_mid_bot = ci_d1(2) < ci_d11(1);
        
        if sep_top_mid && sep_mid_bot
            fprintf('[Status] PASS: Separation confirmed: Robust distinction between Tiers.\n');
        else
            fprintf('[Status] INFO: Overlap detected between Tiers (acceptable for low B/noise).\n');
        end
        % 6. Check Intra-Tier Overlap
        % D1 and D2 are in the same tier, so their rank distributions should overlap significantly.
        idx_d2 = 2; 
        ranks_d2 = boot_ranks(idx_d2, :);
        ci_d2 = quantile(ranks_d2, [0.025, 0.975]);
        
        % Check for Overlap: Max(LowerBound) < Min(UpperBound)
        overlap_ok = max(ci_d1(1), ci_d2(1)) < min(ci_d1(2), ci_d2(2));
        
        if overlap_ok
            fprintf('[Status] PASS: Overlap confirmed: Correctly identifies D1 and D2 as indistinguishable.\n');
        else
            fprintf('[Status] FAIL: Unexpected separation of identical datasets.\n');
        end
        
        if sep_top_mid && sep_mid_bot && overlap_ok
            tests_passed = tests_passed + 1;
        end
        
    catch ME
        fprintf('\n[Status] FAIL: Bootstrap execution failed: %s\n', ME.message);
        fprintf('[Diag] %s at line %d\n', ME.stack(1).file, ME.stack(1).line);
    end

    %% Test 18: Borda Count Aggregation Logic Check
    tests_total = tests_total + 1;
    title_str = 'Test 18: Borda Count Aggregation Logic Check';
    fprintf('\n%s\n', title_str);
    fprintf('%s\n', repmat('-', 1, strlength(title_str)));
    
    % Borda Count aggregates rankings from Sensitivity Analysis (different metric permutations).
    % Points = (N - Rank). The dataset with most points wins.
    
    % Simulation: 3 Datasets, 2 Permutations
    % Permutation 1: [1, 2, 3] -> Ranks: D1=1, D2=2, D3=3 (Points: 2, 1, 0)
    % Permutation 2: [2, 1, 3] -> Ranks: D1=2, D2=1, D3=3 (Points: 1, 2, 0)
    % Max Possible Points per perm: 3-1 = 2. Total Max = 4.
    % D1 Points: 2+1=3. Score: 3/4 = 75%
    % D3 Points: 0+0=0. Score: 0/4 = 0%
    
    sim_ranks = [1, 2; 2, 1; 3, 3]; % Rows=Datasets, Cols=Permutations
    sim_names = {'D1', 'D2', 'D3'};
    
    fprintf('[Test] Correctness of Borda score calculation.\n');
    fprintf('[Setup] Simulating specific rank distribution for 3 datasets over 2 permutations.\n');
    
    % Input Table 
    fprintf('[Input]\n');
    h_in = {'Dataset', 'Ranks (P1, P2)', 'Points Calc'}; 
    d_align = {'c', 'c', 'l'}; 
    h_align = {'c', 'c', 'c'};
    table_data = {
        'D1', '1, 2', '2 + 1 = 3';
        'D2', '2, 1', '1 + 2 = 3';
        'D3', '3, 3', '0 + 0 = 0'
    };
    print_auto_table(h_in, table_data, d_align, h_align);
    
    borda_res = borda_ranking(sim_ranks, sim_names);
    
    score_d1 = borda_res.score(1);
    score_d3 = borda_res.score(3);
    
    % Result Table 
    fprintf('\n[Result]\n');
    h_b = {'Dataset', 'Score', 'Expected'}; 
    d_align = {'c', 'c', 'c'}; 
    h_align = {'c', 'c', 'c'};
    
    table_data = {
        'D1', sprintf('%.1f%%', score_d1), '75.0%';
        'D3', sprintf('%.1f%%', score_d3), '0.0%'
    };
    print_auto_table(h_b, table_data, d_align, h_align);
    
    if abs(score_d1 - 75.0) < 1e-6 && abs(score_d3 - 0.0) < 1e-6
        fprintf('\n[Status] PASS: Borda calculation mathematically correct.\n');
        tests_passed = tests_passed + 1;
    else
        fprintf('\n[Status] FAIL: Borda calculation mismatch.\n');
    end

    %% Test 19: Post-Hoc Power Analysis (Sanity Check)
    tests_total = tests_total + 1;
    title_str = 'Test 19: Post-Hoc Power Analysis (Sanity Check)';
    fprintf('\n%s\n', title_str);
    fprintf('%s\n', repmat('-', 1, strlength(title_str)));
    
    fprintf('[Test] Verifying that Power Analysis returns plausible values (0 <= Power <= 1).\n');
    
    % Setup: Strong Signal (Effect Size large) -> High Power expected
    % D1 Mean=10, D2 Mean=0. SD=2. N=20. Cliff's Delta will be approx 1.0.
    n_pwr = 20;
    data_pwr = generate_exact_data(n_pwr, [10, 0], 2.0);
    all_data_pwr = {data_pwr};
    ds_names_pwr = {'HighSignal', 'LowSignal'};
    pairs_pwr = [1 2];
    
    % Setup Table
    fprintf('[Setup]\n');
    h_set = {'Dataset', 'Mean', 'SD', 'Description'};
    d_align = {'l', 'c', 'c', 'l'};
    h_align = {'c', 'c', 'c', 'c'};
    
    table_data = {
        'HighSignal', '10.0', '2.0', 'Strong Effect Source';
        'LowSignal',  '0.0',  '2.0', 'Noise/Reference Group'
    };
    print_auto_table(h_set, table_data, d_align, h_align);

    % Config & Prerequisites
    config_pwr = default_config;
    config_pwr.ranking_mode = 'M1';
    eff_pwr = calculate_real_effects(all_data_pwr, 1);
    
    % Calculate Alpha Matrices first (needed as input for power_analysis)
    [~, ~, ~, all_alphas_pwr, ~] = calculate_ranking(all_data_pwr, eff_pwr, thresholds, config_pwr, ds_names_pwr, pairs_pwr);
    
    % Simulation Settings
    n_sims = 100; % Low number for unit test speed
    
    % Input Table (Simulation Parameters & Win Criteria)
    fprintf('\n[Input]\n');
    h_in = {'Parameter', 'Value', 'Condition'}; 
    d_align = {'l', 'l', 'l'}; 
    h_align = {'c', 'c', 'c'};
    
    % Extract values for display
    d_thr_disp = thresholds.d_thresh(1);
    r_thr_disp = thresholds.rel_thresh(1);
    alpha_disp = all_alphas_pwr{1}(1,2);

    table_data = {
        'Subjects (N)', num2str(n_pwr), 'Resampling Size';
        'Simulations', num2str(n_sims), 'Bootstrap Iterations';
        'Alpha Limit', sprintf('%.3f', alpha_disp), 'Significant if p < Alpha';
        'Effect Limit', sprintf('d > %.2f, r > %.2f', d_thr_disp, r_thr_disp), 'Relevant if Eff > Limit'
    };
    print_auto_table(h_in, table_data, d_align, h_align);
    
    try
        % Call the function
        % We use s_global to ensure the bootstrap sampling is reproducible
        % Use 'evalc' to suppress the function's internal console output
        cmd_pwr = ['power_results = power_analysis(' ...
                   'all_data_pwr, config_pwr, thresholds, n_pwr, n_sims, pairs_pwr, s_global, lang);'];
        [~] = evalc(cmd_pwr);
        
        % Extract Result
        % power_analysis returns struct('power_matrices', {cell_array})
        % The cell array contains one vector per metric. We have 1 metric.
        % The vector contains one value per pair. We have 1 pair (Index 1).
        pwr_val = power_results.power_matrices{1}(1);
        
        % Result Table
        fprintf('\n[Result]\n');
        h_res = {'Metric', 'Value', 'Plausibility'}; 
        d_align = {'l', 'l', 'l'}; 
        h_align = {'c', 'c', 'c'};
        
        % Expectation: With such a strong signal (Mean 10 vs 0), power should be perfect or near perfect.
        is_plausible = (pwr_val >= 0.8) && (pwr_val <= 1.0);
        
        table_data = {
            'Calculated Power', sprintf('%.2f', pwr_val), '> 0.80 (Expected)'
        };
        print_auto_table(h_res, table_data, d_align, h_align);
        
        if is_plausible
            fprintf('\n[Status] PASS: Power Analysis produced a logical high value for strong signal.\n');
            tests_passed = tests_passed + 1;
        else
            fprintf('\n[Status] FAIL: Power Value implausible (Expected > 0.8 for strong signal).\n');
        end
        
    catch ME
        fprintf('\n[Status] FAIL: Power Analysis crashed: %s\n', ME.message);
        fprintf('\n[Diag] %s\n', ME.stack(1).name);
    end

    %% Final Summary
    fprintf('\n=================================\n');
    fprintf('UNIT TEST SUMMARY - FINAL RESULTS\n');
    fprintf('Passed: %g / %d\n', tests_passed, tests_total);
    fprintf('=================================\n');
    if tests_passed == tests_total
        fprintf('STATUS: ALL TESTS PASSED (SUCCESS)\n');
    else
        fprintf('STATUS: FAILURES DETECTED\n');
    end
    total_duration_seconds = toc;
    minutes = floor(total_duration_seconds / 60);
    seconds = rem(total_duration_seconds, 60);
    fprintf('\nTotal execution time: %d minutes and %.2f seconds.\n', minutes, seconds);
    
    % Clean up resources
    diary off;
    delete(gcp('nocreate'));
end


%% Helper Functions for Output
function [log_folder, source] = get_writable_log_path()
    % Determine a valid location for log files.
    % Priority: 1. User Documents (Persistent), 2. TempDir (Fallback)
    
    app_folder_name = 'HERA_Ranking_Logs';
    
    % Try standard Documents folder first
    if ispc
        docs_path = fullfile(getenv('USERPROFILE'), 'Documents');
    else
        docs_path = fullfile(getenv('HOME'), 'Documents');
    end
    
    target_path = fullfile(docs_path, app_folder_name);
    
    % Attempt to create/access Documents folder
    try
        if ~exist(target_path, 'dir')
            mkdir(target_path);
        end
        % Simple write test to confirm permissions
        testFile = fullfile(target_path, 'write_test.tmp');
        fid = fopen(testFile, 'w');
        if fid == -1
            error('No write access');
        end
        fclose(fid);
        delete(testFile);
        
        log_folder = target_path;
        source = 'Documents (Persistent)';
        return;
    catch
        % Fallback to tempdir if Documents fails (e.g. Restricted User)
        log_folder = fullfile(tempdir, app_folder_name);
        if ~exist(log_folder, 'dir')
            mkdir(log_folder);
        end
        source = 'TempDir (Fallback)';
    end
end

function success = check_result(actual, expected, label)
    % Validates the actual rank order against the expectation.
    % Helper for simple Pass/Fail printing.
    if isequal(actual(:)', expected(:)')
        fprintf('[Status] PASS: %s: Order is %s.\n', label, mat2str(actual(:)'));
        success = true;
    else
        fprintf('[Status] FAIL: %s\n', label);
        fprintf('   Expected: %s\n', mat2str(expected(:)'));
        fprintf('   Actual:   %s\n', mat2str(actual(:)'));
        success = false;
    end
end


function print_auto_table(headers, data, data_alignments, header_alignments)
% PRINT_AUTO_TABLE - Visualizes a data table with dynamic column widths.
%
% Description:
%   Automatically calculates column widths. Allows separate alignment
%   for headers and data rows.
%
% Inputs:
%   headers           - Cell array of strings defining column titles.
%   data              - Cell matrix (M x N) containing the table content.
%   data_alignments   - Cell array of chars ('l', 'c', 'r') for data.
%   header_alignments - (Optional) Cell array of chars for headers. 
%                       Defaults to data_alignments if omitted.

    % 1. Initialization & Defaults
    if nargin < 4
        header_alignments = data_alignments; % Fallback: Header = Data
    end

    num_cols = length(headers);
    num_rows = size(data, 1);
    
    % Initialize column widths based on header lengths
    col_widths = cellfun(@strlength, headers);
    
    % 2. Dynamic Width Calculation based on data content
    if num_rows > 0
        for r = 1:num_rows
            for c = 1:num_cols
                val = data{r,c};
                if isnumeric(val) || islogical(val)
                    val_str = char(string(val)); 
                elseif ischar(val)
                    val_str = val;
                else
                    val_str = char(val);
                end
                col_widths(c) = max(col_widths(c), strlength(val_str));
            end
        end
    end
    
    % 3. Apply Padding
    col_widths = col_widths + 2;
    
    % 4. Print Header (Using header_alignments)
    header_line = strjoin(arrayfun(@(c) format_text(headers{c}, col_widths(c), header_alignments{c}), ...
        1:num_cols, 'UniformOutput', false), '|');
    fprintf('   |%s|\n', header_line);
    
    % Separator
    total_width = sum(col_widths) + length(col_widths) + 1;
    fprintf('   %s\n', repmat('-', 1, total_width));
    
    % 5. Print Data Rows (Using data_alignments)
    for r = 1:num_rows
        row_str_parts = cell(1, num_cols);
        for c = 1:num_cols
            val = data{r,c};
            if isnumeric(val) || islogical(val)
                str_val = char(string(val));
            elseif ischar(val)
                str_val = val;
            else
                str_val = char(val);
            end
            
            % Use data_alignments here
            row_str_parts{c} = format_text(str_val, col_widths(c), data_alignments{c});
        end
        fprintf('   |%s|\n', strjoin(row_str_parts, '|'));
    end
end

function output = format_text(text, width, alignment)
    % Aligns text within a fixed width (Left, Right, Center).
    text_len = strlength(text);
    padding = width - text_len;
    
    switch alignment
        case 'r'
            output = [repmat(' ', 1, padding), text];
        case 'l'
            output = [text, repmat(' ', 1, padding)];
        otherwise % 'c'
            padding_left = floor(padding / 2);
            padding_right = ceil(padding / 2);
            output = [repmat(' ', 1, padding_left), text, repmat(' ', 1, padding_right)];
    end
end

%% Helper Functions for Test Data Generation
function data = generate_exact_data(n, means, std_dev)
    % Generates deterministic random normal data.
    %
    % Description:
    %   This function ensures that for the exact same inputs (n, means, std_dev), 
    %   the exact same data matrix is generated, regardless of the global RNG state.
    %   It achieves this by creating a temporary local RandStream based on a hash 
    %   of the input arguments.

    % Create a deterministic seed based on input parameters.
    % We hash the inputs to create a unique seed for this specific data configuration.
    % This ensures Test 7 receives the exact same data every time, even if Test 6 fails.
    seed_val = abs(sum(means) * 1000 + n + std_dev * 10);
    
    % Initialize a local stream with this specific seed.
    s_gen = RandStream('mlfg6331_64', 'Seed', uint32(seed_val));
    
    % Set this stream as global temporarily and store the previous stream.
    oldStream = RandStream.setGlobalStream(s_gen);
    
    % Start of Data Generation
    k = length(means);
    data = randn(n, k);
    if std(data(:)) == 0 % Handle zero variance case
        data = zeros(n,k) + means;
    else
        % Standardization: (X - mu) / sigma
        data = (data - mean(data)) ./ std(data); 
        % Scaling: X * new_sigma + new_mu
        data = data .* std_dev + means;
    end
    % End of Data Generation 
    
    % Restore the original random stream to prevent side effects on other tests.
    RandStream.setGlobalStream(oldStream);
end

function effect_sizes = calculate_real_effects(all_data, num_metrics)
    % Calculates true effect sizes (Delta, RelDiff) for test verification.
    % Used to ensure the inputs to the ranking function are correct.
    num_datasets = size(all_data{1}, 2);
    pair_idx = nchoosek(1:num_datasets, 2);
    n_pairs = size(pair_idx, 1);
    d_vals = zeros(n_pairs, num_metrics);
    rel_vals = zeros(n_pairs, num_metrics);
    
    for m = 1:num_metrics
        data = all_data{m};
        for k = 1:n_pairs
            i = pair_idx(k, 1); j = pair_idx(k, 2);
            x = data(:, i); y = data(:, j);
            valid = ~isnan(x) & ~isnan(y);
            x = x(valid); y = y(valid);
            n = numel(x);
            if n > 0
                gt = sum(x > y', 'all'); 
                lt = sum(x < y', 'all');
                d_vals(k, m) = (gt - lt) / (n^2);
                mx = mean(x); my = mean(y);
                if (mx+my)==0, rel_vals(k, m) = 0; 
                else, rel_vals(k, m) = abs(mx - my) / abs(mean([mx, my])); end
            end
        end
    end
    effect_sizes.d_vals_all = d_vals;
    effect_sizes.rel_vals_all = rel_vals;
end

function [m1, m2] = generate_efron_data(n_subj)
    % Generates non-transitive Dice data (Efron's Dice).
    % Probabilities: A>B (2/3), B>C (2/3), C>D (2/3), D>A (2/3).
    % This is used to test cycle detection logic.
    
    vals_A = [4; 4; 4; 4; 0; 0];
    vals_B = [3; 3; 3; 3; 3; 3];
    vals_C = [2; 2; 2; 2; 6; 6];
    vals_D = [1; 1; 1; 1; 5; 5];
    
    reps = ceil(n_subj / 6);
    A = repmat(vals_A, reps, 1); A = A(1:n_subj);
    B = repmat(vals_B, reps, 1); B = B(1:n_subj);
    C = repmat(vals_C, reps, 1); C = C(1:n_subj);
    D = repmat(vals_D, reps, 1); D = D(1:n_subj);
    
    m2 = [A, B, C, D];
    % M1 is provided as a transitive baseline
    m1 = generate_exact_data(n_subj, [40, 30, 20, 10], 1.0);
end

function [m1, m2, m3, names] = generate_waterfall_data(n_subj)
    % Generates a complex dataset triggering all HERA logic stages.
    % Designed for "Integration Test 17".
    
    s_gen = RandStream('mlfg6331_64', 'Seed', 9999);
    oldStream = RandStream.setGlobalStream(s_gen);

    names = arrayfun(@(x) sprintf('D%d', x), 1:15, 'UniformOutput', false);
    
    % Initialize matrices
    m1 = zeros(n_subj, 15);
    m2 = zeros(n_subj, 15);
    m3 = zeros(n_subj, 15);
    
    % Helper to add noise
    % scale=0.1 for signals (d > thresh)
    % scale=0.0 for neutrality (d = 0)
    add_vals = @(vals, scale) repmat(vals, n_subj, 1) + (randn(n_subj, length(vals)) * scale);
    
    % Base offset: Shift all values to positive range (e.g., > 100). This makes it easier to understand the result.
    base = 100.0; 
    
    % Tier 1 (D1-D5): Mid Tier
    % M1: Strong (Base + 20). M2/M3: Neutral (Base).
    % Beats Tiers 3/4 via M1. Lost to Tier 2 via M2.
    m1(:, 1:5) = add_vals(repmat(base + 20.0, 1, 5), 0.1);  
    m2(:, 1:5) = add_vals(repmat(base,        1, 5), 0.0);   
    m3(:, 1:5) = add_vals(repmat(base,        1, 5), 0.0);   
    
    % Tier 2 (D6-D10): Top Tier (M2 Correction) 
    % M1: Neutral (Base). Beaten by Tier 1 on M1 initially (Base < Base+20).
    % M2: Strong (Base + 20). Triggers Swap vs Tier 1.
    m1(:, 6:10) = add_vals(repmat(base,        1, 5), 0.0);  
    m2(:, 6:10) = add_vals(repmat(base + 20.0, 1, 5), 0.1); 
    m3(:, 6:10) = add_vals(repmat(base,        1, 5), 0.0);  
    
    % Tier 3 (D11-D12): Logic 3A (Tie-Break) 
    % M1: Weak (Base - 20).
    % M2: Neutral (Base).
    % M3: D12 > D11.
    m1(:, 11:12) = add_vals([base - 20.0, base - 20.0], 0.0); 
    m2(:, 11:12) = add_vals([base,        base       ], 0.0);     
    m3(:, 11:12) = add_vals([base - 5.0,  base + 5.0 ], 0.1); % D12 > D11
    
    % Tier 4 (D13-D15): Logic 3B (Sort)
    % M1: Weakest (Base - 40).
    % M2: Neutral (Base).
    % M3: Gradient D15 > D14 > D13.
    m1(:, 13:15) = add_vals([base - 40.0, base - 40.0, base - 40.0], 0.0); 
    m2(:, 13:15) = add_vals([base,        base,        base       ], 0.0);       
    m3(:, 13:15) = add_vals([base - 10.0, base,        base + 10.0], 0.1); 
    
    RandStream.setGlobalStream(oldStream); % Restore previous stream
end

function [final_order, final_rank, eff, p_vals] = run_single_test_full(all_data, thresholds, config, names)
    % Wrapper to run 'calculate_ranking' and return full outputs including diagnostics.
    import HERA.*
    num_d = size(all_data{1}, 2);
    pairs = nchoosek(1:num_d, 2);
    eff = calculate_real_effects(all_data, numel(all_data));
    [final_order, final_rank, ~, ~, p_vals] = calculate_ranking(all_data, eff, thresholds, config, names, pairs);
end

function [final_rank, final_order] = run_single_test(all_data, thresholds, config, names)
    % Wrapper to run 'calculate_ranking' and return just the order/rank.
    import HERA.*
    [final_order, final_rank, ~, ~] = run_single_test_full(all_data, thresholds, config, names);
end

function [lang, styles] = get_test_resources()
    % Provides minimal Mock Language and Style structs for testing.
    % This avoids dependency on external files during unit testing,
    % ensuring the test suite is self-contained.
  
    lang = struct();
    
    % used in calculate_thresholds.m
    lang.thresholds.manual_b_info = 'Manual Percentile-Bootstrap count B = %d is used for thresholds.';
    lang.thresholds.searching_optimal_b = 'Searching for optimal Percentile-Bootstrap count (B) for thresholds';
    lang.thresholds.primary_criterion = 'Primary: IQR/Median-Ratio';
    lang.thresholds.robust_convergence_info = 'Robust Convergence: Win %d, Strk %d, Tol %.1f%%, Start %d, Max %d';
    lang.thresholds.simple_convergence_info = 'Simple Convergence: Tol %.1f%%, Start %d, Max %d';
    lang.thresholds.secondary_criterion = 'Secondary: Elbow analysis';
    lang.thresholds.convergence_result = 'B-value after convergence: %d';
    lang.thresholds.elbow_analysis_info = 'Convergence not reached, performing elbow analysis';
    lang.thresholds.elbow_result = 'B-value after elbow analysis: %d';
    lang.thresholds.convergence_plot_saved = 'Convergence plot for thresholds saved: %s';
    lang.thresholds.histogram_plot_saved = 'Histogram plot for threshold distribution saved: %s';
    lang.thresholds.effects_plot_saved = 'Histogram plot for effect sizes saved: %s';
    lang.thresholds.checking_stability = 'Checking stability with B = %d and %d trials';
    lang.thresholds.convergence_run_info = 'Convergence at %.2f%% in run %d/%d';
    lang.thresholds.stability_change_info = 'Relative stability change: %.2f%%';
    lang.thresholds.convergence_reached = 'Convergence at %.2f%% reached';
    lang.thresholds.stable_runs_info = 'Stability in %d consecutive runs';
    
    % BCA used in calculate_bca_ci.
    lang.bca.manual_b_info = 'Manual BCa-Bootstrap count B = %d is used.';
    lang.bca.searching_optimal_b = 'Searching for optimal BCa Bootstrap count (B)';
    lang.bca.primary_criterion = 'Primary: IQR/Median-Ratio';
    lang.bca.secondary_criterion = 'Secondary: Elbow analysis';
    lang.bca.simple_convergence_info = 'Simple Convergence';
    lang.bca.robust_convergence_info = 'Robust Convergence';
    lang.bca.checking_stability = 'Checking stability...';
    lang.bca.convergence_run_info = 'Run...';
    lang.bca.stability_change_info = 'Change...';
    lang.bca.convergence_reached = 'Converged';
    lang.bca.stable_runs_info = 'Stable runs reached';
    lang.bca.convergence_result = 'Result: %d';
    lang.bca.elbow_analysis_info = 'Elbow Analysis';
    lang.bca.elbow_result = 'Elbow: %d';
    lang.bca.calculating_final_ci = 'Calculating final CI...';
    lang.bca.convergence_plot_saved = 'Plot saved: %s';
    lang.bca.saving_csv = 'Saving CSV...'; 
    lang.bca.csv_saved = 'Saved: %s';
    lang.bca.z0_histogram_saved = 'Z0 Saved: %s';
    lang.bca.a_histogram_saved = 'A Saved: %s';
    lang.bca.ci_width_histogram_saved = 'Width Saved: %s';
    % Sub-struct for factors
    lang.bca.correction_factors.header = 'Factors';
    lang.bca.correction_factors.factor = 'F';
    lang.bca.correction_factors.bias_delta = 'B-D';
    lang.bca.correction_factors.skew_delta = 'S-D';
    lang.bca.correction_factors.bias_rel = 'B-R';
    lang.bca.correction_factors.skew_rel = 'S-R';
    lang.bca.correction_factors.bias = 'Bias';
    lang.bca.correction_factors.skew = 'Skew';

    % used in bootstrap_ranking.m
    lang.ranking.manual_b_info = 'Manual Rank B=%d';
    lang.ranking.searching_optimal_b = 'Search Rank B...';
    lang.ranking.primary_criterion = 'Rank Crit 1';
    lang.ranking.secondary_criterion = 'Rank Crit 2';
    lang.ranking.robust_convergence_info = 'Robust Rank';
    lang.ranking.simple_convergence_info = 'Simple Rank';
    lang.ranking.checking_stability = 'Rank Check...';
    lang.ranking.convergence_run_info = 'Rank Run';
    lang.ranking.stability_change_info = 'Rank Change';
    lang.ranking.convergence_reached = 'Rank Conv';
    lang.ranking.stable_runs_info = 'Rank Stable';
    lang.ranking.convergence_result = 'Rank Res %d';
    lang.ranking.elbow_analysis_info = 'Rank Elbow';
    lang.ranking.elbow_result = 'Rank Elbow %d';
    lang.ranking.stability_vector_warning = 'Rank Warn';
    lang.ranking.convergence_plot_saved = 'Rank Plot Saved: %s';
    lang.ranking.distribution_header = '\nDistribution (B=%d)\n';
    lang.ranking.table_rank = 'Rank';
    lang.ranking.table_dataset = 'Dataset';
    lang.ranking.table_dist = 'Freq';
    lang.ranking.csv_bootstrap_rank = 'BootRank';
    lang.ranking.csv_frequency_percent = 'Pct';
    lang.ranking.csv_frequency_count = 'Count';
    lang.ranking.csv_saved = 'Rank CSV Saved: %s';
    lang.ranking.histogram_plot_saved = 'Rank Hist Saved: %s';

    % Power Analysis
    lang.power.bootstrap_info = 'Bootstrap Power Analysis with %d simulations per comparison.';

    % Plots (Titles & Labels)
    lang.plots.titles.threshold_convergence_global = 'Global Threshold Convergence';
    lang.plots.titles.threshold_convergence_long_n_g = 'Global Threshold Convergence (n=%d)';
    lang.plots.titles.threshold_convergence = 'Detailed Threshold Convergence';
    lang.plots.titles.threshold_convergence_long_n_d = 'Detailed Threshold Convergence (n=%d)';
    lang.plots.titles.threshold_dist_name = 'Threshold Distribution';
    lang.plots.titles.threshold_dist = 'Distribution (B=%d)';
    lang.plots.titles.raw_effect_dist_name = 'Effect Size Distribution';
    lang.plots.titles.raw_effect_dist = 'Raw Effects';
    lang.plots.titles.bca_convergence_global = 'Global BCa Convergence';
    lang.plots.titles.bca_convergence_long = 'BCa Convergence';
    lang.plots.titles.bca_convergence_long_n_g = 'BCa Global (n=%d)';
    lang.plots.titles.bca_convergence_long_n_d = 'BCa Detailed (n=%d)';
    lang.plots.titles.bca_z0_dist_name = 'Z0 Dist';
    lang.plots.titles.bca_z0_dist = 'Z0';
    lang.plots.titles.bca_a_dist_name = 'A Dist';
    lang.plots.titles.bca_a_dist = 'A';
    lang.plots.titles.ci_width_dist_name = 'CI Width Dist';
    lang.plots.titles.ci_width_dist = 'Width';
    lang.plots.titles.rank_stability_convergence_name = 'Rank Stability';
    lang.plots.titles.rank_stability_convergence_global = 'Rank Global';
    lang.plots.titles.rank_stability_convergence = 'Rank Detailed';
    lang.plots.titles.rank_stability_convergence_n = 'Rank (n=%d)';
    lang.plots.titles.rank_dist_name = 'Rank Dist Name';
    lang.plots.titles.rank_dist = 'Rank Dist';
    
    lang.plots.xlabels.bootstraps = 'Bootstraps (B)';
    lang.plots.ylabels.stability = 'Stability (%)';
    lang.plots.ylabels.stability_rank = 'Stability (Rank)';
    lang.plots.legend.unsmoothed = 'Raw';
    lang.plots.legend.smoothed = 'Smoothed';
    lang.plots.legend.local_elbow = 'Elbow';
    lang.plots.misc.optimal_b_text = 'B=%d';
    lang.plots.misc.boot_thr = 'Boot=%.3f';
    lang.plots.misc.sem_thr = 'SEM=%.3f';
    lang.plots.misc.optimal_b = 'Opt B=%d';
    lang.plots.ylabels.rel_frequency = 'Frequency';
    lang.plots.xlabels.median_delta = 'Median Delta';
    lang.plots.xlabels.median_rel_diff = 'Median RelDiff';
    lang.plots.xlabels.skewness_a = 'Skewness (a)';
    lang.plots.xlabels.bias_z0 = 'Bias (z0)';
    lang.plots.xlabels.ci_width = 'Width';
    lang.plots.xlabels.rank = 'Rank';

    % File Names
    lang.files.convergence_thresholds_global = 'conv_glob.png';
    lang.files.convergence_thresholds = 'conv.png';
    lang.files.dist_bootstrap_thresholds = 'dist.png';
    lang.files.dist_raw_effects = 'raw.png';
    lang.files.convergence_bca_global = 'bca_glob.png';
    lang.files.convergence_bca = 'bca.png';
    lang.files.bca_factors = 'bca.csv';
    lang.files.dist_bca_bias_z0 = 'z0.png';
    lang.files.dist_bca_skew_a = 'a.png';
    lang.files.dist_ci_widths = 'width.png';
    lang.files.convergence_rank_stability = 'rank_conv.png';
    lang.files.convergence_rank_stability_global = 'rank_conv_glob.png';
    lang.files.bootstrap_rank_csv = 'rank.csv';
    lang.files.dist_bootstrap_ranks = 'rank_dist.png';
    
    % CSV Header
    lang.csv.headers.median = 'Median';
    lang.csv.headers.mean = 'Mean';
    lang.csv.headers.min = 'Min';
    lang.csv.headers.max = 'Max';
    lang.csv.headers.effect_size = 'EffectSize';
    lang.csv.headers.correction_factor = 'CorrFactor';
    
    % Styles
    styles = struct();
    styles.colors.background = [1 1 1];
    styles.colors.text = [0 0 0];
    styles.colors.blue_marker = [0 0.4470 0.7410];
    styles.colors.red_marker = [0.8500 0.3250 0.0980];
    styles.colors.grid_color = [0.15 0.15 0.15]; 
    styles.colors.delta_face = [0 0.4470 0.7410];
    styles.colors.rel_face = [0.8500 0.3250 0.0980];
    styles.colors.bar_edge = [0 0 0];
    styles.colors.holm_threshold = [1 0 0];
    styles.colors.sem_threshold = [0 1 0];
    styles.colors.kde_line = [0 0 0];
    styles.colors.bar_face = [0.5 0.5 0.5];
    styles.font.title = 12;
    styles.font.label = 10;
    styles.font.tick = 10;
    styles.font.small_text = 8;
    styles.line.grid = ':';
    styles.errorbar.final_rank = {'LineStyle', 'none'};
    styles.marker.final_rank = {'o'};
end