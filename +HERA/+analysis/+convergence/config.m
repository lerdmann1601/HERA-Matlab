function [n_datasets, modes, scenarios, params, refs, limits, cfg_base, colors] = config(n_sims_per_cond)
% CONFIG - Central configuration for the convergence robustness study.
%
% Syntax:
%   [n_ds, modes, sc, params, refs, lim, cfg, col] = HERA.analysis.convergence.config(n_sims_per_cond)
%
% Description:
%   This function defines the constant parameters, data scenarios, and method configurations
%   used throughout the robustness study. It ensures that all components (simulate, plot)
%   share the exact same definitions.
%
% Inputs:
%   n_sims_per_cond - Not directly used here, but part of the signature for potential dynamic adjustments.
%
% Outputs:
%   n_datasets  - Number of datasets in the scenarios.
%   modes       - Cell array of mode names (e.g., {'Relaxed', 'Default', 'Strict'}).
%   scenarios   - Struct array defining each data scenario (name, N, Dist).
%   params      - Struct containing specific parameter sets for Thr, BCa, Rnk.
%   refs        - Struct defining the high-precision B values for reference calculations.
%   limits      - Struct defining the maximum B limits for the methods.
%   cfg_base    - The standard HERA configuration object initialized with defaults.
%   colors      - RGB matrix for consistent plotting colors.
%
% Author: Lukas von Erdmannsdorff

    %% 1. General Settings
    n_datasets = 6; % Number of datasets (User requested 6 for better scientific rigor)
    modes = {'Relaxed', 'Default', 'Strict'};
    
    %% 2. Data Scenarios
    scenarios = struct();
    scenarios(1).name = 'N = 25 (Normal)';   scenarios(1).N = 25; scenarios(1).Dist = 'Normal'; 
    scenarios(1).DataSummary = sprintf('%d Datasets, Means 10-%d, SD = 2.0', n_datasets, 10 + n_datasets - 1);
    scenarios(2).name = 'N = 50 (Normal)';   scenarios(2).N = 50; scenarios(2).Dist = 'Normal';
    scenarios(2).DataSummary = sprintf('%d Datasets, Means 10-%d, SD = 2.0', n_datasets, 10 + n_datasets - 1);
    scenarios(3).name = 'N = 100 (Normal)';  scenarios(3).N = 100; scenarios(3).Dist = 'Normal';
    scenarios(3).DataSummary = sprintf('%d Datasets, Means 10-%d, SD = 2.0', n_datasets, 10 + n_datasets - 1);
    scenarios(4).name = 'N = 50 (Skewed)';   scenarios(4).N = 50; scenarios(4).Dist = 'LogNormal';
    scenarios(4).DataSummary = sprintf('%d Datasets, Means 2.0-%.1f (Log), SD = 0.4', n_datasets, 2.0 + (n_datasets-1)*0.1);
    scenarios(5).name = 'N = 50 (Likert)';   scenarios(5).N = 50; scenarios(5).Dist = 'Likert';
    scenarios(5).DataSummary = sprintf('%d Datasets, Scale 1-7, Means 3-5, SD = 1.5', n_datasets);
    scenarios(6).name = 'N = 50 (Bimodal)';  scenarios(6).N = 50; scenarios(6).Dist = 'Bimodal';
    scenarios(6).DataSummary = sprintf('%d Datasets, Mix Means 10 & 15, SD = 2.7', n_datasets);
    scenarios(7).name = 'N = 50 (Large Effect)';  scenarios(7).N = 50; scenarios(7).Dist = 'NormalLarge';
    scenarios(7).DataSummary = sprintf('%d Datasets, Means 10-%.1f, SD = 2.0', n_datasets, 10 + (n_datasets - 1) * 2.0);
    
    %% 3. Parameter Sets
    % Thresholds
    p_thr{1} = struct('n', 15, 'sm', 2, 'st', 2, 'tol', 0.01, 'start', 100, 'step', 100, 'end', 10000); 
    p_thr{2} = struct('n', 25, 'sm', 3, 'st', 3, 'tol', 0.01, 'start', 100, 'step', 100, 'end', 10000); 
    p_thr{3} = struct('n', 35, 'sm', 4, 'st', 4, 'tol', 0.01, 'start', 100, 'step', 100, 'end', 10000); 
    
    % BCa Confidence Intervals
    p_bca{1} = struct('n', 20, 'sm', 2, 'st', 2, 'tol', 0.03, 'start', 100, 'step', 200, 'end', 20000); 
    p_bca{2} = struct('n', 30, 'sm', 3, 'st', 3, 'tol', 0.03, 'start', 100, 'step', 200, 'end', 20000);
    p_bca{3} = struct('n', 40, 'sm', 4, 'st', 4, 'tol', 0.03, 'start', 100, 'step', 200, 'end', 20000); 
    
    % Ranking Stability
    p_rank{1} = struct('n', 10, 'sm', 2, 'st', 2, 'tol', 0.005, 'start', 50, 'step', 25, 'end', 2000); 
    p_rank{2} = struct('n', 15, 'sm', 3, 'st', 3, 'tol', 0.005, 'start', 50, 'step', 25, 'end', 2000);
    p_rank{3} = struct('n', 20, 'sm', 4, 'st', 4, 'tol', 0.005, 'start', 50, 'step', 25, 'end', 2000); 
    
    params.thr = p_thr; params.bca = p_bca; params.rnk = p_rank;

    % Reference settings (High B values for "Truth")
    ref_B_thr = 15000; ref_B_bca = 30000; ref_B_rnk = 5000;
    refs.thr = ref_B_thr; refs.bca = ref_B_bca; refs.rnk = ref_B_rnk;
    
    % Upper limits for visualization
    limits.thr = p_thr{1}.end; limits.bca = p_bca{1}.end; limits.rnk = p_rank{1}.end;
    
    %% 4. Base Framework Config
    cfg_base = HERA.default();
    cfg_base.create_reports = false;
    cfg_base.metric_names = {'SimMetric'};
    cfg_base.ranking_mode = 'M1';
    cfg_base.timestamp = 'RobustnessStudy';
    % Dynamic Target Memory Calculation (Inline)
    try
        ram_gb = NaN;
        if ispc
            [~, sysview] = memory;
            ram_gb = sysview.PhysicalMemory.Total / (1024^3);
        elseif ismac
            [stat, out] = system('sysctl hw.memsize');
            if stat == 0
                tok = regexp(out, '\d+', 'match');
                if ~isempty(tok), ram_gb = str2double(tok{1}) / (1024^3); end
            end
        elseif isunix
            [stat, out] = system('grep MemTotal /proc/meminfo');
            if stat == 0
                tok = regexp(out, '\d+', 'match');
                if ~isempty(tok), ram_gb = str2double(tok{1}) / (1024^2); end  % kB to GB
            end
        end
        
        if isnan(ram_gb)
             target_mem = 200; % Fallback
             fprintf('Warning: RAM detection failed. Using fallback: %d MB.\n', target_mem);
        else
             % Formula: 12.5 MB per GB (Conservative for parfeval overhead in this module)
             target_mem = max(200, round(ram_gb * 12.5));
             fprintf('System RAM: %.1f GB. Dynamic target_memory: %d MB.\n', ram_gb, target_mem);
        end
    catch
        target_mem = 200;
        fprintf('Error in RAM detection. Using fallback: %d MB.\n', target_mem);
    end
    cfg_base.system.target_memory = target_mem;
    
    colors = [0.8 0.3 0.3; 0.2 0.5 0.8; 0.3 0.7 0.4]; 
end
