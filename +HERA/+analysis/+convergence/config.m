function [N, modes, scenarios, params, refs, limits, cfg_base, colors, ram_gb] = config(n_sims_per_cond, customConfig)
% CONFIG - Central configuration for the convergence robustness study.
%
% Syntax:
%   [N, modes, sc, params, refs, lim, cfg, col, ram_gb] = HERA.analysis.convergence.config(n_sims_per_cond)
%
% Description:
%   This function defines the parameters and configuration for the convergence robustness study. 
%   It ensures that all components (simulate, plot) share the exact same definitions. 
%   Supports dynamic scaling of statistical parameters based on the number of candidates (N).
%
% Inputs:
%   n_sims_per_cond - Not directly used here, but part of the signature for potential dynamic adjustments.
%   customConfig    - (Optional) Struct loaded from JSON for overriding internal arrays safely.
%                     Supports fields such as:
%                       - `N` (Candidate count, should be between 3 and 15)
%                       - `selected_scenarios` (Array of indices to run, e.g. [1, 3, 5])
%                       - `num_scenarios` (Shortcut: Number of first scenarios to run, 1-8)
%                       - `scenarios` (Array of overrides for n, Step, SD, Base)
%                       - `target_memory`
%                       - `bootstrap_seed_offset`
%                       - `scenario_seed_offset`
%                       - `reference_seed_offset`
%                       - `reference_step_offset`
%                       - `modes` (with sub-structs for 'thr', 'bca', 'rnk')
%
% Outputs:
%   N           - Number of datasets/candidates in the scenarios.
%   modes       - Cell array of mode names (e.g., {'Relaxed', 'Default', 'Strict'}).
%   scenarios   - Struct array defining each data scenario (name, n, Dist).
%   params      - Struct containing specific parameter sets for Thr, BCa, Rnk.
%   refs        - Struct defining the high-precision B values for reference calculations.
%   limits      - Struct defining the maximum B limits for the methods.
%   cfg_base    - The standard HERA configuration object initialized with defaults.
%   colors      - RGB matrix for consistent plotting colors.
%   ram_gb      - Detected system RAM in GB.
%
% Author: Lukas von Erdmannsdorff

    if nargin < 2
        customConfig = struct();
    end

    %% 1. General Settings
    N = 6; % Number of datasets/candidates (Default)
    if isfield(customConfig, 'N') && isnumeric(customConfig.N)
        N = round(customConfig.N);
        if N < 3 && N >= 1
            fprintf('Warning: N = %d is below the recommended minimum of 3. Results may be less representative.\n', N);
        elseif N > 15
            fprintf('Warning: N = %d exceeds the recommended maximum of 15.\n', N);
            fprintf('         This will significantly increase computation time and memory usage.\n');
            fprintf('         Note: The Holm-Bonferroni correction becomes extremely strict at this level.\n');
            fprintf('         Please adjust your parameters accordingly.\n');
        end
    end
    % Hard constraint: N must be at least 2 for pairs calculation (nchoosek)
    if N < 2
        error('HERA:Analysis:InvalidConfig', 'N must be at least 2 (minimum one pair required). Provided N = %d.', N);
    end
    
    modes = {'Relaxed', 'Default', 'Strict'};
    
    %% 2. Data Scenarios
    % Define the 8 core scenarios with their scaling logic defaults.
    % Base: Starting mean/offset, Step: Gap between means, SD: Noise.
    sc_defs = [
        struct('name', '', 'n', 25,  'Dist', 'Normal',       'Base', 10.0, 'Step', 1.0, 'End', NaN,  'Base2', NaN,  'SD', 2.0); % Default: d = 0.5
        struct('name', '', 'n', 50,  'Dist', 'Normal',       'Base', 10.0, 'Step', 1.0, 'End', NaN,  'Base2', NaN,  'SD', 2.0); % Default: d = 0.5
        struct('name', '', 'n', 100, 'Dist', 'Normal',       'Base', 10.0, 'Step', 1.0, 'End', NaN,  'Base2', NaN,  'SD', 2.0); % Default: d = 0.5
        struct('name', '', 'n', 50,  'Dist', 'Skewed',       'Base', 2.0,  'Step', 0.1, 'End', NaN,  'Base2', NaN,  'SD', 0.4);
        struct('name', '', 'n', 50,  'Dist', 'Likert',       'Base', 3.0,  'Step', NaN, 'End', 5.0,  'Base2', NaN,  'SD', 1.5);
        struct('name', '', 'n', 50,  'Dist', 'Bimodal',      'Base', 10.0, 'Step', 1.5, 'End', NaN,  'Base2', 15.0, 'SD', 1.0);
        struct('name', '', 'n', 50,  'Dist', 'Small Effect', 'Base', 10.0, 'Step', 0.4, 'End', NaN,  'Base2', NaN,  'SD', 2.0); % Small: d = 0.2
        struct('name', '', 'n', 50,  'Dist', 'Large Effect', 'Base', 10.0, 'Step', 2.0, 'End', NaN,  'Base2', NaN,  'SD', 2.0)  % Large: d = 1.0
    ];
    
    % Parse JSON Overrides for Scenarios
    if isfield(customConfig, 'scenarios') && isstruct(customConfig.scenarios)
        for i = 1:min(length(sc_defs), length(customConfig.scenarios))
            ov = customConfig.scenarios(i);
            % Loop through fields to allow partial overrides
            f_ov = fieldnames(ov);
            for f = 1:length(f_ov)
                if isfield(sc_defs(i), f_ov{f})
                    sc_defs(i).(f_ov{f}) = ov.(f_ov{f});
                end
            end
        end
    end

    % Scenario Selection Logic ID
    % Default: Run all defined scenarios
    selected_idx = 1:length(sc_defs);
   
    % Selected Scenarios
    if isfield(customConfig, 'selected_scenarios') && isnumeric(customConfig.selected_scenarios)
        selected_idx = customConfig.selected_scenarios;
        % Validate all indices
        if any(selected_idx < 1 | selected_idx > length(sc_defs))
            error('HERA:Analysis:InvalidConfig', ...
                'selected_scenarios contains invalid indices. Must be between 1 and %d.', length(sc_defs));
        end
    elseif isfield(customConfig, 'num_scenarios') && isnumeric(customConfig.num_scenarios)
        num_target = round(customConfig.num_scenarios);
        if num_target >= 1 && num_target <= length(sc_defs)
            selected_idx = 1:num_target;
        else
             error('HERA:Analysis:InvalidConfig', ...
                 'num_scenarios must be between 1 and %d. Provided: %d.', length(sc_defs), num_target);
        end
    end
    
    % Apply Selection
    sc_defs = sc_defs(selected_idx);

    scenarios = struct();
    for i = 1:length(sc_defs)
        def = sc_defs(i);
        
        % Robustness Check for n (Sample size)
        if def.n < 20
            fprintf('Warning (Scenario %d): n = %d is quite small for bootstrap-based accuracy studies.\n', i, def.n);
        end
        if def.n < 5
            fprintf('Error (Scenario %d): n must be at least 5. Setting n = 5.\n', i);
            def.n = 5;
        end

        % Generate Dynamic Name if not overridden (simple n = [val])
        if isempty(def.name)
            scenarios(i).name = sprintf('n = %d', def.n);
        else
            scenarios(i).name = def.name;
        end
        
        scenarios(i).n = def.n;
        scenarios(i).Dist = def.Dist;
        
        % Store scaling parameters and metadata for simulate.m/reports
        scenarios(i).Base  = def.Base;
        scenarios(i).Step  = def.Step;
        scenarios(i).End   = def.End;
        scenarios(i).Base2 = def.Base2;
        scenarios(i).SD    = def.SD;
        
        % Generate dynamic DataSummary
        switch def.Dist
            case {'Normal', 'Large Effect', 'Small Effect'}
                scenarios(i).DataSummary = sprintf('%d Datasets, Means %.1f-%.1f, SD = %.1f', ...
                    N, def.Base, def.Base + (N-1)*def.Step, def.SD);
            case 'Skewed'
                scenarios(i).DataSummary = sprintf('%d Datasets, Means %.1f-%.1f (Log), SD = %.1f (Log)', ...
                    N, def.Base, def.Base + (N-1)*def.Step, def.SD);
            case 'Likert'
                scenarios(i).DataSummary = sprintf('%d Datasets, Scale 1-7, Means %.1f-%.1f, SD = %.1f', ...
                    N, def.Base, def.End, def.SD);
            case 'Bimodal'
                % For Bimodal, we show both means and the internal SD
                scenarios(i).DataSummary = sprintf('%d Datasets, Mix Means %.1f & %.1f, Mode SD = %.1f', ...
                    N, def.Base, def.Base2, def.SD);
        end
    end
    
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
    p_rank{1} = struct('n', 10, 'sm', 2, 'st', 2, 'tol', 0.005, 'start', 50, 'step', 25, 'end', 2500); 
    p_rank{2} = struct('n', 15, 'sm', 3, 'st', 3, 'tol', 0.005, 'start', 50, 'step', 25, 'end', 2500);
    p_rank{3} = struct('n', 20, 'sm', 4, 'st', 4, 'tol', 0.005, 'start', 50, 'step', 25, 'end', 2500); 
    
    params.thr = p_thr; params.bca = p_bca; params.rnk = p_rank;

    % Reference settings (High B values for "Truth")
    ref_B_thr = 25000; ref_B_bca = 50000; ref_B_rnk = 10000;
    
    % Override Reference settings
    if isfield(customConfig, 'refs') && isstruct(customConfig.refs)
        if isfield(customConfig.refs, 'thr'), ref_B_thr = customConfig.refs.thr; end
        if isfield(customConfig.refs, 'bca'), ref_B_bca = customConfig.refs.bca; end
        if isfield(customConfig.refs, 'rnk'), ref_B_rnk = customConfig.refs.rnk; end
    end
    refs.thr = ref_B_thr; refs.bca = ref_B_bca; refs.rnk = ref_B_rnk;
    
    % Upper limits for visualization
    limits.thr = p_thr{1}.end; limits.bca = p_bca{1}.end; limits.rnk = p_rank{1}.end;
    
    %% 4. Base Framework Config
    cfg_base = HERA.default();
    cfg_base.create_reports = false;
    cfg_base.metric_names = {'SimMetric'};
    cfg_base.ranking_mode = 'M1';
    cfg_base.timestamp = 'RobustnessStudy';
    cfg_base.quiet_mode = true;
    cfg_base.create_csvs = false;
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
             % Formula: 25 MB per GB 
             % (Conservative for parfeval overhead in this module but should be safe since we don't store all data and no graphics are generated)
             target_mem = max(200, round(ram_gb * 25));
             fprintf('System RAM: %.1f GB. Dynamic target_memory: %d MB.\n', ram_gb, target_mem);
        end
    catch
        target_mem = 200;
        fprintf('Error in RAM detection. Using fallback: %d MB.\n', target_mem);
    end
    if isfield(customConfig, 'target_memory') && isnumeric(customConfig.target_memory)
        target_mem = customConfig.target_memory;
        fprintf('Overriding memory via Config: %d MB.\n', target_mem);
    end
    cfg_base.system.target_memory = target_mem;
    
    cfg_base.simulation_seed = 123;
    if isfield(customConfig, 'simulation_seed') && isnumeric(customConfig.simulation_seed)
        cfg_base.simulation_seed = customConfig.simulation_seed;
    end
    
    if isfield(customConfig, 'bootstrap_seed_offset') && isnumeric(customConfig.bootstrap_seed_offset)
        cfg_base.bootstrap_seed_offset = customConfig.bootstrap_seed_offset;
    end
    
    cfg_base.scenario_seed_offset = 10000;
    if isfield(customConfig, 'scenario_seed_offset') && isnumeric(customConfig.scenario_seed_offset)
        cfg_base.scenario_seed_offset = customConfig.scenario_seed_offset;
    end
    
    cfg_base.reference_seed_offset = 5000;
    if isfield(customConfig, 'reference_seed_offset') && isnumeric(customConfig.reference_seed_offset)
        cfg_base.reference_seed_offset = customConfig.reference_seed_offset;
    end
    
    cfg_base.reference_step_offset = 1000; % Default gap between ref components (Thr, BCa, Rnk)
    if isfield(customConfig, 'reference_step_offset') && isnumeric(customConfig.reference_step_offset)
        cfg_base.reference_step_offset = customConfig.reference_step_offset;
    end
    
    % Override Params using customConfig
    if isfield(customConfig, 'modes') && isstruct(customConfig.modes)
        fields = {'Relaxed', 'Default', 'Strict'};
        for i=1:length(fields)
            f = fields{i};
            if isfield(customConfig.modes, f) && isstruct(customConfig.modes.(f))
                cu = customConfig.modes.(f);
                
                % Parse 'thr' sub-struct if provided
                if isfield(cu, 'thr') && isstruct(cu.thr)
                    if isfield(cu.thr, 'n'), p_thr{i}.n = cu.thr.n; end
                    if isfield(cu.thr, 'sm'), p_thr{i}.sm = cu.thr.sm; end
                    if isfield(cu.thr, 'st'), p_thr{i}.st = cu.thr.st; end
                    if isfield(cu.thr, 'tol'), p_thr{i}.tol = cu.thr.tol; end
                    if isfield(cu.thr, 'start'), p_thr{i}.start = cu.thr.start; end
                    if isfield(cu.thr, 'step'), p_thr{i}.step = cu.thr.step; end
                    if isfield(cu.thr, 'end'), p_thr{i}.end = cu.thr.end; end
                end
                
                % Parse 'bca' sub-struct if provided
                if isfield(cu, 'bca') && isstruct(cu.bca)
                    if isfield(cu.bca, 'n'), p_bca{i}.n = cu.bca.n; end
                    if isfield(cu.bca, 'sm'), p_bca{i}.sm = cu.bca.sm; end
                    if isfield(cu.bca, 'st'), p_bca{i}.st = cu.bca.st; end
                    if isfield(cu.bca, 'tol'), p_bca{i}.tol = cu.bca.tol; end
                    if isfield(cu.bca, 'start'), p_bca{i}.start = cu.bca.start; end
                    if isfield(cu.bca, 'step'), p_bca{i}.step = cu.bca.step; end
                    if isfield(cu.bca, 'end'), p_bca{i}.end = cu.bca.end; end
                end
                
                % Parse 'rnk' sub-struct if provided
                if isfield(cu, 'rnk') && isstruct(cu.rnk)
                    if isfield(cu.rnk, 'n'), p_rank{i}.n = cu.rnk.n; end
                    if isfield(cu.rnk, 'sm'), p_rank{i}.sm = cu.rnk.sm; end
                    if isfield(cu.rnk, 'st'), p_rank{i}.st = cu.rnk.st; end
                    if isfield(cu.rnk, 'tol'), p_rank{i}.tol = cu.rnk.tol; end
                    if isfield(cu.rnk, 'start'), p_rank{i}.start = cu.rnk.start; end
                    if isfield(cu.rnk, 'step'), p_rank{i}.step = cu.rnk.step; end
                    if isfield(cu.rnk, 'end'), p_rank{i}.end = cu.rnk.end; end
                end
            end
        end
        % Re-assign to params
        params.thr = p_thr; params.bca = p_bca; params.rnk = p_rank;
    end
    
    % Upper limits for visualization (re-evaluating in case it changed)
    limits.thr = p_thr{1}.end; limits.bca = p_bca{1}.end; limits.rnk = p_rank{1}.end;
    
    colors = [0.8 0.3 0.3; 0.2 0.5 0.8; 0.3 0.7 0.4]; 
end
