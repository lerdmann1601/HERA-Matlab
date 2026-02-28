function results = simulate(scenarios, params, n_sims_per_cond, refs, cfg_base, temp_dir, styles, lang, hWait, out_dir, ts_str, final_out_dir, colors, modes, limits, n_datasets)
% SIMULATE - Executes the parallel simulation core for the convergence robustness study.
%
% Syntax:
%   results = HERA.analysis.convergence.simulate(scen, params, n, refs, cfg, tmp, styles, lang, hWait, ...)
%
% Description:
%   This function manages the computationally intensive part of the analysis using an optimized
%   Interleaved Batching approach. It iterates over the defined scenarios and executes the 
%   validations in small batches to minimize memory footprint while maximizing worker utilization.
%
% Workflow:
%   1. Scenario Loop: 
%      Iterates through each configured scenario (Distribution/N).
%   2. Interleaved Batching: 
%      - Breaks the total simulations into small chunks (e.g., 4 sims/batch).
%      - For each batch, data and reference values are generated locally (Just-in-Time).
%   3. Parallel Testing (Worker Balancing): 
%      - Submits all test combinations (Thresholds, BCa, Ranking) for the batch to a shared parallel pool.
%      - Uses fetchNext to process results immediately as workers finish, reducing idle time.
%   4. Memory Management:
%      - Immediately clears simulation data and futures after each batch to prevent RAM growth.
%   5. Reporting:
%      - Updates progress incrementally and saves intermediate plots after each scenario.
%
% Inputs:
%   scenarios        - Struct array defining the data scenarios (N, Distribution).
%   params           - Struct containing the parameter sets for Thr, BCa, Rnk.
%   n_sims_per_cond  - Number of simulations per scenario.
%   refs             - Struct with reference B values.
%   cfg_base         - Base HERA configuration object.
%   temp_dir         - Directory for temporary file artifacts.
%   styles, lang     - HERA design/language structs.
%   hWait            - Handle to the waitbar for progress updates.
%   out_dir          - Base output directory (for potential reports).
%   ts_str           - Timestamp string for file naming.
%   final_out_dir    - Directory for final PDFs.
%   colors           - Color matrix for plots.
%   modes            - Cell array of mode names.
%   limits           - Struct with B limits.
%
% Outputs:
%   results          - Nested struct containing all simulation metrics (err, cost, fail) 
%                      organized by scenario and method.
%
% Author: Lukas von Erdmannsdorff

    %% Main Execution Loop
    num_modes = length(params.thr);
    
    % Memory configuration
    if isfield(cfg_base, 'system') && isfield(cfg_base.system, 'target_memory')
        TARGET_MEMORY = cfg_base.system.target_memory;
    else
        % Should not happen given config.m structure, but safe fallback
        TARGET_MEMORY = 200; 
    end
    
    if isfield(cfg_base, 'num_workers') && isnumeric(cfg_base.num_workers)
        num_workers = cfg_base.num_workers;
    else
        num_workers = feature('numcores');
    end
    effective_memory_mb = TARGET_MEMORY / max(1, num_workers);
    
    % Initialize Result Structure
    scenario_res = repmat(struct('name', '', 'N', 0, 'Dist', '', 'DataSummary', '', ...
                                 'thr', init_storage(n_sims_per_cond, num_modes), ...
                                 'bca', init_storage(n_sims_per_cond, num_modes), ...
                                 'rnk', init_storage(n_sims_per_cond, num_modes)), ...
                          length(scenarios), 1);
                      
    % Global Progress Tracking
    total_tests_global = length(scenarios) * n_sims_per_cond * 3 * num_modes;
    global_tests_completed = 0;

    % Loop over Scenarios
    for sc_idx = 1:length(scenarios)
        sc = scenarios(sc_idx);
        scenario_res(sc_idx).name = sc.name; 
        scenario_res(sc_idx).N = sc.N; 
        scenario_res(sc_idx).Dist = sc.Dist; 
        scenario_res(sc_idx).DataSummary = sc.DataSummary;
        
        fprintf('\n--- Scenario %d/%d: %s ---\n', sc_idx, length(scenarios), sc.name);
        t_scenario = tic;
        
        % Dynamic batch sizing
        n_pairs = nchoosek(n_datasets, 2);
        bytes_per_double = 8;
        bytes_per_sim = (sc.N * n_datasets + n_pairs*4 + n_datasets) * bytes_per_double;
        total_memory_needed = (n_sims_per_cond * bytes_per_sim) / (1024^2);
        
        if total_memory_needed <= effective_memory_mb
            sims_per_batch = n_sims_per_cond;
        else
            sims_per_batch = max(2, floor((effective_memory_mb * 1024^2) / bytes_per_sim));
        end
        sims_per_batch = min(sims_per_batch, num_workers);  % Match batch size to worker count
        
        fprintf('  [Config: sims_per_batch=%d for N=%d]\n', sims_per_batch, sc.N);
        
        %% Interleaved Batch Loop
        % Divide simulations into small batches to keep memory footprint low
        for batch_start = 1:sims_per_batch:n_sims_per_cond
            t_batch = tic;
            batch_end = min(batch_start + sims_per_batch - 1, n_sims_per_cond);
            batch_sims = batch_start:batch_end;
            num_in_batch = length(batch_sims);
            
            fprintf('  Batch %d-%d: Preparing Data... ', batch_start, batch_end);
            
            % 1. Generate Data & References (Parallelized / Coarse-Grained)
            sim_data_batch = cell(num_in_batch, 1);
            
            % Parallelize outer loop to saturate cores during "Preparing Data"
            % Inner functions (calculate_*) will detect single-batch and run serially
            % inside each worker, avoiding nested parallelism overhead.
            parfor i = 1:num_in_batch
                s_idx = batch_sims(i);
                
                % Create unique temp dir for this worker to prevent plot collisions
                % quiet_* functions write to graphics_dir/Threshold_Analysis etc.
                worker_temp_dir = fullfile(temp_dir, sprintf('sim_%d', s_idx));
                if ~exist(worker_temp_dir, 'dir'), mkdir(worker_temp_dir); end
                
                % Bit-Perfect Seeding Constraint
                sim_seed = 123 + (sc_idx-1)*10000 + s_idx;
                ref_seed = sim_seed + 1;
                
                % Generate Data
                dataStream = RandStream('mlfg6331_64', 'Seed', sim_seed);
                d_all = generate_data_vectorized(sc, dataStream, n_datasets);
                all_data = {d_all};
                p_idx_sim = nchoosek(1:n_datasets, 2);
                ds_names = cellstr("D" + (1:n_datasets));
                
                % Calc Real Effects (Fast)
                eff = HERA.test.TestHelper.calculate_real_effects(all_data, 1);
                
                % Reference Calculations (Using High Precision Refs from input)
                
                % Ref: Threshold
                refStream = RandStream('mlfg6331_64', 'Seed', ref_seed);
                [ref_d_t, ref_r_t, ~, ~, d_vals_all, rel_vals_all] = ...
                    quiet_thresholds(all_data, sc.N, cfg_base, worker_temp_dir, refs.thr, refStream, styles, lang);
                ref_thr_struct = struct('d_thresh', ref_d_t, 'rel_thresh', ref_r_t);
                
                % Ref: BCa
                refStream = RandStream('mlfg6331_64', 'Seed', ref_seed);
                [~, ref_ci_d] = quiet_bca_ci(all_data, d_vals_all, rel_vals_all, p_idx_sim, sc.N, ...
                    cfg_base, cfg_base.metric_names, worker_temp_dir, worker_temp_dir, refs.bca, refStream, styles, lang, 'Ref');
                
                % Ref: Ranking
                [~, base_rank] = HERA.calculate_ranking(all_data, eff, ref_thr_struct, cfg_base, ds_names, p_idx_sim);
                refStream = RandStream('mlfg6331_64', 'Seed', ref_seed);
                [boot_r_ref] = quiet_bootstrap_ranking(all_data, ref_thr_struct, cfg_base, ds_names, base_rank, p_idx_sim, sc.N, ...
                    worker_temp_dir, worker_temp_dir, refs.rnk, refStream, styles, lang, 'Ref');
                
                % Store Package (Data transfer back to client)
                sim_data_batch{i} = struct(...
                    'all_data', {all_data}, 'd_vals_all', d_vals_all, 'rel_vals_all', rel_vals_all, ...
                    'p_idx', p_idx_sim, 'ds_names', {ds_names}, 'base_rank', base_rank, ...
                    'ref_thr_struct', ref_thr_struct, 'ref_thr_d', ref_d_t(1), ...
                    'ref_bca_width', ref_ci_d(1,2) - ref_ci_d(1,1), ...
                    'ref_rnk_mean', mean(boot_r_ref(2,:)), 'ref_seed', ref_seed);
                    
                % Cleanup worker temp files
                if exist(worker_temp_dir, 'dir')
                    rmdir(worker_temp_dir, 's');
                end
            end
            fprintf('Done (%.2fs).\n', toc(t_batch));
            
            % 2. Parallel Processing (Worker Balancing)
            futures = parallel.FevalFuture.empty(0, 0);
            
            % Submit All Tasks (Interleaved Methods & Modes)
            % NOTE: Modes are submitted in REVERSE order (Strict → Default → Relaxed)
            % so that slower Strict tasks start first and don't block at batch boundaries.
            for i = 1:num_in_batch
                s_idx = batch_sims(i);
                sd = sim_data_batch{i};
                
                % Thresholds (Method ID = 1) - Strict first
                for m = num_modes:-1:1
                    pp = map_params(params.thr{m});
                    futures(end+1) = parfeval(@run_single_test, 6, ...
                        s_idx, m, 1, sd, pp, sc.N, cfg_base, temp_dir, styles, lang); 
                end
                
                % BCa (Method ID = 2) - Strict first
                for m = num_modes:-1:1
                    pp = map_params(params.bca{m});
                    futures(end+1) = parfeval(@run_single_test, 6, ...
                        s_idx, m, 2, sd, pp, sc.N, cfg_base, temp_dir, styles, lang);
                end
                
                % Ranking (Method ID = 3) - Strict first
                for m = num_modes:-1:1
                    pp = map_params(params.rnk{m});
                    futures(end+1) = parfeval(@run_single_test, 6, ...
                        s_idx, m, 3, sd, pp, sc.N, cfg_base, temp_dir, styles, lang);
                end
            end
            
            % Collect Results (Asynchronous / Load Balanced)
            tests_in_batch = numel(futures);
            tests_done_batch = 0;
            
            % Use loop to fetch exactly as many results as submitted
            for f_k = 1:tests_in_batch
                [~, ret_s_idx, ret_m_idx, ret_method_id, ret_val, ret_cost, ret_fail] = fetchNext(futures);
                
                % Assign to Storage
                switch ret_method_id
                    case 1 % Thr
                        scenario_res(sc_idx).thr.err(ret_s_idx, ret_m_idx) = ret_val;
                        scenario_res(sc_idx).thr.cost(ret_s_idx, ret_m_idx) = ret_cost;
                        scenario_res(sc_idx).thr.fail(ret_s_idx, ret_m_idx) = ret_fail;
                    case 2 % BCa
                        scenario_res(sc_idx).bca.err(ret_s_idx, ret_m_idx) = ret_val;
                        scenario_res(sc_idx).bca.cost(ret_s_idx, ret_m_idx) = ret_cost;
                        scenario_res(sc_idx).bca.fail(ret_s_idx, ret_m_idx) = ret_fail;
                    case 3 % Rnk
                        scenario_res(sc_idx).rnk.err(ret_s_idx, ret_m_idx) = ret_val;
                        scenario_res(sc_idx).rnk.cost(ret_s_idx, ret_m_idx) = ret_cost;
                        scenario_res(sc_idx).rnk.fail(ret_s_idx, ret_m_idx) = ret_fail;
                end
                
                tests_done_batch = tests_done_batch + 1;
                global_tests_completed = global_tests_completed + 1;
                
                % Calculate sub-progress metrics
                sims_done_in_scenario = batch_start - 1 + ceil(tests_done_batch / (3 * num_modes));
                batch_pct = (tests_done_batch / tests_in_batch) * 100;
                
                % Global progress: cap at 99% until truly complete to avoid premature "done" display
                global_pct = global_tests_completed / total_tests_global;
                if global_pct >= 1 && tests_done_batch < tests_in_batch
                    global_pct = 0.99;  % Still processing last batch
                end
                
                % Incremental Waitbar with enhanced info
                if isvalid(hWait)
                     tests_per_sim = 3 * num_modes;  % 9 tests per simulation
                     total_tests_in_scenario = n_sims_per_cond * tests_per_sim;
                     tests_done_in_scenario = (batch_start - 1) * tests_per_sim + tests_done_batch;
                     
                     waitbar(global_pct, hWait, ...
                         sprintf('Scenario %d/%d (%s)\nBatch %d-%d | Sims %d/%d | Test %d/%d | %.0f%%', ...
                         sc_idx, length(scenarios), sc.name, batch_start, batch_end, ...
                         sims_done_in_scenario, n_sims_per_cond, ...
                         tests_done_in_scenario, total_tests_in_scenario, ...
                         global_pct * 100));
                end
            end
            
            % 3. Cleanup & Report
            % Close graphics (excluding waitbar)
            all_figs = findall(0, 'Type', 'figure');
            if isvalid(hWait)
                all_figs(all_figs == hWait) = [];
            end
            close(all_figs);
            clear futures sim_data_batch; % Then release memory
            
            dt = toc(t_batch);
            fprintf('    Batch finished in %.1fs (%.2fs/sim) [Progress: %.0f%%]\n', dt, dt/num_in_batch, global_pct*100);
        end
        
        %% End of Scenario
        t_sc_elapsed = toc(t_scenario);
        fprintf('  Scenario completed in %.1fs (%.1fs/sim)\n', t_sc_elapsed, t_sc_elapsed / n_sims_per_cond);
        
        % Incremental Plotting
        try
            HERA.analysis.convergence.plot('scientific_reports', [scenario_res(sc_idx)], modes, styles, refs, limits, params, final_out_dir, 'Incremental', ts_str);
            
            % Incremental CSV Saving
            % Ensure CSV directory matches what is expected in save_csv
            csv_dir = fullfile(final_out_dir, 'CSV');
            HERA.analysis.convergence.save_csv([scenario_res(sc_idx)], modes, csv_dir, ts_str);
        catch ME
            fprintf('Warning: Incremental reporting/saving failed: %s\n', ME.message);
        end
    end
    
    results = scenario_res;
end

%% Local Worker Function (Executed by parfeval)

function [s_idx, m_idx, method_id, res_val, res_cost, res_fail] = run_single_test(s_idx, m_idx, method_id, sd, param, N, cfg, tmp, styles, lang)
    % Runs a single convergence test and returns ID + Results
    % method_id: 1=Thr, 2=BCa, 3=Rnk
    
    stream = RandStream('mlfg6331_64', 'Seed', sd.ref_seed);
    res_val = 0; res_cost = 0; res_fail = false;
    
    try
        if method_id == 1
            % Thresholds
            c = cfg; c.bootstrap_thresholds = param;
            [d_t, ~, ~, ~, ~, ~, ~, conv_B, stab_data] = ...
                quiet_thresholds(sd.all_data, N, c, tmp, [], stream, styles, lang);
            
            res_val = (d_t(1) - sd.ref_thr_d) / sd.ref_thr_d * 100;
            res_cost = conv_B;
            res_fail = ~stab_data.converged;
            
        elseif method_id == 2
            % BCa
            c = cfg; c.bootstrap_ci = param;
            [conv_B, ci_d, ~, ~, ~, ~, ~, stab_data] = ...
                quiet_bca_ci(sd.all_data, sd.d_vals_all, sd.rel_vals_all, sd.p_idx, N, ...
                c, c.metric_names, tmp, tmp, [], stream, styles, lang, 'Sim');
            
            width = ci_d(1, 2) - ci_d(1, 1);
            res_val = (width - sd.ref_bca_width) / sd.ref_bca_width * 100;
            res_cost = conv_B;
            res_fail = ~stab_data.converged;
            
        elseif method_id == 3
            % Ranking
            c = cfg; c.bootstrap_ranks = param;
            [boot_r, conv_B, stab_data] = ...
                quiet_bootstrap_ranking(sd.all_data, sd.ref_thr_struct, c, sd.ds_names, sd.base_rank, sd.p_idx, N, ...
                tmp, tmp, [], stream, styles, lang, 'Sim');
            
            mean_r = mean(boot_r(2, :));
            res_val = (mean_r - sd.ref_rnk_mean) / sd.ref_rnk_mean * 100;
            res_cost = conv_B;
            res_fail = ~stab_data.converged;
        end
    catch ME
        % In case of worker failure, return NaNs/true so we don't crash main loop
        fprintf('Error in worker (S:%d M:%d ID:%d): %s\n', s_idx, m_idx, method_id, ME.message);
        res_val = NaN; res_cost = 0; res_fail = true;
    end
end

%% Helper Functions

function so = map_params(pi)
    % Maps study params to HERA config
    so.n_trials = pi.n; so.smoothing_window = pi.sm;
    so.convergence_streak_needed = pi.st; so.convergence_tolerance = pi.tol;
    so.B_start = pi.start; so.B_step = pi.step; so.B_end = pi.end;
    so.min_steps_for_convergence_check = 1;
end

function s = init_storage(n, num_modes)
    if nargin < 2, num_modes = 3; end
    s.err = zeros(n, num_modes); s.cost = zeros(n, num_modes); s.fail = false(n, num_modes); 
end

function d_all = generate_data_vectorized(sc, stream, n_datasets)
    % Generates synthetic data based on the scenario configuration
    d_all = zeros(sc.N, n_datasets);
    switch sc.Dist
        case 'Normal'
            means = 10 + (0:n_datasets-1);
            r = randn(stream, sc.N, n_datasets);
            d_all = means + 2.0 * r;
        case 'LogNormal'
            means = 2.0 + (0:n_datasets-1)*0.1;
            r = randn(stream, sc.N, n_datasets);
            d_all = exp(means + 0.4 * r);
        case 'Likert'
            centers = linspace(3, 5, n_datasets);
            r = randn(stream, sc.N, n_datasets);
            raw = centers + 1.5 * r;
            d_all = min(7, max(1, round(raw)));
        case 'Bimodal'
            for k = 1:n_datasets
                is_mode2 = rand(stream, sc.N, 1) > 0.5;
                vals = zeros(sc.N, 1);
                n1 = sum(~is_mode2); n2 = sum(is_mode2);
                vals(~is_mode2) = 10 + randn(stream, n1, 1);
                vals(is_mode2)  = 15 + randn(stream, n2, 1);
                % Adjusted to 1.5 (approx d=0.55) to test robustness, not power
                d_all(:, k) = vals + (k-1)*1.5;
            end
        case 'NormalLarge'
            % Large Effect: Cohen's d = 1.0 (mean diff 2.0 / SD 2.0)
            means = 10 + (0:n_datasets-1) * 2.0;
            r = randn(stream, sc.N, n_datasets);
            d_all = means + 2.0 * r;
    end
end

%% Quiet Wrappers

function [d_t, r_t, rel_thresh_b, min_rel_dynamic, d_vals_all, rel_vals_all, pair_idx_all, conv_B, stab_data] = ...
    quiet_thresholds(all_data, num_probanden, config, graphics_dir, manual_B, s, styles, lang)
    cmd = ['[d_t, r_t, rel_thresh_b, min_rel_dynamic, d_vals_all, rel_vals_all, pair_idx_all, conv_B, stab_data, ' ...
           'h1, h2, h3, h4] = ' ...
           'HERA.calculate_thresholds(all_data, num_probanden, config, graphics_dir, manual_B, s, styles, lang);'];
    evalc(cmd);
    to_close = gobjects(0);
    if exist('h1', 'var'), to_close = [to_close; h1(:)]; end
    if exist('h2', 'var'), to_close = [to_close; h2(:)]; end
    if exist('h3', 'var'), to_close = [to_close; h3(:)]; end
    if exist('h4', 'var'), to_close = [to_close; h4(:)]; end
    if ~isempty(to_close), close(to_close(isgraphics(to_close))); end
end

function [conv_B, ci_d, ci_r, z0_d, a_d, z0_r, a_r, stab_data] = ...
    quiet_bca_ci(all_data, d_vals_all, rel_vals_all, pair_idx_all, num_probanden, ...
                 config, metric_names, graphics_dir, csv_dir, manual_B, s, styles, lang, base_name)
    cmd = ['[conv_B, ci_d, ci_r, z0_d, a_d, z0_r, a_r, stab_data, ' ...
           'h1, h2, h3, h4, h5] = ' ...
           'HERA.calculate_bca_ci(all_data, d_vals_all, rel_vals_all, pair_idx_all, num_probanden, ' ...
           'config, metric_names, graphics_dir, csv_dir, manual_B, s, styles, lang, base_name);'];
    evalc(cmd);
    to_close = gobjects(0);
    if exist('h1', 'var'), to_close = [to_close; h1(:)]; end
    if exist('h2', 'var'), to_close = [to_close; h2(:)]; end
    if exist('h3', 'var'), to_close = [to_close; h3(:)]; end
    if exist('h4', 'var'), to_close = [to_close; h4(:)]; end
    if exist('h5', 'var'), to_close = [to_close; h5(:)]; end
    if ~isempty(to_close), close(to_close(isgraphics(to_close))); end
end

function [boot_r, conv_B, stab_data] = ...
    quiet_bootstrap_ranking(all_data, thresholds, config, dataset_names, final_rank, pair_idx_all, ...
                            num_probanden, graphics_dir, csv_dir, manual_B, s, styles, lang, base_name)
    cmd = ['[boot_r, conv_B, stab_data, h1, h2] = ' ...
           'HERA.bootstrap_ranking(all_data, thresholds, config, dataset_names, final_rank, pair_idx_all, ' ...
           'num_probanden, graphics_dir, csv_dir, manual_B, s, styles, lang, base_name);'];
    evalc(cmd);
    to_close = gobjects(0);
    if exist('h1', 'var'), to_close = [to_close; h1(:)]; end
    if exist('h2', 'var'), to_close = [to_close; h2(:)]; end
    if ~isempty(to_close), close(to_close(isgraphics(to_close))); end
end
