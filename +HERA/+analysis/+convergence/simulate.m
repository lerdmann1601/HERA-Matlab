function results = simulate(scenarios, params, n_sims_per_cond, refs, cfg_base, temp_dir, styles, lang, hWait, out_dir, ts_str, final_out_dir, colors, modes, limits, N)
% SIMULATE - Executes the parallel simulation core for the convergence robustness study.
%
% Syntax:
%   results = HERA.analysis.convergence.simulate(scen, params, n_sims, refs, cfg, tmp, styles, lang, hWait, out_dir, ts, final_out, colors, modes, limits, N)
%
% Description:
%   This function manages the computationally intensive part of the analysis using an optimized
%   Interleaved Batching approach. It iterates over the defined scenarios and executes the 
%   validations in small batches to minimize memory footprint while maximizing worker utilization.
%
% Workflow:
%   1. Scenario Loop: 
%      Iterates through each configured scenario (Distribution/n).
%   2. Interleaved Batching (Hybrid Parallelization): 
%      - Breaks the total simulations into small chunks (matches worker count).
%      - Hybrid Strategy: Automatically switches between coarse-grained (Strategy A) 
%        and fine-grained (Strategy B / Tail Mode) based on batch size AND n-complexity.
%      - For small n (< 40), Strategy A is preferred to avoid dispatch latency.
%      - For large n (>= 40), Strategy B ensures that compute-heavy BCa and Ranking 
%        tasks always saturate the full parallel pool.
%   3. Parallel Testing (Method-Major Scheduling): 
%      - Submits all test combinations to a shared parallel pool using Method-Major ordering.
%      - By scheduling BCa tests across all modes first, we minimize worker idle time 
%        at the end of a batch (Longest Processing Time First).
%      - Uses fetchNext to process results immediately as workers finish, reducing idle time.
%   4. Memory Management:
%      - Immediately clears simulation data and futures after each batch to prevent RAM growth.
%   5. Reporting:
%      - Updates progress incrementally and saves intermediate plots after each scenario.
%
% Inputs:
%   scenarios        - Struct array defining the data scenarios (n, Distribution).
%   params           - Struct containing the parameter sets for Thr, BCa, Rnk.
%   n_sims_per_cond  - Number of simulations per scenario.
%   refs             - Struct with reference B values.
%   cfg_base         - Base HERA configuration object. (e.g. holds system.target_memory, 
%                      simulation_seed, scenario_seed_offset, reference_seed_offset).
%   temp_dir         - Directory for temporary file artifacts.
%   styles, lang     - HERA design/language structs.
%   hWait            - Handle to the waitbar for progress updates.
%   out_dir          - Base output directory (for potential reports).
%   ts_str           - Timestamp string for file naming.
%   final_out_dir    - Directory for final PDFs.
%   colors           - Color matrix for plots.
%   modes            - Cell array of mode names.
%   limits           - Struct with B limits.
%   N                - Number of candidates/datasets.
%
% Outputs:
%   results          - Struct array (one per scenario) containing:
%                      - name, n, Dist, DataSummary : Metadata
%                      - thr, bca, rnk             : Nested metrics (err, cost, fail)
%                      - eff_all, eff_median       : Real effect size statistics (Cliff's d)
%
% Author: Lukas von Erdmannsdorff

    %% Main Execution Loop
    if isfield(params, 'thr')
        num_modes = length(params.thr);
    elseif isfield(params, 'bca')
        num_modes = length(params.bca);
    elseif isfield(params, 'rnk')
        num_modes = length(params.rnk);
    else
        num_modes = 3;
    end
    
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
        cfg_base.num_workers = num_workers; % Store back to prevent redundant calls in workers
    end
    effective_memory_mb = TARGET_MEMORY / max(1, num_workers);
    
    base_seed = 123;
    if isfield(cfg_base, 'simulation_seed') && isnumeric(cfg_base.simulation_seed)
        base_seed = cfg_base.simulation_seed;
    end
    
    % Offset scenario seeds dynamically to prevent overlaps during extensive runs
    scenario_seed_offset = max(10000, n_sims_per_cond * 3);
    if isfield(cfg_base, 'scenario_seed_offset') && isnumeric(cfg_base.scenario_seed_offset)
        scenario_seed_offset = cfg_base.scenario_seed_offset;
    end
    
    reference_seed_offset = 5000;
    if isfield(cfg_base, 'reference_seed_offset') && isnumeric(cfg_base.reference_seed_offset)
        reference_seed_offset = cfg_base.reference_seed_offset;
    end
    
    reference_step_offset = 1000;
    if isfield(cfg_base, 'reference_step_offset') && isnumeric(cfg_base.reference_step_offset)
        reference_step_offset = cfg_base.reference_step_offset;
    end

    
    % Initialize Result Structure
    scenario_res = repmat(struct('name', '', 'n', 0, 'Dist', '', 'DataSummary', '', ...
                                 'thr', init_storage(n_sims_per_cond, num_modes), ...
                                 'bca', init_storage(n_sims_per_cond, num_modes), ...
                                 'rnk', init_storage(n_sims_per_cond, num_modes), ...
                                 'eff_all', zeros(n_sims_per_cond, 1), ...
                                 'eff_median', 0), ...
                          length(scenarios), 1);
                      
    % Global Progress Tracking
    num_selected_global = 0;
    if isfield(params, 'thr'), num_selected_global = num_selected_global + 1; end
    if isfield(params, 'bca'), num_selected_global = num_selected_global + 1; end
    if isfield(params, 'rnk'), num_selected_global = num_selected_global + 1; end
    total_tests_global = length(scenarios) * n_sims_per_cond * num_selected_global * num_modes;
    global_tests_completed = 0;

    % Loop over Scenarios
    for sc_idx = 1:length(scenarios)
        sc = scenarios(sc_idx);
        scenario_res(sc_idx).name = sc.name; 
        scenario_res(sc_idx).n = sc.n; 
        scenario_res(sc_idx).Dist = sc.Dist; 
        scenario_res(sc_idx).DataSummary = sc.DataSummary;
        
        fprintf('\n--- Scenario %d/%d: %s ---\n', sc_idx, length(scenarios), sc.name);
        t_scenario = tic;
        
        % Dynamic batch sizing
        n_pairs = nchoosek(N, 2);
        bytes_per_double = 8;
        bytes_per_int = 4;
        
        % Determine peak B-value from parameters for realistic RAM estimation
        max_B = 0;
        for m_cfg = 1:num_modes
            if isfield(params, 'thr'), max_B = max(max_B, params.thr{m_cfg}.end); end
            if isfield(params, 'bca'), max_B = max(max_B, params.bca{m_cfg}.end); end
            if isfield(params, 'rnk'), max_B = max(max_B, params.rnk{m_cfg}.end); end
        end
        if isfield(refs, 'bca') && isnumeric(refs.bca), max_B = max(max_B, refs.bca); end
        if isfield(refs, 'thr') && isnumeric(refs.thr), max_B = max(max_B, refs.thr); end
        
        % Memory per simulation task (one mode/method at a time)
        % Peak is dominated by the bootstrap result matrix [Pairs x B] and indices [n x B].
        bytes_per_task = (sc.n * N * bytes_per_double) + ...                  % Input Data
                         (n_pairs * max_B * bytes_per_double) + ...           % Bootstrap Matrix peak (one task)
                         (sc.n * max_B * bytes_per_int);                      % Shuffle Indices
        
        % --- DRAS (Dynamic Resource-Aware Scaling) ---
        % Determine how many simulations we can safely process in parallel.
        % We prioritze filling the scenario (e.g. 50 sims) to eliminate the tail effect,
        % as long as the memory per worker is respected for the *active* tasks.
        if (bytes_per_task / 1024^2) <= effective_memory_mb
            % Standard case: Each task's peak RAM (bootstrap matrix) fits in a worker's slice.
            % We allow the batch to cover the entire scenario to saturate the parfeval queue.
            sims_per_batch = n_sims_per_cond;
        else
            % Very low RAM / Huge B: Scale down batch size to fit.
            sims_per_batch = max(2, floor(TARGET_MEMORY / (bytes_per_task / 1024^2)));
        end
        sims_per_batch = min(sims_per_batch, n_sims_per_cond);
        
        fprintf('  [Config: sims_per_batch=%d for n=%d]\n', sims_per_batch, sc.n);
        
        %% Interleaved Batch Loop
        % Divide simulations into small batches to keep memory footprint low
        for batch_start = 1:sims_per_batch:n_sims_per_cond
            t_batch = tic;
            batch_end = min(batch_start + sims_per_batch - 1, n_sims_per_cond);
            batch_sims = batch_start:batch_end;
            num_in_batch = length(batch_sims);
            
            fprintf('  Batch %d-%d: Preparing Data... ', batch_start, batch_end);
            
            % 1. Generate Data & References (Hybrid Parallelization Strategy)
            sim_data_batch = cell(num_in_batch, 1);
            
            % Determine separate parallelization strategies for preparation and testing phases.
            
            % --- Preparation (Reference calculations) ---
            % Determine strategy based on worker saturation and task complexity.
            % Threshold A (Batch Size): If batch fills more than 50% of the pool.
            % Threshold B (Complexity): For small batches, switch to Tail Mode (Strategy B) 
            % only if n is large enough to justify parfor management overhead (Crossover n ~ 40).
            % For small n, coarse-grained Strategy A is preferred to avoid dispatch latency.
            if num_in_batch > (num_workers / 2)
                use_parallel_prep = true;
            else
                use_parallel_prep = sc.n < 40; 
            end
            
            % --- Testing (Thresholds, BCa, Ranking) ---
            % Since each simulation generates 9 tasks (3 methods x 3 modes), we can 
            % saturate the worker pool with individual tasks much longer than for Refs. 
            % Threshold: Only use Tail Mode (Strategy B) if the total number of sub-tasks 
            % in the batch is smaller than the available worker pool.
            num_selected = 0;
            if isfield(params, 'thr'), num_selected = num_selected + 1; end
            if isfield(params, 'bca'), num_selected = num_selected + 1; end
            if isfield(params, 'rnk'), num_selected = num_selected + 1; end
            total_batch_tests = num_in_batch * num_selected * num_modes;
            use_parallel_tests = (total_batch_tests >= num_workers);
            
                if use_parallel_prep
                    % Strategy A: Coarse-Grained (High Saturation)
                    % Parallelize the outer loop, one simulation per worker.
                    parfor i = 1:num_in_batch
                        sim_data_batch{i} = prepare_simulation(i, batch_sims, sc, sc_idx, ...
                            base_seed, scenario_seed_offset, reference_seed_offset, reference_step_offset, ...
                            N, refs, cfg_base, temp_dir, styles, lang);
                    end
                else
                    % Strategy B: Fine-Grained (Low Saturation / Tail Processing)
                    % Execute simulation loop serially to enable internal parallelization on the workers.
                    fprintf(' [Prepare Mode: Serial Outer Loop (%d sim(s) across %d workers)] ', num_in_batch, num_workers);
                    for i = 1:num_in_batch
                        sim_data_batch{i} = prepare_simulation(i, batch_sims, sc, sc_idx, ...
                            base_seed, scenario_seed_offset, reference_seed_offset, reference_step_offset, ...
                            N, refs, cfg_base, temp_dir, styles, lang);
                    end
                end

            
            % Eliminate redundant rmdir calls to maintain filesystem performance
            % (cleanup relies on overwriting and final batch complete instead)
            fprintf('Done (%.2fs).\n', toc(t_batch));
            
            % 2. Parallel Processing (Worker Balancing)
            futures = [];
            
            % Determine strategy for test execution:
            % For small batches, serial-outer execution allows the compute-heavy inner functions 
            % (BCa, Ranking) to utilize the full parallel pool instead of being siloed on a single worker.
            if use_parallel_tests
                % Strategy A: Coarse-Grained (High Saturation)
                % Submit all test combinations to the shared pool using Method-Major scheduling.
                futures = parallel.FevalFuture.empty(0, 0);
                
                % Submit All Tasks (Longest Processing Time First)
                % To minimize worker idle time at the end of a batch, we schedule the most computationally expensive tasks first. 
                % Final complexity order: Method-Major (BCa > Ranking > Thresholds) ensures full core utilization.
                
                % 1. BCa (Method ID = 2) - Most expensive
                if isfield(params, 'bca')
                    for m = num_modes:-1:1
                        for i = 1:num_in_batch
                            s_idx = batch_sims(i);
                            % Minimize IPC overhead: Strip unused variables prior to parfeval
                            task_sd = rmfield(sim_data_batch{i}, {'ref_thr_struct', 'ds_names', 'base_rank', 'ref_thr_d', 'ref_rnk_mean'});
                            pp = map_params(params.bca{m});
                            futures(end+1) = parfeval(@run_single_test, 6, ...
                                s_idx, m, 2, task_sd, pp, sc.n, cfg_base, temp_dir, styles, lang, N);
                        end
                    end
                end
                
                % 2. Ranking (Method ID = 3) - Medium expense
                if isfield(params, 'rnk')
                    for m = num_modes:-1:1
                        for i = 1:num_in_batch
                            s_idx = batch_sims(i);
                            % Minimize IPC overhead: Strip large matrices not needed for ranking
                            task_sd = rmfield(sim_data_batch{i}, {'d_vals_all', 'rel_vals_all', 'ref_thr_d', 'ref_bca_width'});
                            pp = map_params(params.rnk{m});
                            futures(end+1) = parfeval(@run_single_test, 6, ...
                                s_idx, m, 3, task_sd, pp, sc.n, cfg_base, temp_dir, styles, lang, N);
                        end
                    end
                end
                
                % 3. Thresholds (Method ID = 1) - Least expensive
                if isfield(params, 'thr')
                    for m = num_modes:-1:1
                        for i = 1:num_in_batch
                            s_idx = batch_sims(i);
                            % Minimize IPC overhead: Retain core data only for thresholds
                            task_sd = rmfield(sim_data_batch{i}, {'d_vals_all', 'rel_vals_all', 'p_idx', 'ds_names', 'base_rank', 'ref_thr_struct', 'ref_bca_width', 'ref_rnk_mean'});
                            pp = map_params(params.thr{m});
                            futures(end+1) = parfeval(@run_single_test, 6, ...
                                s_idx, m, 1, task_sd, pp, sc.n, cfg_base, temp_dir, styles, lang, N); 
                        end
                    end
                end
                
                % Collect Results (Asynchronous / Load Balanced)
                tests_in_batch = numel(futures);
                tests_done_batch = 0;
                
                % Use loop to fetch exactly as many results as submitted
                for f_k = 1:tests_in_batch
                    [~, ret_s_idx, ret_m_idx, ret_method_id, ret_val, ret_cost, ret_fail] = fetchNext(futures);
                    % Internal Assignment Helper
                    scenario_res = assign_result(scenario_res, sc_idx, ret_method_id, ret_s_idx, ret_m_idx, ret_val, ret_cost, ret_fail, sim_data_batch, batch_start);
                    
                    tests_done_batch = tests_done_batch + 1;
                    global_tests_completed = global_tests_completed + 1;
                    update_progress(hWait, scenarios, sc, sc_idx, batch_start, batch_end, tests_done_batch, tests_in_batch, global_tests_completed, total_tests_global, num_modes, n_sims_per_cond, num_selected);
                end
            else
                % Strategy B: Fine-Grained (Low Saturation / Tail Processing)
                % Execute simulation loops serially on the client to enable internal parfor 
                % in BCa/Ranking to saturate the pool. 
                fprintf(' [Tail Mode: Serial Outer Loop (%d sim(s) across %d workers)] ', num_in_batch, num_workers);
                tests_in_batch = total_batch_tests;
                tests_done_batch = 0;
                
                % Method-Major scheduling consistency (BCa -> Ranking -> Thresholds)
                % 1. BCa (Method ID = 2) - Most expensive
                if isfield(params, 'bca')
                    for m = num_modes:-1:1
                        for i = 1:num_in_batch
                            s_idx = batch_sims(i);
                            % Minimize IPC overhead: Strip unused variables prior to parfeval
                            task_sd = rmfield(sim_data_batch{i}, {'ref_thr_struct', 'ds_names', 'base_rank', 'ref_thr_d', 'ref_rnk_mean'});
                            pp = map_params(params.bca{m});
                            [~, ~, ~, ret_val, ret_cost, ret_fail] = run_single_test(s_idx, m, 2, task_sd, pp, sc.n, cfg_base, temp_dir, styles, lang, N);
                            scenario_res = assign_result(scenario_res, sc_idx, 2, s_idx, m, ret_val, ret_cost, ret_fail, sim_data_batch, batch_start);
                            tests_done_batch = tests_done_batch + 1;
                            global_tests_completed = global_tests_completed + 1;
                            update_progress(hWait, scenarios, sc, sc_idx, batch_start, batch_end, tests_done_batch, tests_in_batch, global_tests_completed, total_tests_global, num_modes, n_sims_per_cond, num_selected);
                        end
                    end
                end
                
                % 2. Ranking (Method ID = 3) - Medium expense
                if isfield(params, 'rnk')
                    for m = num_modes:-1:1
                        for i = 1:num_in_batch
                            s_idx = batch_sims(i);
                            % Minimize IPC overhead: Strip large matrices not needed for ranking
                            task_sd = rmfield(sim_data_batch{i}, {'d_vals_all', 'rel_vals_all', 'ref_thr_d', 'ref_bca_width'});
                            pp = map_params(params.rnk{m});
                            [~, ~, ~, ret_val, ret_cost, ret_fail] = run_single_test(s_idx, m, 3, task_sd, pp, sc.n, cfg_base, temp_dir, styles, lang, N);
                            scenario_res = assign_result(scenario_res, sc_idx, 3, s_idx, m, ret_val, ret_cost, ret_fail, sim_data_batch, batch_start);
                            tests_done_batch = tests_done_batch + 1;
                            global_tests_completed = global_tests_completed + 1;
                            update_progress(hWait, scenarios, sc, sc_idx, batch_start, batch_end, tests_done_batch, tests_in_batch, global_tests_completed, total_tests_global, num_modes, n_sims_per_cond, num_selected);
                        end
                    end
                end
                
                % 3. Thresholds (Method ID = 1) - Least expensive
                if isfield(params, 'thr')
                    for m = num_modes:-1:1
                        for i = 1:num_in_batch
                            s_idx = batch_sims(i);
                            % Minimize IPC overhead: Retain core data only for thresholds
                            task_sd = rmfield(sim_data_batch{i}, {'d_vals_all', 'rel_vals_all', 'p_idx', 'ds_names', 'base_rank', 'ref_thr_struct', 'ref_bca_width', 'ref_rnk_mean'});
                            pp = map_params(params.thr{m});
                            [~, ~, ~, ret_val, ret_cost, ret_fail] = run_single_test(s_idx, m, 1, task_sd, pp, sc.n, cfg_base, temp_dir, styles, lang, N);
                            scenario_res = assign_result(scenario_res, sc_idx, 1, s_idx, m, ret_val, ret_cost, ret_fail, sim_data_batch, batch_start);
                            tests_done_batch = tests_done_batch + 1;
                            global_tests_completed = global_tests_completed + 1;
                            update_progress(hWait, scenarios, sc, sc_idx, batch_start, batch_end, tests_done_batch, tests_in_batch, global_tests_completed, total_tests_global, num_modes, n_sims_per_cond, num_selected);
                        end
                    end
                end
            end
            
            global_pct = global_tests_completed / total_tests_global;
            
            % 3. Cleanup & Report
            % Close graphics (excluding waitbar)
            all_figs = findall(0, 'Type', 'figure');
            if isscalar(hWait) && isgraphics(hWait) && isvalid(hWait)
                all_figs(all_figs == hWait) = [];
            end
            close(all_figs);
            clear futures sim_data_batch; % Then release memory
            
            dt = toc(t_batch);
            fprintf('    Batch finished in %.1fs (%.2fs/sim) [Progress: %.0f%%]\n', dt, dt/num_in_batch, global_pct*100);
        end
        
        %% End of Scenario
        scenario_res(sc_idx).eff_median = median(scenario_res(sc_idx).eff_all);
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
function [s_idx, m_idx, method_id, res_val, res_cost, res_fail] = run_single_test(s_idx, m_idx, method_id, sd, param, sample_size, cfg, tmp, styles, lang, num_candidates)
    % Runs a single convergence test and returns ID + Results
    % method_id: 1=Thr, 2=BCa, 3=Rnk
    
    % Generate independent streams for the simulations to avoid correlating
    % their internal paths directly with the reference computations.
    worker_sim_seed = sd.sim_seed + method_id * 100000;
    stream = RandStream('mlfg6331_64', 'Seed', worker_sim_seed);
    
    % Create isolated temporary directories for file I/O safety across parallel workers
    t_obj = getCurrentTask();
    worker_id = 0; if ~isempty(t_obj), worker_id = t_obj.ID; end
    worker_tmp = fullfile(tmp, sprintf('worker_%d', worker_id));
    if ~exist(worker_tmp, 'dir'), mkdir(worker_tmp); end

    res_val = 0; res_cost = 0; res_fail = false;
    
    try
        if method_id == 1
            % Thresholds
            c = cfg; c.bootstrap_thresholds = param;
            [d_t, ~, ~, ~, ~, ~, ~, conv_B, stab_data, h1, h2, h3, h4] = ...
                HERA.calculate_thresholds(sd.all_data, sample_size, c, worker_tmp, [], stream, styles, lang);
            
            % Memory Management: Close generated figures to prevent leaks in workers
            to_close = [h1(:); h2(:); h3(:); h4(:)];
            close(to_close(isgraphics(to_close)));

            % Threshold Error: Use Absolute Difference (instead of percentage) 
            % Since Cliff's Delta is normalized already, absolute deviation is more representative.
            res_val = (d_t(1) - sd.ref_thr_d); 
            res_cost = conv_B;
            res_fail = ~stab_data.converged;

        elseif method_id == 2
            % BCa
            c = cfg; c.bootstrap_ci = param;
            [conv_B, ci_d, ~, ~, ~, ~, ~, stab_data, h1, h2, h3, h4, h5] = ...
                HERA.calculate_bca_ci(sd.all_data, sd.d_vals_all, sd.rel_vals_all, sd.p_idx, sample_size, ...
                c, c.metric_names, worker_tmp, worker_tmp, [], stream, styles, lang, 'Sim');
            
            % Memory Management
            to_close = [h1(:); h2(:); h3(:); h4(:); h5(:)];
            close(to_close(isgraphics(to_close)));

            width = ci_d(1, 2) - ci_d(1, 1);
            res_val = (width - sd.ref_bca_width) / max(abs(sd.ref_bca_width), 1e-6) * 100;
            res_cost = conv_B;
            res_fail = ~stab_data.converged;
            
        elseif method_id == 3
            % Ranking
            c = cfg; c.bootstrap_ranks = param;
            [boot_r, conv_B, stab_data, h1, h2] = ...
                HERA.bootstrap_ranking(sd.all_data, sd.ref_thr_struct, c, sd.ds_names, sd.base_rank, sd.p_idx, sample_size, ...
                worker_tmp, worker_tmp, [], stream, styles, lang, 'Sim');
            
            % Memory Management
            to_close = [h1(:); h2(:)];
            close(to_close(isgraphics(to_close)));

            mean_r = mean(boot_r(2, :));
            res_val = (mean_r - sd.ref_rnk_mean) / max(abs(sd.ref_rnk_mean), 1e-6) * 100;
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
    s.err = NaN(n, num_modes); s.cost = NaN(n, num_modes); s.fail = false(n, num_modes); 
end

function d_all = generate_data_vectorized(sc, stream, N)
    % Generates synthetic data based on the scenario configuration
    d_all = zeros(sc.n, N);
    switch sc.Dist
        case 'Normal'
            % Medium Effect (Calibrated for Cliff's Delta ~ 0.40)
            % Median Pairwise Difference approx 1.50 / SD 2.0 (Step 0.75)
            means = sc.Base + (0:N-1) * sc.Step;
            r = randn(stream, sc.n, N);
            d_all = means + sc.SD * r;
        case 'Small Effect'
            % Small Effect (Calibrated for Cliff's Delta ~ 0.15)
            % Median Pairwise Difference approx 0.54 / SD 2.0 (Step 0.27)
            means = sc.Base + (0:N-1) * sc.Step;
            r = randn(stream, sc.n, N);
            d_all = means + sc.SD * r;
        case 'Large Effect'
            % Large Effect (Calibrated for Cliff's Delta ~ 0.70)
            % Median Pairwise Difference approx 3.00 / SD 2.0 (Step 1.50)
            means = sc.Base + (0:N-1) * sc.Step;
            r = randn(stream, sc.n, N);
            d_all = means + sc.SD * r;
        case 'Skewed'
            means = sc.Base + (0:N-1) * sc.Step;
            r = randn(stream, sc.n, N);
            d_all = exp(means + sc.SD * r);
        case 'Likert'
            centers = linspace(sc.Base, sc.End, N);
            r = randn(stream, sc.n, N);
            raw = centers + sc.SD * r;
            d_all = min(7, max(1, round(raw)));
        case 'Bimodal'
            for k = 1:N
                is_mode2 = rand(stream, sc.n, 1) > 0.5;
                vals = zeros(sc.n, 1);
                n1 = sum(~is_mode2); n2 = sum(is_mode2);
                vals(~is_mode2) = sc.Base + sc.SD * randn(stream, n1, 1);
                vals(is_mode2)  = sc.Base2 + sc.SD * randn(stream, n2, 1);
                % Adjusted to 1.5 (approx d=0.55) to test robustness, not power
                d_all(:, k) = vals + (k-1) * sc.Step;
            end
    end
end

%% Simulation Case Preparation
function sim_entry = prepare_simulation(i, batch_sims, sc, sc_idx, base_seed, scenario_seed_offset, reference_seed_offset, reference_step_offset, N, refs, cfg_base, temp_dir, styles, lang)
    % PREPARE_SIMULATION - Generates data and reference values for one simulation case.
    % This function is designed to be called either from a parfor (Strategy A) 
    % or a serial for-loop (Strategy B) in simulate.m.
    
    s_idx = batch_sims(i);
    
    % Create a unique temp dir for this worker to prevent collisions.
    t_obj = getCurrentTask();
    if isempty(t_obj)
        worker_id = 0; % Client/serial execution
    else
        worker_id = t_obj.ID;
    end
    
    worker_temp_dir = fullfile(temp_dir, sprintf('worker_%d', worker_id));
    if ~exist(worker_temp_dir, 'dir'), mkdir(worker_temp_dir); end
    
    % Bit-Perfect Seeding Constraint
    sim_seed = base_seed + (sc_idx-1)*scenario_seed_offset + s_idx;
    ref_seed = sim_seed + reference_seed_offset;
    
    % Generate Data
    dataStream = RandStream('mlfg6331_64', 'Seed', sim_seed);
    d_all = generate_data_vectorized(sc, dataStream, N);
    all_data = {d_all};
    p_idx_sim = nchoosek(1:N, 2);
    ds_names = cellstr("D" + (1:N));
    
    % Calc Real Effects (Fast)
    eff = HERA.test.TestHelper.calculate_real_effects(all_data, 1);
    
    % Reference Calculations
    
    % Ref: Thresholds
    ref_seed_thr = ref_seed + 1 * reference_step_offset;
    refStream = RandStream('mlfg6331_64', 'Seed', ref_seed_thr);
    if isstruct(refs.thr)
        c_ref_thr = cfg_base; c_ref_thr.bootstrap_thresholds = map_params(refs.thr);
        man_thr = [];
    else
        c_ref_thr = cfg_base; man_thr = refs.thr;
    end
    [ref_d_t, ref_r_t, ~, ~, d_vals_all, rel_vals_all, ~, ~, ~, h1, h2, h3, h4] = ...
        HERA.calculate_thresholds(all_data, sc.n, c_ref_thr, worker_temp_dir, man_thr, refStream, styles, lang);
    ref_thr_struct = struct('d_thresh', ref_d_t, 'rel_thresh', ref_r_t);
    
    % Memory Management
    close([h1(:); h2(:); h3(:); h4(:)]);
    
    % Ref: BCa (Internal parfor activates if called from client thread)
    ref_seed_bca = ref_seed + 2 * reference_step_offset;
    if isfield(cfg_base.system, 'selected_methods') && ismember('bca', cfg_base.system.selected_methods)
        refStream = RandStream('mlfg6331_64', 'Seed', ref_seed_bca);
        if isstruct(refs.bca)
            c_ref_bca = cfg_base; c_ref_bca.bootstrap_ci = map_params(refs.bca);
            man_bca = [];
        else
            c_ref_bca = cfg_base; man_bca = refs.bca;
        end
        [~, ref_ci_d, ~, ~, ~, ~, ~, ~, h1, h2, h3, h4, h5] = HERA.calculate_bca_ci(all_data, d_vals_all, rel_vals_all, p_idx_sim, sc.n, ...
            c_ref_bca, cfg_base.metric_names, worker_temp_dir, worker_temp_dir, man_bca, refStream, styles, lang, 'Ref');

        % Memory Management
        close([h1(:); h2(:); h3(:); h4(:); h5(:)]);
    else
        ref_ci_d = [0 0];
    end

    % Ref: Ranking
    [~, base_rank] = HERA.calculate_ranking(all_data, eff, ref_thr_struct, cfg_base, ds_names, p_idx_sim, lang);
    ref_seed_rnk = ref_seed + 3 * reference_step_offset;
    if isfield(cfg_base.system, 'selected_methods') && ismember('rnk', cfg_base.system.selected_methods)
        refStream = RandStream('mlfg6331_64', 'Seed', ref_seed_rnk);
        if isstruct(refs.rnk)
            c_ref_rnk = cfg_base; c_ref_rnk.bootstrap_ranks = map_params(refs.rnk);
            man_rnk = [];
        else
            c_ref_rnk = cfg_base; man_rnk = refs.rnk;
        end
        [boot_r_ref, ~, ~, h1, h2] = HERA.bootstrap_ranking(all_data, ref_thr_struct, c_ref_rnk, ds_names, base_rank, p_idx_sim, sc.n, ...
            worker_temp_dir, worker_temp_dir, man_rnk, refStream, styles, lang, 'Ref');

        % Memory Management
        close([h1(:); h2(:)]);
    else
        boot_r_ref = [0 0; 0 0];
    end
    
    % Store Package
    sim_entry = struct(...
        'all_data', {all_data}, 'd_vals_all', d_vals_all, 'rel_vals_all', rel_vals_all, ...
        'p_idx', p_idx_sim, 'ds_names', {ds_names}, 'base_rank', base_rank, ...
        'ref_thr_struct', ref_thr_struct, 'ref_thr_d', ref_d_t(1), ...
        'ref_bca_width', ref_ci_d(1,2) - ref_ci_d(1,1), ...
        'ref_rnk_mean', mean(boot_r_ref(2,:)), 'ref_seed', ref_seed, 'sim_seed', sim_seed, ...
        'eff_median', median(abs(d_vals_all(:))), 'n', sc.n);
end

%% Internal Helpers for Hybrid Execution
function scenario_res = assign_result(scenario_res, sc_idx, ret_method_id, ret_s_idx, ret_m_idx, ret_val, ret_cost, ret_fail, sim_data_batch, batch_start)
    % Assigns worker results to the result structure (Centralized to ensure consistency)
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
    
    % Store scenario-wide constant metadata (once per simulation)
    if ret_method_id == 1 % arbitrary choice
        i_in_batch = ret_s_idx - batch_start + 1;
        scenario_res(sc_idx).eff_all(ret_s_idx) = sim_data_batch{i_in_batch}.eff_median;
    end
end

function update_progress(hWait, scenarios, sc, sc_idx, batch_start, batch_end, tests_done_batch, tests_in_batch, global_tests_completed, total_tests_global, num_modes, n_sims_per_cond, num_selected)
    % Updates the global waitbar and console info
    if isempty(hWait) || ~isgraphics(hWait) || ~isvalid(hWait), return; end
    
    if nargin < 13
        num_selected = 3;
    end
    
    global_pct = global_tests_completed / total_tests_global;
    if global_pct >= 1 && tests_done_batch < tests_in_batch
        global_pct = 0.99;
    end
    
    tests_per_sim = num_selected * num_modes;
    total_tests_in_scenario = n_sims_per_cond * tests_per_sim;
    tests_done_in_scenario = (batch_start - 1) * tests_per_sim + tests_done_batch;
    sims_done_in_scenario = batch_start - 1 + ceil(tests_done_batch / tests_per_sim);
    
    waitbar(global_pct, hWait, ...
        sprintf('Scenario %d/%d (%s)\nBatch %d-%d | Sims %d/%d | Test %d/%d | %.0f%%', ...
        sc_idx, length(scenarios), sc.name, batch_start, batch_end, ...
        sims_done_in_scenario, n_sims_per_cond, ...
        tests_done_in_scenario, total_tests_in_scenario, ...
        global_pct * 100));
end

