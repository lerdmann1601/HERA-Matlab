function [batch_size, num_batches] = get_batch_config(config, num_iterations, bytes_per_iteration)
% GET_BATCH_CONFIG - Calculates optimal batch sizes for parallel processing.
%
% Syntax:
%   [batch_size, num_batches] = HERA.run.get_batch_config(config, num_iterations, bytes_per_iteration)
%
% Description:
%   Calculates the optimal batch configuration for parallel processing to ensure
%   memory stability. It dynamically adjusts batch sizes based on the system's 
%   available resources and the estimated memory footprint per iteration.
%
% Workflow:
%   1. System Configuration:
%      Retrieves memory limits ('target_memory') and batch constraints 
%      ('min_batch_size') from the config struct.
%   2. Resource Detection:
%      Determines the available worker count (via config or feature detection) 
%      to calculate the effective memory available per worker.
%   3. Safety Limit Calculation:
%      Estimates the total memory requirement for the requested iterations.
%   4. Batch Sizing:
%      - If Total Memory <= Effective Memory: Returns a single batch (all items).
%      - If Total Memory > Effective Memory: Splits iterations into chunks 
%        that fit within the effective memory limit per worker.
%   5. Final Validation:
%      Applies hard caps and defensive checks (NaN/Inf) to ensure a valid 
%      integer output for downstream parallel loops.
%
% Inputs:
%   config              - (struct) Configuration struct containing:
%                           .system.target_memory (optional, default: 200 MB)
%                           .system.min_batch_size (optional, default: 100)
%                           .num_workers (optional)
%   num_iterations      - (int) Total number of iterations/samples to process (e.g., B).
%   bytes_per_iteration - (double) Estimated memory usage per single iteration in bytes.
%
% Outputs:
%   batch_size          - (double) Number of iterations to process in one batch.
%   num_batches         - (int32) Total number of batches required.
%
% Example:
%   [b_size, n_batches] = HERA.run.get_batch_config(cfg, 1000, 800);
%
% Author: Lukas von Erdmannsdorff

    arguments
        config (1,1) struct
        num_iterations (1,1) double
        bytes_per_iteration (1,1) double
    end

    % 1. Determine Target Memory (default: 200 MB)
    %    This is the maximum RAM we want to consume per worker "slot".
    target_memory = 200;
    if isfield(config, 'system') && isfield(config.system, 'target_memory')
        sys_mem = config.system.target_memory;
        if ~isempty(sys_mem) && isnumeric(sys_mem)
            target_memory = sys_mem(1);
        end
    end

    % 2. Determine Minimum Batch Size (default: 100)
    %    Prevents overhead from becoming too dominant for tiny batches.
    min_batch_size = 100;
    if isfield(config, 'system') && isfield(config.system, 'min_batch_size')
        sys_min = config.system.min_batch_size;
        if ~isempty(sys_min) && isnumeric(sys_min)
            min_batch_size = sys_min(1);
        end
    end

    % 3. Determine Number of Workers
    %    Used to scale effective memory down (assuming data broadcast overhead).
    if isfield(config, 'num_workers') && isnumeric(config.num_workers) && ~isempty(config.num_workers)
        num_workers = config.num_workers(1);
    else
        % Robust feature detection for core count
        num_workers = feature('numcores');
    end
    
    % Safety: Ensure at least 1 worker to avoid division by zero
    num_workers = max(1, num_workers);

    % 4. Calculate Effective Memory per Batch
    %    Accounts for the fact that parallel workers share system RAM.
    effective_memory = target_memory / num_workers;

    % 5. Estimate Total Memory Needed
    %    How much RAM would satisfying the full request take?
    total_memory_needed = (double(num_iterations) * double(bytes_per_iteration)) / (1024^2); % in MB

    % 6. Determine Batch Size
    %    Logic: If it fits, take all. If not, chunk it.
    if total_memory_needed <= double(effective_memory)
        batch_size = double(num_iterations);
    else
        % Calculate max items that fit in effective_memory
        max_items_fitting = floor((double(effective_memory) * 1024^2) / double(bytes_per_iteration));
        
        % Constrain between min_batch_size and hard cap (e.g., 20,000 to prevent timeout/overhead issues)
        hard_cap = 20000; 
        
        batch_size = max(min_batch_size, min(max_items_fitting, hard_cap));
    end

    % 7. Calculate Number of Batches
    calc_num_batches = double(ceil(double(num_iterations) ./ batch_size));
    
    % Handle potential array outputs (scalar enforcement)
    if numel(calc_num_batches) > 1
        calc_num_batches = calc_num_batches(1);
    end

    % 8. Final Safety Checks (Defense against Inf/NaN)
    if isempty(calc_num_batches) || any(isnan(calc_num_batches)) || any(isinf(calc_num_batches))
        calc_num_batches = 1;
        batch_size = num_iterations; % Fallback: process all at once if logic fails
    end

    % Cast to standard types for output
    num_batches = int32(round(calc_num_batches));
    batch_size = double(batch_size);
    
end
