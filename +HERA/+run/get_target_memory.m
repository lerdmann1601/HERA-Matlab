function [target_mem, ram_gb, status_msg] = get_target_memory()
% GET_TARGET_MEMORY - Automatically calculates the target memory for HERA.
%
% Syntax:
%   [target_mem, ram_gb, status_msg] = HERA.run.get_target_memory()
%
% Description:
%   Determines the available physical RAM on the system and calculates an
%   optimal 'target_memory' value for HERA's batch processing.
%   The formula used is: target_memory = RAM_GB * 25.
%   (Example: 16 GB -> 400 MB target).
%
%   Minimum return value is 200 MB.
%
% Outputs:
%   target_mem - (double) Calculated target memory in MB.
%   ram_gb     - (double) Detected System RAM in GB (or NaN if failed).
%   status_msg - (string) Technical status message (empty if success).
%
% Author: Lukas von Erdmannsdorff

    % Default fallback
    min_mem = 200; 
    status_msg = "";
    ram_gb = NaN;
    
    try
        if ispc
            % Windows: Use MATLAB's built-in memory function
            [~, sysview] = memory;
            ram_bytes = sysview.PhysicalMemory.Total; 
            
        elseif ismac
            % macOS: Use sysctl
            [status, cmdout] = system('sysctl hw.memsize');
            if status == 0
                % cmdout is typically "hw.memsize: 17179869184"
                tokens = regexp(cmdout, '\d+', 'match');
                if ~isempty(tokens)
                    ram_bytes = str2double(tokens{1});
                else
                    error('Could not parse sysctl output.');
                end
            else
                error('sysctl command failed.');
            end
            
        elseif isunix
            % Linux: Read /proc/meminfo
            [status, cmdout] = system('grep MemTotal /proc/meminfo');
            if status == 0
                % Output example: "MemTotal:       16393932 kB"
                tokens = regexp(cmdout, '\d+', 'match');
                if ~isempty(tokens)
                    kb_val = str2double(tokens{1});
                    ram_bytes = kb_val * 1024;
                else
                    error('Could not parse /proc/meminfo.');
                end
            else
                error('Could not read /proc/meminfo.');
            end
        else
            error('Unknown operating system.');
        end
        
        % Convert to GB for the formula
        ram_gb = ram_bytes / (1024^3);
        
        % Formula: 25 MB per GB RAM
        % (16 GB -> 400 MB, 128 GB -> 3200 MB)
        calc_mem = ram_gb * 25;
        
        % Enforce minimum safety floor
        target_mem = max(calc_mem, min_mem);
        
        % Round to nearest integer for cleanliness
        target_mem = round(target_mem);
        
    catch ME
        % Fallback on error (return minimum safe value)
        target_mem = min_mem;
        status_msg = string(ME.message);
    end
end
