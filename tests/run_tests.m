function run_tests()
% RUN_TESTS - Wrapper script to execute HERA unit tests from outside the package.
%
% Syntax:
%   run_tests()
%   matlab -batch "run('tests/run_tests.m')"
%
% Description:
%   This function serves as an entry point for the Continuous Integration (CI) pipeline
%   and manual testing. It ensures that the project root is correctly added to the
%   MATLAB path so that the '+HERA' package is visible. It then triggers the main
%   unit test suite located in 'HERA.run_unit_test'.
%
% Workflow:
%   1.  Path Configuration: Determine the script's location and add the project root to the path.
%   2.  Environment Setup: Create a dedicated log directory for test artifacts.
%   3.  Test Execution: Call the main 'HERA.run_unit_test' function.
%   4.  Error Handling: Catch errors, display stack traces, and manage exit codes for CI environments.
%
% Inputs:
%   None
%
% Outputs:
%   None (Console output and exit codes).
%
% Author:   Lukas von Erdmannsdorff

    %% 1. Path Configuration
    % Get the full path of this script to locate the project root relative to it.
    script_path = fileparts(mfilename('fullpath'));
    
    % The project root is one level up from the 'tests/' folder.
    project_root = fileparts(script_path);
    
    % Add project root to path so the +HERA package is visible to MATLAB.
    addpath(project_root);
    
    fprintf('=======================\n');
    fprintf('Starting HERA CI Runner\n');
    fprintf('=======================\n');
    fprintf('Project Root: %s\n', project_root);
    
    try
        %% 2. Environment Setup
        % Define a log folder inside tests/logs to store test outputs.
        log_dir = fullfile(script_path, 'logs');
        if ~exist(log_dir, 'dir')
            mkdir(log_dir);
        end
        
        %% 3. Test Execution
        % Run the package-based unit test suite.
        % We pass the log directory to ensure logs are stored in the correct location.
        HERA.run_unit_test(log_dir);
        
        fprintf('\nSUCCESS: All tests executed without critical errors.\n');
        
    catch ME
        %% 4. Error Handling
        % Report failure and stack trace for debugging in CI logs.
        fprintf('\nFAILURE: Test execution failed.\n');
        fprintf('Error Message: %s\n', ME.message);
        fprintf('Stack Trace:\n');
        for k = 1:length(ME.stack)
            fprintf('  File: %s, Line: %d, Name: %s\n', ...
                ME.stack(k).file, ME.stack(k).line, ME.stack(k).name);
        end
        
        % Return a non-zero exit code to signal failure to the CI system (e.g., GitHub Actions).
        if ~usejava('desktop') % Only exit if running in batch mode (no GUI)
            exit(1); 
        else
            error('Test execution failed. See details above.');
        end
    end
end
