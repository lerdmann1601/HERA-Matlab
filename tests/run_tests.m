function run_tests()
% RUN_TESTS - Wrapper script to execute HERA unit tests from outside the package.
%
% This script ensures that the parent directory of +HERA is in the MATLAB path
% and then triggers the main unit test suite.
%
% Usage:
%   run_tests()
%   matlab -batch "run('tests/run_tests.m')"

    % Get the full path of this script
    scriptPath = fileparts(mfilename('fullpath'));
    
    % The project root is one level up from 'tests/'
    projectRoot = fileparts(scriptPath);
    
    % Add project root to path so +HERA is visible
    addpath(projectRoot);
    
    fprintf('Running HERA Unit Tests...\n');
    fprintf('Project Root: %s\n', projectRoot);
    
    try
        % Define log folder inside tests/logs
        logDir = fullfile(scriptPath, 'logs');
        if ~exist(logDir, 'dir')
            mkdir(logDir);
        end
        
        % Run the package-based unit test
        % Passing the log directory to the unit test function
        HERA.run_unit_test(logDir);
        
        fprintf('\nSUCCESS: All tests executed without critical errors.\n');
    catch ME
        fprintf('\nFAILURE: Test execution failed.\n');
        fprintf('Error: %s\n', ME.message);
        fprintf('Stack Trace:\n');
        for k = 1:length(ME.stack)
            fprintf('  File: %s, Line: %d, Name: %s\n', ...
                ME.stack(k).file, ME.stack(k).line, ME.stack(k).name);
        end
        exit(1); % Return non-zero exit code for CI
    end
end
