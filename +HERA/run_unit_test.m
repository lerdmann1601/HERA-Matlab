function run_unit_test(log_path_or_mode)
% RUN_UNIT_TEST - Comprehensive Scientific Validation Suite for the HERA Ranking.
%
% Syntax:
%   HERA.run_unit_test()                            1. Auto-Log Mode (Default)
%   HERA.run_unit_test('interactive')               2. Interactive Mode (Select Folder)
%   HERA.run_unit_test('path/to/custom/log/folder') 3. Custom Log Path Mode
%
% Description:
%   This function serves as the master runner for the HERA validation framework.
%   It employs a Clean Code Delegation Architecture to execute:
%     1. Unit Tests (Component Validation)
%     2. Scientific Tests (Algorithm Accuracy, Robustness, Logic)
%
%   Architecture:
%     - Controller: +HERA/+test/run_unit_test.m
%     - Subcontroller for scientific tests: +HERA/+test/ScientificSuite.m
%     - Scientific Test Cases: +HERA/+test/+cases/*.m
%     - Subcontroller for unit tests: +HERA/+test/UnitTests.m
%
% Workflow:
%   1. Environment Validation (Dependencies, Version)
%   2. Logging Setup (Diary persistence)
%   3. Parallel Pool Initialization
%   4. Test Suite Execution (via matlab.unittest)
%   5. Summary Reporting
%
% Inputs:
%   log_path_or_mode - (Optional) "interactive" or path to custom log folder.
%
% Author: Lukas von Erdmannsdorff

arguments
    log_path_or_mode (1,1) string = ""
end

    import HERA.test.TestHelper
    import matlab.unittest.TestSuite
    import matlab.unittest.TestRunner
    import matlab.unittest.plugins.TestReportPlugin
    
    clc;
    fprintf('===============================================\n');
    fprintf('Scientific Validation and Testing Unit for HERA\n');
    fprintf('===============================================\n');
    
    % Display Version Info if available
    try
        v = HERA.get_version();
        fprintf('HERA Version: %s\n', v);
    catch
        fprintf('HERA Version: Unknown (get_version failed)\n');
    end
    fprintf('\n');

    %% 1. Dependency Check
    % Ensure all required toolboxes are present before starting expensive tests.
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
    else
        fprintf('[System] All dependencies verified.\n');
    end

    %% 2. Setup Logging
    [log_folder, path_source] = TestHelper.get_writable_log_path();
    % Inputs:
    % log_path_or_mode - (Optional) "interactive" or path to custom log folder.
    % Handle arguments for custom path
    if log_path_or_mode ~= ""
        if strcmpi(log_path_or_mode, 'interactive')
            selected_dir = uigetdir(pwd, 'Select Log Folder');
            if ~isequal(selected_dir, 0), log_folder = selected_dir; end
        else
            if ~exist(log_path_or_mode, 'dir'), mkdir(log_path_or_mode); end
            log_folder = log_path_or_mode;
        end
    end

    timestamp_str = string(datetime('now'), 'yyyyMMdd_HHmmss');
    log_filename = "HERA_Validation_" + timestamp_str + ".txt";
    log_path = fullfile(log_folder, log_filename);
    
    if exist(log_path, 'file'), delete(log_path); end
    diary(char(log_path));
    
    fprintf('Log File: %s (%s)\n', log_path, path_source);
    
    %% 3. Start Parallel Pool
    pool = gcp('nocreate');
    if isempty(pool)
        try
            fprintf('[System] Starting Parallel Pool...\n');
            % SpmdEnabled=false reduces overhead since HERA only uses parfor.
            parpool('local', 'SpmdEnabled', false);
        catch
            fprintf('WARNING: Parallel Pool failed to start. Running serially.\n');
        end
    else
        fprintf('[System] Parallel Pool active: %d workers.\n', pool.NumWorkers);
    end
    
    %% 4. Execution
    fprintf('\nRunning Test Suites...\n');
    
    try
        % Combine suites
        % Combine suites
        suite = [TestSuite.fromClass(?HERA.test.UnitTests), ...
                 TestSuite.fromClass(?HERA.test.ScientificSuite), ...
                 TestSuite.fromClass(?HERA.test.SystemTests)];

        runner = TestRunner.withTextOutput;
        
        % Add JUnit XML Plugin for CI integration
        xmlFile = fullfile(log_folder, 'testResults.xml');
        runner.addPlugin(matlab.unittest.plugins.XMLPlugin.producingJUnitFormat(xmlFile));
        
        result = runner.run(suite);

        %% 5. Summary
        fprintf('\nTest Results Summary:\n');
        disp(table(result));
        
        % Check for failures and propagate error to caller (critical for CI)
        if any([result.Failed]) || any([result.Incomplete])
             error('HERA:TestsFailed', 'One or more tests failed or were incomplete.');
        end

    catch ME
        fprintf('\nCRITICAL ERROR during Test Execution:\n%s\n', ME.message);
        fprintf('Stack Trace:\n');
        for k = 1:length(ME.stack)
            fprintf('  File: %s, Line: %d, Name: %s\n', ME.stack(k).file, ME.stack(k).line, ME.stack(k).name);
        end
        % Re-throw to ensure the wrapper script sees the failure
        rethrow(ME);
    end
    
    diary off;
    delete(gcp('nocreate'));
end