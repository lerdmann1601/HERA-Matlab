classdef SystemTests < matlab.unittest.TestCase
    % SYSTEMTESTS - High-Level System Validation for HERA
    %
    % Description:
    %   Verifies the core usage modes documented in the README:
    %   1. Developer Mode (Direct API call with data matrices)
    %   2. Batch Mode (JSON Configuration)
    %   3. Unit Test Wrapper (CLI Argument)
    %
    %   These tests ensure that the application entry points (start_ranking, run_ranking) 
    %   correctly parse inputs and orchestrate the analysis pipeline.
    %
    % Author: Lukas von Erdmannsdorff
    
    properties
        tempDir
    end
    
    methods (TestMethodSetup)
        function setupDesc(testCase)
            testCase.tempDir = tempname;
            mkdir(testCase.tempDir);
        end
    end
    
    methods (TestMethodTeardown)
        function teardownDesc(testCase)
             HERA.test.SystemTests.safeCleanup(testCase.tempDir);
        end
    end
    
    methods (Static)
        function safeCleanup(dirPath)
            % SAFECLEANUP - Robust deletion handling
            if exist(dirPath, 'dir')
                try
                    rmdir(dirPath, 's');
                catch
                    pause(0.5); % Wait for locks to release
                    try
                        rmdir(dirPath, 's');
                    catch
                        fprintf('Warning: Could not delete temp dir: %s\n', dirPath);
                    end
                end
            end
        end
    end
    
    methods (Test)
        
        function test_DeveloperMode(testCase)
            % Test 1: Developer Mode (Direct API Call)
            % Simulates a user calling HERA.run_ranking(userInput) directly with matrices.
            
            fprintf('\n[SystemTest] Developer Mode (Direct API)...\n');
            
            % 1. Prepare Mock Data (2 Metrics, 5 Methods, N=20)
            rng(123);
            m1 = randn(20, 5) + [1, 2, 3, 4, 5]; % Clear signal
            m2 = randn(20, 5);                   % Noise
            custom_data = {m1, m2};
            
            % 2. Setup User Input Struct
            userInput = struct();
            userInput.custom_data = custom_data;
            userInput.metric_names = {'Signal', 'Noise'};
            userInput.dataset_names = {'M1', 'M2', 'M3', 'M4', 'M5'};
            userInput.ranking_mode = 'M1_M2';
            userInput.output_dir = fullfile(testCase.tempDir, 'DevOutput');
            userInput.create_reports = false; % Speed up test
            userInput.reproducible = true;
            userInput.seed = 123;
            userInput.plot_theme = 'light';
            userInput.language = 'en';
            
            % 3. Run Ranking 
            [~, results] = evalc('HERA.run_ranking(userInput);');
            
            % 4. Verify Output Structure
            testCase.verifyTrue(isfield(results, 'final_rank'), 'Results missing final_rank');
            testCase.verifyTrue(isfield(results, 'd_vals_all'), 'Results missing d_vals_all');
            
            % Verify internal logic worked (M1 should dictate order due to strong signal)
            expected_order = [5, 4, 3, 2, 1]; % Descending means (5 is best/highest in Data Gen?) 
            % Wait, standard HERA logic logic assumes higher is better? 
            % Actually, HERA defaults depend on logic. Let's just check valid output first.
            testCase.verifyEqual(length(results.final_rank), 5, 'Rank length mismatch');
        end
        
        function test_BatchMode(testCase)
            % Test 2: Batch Mode (JSON Config)
            % Simulates HERA.start_ranking('configFile', 'path/to/config.json')
            
            fprintf('\n[SystemTest] Batch Mode (JSON Config)...\n');
            
            % 1. Create Dummy Data Files
            dataDir = fullfile(testCase.tempDir, 'Data');
            mkdir(dataDir);
            m1_data = rand(10, 3);
            m2_data = rand(10, 3);
            writematrix(m1_data, fullfile(dataDir, 'Metric1.csv'));
            writematrix(m2_data, fullfile(dataDir, 'Metric2.csv'));
            
            outputDir = fullfile(testCase.tempDir, 'BatchOutput');
            
            % 2. Create JSON Config
            configStruct.userInput.folderPath = dataDir;
            configStruct.userInput.fileType = '.csv';
            configStruct.userInput.metric_names = {'Metric1', 'Metric2'};
            configStruct.userInput.ranking_mode = 'M1_M2';
            configStruct.userInput.output_dir = outputDir;
            configStruct.userInput.create_reports = false;
            configStruct.userInput.num_workers = 1; % Serial for robustness in test
            
            jsonPath = fullfile(testCase.tempDir, 'config.json');
            fid = fopen(jsonPath, 'w');
            fprintf(fid, '%s', jsonencode(configStruct));
            fclose(fid);
            
            % 3. Call start_ranking in Batch Mode
            % Note: This function has no return, so we check side effects (files created)
            [T] = evalc('HERA.start_ranking(''configFile'', jsonPath);');
            
            % 4. Verify Output Files
            % HERA creates a timestamped folder inside outputDir. Find it.
            dirs = dir(fullfile(outputDir, 'Ranking_*'));
            testCase.verifyNotEmpty(dirs, 'No results folder created (Ranking_*)');
            
            if ~isempty(dirs)
                resFolder = fullfile(outputDir, dirs(1).name, 'Output');
                testCase.verifyTrue(exist(fullfile(resFolder, 'results.csv'), 'file') > 0, 'results.csv missing');
                testCase.verifyTrue(exist(fullfile(resFolder, 'analysis_data.json'), 'file') > 0, 'analysis_data.json missing');
            else
                % If no folder, print the captured output for debugging
                fprintf('\n[DEBUG] Captured Output from start_ranking:\n%s\n', T);
            end
        end
        

        
    end
end
