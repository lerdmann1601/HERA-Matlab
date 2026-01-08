classdef UnitTests < matlab.unittest.TestCase
    % UNITTESTS - Component Validation Suite for HERA
    %
    % Description:
    %   This class provides rigorous unit testing for the helper logic and
    %   configuration validators within the HERA framework. It specifically targets
    %   the +start/ConfigValidator.m class to ensure robust handling of user inputs.
    %
    % Scope:
    %   - Main Action Selection (Standard, Manual, Load)
    %   - Yes/No Boolean Parsing
    %   - Numeric Input Validation (Integer, Float, Ranges)
    %   - Metric Count Logic
    %   - Ranking Mode Configuration (M1_M2, M1_M3A, etc.)
    %   - Default Configuration Merging (Utils.fill_defaults)
    %   - Input Sanitization & Type Conversion (Utils.clean_struct)
    %   - Bootstrap Summary Formatting (UserInterface.get_bootstrap_string)
    %   - get_batch_config: Batch sizing logic
    %   - get_target_memory: Memory detection and safety limits
    %
    % Author: Lukas von Erdmannsdorff

    properties
        lang
    end
    
    methods (TestClassSetup)
        function setup(testCase)
            import HERA.test.TestHelper
            [testCase.lang, ~] = TestHelper.get_test_resources();
        end
    end
    
    methods (Test)
        
        %% --- Start Package: Validators ---  
        function validate_main_action(testCase)
            import HERA.start.ConfigValidator
            
            % Valid
            [isValid, err, val] = ConfigValidator.validate_main_action('s', testCase.lang);
            testCase.verifyTrue(isValid);
            testCase.verifyEqual(val, 's');
            
            [isValid, err, val] = ConfigValidator.validate_main_action('M', testCase.lang);
            testCase.verifyTrue(isValid);
            testCase.verifyEqual(val, 'm'); % Should normalize to lower
            
            % Invalid
            [isValid, err, ~] = ConfigValidator.validate_main_action('x', testCase.lang);
            testCase.verifyFalse(isValid);
            testCase.verifyNotEmpty(err);
            
            % Empty (Default)
            [isValid, err, val] = ConfigValidator.validate_main_action('', testCase.lang);
            testCase.verifyTrue(isValid);
            testCase.verifyEqual(val, 's'); 
        end
        
        function validate_yes_no(testCase)
            import HERA.start.ConfigValidator
            
            % Standard yes
            [isValid, ~, val] = ConfigValidator.validate_yes_no('y', true, testCase.lang);
            testCase.verifyTrue(isValid);
            testCase.verifyTrue(val);
            
            % Standard no
            [isValid, ~, val] = ConfigValidator.validate_yes_no('n', true, testCase.lang);
            testCase.verifyTrue(isValid);
            testCase.verifyFalse(val); 
            
            % Invalid
            [isValid, err, ~] = ConfigValidator.validate_yes_no('maybe', true, testCase.lang);
            testCase.verifyFalse(isValid);
            testCase.verifyNotEmpty(err);
        end
        
        function validate_numeric(testCase)
            import HERA.start.ConfigValidator
            
            % Valid Integer
            [isValid, ~, val] = ConfigValidator.validate_numeric('10', 5, true, testCase.lang);
            testCase.verifyTrue(isValid);
            testCase.verifyEqual(val, 10);
            
            % Invalid (Float when integer required)
            [isValid, err, ~] = ConfigValidator.validate_numeric('10.5', 5, true, testCase.lang);
            testCase.verifyFalse(isValid);
            testCase.verifyNotEmpty(err);
            
            % Valid Float
            [isValid, ~, val] = ConfigValidator.validate_numeric('0.05', 0.1, false, testCase.lang);
            testCase.verifyTrue(isValid);
            testCase.verifyEqual(val, 0.05);
            
            % Garbage
            [isValid, err, ~] = ConfigValidator.validate_numeric('abc', 5, true, testCase.lang);
            testCase.verifyFalse(isValid);
            testCase.verifyNotEmpty(err);
        end
        
        function validate_metric_count(testCase)
            import HERA.start.ConfigValidator
            
            % Valid: 1, 2, 3
            testCase.verifyTrue(ConfigValidator.validate_metric_count('1', testCase.lang));
            testCase.verifyTrue(ConfigValidator.validate_metric_count('2', testCase.lang));
            testCase.verifyTrue(ConfigValidator.validate_metric_count('3', testCase.lang));
            
            % Invalid: 4
            [isValid, err, ~] = ConfigValidator.validate_metric_count('4', testCase.lang);
            testCase.verifyFalse(isValid);
            testCase.verifyNotEmpty(err);
        end
        
        function validate_ranking_mode(testCase)
            import HERA.start.ConfigValidator
            
            % Mode 2 logic
            [isValid, ~, val] = ConfigValidator.validate_ranking_mode_2_metrics('1', testCase.lang);
            testCase.verifyTrue(isValid);
            testCase.verifyEqual(val, 'M1_M2');
            
            [isValid, ~, val] = ConfigValidator.validate_ranking_mode_2_metrics('2', testCase.lang);
            testCase.verifyTrue(isValid);
            testCase.verifyEqual(val, 'M1_M3A');
            
            % Invalid
            [isValid, err, ~] = ConfigValidator.validate_ranking_mode_2_metrics('3', testCase.lang);
            testCase.verifyFalse(isValid);
        end
        
        function validate_exit_mechanism(testCase)
            import HERA.start.UserInterface
            
            % Should throw error for exit commands
            testCase.verifyError(@() UserInterface.check_exit_command('q', testCase.lang), 'HERA:UserExit');
            testCase.verifyError(@() UserInterface.check_exit_command('exit', testCase.lang), 'HERA:UserExit');
            testCase.verifyError(@() UserInterface.check_exit_command('QUIT', testCase.lang), 'HERA:UserExit');
            
            % Should NOT throw error for normal input
            try
                UserInterface.check_exit_command('y', testCase.lang);
                UserInterface.check_exit_command('123', testCase.lang);
                UserInterface.check_exit_command('', testCase.lang);
            catch
                testCase.verifyFail('Check exit command threw error for valid input.');
            end
        end

        %% --- Start Package: Utilities & Logic ---
        function test_fill_defaults_basic(testCase)
            % Case: Basic flat structure merging
            % Description: Verifies that missing fields are filled with defaults.
            
            import HERA.start.Utils
            
            defaults.a = 1;
            defaults.b = 2;
            
            % Input missing 'b'
            inputStruct.a = 10;
            
            result = Utils.fill_defaults(inputStruct, defaults);
            
            testCase.verifyEqual(result.a, 10, 'Existing value should be preserved');
            testCase.verifyEqual(result.b, 2, 'Missing value should be filled from defaults');
        end
        
        function test_fill_defaults_recursive(testCase)
            % Case: Nested structure merging
            % Description: Verifies traversal into sub-structs.
            
            import HERA.start.Utils
            
            defaults.sub.x = 100;
            defaults.sub.y = 200;
            
            % Input has sub struct but is missing 'y'
            inputStruct.sub.x = 999;
            
            result = Utils.fill_defaults(inputStruct, defaults);
            
            testCase.verifyEqual(result.sub.x, 999);
            testCase.verifyEqual(result.sub.y, 200);
        end
        
        function test_fill_defaults_dynamic_fields(testCase)
            % Case: Dynamic fields (Paths, special configs)
            % Description: 'folderPath', 'output_dir' etc. should NOT be auto-filled if empty,
            % as they require explicit user interaction usually.
            
            import HERA.start.Utils
            
            defaults.folderPath = 'C:/Default/Path';
            defaults.output_dir = 'C:/Default/Out';
            defaults.normal_field = 'default_val';
            
            % Input has empty fields
            inputStruct.folderPath = '';
            inputStruct.output_dir = '';
            inputStruct.normal_field = '';
            
            result = Utils.fill_defaults(inputStruct, defaults);
            
            testCase.verifyEmpty(result.folderPath, 'Dynamic path field should remain empty');
            testCase.verifyEmpty(result.output_dir, 'Dynamic output field should remain empty');
            testCase.verifyEqual(result.normal_field, 'default_val', 'Standard field should be filled');
        end
        
        function test_clean_struct_numeric_conversion(testCase)
            % Case: Numeric strings to Double conversion
            % Description: "3" should become 3. "3.5" should become 3.5.
            
            import HERA.start.Utils
            
            inputStruct.val1 = '3';
            inputStruct.val2 = '4.5';
            inputStruct.val3 = 'not_a_number';
            
            result = Utils.clean_struct(inputStruct);
            
            testCase.verifyEqual(result.val1, 3);
            testCase.verifyClass(result.val1, 'double');
            
            testCase.verifyEqual(result.val2, 4.5);
            
            testCase.verifyEqual(result.val3, 'not_a_number', 'Non-numeric strings should be untouched');
        end
        
        function test_clean_struct_recursive_cell(testCase)
            % Case: Nested Cells and Structs
            % Description: Should deep clean complex data structures.
            
            import HERA.start.Utils
            
            % Nested: Struct -> Cell -> Struct
            inputStruct.list = {struct('id', '10'), struct('id', '20')};
            
            result = Utils.clean_struct(inputStruct);
            
            testCase.verifyEqual(result.list{1}.id, 10);
            testCase.verifyEqual(result.list{2}.id, 20);
            testCase.verifyClass(result.list{1}.id, 'double');
        end
        
        function test_bootstrap_summary_strings(testCase)
            % Case: Bootstrap Configuration Formatting
            % Description: Verifies that the summary string correctly reflects the config state.
            
            import HERA.start.UserInterface
            
            % 1. Manual Configuration
            configManual.manual_B_thr = 500;
            str = UserInterface.get_bootstrap_string(configManual, testCase.lang);
            % Expecting something containing the number 500
            testCase.verifySubstring(str, '500');
            
            % 2. Robust (Auto) Configuration
            configRobust.bootstrap_thresholds.smoothing_window = 5; % implies robust
            str = UserInterface.get_bootstrap_string(configRobust, testCase.lang);
            % Expecting 'Robust' (or localized equivalent, but we check key phrase)
            % The lang file usually has "Robust" or similar.
            % We will check that it is NOT empty.
            testCase.verifyNotEmpty(str);
            
            % 3. Simple (Auto) Configuration
            configSimple.bootstrap_thresholds.smoothing_window = []; % implies simple
            % It needs to be a struct though
            configSimple.bootstrap_thresholds.convergence_tolerance = 0.01;
            
            str = UserInterface.get_bootstrap_string(configSimple, testCase.lang);
            testCase.verifyNotEmpty(str);
        end


        %% --- Run Package: Execution Logic ---      
        function test_batch_fits_memory(testCase)
            % Case: Total memory needed is well within effective memory
            % Expectation: 1 batch, batch_size = num_iterations
            
            import HERA.run.get_batch_config
            
            config.system.target_memory = 1000; % 1000 MB
            config.num_workers = 1;
            
            iterations = 100; 
            bytes_per_iter = 1024; % 1 KB per iter -> Total ~0.1 MB
            
            [batch_size, num_batches] = get_batch_config(config, iterations, bytes_per_iter);
            
            testCase.verifyEqual(num_batches, int32(1));
            testCase.verifyEqual(batch_size, double(iterations));
        end
        
        function test_batch_splitting(testCase)
            % Case: Total memory needed exceeds effective memory
            % Expectation: >1 batches, batch_size < iterations
            
            import HERA.run.get_batch_config
            
            target_mem = 10; % 10 MB restricted memory
            config.system.target_memory = target_mem;
            config.num_workers = 1;
            
            % Request ~20 MB total (double the capacity)
            % 20,000 iters * 1024 bytes = ~20 MB
            iterations = 20000;
            bytes_per_iter = 1024; 
            
            [batch_size, num_batches] = get_batch_config(config, iterations, bytes_per_iter);
            
            testCase.verifyGreaterThan(num_batches, int32(1), 'Should require multiple batches');
            testCase.verifyLessThan(batch_size, double(iterations), 'Batch size should be reduced');
            
            % Verify calculated batch size roughly fits in 10MB
            % 10 MB * 1024^2 / 1024 bytes = 10,240 items max
            testCase.verifyLessThanOrEqual(batch_size, 10240, 'Batch size exceeds memory limit calculation');
        end
        
        function test_min_batch_size(testCase)
            % Case: Calculated fit is very small, should respect minimum
            
            import HERA.run.get_batch_config
            
            config.system.target_memory = 1; % 1 MB tiny memory
            config.system.min_batch_size = 50; 
            config.num_workers = 1;
            
            % Huge per-iter size to force small batches
            % 100 KB per iter -> 10 iters = 1 MB. 
            % 'max_items_fitting' would be ~10.
            bytes_per_iter = 100 * 1024; 
            iterations = 100;
            
            [batch_size, ~] = get_batch_config(config, iterations, bytes_per_iter);
            
            testCase.verifyEqual(batch_size, double(50), 'Should clamp to min_batch_size even if it exceeds memory target');
        end
        
        function test_worker_scaling(testCase)
            % Case: More workers reduce effective memory per worker
            
            import HERA.run.get_batch_config
            
            config.system.target_memory = 100; % 100 MB total system RAM target
            
            iterations = 10000;
            bytes_per_iter = 1024; % total 10 MB needed
            
            % 1 Worker: 100 MB effective > 10 MB needed -> 1 batch
            config.num_workers = 1;
            [~, n_b1] = get_batch_config(config, iterations, bytes_per_iter);
            
            % 20 Workers: 5 MB effective < 10 MB needed -> splitting needed
            config.num_workers = 20;
            [~, n_b2] = get_batch_config(config, iterations, bytes_per_iter);
            
            testCase.verifyEqual(n_b1, int32(1));
            testCase.verifyGreaterThan(n_b2, int32(1));
        end
        
        function test_defaults(testCase)
            % Case: Empty config, should use internal defaults (200MB, 100 min)
            
            import HERA.run.get_batch_config
            config = struct();
            
            % Small request should fit
            [bs, nb] = get_batch_config(config, 50, 100);
            
            testCase.verifyEqual(nb, int32(1));
            testCase.verifyEqual(bs, 50);
        end
        
        function test_get_target_memory_execution(testCase)
            % Just ensure it runs on the current platform without crashing
            import HERA.run.get_target_memory
            
            [target, ram, status] = get_target_memory();
            
            testCase.verifyNotEmpty(target);
            testCase.verifyClass(target, 'double');
            
            if isnan(ram)
                % If detection failed, check we have a status message
                testCase.verifyNotEmpty(status);
            else
                testCase.verifyGreaterThan(ram, 0);
            end
        end
        
        function test_target_memory_constraints(testCase)
            % Must always be >= 200 MB
            import HERA.run.get_target_memory
            
            target = get_target_memory();
            
            testCase.verifyGreaterThanOrEqual(target, 200);
        end
        
    end
end
