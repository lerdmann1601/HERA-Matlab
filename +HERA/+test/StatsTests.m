classdef StatsTests < matlab.unittest.TestCase
    % STATSTESTS - Validation of low-level statistical functions and performance heuristics
    %
    % Description:
    %   This class verifies that the core statistical functions (jackknife, cliffs_delta)
    %   correctly accept and respect the new system configuration parameters
    %   (vec_limit, mat_limit). It ensures that limiting the sample size forces
    %   the algorithm to switch between vectorized and loop-based implementations,
    %   while producing bit-identical results.
    %
    % Scope:
    %   - Jackknife: Verifies switch between Vectorized (O(N^2)) and Loop (O(N)) modes.
    %   - Cliff's Delta: Verifies switch between Matrix (O(N^2)) and Rank-based (O(N log N)) modes.
    %   - Parameter Propagation: Ensures vec_limit and mat_limit are respected.
    %   - BCA Parallel Logic: Verifies that switching between serial and parallel (parfor) pre-computation
    %     in calculate_bca_ci produces identical results.
    %   - Fallback Logic: Verifies that missing configuration fields trigger correct internal defaults.
    %
    % Author: Lukas von Erdmannsdorff
    
    properties
        tempDir
        lang
        styles
    end
    
    methods (TestClassSetup)
        function setup(testCase)
            import HERA.test.TestHelper
            
            % Setup temp directory
            root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            testCase.tempDir = fullfile(root, 'tests', 'TestOutput_Stats');
            if exist(testCase.tempDir, 'dir')
               rmdir(testCase.tempDir, 's');
            end
            mkdir(testCase.tempDir);
            
            % Load resources
            [testCase.lang, testCase.styles] = TestHelper.get_test_resources();
        end
    end
    
    methods (TestClassTeardown)
        function teardown(testCase)
            if exist(testCase.tempDir, 'dir')
                try
                    rmdir(testCase.tempDir, 's');
                catch
                    % Ignore file locks during test
                end
            end
        end
    end
    
    methods (Test)
        
        function test_jackknife_limits(testCase)
            % Verify that jackknife accepts vec_limit and produces consistent results
            
            % 1. Create Data (N=50)
            rng(123);
            x = randn(50, 1);
            y = randn(50, 1) + 0.5;
            
            % 2. Run with default limit (150) -> Should be Vectorized (N < 150)
            [val_def, a_def] = HERA.stats.jackknife(x, y, 'delta');
            
            % 3. Run with low limit (10) -> Should force Loop Mode (N > 10)
            [val_loop, a_loop] = HERA.stats.jackknife(x, y, 'delta', 10);
            
            % 4. Run with high limit (1000) -> Should be Vectorized (N < 1000)
            [val_vec, a_vec] = HERA.stats.jackknife(x, y, 'delta', 1000);
            
            % 5. Compare Results (Must be identical)
            testCase.verifyEqual(val_loop, val_def, 'AbsTol', 1e-10, 'Loop result differs from default');
            testCase.verifyEqual(val_vec, val_loop, 'AbsTol', 1e-10, 'Vectorized result differs from loop');
            
            % Double check acceleration factor 'a'
            testCase.verifyEqual(a_vec, a_loop, 'AbsTol', 1e-10, 'Acceleration factor differs');
        end
        
        function test_cliffs_delta_limits(testCase)
            % Verify that cliffs_delta accepts mat_limit and produces consistent results
            
            % 1. Create Data (nx=30, ny=30 -> prod=900)
            rng(456);
            x = randn(30, 1);
            y = randn(30, 1) + 0.2;
            
            % 2. Run with default limit (30000) -> Should be Matrix Mode (900 < 30000)
            d_def = HERA.stats.cliffs_delta(x, y);
            
            % 3. Run with low limit (100) -> Should force Rank Mode (900 > 100)
            d_rank = HERA.stats.cliffs_delta(x, y, 100);
            
            % 4. Run with high limit (50000) -> Should be Matrix Mode (900 < 50000)
            d_mat = HERA.stats.cliffs_delta(x, y, 50000);
            
            % 5. Compare Results (Must be identical)
            testCase.verifyEqual(d_rank, d_def, 'AbsTol', 1e-10, 'Rank result differs from default');
            testCase.verifyEqual(d_mat, d_rank, 'AbsTol', 1e-10, 'Matrix result differs from rank');
        end
        
        function test_parfor_switch(testCase)
             % Verify that calculate_bca_ci produces identical results whether
             % executing serially (N < limit) or parallel (N > limit).
             % This ensures both branches in the optimization logic are correct.
             
             import HERA.calculate_bca_ci
             
             % 1. Setup minimal input
             % Use N=20 subjects, 2 methods (1 pair)
             % This is small enough to be fast, but we can force it to be "large" by setting thr=0.
             n = 20;
             rng(999);
             d1 = randn(n, 1);
             d2 = randn(n, 1) + 1;
             all_data = { [d1, d2] }; % 1 Metric, 2 Methods
             
             % Required dummy args
             d_vals_all = [0.5];         % Dummy
             rel_vals_all = [0.1];       % Dummy
             pair_idx_all = [1, 2];      % Compare Col 1 vs Col 2
             metric_names = {'M1'};
             graphics_dir = testCase.tempDir;
             csv_dir = testCase.tempDir;
             manual_B = 100; % Small fixed B for speed
             s = RandStream('mlfg6331_64', 'Seed', 12345);
             base_name = 'Test';
             
             % 2. Config A: Force Parallel (Limit = 0 < N)
             cfg_par = HERA.default();
             cfg_par.timestamp = '20250101_000000';
             % Important: Hardcode numeric memory because this unit test calls calculate_bca_ci directly,
             % bypassing 'setup_environment' where the string "auto" is resolved to a number.
             cfg_par.system.target_memory = 2000; 
             cfg_par.system.jack_parfor_thr = 0; % Force parfor
             cfg_par.bootstrap_ci.n_trials = 2;   % Minimal
             cfg_par.bootstrap_ci.B_start = 100;
             cfg_par.bootstrap_ci.B_end = 100;
             
             % 4. Assert Equivalence
             % Results must be IDENTICAL because the Random Stream 's' is controlled 
             % and the logic should be mathematically equivalent.
             
             % Wrap in evalc to suppress output during testing
             % We execute the call in the caller workspace context
             T_par = evalc(['[~, ci_d_par, ~, z0_d_par, a_d_par] = HERA.calculate_bca_ci(' ...
                 'all_data, d_vals_all, rel_vals_all, pair_idx_all, n, ' ...
                 'cfg_par, metric_names, graphics_dir, csv_dir, manual_B, s, testCase.styles, testCase.lang, base_name);']);
                 
             % 3. Config B: Force Serial (Limit = 100 > N)
             cfg_ser = cfg_par;
             cfg_ser.system.jack_parfor_thr = 100; % Force serial
             
             % Reset RNG to ensure bootstrap sampling is identical
             s.reset(); 
             
             T_ser = evalc(['[~, ci_d_ser, ~, z0_d_ser, a_d_ser] = HERA.calculate_bca_ci(' ...
                 'all_data, d_vals_all, rel_vals_all, pair_idx_all, n, ' ...
                 'cfg_ser, metric_names, graphics_dir, csv_dir, manual_B, s, testCase.styles, testCase.lang, base_name);']);
             
             testCase.verifyEqual(ci_d_par, ci_d_ser, 'AbsTol', 1e-10, 'Parallel CI differs from Serial CI');
             testCase.verifyEqual(z0_d_par, z0_d_ser, 'AbsTol', 1e-10, 'Parallel z0 differs from Serial z0');
             testCase.verifyEqual(a_d_par, a_d_ser, 'AbsTol', 1e-10, 'Parallel acceleration differs from Serial acceleration');
        end
        
        function test_defaults_fallback(testCase)
            % Verify that missing configuration fields trigger correct fallbacks.
            % This ensures backward compatibility and robustness against partial configs.
            
            import HERA.calculate_bca_ci
            
            % Setup minimal data
            n = 10;
            rng(123);
            all_data = {[randn(n,1), randn(n,1)]};
             d_vals_all = [0.5];
             rel_vals_all = [0.1];
             pair_idx_all = [1, 2];
             metric_names = {'TestMetric'};
             graphics_dir = testCase.tempDir;
             csv_dir = testCase.tempDir;
             manual_B = 10; % Minimal B
             s = RandStream('mlfg6331_64', 'Seed', 123);
             base_name = 'FallbackTest';
             
             % Create a config STRIPPED of 'system' field entirely
             % simulating an old config or incomplete setup.
             % We MUST provide target_memory because calculate_bca_ci relies on it for Batch Size,
             % but we omit the new fields (jack_parfor_thr, limits).
             
             % FIX: Start with full default to get ci_level etc.
             cfg_fallback = HERA.default();
             cfg_fallback.timestamp = '20250101_000000';
             
             % Overwrite 'system' with a minimal struct containing ONLY target_memory
             % This effectively removes 'jack_parfor_thr', 'jack_vec_limit' etc.
             cfg_fallback.system = struct('target_memory', 1000); 
             
             % Execute
             % Should NOT error. Should use internal default (300) for parfor threshold.
             % Since N=10 < 300, it should run serially.
             try
                 [~, ~] = evalc(['HERA.calculate_bca_ci(' ...
                     'all_data, d_vals_all, rel_vals_all, pair_idx_all, n, ' ...
                     'cfg_fallback, metric_names, graphics_dir, csv_dir, manual_B, s, testCase.styles, testCase.lang, base_name);']);
             catch ME
                 testCase.verifyFail(['Fallback logic failed: ' ME.message]);
             end
        end
        
    end
end
