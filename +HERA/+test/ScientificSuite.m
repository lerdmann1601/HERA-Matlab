classdef ScientificSuite < matlab.unittest.TestCase
    % SCIENTIFICSUITE - Master Test Suite for HERA Scientific Validation
    %
    % Description:
    %   This class orchestrates the execution of 19 complex scientific validation tests. 
    %   Each test corresponds to a specific logic requirement or math boundary case in the HERA.
    %
    % Architecture:
    %   - Delegates internal logic to: +HERA/+test/+cases/tXX_*.m
    %   - Manages shared environment (RNG, Config, Resources)
    %
    % List of Tests:
    %   T01: Small Sample (Wilcoxon Exactness)
    %   T02: Holm-Bonferroni Correction
    %   T03: Outlier Robustness
    %   T04: Logical Math Validations (SEM, BCa)
    %   T05: Zero-Variance Stability
    %   T06: Bootstrap Convergence
    %   T07: Initial Sort logic (Dominance/Mean)
    %   T08: Mode M1_M2 (Correction)
    %   T09: Mode M1_M3A (Swap)
    %   T10: Mode M1_M2_M3 (Sort)
    %   T11: Missing Data (NaN)
    %   T12: Cycle Detection (Artificial)
    %   T13: Pathological Data (Real Cycle)
    %   T14: Dynamic Cycle (Pipeline)
    %   T15: System Test (Static)
    %   T16: Full Integration (Dynamic)
    %   T17: Bootstrap Stability
    %   T18: Borda Count
    %   T19: Power Analysis
    %
    % Author: Lukas von Erdmannsdorff

    properties
        config
        thresholds
        n_subj
        styles
        lang
        s_global
    end
    
    methods (TestClassSetup)
        function setupEnvironment(testCase)
            import HERA.test.TestHelper
            
            % 1. Global RNG
            global_seed = 123;
            testCase.s_global = RandStream('mlfg6331_64', 'Seed', global_seed);
            RandStream.setGlobalStream(testCase.s_global);
            
            % 2. Standard Parameters
            testCase.n_subj = 50;
            timestamp_str = string(datetime('now'), 'yyyyMMdd_HHmmss');
            
            % 3. Default Configuration
            config = struct();
            config.alphas = [0.05, 0.05, 0.05]; 
            config.ci_level = 0.95; 
            config.timestamp = char(timestamp_str);
            config.metric_names = {'M1', 'M2', 'M3'}; 
            config.ranking_mode = 'M1_M2_M3'; 
            
            % Thresholds
            thresholds = struct();
            thresholds.d_thresh = [0.15, 0.15, 0.15];       
            thresholds.rel_thresh = [0.05, 0.05, 0.05];     
            thresholds.rel_thresh_b = [0.05, 0.05, 0.05];   
            thresholds.min_rel_thresh = [0, 0, 0];          
            
            % Bootstrap Config
            bs_config = struct();
            bs_config.B_start = 50; bs_config.B_step = 50; bs_config.B_end = 150;
            bs_config.n_trials = 3; bs_config.convergence_tolerance = 0.05;
            bs_config.min_steps_for_convergence_check = 2; bs_config.smoothing_window = 3;
            bs_config.convergence_streak_needed = 2;
            
            config.bootstrap_thresholds = bs_config;
            config.bootstrap_ci = bs_config;
            config.bootstrap_ranks = bs_config;
            
            testCase.config = config;
            testCase.thresholds = thresholds;
            
            % 4. Resources
            [testCase.lang, testCase.styles] = TestHelper.get_test_resources();
        end
    end

    methods (TestMethodTeardown)
        function cleanupFigures(testCase)
            % Ensure no figures from direct function calls persist between tests.
            close all force;
        end
    end
    
    methods (Test)
        
        function t01_SmallSample(testCase)
            import HERA.test.cases.t01_SmallSample
            passed = t01_SmallSample(testCase.config, testCase.thresholds, testCase.n_subj, testCase.styles, testCase.lang);
            testCase.verifyTrue(passed, "Test 01 (Small Sample) Failed");
        end
        
        function t02_HolmBonferroni(testCase)
            import HERA.test.cases.t02_HolmBonferroni
            passed = t02_HolmBonferroni(testCase.config, testCase.thresholds, testCase.n_subj, testCase.styles, testCase.lang);
            testCase.verifyTrue(passed, "Test 02 (Holm-Bonferroni) Failed");
        end
        
        function t03_OutlierRobustness(testCase)
            import HERA.test.cases.t03_OutlierRobustness
            passed = t03_OutlierRobustness(testCase.config, testCase.thresholds, testCase.n_subj, testCase.styles, testCase.lang);
            testCase.verifyTrue(passed, "Test 03 (Outlier Robustness) Failed");
        end
        
        function t04_MathValidation(testCase)
            import HERA.test.cases.t04_MathValidation
            passed = t04_MathValidation(testCase.config, testCase.thresholds, testCase.n_subj, testCase.styles, testCase.lang);
            testCase.verifyTrue(passed, "Test 04 (Math Validation) Failed");
        end
        
        function t05_Stability(testCase)
            import HERA.test.cases.t05_Stability
            passed = t05_Stability(testCase.config, testCase.thresholds, testCase.n_subj, testCase.styles, testCase.lang);
            testCase.verifyTrue(passed, "Test 05 (Stability) Failed");
        end
        
        function t06_BootstrapConvergence(testCase)
            import HERA.test.cases.t06_BootstrapConvergence
            passed = t06_BootstrapConvergence(testCase.config, testCase.thresholds, testCase.n_subj, testCase.styles, testCase.lang);
            testCase.verifyTrue(passed, "Test 06 (Convergence) Failed");
        end
        
        function t07_InitialSorting(testCase)
            import HERA.test.cases.t07_InitialSorting
            passed = t07_InitialSorting(testCase.config, testCase.thresholds, testCase.n_subj, testCase.styles, testCase.lang);
            testCase.verifyTrue(passed, "Test 07 (Initial Sorting) Failed");
        end
        
        function t08_ModeM1M2(testCase)
            import HERA.test.cases.t08_ModeM1M2
            passed = t08_ModeM1M2(testCase.config, testCase.thresholds, testCase.n_subj, testCase.styles, testCase.lang);
            testCase.verifyTrue(passed, "Test 08 (M1_M2) Failed");
        end
        
        function t09_ModeM1M3A(testCase)
            import HERA.test.cases.t09_ModeM1M3A
            passed = t09_ModeM1M3A(testCase.config, testCase.thresholds, testCase.n_subj, testCase.styles, testCase.lang);
            testCase.verifyTrue(passed, "Test 09 (M1_M3A) Failed");
        end
        
        function t10_ModeM1M2M3(testCase)
            import HERA.test.cases.t10_ModeM1M2M3
            passed = t10_ModeM1M2M3(testCase.config, testCase.thresholds, testCase.n_subj, testCase.styles, testCase.lang);
            testCase.verifyTrue(passed, "Test 10 (M1_M2_M3) Failed");
        end
        
        function t11_MissingData(testCase)
            import HERA.test.cases.t11_MissingData
            passed = t11_MissingData(testCase.config, testCase.thresholds, testCase.n_subj, testCase.styles, testCase.lang);
            testCase.verifyTrue(passed, "Test 11 (Missing Data) Failed");
        end
        
        function t12_CycleDetection(testCase)
            import HERA.test.cases.t12_CycleDetection
            passed = t12_CycleDetection(testCase.config, testCase.thresholds, testCase.n_subj, testCase.styles, testCase.lang);
            testCase.verifyTrue(passed, "Test 12 (Cycle Detection) Failed");
        end
        
        function t13_PathologicalData(testCase)
            import HERA.test.cases.t13_PathologicalData
            passed = t13_PathologicalData(testCase.config, testCase.thresholds, testCase.n_subj, testCase.styles, testCase.lang);
            testCase.verifyTrue(passed, "Test 13 (Pathological Data) Failed");
        end
        
        function t14_DynamicCycle(testCase)
            import HERA.test.cases.t14_DynamicCycle
            passed = t14_DynamicCycle(testCase.config, testCase.thresholds, testCase.n_subj, testCase.styles, testCase.lang);
            testCase.verifyTrue(passed, "Test 14 (Dynamic Cycle) Failed");
        end
        
        function t15_SystemTest(testCase)
            import HERA.test.cases.t15_SystemTest
            passed = t15_SystemTest(testCase.config, testCase.thresholds, testCase.n_subj, testCase.styles, testCase.lang);
            testCase.verifyTrue(passed, "Test 15 (System Test) Failed");
        end
        
        function t16_FullIntegration(testCase)
            import HERA.test.cases.t16_FullIntegration
            passed = t16_FullIntegration(testCase.config, testCase.thresholds, testCase.n_subj, testCase.styles, testCase.lang);
            testCase.verifyTrue(passed, "Test 16 (Full Integration) Failed");
        end
        
        function t17_BootstrapStability(testCase)
            import HERA.test.cases.t17_BootstrapStability
            passed = t17_BootstrapStability(testCase.config, testCase.thresholds, testCase.n_subj, testCase.styles, testCase.lang);
            testCase.verifyTrue(passed, "Test 17 (Bootstrap Stability) Failed");
        end
        
        function t18_BordaCount(testCase)
            import HERA.test.cases.t18_BordaCount
            passed = t18_BordaCount(testCase.config, testCase.thresholds, testCase.n_subj, testCase.styles, testCase.lang);
            testCase.verifyTrue(passed, "Test 18 (Borda Count) Failed");
        end
        
        function t19_PowerAnalysis(testCase)
            import HERA.test.cases.t19_PowerAnalysis
            passed = t19_PowerAnalysis(testCase.config, testCase.thresholds, testCase.n_subj, testCase.styles, testCase.lang);
            testCase.verifyTrue(passed, "Test 19 (Power Analysis) Failed");
        end
        
    end
end
