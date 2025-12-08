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
        
    end
end
