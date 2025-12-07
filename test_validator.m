% Test script for HERA.start.ConfigValidator
% Run this to verify validation logic

% Import necessary classes
import HERA.start.ConfigValidator;
import HERA.start.Utils;

fprintf('Running ConfigValidator Tests...\n');

% Load the real language file
try
    lang = Utils.language_code('en');
catch ME
    error('Failed to load language file via Utils: %s', ME.message);
end

% 1. Main Action
[ok, ~, val] = ConfigValidator.validate_main_action('s', lang); assert(ok && strcmp(val, 's'));
[ok, ~, val] = ConfigValidator.validate_main_action('', lang); assert(ok && strcmp(val, 's'));
[ok, msg, ~] = ConfigValidator.validate_main_action('x', lang); assert(~ok); 

% 2. Yes/No
[ok, ~, val] = ConfigValidator.validate_yes_no('y', false, lang); assert(ok && val == true);
[ok, ~, val] = ConfigValidator.validate_yes_no('n', true, lang); assert(ok && val == false);
[ok, ~, val] = ConfigValidator.validate_yes_no('', true, lang); assert(ok && val == true);
[ok, msg, ~] = ConfigValidator.validate_yes_no('x', true, lang); assert(~ok);

% 3. Numeric
[ok, ~, val] = ConfigValidator.validate_numeric('10', 0, true, lang); assert(ok && val==10);
[ok, ~, val] = ConfigValidator.validate_numeric('10.5', 0, false, lang); assert(ok && val==10.5);
[ok, msg, ~] = ConfigValidator.validate_numeric('10.5', 0, true, lang); assert(~ok); % Not integer
[ok, msg, ~] = ConfigValidator.validate_numeric('abc', 0, true, lang); assert(~ok);

% 4. Metric Count
[ok, ~, val] = ConfigValidator.validate_metric_count('1', lang); assert(ok && val==1);
[ok, ~, val] = ConfigValidator.validate_metric_count('3', lang); assert(ok && val==3);
[ok, msg, ~] = ConfigValidator.validate_metric_count('4', lang); assert(~ok);

% 5. Metric Order
[ok, ~, val] = ConfigValidator.validate_metric_order('1 2 3', 3, 3, lang); assert(ok && isequal(val, [1 2 3]));
[ok, msg, ~] = ConfigValidator.validate_metric_order('1 2 2', 3, 3, lang); assert(~ok); % Duplicate
[ok, msg, val] = ConfigValidator.validate_metric_order('1 5', 2, 5, lang); assert(ok && isequal(val, [1 5]));
[ok, msg, ~] = ConfigValidator.validate_metric_order('1 6', 2, 5, lang); assert(~ok); % Out of bound

% 6. Logic 2 Metrics
[ok, ~, val] = ConfigValidator.validate_ranking_mode_2_metrics('1', lang); assert(ok && strcmp(val, 'M1_M2'));
[ok, ~, val] = ConfigValidator.validate_ranking_mode_2_metrics('2', lang); assert(ok && strcmp(val, 'M1_M3A'));
[ok, msg, ~] = ConfigValidator.validate_ranking_mode_2_metrics('3', lang); assert(~ok);

fprintf('All tests passed!\n');
exit(0);
