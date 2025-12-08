classdef TestHelper
    % TESTHELPER - Static Utility Class for HERA Scientific Validation
    %
    % Description:
    %   Provides a centralized repository of shared helper functions used
    %   across all test cases. This ensures DRY (Don't Repeat Yourself)
    %   principles and standardizes data generation, formatting, and effect
    %   size calculation for testing.
    %
    % Categories:
    %   1. Data Generation (Exact, Efron, Waterfall)
    %   2. Wrapper Wrappers (Calculations)
    %   3. Output Formatting (Tables)
    %   4. System Utils (Paths, Resources)
    %
    % Author: Lukas von Erdmannsdorff
    
    methods (Static)
        
        function [log_folder, source] = get_writable_log_path()
            % Determine a valid location for log files.
            % Priority: 1. User Documents (Persistent), 2. TempDir (Fallback)
            
            app_folder_name = 'HERA_Ranking_Logs';
            
            % Try standard Documents folder first
            if ispc
                docs_path = fullfile(getenv('USERPROFILE'), 'Documents');
            else
                docs_path = fullfile(getenv('HOME'), 'Documents');
            end
            
            target_path = fullfile(docs_path, app_folder_name);
            
            % Attempt to create/access Documents folder
            try
                if ~exist(target_path, 'dir')
                    mkdir(target_path);
                end
                % Simple write test to confirm permissions
                testFile = fullfile(target_path, 'write_test.tmp');
                fid = fopen(testFile, 'w');
                if fid == -1
                    error('No write access');
                end
                fclose(fid);
                delete(testFile);
                
                log_folder = target_path;
                source = 'Documents (Persistent)';
                return;
            catch
                % Fallback to tempdir if Documents fails (e.g. Restricted User)
                log_folder = fullfile(tempdir, app_folder_name);
                if ~exist(log_folder, 'dir')
                    mkdir(log_folder);
                end
                source = 'TempDir (Fallback)';
            end
        end

        function success = check_result(actual, expected, label)
            % Validates the actual rank order against the expectation.
            % Helper for simple Pass/Fail printing.
            if isequal(actual(:)', expected(:)')
                fprintf('[Status] PASS: %s: Order is %s.\n', label, mat2str(actual(:)'));
                success = true;
            else
                fprintf('[Status] FAIL: %s\n', label);
                fprintf('   Expected: %s\n', mat2str(expected(:)'));
                fprintf('   Actual:   %s\n', mat2str(actual(:)'));
                success = false;
            end
        end
        
        function print_auto_table(headers, data, data_alignments, header_alignments)
            % PRINT_AUTO_TABLE - Visualizes a data table with dynamic column widths.
            
            % 1. Initialization & Defaults
            import HERA.test.TestHelper
            
            if nargin < 4
                header_alignments = data_alignments; % Fallback: Header = Data
            end

            num_cols = length(headers);
            num_rows = size(data, 1);
            
            % Initialize column widths based on header lengths
            col_widths = cellfun(@strlength, headers);
            
            % 2. Dynamic Width Calculation based on data content
            if num_rows > 0
                for r = 1:num_rows
                    for c = 1:num_cols
                        val = data{r,c};
                        if isnumeric(val) || islogical(val)
                            val_str = char(string(val)); 
                        elseif ischar(val)
                            val_str = val;
                        else
                            val_str = char(val);
                        end
                        col_widths(c) = max(col_widths(c), strlength(val_str));
                    end
                end
            end
            
            % 3. Apply Padding
            col_widths = col_widths + 2;
            
            % 4. Print Header (Using header_alignments)
            header_line = strjoin(arrayfun(@(c) TestHelper.format_text(headers{c}, col_widths(c), header_alignments{c}), ...
                1:num_cols, 'UniformOutput', false), '|');
            fprintf('   |%s|\n', header_line);
            
            % Separator
            total_width = sum(col_widths) + length(col_widths) + 1;
            fprintf('   %s\n', repmat('-', 1, total_width));
            
            % 5. Print Data Rows (Using data_alignments)
            for r = 1:num_rows
                row_str_parts = cell(1, num_cols);
                for c = 1:num_cols
                    val = data{r,c};
                    if isnumeric(val) || islogical(val)
                        str_val = char(string(val));
                    elseif ischar(val)
                        str_val = val;
                    else
                        str_val = char(val);
                    end
                    
                    % Use data_alignments here
                    row_str_parts{c} = TestHelper.format_text(str_val, col_widths(c), data_alignments{c});
                end
                fprintf('   |%s|\n', strjoin(row_str_parts, '|'));
            end
        end

        function output = format_text(text, width, alignment)
            % Aligns text within a fixed width (Left, Right, Center).
            text_len = strlength(text);
            padding = width - text_len;
            
            switch alignment
                case 'r'
                    output = [repmat(' ', 1, padding), text];
                case 'l'
                    output = [text, repmat(' ', 1, padding)];
                otherwise % 'c'
                    padding_left = floor(padding / 2);
                    padding_right = ceil(padding / 2);
                    output = [repmat(' ', 1, padding_left), text, repmat(' ', 1, padding_right)];
            end
        end
        
        %% Data Generation
        
        function data = generate_exact_data(n, means, std_dev)
            % Generates deterministic random normal data.
            % Uses input-based seeding to ensure reproducibility.
            
            seed_val = abs(sum(means) * 1000 + n + std_dev * 10);
            
            % Initialize a local stream with this specific seed.
            s_gen = RandStream('mlfg6331_64', 'Seed', uint32(seed_val));
            
            % Set this stream as global temporarily and store the previous stream.
            oldStream = RandStream.setGlobalStream(s_gen);
            
            % Start of Data Generation
            k = length(means);
            data = randn(n, k);
            if std(data(:)) == 0 % Handle zero variance case
                data = zeros(n,k) + means;
            else
                % Standardization: (X - mu) / sigma
                data = (data - mean(data)) ./ std(data); 
                % Scaling: X * new_sigma + new_mu
                data = data .* std_dev + means;
            end
            
            % Restore the original random stream
            RandStream.setGlobalStream(oldStream);
        end
        
        function [m1, m2] = generate_efron_data(n_subj)
            % Generates non-transitive Dice data (Efron's Dice).
            import HERA.test.TestHelper
            
            vals_A = [4; 4; 4; 4; 0; 0];
            vals_B = [3; 3; 3; 3; 3; 3];
            vals_C = [2; 2; 2; 2; 6; 6];
            vals_D = [1; 1; 1; 1; 5; 5];
            
            reps = ceil(n_subj / 6);
            A = repmat(vals_A, reps, 1); A = A(1:n_subj);
            B = repmat(vals_B, reps, 1); B = B(1:n_subj);
            C = repmat(vals_C, reps, 1); C = C(1:n_subj);
            D = repmat(vals_D, reps, 1); D = D(1:n_subj);
            
            m2 = [A, B, C, D];
            % M1 is provided as a transitive baseline
            m1 = TestHelper.generate_exact_data(n_subj, [40, 30, 20, 10], 1.0);
        end
        
        function [m1, m2, m3, names] = generate_waterfall_data(n_subj)
            % Generates a complex dataset triggering all HERA logic stages.
            % Designed for "Integration Test 17".
            
            s_gen = RandStream('mlfg6331_64', 'Seed', 9999);
            oldStream = RandStream.setGlobalStream(s_gen);

            names = arrayfun(@(x) sprintf('D%d', x), 1:15, 'UniformOutput', false);
            
            % Initialize matrices
            m1 = zeros(n_subj, 15);
            m2 = zeros(n_subj, 15);
            m3 = zeros(n_subj, 15);
            
            % Helper to add noise
            % scale=0.1 for signals (d > thresh)
            % scale=0.0 for neutrality (d = 0)
            add_vals = @(vals, scale) repmat(vals, n_subj, 1) + (randn(n_subj, length(vals)) * scale);
            
            % Base offset
            base = 100.0; 
            
            % Tier 1 (D1-D5): Mid Tier
            m1(:, 1:5) = add_vals(repmat(base + 20.0, 1, 5), 0.1);  
            m2(:, 1:5) = add_vals(repmat(base,        1, 5), 0.0);   
            m3(:, 1:5) = add_vals(repmat(base,        1, 5), 0.0);   
            
            % Tier 2 (D6-D10): Top Tier (M2 Correction) 
            m1(:, 6:10) = add_vals(repmat(base,        1, 5), 0.0);  
            m2(:, 6:10) = add_vals(repmat(base + 20.0, 1, 5), 0.1); 
            m3(:, 6:10) = add_vals(repmat(base,        1, 5), 0.0);  
            
            % Tier 3 (D11-D12): Logic 3A (Tie-Break) 
            m1(:, 11:12) = add_vals([base - 20.0, base - 20.0], 0.0); 
            m2(:, 11:12) = add_vals([base,        base       ], 0.0);     
            m3(:, 11:12) = add_vals([base - 5.0,  base + 5.0 ], 0.1); 
            
            % Tier 4 (D13-D15): Logic 3B (Sort)
            m1(:, 13:15) = add_vals([base - 40.0, base - 40.0, base - 40.0], 0.0); 
            m2(:, 13:15) = add_vals([base,        base,        base       ], 0.0);       
            m3(:, 13:15) = add_vals([base - 10.0, base,        base + 10.0], 0.1); 
            
            RandStream.setGlobalStream(oldStream); 
        end
        
        %% Calculations & Wrappers
        
        function effect_sizes = calculate_real_effects(all_data, num_metrics)
            % Calculates true effect sizes (Delta, RelDiff) for test verification.
            num_datasets = size(all_data{1}, 2);
            pair_idx = nchoosek(1:num_datasets, 2);
            n_pairs = size(pair_idx, 1);
            d_vals = zeros(n_pairs, num_metrics);
            rel_vals = zeros(n_pairs, num_metrics);
            
            for m = 1:num_metrics
                data = all_data{m};
                for k = 1:n_pairs
                    i = pair_idx(k, 1); j = pair_idx(k, 2);
                    x = data(:, i); y = data(:, j);
                    valid = ~isnan(x) & ~isnan(y);
                    x = x(valid); y = y(valid);
                    n = numel(x);
                    if n > 0
                        d_vals(k, m) = HERA.stats.cliffs_delta(x, y);
                        rel_vals(k, m) = HERA.stats.relative_difference(x, y);
                    end
                end
            end
            effect_sizes.d_vals_all = d_vals;
            effect_sizes.rel_vals_all = rel_vals;
        end
        
        function [final_order, final_rank, eff, p_vals] = run_single_test_full(all_data, thresholds, config, names)
            % Wrapper to run 'calculate_ranking' and return full outputs including diagnostics.
            import HERA.test.TestHelper
            
            % Check dimensions
            num_d = size(all_data{1}, 2);
            pairs = nchoosek(1:num_d, 2);
            
            eff = TestHelper.calculate_real_effects(all_data, numel(all_data));
            [final_order, final_rank, ~, ~, p_vals] = calculate_ranking(all_data, eff, thresholds, config, names, pairs);
        end
        
        function [final_rank, final_order] = run_single_test(all_data, thresholds, config, names)
            % Wrapper to run 'calculate_ranking' and return just the order/rank.
            import HERA.test.TestHelper
            [final_order, final_rank, ~, ~] = TestHelper.run_single_test_full(all_data, thresholds, config, names);
        end
        
        function [lang, styles] = get_test_resources()
            % GET_TEST_RESOURCES - Loads production resources for testing.
            %
            % Description:
            %   Loads the actual language pack and design styles used in the application.
            %   This acts as an implicit unit test for the resource loaders themselves.
            %   If this fails, the test suite aborts, which is the intended behavior
            %   as the application cannot run without these resources.
            %   
            %   I know this is a bit of a hack, but printing here all the lang and styles would be
            %   a bit of a mess. Also I want to check that all functions are working as expected when called
            %   directly in the test suite. This is the best way to test that all called functions are using
            %   correctly all the subfunctions.
            %
            % LvE
            
            import HERA.get_language
            import HERA.design
            
            % 1. Load Language (Fatal error if fails)
            lang = get_language();
            
            % 2. Load Styles (Fatal error if fails)
            % Use standard Parameters: Theme='light', N=5 (Arbitrary), Power=true
            styles = design('light', 5, true);
        end
        
    end
end
