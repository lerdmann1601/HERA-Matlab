function [selected_B, elbow_indices] = find_elbow_point(x_values, y_values)
% FIND_ELBOW_POINT - Vectorized "Knee/Elbow" detection using Max Distance method.
%
% Syntax:
%   [selected_B, elbow_indices] = HERA.stats.find_elbow_point(x_values, y_values)
%
% Description:
%   Determines the optimal cutoff point (Elbow) for one or multiple stability curves.
%   It calculates the point of maximum perpendicular distance to the line connecting
%   the start and end of the curve.
%
%   Efficiency Note:
%   This implementation is fully vectorized using matrix operations. It processes
%   multiple curves (columns) simultaneously without 'for' loops, which is 
%   significantly faster for high-dimensional stability checks.
%
% Inputs:
%   x_values      - Vector of x-axis values (B-steps). [N x 1]
%   y_values      - Stability metrics. Can be Vector [N x 1], Matrix [N x M], or Cell.
%
% Outputs:
%   selected_B    - The single optimal B (Maximum of all individual elbows).
%   elbow_indices - Indices of the elbow points for each curve [1 x M].
%
% Author: Lukas von Erdmannsdorff

    %% 1. Input Standardization (Ensure Matrix Format)
    x_values = x_values(:); % Ensure Column [N x 1]
    n_points = length(x_values);
    
    % Handle different input types -> Convert to Matrix [N x NumCurves]
    if iscell(y_values)
        % Convert cell array of vectors to matrix (assuming equal length)
        try
            Y = [y_values{:}];
        catch
            % Fallback if lengths differ (rare edge case in HERA) no lang error message jet..
            error('HERA:Elbow:DimensionMismatch', 'Curve lengths in cell array must match X.');
        end
    elseif ismatrix(y_values) && ~isvector(y_values)
        % Matrix: Auto-detect orientation
        [r, c] = size(y_values);
        if r == n_points
            Y = y_values;       % Already [N x M]
        elseif c == n_points
            Y = y_values';      % Transpose to [N x M]
        else
            Y = y_values(:);    % Fallback: Flatten
        end
    else
        % Single Vector
        Y = y_values(:);
    end
    
    % Check for empty or too short data
    if size(Y, 1) < 2
        selected_B = max(x_values);
        elbow_indices = ones(1, size(Y, 2));
        return;
    end

    %% 2. Vectorized Normalization (Min-Max Scaling)
    % We normalize columns independently to [0, 1] range to make distance geometric.
    
    min_vals = min(Y, [], 1);
    max_vals = max(Y, [], 1);
    ranges   = max_vals - min_vals;
    
    % Handle Flat Lines (Range = 0)
    % If a curve is perfectly flat, the "elbow" is conventionally the start (Index 1).
    % We set range to 1 to avoid Div/0, result will be 0 everywhere.
    is_flat = ranges == 0;
    ranges(is_flat) = 1; 
    
    % Normalize Y matrix: [N x M]
    Y_norm = (Y - min_vals) ./ ranges;
    
    % Normalize X vector: [N x 1] -> Replicate to [N x M] via implicit expansion later
    x_min = min(x_values);
    x_range = max(x_values) - x_min;
    if x_range == 0, x_range = 1; end
    X_norm = (x_values - x_min) / x_range;
    
    %% 3. Vectorized Geometric Distance Calculation
    % Line P1 -> P_end
    % Start Points (Normalized): (0, Y_norm(1,:))
    % End Points   (Normalized): (1, Y_norm(end,:))
    
    % Vector connecting Start to End for each curve: [1 x M] (dX is always 1)
    % vec_line = [dX; dY]
    line_vec_x = 1; 
    line_vec_y = Y_norm(end, :) - Y_norm(1, :); % [1 x M]
    
    % Norm of the line vector (Hypotenuse length): sqrt(1^2 + dy^2)
    line_norm = sqrt(line_vec_x^2 + line_vec_y.^2); % [1 x M]
    
    % Vectors from Start Point to Every Point P_i
    % P_i = (X_norm(i), Y_norm(i,:))
    % Vec_P = [X_norm(i) - 0; Y_norm(i,:) - Y_norm(1,:)]
    vec_p_x = X_norm;                     % [N x 1]
    vec_p_y = Y_norm - Y_norm(1, :);      % [N x M]
    
    % 2D Cross Product (Determinant Area)
    % |x1*y2 - x2*y1|
    % Here: vec_p_x * line_vec_y - vec_p_y * line_vec_x
    % Broadcasting: [N x 1] * [1 x M] - [N x M] * 1
    cross_prod = abs(vec_p_x .* line_vec_y - vec_p_y .* line_vec_x); % [N x M]
    
    % Perpendicular Distance = Area / Base_Length
    distances = cross_prod ./ line_norm; % [N x M] ./ [1 x M]
    
    %% 4. Selection
    % Find index of max distance for each column
    [~, elbow_indices] = max(distances, [], 1); % [1 x M]
    
    % Correction for flat lines: If curve was flat, force index 1
    elbow_indices(is_flat) = 1;
    
    % Determine the conservative global B (Maximum of all elbows)
    % We wait until the SLOWEST metric has stabilized.
    if isempty(elbow_indices)
        idx_final = 1;
    else
        idx_final = max(elbow_indices);
    end
    
    selected_B = x_values(idx_final);
end
