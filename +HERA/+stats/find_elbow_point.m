function [selected_B, elbow_indices] = find_elbow_point(x_values, y_values)
% FIND_ELBOW_POINT - Uses the "Knee/Elbow Point" method to find the optimal B.
%
% Syntax:
%   [selected_B, elbow_indices] = find_elbow_point(x_values, y_values)
%
% Description:
%   This function implements the "Maximum Distance to Line" method (often called the Elbow 
%   or Knee method) to determine the optimal cutoff point in a curve.
%   It calculates the perpendicular distance from each point on the curve to the line 
%   connecting the first and last points. The point with the maximum distance is selected.
%
%   This is useful for BCa bootstrapping to find a B where the stability gain diminishes.
%   The function handles single curves (vector) or multiple curves (matrix/cell) simultaneously.
%   If multiple curves are provided, it calculates the elbow for each and returns the 
%   maximum B among them (conservative approach).
%
% Inputs:
%   x_values      - Vector of x-axis values (e.g., the tested B steps).
%   y_values      - The dependent values (e.g., stability metrics). Can be:
%                   * Single vector.
%                   * Matrix (Rows = X-steps, Cols = Curves OR Cols = X-steps, Rows = Curves). 
%                     (Function attempts to auto-detect orientation).
%                   * Cell Array of vectors (each cell is a curve).
%
% Outputs:
%   selected_B    - The single optimal x-value (B) selected. Max of all individual elbows.
%   elbow_indices - Vector containing the index of the elbow point for each curve analyzed.
%                   Useful for plotting/visualization.
%
% Robustness:
%   - Handles NaNs by forcing normalized values to 0.
%   - handles "Flat" or "No Variance" curves by selecting the first index.
%   - Normalizes both X and Y axes to [0,1] to ensure geometric distance is scale-invariant.
%
% Author: Lukas von Erdmannsdorff

    x_values = x_values(:); % Ensure x is a column vector
    
    %% 1. Input Standardization
    % Convert all valid input formats (Vector, Matrix, Cell) into a Cell Array of column vectors
    if iscell(y_values)
        curves = y_values;
    elseif ismatrix(y_values) && ~isvector(y_values)
        % Matrix input: Determine orientation based on x_values length
        nx = length(x_values);
        [r, c] = size(y_values);
        
        if r == nx
            % Rows match X: Columns are curves
            curves = num2cell(y_values, 1);
        elseif c == nx
             % Columns match X: Rows are curves -> Transpose
             curves = num2cell(y_values', 1);
        else
             % Ambiguous or Mismatch: Treat as single flattened curve (Fallback)
             % This mimics original behavior if dimension check failed somehow
             curves = {y_values(:)}; 
        end
    else
        % Single vector input
        curves = {y_values(:)};
    end
    
    num_curves = numel(curves);
    elbow_indices = zeros(1, num_curves);
    
    %% 2. Elbow Detection Loop
    for k = 1:num_curves
        y = curves{k};
        
        % Safety Check: Skip if curve has no valid data or no variation (flat line)
        % If flat, practically 'stable' from the start -> index 1
        if any(isnan(y)) || numel(unique(y)) < 2
            elbow_indices(k) = 1;
            continue;
        end
        
        % Normalization to [0, 1] range
        % This is crucial because X (e.g., 1000s) and Y (e.g., 0.01s) have vastly different scales.
        x_norm = (x_values - min(x_values)) / (max(x_values) - min(x_values));
        y_norm = (y - min(y)) / (max(y) - min(y));
        
        % Robust handling for NaN results from normalization (e.g., if max == min despite unique check)
        if all(isnan(x_norm)), x_norm(:) = 0; end
        if all(isnan(y_norm)), y_norm(:) = 0; end
        
        % Vector Geometry: Distance from point P to line AB (Start to End)
        % Line Vector (Vector connecting first and last point)
        line_vec = [x_norm(end) - x_norm(1); y_norm(end) - y_norm(1)];
        
        % Vector from First Point to every other point
        vec_from_first = [x_norm - x_norm(1), y_norm - y_norm(1)];
        
        % 2D Cross Product (Determinant) represents the area of parallelogram.
        % Area / Base_Length = Height (Perpendicular Distance)
        cross_prod = vec_from_first(:, 1) * line_vec(2) - vec_from_first(:, 2) * line_vec(1);
        
        % Find index of maximum distance
        [~, idx] = max(abs(cross_prod) / norm(line_vec));
        
        % Fallback: If calculation fails (e.g., empty), default to last index
        if isempty(idx), idx = numel(x_values); end
        
        elbow_indices(k) = idx;
    end
    
    %% 3. Result Selection
    % We choose the MAXIMUM B among all curves.
    % Conservative strategy: We wait until the SLOWEST metric has stabilized.
    selected_B = max(x_values(elbow_indices));
end
