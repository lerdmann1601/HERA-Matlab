function results = run_ranking(varargin)
    % Wrapper to expose HERA.run_ranking at the top level of the Python package
    if nargout > 0
        results = HERA.run_ranking(varargin{:});
    else
        HERA.run_ranking(varargin{:});
    end
end
