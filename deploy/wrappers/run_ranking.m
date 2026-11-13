function results = run_ranking(varargin)
    % Wrapper to expose HERA.run_ranking at the top level of the Python package.
    % Converts non-vector cell arrays in the returned struct (which matlab_pysdk
    % cannot marshal across the C++/Python boundary) into nested 1-by-N cells.
    if nargout > 0
        results = HERA.run_ranking(varargin{:});
        results = sanitize_for_python(results);
    else
        HERA.run_ranking(varargin{:});
    end
end

function val = sanitize_for_python(val)
    if isstruct(val)
        fnames = fieldnames(val);
        for idx = 1:numel(val)
            for f = 1:numel(fnames)
                fn = fnames{f};
                val(idx).(fn) = sanitize_for_python(val(idx).(fn));
            end
        end
    elseif iscell(val)
        for idx = 1:numel(val)
            val{idx} = sanitize_for_python(val{idx});
        end
        sz = size(val);
        if ndims(val) > 2
            nRows = sz(1);
            out = cell(1, nRows);
            otherDims = sz(2:end);
            for r = 1:nRows
                slice = reshape(val(r, :), otherDims);
                out{r} = sanitize_for_python(slice);
            end
            val = out;
        elseif sz(1) ~= 1 && sz(2) ~= 1
            if sz(1) == 0 || sz(2) == 0
                val = cell(1, 0);
            else
                nRows = sz(1);
                out = cell(1, nRows);
                for r = 1:nRows
                    out{r} = val(r, :);
                end
                val = out;
            end
        end
    end
end
