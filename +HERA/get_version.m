function version_str = get_version()
% GET_VERSION - Retrieves the current version of the HERA toolbox.
%
%   version_str = HERA.get_version()
%
% Description:
%   Attempts to retrieve the version string from the Git repository using
%   'git describe --tags --always --dirty'. If Git is not available or
%   execution fails, it falls back to a default version string.
%
% Outputs:
%   version_str - (char) The version string (e.g., 'v1.0.0-27-gec99cb4-dirty' or '2.1.0-dev').
%
% Author: Lukas von Erdmannsdorff

    % Default fallback version if Git detection fails
    fallback_version = '1.0.0'; 

    try
        % Attempt to run git describe
        [status, cmdout] = system('git describe --tags --always --dirty');
        
        if status == 0
            version_str = strtrim(cmdout);
        else
            % Git failed (nonzero exit code), use fallback
            version_str = fallback_version;
        end
    catch
        % System command failed (e.g., permissions), use fallback
        version_str = fallback_version;
    end
end
