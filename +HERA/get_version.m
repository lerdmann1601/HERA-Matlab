function version_str = get_version()
% GET_VERSION - Retrieves the current version of the HERA toolbox.
%
% Syntax:
%   version_str = HERA.get_version()
%
% Description:
%   Attempts to retrieve the version string from the Git repository using 'git describe --tags --always --dirty'. 
%   This provides improved traceability by including the latest tag, commit offset, and dirty state.
%
%   If the Git command fails (e.g., Git not installed, not a repository, 
%   or no tags found), it falls back to defined default version. 
%   The fallback version is currently defined as 'v1.1.1'. 
%   Please make sure to update this version if a new tag is created!
%
% Outputs:
%   version_str - (char) The version string (e.g., 'v1.1.1').
%
% Author: Lukas von Erdmannsdorff

    % Default fallback version used if dynamic retrieval fails
    fallback_version = 'v1.1.1'; 

    try
        % Attempt to run git describe.
        % We use '2>&1' to redirect stderr to stdout. 
        % This prevents Git error messages (like "fatal: not a git repository") from appearing in the 
        % MATLAB Command Window and allows us to handle them logic-wise.
        [status, cmdout] = system('git describe --tags --abbrev=0 2>&1');
        
        % Check if the command was successful (status 0) and outcome is valid
        if status == 0 && ~isempty(cmdout)
            version_str = strtrim(cmdout);
        else
            % Fallback if Git returns non-zero status or empty output
            version_str = fallback_version;
        end
    catch
        % Fallback if the system command itself fails (e.g., permissions)
        version_str = fallback_version;
    end
end
