function check_updates(lang)
% CHECK_UPDATES - Verifies if a new HERA version is available on GitHub.
%
% Syntax:
%   HERA.start.check_updates(lang)
%
% Description:
%   Queries the GitHub API for the latest release tag and compares it with 
%   the local version. Notifies the user in the console if an update exists.
%   It runs only once per session, skips on clusters, and never blocks due to a 2-second timeout.
%
% Inputs:
%   lang - (struct) Language structure for localized notifications.
%
% Author: Lukas von Erdmannsdorff

    % 1. Persistent flag to prevent multiple checks in one session
    persistent hasChecked;
    if ~isempty(hasChecked) && hasChecked
        return;
    end
    hasChecked = true;

    % 2. Safety Checks: Skip if not appropriate
    % Skip if we are running as a parallel worker or in a batch job to avoid overhead
    try
        if ~isempty(getCurrentWorker()) || ~isempty(getCurrentJob())
            return;
        end
    catch
        % Parallel Computing Toolbox not installed or configured
    end

    % Check if we are in a non-interactive/headless environment that doesn't 
    % need update notifications (e.g. CI/CD)
    if ~feature('ShowFigureWindows') || ~isempty(getenv('GITHUB_ACTIONS'))
        return;
    end

    % 3. Check for Internet and Fetch Latest Version
    try
        % Define GitHub API URL
        url = 'https://api.github.com/repos/lerdmann1601/HERA-Matlab/releases/latest';
        
        % Set a strict timeout to ensure we never block the user
        options = weboptions('Timeout', 2); 
        
        % Fetch JSON
        data = webread(url, options);
        
        if isfield(data, 'tag_name')
            latest_version = data.tag_name;
            current_version = HERA.get_version();
            
            % Simple comparison (exact match). 
            % Note: For more complex semver, we could parse the numbers, 
            % but exact mismatch is usually sufficient to trigger a notice.
            if ~strcmp(latest_version, current_version)
                fprintf('\n--------------------------------------------------\n');
                fprintf('  %s\n', lang.update_check.notice);
                fprintf('  %s\n', sprintf(lang.update_check.current, current_version, latest_version));
                fprintf('  %s\n', lang.update_check.download);
                fprintf('--------------------------------------------------\n\n');
            end
        end
    catch
        % Silently fail on network issues, timeouts, or API limits.
        % The update check is a "nice-to-have" and should never cause an error.
    end
end
