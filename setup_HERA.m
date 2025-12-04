function setup_HERA()
% SETUP_HERA - Sets up the MATLAB path for the HERA Toolbox.
%
% Description:
%   This script adds the current project directory to the MATLAB path.
%   It ensures that the '+HERA' package is accessible from anywhere.
%   Run this script once after cloning the repository.
%
% Usage:
%   setup_HERA()
%
% Author: Lukas von Erdmannsdorff

    fprintf('Initializing HERA Toolbox...\n');

    % 1. Identify the project root (directory of this script)
    projectRoot = fileparts(mfilename('fullpath'));

    % 2. Add the project root to the MATLAB path
    % We do not use genpath() to avoid adding .git or test folders to the path, which prevents namespace pollution. 
    % MATLAB finds '+HERA' automatically in root.
    addpath(projectRoot);

    % 3. Save the path (optional, user might not have write permissions)
    try
        savepath;
        fprintf(' -> HERA added to path and saved permanently.\n');
    catch
        fprintf(' -> HERA added to path for this session.\n');
        fprintf('    (Could not save path permanently. You may need to run this again next time.)\n');
    end

    fprintf('Setup complete. You can now run "HERA.start_ranking" from anywhere.\n');
end
