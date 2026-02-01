function package_HERA_toolbox()
% PACKAGE_HERA_TOOLBOX - Packages HERA as a MATLAB Toolbox (.mltbx).
%
% Syntax:
%   package_HERA_toolbox()
%
% Description:
%   This script automates the creation of a MATLAB Toolbox file (.mltbx)
%   for easy installation and distribution via the File Exchange.
%   It uses 'matlab.addons.toolbox.packageToolbox'.
%
% Workflow:
%   1.  Initialization: Detects the project root and sets up output directories.
%   2.  Configuration: Defines toolbox metadata (Version, Name, Author, etc.) and included files.
%   3.  Packaging: Generates the .mltbx file in 'deploy/output/toolbox'.
%
% Inputs:
%   None
%
% Outputs:
%   The .mltbx file is saved in 'deploy/output/toolbox'.
%   Console output indicates packaging progress.
%
% Author: Lukas von Erdmannsdorff

    clc
    %% 1. Initialization and Path Detection
    fprintf('Initializing Toolbox Packaging...\n');

    % Get the full path of this script
    scriptPath = mfilename('fullpath');
    
    % Get the directory containing this script (the 'deploy' folder)
    deployDir = fileparts(scriptPath);
    
    % Get the parent directory (the Project Root)
    projectRoot = fileparts(deployDir);
    
    % Verify that the +HERA package exists in the root
    if ~exist(fullfile(projectRoot, '+HERA'), 'dir')
        error('Could not locate +HERA package. Expected at: %s', fullfile(projectRoot, '+HERA'));
    end
    
    % Define and create the output directory.
    outputDir = fullfile(projectRoot, 'deploy', 'output', 'toolbox');
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end

    % Get Version
    version_str = HERA.get_version();
    fprintf('Detected Version: %s\n', version_str);

    % Clean up problematic files before packaging
    % These files can cause File Exchange upload validation to fail
    fprintf('Cleaning up temporary/unwanted files...\n');
    
    % Patterns to delete (relative to projectRoot)
    cleanupPatterns = {'**/__pycache__', '**/*.pyc', '**/.DS_Store'};
    
    totalDeleted = 0;
    for i = 1:length(cleanupPatterns)
        pattern = cleanupPatterns{i};
        files = dir(fullfile(projectRoot, pattern));
        for j = 1:length(files)
            itemPath = fullfile(files(j).folder, files(j).name);
            
            % Skip .venv, .git, deploy/output, deploy/dist, and release folders
            if contains(itemPath, fullfile(projectRoot, '.venv')) || ...
               contains(itemPath, fullfile(projectRoot, '.git')) || ...
               contains(itemPath, fullfile(projectRoot, 'deploy', 'output')) || ...
               contains(itemPath, fullfile(projectRoot, 'deploy', 'dist')) || ...
               contains(itemPath, fullfile(projectRoot, 'release'))
               continue;
            end

            try
                if files(j).isdir
                    rmdir(itemPath, 's');
                else
                    delete(itemPath);
                end
                totalDeleted = totalDeleted + 1;
                % Only print if it's a file inside project sources to reduce noise
                if ~contains(itemPath, '.venv')
                    fprintf('  Removed: %s\n', strrep(itemPath, projectRoot, ''));
                end
            catch
                % Suppress warnings for expected conflicts
            end
        end
    end
    fprintf('Cleanup complete. Removed %d items.\n', totalDeleted);

    %% 2. Toolbox Configuration
    fprintf('Configuring Toolbox Options...\n');
    
    % Define Toolbox Name
    toolboxName = 'HERA';
    
    % Create ToolboxOptions object
    % syntax: ToolboxOptions(toolboxPath, toolboxName)
    % We point it to projectRoot as the source of files
    opts = matlab.addons.toolbox.ToolboxOptions(projectRoot, toolboxName);
    
    % Metadata
    opts.ToolboxName = 'HERA'; % The display name
    opts.ToolboxVersion = replace(version_str, 'v', '');
    opts.AuthorName = 'Lukas von Erdmannsdorff';
    opts.AuthorEmail = ''; % Optional: Leave empty or fill if known
    
    % Description
    opts.Summary = 'HERA: A High-Efficiency Ranking Tool for MCDA';
    opts.Description = 'HERA provides tools for Multi-Criteria Decision Analysis, including ranking, visualization, and validation.';
    
    % Define essential files and directories to include in the toolbox.
    % Explicitly listing items prevents the inclusion of development artifacts (e.g., .git, tests).
    % Note: ToolboxOptions automatically handles +NAmespace folders (like +HERA) if listed.
    
    includedPaths = { ...
        fullfile(projectRoot, '+HERA'), ...
        fullfile(projectRoot, 'assets'), ...
        fullfile(projectRoot, 'data'), ...
        fullfile(projectRoot, 'deploy', 'build_HERA_matlab.m'), ...
        fullfile(projectRoot, 'deploy', 'build_HERA_python.m'), ...
        fullfile(projectRoot, 'deploy', 'build_and_prep_pypi.sh'), ...
        fullfile(projectRoot, 'deploy', 'package_HERA_toolbox.m'), ...
        fullfile(projectRoot, 'deploy', 'readme.txt'), ...
        fullfile(projectRoot, 'deploy', 'dependencies'), ...
        fullfile(projectRoot, 'deploy', 'python_assets'), ...
        fullfile(projectRoot, 'deploy', 'wrappers'), ...
        fullfile(projectRoot, 'docs'), ...
        fullfile(projectRoot, 'paper'), ... 
        fullfile(projectRoot, 'tests', 'run_HERA_tests.m'), ...
        fullfile(projectRoot, 'tests', 'HERA_Validation_Example.txt'), ...
        fullfile(projectRoot, 'CITATION.cff'), ...
        fullfile(projectRoot, 'CODE_OF_CONDUCT.md'), ...
        fullfile(projectRoot, 'CONTRIBUTING.md'), ...
        fullfile(projectRoot, 'license.txt'), ...
        fullfile(projectRoot, 'README.md'), ...
        fullfile(projectRoot, 'setup_HERA.m') ...
    };

    % Filter to ensuring they exist before adding
    validFiles = {};
    for i = 1:length(includedPaths)
        if exist(includedPaths{i}, 'file') || exist(includedPaths{i}, 'dir')
             validFiles{end+1} = includedPaths{i};
        else
             fprintf('Warning: Resource not found, skipping: %s\n', includedPaths{i});
        end
    end
    
    if isempty(validFiles)
        error('No valid files found to include in the toolbox.');
    end
    
    opts.ToolboxFiles = validFiles;

    % MATLAB Path Setup
    % By default, the root of the installed toolbox is added to the path.
    % That is sufficient for +HERA to work.
    
    opts.OutputFile = fullfile(outputDir, [toolboxName '_' version_str '.mltbx']);

    %% 3. Packaging
    fprintf('==================================\n');
    fprintf('Starting HERA Toolbox Packaging...\n');
    fprintf('==================================\n');

    try
        fprintf('Creating .mltbx file...\n');
        matlab.addons.toolbox.packageToolbox(opts);
        
        fprintf('\n========================================\n');
        fprintf('SUCCESS!\n');
        fprintf('Toolbox File: %s\n', opts.OutputFile);
        fprintf('----------------------------------------\n');
        fprintf('Ready for Release! Upload the .mltbx file to GitHub/FileExchange.\n');
        fprintf('========================================\n');
        
    catch ME
        % Error Handling
        fprintf('\n========================================\n');
        fprintf('Packaging Failed:\n%s\n', ME.message);
        fprintf('Stack Trace:\n');
        for k = 1:length(ME.stack)
            fprintf('  File: %s\n  Name: %s\n  Line: %d\n', ME.stack(k).file, ME.stack(k).name, ME.stack(k).line);
        end
        fprintf('========================================\n');
        rethrow(ME);
    end
end
