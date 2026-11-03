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

% Get Version and validate (Strict for CI/CD)
is_cicd = ~isempty(getenv('GITHUB_ACTIONS'));
if is_cicd
    version_str = getenv('GITHUB_REF_NAME');
    if isempty(version_str) || ~startsWith(version_str, 'v')
        error('Error: Could not determine a valid version for the build. In CI/CD environments, a Git Tag (e.g., v1.4.8) must be set via GITHUB_REF_NAME.');
    end
else
    % Local build: Use get_version() helper
    version_str = HERA.get_version();
end
fprintf('Detected Version: %s\n', version_str);

% Clean up problematic files before packaging
% These files can cause File Exchange upload validation to fail
fprintf('Cleaning up temporary/unwanted files...\n');

% Clean up temporary Google Drive sync directories if present
driveTempDirs = {fullfile(projectRoot, '.tmp.driveupload'), fullfile(projectRoot, '.tmp.drivedownload')};
for k = 1:length(driveTempDirs)
    if exist(driveTempDirs{k}, 'dir')
        try
            rmdir(driveTempDirs{k}, 's');
        catch
        end
    end
end

% Patterns and target directories to clean (avoids scanning .git, .venv, release, or data/results)
cleanupPatterns = {'**/__pycache__', '**/*.pyc', '**/.DS_Store'};
cleanupTargets = { ...
    fullfile(projectRoot, '+HERA'), ...
    fullfile(projectRoot, 'data', 'examples'), ...
    fullfile(projectRoot, 'tests'), ...
    fullfile(projectRoot, 'deploy', 'wrappers'), ...
    fullfile(projectRoot, 'deploy', 'python_assets'), ...
    fullfile(projectRoot, '.github') ...
    };

totalDeleted = 0;

% Clean root-level .DS_Store files directly
rootDsStores = {fullfile(projectRoot, '.DS_Store'), fullfile(projectRoot, 'data', '.DS_Store')};
for k = 1:length(rootDsStores)
    if exist(rootDsStores{k}, 'file')
        try
            delete(rootDsStores{k});
            totalDeleted = totalDeleted + 1;
            fprintf('  Removed: %s\n', strrep(rootDsStores{k}, projectRoot, ''));
        catch
        end
    end
end

% Clean target directories recursively
for t = 1:length(cleanupTargets)
    targetDir = cleanupTargets{t};
    if ~exist(targetDir, 'dir')
        continue;
    end
    for i = 1:length(cleanupPatterns)
        pattern = cleanupPatterns{i};
        files = dir(fullfile(targetDir, pattern));
        for j = 1:length(files)
            itemPath = fullfile(files(j).folder, files(j).name);
            try
                if exist(itemPath, 'dir')
                    rmdir(itemPath, 's');
                    totalDeleted = totalDeleted + 1;
                    fprintf('  Removed: %s\n', strrep(itemPath, projectRoot, ''));
                elseif exist(itemPath, 'file')
                    delete(itemPath);
                    totalDeleted = totalDeleted + 1;
                    fprintf('  Removed: %s\n', strrep(itemPath, projectRoot, ''));
                end
            catch
                % Suppress warnings for expected conflicts
            end
        end
    end
end
fprintf('Cleanup complete. Removed %d items.\n', totalDeleted);

%% 2. Toolbox Configuration
fprintf('Configuring Toolbox Options...\n');

% Define Toolbox Name
toolboxName = 'HERA';

% Define essential files and directories to include in the toolbox (relative to projectRoot).
% Explicitly listing items prevents the inclusion of development artifacts (e.g., .git, tests).
% Note: ToolboxOptions automatically handles +Namespace folders (like +HERA) if listed.

% NOTE: Exclude 'docs', 'assets', 'deploy', 'data/results', and 'data/utils' to keep the toolbox lean.
% Internal deployment scripts, raw test files, benchmark results ('data/results'), and development utilities ('data/utils') are excluded to ensure a professional production distribution.
relativePaths = { ...
    '+HERA', ...
    fullfile('data', 'README.md'), ...
    fullfile('data', 'examples'), ...
    fullfile('tests', 'HERA_Validation_Example.txt'), ...
    'CITATION.cff', ...
    'CODE_OF_CONDUCT.md', ...
    'CONTRIBUTING.md', ...
    'license.txt', ...
    'README.md', ...
    'setup_HERA.m' ...
    };

% Stage valid toolbox files in a clean temporary directory.
% Pointing ToolboxOptions to an isolated staging directory prevents MATLAB from
% recursively scanning heavy non-toolbox folders (.git, .venv, release, data/results)
% or stalling on cloud-synced workspace folders.
stagingDir = tempname;
mkdir(stagingDir);
cleanupObj = onCleanup(@() remove_staging_dir(stagingDir));

validFiles = {};
for i = 1:length(relativePaths)
    relPath = relativePaths{i};
    srcPath = fullfile(projectRoot, relPath);
    dstPath = fullfile(stagingDir, relPath);

    % Explicit guard: ensure 'data/results' and 'data/utils' are always ignored
    if contains(relPath, fullfile('data', 'results')) || contains(relPath, fullfile('data', 'utils'))
        continue;
    end

    if exist(srcPath, 'dir')
        if ~exist(dstPath, 'dir')
            mkdir(dstPath);
        end
        copyfile(srcPath, dstPath);
        validFiles{end+1} = dstPath;
    elseif exist(srcPath, 'file')
        dstParent = fileparts(dstPath);
        if ~exist(dstParent, 'dir')
            mkdir(dstParent);
        end
        copyfile(srcPath, dstPath);
        validFiles{end+1} = dstPath;
    else
        fprintf('Warning: Resource not found, skipping: %s\n', srcPath);
    end
end

if isempty(validFiles)
    error('No valid files found to include in the toolbox.');
end

% Ensure no residual temporary files (.DS_Store, __pycache__, .pyc) exist in the staging directory
for i = 1:length(cleanupPatterns)
    stagedCleanup = dir(fullfile(stagingDir, cleanupPatterns{i}));
    for j = 1:length(stagedCleanup)
        stagedItem = fullfile(stagedCleanup(j).folder, stagedCleanup(j).name);
        try
            if exist(stagedItem, 'dir')
                rmdir(stagedItem, 's');
            elseif exist(stagedItem, 'file')
                delete(stagedItem);
            end
        catch
        end
    end
end

% Create ToolboxOptions object from the clean staging directory
% syntax: ToolboxOptions(toolboxPath, toolboxUUID)
% UUID that identifies this toolbox on File Exchange
packageId = "991899d7-67cb-40a5-a684-e51f9681949b";
opts = matlab.addons.toolbox.ToolboxOptions(stagingDir, packageId);

% Metadata
opts.ToolboxName = 'HERA'; % The display name
opts.ToolboxVersion = replace(version_str, 'v', '');
opts.AuthorName = 'Lukas von Erdmannsdorff';
opts.AuthorEmail = ''; % Optional: Leave empty or fill if known

% Description
opts.Summary = 'A Scientific Ranking Framework for Paired Benchmarking';
opts.Description = ['HERA (Hierarchical-Compensatory, Effect-Size-Driven Ranking Algorithm) is a MATLAB toolbox designed to automate the objective comparison of algorithms, experimental conditions, or datasets across up to three quality metrics. ' ...
    'It implements a hierarchical-compensatory logic integrating non-parametric significance testing, robust effect size estimation (Cliff''s Delta and Relative Mean Difference), and bootstrapping to produce data-driven thresholds and statistically robust rankings. ' ...
    'For more information, visit the project website: https://lerdmann1601.github.io/HERA-Matlab/'];

% Compatibility (Minimum MATLAB Release)
opts.MinimumMatlabRelease = 'R2020a';
opts.MaximumMatlabRelease = '';

opts.ToolboxFiles = validFiles;

% Defensive guard: ensure no files from 'data/results' or 'data/utils' can be packaged
opts.ToolboxFiles(contains(opts.ToolboxFiles, fullfile('data', 'results')) | ...
                  contains(opts.ToolboxFiles, fullfile('data', 'utils'))) = [];

% Set Toolbox Image (for File Exchange and Add-On Explorer thumbnail)
logoSrc = fullfile(projectRoot, 'release', 'hera_logo_submission.png');
if exist(logoSrc, 'file')
    logoDst = fullfile(stagingDir, 'hera_logo_submission.png');
    copyfile(logoSrc, logoDst);
    opts.ToolboxImageFile = logoDst;
    fprintf('Configured Toolbox Image: %s\n', logoSrc);
else
    fprintf('Warning: Toolbox logo image not found at: %s\n', logoSrc);
end

% MATLAB Path Setup
% By default, only the root of the installed toolbox is added to the path.
% That is sufficient for +HERA to work and keeps the user's MATLAB path clean.
opts.ToolboxMatlabPath = {stagingDir};

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

%% Local Helper Functions
function remove_staging_dir(stagingDir)
% REMOVE_STAGING_DIR - Safely removes the temporary staging directory.
if exist(stagingDir, 'dir')
    try
        rmdir(stagingDir, 's');
    catch
        % Suppress warnings if OS still holds a transient lock
    end
end
end
