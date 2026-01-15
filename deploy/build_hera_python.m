function build_HHERA_python()
% BUILD_HERA_PYTHON - Compiles HERA into a Python Package.
%
% Syntax:
%   build_hera_python()
%
% Description:
%   This script automates the creation of a Python package for HERA using 
%   the 'compiler.build.pythonPackage' workflow.
%   It allows HERA to be installed via pip and used in Python environments,
%   MATLAB Runtime should be installed manually via 
%   https://www.mathworks.com/products/compiler/matlab-runtime.html.
%
% Workflow:
%   1.  Initialization: Detects the project root and sets up output directories.
%   2.  Configuration: Defines build options, including packaging name and resources.
%   3.  Compilation: Builds the Python package using 'compiler.build.pythonPackage'.
%
% Inputs:
%   None
%
% Outputs:
%   The compiled Python package is saved in 'deploy/output/python'.
%
% Author: Lukas von Erdmannsdorff
%

clc
    %% 1. Initialization and Path Detection
    % Determine the project root directory relative to this script.
    % This script is located in <ProjectRoot>/deploy/build_hera_python.m
    
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
    outputDir = fullfile(projectRoot, 'deploy', 'output', 'python');
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end

    % Define Package Name (This becomes the python import name)
    pkgName = 'hera_matlab';
    
    % Get Version
    version_str = replace(HERA.get_version(), 'v', '');
    fprintf('Detected Version: %s\n', version_str);

    % Define Resources to Include
    % Note: +HERA package is automatically analyzed by the compiler, 
    % but explicit assets/languages need to be added.
    additionalFiles = [ ...
        string(fullfile(projectRoot, 'assets')), ...
        string(fullfile(projectRoot, '+HERA', 'language')), ... 
        string(fullfile(projectRoot, 'paper')) ...
    ];

    %% 2. Build Configuration
    fprintf('Configuring Build Options for Python Package...\n');
    
    % Define the primary entry points for export.
    % Explicitly selecting functions such as 'start_ranking' and 'run_ranking' ensures
    % a clean and defined API for the generated Python package.
    % This method is preferred over exporting the entire '+HERA' namespace folder directly.
    % We likely want users to call `hera_matlab.start_ranking()` or `hera_matlab.HERA`.
    exportedFunctions = [ ...
        string(fullfile(projectRoot, '+HERA', 'start_ranking.m')), ...
        string(fullfile(projectRoot, '+HERA', 'run_ranking.m')) ...
    ];
    
    % Initialize options
    buildOpts = compiler.build.PythonPackageOptions(exportedFunctions);
    buildOpts.PackageName = pkgName;
    buildOpts.OutputDir = outputDir;
    buildOpts.Verbose = true;
    
    % Add resource folders
    buildOpts.AdditionalFiles = additionalFiles;

    %% 3. Compilation
    fprintf('========================================\n');
    fprintf('Starting HERA Python Package Build...\n');
    fprintf('Package Name: %s\n', pkgName);
    fprintf('========================================\n');

    try
        % Step 1: Compile the Package
        fprintf('1. Compiling Python Package (may take a few minutes)...\n');
        buildResults = compiler.build.pythonPackage(buildOpts);
        
    % Step 2: Create the Installer (Optional - Currently Disabled)
    % We focus on PyPI distribution where the user installs the Runtime manually.
    % fprintf('\n2. Creating Installer with Auto-Runtime Download...\n');
    % installerName = [pkgName '_Installer_' version_str];
    % compiler.package.installer(buildResults, ...
    %     'InstallerName', installerName, ...
    %     'OutputDir', outputDir, ...
    %     'RuntimeDelivery', 'web');

    % Success Message
    fprintf('========================================\n');
    fprintf('SUCCESS!\n');
    % fprintf('1. Python Package: %s\n', fullfile(outputDir, pkgName));
    % fprintf('2. Installer:      %s\n', fullfile(outputDir, installerName));
    fprintf('Python Package generated at:\n%s\n', fullfile(outputDir, pkgName));
    fprintf('----------------------------------------\n');
    fprintf('To release:\n');
    fprintf('1. Verify: cd %s && pip install .\n', fullfile(outputDir, pkgName));
    fprintf('2. Upload: twine upload dist/* (inside the package folder)\n');
    fprintf('========================================\n');
        
    catch ME
        % Error Handling
        fprintf('\n========================================\n');
        fprintf('Build Failed:\n%s\n', ME.message);
        fprintf('Stack Trace:\n');
        for k = 1:length(ME.stack)
            fprintf('  File: %s\n  Name: %s\n  Line: %d\n', ME.stack(k).file, ME.stack(k).name, ME.stack(k).line);
        end
        fprintf('========================================\n');
        
        % Rethrow the error
        rethrow(ME);
    end
end
