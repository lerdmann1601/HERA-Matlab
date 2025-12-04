function build_hera()
% BUILD_HERA - Compiles HERA and creates an Installer with auto-runtime download.
%
% Syntax:
%   build_hera()
%
% Description:
%   This script automates the deployment process using the 'compiler.build' workflow.
%   It performs two main tasks:
%   1.  Compiles the HERA toolbox into a standalone executable.
%   2.  Packages the executable into a web-based installer that automatically downloads the required MATLAB Runtime during installation.
%
% Workflow:
%   1.  Initialization: Detects the project root and sets up output directories.
%   2.  Configuration: Defines build options, including the main entry point and resource folders.
%   3.  Compilation: Builds the standalone application using 'compiler.build.standaloneApplication'.
%   4.  Packaging: Creates an installer using 'compiler.package.installer' with 'RuntimeDelivery' set to 'web'.
%
% Inputs:
%   None
%
% Outputs:
%   The compiled application and installer are saved in 'deploy/output'.
%   Console output indicates build progress and status.
%
% Author: Lukas von Erdmannsdorff

    %% 1. Initialization and Path Detection
    % Determine the project root directory relative to this script.
    % This script is located in <ProjectRoot>/deploy/build_hera.m
    
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
    outputDir = fullfile(projectRoot, 'deploy', 'output');
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end

    % Define Main File
    mainFile = fullfile(projectRoot, '+HERA', 'start_ranking.m');
    
    % Define Application Name
    appName = 'HERA_Runtime';

    % Define Resources to Include
    % Note: +HERA package is automatically analyzed by the compiler, 
    % but explicit assets/languages need to be added.
    additionalFiles = [ ...
        string(fullfile(projectRoot, 'assets')), ...
        string(fullfile(projectRoot, '+HERA', 'language')), ... 
        string(fullfile(projectRoot, 'paper')) ...
    ];

    %% 2. Build Configuration
    fprintf('Configuring Build Options...\n');
    
    % Initialize options for standalone application
    buildOpts = compiler.build.StandaloneApplicationOptions(mainFile);
    buildOpts.ExecutableName = appName;
    buildOpts.OutputDir = outputDir;
    buildOpts.Verbose = true;
    buildOpts.TreatInputsAsNumeric = false; % Adjust if arguments are passed from CLI
    
    % Add resource folders
    buildOpts.AdditionalFiles = additionalFiles;

    %% 3. Compilation and Packaging
    fprintf('========================================\n');
    fprintf('Starting HERA Build for Platform: %s\n', computer);
    fprintf('========================================\n');

    try
        % Step 1: Compile the Application
        fprintf('1. Compiling Application (may take a few minutes)...\n');
        buildResults = compiler.build.standaloneApplication(buildOpts);
        
        % Step 2: Create the Installer (Packaging)
        fprintf('2. Creating Installer with Auto-Runtime Download...\n');
        
        % Use Name-Value pairs directly to avoid potential object handling issues
        installerName = [appName '_Installer'];
        
        compiler.package.installer(buildResults, ...
            'InstallerName', installerName, ...
            'OutputDir', outputDir, ...
            'RuntimeDelivery', 'web');
        
        % Step 3: Copy Launcher Scripts to Output (macOS only)
        if ismac
            fprintf('3. Copying Launcher Scripts...\n');
            launcherSrc = fullfile(projectRoot, 'release', 'macos', 'HERA_Launcher.command');
            launcherDst = fullfile(outputDir, 'HERA_Launcher.command');
            copyfile(launcherSrc, launcherDst);
        end
        
        % Success Message
        fprintf('\n========================================\n');
        fprintf('SUCCESS!\n');
        fprintf('1. Standalone App: %s\n', fullfile(outputDir, appName));
        fprintf('2. Installer:      %s\n', fullfile(outputDir, [appName '_Installer']));
        if ismac
            fprintf('3. Launcher:       %s\n', launcherDst);
        end
        fprintf('----------------------------------------\n');
        fprintf('Share the Installer with users.\n');
        if ismac
            fprintf(' Also share the Launcher script (needed for Terminal on macOS).\n');
        end
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
        if exist('buildResults', 'var')
            disp('Build Results:');
            disp(buildResults);
        end
        
        % Rethrow the error so the CI/CD pipeline knows the build failed
        rethrow(ME);
    end
end