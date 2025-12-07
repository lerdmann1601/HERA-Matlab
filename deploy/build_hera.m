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
    
    % Get Version
    version_str = HERA.get_version();
    fprintf('Detected Version: %s\n', version_str);

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
        installerName = [appName '_Installer_' version_str];
        
        compiler.package.installer(buildResults, ...
            'InstallerName', installerName, ...
            'OutputDir', outputDir, ...
            'RuntimeDelivery', 'web');
        
        % Step 3: Copy Launcher Scripts to Output (macOS only)
        if ismac
            fprintf('3. Copying Launcher Scripts...\n');
            launcherSrc = fullfile(projectRoot, 'deploy', 'dependencies', 'HERA_Launcher.command');
            launcherDst = fullfile(outputDir, 'HERA_Launcher.command');
            copyfile(launcherSrc, launcherDst);
            % Ensure it is executable
            system(['chmod +x "' launcherDst '"']);
        end
        
        
        % Step 4: Compress Artifacts (ZIP)
        fprintf('4. Compressing Artifacts into ZIP...\n');
        
        % Define Zip Name (e.g., HERA_Runtime_v1.0.1_maci64.zip)
        % Sanitize version string for filename (remove potential illegal chars)
        safe_ver = regexprep(version_str, '[^a-zA-Z0-9_\-\.]', '_');
        zipName = sprintf('%s_%s_%s.zip', appName, safe_ver, computer('arch'));
        
        % Change to output directory to create a clean zip structure (relative paths)
        owd = pwd;
        cd(outputDir);
        cleanupObj = onCleanup(@() cd(owd)); % Ensure we return to original dir
        
        filesToZip = {};
        
        % Add Installer
        if ismac
             filesToZip{end+1} = [installerName '.app'];
             filesToZip{end+1} = 'HERA_Launcher.command';
        elseif ispc
             filesToZip{end+1} = [installerName '.exe'];
        else
             filesToZip{end+1} = [installerName '.install'];
        end
        
        % Create Zip
        zip(zipName, filesToZip);
        
        % Success Message
        fprintf('\n========================================\n');
        fprintf('SUCCESS!\n');
        fprintf('1. Standalone App: %s\n', fullfile(outputDir, appName));
        fprintf('2. Installer:      %s\n', fullfile(outputDir, installerName));
        fprintf('3. ZIP Archive:    %s\n', fullfile(outputDir, zipName));
        if ismac
             fprintf('Includes Launcher: %s\n', fullfile(outputDir, 'HERA_Launcher.command'));
        end
        fprintf('----------------------------------------\n');
        fprintf('Ready for Release! Upload the ZIP file to GitHub.\n');
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