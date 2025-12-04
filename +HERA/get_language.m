function lang = get_language()
% GET_LANGUAGE - Loads the default English language structure from an external JSON file.
%
% Syntax:
%   lang = get_language()
%
% Description:
%   This function loads the default language file ('en.json') from the 'language' subfolder. 
%   It is designed to work seamlessly in both the MATLAB development environment and in a compiled runtime application. 
%   By centralizing text strings in a JSON file, this approach ensures the application is easy to maintain. 
%   It can be extended to support other languages in the future.
%
% Workflow:
%   1.  Set the default language code to 'en'.
%   2.  Determine the application's base path to locate the 'language' folder.
%   3.  Construct the full path to the 'en.json' file.
%   4.  Validate the file's existence.
%   5.  Read and parse the JSON file into a MATLAB struct.
%
% Outputs:
%   lang - A struct containing all localized text strings for the English language.
%
% Author: Lukas von Erdmannsdorff

%% 1. Initialization
    % Sets the language code directly to English.
    language_code = 'en';

%% 2. Determine Application Base Path
    % This block identifies the root directory of the application, which is
    % crucial for locating the 'language' folder reliably.
    if isdeployed
        % WHEN COMPILED: Get the root directory of the deployed application.
        base_path = mcr.runtime.getApplicationRoot;
    else
        % WHEN RUNNING IN MATLAB IDE: Get the directory of this .m file.
        base_path = fileparts(mfilename('fullpath'));
    end

%% 3. Construct Path to Language File
    % Builds the full path to the folder and the specific 'en.json' file.
    language_folder = fullfile(base_path, 'language');
    file_path = fullfile(language_folder, [language_code, '.json']);

%% 4. Validate File Existence and Load JSON
    % Checks if the default language file exists. If not, it throws a fatal error.
    if ~exist(file_path, 'file')
        error('Default language file could not be found: "%s". Application cannot continue.', file_path);
    end
    
    % Reads the JSON file content and parses it into a MATLAB struct.
    % A try-catch block ensures robust error handling if the file is corrupted or improperly formatted.
    try
        json_text = fileread(file_path);
        lang = jsondecode(json_text);
    catch ME
        error('Failed to load or parse language file "%s": %s', file_path, ME.message);
    end
end