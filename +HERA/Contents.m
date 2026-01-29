% HERA Toolbox
% Version 1.1.2 26-Jan-2026
%
% HERA: Hierarchical-Compensatory, Effect-Size-Driven and Non-Parametric Ranking Algorithm
%
% Files
%   bootstrap_ranking    - Performs bootstrap analysis for ranking stability.
%   borda_ranking        - Calculates Borda count rankings.
%   calculate_bca_ci     - Calculates Bias-Corrected and Accelerated (BCa) confidence intervals.
%   calculate_ranking    - Core function to calculate rankings based on input data.
%   calculate_thresholds - Computes thresholds for ranking criteria.
%   default              - Returns default settings and parameters for HERA.
%   design               - Manages design aspects/colors for HERA visualizations.
%   generate_output      - Handles output generation (exports, reports).
%   generate_plots       - Generates visualizations of ranking results.
%   get_language         - Retrieves language settings/localized strings.
%   get_version          - Returns the current version of the toolbox.
%   load_data            - Loads and preprocesses data for analysis.
%   power_analysis       - Conducts power analysis for the ranking results.
%   run_ranking          - Wrapper to execute the full ranking pipeline.
%   run_unit_test        - Runs unit tests for the toolbox.
%   start_ranking        - Main entry point for the HERA ranking CLI.
%
% Packages
%   analysis             - Convergence analysis (robust mode default values).
%   language             - Localization files (JSON).
%   output               - Output generation modules.
%   plot                 - Plotting functions.
%   run                  - Execution control modules.
%   start                - Initialization and startup routines.
%   stats                - Statistical utility functions.
%   test                 - Unit and integration tests.
%
%   license.txt          - License information.
%   README.md            - Toolbox documentation and usage instructions.
%
% See also HERA.start_ranking
