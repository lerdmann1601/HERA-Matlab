function styles = design(plot_theme, num_datasets, has_power_results)
% DESIGN - Defines and returns a structure with all style settings for graphics.
%
% Syntax:
%   styles = design(plot_theme, num_datasets, has_power_results)
%
% Description:
%   This function acts as a central style repository. 
%   It generates a struct containing all color definitions, font sizes, and line/marker styles required for generating plots. 
%   It supports both a 'light' and a 'dark' theme, which is selected via the input argument.
%   Future designs like color blind will be added here.
%
% Workflow:
%   1. Theme Selection: 
%      Sets base colors (background, text, grid) based on whether 'plot_theme' is 'dark' or 'light'.
%   2. General Design: 
%      Defines global styles for lines, markers, fonts and error bars used across multiple plots.
%   3. Specific Settings: 
%      Defines unique color palettes and styles for plots like Win-Loss-Matrix (adjusting for power analysis) and Sankey diagram (creating an 'hsv' palette).
%
% Inputs:
%   plot_theme        - String defining the theme ('light' or 'dark').
%   num_datasets      - Scalar, the number of datasets (required for color palettes).
%   has_power_results - Boolean, true if power analysis results are available.
%
% Outputs:
%   styles            - Structure containing all color, line, marker, and font definitions.
%
% Author: Lukas von Erdmannsdorff

styles = struct();
is_dark_theme = contains(plot_theme, 'dark');
is_colourblind = contains(plot_theme, 'colourblind');

%% Dark Theme
if is_dark_theme
    % Dark Mode
    styles.colors.background = [0 0 0];
    styles.colors.text = [1 1 1];
    styles.colors.convergence = [0.8500 0.3250 0.0980];
    styles.colors.sem_threshold = [0.98, 0.75, 0.2];
    %  Colors 
    styles.colors.p_significant = [0.0, 0.45, 0.74];     % Rich Blue
    styles.colors.p_nonsignificant = [0.85, 0.85, 0.85]; % Light Gray
    styles.colors.holm_threshold = [0.8, 0, 0];          % Strong Red
    styles.colors.intermediate_ci = [0.6, 0.6, 0.6];     % Gray 
    styles.colors.win = [0.1, 0.65, 0.1];
    styles.colors.loss = [0.85, 0.1, 0.1];
    styles.colors.grid_color = [0.45, 0.45, 0.45];
    styles.colors.marker_face = [1 1 1];
    styles.colors.marker_edge = [1 1 1];
    styles.colors.delta_face = [0.3, 0.6, 1.0];
    styles.colors.rel_face = [1.0, 0.5, 0.2];
    styles.colors.bar_face = [0.3, 0.6, 1.0];
    styles.colors.bar_edge = [0.8, 0.8, 0.8];
    styles.colors.kde_line = [0.9, 0.9, 0.9];
    styles.colors.red_marker = [1 0.1 0.1];
    styles.colors.blue_marker = [0.3 0.6 1.0];
    styles.colors.borda = [0.8, 0.8, 0.8];

%% Light Theme
else
    % Light Mode
    styles.colors.background = [1 1 1];
    styles.colors.text = [0 0 0];
    styles.colors.convergence = [0.8500 0.3250 0.0980];
    styles.colors.sem_threshold = [0.9290 0.6940 0.1250];
    styles.colors.p_significant = [0.0, 0.45, 0.74];     % Rich Blue
    styles.colors.p_nonsignificant = [0.85, 0.85, 0.85]; % Light Gray
    styles.colors.holm_threshold = [0.8, 0, 0];          % Strong Red
    styles.colors.intermediate_ci = [0.6, 0.6, 0.6];     % Gray 
    styles.colors.win = [0.1, 0.65, 0.1];
    styles.colors.loss = [0.85, 0.1, 0.1];
    styles.colors.grid_color = [0.5, 0.5, 0.5];
    styles.colors.marker_face = [0 0 0];
    styles.colors.marker_edge = [0 0 0];
    styles.colors.delta_face = [0 0.4470 0.7410];
    styles.colors.rel_face = [0.8500 0.3250 0.0980];
    styles.colors.bar_face = [0 0.4470 0.7410];
    styles.colors.bar_edge = [0.15 0.15 0.15];
    styles.colors.kde_line = [0.2 0.2 0.2];
    styles.colors.red_marker = [1 0 0];
    styles.colors.blue_marker = [0 0.4470 0.7410];
    styles.colors.borda = [0.8, 0.8, 0.8];
end

%% General Design
styles.colors.neutral = styles.colors.p_nonsignificant;

% Lines and Markers
styles.line.grid = ':';
styles.line.threshold = {styles.line.grid, 'Color', styles.colors.holm_threshold, 'LineWidth', 1.5};
styles.marker.final_rank = {'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 5};
styles.marker.holm = {'o', 'MarkerFaceColor', styles.colors.holm_threshold, 'MarkerEdgeColor', 'none', 'MarkerSize', 5};

% Error Bars 
styles.errorbar.intermediate = {'horizontal', 'LineStyle', 'none', 'Color', styles.colors.intermediate_ci, 'CapSize', 5, 'LineWidth', 1.5};
styles.errorbar.final_rank = {'horizontal', 'LineStyle', '-', 'Color', styles.colors.intermediate_ci, 'LineWidth', 1, 'CapSize', 6};
styles.errorbar.detail_ci = {'horizontal', 'LineStyle', 'none', 'Color', styles.colors.intermediate_ci, 'LineWidth', 1, 'CapSize', 6}; 

% Font sizes remain global
styles.font.title = 12; styles.font.label = 11;
styles.font.tick = 10; styles.font.text = 10; styles.font.small_text = 9;

% Specific settings for Win-Loss-Matrix 
styles.win_loss.intensities = [0.3, 0.45, 0.7, 1.0]; % For good contrast
styles.win_loss.gray_values = [0.15, 0.3, 0.6, 0.9];
styles.win_loss.diag_color  = [0.95, 0.95, 0.95];
styles.win_loss.legend_intensities = styles.win_loss.intensities;
styles.win_loss.legend_gray_values = styles.win_loss.gray_values;

% Override colors if no power analysis was performed
if ~has_power_results
    styles.win_loss.intensities = [1, 1, 1, 1]; % Full colors
    styles.win_loss.gray_values = [0.0, 0.0, 0.0, 0.0]; % Neutral becomes black
end

% Specific settings for Sankey diagram
styles.sankey.colors = hsv(num_datasets);

%% Colourblind Theme Override
styles.is_colourblind = is_colourblind;

if is_colourblind
    if is_dark_theme
        % Override using Okabe-Ito colorblind-friendly palette (Optimized for Dark Background)
        styles.colors.p_significant = [0.34, 0.71, 0.91];     % Sky Blue
        styles.colors.holm_threshold = [0.90, 0.62, 0.0];     % Orange
        styles.colors.win = [0.34, 0.71, 0.91];               % Sky Blue
        styles.colors.loss = [0.90, 0.62, 0.0];               % Orange
        styles.colors.delta_face = [0.34, 0.71, 0.91];        
        styles.colors.bar_face = [0.34, 0.71, 0.91];          
        styles.colors.red_marker = [0.90, 0.62, 0.0];         
        styles.colors.blue_marker = [0.34, 0.71, 0.91];       
        styles.colors.convergence = [0.95, 0.90, 0.25];       % Yellow
        styles.colors.rel_face = [0.95, 0.90, 0.25];          
        
        styles.win_loss.intensities = [0.5, 0.7, 0.85, 1.0];
        styles.win_loss.gray_values = [0.3, 0.5, 0.7, 0.9];
    else
        % Override using Okabe-Ito colorblind-friendly palette (Optimized for Light Background)
        styles.colors.p_significant = [0.0, 0.45, 0.70];      % Blue
        styles.colors.holm_threshold = [0.84, 0.37, 0.0];     % Vermilion
        styles.colors.win = [0.0, 0.45, 0.70];                % Blue
        styles.colors.loss = [0.84, 0.37, 0.0];               % Vermilion
        styles.colors.delta_face = [0.0, 0.45, 0.70];         
        styles.colors.bar_face = [0.0, 0.45, 0.70];           
        styles.colors.red_marker = [0.84, 0.37, 0.0];         
        styles.colors.blue_marker = [0.0, 0.45, 0.70];        
        styles.colors.convergence = [0.90, 0.62, 0.0];        % Orange
        styles.colors.rel_face = [0.90, 0.62, 0.0];           
        
        styles.win_loss.intensities = [0.4, 0.6, 0.8, 1.0];
        styles.win_loss.gray_values = [0.2, 0.4, 0.6, 0.8];
    end
    
    % Adjust contrast for win-loss matrix to maintain optimal readability
    if ~has_power_results
        styles.win_loss.intensities = [1, 1, 1, 1];
        styles.win_loss.gray_values = [0, 0, 0, 0];
    end
    styles.win_loss.legend_intensities = styles.win_loss.intensities;
    styles.win_loss.legend_gray_values = styles.win_loss.gray_values;
end

end