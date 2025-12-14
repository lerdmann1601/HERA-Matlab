function generate_plots(results, thresholds, shared_info, styles)
% GENERATE_PLOTS - Creates and saves all ranking graphics.
%
% Syntax:
%   generate_plots(results, thresholds, shared_info, styles)
%
% Description:
%   This function visualizes the results of the ranking algorithm by creating a series of graphics.
%   All figures are automatically saved as high-resolution PNG files.
%   Additionally, the graphics are grouped thematically and compiled into several PDF reports for easy sharing.
%
% Workflow:
%   1.  Creation of intermediate ranking plots: Visualizes the rankings after sorting by Metric 1.
%       Plots for Metric 2 and 3 are created conditionally based on the 'ranking_mode'.
%   2.  Creation of detailed comparison plots: Displays all pairwise comparisons (p-value, delta, rel. diff.) for each loaded metric.
%   3.  Creation of a win-loss matrix to visualize direct duels, dynamically showing 1, 2, or 3 matrices.
%   4.  Creation of a Sankey diagram to visualize rank shifts across the executed stages (dynamically 1, 2, or 3 stages).
%   5.  Creation of the final summary plot: Combines the final ranking (with bootstrap CI) with a table of metric values (dynamically 1-3 columns).
%   6.  Creation of the Borda plot (optional): Visualizes the result of the sensitivity analysis as a consensus ranking.
%   7.  Compilation into PDF reports: Saves the created graphics in thematic PDF files (e.g., Final_Report.pdf, Ranking_Report.pdf).
%
% Inputs:
%   results     - (struct) Contains all calculated results (ranks, CIs, p-values, etc.).
%   thresholds  - (struct) Contains the calculated statistical thresholds.
%   shared_info - (struct) Contains general data (names, paths) and configurations (incl. 'config' struct with 'ranking_mode').
%   styles      - (struct) Contains all color and style settings from the 'design' function.
%
% Outputs:
%   The function has no direct return values. 
%   Instead, it generates files:
%   Ranking plots           - PNG files of intermediate and final rankings.
%   Detail comparison plots - PNG files of pairwise comparisons for each metric.
%   Win-Loss Matrix         - PNG file of the duel overview.
%   Sankey Diagram          - PNG file of the rank shifts.
%   Borda Plot              - PNG file of the consensus ranking (if calculated).
%   PDF reports             - Several PDF files that bundle the graphics thematically.
%
% Author: Lukas von Erdmannsdorff

arguments
    results (1,1) struct
    thresholds (1,1) struct
    shared_info (1,1) struct
    styles (1,1) struct
end

%% 1. Initialization and Style Definition
% Unpack the passed structures into local variables for easier access.
metric_names = shared_info.metric_names;
dataset_names = shared_info.dataset_names;
num_datasets = shared_info.num_datasets;
pair_idx_all = shared_info.pair_idx_all;
d_thresh = thresholds.d_thresh;
rel_thresh = thresholds.rel_thresh;
output_dir = shared_info.output_dir;
graphics_dir = shared_info.graphics_dir;
pdf_dir      = shared_info.pdf_dir;
plot_theme = shared_info.plot_theme;
lang = shared_info.lang; % Load the full language pack.
ts = shared_info.config.timestamp; 

% Get dynamic metric count and ranking logic
num_metrics = numel(metric_names);
ranking_mode = shared_info.config.ranking_mode;

% Create subfolders for detailed comparison and ranking plots.
subfolder_detail_comp = fullfile(graphics_dir, 'Detail_Comparison');
if ~exist(subfolder_detail_comp, 'dir')
    mkdir(subfolder_detail_comp);
end
subfolder_ranking = fullfile(graphics_dir, 'Ranking');
if ~exist(subfolder_ranking, 'dir')
    mkdir(subfolder_ranking);
end

% Initialize arrays to collect all figure handles for PDF export.
final_report_handles = gobjects(0);
ranking_report_handles = gobjects(0);
convergence_report_handles = gobjects(0);

% Log the start of the plot creation process.
fprintf(['\n' lang.plots.log_messages.start_creation '\n'], upper(plot_theme));

% Set global background and font for consistency across all plots.
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');
set(0, 'DefaultAxesColor', styles.colors.background);

%% 2. Plot: Ranking by Metric 1
% Generate the initial ranking plot based on the primary metric.
% This is always created.
h_fig = HERA.plot.ranking_bar(1, ...
    results.intermediate_orders.after_metric1, ...
    lang.plots.figure_names.initial_ranking, ...
    lang.plots.ylabels.rank_metric_1, ...
    lang.plots.titles.initial_ranking, ...
    lang.files.initial_ranking_metric1, ...
    shared_info, styles, lang, subfolder_detail_comp);
ranking_report_handles(end+1) = h_fig;

%% 3. Plot: Detailed Pairwise Comparisons for Metric 1
% Generate the detailed comparison plots for the primary metric.
% This is always created.
handles = HERA.plot.detail_comparison(1, results, thresholds, shared_info, styles, lang, subfolder_detail_comp);
ranking_report_handles = [ranking_report_handles, handles];

%% 4. Plot: Ranking by Metric 2
% Conditional plotting based on logic
if num_metrics >= 2 && (strcmp(ranking_mode, 'M1_M2') || strcmp(ranking_mode, 'M1_M2_M3'))
    % Generate the ranking plot after corrections based on the second metric.
    h_fig = HERA.plot.ranking_bar(2, ...
        results.intermediate_orders.after_metric2, ...
        lang.plots.figure_names.rank_correction_m2, ...
        lang.plots.ylabels.rank_metric_2, ...
        sprintf(lang.plots.titles.ranking_after_correction_m2, metric_names{2}), ... 
        lang.files.correction_metric2, ...
        shared_info, styles, lang, subfolder_detail_comp);
    ranking_report_handles(end+1) = h_fig;
end

%% 5. Plot: Detailed Pairwise Comparisons for Metric 2
% Conditional plotting (if metric 2 exists, but NOT for M1_M3A, which is handled in Sec 6)
if num_metrics >= 2 && ~strcmp(ranking_mode, 'M1_M3A')
    % Generate the detailed comparison plots for the second metric.
    handles = HERA.plot.detail_comparison(2, results, thresholds, shared_info, styles, lang, subfolder_detail_comp);
    ranking_report_handles = [ranking_report_handles, handles];
end

%% 6. Plot: Ranking by Metric 3
% Conditional plotting based on logic
if num_metrics == 2 && strcmp(ranking_mode, 'M1_M3A')
    % Plot the final ranking after M3A logic (which used M2 data)
    h_fig = HERA.plot.ranking_bar(2, ... % Uses M2 data
        results.intermediate_orders.after_metric3, ... % Uses final order
        lang.plots.figure_names.rank_correction_m3, ...
        lang.plots.ylabels.rank_metric_2, ...
        sprintf(lang.plots.titles.ranking_after_correction_m3a, metric_names{2}), ... 
        lang.files.correction_metric3, ... % Saves as M3 correction
        shared_info, styles, lang, subfolder_detail_comp);
    ranking_report_handles(end+1) = h_fig;
    
    handles = HERA.plot.detail_comparison(2, results, thresholds, shared_info, styles, lang, subfolder_detail_comp);
    ranking_report_handles = [ranking_report_handles, handles];
        
elseif num_metrics == 3 && strcmp(ranking_mode, 'M1_M2_M3')
    % Plot the final ranking after M3 (A+B) logic
    % Generate the ranking plot after corrections based on the third metric.
    h_fig = HERA.plot.ranking_bar(3, ...
        results.intermediate_orders.after_metric3, ...
        lang.plots.figure_names.rank_correction_m3, ...
        lang.plots.ylabels.rank_metric_3, ...
        sprintf(lang.plots.titles.ranking_after_correction_m3, metric_names{3}), ... 
        lang.files.correction_metric3, ...
        shared_info, styles, lang, subfolder_detail_comp);
    ranking_report_handles(end+1) = h_fig;
end

%% 7. Plot: Detailed Pairwise Comparisons for Metric 3
% Conditional plotting (if metric 3 exists, regardless of logic)
if num_metrics == 3
    % Generate the detailed comparison plots for the third metric.
    handles = HERA.plot.detail_comparison(3, results, thresholds, shared_info, styles, lang, subfolder_detail_comp);
    ranking_report_handles = [ranking_report_handles, handles];
end

%% 8. Plot: Win-Loss-Matrix
h_fig_winloss = HERA.plot.win_loss_matrix(results, shared_info, styles, lang, subfolder_ranking);

%% 9. Plot: Sankey-Diagram of Rank Shifts
h_fig_sankey = HERA.plot.sankey_diagram(results, shared_info, styles, lang, subfolder_ranking);

%% 10. Plot: Final Summary (Ranking + Metrics Table)
h_fig_summary = HERA.plot.final_summary(results, shared_info, styles, lang, output_dir);

%% 11. Borda-Plot (if sensitivity analysis was performed)
h_fig_borda = HERA.plot.borda_consensus(results, shared_info, styles, lang, subfolder_ranking);

%% 12. Compile all graphics into PDF files
% Collect handles for the Convergence Report PDF.
h_fig_thr_global = shared_info.h_fig_thr_global;
h_fig_thr_detailed = shared_info.h_fig_thr_detailed;
h_fig_bca_global = shared_info.h_fig_bca_global;
h_fig_bca_detailed = shared_info.h_fig_bca_detailed;
h_fig_rank_conv = shared_info.h_fig_rank_conv;
convergence_report_handles = [h_fig_thr_global, h_fig_thr_detailed, h_fig_bca_global, h_fig_bca_detailed, h_fig_rank_conv];

% Collect handles for the Bootstrap Report PDF.
h_fig_hist_thr = shared_info.h_fig_hist_thr;
h_fig_hist_raw = shared_info.h_fig_hist_raw;
h_fig_hist_z0 = shared_info.h_fig_hist_z0;
h_fig_hist_a = shared_info.h_fig_hist_a;
h_fig_hist_widths = shared_info.h_fig_hist_widths;
h_fig_hist_rank = shared_info.h_fig_hist_rank;
bootstrap_report_handles = [h_fig_hist_thr, h_fig_hist_raw, h_fig_hist_z0, h_fig_hist_a, h_fig_hist_widths];

% Collect handles for the Final Report PDF.
final_report = gobjects(0);
final_report(end+1) = h_fig_summary;
if isgraphics(h_fig_borda) % Check if borda plot was created
    final_report(end+1) = h_fig_borda;
end
final_report(end+1) = h_fig_winloss;
if isgraphics(h_fig_sankey) % Check if sankey plot was created
    final_report(end+1) = h_fig_sankey;
end
final_report(end+1) = h_fig_hist_rank;

% Collect all figure handles to close them after saving.
all_handles_to_close = unique([ranking_report_handles, final_report, convergence_report_handles, bootstrap_report_handles]);

% Define a local helper function to save an array of figures to a single PDF.
save_pdf_report = @(handles, filename) ...
    arrayfun(@(i) exportgraphics(handles(i), filename, 'ContentType', 'vector', 'Append', i > 1, 'BackgroundColor', styles.colors.background), 1:numel(handles));

% Generate Convergence_Report.pdf.
[~, fName, fExt] = fileparts(lang.files.convergence_report);
pdf_filename_conv = fullfile(pdf_dir, [fName, '_', ts, fExt]);
if exist(pdf_filename_conv, 'file'), delete(pdf_filename_conv); end % Delete old file if it exists.
if ~isempty(convergence_report_handles)
    save_pdf_report(convergence_report_handles, pdf_filename_conv);
    fprintf(['  ' lang.plots.log_messages.pdf_saved '\n'], pdf_filename_conv);
end

% Generate Bootstrap_Report.pdf.
[~, fName, fExt] = fileparts(lang.files.bootstrap_report);
pdf_filename_boot = fullfile(pdf_dir, [fName, '_', ts, fExt]);
if exist(pdf_filename_boot, 'file'), delete(pdf_filename_boot); end
if ~isempty(bootstrap_report_handles)
    save_pdf_report(bootstrap_report_handles, pdf_filename_boot);
    fprintf(['  ' lang.plots.log_messages.pdf_saved '\n'], pdf_filename_boot);
end

% Generate Ranking_Report.pdf.
[~, fName, fExt] = fileparts(lang.files.ranking_report);
pdf_filename_rank = fullfile(pdf_dir, [fName, '_', ts, fExt]);
if exist(pdf_filename_rank, 'file'), delete(pdf_filename_rank); end
if ~isempty(ranking_report_handles)
    save_pdf_report(ranking_report_handles, pdf_filename_rank);
    fprintf(['  ' lang.plots.log_messages.pdf_saved '\n'], pdf_filename_rank);
end

% Generate Final_Report.pdf.
[~, fName, fExt] = fileparts(lang.files.final_report);
pdf_filename_final = fullfile(output_dir, [fName, '_', ts, fExt]);
if exist(pdf_filename_final, 'file'), delete(pdf_filename_final); end
if ~isempty(final_report)
    save_pdf_report(final_report, pdf_filename_final);
    fprintf(['  ' lang.plots.log_messages.pdf_saved '\n'], pdf_filename_final);
end

% Close all generated figures to free up memory.
close(all_handles_to_close);

fprintf([lang.plots.log_messages.creation_finished '\n']);
end