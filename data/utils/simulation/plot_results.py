"""
Plotting Suite for HERA Methodological Ground-Truth Validation.

Reads Monte Carlo benchmark results and generates publication-grade,
color-blind compliant visualizations comparing HERA against standard
Multi-Criteria Decision Making (MCDM) baselines (TOPSIS, Borda Count, Copeland).

Visualizations generated in <data_folder>/Simulation_Output_<timestamp>/Graphics/:
1. validation_overview_<timestamp>.png      : Master 8-panel publication overview.
2. top_choice_recovery_vs_noise_<ts>.png   : Top-choice recovery rates faceted across candidate scale (N=10, 12, 14).
3. complete_rank_recovery_vs_noise_<ts>.png : Complete-rank recovery rates faceted across candidate scale.
4. rank_displacement_vs_noise_<ts>.png     : Mean rank displacement across noise levels.
5. false_superiority_vs_noise_<ts>.png     : Pairwise false superiority rates across noise levels and candidate scales.
6. rank_correlation_vs_sample_<ts>.png     : Kendall's Tau and Spearman's Rho rank correlation vs sample size.
7. effect_magnitude_sensitivity_<ts>.png   : Assigned Rank vs Ground Truth across Small, Medium, and Large effects.
8. correlation_sensitivity_<ts>.png        : Robustness under inter-metric collinearity (rho = 0.0 vs rho = 0.5).
9. compensatory_failure_vs_noise_<ts>.png  : Non-compensatory diagnostic stress test failure rates.
10. candidate_rank_stability_<ts>.png      : Candidate-level middle-field rank stability across scales.
11. pooled_core_distributions_<ts>.png     : Pooled non-parametric distributions across core conditions.

All plots are rendered at 300 DPI using the Wong/Okabe-Ito colorblind palette 
with distinct marker shapes and linestyles to ensure accessibility in print and digital media.

Usage:
    python3 plot_results.py

Author: Lukas von Erdmannsdorff
"""

import os
import shutil
from pathlib import Path
from typing import Tuple, Dict, List, Optional, Any
from datetime import datetime

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter, MultipleLocator
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.image as mpimg
import numpy as np
import pandas as pd
import seaborn as sns

# Import simulation configuration parameters
try:
    from .config import (
        DEFAULT_CANDIDATES, DEFAULT_SAMPLE_SIZE, DEFAULT_NOISE,
        DEFAULT_EFFECT, DEFAULT_CORRELATION, NUM_CANDIDATES_LIST,
        SUITE_CORE, SUITE_SAMPLE_SIZE, SUITE_EFFECT, SUITE_CORRELATION
    )
    from .reporting import compute_superordinate_summary
except ImportError:
    try:
        from simulation.config import (
            DEFAULT_CANDIDATES, DEFAULT_SAMPLE_SIZE, DEFAULT_NOISE,
            DEFAULT_EFFECT, DEFAULT_CORRELATION, NUM_CANDIDATES_LIST,
            SUITE_CORE, SUITE_SAMPLE_SIZE, SUITE_EFFECT, SUITE_CORRELATION
        )
        from simulation.reporting import compute_superordinate_summary
    except ImportError:
        DEFAULT_CANDIDATES = 12
        DEFAULT_SAMPLE_SIZE = 50
        DEFAULT_NOISE = 4.0
        DEFAULT_EFFECT = "Medium"
        DEFAULT_CORRELATION = 0.0
        NUM_CANDIDATES_LIST = [10, 12, 14]
        SUITE_CORE = "Core"
        compute_superordinate_summary = None

# --- Directory Configuration ---
BASE_DIR = Path(__file__).parent.resolve()
if BASE_DIR.name == "simulation":
    UTILS_DIR = BASE_DIR.parent
    DATA_DIR = UTILS_DIR.parent
else:
    UTILS_DIR = BASE_DIR
    DATA_DIR = BASE_DIR.parent
RESULTS_FILE = DATA_DIR / "simulation_results.csv"
SUMMARY_FILE = DATA_DIR / "global_summary.csv"
OUTPUT_DIR = DATA_DIR / "plots"

# --- Visual Styling (Publication Standards) ---
sns.set_theme(style="whitegrid", context="paper", font_scale=1.2)
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["axes.edgecolor"] = "#333333"
plt.rcParams["axes.linewidth"] = 0.8
plt.rcParams["grid.color"] = "#E0E0E0"
plt.rcParams["grid.linestyle"] = "--"
plt.rcParams["grid.alpha"] = 0.7

# Colorblind-safe palette (Okabe-Ito / Seaborn Colorblind)
COLORS = {
    "HERA": "#0173b2",              # Accessible Blue
    "TOPSIS": "#de8f05",            # High-contrast Amber
    "Borda": "#029e73",             # Jade Green
    "Wilcoxon-Copeland": "#cc78bc", # Soft Purple
    "Copeland": "#cc78bc"           # Backwards compatibility alias
}

# Concentric marker and hierarchical linestyle specifications
# Guarantees that when methods achieve identical values (e.g. 1.0 or 0.0),
# all 4 methods remain distinctly visible as a nested 'bullseye' target
# and through interleaved dashed/dotted line styles.
METHOD_SPECS = {
    "HERA": {
        "color": "#0173b2",       # Blue
        "linestyle": "-",          # Solid
        "linewidth": 2.8,
        "marker": "o",             # Large filled circle
        "markersize": 9.5,
        "markerfacecolor": "#0173b2",
        "markeredgecolor": "white",
        "markeredgewidth": 1.2,
        "alpha": 0.85,
        "zorder": 10
    },
    "TOPSIS": {
        "color": "#de8f05",       # Amber
        "linestyle": "--",         # Dashed
        "linewidth": 2.2,
        "dashes": (4, 2),
        "marker": "s",             # Open square
        "markersize": 8.0,
        "markerfacecolor": "none",
        "markeredgecolor": "#de8f05",
        "markeredgewidth": 2.0,
        "alpha": 0.90,
        "zorder": 9
    },
    "Borda": {
        "color": "#029e73",       # Green
        "linestyle": ":",          # Dotted
        "linewidth": 2.2,
        "dashes": (1.5, 2),
        "marker": "^",             # Open triangle
        "markersize": 6.5,
        "markerfacecolor": "none",
        "markeredgecolor": "#029e73",
        "markeredgewidth": 2.0,
        "alpha": 0.95,
        "zorder": 8
    },
    "Wilcoxon-Copeland": {
        "color": "#cc78bc",       # Purple
        "linestyle": "-.",         # Dash-dot
        "linewidth": 1.8,
        "dashes": (5, 2, 1.5, 2),
        "marker": "D",             # Small filled diamond
        "markersize": 4.5,
        "markerfacecolor": "#cc78bc",
        "markeredgecolor": "white",
        "markeredgewidth": 0.8,
        "alpha": 0.95,
        "zorder": 7
    },
    "Copeland": {
        "color": "#cc78bc",       # Purple
        "linestyle": "-.",         # Dash-dot
        "linewidth": 1.8,
        "dashes": (5, 2, 1.5, 2),
        "marker": "D",             # Small filled diamond
        "markersize": 4.5,
        "markerfacecolor": "#cc78bc",
        "markeredgecolor": "white",
        "markeredgewidth": 0.8,
        "alpha": 0.95,
        "zorder": 7
    }
}

__all__ = [
    "COLORS",
    "METHOD_SPECS",
    "get_runs_per_point",
    "plot_method_lines",
    "plot_master_overview",
    "plot_faceted_top_choice",
    "plot_faceted_complete_rank",
    "plot_faceted_rank_displacement",
    "plot_faceted_regret",
    "plot_faceted_false_superiority",
    "plot_faceted_rank_correlation",
    "plot_effect_sensitivity",
    "plot_correlation_sensitivity",
    "plot_compensatory_failure",
    "plot_candidate_rank_stability",
    "plot_pooled_core_distributions",
    "plot_pooled_core_marginal_scales",
    "plot_pooled_compensatory_stress",
    "plot_pooled_sensitivity_summary",
    "plot_superordinate_table_figure",
    "generate_global_summary_pdf",
    "generate_all_plots",
    "main"
]


def get_runs_per_point(df: Optional[pd.DataFrame], default: int = 10) -> int:
    """Detects number of Monte Carlo iterations/runs per condition in dataset.

    Syntax:
        K = get_runs_per_point(df, default=10)

    Description:
        Inspects the 'Iteration' column of the provided results DataFrame to compute
        the exact empirical number of Monte Carlo replications per experimental point.
        Falls back to the provided default if the DataFrame is empty or unindexed.

    Parameters:
        df (Optional[pd.DataFrame]): Input results DataFrame.
        default (int): Fallback number of runs if detection fails.

    Returns:
        int: Number of runs per experimental condition.

    Author:
        Lukas von Erdmannsdorff
    """
    if df is not None and not df.empty and "Iteration" in df.columns:
        valid_iters = df["Iteration"].dropna().unique()
        if len(valid_iters) > 0:
            return int(len(valid_iters))
    return default


def plot_method_lines(
    ax: plt.Axes, 
    df_data: pd.DataFrame, 
    x_col: str, 
    y_col: str, 
    show_ci: bool = True,
    estimator: Any = "mean",
    errorbar: Any = ("ci", 95)
) -> None:
    """Plots lines with concentric markers and hierarchical linestyles for all methods.

    Syntax:
        plot_method_lines(ax, df_data, x_col, y_col, show_ci=True, estimator='mean', errorbar=('ci', 95))

    Description:
        Renders method-specific curves across an experimental sweep with publication-grade
        styling. Enforces distinct visual markers, colors, line widths, and z-orders
        calibrated for colorblind accessibility and black-and-white print legibility.

    Parameters:
        ax (plt.Axes): Target matplotlib Axes instance.
        df_data (pd.DataFrame): Dataset containing experimental sweep data.
        x_col (str): Column name for the independent variable (X-axis).
        y_col (str): Column name for the dependent metric (Y-axis).
        show_ci (bool): Whether to calculate and display confidence intervals / error bands.
        estimator (Any): Aggregation function (e.g., 'mean', np.median).
        errorbar (Any): Error band specification for Seaborn (e.g., ('ci', 95), ('pi', 50)).

    Returns:
        None

    Author:
        Lukas von Erdmannsdorff
    """
    target_methods = ["HERA", "TOPSIS", "Borda", "Wilcoxon-Copeland", "Copeland"]
    methods = [m for m in target_methods if m in df_data["Method"].unique()]
    for m in methods:
        sub = df_data[df_data["Method"] == m]
        if sub.empty:
            continue
        spec = METHOD_SPECS[m]
        line_kws = {
            "linestyle": spec["linestyle"],
            "linewidth": spec["linewidth"],
            "marker": spec["marker"],
            "markersize": spec["markersize"],
            "markerfacecolor": spec["markerfacecolor"],
            "markeredgecolor": spec["markeredgecolor"],
            "markeredgewidth": spec["markeredgewidth"],
            "zorder": spec["zorder"],
            "alpha": spec["alpha"]
        }
        if "dashes" in spec and spec["dashes"]:
            line_kws["dashes"] = spec["dashes"]

        sns.lineplot(
            data=sub,
            x=x_col,
            y=y_col,
            ax=ax,
            color=spec["color"],
            label=m,
            estimator=estimator,
            errorbar=errorbar if show_ci else None,
            **line_kws
        )


def add_external_legend(
    fig: plt.Figure, 
    ax_source: plt.Axes, 
    ncol: int = 4, 
    y_pos: float = -0.05, 
    title: Optional[str] = None
) -> Optional[matplotlib.legend.Legend]:
    """Attaches a consolidated, single external legend centered below subplots.

    Syntax:
        leg = add_external_legend(fig, ax_source, ncol=4, y_pos=-0.05, title=None)

    Description:
        Extracts method handles and labels from a source subplot axis, de-duplicates
        them while preserving canonical order, and renders an elegant, unified
        legend beneath the figure grid.

    Parameters:
        fig (plt.Figure): Matplotlib Figure canvas.
        ax_source (plt.Axes): Subplot axes containing method line handles.
        ncol (int): Number of horizontal columns in the legend.
        y_pos (float): Vertical Y-coordinate in figure coordinates (typically negative).
        title (Optional[str]): Optional legend title string.

    Returns:
        Optional[matplotlib.legend.Legend]: The created legend instance.

    Author:
        Lukas von Erdmannsdorff
    """
    handles, labels = ax_source.get_legend_handles_labels()
    # Filter to unique labels preserving order
    by_label = dict(zip(labels, handles))
    leg = fig.legend(
        by_label.values(),
        by_label.keys(),
        loc="lower center",
        ncol=ncol,
        title=title,
        frameon=True,
        facecolor="white",
        edgecolor="#CCCCCC",
        framealpha=0.95,
        fontsize=10.5,
        title_fontsize=11,
        bbox_to_anchor=(0.5, y_pos)
    )
    return leg


def save_figure(fig: plt.Figure, png_path: Path, dpi: int = 300) -> None:
    """Saves a matplotlib figure as high-resolution publication PNG.

    Syntax:
        save_figure(fig, png_path, dpi=300)

    Description:
        Ensures target parent directory exists, applies tight bounding box clipping,
        and saves the figure canvas at 300 DPI.

    Parameters:
        fig (plt.Figure): Matplotlib Figure to export.
        png_path (Path): Destination filesystem path.
        dpi (int): Resolution in dots per inch (default: 300).

    Returns:
        None

    Author:
        Lukas von Erdmannsdorff
    """
    png_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png_path, dpi=dpi, bbox_inches="tight")


def save_figure_dual(
    fig: plt.Figure, 
    png_path: Path, 
    pdf_path: Optional[Path] = None, 
    dpi: int = 300
) -> None:
    """Saves figure as publication PNG (consolidated PDF reports bundle all figures).

    Syntax:
        save_figure_dual(fig, png_path, pdf_path=None, dpi=300)

    Description:
        Wrapper maintaining the dual-format save API. Saves 300 DPI PNG; loose PDFs
        are superseded by the consolidated multi-page PDF reports in Reports/.

    Parameters:
        fig (plt.Figure): Matplotlib Figure to export.
        png_path (Path): Destination filesystem path for PNG.
        pdf_path (Optional[Path]): Optional target path for loose PDF (ignored).
        dpi (int): Resolution in dots per inch (default: 300).

    Returns:
        None

    Author:
        Lukas von Erdmannsdorff
    """
    save_figure(fig, png_path, dpi=dpi)


def load_data(
    results_file: Optional[Path] = None, 
    summary_file: Optional[Path] = None
) -> Tuple[pd.DataFrame, Optional[pd.DataFrame]]:
    """Loads and validates the simulation raw results dataset and aggregate summary.

    Syntax:
        df_raw, df_summary = load_data(results_file=None, summary_file=None)

    Description:
        Locates and loads raw simulation CSV results. Automatically searches fallback
        paths and parent directories, augments TopChoiceRankRegret from candidate ranks
        if missing, and pairs the raw dataset with its corresponding global summary.

    Parameters:
        results_file (Optional[Path]): Explicit path to simulation_results_*.csv.
        summary_file (Optional[Path]): Explicit path to global_summary_*.csv.

    Returns:
        Tuple[pd.DataFrame, Optional[pd.DataFrame]]: Loaded raw results and summary DataFrames.

    Author:
        Lukas von Erdmannsdorff
    """
    res_path = results_file if results_file is not None else RESULTS_FILE
    sum_path = summary_file if summary_file is not None else SUMMARY_FILE

    if not res_path.exists():
        # Fallback: search for latest Simulation_Output_* directory in CSVs/ or root
        candidates = sorted(DATA_DIR.glob("Simulation_Output_*/CSVs/simulation_results*.csv"))
        if not candidates:
            candidates = sorted(DATA_DIR.glob("Simulation_Output_*/simulation_results*.csv"))
        if candidates:
            res_path = candidates[-1]
            sum_path = None
        else:
            raise FileNotFoundError(f"Results file not found: {res_path}\nPlease run run_simulation.py first.")

    # Search for matching summary CSV if not found directly
    if sum_path is None or not sum_path.exists():
        for s_dir in [res_path.parent, res_path.parent.parent / "CSVs", res_path.parent / "CSVs"]:
            if s_dir.exists():
                cands = sorted(s_dir.glob("global_summary*.csv"))
                if cands and os.path.getsize(cands[-1]) > 0:
                    sum_path = cands[-1]
                    break
    
    df_raw = pd.read_csv(res_path)
    if "TopChoiceRankRegret" not in df_raw.columns:
        cr_path = None
        search_dirs = [res_path.parent, res_path.parent.parent / "CSVs", res_path.parent / "CSVs", res_path.parent.parent]
        for s_dir in search_dirs:
            if s_dir.exists():
                cands = sorted(s_dir.glob("candidate_ranks*.csv"))
                if cands and os.path.getsize(cands[-1]) > 0:
                    cr_path = cands[-1]
                    break
        if cr_path and cr_path.exists():
            df_cr = pd.read_csv(cr_path)
            top_picks = df_cr[df_cr["PredictedRank"] == 1].copy()
            top_picks["TopChoiceRankRegret"] = (top_picks["TrueRank"] - 1).astype(float)
            df_raw = pd.merge(df_raw, top_picks[["ScenarioID", "Iteration", "Method", "TopChoiceRankRegret"]], on=["ScenarioID", "Iteration", "Method"], how="left")

    df_summary = pd.read_csv(sum_path) if (sum_path and sum_path.exists()) else None
    return df_raw, df_summary


def plot_master_overview(
    df: pd.DataFrame, 
    output_path: Path, 
    pdf_path: Optional[Path] = None,
    df_ranks: Optional[pd.DataFrame] = None
) -> None:
    r"""Generates the comprehensive 8-panel Master Overview figure directly addressing Reviewer Point 6.

    Syntax:
        plot_master_overview(df, output_path, pdf_path=None, df_ranks=None)

    Description:
        Creates an 8-panel multi-dimensional validation figure summarizing core
        performance metrics across noise levels and sensitivity sweeps.
        
        Layout (2 rows x 4 columns):
          Top Row (Core Decision Metrics vs. Noise Level sigma):
            (A) Top-Choice Recovery Rate vs. Noise (sigma)
            (B) Complete-Rank Recovery Rate vs. Noise (sigma)
            (C) Pairwise Inversion Rate (False Superiority Rate / Kendall Distance) vs. Noise (sigma)
            (D) Monotonic Rank Alignment (Spearman's \rho) vs. Noise (\sigma)
          Bottom Row (Multi-Dimensional Robustness Dimensions):
            (E) Rank Fidelity vs. Sample Size (n = 25, 50, 100)
            (F) Scalability across Candidate Scale (N = 10, 12, 14)
            (G) Rank Fidelity vs. Effect Size Magnitude (Cliff's d \approx 0.25, 0.50, 0.80)
            (H) Robustness under Inter-Metric Collinearity (\rho = 0.0 vs. 0.5)

    Parameters:
        df (pd.DataFrame): Raw Monte Carlo simulation results dataset.
        output_path (Path): Target file path for the 300 DPI PNG figure.
        pdf_path (Optional[Path]): Optional target path for vector PDF export.
        df_ranks (Optional[pd.DataFrame]): Optional candidate ranks dataframe for rank fidelity.

    Returns:
        None

    Author:
        Lukas von Erdmannsdorff
    """
    ref_N = DEFAULT_CANDIDATES if DEFAULT_CANDIDATES in df["Candidates"].values else df["Candidates"].iloc[0]
    target_methods = ["HERA", "TOPSIS", "Borda", "Wilcoxon-Copeland"]

    fig, axes = plt.subplots(2, 4, figsize=(18.5, 9.6))
    title_fontsize = 13.5
    subplot_title_fontsize = 11.5
    fig.suptitle("Benchmark Validation Overview", fontweight="bold", fontsize=title_fontsize, y=0.985)
    K_overview = get_runs_per_point(df)
    fig.text(
        0.5, 0.948,
        rf"Default parameters (unless otherwise specified): $\mathbf{{N = {ref_N}}}$, $\mathbf{{n = 50}}$, Cliff's $\mathbf{{d \approx 0.50}}$, $\boldsymbol{{\rho = 0.0}}$, $\mathbf{{K = {K_overview}}}$ runs per point",
        ha="center", va="top", fontsize=subplot_title_fontsize, fontweight="bold", color="#222222"
    )

    # Filter to reference condition for panels A, B, C, D (varying noise at standard n=50, d=0.50, rho=0.0)
    df_core = df[
        (df["Candidates"] == ref_N) & 
        (df["SampleSize"] == DEFAULT_SAMPLE_SIZE) & 
        (df["EffectMagnitude"] == DEFAULT_EFFECT) & 
        (df["Correlation"] == DEFAULT_CORRELATION)
    ]
    if df_core.empty:
        df_core = df[(df["Candidates"] == ref_N) & (df["EffectMagnitude"] == DEFAULT_EFFECT) & (df["Correlation"] == DEFAULT_CORRELATION)]

    # (A) Top-Choice Recovery vs Noise
    ax = axes[0, 0]
    plot_method_lines(ax, df_core, "Noise", "TopChoice")
    ax.set_title("(A) Top-Choice Recovery vs. Noise", fontweight="bold", loc="left")
    ax.set_xlabel(r"Noise Level ($\sigma$, %)")
    ax.set_ylabel("Top-Choice Recovery (Mean [95% CI], %)")
    ax.tick_params(labelleft=True)
    ax.set_ylim(-0.02, 1.05)
    ax.yaxis.set_major_locator(MultipleLocator(0.20))
    ax.yaxis.set_major_formatter(PercentFormatter(1.0, decimals=0))

    # (B) Total Rank Recovery vs Noise
    ax = axes[0, 1]
    plot_method_lines(ax, df_core, "Noise", "CompleteRank")
    ax.set_title("(B) Total Rank Recovery vs. Noise", fontweight="bold", loc="left")
    ax.set_xlabel(r"Noise Level ($\sigma$, %)")
    ax.set_ylabel("Total Rank Recovery (Mean [95% CI], %)")
    ax.tick_params(labelleft=True)
    ax.set_ylim(-0.02, 1.05)
    ax.yaxis.set_major_locator(MultipleLocator(0.20))
    ax.yaxis.set_major_formatter(PercentFormatter(1.0, decimals=0))

    # (C) False Superiority vs Noise
    ax = axes[0, 2]
    plot_method_lines(ax, df_core, "Noise", "FalseSuperiority", estimator=np.median, errorbar=("pi", 50))
    ax.set_title("(C) Pairwise Inversion Rate (FSR) vs. Noise", fontweight="bold", loc="left")
    ax.set_xlabel(r"Noise Level ($\sigma$, %)")
    ax.set_ylabel("False Superiority (Median [IQR], %)")
    ax.tick_params(labelleft=True)
    ax.set_ylim(-0.005, 0.32)
    ax.yaxis.set_major_locator(MultipleLocator(0.05))
    ax.yaxis.set_major_formatter(PercentFormatter(1.0, decimals=0))

    # (D) Monotonic Alignment (Spearman's rho) vs Noise
    ax = axes[0, 3]
    plot_method_lines(ax, df_core, "Noise", "SpearmanRho", estimator=np.median, errorbar=("pi", 50))
    ax.set_title(r"(D) Monotonic Alignment (Spearman's $\boldsymbol{\rho}$) vs. Noise", fontweight="bold", loc="left")
    ax.set_xlabel(r"Noise Level ($\sigma$, %)")
    ax.set_ylabel(r"Spearman's $\rho$ (Median [IQR])")
    ax.tick_params(labelleft=True)
    ax.set_ylim(0.65, 1.02)
    ax.yaxis.set_major_locator(MultipleLocator(0.10))

    # (E) Kendall Tau vs Sample Size (n = 25, 50, 100)
    ax = axes[1, 0]
    df_sample = df[
        (df["Candidates"] == ref_N) & 
        (df["Noise"] == DEFAULT_NOISE) & 
        (df["EffectMagnitude"] == DEFAULT_EFFECT) & 
        (df["Correlation"] == DEFAULT_CORRELATION)
    ]
    if df_sample.empty:
        df_sample = df_core
    plot_method_lines(ax, df_sample, "SampleSize", "KendallTau", estimator=np.median, errorbar=("pi", 50))
    ax.set_title("(E) Rank Fidelity vs. Sample Size", fontweight="bold", loc="left")
    ax.set_xlabel(r"Cohort Sample Size ($n$)")
    ax.set_ylabel(r"Kendall's $\tau$ (Median [IQR])")
    ax.tick_params(labelleft=True)
    ax.set_ylim(0.45, 1.02)
    ax.yaxis.set_major_locator(MultipleLocator(0.10))
    sample_ticks = sorted(df_sample["SampleSize"].dropna().unique())
    if sample_ticks:
        ax.set_xticks(sample_ticks)

    # (F) Candidate Scalability across Candidates N (N = 10, 12, 14)
    ax = axes[1, 1]
    df_cands = df[
        (df["SampleSize"] == DEFAULT_SAMPLE_SIZE) & 
        (df["Noise"] == DEFAULT_NOISE) & 
        (df["EffectMagnitude"] == DEFAULT_EFFECT) & 
        (df["Correlation"] == DEFAULT_CORRELATION)
    ]
    if df_cands.empty:
        df_cands = df_core
    plot_method_lines(ax, df_cands, "Candidates", "KendallTau", estimator=np.median, errorbar=("pi", 50))
    ax.set_title(r"(F) Scalability across Candidate Scale ($\mathbf{N}$)", fontweight="bold", loc="left")
    ax.set_xlabel(r"Number of Candidates ($N$)")
    ax.set_ylabel(r"Kendall's $\tau$ (Median [IQR])")
    ax.tick_params(labelleft=True)
    ax.set_ylim(0.45, 1.02)
    ax.yaxis.set_major_locator(MultipleLocator(0.10))
    cand_ticks = sorted(df_cands["Candidates"].dropna().unique())
    if cand_ticks:
        ax.set_xticks(cand_ticks)

    # (G) Effect Magnitude Sensitivity (Small, Medium, Large)
    ax = axes[1, 2]
    df_eff = df[
        (df["Candidates"] == ref_N) & 
        (df["SampleSize"] == DEFAULT_SAMPLE_SIZE) & 
        (df["Noise"] == DEFAULT_NOISE) & 
        (df["Correlation"] == DEFAULT_CORRELATION)
    ].copy()
    eff_labels = {
        "Small": "Small\n(d ≈ 0.25)",
        "Medium": "Medium\n(d ≈ 0.50)",
        "Large": "Large\n(d ≈ 0.80)"
    }
    df_eff["EffDisp"] = df_eff["EffectMagnitude"].map(eff_labels)
    eff_order = [eff_labels["Small"], eff_labels["Medium"], eff_labels["Large"]]
    for m in target_methods:
        sub_m = df_eff[df_eff["Method"] == m]
        if sub_m.empty:
            continue
        spec = METHOD_SPECS[m]
        meds = [sub_m[sub_m["EffDisp"] == el]["KendallTau"].median() for el in eff_order]
        ax.plot(
            eff_order, meds, label=m, color=spec["color"], linestyle=spec["linestyle"],
            linewidth=spec["linewidth"], marker=spec["marker"], markersize=spec["markersize"],
            markerfacecolor=spec["markerfacecolor"], markeredgecolor=spec["markeredgecolor"],
            markeredgewidth=spec["markeredgewidth"], zorder=spec["zorder"]
        )
    ax.set_title("(G) Rank Fidelity vs. Effect Size Magnitude", fontweight="bold", loc="left")
    ax.set_xlabel(r"Effect Size (Cliff's $d$)")
    ax.set_ylabel(r"Kendall's $\tau$ (Median)")
    ax.tick_params(labelleft=True)
    ax.set_ylim(0.45, 1.02)
    ax.yaxis.set_major_locator(MultipleLocator(0.10))

    # (H) Correlation Sensitivity (rho = 0.0 vs 0.5)
    ax = axes[1, 3]
    df_corr = df[
        (df["Candidates"] == ref_N) & 
        (df["SampleSize"] == DEFAULT_SAMPLE_SIZE) & 
        (df["Noise"] == DEFAULT_NOISE) & 
        (df["EffectMagnitude"] == DEFAULT_EFFECT)
    ].copy()
    corr_labels = {0.0: "Independent\n(" + r"$\rho = 0.0$" + ")", 0.5: "Correlated\n(" + r"$\rho = 0.5$" + ")"}
    df_corr["CorrDisp"] = df_corr["Correlation"].map(corr_labels)
    corr_order = [corr_labels[0.0], corr_labels[0.5]]
    for m in target_methods:
        sub_m = df_corr[df_corr["Method"] == m]
        if sub_m.empty:
            continue
        spec = METHOD_SPECS[m]
        meds = [sub_m[sub_m["CorrDisp"] == cl]["KendallTau"].median() for cl in corr_order]
        ax.plot(
            corr_order, meds, label=m, color=spec["color"], linestyle=spec["linestyle"],
            linewidth=spec["linewidth"], marker=spec["marker"], markersize=spec["markersize"],
            markerfacecolor=spec["markerfacecolor"], markeredgecolor=spec["markeredgecolor"],
            markeredgewidth=spec["markeredgewidth"], zorder=spec["zorder"]
        )
    ax.set_title(r"(H) Robustness under Metric Collinearity", fontweight="bold", loc="left")
    ax.set_xlabel(r"Inter-Metric Correlation ($\rho$)")
    ax.set_ylabel(r"Kendall's $\tau$ (Median)")
    ax.tick_params(labelleft=True)
    ax.set_ylim(0.45, 1.02)
    ax.yaxis.set_major_locator(MultipleLocator(0.10))

    # Strip internal subplot legends to prevent occluding data
    for a in axes.flat:
        if a.get_legend():
            a.get_legend().remove()

    # Unified external legend at bottom
    add_external_legend(fig, axes[0, 0], ncol=4, y_pos=-0.04)

    plt.tight_layout(rect=[0, 0.05, 1, 0.91])
    save_figure_dual(fig, output_path, pdf_path)
    plt.close(fig)


def plot_faceted_top_choice(
    df: pd.DataFrame, 
    output_path: Path, 
    pdf_path: Optional[Path] = None
) -> None:
    """Plots Top-Choice Recovery Rate faceted across candidate scales (N = 10, 12, 14) vs. Noise.

    Syntax:
        plot_faceted_top_choice(df, output_path, pdf_path=None)

    Description:
        Generates a 3-panel publication figure displaying Top-Choice Recovery Rate
        as a function of measurement noise (sigma in {2, 4, 6, 8, 10}%), faceted
        across candidate scales (N=10, 12, 14). Evaluates the frequency with which
        the true global optimum (C1) is correctly placed in rank 1.

    Parameters:
        df (pd.DataFrame): Raw Monte Carlo simulation results dataset.
        output_path (Path): Target file path for the 300 DPI PNG figure.
        pdf_path (Optional[Path]): Optional target path for vector PDF export.

    Returns:
        None

    Author:
        Lukas von Erdmannsdorff
    """
    df_core = df[(df["EffectMagnitude"] == "Medium") & (df["Correlation"] == 0.0)]
    candidates = sorted(df_core["Candidates"].unique())
    if not candidates:
        candidates = [10, 12, 14]
    
    num_cols = len(candidates)
    fig, axes = plt.subplots(1, num_cols, figsize=(4.8 * num_cols, 5.0), sharey=False)
    if num_cols == 1:
        axes = [axes]

    for idx, N in enumerate(candidates):
        ax = axes[idx]
        sub = df_core[df_core["Candidates"] == N]
        plot_method_lines(ax, sub, "Noise", "TopChoice")
        ax.set_title(f"N = {N} Candidates", fontweight="bold", fontsize=11.5)
        ax.set_xlabel(r"Noise Level ($\sigma$, %)", fontsize=11)
        ax.set_ylabel("Top-Choice Recovery (Mean [95% CI], %)", fontsize=11)
        ax.tick_params(labelleft=True)
        if ax.get_legend():
            ax.get_legend().remove()
        ax.set_ylim(-0.02, 1.05)
        ax.yaxis.set_major_locator(MultipleLocator(0.20))
        ax.yaxis.set_major_formatter(PercentFormatter(1.0, decimals=0))

    add_external_legend(fig, axes[0], ncol=4, y_pos=-0.06)
    cand_str = ", ".join(map(str, candidates))
    K = get_runs_per_point(df_core)
    fig.suptitle("Top-Choice Recovery Rate vs. Noise Level", fontweight="bold", fontsize=13.0, y=0.985)
    fig.text(
        0.5, 0.925,
        rf"$\mathbf{{N \in \{{{cand_str}\}}}}$, $\mathbf{{n = 50}}$, $\boldsymbol{{\rho = 0.0}}$, $\mathbf{{K = {K}}}$ runs per point",
        ha="center", va="top", fontsize=11.5, fontweight="bold", color="#222222"
    )
    plt.tight_layout(rect=[0, 0.06, 1, 0.86])
    save_figure_dual(fig, output_path, pdf_path)
    plt.close(fig)


def plot_faceted_complete_rank(
    df: pd.DataFrame, 
    output_path: Path, 
    pdf_path: Optional[Path] = None
) -> None:
    """Plots Complete-Rank Recovery Rate faceted across candidate scales (N = 10, 12, 14) vs. Noise.

    Syntax:
        plot_faceted_complete_rank(df, output_path, pdf_path=None)

    Description:
        Generates a 3-panel publication figure displaying Complete-Rank Recovery Rate
        (exact permutation match with ground truth across all N positions) as a function
        of noise level (sigma), faceted across candidate scales (N=10, 12, 14).

    Parameters:
        df (pd.DataFrame): Raw Monte Carlo simulation results dataset.
        output_path (Path): Target file path for the 300 DPI PNG figure.
        pdf_path (Optional[Path]): Optional target path for vector PDF export.

    Returns:
        None

    Author:
        Lukas von Erdmannsdorff
    """
    df_core = df[(df["EffectMagnitude"] == "Medium") & (df["Correlation"] == 0.0)]
    candidates = sorted(df_core["Candidates"].unique())
    if not candidates:
        candidates = [10, 12, 14]
    
    num_cols = len(candidates)
    fig, axes = plt.subplots(1, num_cols, figsize=(4.8 * num_cols, 5.0), sharey=False)
    if num_cols == 1:
        axes = [axes]

    for idx, N in enumerate(candidates):
        ax = axes[idx]
        sub = df_core[df_core["Candidates"] == N]
        plot_method_lines(ax, sub, "Noise", "CompleteRank")
        ax.set_title(f"N = {N} Candidates", fontweight="bold", fontsize=11.5)
        ax.set_xlabel(r"Noise Level ($\sigma$, %)", fontsize=11)
        ax.set_ylabel("Total Rank Recovery (Mean [95% CI], %)", fontsize=11)
        ax.tick_params(labelleft=True)
        if ax.get_legend():
            ax.get_legend().remove()
        ax.set_ylim(-0.02, 1.05)
        ax.yaxis.set_major_locator(MultipleLocator(0.20))
        ax.yaxis.set_major_formatter(PercentFormatter(1.0, decimals=0))

    add_external_legend(fig, axes[0], ncol=4, y_pos=-0.06)
    cand_str = ", ".join(map(str, candidates))
    K = get_runs_per_point(df_core)
    fig.suptitle("Total Rank Recovery vs. Noise Level", fontweight="bold", fontsize=13.0, y=0.985)
    fig.text(
        0.5, 0.925,
        rf"$\mathbf{{N \in \{{{cand_str}\}}}}$, $\mathbf{{n = 50}}$, $\boldsymbol{{\rho = 0.0}}$, $\mathbf{{K = {K}}}$ runs per point",
        ha="center", va="top", fontsize=11.5, fontweight="bold", color="#222222"
    )
    plt.tight_layout(rect=[0, 0.06, 1, 0.86])
    save_figure_dual(fig, output_path, pdf_path)
    plt.close(fig)


def plot_faceted_rank_displacement(
    df: pd.DataFrame, 
    output_path: Path, 
    pdf_path: Optional[Path] = None
) -> None:
    """Plots Mean Absolute Rank Displacement (|Rank - TrueRank|) faceted across candidate scales vs. Noise.

    Syntax:
        plot_faceted_rank_displacement(df, output_path, pdf_path=None)

    Description:
        Generates a 3-panel publication figure displaying the Mean Absolute Rank
        Displacement across all candidate positions from their true ground-truth ranks
        as a function of measurement noise (sigma), faceted across candidate scales (N=10, 12, 14).

    Parameters:
        df (pd.DataFrame): Raw Monte Carlo simulation results dataset.
        output_path (Path): Target file path for the 300 DPI PNG figure.
        pdf_path (Optional[Path]): Optional target path for vector PDF export.

    Returns:
        None

    Author:
        Lukas von Erdmannsdorff
    """
    df_core = df[(df["EffectMagnitude"] == "Medium") & (df["Correlation"] == 0.0)]
    candidates = sorted(df_core["Candidates"].unique())
    if not candidates:
        candidates = [10, 12, 14]

    num_cols = len(candidates)
    fig, axes = plt.subplots(1, num_cols, figsize=(4.8 * num_cols, 5.0), sharey=False)
    if num_cols == 1:
        axes = [axes]

    for idx, N in enumerate(candidates):
        ax = axes[idx]
        sub = df_core[df_core["Candidates"] == N]
        y_metric = "RankDisplacement" if "RankDisplacement" in sub.columns else "KendallTau"
        y_label = "Rank Displacement (Mean [95% CI])" if y_metric == "RankDisplacement" else r"Kendall's $\tau$ (Mean [95% CI])"
        plot_method_lines(ax, sub, "Noise", y_metric)
        ax.set_title(f"N = {N} Candidates", fontweight="bold", fontsize=11.5)
        ax.set_xlabel(r"Noise Level ($\sigma$, %)", fontsize=11)
        ax.set_ylabel(y_label, fontsize=11)
        ax.tick_params(labelleft=True)
        if y_metric == "RankDisplacement":
            ax.set_ylim(-0.05, 3.2)
            ax.yaxis.set_major_locator(MultipleLocator(0.5))
        if ax.get_legend():
            ax.get_legend().remove()

    add_external_legend(fig, axes[0], ncol=4, y_pos=-0.06)
    cand_str = ", ".join(map(str, candidates))
    K = get_runs_per_point(df_core)
    fig.suptitle("Mean Rank Displacement vs. Noise Level", fontweight="bold", fontsize=13.0, y=0.985)
    fig.text(
        0.5, 0.925,
        rf"$\mathbf{{N \in \{{{cand_str}\}}}}$, $\mathbf{{n = 50}}$, $\boldsymbol{{\rho = 0.0}}$, $\mathbf{{K = {K}}}$ runs per point",
        ha="center", va="top", fontsize=11.5, fontweight="bold", color="#222222"
    )
    plt.tight_layout(rect=[0, 0.06, 1, 0.86])
    save_figure_dual(fig, output_path, pdf_path)
    plt.close(fig)


def plot_faceted_regret(
    df: pd.DataFrame, 
    output_path: Path, 
    pdf_path: Optional[Path] = None
) -> None:
    """
    Backwards-compatibility routing to plot_faceted_rank_displacement.

    Syntax:
        plot_faceted_regret(df, output_path, pdf_path=None)

    Description:
        Provides backwards compatibility for legacy callers requesting regret plots
        by forwarding execution to `plot_faceted_rank_displacement`.

    Parameters:
        df (pd.DataFrame): Simulation results dataset.
        output_path (Path): Destination PNG file path.
        pdf_path (Optional[Path]): Optional destination vector PDF file path.

    Returns:
        None

    Author:
        Lukas von Erdmannsdorff
    """
    plot_faceted_rank_displacement(df, output_path, pdf_path)


def plot_faceted_false_superiority(
    df: pd.DataFrame, 
    output_path: Path, 
    pdf_path: Optional[Path] = None
) -> None:
    """Plots False Superiority Rate (Pairwise Inversions) faceted across candidate scales vs. Noise.

    Syntax:
        plot_faceted_false_superiority(df, output_path, pdf_path=None)

    Description:
        Generates a 3-panel publication figure displaying False Superiority Rate
        (the proportion of pairwise order inversions: Inversions / [N*(N-1)/2]) as a function
        of measurement noise (sigma), faceted across candidate scales (N=10, 12, 14).

    Parameters:
        df (pd.DataFrame): Raw Monte Carlo simulation results dataset.
        output_path (Path): Target file path for the 300 DPI PNG figure.
        pdf_path (Optional[Path]): Optional target path for vector PDF export.

    Returns:
        None

    Author:
        Lukas von Erdmannsdorff
    """
    df_core = df[(df["EffectMagnitude"] == "Medium") & (df["Correlation"] == 0.0)]
    candidates = sorted(df_core["Candidates"].unique())
    if not candidates:
        candidates = [10, 12, 14]
    
    num_cols = len(candidates)
    fig, axes = plt.subplots(1, num_cols, figsize=(4.8 * num_cols, 5.0), sharey=False)
    if num_cols == 1:
        axes = [axes]

    for idx, N in enumerate(candidates):
        ax = axes[idx]
        sub = df_core[df_core["Candidates"] == N]
        plot_method_lines(ax, sub, "Noise", "FalseSuperiority", estimator=np.median, errorbar=("pi", 50))
        ax.set_title(f"N = {N} Candidates", fontweight="bold", fontsize=11.5)
        ax.set_xlabel(r"Noise Level ($\sigma$, %)", fontsize=11)
        ax.set_ylabel("False Superiority (Median [IQR], %)", fontsize=11)
        ax.tick_params(labelleft=True)
        if ax.get_legend():
            ax.get_legend().remove()
        ax.set_ylim(-0.005, 0.32)
        ax.yaxis.set_major_locator(MultipleLocator(0.05))
        ax.yaxis.set_major_formatter(PercentFormatter(1.0, decimals=0))

    add_external_legend(fig, axes[0], ncol=4, y_pos=-0.06)
    cand_str = ", ".join(map(str, candidates))
    K = get_runs_per_point(df_core)
    fig.suptitle("Pairwise Inversion Rate (FSR) vs. Noise Level", fontweight="bold", fontsize=13.0, y=0.985)
    fig.text(
        0.5, 0.925,
        rf"$\mathbf{{N \in \{{{cand_str}\}}}}$, $\mathbf{{n = 50}}$, $\boldsymbol{{\rho = 0.0}}$, $\mathbf{{K = {K}}}$ runs per point",
        ha="center", va="top", fontsize=11.5, fontweight="bold", color="#222222"
    )
    plt.tight_layout(rect=[0, 0.06, 1, 0.86])
    save_figure_dual(fig, output_path, pdf_path)
    plt.close(fig)


def plot_faceted_rank_correlation(
    df: pd.DataFrame, 
    output_path: Path, 
    pdf_path: Optional[Path] = None
) -> None:
    """Plots Kendall's Tau rank correlation faceted across candidate scales vs. Sample Size.

    Syntax:
        plot_faceted_rank_correlation(df, output_path, pdf_path=None)

    Description:
        Generates a 3-panel publication figure displaying Kendall's Tau rank correlation
        as a function of sample size (n in {25, 50, 100}), faceted across
        candidate scales (N=10, 12, 14) at fixed baseline noise (sigma = 4.0%).

    Parameters:
        df (pd.DataFrame): Raw Monte Carlo simulation results dataset.
        output_path (Path): Target file path for the 300 DPI PNG figure.
        pdf_path (Optional[Path]): Optional target path for vector PDF export.

    Returns:
        None

    Author:
        Lukas von Erdmannsdorff
    """
    df_sample = df[
        (df["Noise"] == DEFAULT_NOISE) & 
        (df["EffectMagnitude"] == DEFAULT_EFFECT) & 
        (df["Correlation"] == DEFAULT_CORRELATION)
    ]
    candidates = sorted(df_sample["Candidates"].unique())
    if not candidates:
        candidates = NUM_CANDIDATES_LIST
    
    num_cols = len(candidates)
    fig, axes = plt.subplots(1, num_cols, figsize=(4.8 * num_cols, 5.0), sharey=False)
    if num_cols == 1:
        axes = [axes]

    for idx, N in enumerate(candidates):
        ax = axes[idx]
        sub = df_sample[df_sample["Candidates"] == N]
        plot_method_lines(ax, sub, "SampleSize", "KendallTau", estimator=np.median, errorbar=("pi", 50))
        ax.set_title(f"N = {N} Candidates", fontweight="bold", fontsize=11.5)
        ax.set_xlabel("Sample Size ($n$)", fontsize=11)
        ax.set_ylabel(r"Kendall's $\tau$ (Median [IQR])", fontsize=11)
        ax.tick_params(labelleft=True)
        sample_ticks = sorted(sub["SampleSize"].dropna().unique())
        if sample_ticks:
            ax.set_xticks(sample_ticks)
        if ax.get_legend():
            ax.get_legend().remove()
        ax.set_ylim(-0.05, 1.05)
        ax.yaxis.set_major_locator(MultipleLocator(0.20))

    add_external_legend(fig, axes[0], ncol=4, y_pos=-0.06)
    cand_str = ", ".join(map(str, candidates))
    K = get_runs_per_point(df_sample)
    fig.suptitle(r"Rank Correlation (Kendall's $\boldsymbol{\tau}$) vs. Sample Size", fontweight="bold", fontsize=13.0, y=0.985)
    fig.text(
        0.5, 0.925,
        rf"$\mathbf{{N \in \{{{cand_str}\}}}}$, $\boldsymbol{{\sigma = 4\%}}$, $\boldsymbol{{\rho = 0.0}}$, $\mathbf{{K = {K}}}$ runs per point",
        ha="center", va="top", fontsize=11.5, fontweight="bold", color="#222222"
    )
    plt.tight_layout(rect=[0, 0.06, 1, 0.86])
    save_figure_dual(fig, output_path, pdf_path)
    plt.close(fig)


def plot_effect_sensitivity(
    df: pd.DataFrame, 
    output_path: Path, 
    pdf_path: Optional[Path] = None,
    rank_data_path: Optional[Path] = None
) -> None:
    """Plots method robustness under Small, Medium, and Large Effect Magnitudes.

    Syntax:
        plot_effect_sensitivity(df, output_path, pdf_path=None, rank_data_path=None)

    Description:
        Generates a 6-panel publication figure (2 rows x 3 columns) visualizing:
          - Top Row (A1-A3): Candidate-Level Rank Dispersion (IQR = Q75 - Q25) across methods.
          - Bottom Row (B1-B3): Observed Assigned Rank (Median ± IQR) vs. Ground Truth along y = x.
        Panels evaluate performance across three standardized effect magnitude tiers:
          (1) Small Effect (Delta = 2.5%, Cliff's d ≈ 0.25 [0.20 – 0.30])
          (2) Medium Effect (Delta = 5.0%, Cliff's d ≈ 0.50 [0.45 – 0.60])
          (3) Large Effect (Delta = 8.0%, Cliff's d ≈ 0.80 [0.75 – 0.90])

    Parameters:
        df (pd.DataFrame): Raw Monte Carlo simulation results dataset.
        output_path (Path): Target file path for the 300 DPI PNG figure.
        pdf_path (Optional[Path]): Optional target path for vector PDF export.
        rank_data_path (Optional[Path]): Optional direct path to candidate_ranks.csv.

    Returns:
        None

    Author:
        Lukas von Erdmannsdorff
    """
    ref_N = DEFAULT_CANDIDATES if DEFAULT_CANDIDATES in df["Candidates"].values else df["Candidates"].iloc[0]

    # Search for candidate_ranks.csv / candidate_ranks_<timestamp>.csv
    csv_candidates = [rank_data_path] if rank_data_path else []
    for base in [output_path.parent, output_path.parent.parent, output_path.parent.parent / "CSVs", output_path.parent / "CSVs", DATA_DIR, DATA_DIR / "Simulation_Output", DATA_DIR / "Simulation_Output" / "CSVs"]:
        if base.exists():
            csv_candidates.extend(sorted(base.glob("CSVs/candidate_ranks*.csv")))
            csv_candidates.extend(sorted(base.glob("candidate_ranks*.csv")))
    target_csv = None
    for p in csv_candidates:
        if p and Path(p).exists() and os.path.getsize(p) > 0:
            target_csv = Path(p)
            break

    if target_csv:
        df_ranks = pd.read_csv(target_csv)
        sub = df_ranks[
            (df_ranks["Candidates"] == ref_N) & 
            (df_ranks["SampleSize"] == DEFAULT_SAMPLE_SIZE) & 
            (df_ranks["Noise"] == DEFAULT_NOISE) & 
            (df_ranks["Correlation"] == DEFAULT_CORRELATION)
        ]
        if not sub.empty:
            effects = [
                ("Small", "Small", "2.5", "0.25"),
                ("Medium", "Medium", "5.0", "0.50"),
                ("Large", "Large", "8.0", "0.80")
            ]

            fig, axes = plt.subplots(2, 3, figsize=(16.5, 9.6), squeeze=False)
            dodge = {"HERA": -0.22, "TOPSIS": -0.07, "Borda": 0.08, "Wilcoxon-Copeland": 0.23}
            target_methods = ["HERA", "TOPSIS", "Borda", "Wilcoxon-Copeland"]
            N = int(ref_N)
            candidates = [f"C{i}" for i in range(1, N + 1)]

            for col_idx, (eff_key, eff_name, delta, d_val) in enumerate(effects):
                sub_eff = sub[sub["EffectMagnitude"] == eff_key]

                # --- Top Row: Non-Parametric Rank Dispersion (IQR = Q75 - Q25) ---
                ax_top = axes[0, col_idx]
                iqr_records = []
                for m in target_methods:
                    sub_m = sub_eff[sub_eff["Method"] == m]
                    for c_idx, c in enumerate(candidates, start=1):
                        c_vals = sub_m[sub_m["Candidate"] == c]["PredictedRank"].dropna().values
                        if len(c_vals) > 0:
                            q25 = float(np.percentile(c_vals, 25))
                            q75 = float(np.percentile(c_vals, 75))
                            iqr = q75 - q25
                        else:
                            iqr = 0.0
                        iqr_records.append({
                            "Method": m,
                            "Candidate": c,
                            "TrueRank": c_idx,
                            "IQR": iqr
                        })

                df_iqr = pd.DataFrame(iqr_records)
                sns.barplot(
                    data=df_iqr, x="Candidate", y="IQR", hue="Method",
                    palette=COLORS, ax=ax_top, edgecolor="#333333", linewidth=0.6
                )
                ax_top.set_title(f"(A{col_idx+1}) Rank Dispersion: {eff_name} (Δ = {delta}%, d ≈ {d_val})", fontweight="bold", fontsize=10.5, loc="left")
                ax_top.set_xlabel(f"Candidate ($C_1 \\dots C_{{{N}}}$)", fontsize=10.5)
                ax_top.set_ylabel("Rank Dispersion (IQR)", fontsize=10.5)
                ax_top.tick_params(labelleft=True, labelsize=9.5)
                max_iqr = df_iqr["IQR"].max() if not df_iqr.empty else 3.0
                ax_top.set_ylim(-0.05, max(3.5, max_iqr + 0.5))
                ax_top.yaxis.set_major_locator(MultipleLocator(1.0))
                if ax_top.get_legend():
                    ax_top.get_legend().remove()

                # --- Bottom Row: Observed Assigned Rank (Median ± IQR) vs Ground Truth ---
                ax_bot = axes[1, col_idx]
                ax_bot.plot([1, N], [1, N], ls="--", color="#888888", lw=1.3, zorder=1)

                for m in target_methods:
                    sub_m = sub_eff[sub_eff["Method"] == m]
                    if sub_m.empty:
                        continue
                    spec = METHOD_SPECS[m]
                    x_vals, medians, yerr_l, yerr_u = [], [], [], []
                    for c_idx in range(1, N + 1):
                        c = f"C{c_idx}"
                        vals = sub_m[sub_m["Candidate"] == c]["PredictedRank"].dropna().values
                        if len(vals) > 0:
                            med = float(np.median(vals))
                            q25 = float(np.percentile(vals, 25))
                            q75 = float(np.percentile(vals, 75))
                        else:
                            med, q25, q75 = c_idx, c_idx, c_idx
                        x_vals.append(c_idx + dodge[m])
                        medians.append(med)
                        yerr_l.append(med - q25)
                        yerr_u.append(q75 - med)

                    ax_bot.errorbar(
                        x_vals, medians, yerr=[yerr_l, yerr_u],
                        fmt=spec["marker"], color=spec["color"], label=m,
                        capsize=3.5, capthick=1.0, elinewidth=1.2, markersize=spec["markersize"],
                        markeredgecolor=spec["markeredgecolor"], markeredgewidth=1.0,
                        markerfacecolor=spec["markerfacecolor"], zorder=spec["zorder"]
                    )

                ax_bot.set_title(f"(B{col_idx+1}) Assigned vs. True Rank: {eff_name} (Δ = {delta}%, d ≈ {d_val})", fontweight="bold", fontsize=10.5, loc="left")
                ax_bot.set_xlabel(f"Ground-Truth Rank ($1 \\dots {N}$)", fontsize=10.5)
                ax_bot.set_ylabel(r"Assigned Rank (Median $\pm$ IQR)", fontsize=10.5)
                ax_bot.set_xticks(range(1, N + 1))
                ax_bot.set_yticks(range(1, N + 1))
                ax_bot.set_xlim(0.4, N + 0.6)
                ax_bot.set_ylim(0.4, N + 0.6)
                ax_bot.tick_params(labelleft=True, labelsize=9.5)
                if ax_bot.get_legend():
                    ax_bot.get_legend().remove()

            add_external_legend(fig, axes[0, 0], ncol=4, y_pos=-0.04)
            K = get_runs_per_point(sub)
            fig.suptitle("Effect Size Sensitivity: Candidate-Level Rank Dispersion & Fidelity", fontweight="bold", fontsize=13.5, y=0.985)
            fig.text(
                0.5, 0.950,
                rf"$\mathbf{{N = {N}}}$, $\mathbf{{n = 50}}$, $\boldsymbol{{\sigma = 4\%}}$, $\boldsymbol{{\rho = 0.0}}$, $\mathbf{{K = {K}}}$ runs per condition",
                ha="center", va="top", fontsize=11.5, fontweight="bold", color="#222222"
            )
            plt.tight_layout(rect=[0, 0.05, 1, 0.91])
            save_figure_dual(fig, output_path, pdf_path)
            plt.close(fig)
            return

    # Fallback: TopChoice and RankDisplacement barplots if candidate_ranks.csv unavailable
    sub_df = df[
        (df["Candidates"] == ref_N) & 
        (df["SampleSize"] == DEFAULT_SAMPLE_SIZE) & 
        (df["Noise"] == DEFAULT_NOISE) & 
        (df["Correlation"] == DEFAULT_CORRELATION)
    ]
    if sub_df.empty:
        print(f"    [Notice] Dataset subset for Effect Sensitivity (N={ref_N}, n={DEFAULT_SAMPLE_SIZE}, sigma={DEFAULT_NOISE}) not present. Skipped.")
        return

    effect_labels = {
        "Small": "Small\n(d ≈ 0.25 [0.20 - 0.30])",
        "Medium": "Medium\n(d ≈ 0.50 [0.45 - 0.60])",
        "Large": "Large\n(d ≈ 0.80 [0.75 - 0.90])"
    }
    sub_df = sub_df.copy()
    sub_df["EffectDisplay"] = sub_df["EffectMagnitude"].map(effect_labels)
    order_display = [effect_labels[k] for k in ["Small", "Medium", "Large"] if k in effect_labels]

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.0))
    
    # Panel 1: Top-Choice Recovery
    sns.barplot(
        data=sub_df, x="EffectDisplay", y="TopChoice", hue="Method", 
        palette=COLORS, order=order_display, ax=axes[0], errorbar=("ci", 95), capsize=0.08
    )
    axes[0].set_title(r"(A) Top-Choice Recovery Rate", fontweight="bold", loc="left")
    axes[0].set_xlabel("Effect Size (Median Cliff's d [Range])", fontsize=10.5)
    axes[0].set_ylabel("Top-Choice Recovery (Mean [95% CI], %)")
    axes[0].tick_params(labelleft=True)
    axes[0].set_ylim(0, 1.05)
    axes[0].yaxis.set_major_locator(MultipleLocator(0.20))
    axes[0].yaxis.set_major_formatter(PercentFormatter(1.0, decimals=0))
    if axes[0].get_legend():
        axes[0].get_legend().remove()

    # Panel 2: Rank Displacement
    y_metric = "RankDisplacement" if "RankDisplacement" in sub_df.columns else "KendallTau"
    y_label = "Rank Displacement (Mean [95% CI])" if y_metric == "RankDisplacement" else r"Kendall's $\tau$ (Mean [95% CI])"
    y_title = r"(B) Mean Rank Displacement (|Rank - TrueRank|)" if y_metric == "RankDisplacement" else r"(B) Rank Fidelity (Kendall's $\tau$)"
    sns.barplot(
        data=sub_df, x="EffectDisplay", y=y_metric, hue="Method", 
        palette=COLORS, order=order_display, ax=axes[1], errorbar=("ci", 95), capsize=0.08
    )
    axes[1].set_title(y_title, fontweight="bold", loc="left")
    axes[1].set_xlabel("Effect Size (Median Cliff's d [Range])", fontsize=10.5)
    axes[1].set_ylabel(y_label)
    axes[1].tick_params(labelleft=True)
    if y_metric == "RankDisplacement":
        axes[1].set_ylim(0, 2.5)
        axes[1].yaxis.set_major_locator(MultipleLocator(0.5))
    if axes[1].get_legend():
        axes[1].get_legend().remove()

    add_external_legend(fig, axes[0], ncol=4, y_pos=-0.07)
    K = get_runs_per_point(sub_df)
    fig.suptitle("Effect Size Sensitivity", fontweight="bold", fontsize=13.0, y=0.985)
    fig.text(
        0.5, 0.925,
        rf"$\mathbf{{N = {ref_N}}}$, $\mathbf{{n = 50}}$, $\boldsymbol{{\sigma = 4\%}}$, $\boldsymbol{{\rho = 0.0}}$, $\mathbf{{K = {K}}}$ runs per point",
        ha="center", va="top", fontsize=11.5, fontweight="bold", color="#222222"
    )
    plt.tight_layout(rect=[0, 0.08, 1, 0.86])
    save_figure_dual(fig, output_path, pdf_path)
    plt.close(fig)


def plot_correlation_sensitivity(
    df: pd.DataFrame, 
    output_path: Path, 
    pdf_path: Optional[Path] = None,
    rank_data_path: Optional[Path] = None
) -> None:
    """Plots method robustness under Inter-Metric Collinearity (rho = 0.0 vs. 0.5).

    Syntax:
        plot_correlation_sensitivity(df, output_path, pdf_path=None, rank_data_path=None)

    Description:
        Generates a 4-panel publication figure (2 rows x 2 columns) visualizing:
          - Top Row (A1-A2): Candidate-Level Rank Dispersion (IQR = Q75 - Q25) across methods.
          - Bottom Row (B1-B2): Observed Assigned Rank (Median ± IQR) vs. Ground Truth along y = x.
        Panels evaluate performance under two correlation structures:
          (1) Independent Criteria (rho = 0.0, Orthogonal Multi-Criteria Benchmark)
          (2) Correlated Criteria (rho = 0.5, Moderate Collinearity)

    Parameters:
        df (pd.DataFrame): Raw Monte Carlo simulation results dataset.
        output_path (Path): Target file path for the 300 DPI PNG figure.
        pdf_path (Optional[Path]): Optional target path for vector PDF export.
        rank_data_path (Optional[Path]): Optional direct path to candidate_ranks.csv.

    Returns:
        None

    Author:
        Lukas von Erdmannsdorff
    """
    ref_N = DEFAULT_CANDIDATES if DEFAULT_CANDIDATES in df["Candidates"].values else df["Candidates"].iloc[0]

    # Search for candidate_ranks.csv / candidate_ranks_<timestamp>.csv
    csv_candidates = [rank_data_path] if rank_data_path else []
    for base in [output_path.parent, output_path.parent.parent, output_path.parent.parent / "CSVs", output_path.parent / "CSVs", DATA_DIR, DATA_DIR / "Simulation_Output", DATA_DIR / "Simulation_Output" / "CSVs"]:
        if base.exists():
            csv_candidates.extend(sorted(base.glob("CSVs/candidate_ranks*.csv")))
            csv_candidates.extend(sorted(base.glob("candidate_ranks*.csv")))
    target_csv = None
    for p in csv_candidates:
        if p and Path(p).exists() and os.path.getsize(p) > 0:
            target_csv = Path(p)
            break

    if target_csv:
        df_ranks = pd.read_csv(target_csv)
        sub = df_ranks[
            (df_ranks["Candidates"] == ref_N) & 
            (df_ranks["SampleSize"] == DEFAULT_SAMPLE_SIZE) & 
            (df_ranks["Noise"] == DEFAULT_NOISE) & 
            (df_ranks["EffectMagnitude"] == DEFAULT_EFFECT)
        ]
        if not sub.empty and len(sub["Correlation"].dropna().unique()) > 1:
            corr_conditions = [
                (0.0, r"Independent Criteria ($\boldsymbol{\rho = 0.0}$)"),
                (0.5, r"Correlated Criteria ($\boldsymbol{\rho = 0.5}$)")
            ]

            fig, axes = plt.subplots(2, 2, figsize=(11.8, 9.6), squeeze=False)
            dodge = {"HERA": -0.22, "TOPSIS": -0.07, "Borda": 0.08, "Wilcoxon-Copeland": 0.23}
            target_methods = ["HERA", "TOPSIS", "Borda", "Wilcoxon-Copeland"]
            N = int(ref_N)
            candidates = [f"C{i}" for i in range(1, N + 1)]

            for col_idx, (rho_val, title_condition) in enumerate(corr_conditions):
                sub_corr = sub[sub["Correlation"] == rho_val]

                # --- Top Row: Non-Parametric Rank Dispersion (IQR = Q75 - Q25) ---
                ax_top = axes[0, col_idx]
                iqr_records = []
                for m in target_methods:
                    sub_m = sub_corr[sub_corr["Method"] == m]
                    for c_idx, c in enumerate(candidates, start=1):
                        c_vals = sub_m[sub_m["Candidate"] == c]["PredictedRank"].dropna().values
                        if len(c_vals) > 0:
                            q25 = float(np.percentile(c_vals, 25))
                            q75 = float(np.percentile(c_vals, 75))
                            iqr = q75 - q25
                        else:
                            iqr = 0.0
                        iqr_records.append({
                            "Method": m,
                            "Candidate": c,
                            "TrueRank": c_idx,
                            "IQR": iqr
                        })

                df_iqr = pd.DataFrame(iqr_records)
                sns.barplot(
                    data=df_iqr, x="Candidate", y="IQR", hue="Method",
                    palette=COLORS, ax=ax_top, edgecolor="#333333", linewidth=0.6
                )
                ax_top.set_title(f"(A{col_idx+1}) Rank Dispersion: {title_condition}", fontweight="bold", fontsize=11.0, loc="left")
                ax_top.set_xlabel(f"Candidate ($C_1 \\dots C_{{{N}}}$)", fontsize=10.5)
                ax_top.set_ylabel("Rank Dispersion (IQR)", fontsize=10.5)
                ax_top.tick_params(labelleft=True, labelsize=9.5)
                max_iqr = df_iqr["IQR"].max() if not df_iqr.empty else 3.0
                ax_top.set_ylim(-0.05, max(3.5, max_iqr + 0.5))
                ax_top.yaxis.set_major_locator(MultipleLocator(1.0))
                if ax_top.get_legend():
                    ax_top.get_legend().remove()

                # --- Bottom Row: Observed Assigned Rank (Median ± IQR) vs Ground Truth ---
                ax_bot = axes[1, col_idx]
                ax_bot.plot([1, N], [1, N], ls="--", color="#888888", lw=1.3, zorder=1)

                for m in target_methods:
                    sub_m = sub_corr[sub_corr["Method"] == m]
                    if sub_m.empty:
                        continue
                    spec = METHOD_SPECS[m]
                    x_vals, medians, yerr_l, yerr_u = [], [], [], []
                    for c_idx in range(1, N + 1):
                        c = f"C{c_idx}"
                        vals = sub_m[sub_m["Candidate"] == c]["PredictedRank"].dropna().values
                        if len(vals) > 0:
                            med = float(np.median(vals))
                            q25 = float(np.percentile(vals, 25))
                            q75 = float(np.percentile(vals, 75))
                        else:
                            med, q25, q75 = c_idx, c_idx, c_idx
                        x_vals.append(c_idx + dodge[m])
                        medians.append(med)
                        yerr_l.append(med - q25)
                        yerr_u.append(q75 - med)

                    ax_bot.errorbar(
                        x_vals, medians, yerr=[yerr_l, yerr_u],
                        fmt=spec["marker"], color=spec["color"], label=m,
                        capsize=3.5, capthick=1.0, elinewidth=1.2, markersize=spec["markersize"],
                        markeredgecolor=spec["markeredgecolor"], markeredgewidth=1.0,
                        markerfacecolor=spec["markerfacecolor"], zorder=spec["zorder"]
                    )

                ax_bot.set_title(f"(B{col_idx+1}) Assigned Rank vs. Ground Truth: {title_condition}", fontweight="bold", fontsize=11.0, loc="left")
                ax_bot.set_xlabel(f"Ground-Truth Rank ($1 \\dots {N}$)", fontsize=10.5)
                ax_bot.set_ylabel(r"Assigned Rank (Median $\pm$ IQR)", fontsize=10.5)
                ax_bot.set_xticks(range(1, N + 1))
                ax_bot.set_yticks(range(1, N + 1))
                ax_bot.set_xlim(0.4, N + 0.6)
                ax_bot.set_ylim(0.4, N + 0.6)
                ax_bot.tick_params(labelleft=True, labelsize=9.5)
                if ax_bot.get_legend():
                    ax_bot.get_legend().remove()

            add_external_legend(fig, axes[0, 0], ncol=4, y_pos=-0.04)
            K = get_runs_per_point(sub)
            fig.suptitle("Inter-Metric Collinearity Sensitivity: Rank Dispersion & Fidelity", fontweight="bold", fontsize=13.5, y=0.985)
            fig.text(
                0.5, 0.950,
                rf"$\mathbf{{N = {N}}}$, $\mathbf{{n = 50}}$, $\boldsymbol{{\sigma = 4\%}}$, Cliff's $\mathbf{{d \approx 0.50}}$, $\mathbf{{K = {K}}}$ runs per condition",
                ha="center", va="top", fontsize=11.5, fontweight="bold", color="#222222"
            )
            plt.tight_layout(rect=[0, 0.05, 1, 0.91])
            save_figure_dual(fig, output_path, pdf_path)
            plt.close(fig)
            return

    # Fallback: TopChoice barplot if candidate_ranks.csv unavailable
    sub = df[
        (df["Candidates"] == ref_N) & 
        (df["SampleSize"] == DEFAULT_SAMPLE_SIZE) & 
        (df["Noise"] == DEFAULT_NOISE) & 
        (df["EffectMagnitude"] == DEFAULT_EFFECT)
    ]
    if sub.empty or "Correlation" not in sub.columns:
        print(f"    [Notice] Dataset subset for Correlation Sensitivity (N={ref_N}, n={DEFAULT_SAMPLE_SIZE}, sigma={DEFAULT_NOISE}) not present. Skipped.")
        return

    sub = sub.copy()
    label_map = {0.0: r"Independent ($\rho = 0.0$)", 0.5: r"Correlated ($\rho = 0.5$)"}
    sub["CorrelationLabel"] = sub["Correlation"].map(label_map)
    order = [r"Independent ($\rho = 0.0$)", r"Correlated ($\rho = 0.5$)"]

    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    sns.barplot(
        data=sub, x="CorrelationLabel", y="TopChoice", hue="Method", 
        order=order, palette=COLORS, ax=ax, errorbar=("ci", 95), capsize=0.08
    )
    ax.set_xlabel("Inter-Metric Correlation Condition")
    ax.set_ylabel("Top-Choice Recovery (Mean [95% CI], %)")
    ax.tick_params(labelleft=True)
    ax.set_ylim(0, 1.05)
    ax.yaxis.set_major_locator(MultipleLocator(0.20))
    ax.yaxis.set_major_formatter(PercentFormatter(1.0, decimals=0))
    if ax.get_legend():
        ax.get_legend().remove()

    add_external_legend(fig, ax, ncol=4, y_pos=-0.08)
    K = get_runs_per_point(sub)
    fig.suptitle("Inter-Metric Collinearity Sensitivity", fontweight="bold", fontsize=13.0, y=0.985)
    fig.text(
        0.5, 0.925,
        rf"$\mathbf{{N = {ref_N}}}$, $\mathbf{{n = 50}}$, $\boldsymbol{{\sigma = 4\%}}$, $\mathbf{{K = {K}}}$ runs per condition",
        ha="center", va="top", fontsize=11.5, fontweight="bold", color="#222222"
    )
    plt.tight_layout(rect=[0, 0.08, 1, 0.86])
    save_figure_dual(fig, output_path, pdf_path)
    plt.close(fig)


def plot_compensatory_failure(
    df: pd.DataFrame, 
    output_path: Path, 
    pdf_path: Optional[Path] = None
) -> None:
    """
    Plots Compensatory Pathology Stress Test (Protection vs. MCDM Failure).

    Syntax:
        plot_compensatory_failure(df, output_path, pdf_path)

    Description:
        Evaluates candidate ranking behavior under deliberate metric defect injection
        (non-compensatory stress test). Specifically contrasts safe model recovery
        (C1 as Rank 1), compensatory failure rate (selection of flawed model C_N into
        the upper half despite critical failure on a gatekeeper metric), and rank displacement.

    Parameters:
        df (pd.DataFrame): Simulation trial dataset containing compensatory benchmark records.
        output_path (Path): Destination filesystem path for the output PNG graphic.
        pdf_path (Optional[Path]): Optional destination path for standalone vector PDF.

    Returns:
        None

    Author:
        Lukas von Erdmannsdorff
    """
    if "Suite" not in df.columns or not any(df["Suite"].astype(str).str.contains("Compensatory", case=False)):
        print("    [Notice] No Compensatory suite present in dataset. Skipped Compensatory plot.")
        return

    ref_N = DEFAULT_CANDIDATES if DEFAULT_CANDIDATES in df["Candidates"].values else df["Candidates"].iloc[0]
    sub = df[
        (df["Suite"].astype(str).str.contains("Compensatory", case=False)) &
        (df["Candidates"] == ref_N) & 
        (df["SampleSize"] == DEFAULT_SAMPLE_SIZE) & 
        (df["EffectMagnitude"] == DEFAULT_EFFECT) & 
        (df["Correlation"] == DEFAULT_CORRELATION)
    ]
    if sub.empty:
        print("    [Notice] Compensatory data subset not present. Skipped Compensatory plot.")
        return

    ref_N = sub["Candidates"].iloc[0] if "Candidates" in sub.columns else 11

    fig, axes = plt.subplots(1, 3, figsize=(15.5, 5.0))

    # Panel 1: Top-Choice Recovery Rate (Recovery of True Safe Model C1)
    ax = axes[0]
    plot_method_lines(ax, sub, "Noise", "TopChoice")
    ax.set_title(r"(A) Safe Model Recovery ($\mathbf{C_1}$ as Rank 1)", fontweight="bold", loc="left")
    ax.set_xlabel(r"Noise Level ($\sigma$, %)")
    ax.set_ylabel("Top-Choice Recovery (Mean [95% CI], %)")
    ax.tick_params(labelleft=True)
    ax.set_ylim(-0.02, 1.05)
    ax.yaxis.set_major_locator(MultipleLocator(0.20))
    ax.yaxis.set_major_formatter(PercentFormatter(1.0, decimals=0))
    if ax.get_legend():
        ax.get_legend().remove()

    # Panel 2: Compensatory Error Rate (selection of Flawed Model C_N in top half)
    ax = axes[1]
    plot_method_lines(ax, sub, "Noise", "CompensatoryError")
    ax.set_title(r"(B) Compensatory Failure Rate ($\mathbf{C}_{\mathbf{flawed}}$ in Top Half)", fontweight="bold", loc="left")
    ax.set_xlabel(r"Noise Level ($\sigma$, %)")
    ax.set_ylabel("Compensatory Failure (Mean [95% CI], %)")
    ax.tick_params(labelleft=True)
    ax.set_ylim(-0.02, 1.05)
    ax.yaxis.set_major_locator(MultipleLocator(0.20))
    ax.yaxis.set_major_formatter(PercentFormatter(1.0, decimals=0))
    if ax.get_legend():
        ax.get_legend().remove()

    # Panel 3: Mean Rank Displacement (|Rank - TrueRank|)
    ax = axes[2]
    y_metric = "RankDisplacement" if "RankDisplacement" in sub.columns else "KendallTau"
    y_title = r"(C) Mean Rank Displacement ($\mathbf{|Rank - TrueRank|}$)" if y_metric == "RankDisplacement" else r"(C) Kendall's $\boldsymbol{\tau}$"
    y_label = "Rank Displacement (Mean [95% CI])" if y_metric == "RankDisplacement" else r"Kendall's $\tau$ (Mean [95% CI])"
    plot_method_lines(ax, sub, "Noise", y_metric)
    ax.set_title(y_title, fontweight="bold", loc="left")
    ax.set_xlabel(r"Noise Level ($\sigma$, %)")
    ax.set_ylabel(y_label)
    ax.tick_params(labelleft=True)
    if y_metric == "RankDisplacement":
        ax.set_ylim(-0.05, 3.2)
        ax.yaxis.set_major_locator(MultipleLocator(0.5))
    if ax.get_legend():
        ax.get_legend().remove()

    add_external_legend(fig, axes[0], ncol=4, y_pos=-0.06)
    K = get_runs_per_point(sub)
    fig.suptitle("Non-Compensatory Stress Test", fontweight="bold", fontsize=13.0, y=0.985)
    fig.text(
        0.5, 0.925,
        rf"$\mathbf{{N = {ref_N}}}$, $\mathbf{{n = 50}}$, $\boldsymbol{{\rho = 0.0}}$, $\mathbf{{K = {K}}}$ runs per point",
        ha="center", va="top", fontsize=11.5, fontweight="bold", color="#222222"
    )
    plt.tight_layout(rect=[0, 0.06, 1, 0.86])
    save_figure_dual(fig, output_path, pdf_path)
    plt.close(fig)


def plot_candidate_rank_stability(
    output_path: Path, 
    rank_data_path: Optional[Path] = None, 
    pdf_path: Optional[Path] = None
) -> None:
    """
    Plots Middle-Field Rank Stability across candidate scales using Non-Parametric Median and IQR.

    Syntax:
        plot_candidate_rank_stability(output_path, rank_data_path, pdf_path)

    Description:
        Renders candidate-level rank dispersion and assigned rank fidelity across candidate
        scales (N in {10, 12, 14}).
        Row 1 shows rank dispersion as Interquartile Range (IQR = Q75 - Q25), exposing
        rank volatility in the crowded middle field.
        Row 2 shows observed assigned ranks (Median +/- IQR error bars) against the ground-truth
        diagonal (y = x).

    Parameters:
        output_path (Path): Destination filesystem path for the output PNG graphic.
        rank_data_path (Optional[Path]): Path to candidate_ranks.csv recording trial-level rankings.
        pdf_path (Optional[Path]): Optional destination path for standalone vector PDF.

    Returns:
        None

    Author:
        Lukas von Erdmannsdorff
    """
    csv_candidates = [
        rank_data_path,
        output_path.parent / "candidate_ranks.csv",
        output_path.parent / "CSVs" / "candidate_ranks.csv",
        output_path.parent.parent / "CSVs" / "candidate_ranks.csv",
        output_path.parent.parent / "candidate_ranks.csv",
        output_path.parent / "candidate_rank_stability.csv",
        output_path.parent.parent / "candidate_rank_stability.csv",
        DATA_DIR / "candidate_ranks.csv",
        DATA_DIR / "CSVs" / "candidate_ranks.csv",
        DATA_DIR / "Simulation_Output" / "candidate_ranks.csv",
        DATA_DIR / "Simulation_Output" / "CSVs" / "candidate_ranks.csv"
    ]
    for base in [output_path.parent, output_path.parent.parent]:
        if base.exists():
            csv_candidates.extend(sorted(base.glob("CSVs/candidate_ranks*.csv")))
            csv_candidates.extend(sorted(base.glob("candidate_ranks*.csv")))

    target_csv = None
    for p in csv_candidates:
        if p and Path(p).exists() and os.path.getsize(p) > 0:
            target_csv = Path(p)
            break

    if not target_csv:
        print("    [Notice] candidate_ranks.csv / candidate_rank_stability.csv not found. Skipped Figure 10.")
        return

    df_ranks = pd.read_csv(target_csv)
    if df_ranks.empty or "TrueRank" not in df_ranks.columns or df_ranks["TrueRank"].dropna().empty:
        print("    [Notice] candidate_ranks.csv contains no valid rank data. Skipped Figure 10.")
        return

    target_methods = ["HERA", "TOPSIS", "Borda", "Wilcoxon-Copeland"]
    df_ranks = df_ranks[df_ranks["Method"].isin(target_methods)]
    if df_ranks.empty:
        print("    [Notice] candidate_ranks.csv contains none of the target methods. Skipped Figure 10.")
        return

    # Detect candidate scales present in dataset
    scales = sorted(df_ranks["Candidates"].unique()) if "Candidates" in df_ranks.columns else []
    if not scales:
        max_rank = int(df_ranks["TrueRank"].max())
        scales = [max_rank]

    # Display target scales (10, 12, 14 if available, otherwise existing scales)
    desired_scales = [s for s in [10, 12, 14] if s in scales]
    if not desired_scales:
        desired_scales = scales[:3] if len(scales) >= 3 else scales

    num_cols = len(desired_scales)
    fig, axes = plt.subplots(2, num_cols, figsize=(5.5 * num_cols, 9.6), squeeze=False)

    dodge_offsets = {"HERA": -0.22, "TOPSIS": -0.07, "Borda": 0.08, "Wilcoxon-Copeland": 0.23}

    for col_idx, N in enumerate(desired_scales):
        sub_scale = df_ranks[df_ranks["Candidates"] == N] if "Candidates" in df_ranks.columns else df_ranks
        candidates = [f"C{i}" for i in range(1, N + 1)]

        # --- Top Row: Non-Parametric Rank Dispersion (IQR = Q75 - Q25) ---
        ax_top = axes[0, col_idx]
        
        # Calculate IQR per candidate and method
        iqr_records = []
        for m in target_methods:
            sub_m = sub_scale[sub_scale["Method"] == m]
            for c_idx, c in enumerate(candidates, start=1):
                c_vals = sub_m[sub_m["Candidate"] == c]["PredictedRank"].dropna().values
                if len(c_vals) > 0:
                    q25 = float(np.percentile(c_vals, 25))
                    q75 = float(np.percentile(c_vals, 75))
                    iqr = q75 - q25
                else:
                    iqr = 0.0
                iqr_records.append({
                    "Method": m,
                    "Candidate": c,
                    "TrueRank": c_idx,
                    "IQR": iqr
                })
        
        df_iqr = pd.DataFrame(iqr_records)
        sns.barplot(
            data=df_iqr, x="Candidate", y="IQR", hue="Method",
            palette=COLORS, ax=ax_top, edgecolor="#333333", linewidth=0.6
        )
        ax_top.set_title(f"(A{col_idx+1}) Rank Dispersion (IQR, N = {N})", fontweight="bold", fontsize=11.5, loc="left")
        ax_top.set_xlabel(f"Candidate ($C_1 \\dots C_{{{N}}}$)", fontsize=11)
        ax_top.set_ylabel("Rank Dispersion (IQR)", fontsize=11)
        ax_top.tick_params(labelleft=True, labelsize=9.5)
        max_iqr = df_iqr["IQR"].max() if not df_iqr.empty else 3.0
        ax_top.set_ylim(-0.05, max(3.5, max_iqr + 0.5))
        ax_top.yaxis.set_major_locator(MultipleLocator(1.0))
        if ax_top.get_legend():
            ax_top.get_legend().remove()

        # --- Bottom Row: Observed Assigned Rank (Median ± IQR Error Bars) vs Ground Truth ---
        ax_bot = axes[1, col_idx]
        # Diagonal ground-truth reference line
        ax_bot.plot([1, N], [1, N], ls="--", color="#888888", lw=1.3, zorder=1, label="Ground Truth Rank ($y = x$)")

        for m in target_methods:
            sub_m = sub_scale[sub_scale["Method"] == m]
            if sub_m.empty:
                continue
            spec = METHOD_SPECS[m]
            x_vals = []
            medians = []
            yerr_lowers = []
            yerr_uppers = []

            for c_idx, c in enumerate(candidates, start=1):
                c_vals = sub_m[sub_m["Candidate"] == c]["PredictedRank"].dropna().values
                if len(c_vals) > 0:
                    med = float(np.median(c_vals))
                    q25 = float(np.percentile(c_vals, 25))
                    q75 = float(np.percentile(c_vals, 75))
                else:
                    med = float(c_idx)
                    q25 = float(c_idx)
                    q75 = float(c_idx)

                x_vals.append(c_idx + dodge_offsets[m])
                medians.append(med)
                yerr_lowers.append(med - q25)
                yerr_uppers.append(q75 - med)

            ax_bot.errorbar(
                x_vals, medians, yerr=[yerr_lowers, yerr_uppers],
                fmt=spec["marker"], color=spec["color"], label=m,
                capsize=3.5, capthick=1.0, elinewidth=1.2, markersize=spec["markersize"],
                markeredgecolor=spec["markeredgecolor"], markeredgewidth=1.0,
                markerfacecolor=spec["markerfacecolor"], zorder=spec["zorder"]
            )

        ax_bot.set_title(f"(B{col_idx+1}) Assigned Rank vs. Ground Truth (N = {N})", fontweight="bold", fontsize=11.5, loc="left")
        ax_bot.set_xlabel(f"Ground-Truth Rank ($1 \\dots {N}$)", fontsize=11)
        ax_bot.set_ylabel("Assigned Rank (Median $\\pm$ IQR)", fontsize=11)
        ax_bot.set_xticks(range(1, N + 1))
        ax_bot.set_yticks(range(1, N + 1))
        ax_bot.set_xlim(0.4, N + 0.6)
        ax_bot.set_ylim(0.4, N + 0.6)
        ax_bot.tick_params(labelleft=True, labelsize=9.5)
        if ax_bot.get_legend():
            ax_bot.get_legend().remove()

    add_external_legend(fig, axes[0, 0], ncol=4, y_pos=-0.04)
    cand_str = ", ".join(map(str, desired_scales))
    K = get_runs_per_point(df_ranks)
    fig.suptitle("Candidate-Level Rank Stability Across Scales", fontweight="bold", fontsize=13.5, y=0.985)
    fig.text(
        0.5, 0.950,
        rf"$\mathbf{{N \in \{{{cand_str}\}}}}$, $\mathbf{{n = 50}}$, $\boldsymbol{{\sigma = 4\%}}$, $\boldsymbol{{\rho = 0.0}}$, $\mathbf{{K = {K}}}$ runs per point",
        ha="center", va="top", fontsize=11.5, fontweight="bold", color="#222222"
    )
    plt.tight_layout(rect=[0, 0.05, 1, 0.91])
    save_figure_dual(fig, output_path, pdf_path)
    plt.close(fig)


def plot_pooled_core_distributions(
    df: pd.DataFrame, 
    output_path: Path, 
    pdf_path: Optional[Path] = None,
    df_ranks: Optional[pd.DataFrame] = None
) -> Optional[plt.Figure]:
    """
    Plots master non-parametric pooled distributions across the Core benchmark conditions.

    Syntax:
        fig = plot_pooled_core_distributions(df, output_path, pdf_path, df_ranks)

    Description:
        Constructs a comprehensive 2x3 panel layout summarizing overall performance across
        all core benchmark conditions:
          (A) Top-Choice Recovery Rate (%) [Barplot with 95% CI]
          (B) Complete-Rank Recovery Rate (%) [Barplot with 95% CI]
          (C) Pairwise Inversion Rate / False Superiority Rate (%) [Boxplot]
          (D) Rank Fidelity (Kendall's tau) [Boxplot]
          (E) Empirical Rank Stability across Repeated Runs (Candidate-IQR) [Boxplot]
          (F) Monotonic Rank Alignment (Spearman's rho) [Boxplot]

    Parameters:
        df (pd.DataFrame): Simulation results dataset.
        output_path (Path): Destination filesystem path for the output PNG graphic.
        pdf_path (Optional[Path]): Optional destination path for standalone vector PDF.
        df_ranks (Optional[pd.DataFrame]): Optional candidate ranks DataFrame for Panel E.

    Returns:
        Optional[plt.Figure]: Rendered matplotlib Figure instance (closed before return).

    Author:
        Lukas von Erdmannsdorff
    """
    df_core = df[df["Suite"].isin(["Core", SUITE_CORE])]
    if df_core.empty:
        df_core = df[df["Suite"].astype(str).str.lower() == "core"]
    if df_core.empty:
        df_core = df
    
    target_methods = ["HERA", "TOPSIS", "Borda", "Wilcoxon-Copeland", "Copeland"]
    methods = [m for m in target_methods if m in df_core["Method"].unique()]

    # Resolve candidate ranks if needed for Panel (E) Candidate-IQR
    if df_ranks is None or df_ranks.empty:
        cand_files = []
        for base in [output_path.parent, output_path.parent.parent, output_path.parent.parent / "CSVs", output_path.parent / "CSVs", DATA_DIR, DATA_DIR / "Simulation_Output", DATA_DIR / "Simulation_Output" / "CSVs"]:
            if base.exists():
                cand_files.extend(sorted(base.glob("candidate_ranks*.csv")))
                cand_files.extend(sorted(base.glob("CSVs/candidate_ranks*.csv")))
        for p_cand in cand_files:
            if p_cand.exists() and os.path.getsize(p_cand) > 0:
                try:
                    df_ranks = pd.read_csv(p_cand)
                    break
                except Exception:
                    pass

    df_iqr = None
    if df_ranks is not None and not df_ranks.empty and "PredictedRank" in df_ranks.columns:
        sub_ranks = df_ranks[df_ranks["Suite"].isin(["Core", SUITE_CORE])] if "Suite" in df_ranks.columns else df_ranks
        if sub_ranks.empty:
            sub_ranks = df_ranks
        group_cols = [c for c in ["ScenarioID", "Scenario", "Candidates", "Noise", "Method", "Candidate"] if c in sub_ranks.columns]
        if "Method" in group_cols and "Candidate" in group_cols:
            df_iqr = sub_ranks.groupby(group_cols)["PredictedRank"].apply(
                lambda x: float(np.percentile(x, 75) - np.percentile(x, 25))
            ).reset_index(name="CandidateIQR")
    
    scen_col = "Scenario" if "Scenario" in df_core.columns else "ScenarioID" if "ScenarioID" in df_core.columns else None
    n_scen = len(df_core[scen_col].unique()) if scen_col else 15

    cands = sorted(df_core["Candidates"].unique())
    cand_str = ", ".join(map(str, cands)) if cands else "10, 12, 14"
    noises = sorted(df_core["Noise"].unique())
    noise_str = ", ".join([f"{int(x) if x == int(x) else x}\\%" for x in noises]) if noises else r"2\%, 4\%, 6\%, 8\%, 10\%"

    fig, axes = plt.subplots(2, 3, figsize=(16, 9.6))
    fig.suptitle("Pooled Distribution across Core Conditions", fontweight="bold", fontsize=13.5, y=0.985)
    m0 = methods[0] if methods else "HERA"
    runs_per_method = len(df_core[df_core["Method"] == m0]) if m0 in df_core["Method"].values else len(df_core)
    fig.text(
        0.5, 0.950,
        rf"$\mathbf{{N \in \{{{cand_str}\}}}}$, $\boldsymbol{{\sigma}} \mathbf{{\in \{{{noise_str}\}}}}$, $\mathbf{{n = 50}}$, $\boldsymbol{{\rho = 0.0}}$; $\mathbf{{K = {runs_per_method}}}$ runs per method",
        ha="center", va="top", fontsize=11.5, fontweight="bold", color="#222222"
    )

    boxprops = dict(linewidth=1.2, edgecolor="#333333")
    medianprops = dict(linewidth=2.0, color="#111111")
    whiskerprops = dict(linewidth=1.1, color="#333333")
    capprops = dict(linewidth=1.1, color="#333333")
    flierprops = dict(marker=".", markersize=3, alpha=0.35, markeredgecolor="none")

    # (A) Top-Choice Recovery Rate (Binary Success Rate: Mean with 95% CI)
    ax = axes[0, 0]
    sns.barplot(
        data=df_core, x="Method", y="TopChoice", order=methods, palette=COLORS, hue="Method", legend=False,
        ax=ax, errorbar=("ci", 95), capsize=0.08, edgecolor="#333333", linewidth=1.1
    )
    for p in ax.patches:
        x = p.get_x() + p.get_width() / 2.
        h = p.get_height()
        if pd.notna(h) and h > 0:
            y_top = h
            for line in ax.lines:
                xs = line.get_xdata()
                if any(abs(xi - x) < 0.1 for xi in xs if not np.isnan(xi)):
                    ys = line.get_ydata()
                    if len(ys) > 0 and not np.all(np.isnan(ys)):
                        y_top = max(y_top, float(np.nanmax(ys)))
            ax.annotate(f"{h*100:.1f}%", (x, y_top),
                        ha='center', va='bottom', fontsize=9.0, xytext=(0, 4),
                        textcoords='offset points', fontweight='bold', color="#222222")
    ax.set_title("(A) Top-Choice Recovery Rate", fontweight="bold", loc="left")
    ax.set_xlabel("")
    ax.set_ylabel("Top-Choice Recovery (Mean [95% CI], %)")
    ax.set_ylim(0, 1.18)
    ax.yaxis.set_major_locator(MultipleLocator(0.20))
    ax.yaxis.set_major_formatter(PercentFormatter(1.0, decimals=0))
    ax.tick_params(labelleft=True)

    # (B) Complete-Rank Recovery Rate (Binary Success Rate: Mean with 95% CI)
    ax = axes[0, 1]
    sns.barplot(
        data=df_core, x="Method", y="CompleteRank", order=methods, palette=COLORS, hue="Method", legend=False,
        ax=ax, errorbar=("ci", 95), capsize=0.08, edgecolor="#333333", linewidth=1.1
    )
    for p in ax.patches:
        x = p.get_x() + p.get_width() / 2.
        h = p.get_height()
        if pd.notna(h) and h > 0:
            y_top = h
            for line in ax.lines:
                xs = line.get_xdata()
                if any(abs(xi - x) < 0.1 for xi in xs if not np.isnan(xi)):
                    ys = line.get_ydata()
                    if len(ys) > 0 and not np.all(np.isnan(ys)):
                        y_top = max(y_top, float(np.nanmax(ys)))
            ax.annotate(f"{h*100:.1f}%", (x, y_top),
                        ha='center', va='bottom', fontsize=9.0, xytext=(0, 4),
                        textcoords='offset points', fontweight='bold', color="#222222")
    ax.set_title("(B) Total Rank Recovery", fontweight="bold", loc="left")
    ax.set_xlabel("")
    ax.set_ylabel("Total Rank Recovery (Mean [95% CI], %)")
    cr_max = df_core.groupby("Method")["CompleteRank"].mean().max() if not df_core.empty else 0.3
    y_top_limit = max(0.48, min(1.15, cr_max * 1.65))
    ax.set_ylim(0, y_top_limit)
    ax.yaxis.set_major_locator(MultipleLocator(0.10 if y_top_limit <= 0.60 else 0.20))
    ax.yaxis.set_major_formatter(PercentFormatter(1.0, decimals=0))
    ax.tick_params(labelleft=True)

    # (C) Pairwise Inversion Rate (False Superiority Rate, %)
    ax = axes[0, 2]
    sns.boxplot(
        data=df_core, x="Method", y="FalseSuperiority", order=methods, palette=COLORS, hue="Method", legend=False,
        ax=ax, boxprops=boxprops, medianprops=medianprops, whiskerprops=whiskerprops,
        capprops=capprops, flierprops=flierprops
    )
    ax.set_title("(C) Pairwise Inversion Rate (False Superiority)", fontweight="bold", loc="left")
    ax.set_xlabel("")
    ax.set_ylabel("False Superiority (Median [IQR], %)")

    # Dynamic upper bound ensuring upper whiskers and outlier fliers are not clipped
    fsr_vals = df_core["FalseSuperiority"].dropna()
    fsr_max_val = float(fsr_vals.max()) if not fsr_vals.empty else 0.35
    upper_whiskers = []
    for _, g in df_core.groupby("Method"):
        v = g["FalseSuperiority"].dropna()
        if len(v) > 0:
            q1, q3 = np.percentile(v, 25), np.percentile(v, 75)
            upper_whiskers.append(min(float(v.max()), float(q3 + 1.5 * (q3 - q1))))
    max_whisker = max(upper_whiskers) if upper_whiskers else fsr_max_val
    y_top_target = max(0.40, max_whisker + 0.03, fsr_max_val + 0.02)
    y_top_fsr = float(np.ceil(y_top_target * 20.0) / 20.0)
    ax.set_ylim(-0.005, y_top_fsr)
    ax.yaxis.set_major_locator(MultipleLocator(0.05 if y_top_fsr <= 0.45 else 0.10))
    ax.yaxis.set_major_formatter(PercentFormatter(1.0, decimals=0))
    ax.tick_params(labelleft=True)

    # (D) Kendall's Tau
    ax = axes[1, 0]
    sns.boxplot(
        data=df_core, x="Method", y="KendallTau", order=methods, palette=COLORS, hue="Method", legend=False,
        ax=ax, boxprops=boxprops, medianprops=medianprops, whiskerprops=whiskerprops,
        capprops=capprops, flierprops=flierprops
    )
    ax.set_title(r"(D) Rank Correlation Fidelity (Kendall's $\boldsymbol{\tau}$)", fontweight="bold", loc="left")
    ax.set_xlabel("")
    ax.set_ylabel(r"Kendall's $\tau$ (Median [IQR])")
    ax.set_ylim(-0.05, 1.05)
    ax.yaxis.set_major_locator(MultipleLocator(0.20))
    ax.tick_params(labelleft=True)

    # (E) Empirical Rank Stability across Repeated Runs (Candidate-IQR)
    ax = axes[1, 1]
    if df_iqr is not None and not df_iqr.empty and "CandidateIQR" in df_iqr.columns:
        sns.boxplot(
            data=df_iqr, x="Method", y="CandidateIQR", order=methods, palette=COLORS, hue="Method", legend=False,
            ax=ax, boxprops=boxprops, medianprops=medianprops, whiskerprops=whiskerprops,
            capprops=capprops, flierprops=flierprops
        )
        ax.set_title(r"(E) Empirical Rank Stability (Candidate-IQR)", fontweight="bold", loc="left")
        ax.set_xlabel("")
        ax.set_ylabel("Rank Dispersion (IQR)")
        max_iqr = float(df_iqr["CandidateIQR"].max())
        ax.set_ylim(-0.2, max(4.0, max_iqr + 0.5))
        ax.yaxis.set_major_locator(MultipleLocator(1.0))
        ax.tick_params(labelleft=True)
    else:
        sns.boxplot(
            data=df_core, x="Method", y="RankDisplacement", order=methods, palette=COLORS, hue="Method", legend=False,
            ax=ax, boxprops=boxprops, medianprops=medianprops, whiskerprops=whiskerprops,
            capprops=capprops, flierprops=flierprops
        )
        ax.set_title(r"(E) Mean Rank Displacement ($\mathbf{|Rank - TrueRank|}$)", fontweight="bold", loc="left")
        ax.set_xlabel("")
        ax.set_ylabel("Rank Displacement (Median [IQR])")
        ax.set_ylim(-0.05, 3.2)
        ax.yaxis.set_major_locator(MultipleLocator(0.5))
        ax.tick_params(labelleft=True)

    # (F) Spearman's Rho
    ax = axes[1, 2]
    rho_col = "SpearmanRho" if "SpearmanRho" in df_core.columns else "KendallTau"
    sns.boxplot(
        data=df_core, x="Method", y=rho_col, order=methods, palette=COLORS, hue="Method", legend=False,
        ax=ax, boxprops=boxprops, medianprops=medianprops, whiskerprops=whiskerprops,
        capprops=capprops, flierprops=flierprops
    )
    ax.set_title(r"(F) Monotonic Alignment (Spearman's $\boldsymbol{\rho}$)", fontweight="bold", loc="left")
    ax.set_xlabel("")
    ax.set_ylabel(r"Spearman's $\rho$ (Median [IQR])")
    ax.set_ylim(-0.05, 1.05)
    ax.yaxis.set_major_locator(MultipleLocator(0.20))
    ax.tick_params(labelleft=True)

    plt.tight_layout(rect=[0, 0.03, 1, 0.91])
    save_figure_dual(fig, output_path, pdf_path)
    plt.close(fig)
    return fig


def plot_pooled_core_marginal_scales(
    df: pd.DataFrame, 
    output_path: Path, 
    pdf_path: Optional[Path] = None
) -> Optional[plt.Figure]:
    """
    Plots marginal pooled performance across Candidate Scales N and Cohort Sample Sizes n.

    Syntax:
        fig = plot_pooled_core_marginal_scales(df, output_path, pdf_path)

    Description:
        Evaluates marginal performance scaling across the two primary structural axes:
        (A1) Top-choice recovery vs. candidate scale N.
        (A2) Kendall's tau vs. candidate scale N.
        (B1) Top-choice recovery vs. cohort sample size n.
        (B2) Kendall's tau vs. cohort sample size n.

    Parameters:
        df (pd.DataFrame): Simulation results dataset.
        output_path (Path): Destination filesystem path for the output PNG graphic.
        pdf_path (Optional[Path]): Optional destination path for standalone vector PDF.

    Returns:
        Optional[plt.Figure]: Rendered matplotlib Figure instance (closed before return).

    Author:
        Lukas von Erdmannsdorff
    """
    df_core = df[df["Suite"] == "Core"]
    if df_core.empty:
        df_core = df

    fig, axes = plt.subplots(2, 2, figsize=(13, 9.6))
    fig.suptitle("Candidate Scalability and Sample Size Margins", fontweight="bold", fontsize=13.5, y=0.985)
    K = get_runs_per_point(df_core)
    n_noises = len(df_core["Noise"].dropna().unique()) if "Noise" in df_core.columns else 5
    K_margin = K * max(1, n_noises)
    fig.text(
        0.5, 0.950,
        rf"$\mathbf{{K = {K_margin}}}$ runs per point",
        ha="center", va="top", fontsize=11.5, fontweight="bold", color="#222222"
    )

    # (A1) TopChoice vs Candidates N
    ax = axes[0, 0]
    plot_method_lines(ax, df_core, "Candidates", "TopChoice")
    ax.set_title(r"(A1) Top-Choice Recovery vs. Candidate Scale ($\mathbf{N}$)", fontweight="bold", loc="left")
    ax.set_xlabel(r"Number of Candidates ($N$)")
    ax.set_ylabel("Top-Choice Recovery (Mean [95% CI], %)")
    ax.set_ylim(-0.02, 1.05)
    ax.yaxis.set_major_locator(MultipleLocator(0.20))
    ax.yaxis.set_major_formatter(PercentFormatter(1.0, decimals=0))
    cands = sorted(df_core["Candidates"].dropna().unique())
    if cands:
        ax.set_xticks(cands)
    if ax.get_legend():
        ax.get_legend().remove()

    # (A2) Kendall Tau vs Candidates N
    ax = axes[0, 1]
    plot_method_lines(ax, df_core, "Candidates", "KendallTau", estimator=np.median, errorbar=("pi", 50))
    ax.set_title(r"(A2) Rank Correlation vs. Candidate Scale ($\mathbf{N}$)", fontweight="bold", loc="left")
    ax.set_xlabel(r"Number of Candidates ($N$)")
    ax.set_ylabel(r"Kendall's $\tau$ (Median [IQR])")
    ax.set_ylim(-0.05, 1.05)
    ax.yaxis.set_major_locator(MultipleLocator(0.20))
    if cands:
        ax.set_xticks(cands)
    if ax.get_legend():
        ax.get_legend().remove()

    # (B1) TopChoice vs Sample Size n
    ax = axes[1, 0]
    plot_method_lines(ax, df_core, "SampleSize", "TopChoice")
    ax.set_title(r"(B1) Top-Choice Recovery vs. Sample Size ($\mathbf{n}$)", fontweight="bold", loc="left")
    ax.set_xlabel(r"Sample Size ($n$)")
    ax.set_ylabel("Top-Choice Recovery (Mean [95% CI], %)")
    ax.set_ylim(-0.02, 1.05)
    ax.yaxis.set_major_locator(MultipleLocator(0.20))
    ax.yaxis.set_major_formatter(PercentFormatter(1.0, decimals=0))
    samples = sorted(df_core["SampleSize"].dropna().unique())
    if samples:
        ax.set_xticks(samples)
    if ax.get_legend():
        ax.get_legend().remove()

    # (B2) Kendall Tau vs Sample Size n
    ax = axes[1, 1]
    plot_method_lines(ax, df_core, "SampleSize", "KendallTau", estimator=np.median, errorbar=("pi", 50))
    ax.set_title(r"(B2) Rank Correlation vs. Sample Size ($\mathbf{n}$)", fontweight="bold", loc="left")
    ax.set_xlabel(r"Sample Size ($n$)")
    ax.set_ylabel(r"Kendall's $\tau$ (Median [IQR])")
    ax.set_ylim(-0.05, 1.05)
    ax.yaxis.set_major_locator(MultipleLocator(0.20))
    if samples:
        ax.set_xticks(samples)
    if ax.get_legend():
        ax.get_legend().remove()

    add_external_legend(fig, axes[0, 0], ncol=4, y_pos=0.015)
    plt.tight_layout(rect=[0, 0.05, 1, 0.91])
    save_figure_dual(fig, output_path, pdf_path)
    plt.close(fig)
    return fig


def plot_pooled_compensatory_stress(
    df: pd.DataFrame, 
    output_path: Path, 
    pdf_path: Optional[Path] = None
) -> Optional[plt.Figure]:
    """
    Plots pooled non-compensatory diagnostic stress test evaluating flawed model protection.

    Syntax:
        fig = plot_pooled_compensatory_stress(df, output_path, pdf_path)

    Description:
        Displays pooled diagnostic metrics evaluating resistance against compensatory error:
        (A) Safe Model Recovery Rate (C1 as Rank 1)
        (B) Compensatory Failure Rate (C_flawed ranked in top half)
        (C) Mean Rank Displacement (|Rank - TrueRank|)

    Parameters:
        df (pd.DataFrame): Simulation results dataset.
        output_path (Path): Destination filesystem path for the output PNG graphic.
        pdf_path (Optional[Path]): Optional destination path for standalone vector PDF.

    Returns:
        Optional[plt.Figure]: Rendered matplotlib Figure instance (closed before return).

    Author:
        Lukas von Erdmannsdorff
    """
    sub = df[df["Suite"] == "Compensatory"]
    if sub.empty:
        sub = df
    
    target_methods = ["HERA", "TOPSIS", "Borda", "Wilcoxon-Copeland", "Copeland"]
    methods = [m for m in target_methods if m in sub["Method"].unique()]

    fig, axes = plt.subplots(1, 3, figsize=(15.5, 5.0))
    fig.suptitle("Non-Compensatory Diagnostic Stress Test", fontweight="bold", fontsize=13.0, y=0.985)
    K = get_runs_per_point(sub)
    fig.text(
        0.5, 0.925,
        rf"$\mathbf{{K = {K}}}$ runs per condition",
        ha="center", va="top", fontsize=11.5, fontweight="bold", color="#222222"
    )

    boxprops = dict(linewidth=1.2, edgecolor="#333333")
    medianprops = dict(linewidth=2.0, color="#111111")
    whiskerprops = dict(linewidth=1.1, color="#333333")
    capprops = dict(linewidth=1.1, color="#333333")
    flierprops = dict(marker=".", markersize=3, alpha=0.35, markeredgecolor="none")

    # (A) Safe Model Recovery Rate (Binary Rate: Mean with 95% CI)
    ax = axes[0]
    sns.barplot(
        data=sub, x="Method", y="TopChoice", order=methods, palette=COLORS, hue="Method", legend=False,
        ax=ax, errorbar=("ci", 95), capsize=0.08, edgecolor="#333333", linewidth=1.1
    )
    for p in ax.patches:
        x = p.get_x() + p.get_width() / 2.
        h = p.get_height()
        if pd.notna(h) and h > 0:
            y_top = h
            for line in ax.lines:
                xs = line.get_xdata()
                if any(abs(xi - x) < 0.1 for xi in xs if not np.isnan(xi)):
                    ys = line.get_ydata()
                    if len(ys) > 0 and not np.all(np.isnan(ys)):
                        y_top = max(y_top, float(np.nanmax(ys)))
            ax.annotate(f"{h*100:.1f}%", (x, y_top),
                        ha='center', va='bottom', fontsize=9.0, xytext=(0, 4),
                        textcoords='offset points', fontweight='bold', color="#222222")
    ax.set_title(r"(A) Safe Model Recovery ($\mathbf{C_1}$ as Rank 1)", fontweight="bold", loc="left")
    ax.set_xlabel("")
    ax.set_ylabel("Safe Model Recovery (Mean [95% CI], %)")
    ax.set_ylim(0, 1.18)
    ax.yaxis.set_major_locator(MultipleLocator(0.20))
    ax.yaxis.set_major_formatter(PercentFormatter(1.0, decimals=0))
    ax.tick_params(labelleft=True)

    # (B) Compensatory Failure Rate (Binary Rate: Mean with 95% CI)
    ax = axes[1]
    sns.barplot(
        data=sub, x="Method", y="CompensatoryError", order=methods, palette=COLORS, hue="Method", legend=False,
        ax=ax, errorbar=("ci", 95), capsize=0.08, edgecolor="#333333", linewidth=1.1
    )
    for p in ax.patches:
        x = p.get_x() + p.get_width() / 2.
        h = p.get_height()
        if pd.notna(h) and h > 0:
            y_top = h
            for line in ax.lines:
                xs = line.get_xdata()
                if any(abs(xi - x) < 0.1 for xi in xs if not np.isnan(xi)):
                    ys = line.get_ydata()
                    if len(ys) > 0 and not np.all(np.isnan(ys)):
                        y_top = max(y_top, float(np.nanmax(ys)))
            ax.annotate(f"{h*100:.1f}%", (x, y_top),
                        ha='center', va='bottom', fontsize=9.0, xytext=(0, 4),
                        textcoords='offset points', fontweight='bold', color="#222222")
    ax.set_title(r"(B) Compensatory Failure Rate ($\mathbf{C}_{\mathbf{flawed}}$ in Top Ranks)", fontweight="bold", loc="left")
    ax.set_xlabel("")
    ax.set_ylabel("Compensatory Failure (Mean [95% CI], %)")
    ax.set_ylim(0, 1.18)
    ax.yaxis.set_major_locator(MultipleLocator(0.20))
    ax.yaxis.set_major_formatter(PercentFormatter(1.0, decimals=0))
    ax.tick_params(labelleft=True)

    # (C) Mean Rank Displacement
    ax = axes[2]
    sns.boxplot(
        data=sub, x="Method", y="RankDisplacement", order=methods, palette=COLORS, hue="Method", legend=False,
        ax=ax, boxprops=boxprops, medianprops=medianprops, whiskerprops=whiskerprops,
        capprops=capprops, flierprops=flierprops
    )
    ax.set_title(r"(C) Mean Rank Displacement ($\mathbf{|Rank - TrueRank|}$)", fontweight="bold", loc="left")
    ax.set_xlabel("")
    ax.set_ylabel("Rank Displacement (Median [IQR])")
    ax.set_ylim(-0.05, 3.2)
    ax.yaxis.set_major_locator(MultipleLocator(0.5))
    ax.tick_params(labelleft=True)

    plt.tight_layout(rect=[0, 0.04, 1, 0.86])
    save_figure_dual(fig, output_path, pdf_path)
    plt.close(fig)
    return fig


def plot_pooled_sensitivity_summary(
    df: pd.DataFrame, 
    output_path: Path, 
    pdf_path: Optional[Path] = None
) -> Optional[plt.Figure]:
    """
    Plots pooled summary across Effect Magnitude and Metric Correlation sensitivity suites.

    Syntax:
        fig = plot_pooled_sensitivity_summary(df, output_path, pdf_path)

    Description:
        Displays side-by-side sensitivity comparisons across experimental suites:
        (A) Effect Magnitude Sensitivity (Small, Medium, Large separation).
        (B) Collinearity Robustness (Independent rho = 0.0 vs. Correlated rho = 0.5).

    Parameters:
        df (pd.DataFrame): Simulation results dataset.
        output_path (Path): Destination filesystem path for the output PNG graphic.
        pdf_path (Optional[Path]): Optional destination path for standalone vector PDF.

    Returns:
        Optional[plt.Figure]: Rendered matplotlib Figure instance (closed before return).

    Author:
        Lukas von Erdmannsdorff
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.0))
    fig.suptitle("Sensitivity Benchmarks: Effect Size and Collinearity", fontweight="bold", fontsize=13.0, y=0.985)
    K = get_runs_per_point(df)
    fig.text(
        0.5, 0.925,
        rf"$\mathbf{{K = {K}}}$ runs per point",
        ha="center", va="top", fontsize=11.5, fontweight="bold", color="#222222"
    )

    ref_N = DEFAULT_CANDIDATES if DEFAULT_CANDIDATES in df["Candidates"].values else (
        df["Candidates"].iloc[0] if "Candidates" in df.columns and not df.empty else 12
    )

    # Panel 1: Effect Magnitude (Small, Medium, Large)
    order_eff = ["Small", "Medium", "Large"]
    sub_eff = df[
        ((df["Suite"] == "EffectSensitivity") & (df["EffectMagnitude"].isin(order_eff))) |
        ((df["Suite"] == "Core") & (df["Candidates"] == ref_N) & (df["SampleSize"] == DEFAULT_SAMPLE_SIZE) & (df["Noise"] == DEFAULT_NOISE))
    ]
    if not sub_eff.empty:
        sns.barplot(
            data=sub_eff, x="EffectMagnitude", y="TopChoice", hue="Method",
            order=order_eff, palette=COLORS, ax=axes[0], errorbar=("ci", 95), capsize=0.08
        )
        axes[0].set_title(r"(A) Effect Size Sensitivity", fontweight="bold", loc="left")
        axes[0].set_xlabel("Effect Size (Median Cliff's d [Range])")
        axes[0].set_ylabel("Top-Choice Recovery (Mean [95% CI], %)")
        axes[0].set_ylim(0, 1.05)
        axes[0].yaxis.set_major_locator(MultipleLocator(0.20))
        axes[0].yaxis.set_major_formatter(PercentFormatter(1.0, decimals=0))
        if axes[0].get_legend():
            axes[0].get_legend().remove()

    # Panel 2: Correlation Sensitivity (rho = 0.0 vs rho = 0.5)
    sub_corr = df[
        ((df["Suite"] == "CorrelationSensitivity") & (df["Correlation"].isin([0.0, 0.5]))) |
        ((df["Suite"] == "Core") & (df["Candidates"] == ref_N) & (df["SampleSize"] == DEFAULT_SAMPLE_SIZE) & (df["Noise"] == DEFAULT_NOISE))
    ].copy()
    if not sub_corr.empty:
        sub_corr["CorrelationLabel"] = sub_corr["Correlation"].map({0.0: r"Independent ($\rho = 0.0$)", 0.5: r"Correlated ($\rho = 0.5$)"})
        order_corr = [r"Independent ($\rho = 0.0$)", r"Correlated ($\rho = 0.5$)"]
        sns.barplot(
            data=sub_corr, x="CorrelationLabel", y="TopChoice", hue="Method",
            order=order_corr, palette=COLORS, ax=axes[1], errorbar=("ci", 95), capsize=0.08
        )
        axes[1].set_title(r"(B) Collinearity Robustness ($\boldsymbol{\rho} = \mathbf{0.0}$ vs. $\boldsymbol{\rho} = \mathbf{0.5}$)", fontweight="bold", loc="left")
        axes[1].set_xlabel("Inter-Metric Correlation Condition")
        axes[1].set_ylabel("Top-Choice Recovery (Mean [95% CI], %)")
        axes[1].set_ylim(0, 1.05)
        axes[1].yaxis.set_major_locator(MultipleLocator(0.20))
        axes[1].yaxis.set_major_formatter(PercentFormatter(1.0, decimals=0))
        if axes[1].get_legend():
            axes[1].get_legend().remove()

    add_external_legend(fig, axes[0], ncol=4, y_pos=-0.06)
    plt.tight_layout(rect=[0, 0.06, 1, 0.86])
    save_figure_dual(fig, output_path, pdf_path)
    plt.close(fig)
    return fig


def plot_superordinate_table_figure(
    df_super: pd.DataFrame, 
    output_path: Optional[Path] = None, 
    pdf_path: Optional[Path] = None
) -> Optional[plt.Figure]:
    """
    Renders the executive superordinate summary table as a publication graphic.

    Syntax:
        fig = plot_superordinate_table_figure(df_super, output_path, pdf_path)

    Description:
        Formats and renders the executive summary table summarizing Grand Pooled
        benchmarks, structural margins, and sensitivity suites as a clean publication figure.

    Parameters:
        df_super (pd.DataFrame): Executive superordinate summary DataFrame.
        output_path (Optional[Path]): Optional destination filesystem path for PNG graphic.
        pdf_path (Optional[Path]): Optional destination path for standalone vector PDF.

    Returns:
        Optional[plt.Figure]: Rendered matplotlib Figure instance (closed before return).

    Author:
        Lukas von Erdmannsdorff
    """
    fig, ax = plt.subplots(figsize=(15.5, 10.5))
    ax.axis("off")
    fig.patch.set_facecolor("white")
    
    fig.suptitle("Superordinate Benchmark Summary", fontweight="bold", fontsize=13.5, y=0.985)
    fig.text(
        0.5, 0.955,
        "Grand Pooled Statistics, Structural Margins, and Sensitivity Suites (Median [IQR])",
        ha="center", va="top", fontsize=11.5, fontweight="bold", color="#222222"
    )

    if df_super.empty:
        ax.text(0.5, 0.5, "No Superordinate Summary Data Available", ha="center", va="center", fontsize=12)
        if output_path:
            save_figure_dual(fig, output_path, pdf_path)
        plt.close(fig)
        return fig

    # Prepare formatted display table
    display_rows = []
    headers = [
        "Superordinate Scenario", "Method", "N",
        "Top-Choice (%)", "Total Rank (%)", "Kendall τ",
        "False Sup. (%)", "Rank Displ. (Ranks)", "Comp. Error (%)"
    ]

    for _, row in df_super.iterrows():
        # TopChoice: Rate (%) or Median [IQR]
        if "TopChoice_Rate (%)" in row and pd.notna(row["TopChoice_Rate (%)"]):
            tc_str = f"{row['TopChoice_Rate (%)']:.1f}%"
        elif "TopChoice_Median (%)" in row and pd.notna(row["TopChoice_Median (%)"]):
            tc_str = f"{row['TopChoice_Median (%)']:.1f} [{row['TopChoice_IQR (%)']:.1f}]"
        else:
            tc_str = "N/A"

        # CompleteRank: Rate (%) or Median [IQR]
        if "CompleteRank_Rate (%)" in row and pd.notna(row["CompleteRank_Rate (%)"]):
            cr_str = f"{row['CompleteRank_Rate (%)']:.1f}%"
        elif "CompleteRank_Median (%)" in row and pd.notna(row["CompleteRank_Median (%)"]):
            cr_str = f"{row['CompleteRank_Median (%)']:.1f} [{row['CompleteRank_IQR (%)']:.1f}]"
        else:
            cr_str = "N/A"

        tau_str = f"{row['KendallTau_Median']:.3f} [{row['KendallTau_IQR']:.3f}]" if "KendallTau_Median" in row and pd.notna(row['KendallTau_Median']) else "N/A"
        fsr_str = f"{row['FalseSuperiority_Median (%)']:.2f} [{row['FalseSuperiority_IQR (%)']:.2f}]" if "FalseSuperiority_Median (%)" in row and pd.notna(row['FalseSuperiority_Median (%)']) else "N/A"
        disp_str = f"{row['RankDisplacement_Median (Ranks)']:.2f} [{row['RankDisplacement_IQR (Ranks)']:.2f}]" if "RankDisplacement_Median (Ranks)" in row and pd.notna(row['RankDisplacement_Median (Ranks)']) else (
            f"{row['ExpectedRegret_Median (%)']:.2f} [{row['ExpectedRegret_IQR (%)']:.2f}]" if "ExpectedRegret_Median (%)" in row and pd.notna(row['ExpectedRegret_Median (%)']) else "N/A"
        )

        # CompensatoryError: Rate (%) or Median [IQR]
        if "CompensatoryError_Rate (%)" in row and pd.notna(row["CompensatoryError_Rate (%)"]):
            cer_str = f"{row['CompensatoryError_Rate (%)']:.1f}%"
        elif "CompensatoryError_Median (%)" in row and pd.notna(row["CompensatoryError_Median (%)"]):
            cer_str = f"{row['CompensatoryError_Median (%)']:.1f} [{row['CompensatoryError_IQR (%)']:.1f}]"
        else:
            cer_str = "N/A"

        scen_name = row.get("Benchmark_Suite", row.get("Superordinate_Scenario", "Grand Pooled"))
        display_rows.append([
            str(scen_name),
            str(row["Method"]),
            f"{int(row['N_Trials'])}",
            tc_str,
            cr_str,
            tau_str,
            fsr_str,
            disp_str,
            cer_str
        ])

    table = ax.table(
        cellText=display_rows,
        colLabels=headers,
        cellLoc="center",
        loc="center"
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8.5)
    table.scale(1.0, 1.45)

    # Style header and alternating row colors
    for (r, c), cell in table.get_celld().items():
        if r == 0:
            cell.set_facecolor("#2B3E50")
            cell.set_text_props(color="white", fontweight="bold", fontsize=9.0)
            cell.set_height(0.042)
        else:
            if r % 2 == 1:
                cell.set_facecolor("#F8F9FA")
            else:
                cell.set_facecolor("white")
            cell.set_edgecolor("#E0E0E0")
            cell.set_linewidth(0.6)
            if c == 0:
                cell.set_text_props(ha="left")

    plt.tight_layout(rect=[0.02, 0.02, 0.98, 0.93])
    if output_path:
        save_figure_dual(fig, output_path, pdf_path)
    plt.close(fig)
    return fig


def generate_consolidated_pdf_reports(
    run_dir: Path,
    graphics_dir: Path,
    reports_dir: Path,
    timestamp: str
) -> Tuple[Path, Path]:
    """
    Compiles publication-grade multi-page PDF reports.

    Syntax:
        pdf_master, pdf_executive = generate_consolidated_pdf_reports(run_dir, graphics_dir, reports_dir, timestamp)

    Description:
        Assembles high-resolution publication PNG graphics into consolidated multi-page PDF reports:
          1. Global_Summary_<timestamp>.pdf:
             Master publication report bundling all primary figures in logical sequence.
             Mirrored to root run_dir for immediate top-level access.
          2. Executive_Summary_<timestamp>.pdf:
             Focused briefing report bundling key executive findings.

    Parameters:
        run_dir (Path): Root directory of the simulation run.
        graphics_dir (Path): Directory containing generated publication PNG figures.
        reports_dir (Path): Output directory for consolidated PDF reports.
        timestamp (str): Execution timestamp string (YYYYMMDD_HHMMSS).

    Returns:
        Tuple[Path, Path]: Tuple containing (master_pdf_path, executive_pdf_path).

    Author:
        Lukas von Erdmannsdorff
    """
    reports_dir.mkdir(parents=True, exist_ok=True)
    pdf_master = run_dir / f"Global_Summary_{timestamp}.pdf"
    pdf_executive = reports_dir / f"Executive_Summary_{timestamp}.pdf"

    # Remove any stale redundant duplicate in reports_dir if present
    stale_reports_pdf = reports_dir / f"Global_Summary_{timestamp}.pdf"
    if stale_reports_pdf.exists() and stale_reports_pdf != pdf_master:
        try:
            stale_reports_pdf.unlink()
        except Exception:
            pass

    def append_image_to_pdf(pdf: PdfPages, img_path: Path) -> None:
        """
        Loads a PNG figure, calculates aspect-ratio scaling, and appends it to a PDF report.
        """
        if not img_path.exists() or os.path.getsize(img_path) == 0:
            return
        try:
            img = mpimg.imread(str(img_path))
            h, w = img.shape[0], img.shape[1]
            aspect = w / max(1, h)
            fig_w = 12.0
            fig_h = max(5.0, min(14.0, fig_w / aspect))
            fig, ax = plt.subplots(figsize=(fig_w, fig_h))
            ax.imshow(img)
            ax.axis("off")
            plt.tight_layout(pad=0.15)
            pdf.savefig(fig, dpi=300, bbox_inches="tight")
            plt.close(fig)
        except Exception as e:
            print(f"    [Warning] Could not append image {img_path.name} to PDF: {e}")

    # Search for publication graphics in graphics_dir, run_dir, or CSVs/
    def find_figure_image(name_pattern: str) -> Optional[Path]:
        """
        Discovers the active figure image path across candidate directories by timestamp or glob pattern.
        """
        for search_folder in [graphics_dir, run_dir, run_dir / "CSVs"]:
            if search_folder.exists():
                direct = search_folder / f"{name_pattern}_{timestamp}.png"
                if direct.exists() and os.path.getsize(direct) > 0:
                    return direct
                cands = sorted(search_folder.glob(f"{name_pattern}_*.png"))
                if cands and os.path.getsize(cands[-1]) > 0:
                    return cands[-1]
                simple = search_folder / f"{name_pattern}.png"
                if simple.exists() and os.path.getsize(simple) > 0:
                    return simple
        return None

    # Ordered list of publication graphics for master report
    master_figure_sequence = [
        find_figure_image("method_parameters"),
        find_figure_image("scenario_configuration"),
        find_figure_image("candidate_table"),
        find_figure_image("candidate_profiles"),
        find_figure_image("validation_overview"),
        find_figure_image("pooled_core_distributions"),
        find_figure_image("top_choice_recovery_vs_noise"),
        find_figure_image("complete_rank_recovery_vs_noise"),
        find_figure_image("rank_displacement_vs_noise"),
        find_figure_image("false_superiority_vs_noise"),
        find_figure_image("rank_correlation_vs_sample"),
        find_figure_image("effect_magnitude_sensitivity"),
        find_figure_image("correlation_sensitivity"),
    ]
    comp_img = find_figure_image("compensatory_failure_vs_noise")
    if comp_img:
        master_figure_sequence.append(comp_img)
    master_figure_sequence.extend([
        find_figure_image("candidate_rank_stability"),
    ])

    # --- 1. Compile Master Global Summary Report directly in run_dir ---
    print(f" -> Compiling Master Global Summary PDF (Root): {pdf_master.name}...")
    with PdfPages(pdf_master) as pdf_m:
        for img_path in master_figure_sequence:
            if img_path and Path(img_path).exists():
                append_image_to_pdf(pdf_m, Path(img_path))

    # --- 2. Compile Focused Executive Summary Report in reports_dir ---
    executive_sequence = [
        find_figure_image("validation_overview"),
        find_figure_image("pooled_core_distributions"),
        find_figure_image("candidate_rank_stability"),
        find_figure_image("effect_magnitude_sensitivity"),
        find_figure_image("correlation_sensitivity")
    ]
    if comp_img:
        executive_sequence.insert(3, comp_img)
    print(f" -> Compiling Focused Executive Summary PDF: {pdf_executive.name}...")
    with PdfPages(pdf_executive) as pdf_e:
        for img_path in executive_sequence:
            if img_path and Path(img_path).exists():
                append_image_to_pdf(pdf_e, Path(img_path))

    return pdf_master, pdf_executive


def generate_global_summary_pdf(
    run_dir: Path,
    plots_dir: Path,
    df_raw: pd.DataFrame,
    timestamp: str,
    output_pdf: Optional[Path] = None
) -> Tuple[Path, Path]:
    """
    Backwards-compatibility wrapper routing to generate_consolidated_pdf_reports.

    Syntax:
        pdf_master, pdf_exec = generate_global_summary_pdf(run_dir, plots_dir, df_raw, timestamp, output_pdf)

    Description:
        Provides backwards-compatible interface for multi-page global summary PDF compilation,
        delegating directly to `generate_consolidated_pdf_reports`.

    Parameters:
        run_dir (Path): Root directory of the simulation run.
        plots_dir (Path): Directory containing generated graphics.
        df_raw (pd.DataFrame): Simulation raw results DataFrame.
        timestamp (str): Execution timestamp string.
        output_pdf (Optional[Path]): Optional explicit output PDF path (unused).

    Returns:
        Tuple[Path, Path]: Tuple of (master_pdf_path, executive_pdf_path).

    Author:
        Lukas von Erdmannsdorff
    """
    reports_dir = run_dir / "Reports" if (run_dir / "Reports").exists() else plots_dir
    return generate_consolidated_pdf_reports(
        run_dir=run_dir,
        graphics_dir=plots_dir,
        reports_dir=reports_dir,
        timestamp=timestamp
    )


def generate_all_plots(
    results_file: Optional[Path] = None,
    summary_file: Optional[Path] = None,
    graphics_dir: Optional[Path] = None,
    reports_dir: Optional[Path] = None,
    output_dir: Optional[Path] = None,
    timestamp: Optional[str] = None
) -> Tuple[Path, Path]:
    """
    Generates all publication figures and compiles consolidated PDF reports.

    Syntax:
        graphics_dir, reports_dir = generate_all_plots(results_file, summary_file, graphics_dir, reports_dir, output_dir, timestamp)

    Description:
        Coordinates the complete visualization and reporting pipeline:
          - Ingests raw simulation and rank distribution records.
          - Generates publication graphics (PNG, 300 DPI) with timestamp suffixes.
          - Compiles multi-page Master Global Summary and Executive Summary PDF reports.
          - Mirrors the Master PDF report to the root run directory for rapid access.

    Parameters:
        results_file (Optional[Path]): Path to simulation_results.csv.
        summary_file (Optional[Path]): Path to global_summary.csv.
        graphics_dir (Optional[Path]): Output directory for PNG figures.
        reports_dir (Optional[Path]): Output directory for consolidated PDF reports.
        output_dir (Optional[Path]): Legacy alias for graphics directory.
        timestamp (Optional[str]): Optional custom timestamp identifier.

    Returns:
        Tuple[Path, Path]: Tuple of (graphics_dir, reports_dir).

    Author:
        Lukas von Erdmannsdorff
    """
    print("[Visualization] Commencing generation of publication plots...")
    df_raw, _ = load_data(results_file, summary_file)
    res_path = results_file if results_file is not None else RESULTS_FILE
    if not res_path.exists():
        candidates = sorted(DATA_DIR.glob("Simulation_Output_*/CSVs/simulation_results*.csv"))
        if not candidates:
            candidates = sorted(DATA_DIR.glob("Simulation_Output_*/simulation_results*.csv"))
        if candidates:
            res_path = candidates[-1]

    # Resolve root run_dir
    if res_path.exists():
        run_dir = res_path.parent.parent if res_path.parent.name == "CSVs" else res_path.parent
    elif output_dir is not None:
        run_dir = output_dir.parent
    else:
        run_dir = DATA_DIR

    # Search for candidate ranks DataFrame
    df_ranks = None
    for s_dir in [res_path.parent, res_path.parent.parent / "CSVs", run_dir / "CSVs", run_dir]:
        if s_dir.exists():
            cands = sorted(s_dir.glob("candidate_ranks*.csv"))
            for c in reversed(cands):
                if c.exists() and os.path.getsize(c) > 0:
                    try:
                        df_ranks = pd.read_csv(c)
                        break
                    except Exception:
                        pass
            if df_ranks is not None:
                break

    # Resolve graphics_dir (default: run_dir / "Graphics")
    if graphics_dir is None:
        if output_dir is not None:
            graphics_dir = output_dir
        else:
            graphics_dir = run_dir / "Graphics"
    graphics_dir.mkdir(parents=True, exist_ok=True)

    # Resolve reports_dir (default: run_dir / "Reports")
    if reports_dir is None:
        reports_dir = run_dir / "Reports"
    reports_dir.mkdir(parents=True, exist_ok=True)

    # Detect or extract timestamp
    if timestamp is None:
        import re
        ts_match = re.search(r"(\d{8}_\d{6})", res_path.name)
        if not ts_match and hasattr(run_dir, "name"):
            ts_match = re.search(r"(\d{8}_\d{6})", run_dir.name)
        timestamp = ts_match.group(1) if ts_match else datetime.now().strftime("%Y%m%d_%H%M%S")

    # --- Generate Publication Figures (PNG only, 300 DPI, with timestamp suffix) ---
    print(" -> Generating: Validation Overview 8-Panel Benchmark...")
    f1_path = graphics_dir / f"validation_overview_{timestamp}.png"
    plot_master_overview(df_raw, f1_path, df_ranks=df_ranks)

    print(" -> Generating: Faceted Top-Choice Recovery Across Candidate Scales...")
    f2_path = graphics_dir / f"top_choice_recovery_vs_noise_{timestamp}.png"
    plot_faceted_top_choice(df_raw, f2_path)

    print(" -> Generating: Faceted Complete-Rank Recovery Across Candidate Scales...")
    f3_path = graphics_dir / f"complete_rank_recovery_vs_noise_{timestamp}.png"
    plot_faceted_complete_rank(df_raw, f3_path)

    print(" -> Generating: Faceted Rank Displacement Across Candidate Scales...")
    f4_path = graphics_dir / f"rank_displacement_vs_noise_{timestamp}.png"
    plot_faceted_rank_displacement(df_raw, f4_path)

    print(" -> Generating: Faceted False Superiority Rate Across Candidate Scales...")
    f5_path = graphics_dir / f"false_superiority_vs_noise_{timestamp}.png"
    plot_faceted_false_superiority(df_raw, f5_path)

    print(" -> Generating: Faceted Rank Correlation vs. Cohort Sample Size...")
    f6_path = graphics_dir / f"rank_correlation_vs_sample_{timestamp}.png"
    plot_faceted_rank_correlation(df_raw, f6_path)

    # Search for candidate_ranks*.csv for candidate-level plots
    cand_csv = None
    for s_dir in [run_dir / "CSVs", run_dir, res_path.parent, res_path.parent.parent / "CSVs"]:
        if s_dir.exists():
            cands = sorted(s_dir.glob("candidate_ranks*.csv"))
            for c in reversed(cands):
                if c.exists() and os.path.getsize(c) > 0:
                    cand_csv = c
                    break
            if cand_csv:
                break

    print(" -> Generating: Effect Size Sensitivity Benchmark...")
    f7_path = graphics_dir / f"effect_magnitude_sensitivity_{timestamp}.png"
    plot_effect_sensitivity(df_raw, f7_path, rank_data_path=cand_csv)

    print(" -> Generating: Inter-Metric Collinearity Sensitivity Benchmark...")
    f8_path = graphics_dir / f"correlation_sensitivity_{timestamp}.png"
    plot_correlation_sensitivity(df_raw, f8_path, rank_data_path=cand_csv)

    # Compensatory Diagnostic Stress Test (only if Compensatory suite present in dataset)
    if "Suite" in df_raw.columns and any(df_raw["Suite"].astype(str).str.contains("Compensatory", case=False)):
        print(" -> Generating: Compensatory Diagnostic Stress Test...")
        f9_path = graphics_dir / f"compensatory_failure_vs_noise_{timestamp}.png"
        plot_compensatory_failure(df_raw, f9_path)

    print(" -> Generating: Candidate-Level Rank Stability Across Scales...")
    f10_path = graphics_dir / f"candidate_rank_stability_{timestamp}.png"
    plot_candidate_rank_stability(f10_path, cand_csv)

    print(" -> Generating: Pooled Distribution...")
    f11_path = graphics_dir / f"pooled_core_distributions_{timestamp}.png"
    plot_pooled_core_distributions(df_raw, f11_path, df_ranks=df_ranks)

    # --- Compile Consolidated Multi-Page PDF Reports ---
    try:
        generate_consolidated_pdf_reports(
            run_dir=run_dir,
            graphics_dir=graphics_dir,
            reports_dir=reports_dir,
            timestamp=timestamp
        )
    except Exception as e:
        print(f"    [Warning] Could not compile consolidated summary PDF: {e}")

    print(f"[Success] All publication figures generated in: {graphics_dir}")
    print(f"[Success] Master global summary PDF generated in: {run_dir / f'Global_Summary_{timestamp}.pdf'}")
    print(f"[Success] Focused executive summary PDF generated in: {reports_dir}")
    return graphics_dir, reports_dir


def main() -> None:
    """
    CLI entry point for publication visualization and report generation.

    Syntax:
        main()

    Description:
        Parses command-line arguments, locates the target simulation run,
        and invokes generate_all_plots with custom or default output directories.

    Author:
        Lukas von Erdmannsdorff
    """
    import argparse
    import re

    parser = argparse.ArgumentParser(
        description="Generate publication plots and multi-page PDF reports from simulation results."
    )
    parser.add_argument(
        "target",
        nargs="?",
        default=None,
        help="Optional path to a simulation directory (e.g. Simulation_Output_*) or a simulation_results_*.csv file."
    )
    parser.add_argument(
        "-o", "--output", "--graphics-dir",
        dest="graphics_dir",
        default=None,
        help="Custom output directory for generated graphics (default: Graphics/ inside simulation folder)"
    )
    parser.add_argument(
        "-r", "--reports-dir",
        dest="reports_dir",
        default=None,
        help="Custom output directory for consolidated PDF reports (default: Reports/ inside simulation folder)"
    )
    args, _ = parser.parse_known_args()

    latest_res = None
    run_dir = None

    if args.target:
        target_path = Path(args.target).expanduser().resolve()
        if target_path.is_dir():
            run_dir = target_path
            res_candidates = sorted(run_dir.glob("CSVs/simulation_results*.csv")) or sorted(run_dir.glob("simulation_results*.csv"))
            if res_candidates:
                latest_res = res_candidates[-1]
        elif target_path.is_file():
            latest_res = target_path
            run_dir = latest_res.parent.parent if latest_res.parent.name == "CSVs" else latest_res.parent

    # If not provided or not resolved, detect latest Simulation_Output_* directory in CSVs/ or root
    if not latest_res:
        candidates = sorted(DATA_DIR.glob("Simulation_Output_*/CSVs/simulation_results*.csv"))
        if not candidates:
            candidates = sorted(DATA_DIR.glob("Simulation_Output_*/simulation_results*.csv"))
        if candidates:
            latest_res = candidates[-1]
            run_dir = latest_res.parent.parent if latest_res.parent.name == "CSVs" else latest_res.parent

    if latest_res and run_dir:
        # Detect timestamp
        ts_match = re.search(r"(\d{8}_\d{6})", latest_res.name)
        if not ts_match and hasattr(run_dir, "name"):
            ts_match = re.search(r"(\d{8}_\d{6})", run_dir.name)
        timestamp = ts_match.group(1) if ts_match else None

        latest_sum = None
        for s_dir in [run_dir / "CSVs", run_dir]:
            if s_dir.exists():
                if timestamp:
                    ts_sum = s_dir / f"global_summary_{timestamp}.csv"
                    if ts_sum.exists():
                        latest_sum = ts_sum
                        break
                sums = sorted(s_dir.glob("global_summary*.csv"))
                if sums:
                    latest_sum = sums[-1]
                    break

        graphics_dir = Path(args.graphics_dir).expanduser().resolve() if args.graphics_dir else run_dir / "Graphics"
        reports_dir = Path(args.reports_dir).expanduser().resolve() if args.reports_dir else run_dir / "Reports"

        generate_all_plots(
            results_file=latest_res,
            summary_file=latest_sum,
            graphics_dir=graphics_dir,
            reports_dir=reports_dir,
            timestamp=timestamp
        )
    else:
        generate_all_plots()


if __name__ == "__main__":
    main()

