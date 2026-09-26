"""
Reporting, Statistical Aggregation & FAIR Publication Graphics.

Handles incremental streaming CSV export, checkpointing, grand pooled summaries
(matching calc_pooled_csv.m), and 300 DPI publication configuration tables.

Author: Lukas von Erdmannsdorff
"""

import json
from pathlib import Path
from typing import Dict, List, Tuple, Any, Optional
import numpy as np
import pandas as pd
from scipy import stats

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from .config import (
    RESULT_COLUMNS, METRIC_COLS, STATISTICAL_EFFECT_COLS,
    EFFECT_CALIBRATION, HERA_CONFIG,
    DEFAULT_CANDIDATES, DEFAULT_SAMPLE_SIZE, DEFAULT_NOISE
)


def append_results_to_csv(
    csv_path: Path, 
    records: List[Dict[str, Any]], 
    expected_cols: Optional[List[str]] = None
) -> None:
    """
    Appends simulation records directly to CSV and enforces OS disk flush.

    Syntax:
        append_results_to_csv(csv_path, records, expected_cols=expected_cols)

    Description:
        Writes a list of trial result dictionaries to the designated CSV file in append mode.
        Filters columns against expected_cols to maintain strictly consistent tabular schemas.

    Parameters:
        csv_path (Path): Path to the target CSV file.
        records (List[Dict[str, Any]]): List of row dictionaries to write.
        expected_cols (Optional[List[str]]): Optional list of permitted column names.

    Returns:
        None

    Author:
        Lukas von Erdmannsdorff
    """
    if not records:
        return
    df_chunk = pd.DataFrame(records)
    cols = expected_cols if expected_cols is not None else (RESULT_COLUMNS if "TopChoice" in df_chunk.columns else list(df_chunk.columns))
    valid_cols = [c for c in cols if c in df_chunk.columns]
    df_chunk = df_chunk[valid_cols]
    df_chunk.to_csv(csv_path, mode="a", header=False, index=False)


def compute_aggregate_summary(df_results: pd.DataFrame) -> pd.DataFrame:
    """
    Calculates pooled descriptive statistics matching HERA's +analysis standard.

    Syntax:
        df_summary = compute_aggregate_summary(df_results)

    Description:
        Aggregates raw Monte Carlo trial results across experimental conditions
        (ScenarioID, Suite, Candidates, SampleSize, Noise, EffectMagnitude, Delta, Correlation, Method).
        Computes Mean, Std, Median, IQR, Q1, Q3, 95% Bootstrap/Empirical CIs, Min, and Max
        for every evaluated metric.

    Parameters:
        df_results (pd.DataFrame): Raw trial-level simulation results DataFrame.

    Returns:
        pd.DataFrame: Aggregated summary statistics DataFrame.

    Author:
        Lukas von Erdmannsdorff
    """
    if df_results.empty:
        return pd.DataFrame()

    raw_group_cols = ["ScenarioID", "Suite", "Candidates", "SampleSize", "Noise", "EffectMagnitude", "Delta", "Calibrated_Median_Cliffs_d", "Correlation", "Method"]
    group_cols = [c for c in raw_group_cols if c in df_results.columns]
    summary_rows = []

    active_metric_cols = [c for c in (METRIC_COLS + STATISTICAL_EFFECT_COLS) if c in df_results.columns]

    for keys, group in df_results.groupby(group_cols):
        row = dict(zip(group_cols, keys))
        n_trials = len(group)
        row["N_Trials"] = n_trials

        for m in active_metric_cols:
            vals = group[m].dropna().values
            if len(vals) == 0:
                continue
            row[f"{m}_Mean"] = float(np.mean(vals))
            row[f"{m}_Std"] = float(np.std(vals, ddof=1)) if len(vals) > 1 else 0.0
            row[f"{m}_Median"] = float(np.median(vals))
            row[f"{m}_IQR"] = float(stats.iqr(vals))
            row[f"{m}_Q1"] = float(np.percentile(vals, 25))
            row[f"{m}_Q3"] = float(np.percentile(vals, 75))
            row[f"{m}_CI95_Lower"] = float(np.percentile(vals, 2.5))
            row[f"{m}_CI95_Upper"] = float(np.percentile(vals, 97.5))
            row[f"{m}_Min"] = float(np.min(vals))
            row[f"{m}_Max"] = float(np.max(vals))

        row["Failure_Rate_Percent"] = 0.0
        summary_rows.append(row)

    return pd.DataFrame(summary_rows)


def compute_benchmark_summary(df_results: pd.DataFrame) -> pd.DataFrame:
    """
    Calculates an executive benchmark summary table across primary benchmark dimensions.

    Syntax:
        df_bench = compute_benchmark_summary(df_results)

    Description:
        Provides an immediately readable, presentation-grade tabular overview where all
        overarching benchmark suites and MCDM algorithms are directly comparable side-by-side
        without unwieldy horizontal column overflow.
        
        Benchmark Suites Evaluated:
          1. Core Benchmark (Grand Pooled across standard conditions)
          2. Core Benchmark Scale Margins (Stratified by N in {10, 12, 14})
          3. Cohort Sample Size Margins (Stratified by n in {25, 50, 100})
          4. Effect Magnitude Sensitivity (Small, Medium, Large Effect sizes)
          5. Correlation Sensitivity (Independent vs. Correlated)
          6. Non-Compensatory Diagnostic Stress Test (Pooled across conditions)

    Parameters:
        df_results (pd.DataFrame): Raw or aggregated simulation DataFrame containing Monte Carlo iterations.

    Returns:
        pd.DataFrame: Structured DataFrame with human-readable scenario names, sample sizes, and
            non-parametric performance metrics (Median and IQR).

    Author:
        Lukas von Erdmannsdorff
    """
    if df_results.empty:
        return pd.DataFrame()

    target_methods = ["HERA", "TOPSIS", "Borda", "Wilcoxon-Copeland", "Copeland"]
    available_methods = [m for m in target_methods if m in df_results["Method"].unique()]
    if not available_methods:
        available_methods = list(df_results["Method"].unique())

    slices = []
    
    # 1. Core Grand Pooled
    slices.append((
        "Core Benchmark (Grand Pooled)",
        "Core",
        "All Standard Conditions (N: 8-14, n: 25-100, noise: 2-10%)",
        lambda df: df["Suite"] == "Core"
    ))
    
    # 2. Core Candidate Scale Margins
    if "Candidates" in df_results.columns:
        core_cands = sorted(df_results[df_results["Suite"] == "Core"]["Candidates"].dropna().unique())
        for N in core_cands:
            slices.append((
                f"Core Scale (N = {int(N)} Candidates)",
                "Core",
                f"Candidate Set N = {int(N)}",
                lambda df, n_val=N: (df["Suite"] == "Core") & (df["Candidates"] == n_val)
            ))
            
    # 3. Core Cohort Sample Size Margins (Evaluated across sample sizes n in {25, 50, 100})
    for n_val in [25, 50, 100]:
        slices.append((
            f"Cohort Sample Size (n = {n_val})",
            "Sample Size Sensitivity",
            f"Cohort Size n = {n_val} (at N = {DEFAULT_CANDIDATES})",
            lambda df, sv=n_val: (
                (df["Suite"].isin(["SampleSizeSensitivity", "Sample Size Sensitivity"]) & (df["SampleSize"] == sv)) |
                ((sv == DEFAULT_SAMPLE_SIZE) & (df["Suite"] == "Core") & (df["Candidates"] == DEFAULT_CANDIDATES) & (df["Noise"] == DEFAULT_NOISE))
            )
        ))
            
    # 4. Effect Magnitude Sensitivity Suite (Small, Medium, Large)
    if "EffectMagnitude" in df_results.columns:
        for eff, label in [("Small", "Small Effect (Delta = 2.5%, d ~ 0.25)"),
                           ("Medium", "Medium Effect (Delta = 5.0%, d ~ 0.50)"),
                           ("Large", "Large Effect (Delta = 8.0%, d ~ 0.80)")]:
            slices.append((
                f"Effect Sensitivity ({eff})",
                "Effect Magnitude Sensitivity",
                label,
                lambda df, e_val=eff: (
                    (df["Suite"].isin(["EffectSensitivity", "Effect Magnitude Sensitivity"]) & (df["EffectMagnitude"] == e_val)) |
                    ((e_val == "Medium") & (df["Suite"] == "Core") & (df["Candidates"] == DEFAULT_CANDIDATES) & (df["SampleSize"] == DEFAULT_SAMPLE_SIZE) & (df["Noise"] == DEFAULT_NOISE))
                )
            ))
            
    # 5. Correlation Sensitivity Suite (Independent vs Correlated)
    if "Correlation" in df_results.columns:
        for corr, label in [(0.0, "Independent Criteria (rho = 0.0)"),
                            (0.5, "Correlated Criteria (rho = 0.5)")]:
            slices.append((
                f"Correlation Sensitivity (rho = {corr})",
                "Correlation Sensitivity",
                label,
                lambda df, c_val=corr: (
                    (df["Suite"].isin(["CorrelationSensitivity", "Correlation Sensitivity"]) & (df["Correlation"] == c_val)) |
                    ((c_val == 0.0) & (df["Suite"] == "Core") & (df["Candidates"] == DEFAULT_CANDIDATES) & (df["SampleSize"] == DEFAULT_SAMPLE_SIZE) & (df["Noise"] == DEFAULT_NOISE))
                )
            ))
            
    # 6. Non-Compensatory Stress Test (only if Compensatory suite exists in dataset)
    if "Compensatory" in df_results["Suite"].values:
        slices.append((
            "Non-Compensatory Stress Test (Pooled)",
            "Compensatory",
            "Pathological Shortcut Model Trap (CN in top half)",
            lambda df: (df["Suite"] == "Compensatory")
        ))

    rows = []
    for scen_title, suite, cond_label, mask_fn in slices:
        sub = df_results[mask_fn(df_results)]
        if sub.empty:
            continue
            
        for method in available_methods:
            m_sub = sub[sub["Method"] == method]
            if m_sub.empty:
                continue
                
            n_trials = len(m_sub)
            
            def get_stat(col: str, multiplier: float = 1.0, decimals: int = 1) -> Tuple[float, float]:
                if col not in m_sub.columns:
                    return np.nan, np.nan
                vals = m_sub[col].dropna().values * multiplier
                if len(vals) == 0:
                    return np.nan, np.nan
                med = round(float(np.median(vals)), decimals)
                iqr_val = round(float(stats.iqr(vals)), decimals)
                return med, iqr_val

            def get_rate(col: str) -> float:
                if col not in m_sub.columns or m_sub[col].dropna().empty:
                    return np.nan
                return round(float(np.mean(m_sub[col].dropna().values)) * 100.0, 1)

            top_rate = get_rate("TopChoice")
            comp_rate = get_rate("CompleteRank")
            tau_med, tau_iqr = get_stat("KendallTau", multiplier=1.0, decimals=3)
            rho_med, rho_iqr = get_stat("SpearmanRho", multiplier=1.0, decimals=3)
            fsr_med, fsr_iqr = get_stat("FalseSuperiority", multiplier=100.0, decimals=2)
            disp_med, disp_iqr = get_stat("RankDisplacement", multiplier=1.0, decimals=2)
            regret_med, regret_iqr = get_stat("TopChoiceRankRegret", multiplier=1.0, decimals=2)

            row = {
                "Benchmark_Suite": scen_title,
                "Suite": suite,
                "Stratification": cond_label,
                "Method": method,
                "N_Trials": n_trials,
                "TopChoice_Rate (%)": top_rate,
                "CompleteRank_Rate (%)": comp_rate,
                "KendallTau_Median": tau_med,
                "KendallTau_IQR": tau_iqr,
                "SpearmanRho_Median": rho_med,
                "SpearmanRho_IQR": rho_iqr,
                "FalseSuperiority_Median (%)": fsr_med,
                "FalseSuperiority_IQR (%)": fsr_iqr,
                "RankDisplacement_Median (Ranks)": disp_med,
                "RankDisplacement_IQR (Ranks)": disp_iqr,
                "TopChoiceRankRegret_Median (Ranks)": regret_med,
                "TopChoiceRankRegret_IQR (Ranks)": regret_iqr
            }
            if "CompensatoryError" in m_sub.columns and m_sub["CompensatoryError"].dropna().any():
                row["CompensatoryError_Rate (%)"] = get_rate("CompensatoryError")
            rows.append(row)

    return pd.DataFrame(rows)


# Backwards compatibility alias
compute_superordinate_summary = compute_benchmark_summary


def compute_detailed_pooled_statistics(df_results: pd.DataFrame) -> pd.DataFrame:
    """
    Calculates granular non-parametric pooled statistics matching calc_pooled_csv.m schema.

    Syntax:
        df_detailed = compute_detailed_pooled_statistics(df_results)

    Description:
        Produces a standardized long-format dataset with complete descriptive properties:
        Scenario;Metric;Method;N_Trials;Median;IQR;Q1;Q3;CI95_Lower;CI95_Upper;Min;Max;Failure_Rate_Percent.

    Parameters:
        df_results (pd.DataFrame): Simulation results dataset.

    Returns:
        pd.DataFrame: Long-format detailed statistics DataFrame.

    Author:
        Lukas von Erdmannsdorff
    """
    if df_results.empty:
        return pd.DataFrame()

    active_metric_cols = [c for c in (METRIC_COLS + STATISTICAL_EFFECT_COLS) if c in df_results.columns]
    
    target_methods = ["HERA", "TOPSIS", "Borda", "Wilcoxon-Copeland", "Copeland"]
    available_methods = [m for m in target_methods if m in df_results["Method"].unique()]
    if not available_methods:
        available_methods = list(df_results["Method"].unique())

    # Superordinate grouping conditions
    conditions = [
        ("Core_Grand_Pooled", lambda df: df["Suite"] == "Core"),
        ("Compensatory_Stress_Pooled", lambda df: df["Suite"] == "Compensatory")
    ]
    if "Candidates" in df_results.columns:
        for N in sorted(df_results[df_results["Suite"] == "Core"]["Candidates"].dropna().unique()):
            conditions.append((f"Core_N_{int(N)}", lambda df, n_val=N: (df["Suite"] == "Core") & (df["Candidates"] == n_val)))
    if "SampleSize" in df_results.columns:
        for n_val in [25, 50, 100]:
            conditions.append((f"Cohort_n_{n_val}", lambda df, sv=n_val: (
                (df["Suite"].isin(["SampleSizeSensitivity", "Sample Size Sensitivity"]) & (df["SampleSize"] == sv)) |
                ((sv == DEFAULT_SAMPLE_SIZE) & (df["Suite"] == "Core") & (df["Candidates"] == DEFAULT_CANDIDATES) & (df["Noise"] == DEFAULT_NOISE))
            )))
    if "EffectMagnitude" in df_results.columns:
        for eff in ["Small", "Medium", "Large"]:
            conditions.append((f"Effect_{eff}", lambda df, e_val=eff: (
                (df["Suite"].isin(["EffectSensitivity", "Effect Magnitude Sensitivity"]) & (df["EffectMagnitude"] == e_val)) |
                ((e_val == "Medium") & (df["Suite"] == "Core") & (df["Candidates"] == DEFAULT_CANDIDATES) & (df["SampleSize"] == DEFAULT_SAMPLE_SIZE) & (df["Noise"] == DEFAULT_NOISE))
            )))
    if "Correlation" in df_results.columns:
        for corr in [0.0, 0.5]:
            c_tag = "Independent" if corr == 0.0 else "Correlated"
            conditions.append((f"Correlation_{c_tag}", lambda df, c_val=corr: (
                (df["Suite"].isin(["CorrelationSensitivity", "Correlation Sensitivity"]) & (df["Correlation"] == c_val)) |
                ((c_val == 0.0) & (df["Suite"] == "Core") & (df["Candidates"] == DEFAULT_CANDIDATES) & (df["SampleSize"] == DEFAULT_SAMPLE_SIZE) & (df["Noise"] == DEFAULT_NOISE))
            )))

    rows = []
    for scen_name, mask_fn in conditions:
        sub = df_results[mask_fn(df_results)]
        if sub.empty:
            continue
            
        for method in available_methods:
            m_sub = sub[sub["Method"] == method]
            if m_sub.empty:
                continue
                
            n_trials = len(m_sub)
            for m in active_metric_cols:
                vals = m_sub[m].dropna().values
                if len(vals) == 0:
                    continue
                rows.append({
                    "Scenario": scen_name,
                    "Metric": m,
                    "Method": method,
                    "N_Trials": n_trials,
                    "Mean": round(float(np.mean(vals)), 4),
                    "Std": round(float(np.std(vals, ddof=1)), 4) if len(vals) > 1 else 0.0,
                    "Median": round(float(np.median(vals)), 4),
                    "IQR": round(float(stats.iqr(vals)), 4),
                    "Q1": round(float(np.percentile(vals, 25)), 4),
                    "Q3": round(float(np.percentile(vals, 75)), 4),
                    "CI95_Lower": round(float(np.percentile(vals, 2.5)), 4),
                    "CI95_Upper": round(float(np.percentile(vals, 97.5)), 4),
                    "Min": round(float(np.min(vals)), 4),
                    "Max": round(float(np.max(vals)), 4),
                    "Failure_Rate_Percent": 0.0
                })

    return pd.DataFrame(rows)


def compute_pooled_summary(df_results: pd.DataFrame) -> pd.DataFrame:
    """
    Calculates grand pooled convergence statistics matching calc_pooled_csv.m.

    Syntax:
        df_pooled = compute_pooled_summary(df_results)

    Description:
        Aggregates metrics both overall (All_Pooled) and per experimental suite (Core,
        Sensitivity suites), calculating full descriptive metrics (Mean, Std, Median, IQR,
        Q1, Q3, CIs, Min, Max).

    Parameters:
        df_results (pd.DataFrame): Simulation results dataset.

    Returns:
        pd.DataFrame: Pooled summary statistics DataFrame.

    Author:
        Lukas von Erdmannsdorff
    """
    if df_results.empty:
        return pd.DataFrame()

    active_metric_cols = [c for c in (METRIC_COLS + STATISTICAL_EFFECT_COLS) if c in df_results.columns]
    summary_rows = []

    # 1. Overall Pooling by Method
    for method, group in df_results.groupby("Method"):
        row = {"Suite": "All_Pooled", "Method": method, "N_Trials": len(group)}
        for m in active_metric_cols:
            vals = group[m].dropna().values
            if len(vals) == 0:
                continue
            row[f"{m}_Mean"] = float(np.mean(vals))
            row[f"{m}_Std"] = float(np.std(vals, ddof=1)) if len(vals) > 1 else 0.0
            row[f"{m}_Median"] = float(np.median(vals))
            row[f"{m}_IQR"] = float(stats.iqr(vals))
            row[f"{m}_Q1"] = float(np.percentile(vals, 25))
            row[f"{m}_Q3"] = float(np.percentile(vals, 75))
            row[f"{m}_CI95_Lower"] = float(np.percentile(vals, 2.5))
            row[f"{m}_CI95_Upper"] = float(np.percentile(vals, 97.5))
            row[f"{m}_Min"] = float(np.min(vals))
            row[f"{m}_Max"] = float(np.max(vals))
        row["Failure_Rate_Percent"] = 0.0
        summary_rows.append(row)

    # 2. Granular Pooling by Suite and Method
    for (suite, method), group in df_results.groupby(["Suite", "Method"]):
        row = {"Suite": suite, "Method": method, "N_Trials": len(group)}
        for m in active_metric_cols:
            vals = group[m].dropna().values
            if len(vals) == 0:
                continue
            row[f"{m}_Mean"] = float(np.mean(vals))
            row[f"{m}_Std"] = float(np.std(vals, ddof=1)) if len(vals) > 1 else 0.0
            row[f"{m}_Median"] = float(np.median(vals))
            row[f"{m}_IQR"] = float(stats.iqr(vals))
            row[f"{m}_Q1"] = float(np.percentile(vals, 25))
            row[f"{m}_Q3"] = float(np.percentile(vals, 75))
            row[f"{m}_CI95_Lower"] = float(np.percentile(vals, 2.5))
            row[f"{m}_CI95_Upper"] = float(np.percentile(vals, 97.5))
            row[f"{m}_Min"] = float(np.min(vals))
            row[f"{m}_Max"] = float(np.max(vals))
        row["Failure_Rate_Percent"] = 0.0
        summary_rows.append(row)

    return pd.DataFrame(summary_rows)


def export_summary_csvs(
    df_results: pd.DataFrame, 
    out_dir: Path, 
    timestamp: str,
    root_dir: Optional[Path] = None
) -> Dict[str, Path]:
    """
    Generates and exports executive benchmark suite summary and detailed pooled CSVs.

    Syntax:
        paths = export_summary_csvs(df_results, out_dir, timestamp, root_dir=root_dir)

    Description:
        Exports standardized tabular results:
          1. benchmark_summary_<timestamp>.csv: Executive table with directly readable key metrics
             (stored directly in the main output folder root_dir for immediate top-level inspection).
          2. pooled_results_detailed_<timestamp>.csv: Comprehensive non-parametric breakdown (CSVs/).
          3. pooled_results_<timestamp>.csv: Wide-format grand pooled dataset (CSVs/).

    Parameters:
        df_results (pd.DataFrame): Simulation results DataFrame or path.
        out_dir (Path): Output directory for CSV files (CSVs/).
        timestamp (str): Execution timestamp string.
        root_dir (Optional[Path]): Optional top-level run directory for main output placement.

    Returns:
        Dict[str, Path]: Mapping of export keys to output file Paths.

    Author:
        Lukas von Erdmannsdorff
    """
    out_dir = Path(out_dir)
    if root_dir is not None:
        root_dir = Path(root_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    if root_dir is not None:
        root_dir.mkdir(parents=True, exist_ok=True)
    paths = {}

    if isinstance(df_results, (str, Path)):
        df_results = pd.read_csv(df_results)

    # 1. Executive Benchmark Suite Summary CSV (Stored solely in main output root_dir if specified)
    df_bench = compute_benchmark_summary(df_results)
    if not df_bench.empty:
        target_bench_dir = root_dir if root_dir is not None else out_dir
        bench_csv = target_bench_dir / f"benchmark_summary_{timestamp}.csv"
        df_bench.to_csv(bench_csv, index=False)
        paths["benchmark_summary"] = bench_csv
        paths["superordinate_summary"] = bench_csv

        # Remove stale duplicate in out_dir if it was previously written there
        if root_dir is not None and root_dir != out_dir:
            stale_csv = out_dir / f"benchmark_summary_{timestamp}.csv"
            if stale_csv.exists() and stale_csv != bench_csv:
                try:
                    stale_csv.unlink()
                except Exception:
                    pass

    # 2. Detailed Pooled Statistics CSV (matching calc_pooled_csv.m)
    df_detailed = compute_detailed_pooled_statistics(df_results)
    if not df_detailed.empty:
        detailed_csv = out_dir / f"pooled_results_detailed_{timestamp}.csv"
        df_detailed.to_csv(detailed_csv, index=False)
        paths["pooled_detailed"] = detailed_csv

    # 3. Grand Pooled Summary CSV
    df_pooled = compute_pooled_summary(df_results)
    if not df_pooled.empty:
        pooled_csv = out_dir / f"pooled_results_{timestamp}.csv"
        df_pooled.to_csv(pooled_csv, index=False)
        paths["pooled_summary"] = pooled_csv

    return paths


# Backwards compatibility alias
export_superordinate_csvs = export_summary_csvs


def update_global_summary_file(summary_csv: Path, df_all_results: pd.DataFrame) -> None:
    """
    Recomputes and overwrites global_summary.csv atomically with all available results.

    Syntax:
        update_global_summary_file(summary_csv, df_all_results)

    Description:
        Aggregates current in-memory results and writes to a temporary file before atomically
        replacing summary_csv, guaranteeing robust file integrity even during parallel execution.

    Parameters:
        summary_csv (Path): Destination CSV path.
        df_all_results (pd.DataFrame): Complete cumulative results DataFrame.

    Returns:
        None

    Author:
        Lukas von Erdmannsdorff
    """
    if df_all_results.empty:
        return
    df_summary = compute_aggregate_summary(df_all_results)
    temp_summary = summary_csv.with_suffix(".tmp")
    df_summary.to_csv(temp_summary, index=False)
    temp_summary.replace(summary_csv)


def save_checkpoint(checkpoint_file: Path, state_data: Dict[str, Any]) -> None:
    """
    Saves execution checkpoint to JSON atomically.

    Syntax:
        save_checkpoint(checkpoint_file, state_data)

    Description:
        Writes the current run progress state to a temporary file and atomically renames
        it to checkpoint_file, preventing corruption on unexpected interrupts.

    Parameters:
        checkpoint_file (Path): Destination checkpoint JSON file path.
        state_data (Dict[str, Any]): State metadata dictionary.

    Returns:
        None

    Author:
        Lukas von Erdmannsdorff
    """
    temp_ckpt = checkpoint_file.with_suffix(".tmp")
    with open(temp_ckpt, "w", encoding="utf-8") as f:
        json.dump(state_data, f, indent=2)
    temp_ckpt.replace(checkpoint_file)


def save_scenario_configuration(
    scenarios: List[Dict[str, Any]], 
    out_dir: Path, 
    timestamp: str,
    graphics_dir: Optional[Path] = None
) -> Tuple[Path, Path]:
    """
    Saves experimental scenario definitions to CSV and renders a 300 DPI publication graphic.

    Syntax:
        csv_path, png_path = save_scenario_configuration(scenarios, out_dir, timestamp, graphics_dir=graphics_dir)

    Description:
        Exports the complete experimental grid to CSV and renders an executive 300 DPI
        publication table summarizing each benchmark suite's parameter sweeps.

    Parameters:
        scenarios (List[Dict[str, Any]]): List of scenario configuration dictionaries.
        out_dir (Path): Output directory for CSV table.
        timestamp (str): Execution timestamp string.
        graphics_dir (Optional[Path]): Directory for rendered publication graphic.

    Returns:
        Tuple[Path, Path]: Tuple of (csv_path, png_path).

    Author:
        Lukas von Erdmannsdorff
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / f"scenarios_{timestamp}.csv"
    df_sc = pd.DataFrame(scenarios)
    df_sc.to_csv(csv_path, index=False)

    # Generate Publication Figure of Scenario Grid in graphics_dir
    target_graphics_dir = graphics_dir if graphics_dir is not None else out_dir
    target_graphics_dir.mkdir(parents=True, exist_ok=True)
    png_path = target_graphics_dir / f"scenario_configuration_{timestamp}.png"

    suite_summary = []
    for suite_name, group in df_sc.groupby("Suite", sort=False):
        cands = sorted(group["Candidates"].unique().tolist())
        samples = sorted(group["SampleSize"].unique().tolist())
        noises = sorted(group["Noise"].unique().tolist())
        deltas = sorted(group["Delta"].unique().tolist())
        rhos = sorted(group["Correlation"].unique().tolist())

        eff_parts = []
        cliffs_parts = []
        for d in deltas:
            calib = None
            for c_name, c_info in EFFECT_CALIBRATION.items():
                if abs(c_info["delta"] - d) < 1e-4:
                    calib = c_info
                    break
            if calib:
                eff_parts.append(f"{d:g}% ({calib['label'].split()[0]})")
                cliffs_parts.append(f"d ≈ {calib['median_cliffs_d']:.2f} [{calib['cliffs_d_range']}]")
            else:
                eff_parts.append(f"{d:g}%")
                cliffs_parts.append("N/A")

        # Format benchmark suite label with proper word spacing
        display_suite = str(suite_name).replace("Sensitivity", " Sensitivity").replace("  ", " ").strip()
        suite_summary.append({
            "Benchmark Suite": display_suite,
            "N Conditions": len(group),
            "Candidates (N)": ", ".join(map(str, cands)),
            "Sample Size (n)": ", ".join(map(str, samples)),
            r"Noise ($\sigma$, %)": ", ".join([f"{x:g}%" for x in noises]),
            r"Effect Size (Cliff's $d$ [Range])": ", ".join(cliffs_parts),
            r"Correlation ($\rho$)": ", ".join([f"{x:g}" for x in rhos])
        })
    df_summary = pd.DataFrame(suite_summary)

    fig, ax = plt.subplots(figsize=(16.5, 4.8), dpi=300)
    ax.axis("off")

    col_widths = [0.18, 0.08, 0.12, 0.12, 0.16, 0.24, 0.10]

    table = ax.table(
        cellText=df_summary.values,
        colLabels=df_summary.columns,
        colWidths=col_widths,
        cellLoc="center",
        loc="center"
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9.5)
    table.scale(1.0, 2.0)

    # Style table header and cells
    for (row, col), cell in table.get_celld().items():
        cell.set_edgecolor("#CCCCCC")
        cell.set_linewidth(0.7)
        if row == 0:
            cell.set_facecolor("#2B3E50")  # Standard slate navy header
            cell.set_text_props(color="white", weight="bold")
        else:
            if row % 2 == 1:
                cell.set_facecolor("#F8F9FA")
            else:
                cell.set_facecolor("#FFFFFF")

    fig.suptitle(
        f"Experimental Benchmark Scenarios ({len(scenarios)} Conditions)",
        fontsize=13.0,
        fontweight="bold",
        y=0.92
    )
    plt.tight_layout()
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    return csv_path, png_path


def save_method_configuration(
    out_dir: Path, 
    timestamp: str,
    graphics_dir: Optional[Path] = None
) -> Tuple[Path, Path]:
    """
    Saves streamlined method parameters to CSV and renders a 300 DPI publication graphic.

    Syntax:
        csv_path, png_path = save_method_configuration(out_dir, timestamp, graphics_dir=graphics_dir)

    Description:
        Provides an essential summary of algorithmic specifications and calibration parameters
        for HERA, TOPSIS, Borda Count, and Wilcoxon-Copeland, exporting both CSV and high-res graphic.

    Parameters:
        out_dir (Path): Output directory for CSV table.
        timestamp (str): Execution timestamp string.
        graphics_dir (Optional[Path]): Directory for rendered publication graphic.

    Returns:
        Tuple[Path, Path]: Tuple of (csv_path, png_path).

    Author:
        Lukas von Erdmannsdorff
    """
    records = [
        {
            "Framework": "HERA", 
            "Core Parameters": "Threshold Bootstrap\n(B = 2000)", 
            "Specification / Baseline Setup": "Sequential Hierarchical (M1 -> M2 -> M3);\nHolm-Bonferroni FWER control (alpha = 0.05)"
        },
        {
            "Framework": "TOPSIS", 
            "Core Parameters": "Criteria Weights:\nw = [0.40, 0.40, 0.20]", 
            "Specification / Baseline Setup": "Vector Normalization; Euclidean distance to\nideal solutions; relative closeness ranking"
        },
        {
            "Framework": "Borda Count", 
            "Core Parameters": "Positional Rank Sum\n(Equal weights across M1-M3)", 
            "Specification / Baseline Setup": "Sum of individual metric ranks;\nconsensus ranking by lowest rank sum"
        },
        {
            "Framework": "Wilcoxon-Copeland", 
            "Core Parameters": "Pairwise Tournament\n(Net wins minus losses)", 
            "Specification / Baseline Setup": "Paired Wilcoxon test with step-down\nHolm-Bonferroni FWER control (alpha = 0.05)"
        }
    ]
    df_methods = pd.DataFrame(records)
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / f"method_parameters_{timestamp}.csv"
    df_methods.to_csv(csv_path, index=False)

    # Render Publication Graphic in target_graphics_dir
    target_graphics_dir = graphics_dir if graphics_dir is not None else out_dir
    target_graphics_dir.mkdir(parents=True, exist_ok=True)
    png_path = target_graphics_dir / f"method_parameters_{timestamp}.png"
    pdf_path = target_graphics_dir / f"method_parameters_{timestamp}.pdf"

    fig, ax = plt.subplots(figsize=(11.5, 3.8), dpi=300)
    ax.axis("off")

    col_widths = [0.18, 0.32, 0.50]

    table = ax.table(
        cellText=df_methods.values,
        colLabels=df_methods.columns,
        colWidths=col_widths,
        cellLoc="left",
        loc="center"
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9.5)
    table.scale(1.0, 2.1)

    # Color code method frameworks
    framework_colors = {
        "HERA": "#E8F0FE",               # Soft Blue tint
        "TOPSIS": "#FEF3E2",             # Soft Amber tint
        "Borda Count": "#E6F4EA",        # Soft Green tint
        "Wilcoxon-Copeland": "#F3E8FD"   # Soft Purple tint
    }

    for (row, col), cell in table.get_celld().items():
        cell.set_edgecolor("#CCCCCC")
        cell.set_linewidth(0.7)
        if row == 0:
            cell.set_facecolor("#2B3E50")
            cell.set_text_props(color="white", weight="bold", fontsize=9.5)
        else:
            fw = df_methods.iloc[row - 1]["Framework"]
            cell.set_facecolor(framework_colors.get(fw, "#FFFFFF"))

    fig.suptitle(
        "Benchmark Method Specifications",
        fontsize=13.0,
        fontweight="bold",
        y=0.94
    )
    plt.tight_layout()
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    return csv_path, png_path


def save_candidate_configuration(
    out_dir: Path, 
    timestamp: str,
    graphics_dir: Optional[Path] = None,
    delta: float = 5.0,
    noise: float = 4.0,
    sample_size: int = 50,
    n_sims: int = 1000,
    seed: int = 123
) -> Tuple[Path, Path]:
    """
    Exports ground-truth candidate benchmark profiles to CSV and renders 300 DPI publication graphics.

    Syntax:
        csv_path, png_path = save_candidate_configuration(out_dir, timestamp, graphics_dir, delta, noise, sample_size, n_sims, seed)

    Description:
        Simulates empirical distributions for all candidate models across scales N in {10, 12, 14},
        details individual candidate roles and metrics, and exports both tabular CSVs and high-res
        publication profile curves and taxonomy tables.

    Parameters:
        out_dir (Path): Output directory for CSV table.
        timestamp (str): Execution timestamp string.
        graphics_dir (Optional[Path]): Directory for rendered publication graphics.
        delta (float): Ground-truth effect delta (%) [default: 5.0].
        noise (float): Gaussian noise sigma (%) [default: 4.0].
        sample_size (int): Cohort sample size n [default: 50].
        n_sims (int): Number of simulated trials for distribution estimation [default: 1000].
        seed (int): PRNG seed for deterministic profile generation [default: 123].

    Returns:
        Tuple[Path, Path]: Tuple of (csv_path, png_path).

    Author:
        Lukas von Erdmannsdorff
    """
    from .synthesis import get_ground_truth_means

    candidate_scales = [10, 12, 14]
    records = []
    
    for N in candidate_scales:
        means, order = get_ground_truth_means("Core", N, delta)
        rng = np.random.default_rng(seed + N)
        
        for rank_idx, cand in enumerate(order, 1):
            mu = means[cand]
            # Simulate n_sims trials of sample_size observations with noise
            trials = rng.normal(loc=mu, scale=noise, size=(n_sims, sample_size, 3))
            all_obs = trials.reshape(-1, 3)
            
            # Benchmark Ground-Truth Role Rationale
            if cand == "C1":
                role = "True Rank 1: Multi-Criterion Leader (Decisive M2 advantage over C2 by Delta+7.0%)"
            elif cand == "C2":
                role = "True Rank 2: Primary Metric Leader (Highest M1=80.0%, baseline M2/M3)"
            elif cand == "C3":
                role = "True Rank 3: Secondary Metric Step Jumper (M2 advantage over C4/C5 by Delta)"
            elif cand == "C4":
                role = "True Rank 4: Tertiary Metric Tie-Break Winner (M3 beats C5 by Delta)"
            elif cand == "C5":
                role = "True Rank 5: Tertiary Metric Tie-Break Baseline (Baseline M1-M3)"
            elif rank_idx == N:
                role = f"True Rank {N}: Compensatory Trap Model (Severe M2 deficit: 60.0%, inflated M3: 99.0%)"
            else:
                role = f"True Rank {rank_idx}: Monotonic Baseline Step (Descending M1 gradient, baseline M2/M3)"
                
            records.append({
                "Candidates_N": N,
                "Candidate": cand,
                "True_Rank": rank_idx,
                "Ranking_Role": role,
                "M1_True_Mean": round(float(mu[0]), 2),
                "M1_Simulated_Median": round(float(np.median(all_obs[:, 0])), 2),
                "M1_Simulated_IQR": round(float(np.percentile(all_obs[:, 0], 75) - np.percentile(all_obs[:, 0], 25)), 2),
                "M1_Simulated_Q1": round(float(np.percentile(all_obs[:, 0], 25)), 2),
                "M1_Simulated_Q3": round(float(np.percentile(all_obs[:, 0], 75)), 2),
                "M2_True_Mean": round(float(mu[1]), 2),
                "M2_Simulated_Median": round(float(np.median(all_obs[:, 1])), 2),
                "M2_Simulated_IQR": round(float(np.percentile(all_obs[:, 1], 75) - np.percentile(all_obs[:, 1], 25)), 2),
                "M2_Simulated_Q1": round(float(np.percentile(all_obs[:, 1], 25)), 2),
                "M2_Simulated_Q3": round(float(np.percentile(all_obs[:, 1], 75)), 2),
                "M3_True_Mean": round(float(mu[2]), 2),
                "M3_Simulated_Median": round(float(np.median(all_obs[:, 2])), 2),
                "M3_Simulated_IQR": round(float(np.percentile(all_obs[:, 2], 75) - np.percentile(all_obs[:, 2], 25)), 2),
                "M3_Simulated_Q1": round(float(np.percentile(all_obs[:, 2], 25)), 2),
                "M3_Simulated_Q3": round(float(np.percentile(all_obs[:, 2], 75)), 2),
            })
            
    df_cands = pd.DataFrame(records)
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / f"candidate_profiles_{timestamp}.csv"
    df_cands.to_csv(csv_path, index=False)
    
    # Render Publication Graphic
    target_graphics_dir = graphics_dir if graphics_dir is not None else out_dir
    target_graphics_dir.mkdir(parents=True, exist_ok=True)
    png_path = target_graphics_dir / f"candidate_profiles_{timestamp}.png"
    
    fig, axes = plt.subplots(1, 3, figsize=(19, 6.8), dpi=300, sharey=True)
    
    # Colorblind-safe palette (Wong / Okabe-Ito compliant)
    # M1: Accessible Blue (#0173b2), Circle, Solid
    # M2: High-contrast Amber (#de8f05), Square, Dashed
    # M3: Soft Purple (#cc78bc), Triangle, Dash-dot
    # Guarantees complete accessibility across protanopia, deuteranopia, tritanopia, and monochrome print
    metric_colors = {
        "M1": "#0173b2",  # Accessible Blue
        "M2": "#de8f05",  # High-contrast Amber
        "M3": "#cc78bc",  # Soft Purple
    }
    metric_labels = {
        "M1": "M1 (Primary Metric)",
        "M2": "M2 (Secondary Metric)",
        "M3": "M3 (Tertiary Metric)",
    }
    metric_styles = {
        "M1": {"marker": "o", "linestyle": "-", "linewidth": 2.4, "markersize": 7.0},
        "M2": {"marker": "s", "linestyle": "--", "linewidth": 2.0, "markersize": 6.5},
        "M3": {"marker": "^", "linestyle": "-.", "linewidth": 2.0, "markersize": 6.5},
    }
    
    for ax_idx, N in enumerate(candidate_scales):
        ax = axes[ax_idx]
        df_sub = df_cands[df_cands["Candidates_N"] == N]
        order = df_sub["Candidate"].tolist()
        x_indices = np.arange(len(order))
        
        for m_key in ["M1", "M2", "M3"]:
            col = metric_colors[m_key]
            lbl = metric_labels[m_key] if ax_idx == 0 else None
            meds = df_sub[f"{m_key}_Simulated_Median"].values
            q1s = df_sub[f"{m_key}_Simulated_Q1"].values
            q3s = df_sub[f"{m_key}_Simulated_Q3"].values
            yerr = [meds - q1s, q3s - meds]
            m_style = metric_styles[m_key]
            
            ax.errorbar(
                x_indices, meds, yerr=yerr,
                fmt=m_style["marker"],
                color=col,
                linestyle=m_style["linestyle"],
                label=lbl,
                linewidth=m_style["linewidth"],
                markersize=m_style["markersize"],
                capsize=3.5,
                capthick=1.1,
                elinewidth=1.1,
                alpha=0.92,
                zorder=3
            )
            ax.fill_between(x_indices, q1s, q3s, color=col, alpha=0.12, zorder=2)
            
        ax.set_title(f"Candidate Profiles (N = {N})", fontsize=11.5, fontweight="bold", pad=10)
        ax.set_xticks(x_indices)
        x_tick_labels = df_sub["Candidate"].tolist()
        ax.set_xticklabels(x_tick_labels, fontsize=9.5)
        ax.set_xlabel("Candidates in Ground-Truth Order (Rank 1 to N)", fontsize=11.0, fontweight="bold", labelpad=8)
        ax.set_ylabel("Metric Level (%)", fontsize=11.0, fontweight="bold", labelpad=8)
        ax.tick_params(labelleft=True)
            
        ax.set_ylim(52, 105)
        ax.grid(True, linestyle="--", alpha=0.45, zorder=1)
        ax.set_facecolor("#FAFAFA")
        
    fig.legend(
        loc="lower center",
        bbox_to_anchor=(0.5, -0.04),
        ncol=3,
        frameon=True,
        facecolor="white",
        edgecolor="#CCCCCC",
        fontsize=10.5
    )
    fig.suptitle(
        "Ground-Truth Candidate Profiles",
        fontsize=13.5,
        fontweight="bold",
        y=0.985
    )
    fig.text(
        0.5, 0.945,
        "Candidate input presentation order was randomly permuted across trials to eliminate positional bias",
        ha="center", va="top", fontsize=11.5, fontweight="bold", color="#222222"
    )
    plt.tight_layout()
    fig.subplots_adjust(bottom=0.15, top=0.845)
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    # Render Publication Candidate Table Graphic (direct visual taxonomy across N=10, 12, 14)
    table_png_path = target_graphics_dir / f"candidate_table_{timestamp}.png"

    fig_tab, axes_tab = plt.subplots(1, 3, figsize=(20, 7.2), dpi=300)

    for ax_idx, N in enumerate(candidate_scales):
        ax_t = axes_tab[ax_idx]
        ax_t.axis("off")
        means_t, order_t = get_ground_truth_means("Core", N, delta)

        table_data = []
        for rank_idx, c in enumerate(order_t, 1):
            mu_t = means_t[c]
            if c == "C1":
                role_short = "Multi-Criterion Leader (M2 Advantage)"
            elif c == "C2":
                role_short = "Primary Metric M1 Leader"
            elif c == "C3":
                role_short = "Secondary M2 Step Jumper"
            elif c == "C4":
                role_short = "Tertiary M3 Tie-Break Winner"
            elif c == "C5":
                role_short = "Tertiary M3 Tie-Break Baseline"
            elif rank_idx == N:
                role_short = "Compensatory Trap (Low M2, High M3)"
            else:
                role_short = "Monotonic Baseline Step"

            table_data.append([
                f"{c} (Rank {rank_idx})",
                f"{mu_t[0]:.1f}%",
                f"{mu_t[1]:.1f}%",
                f"{mu_t[2]:.1f}%",
                role_short
            ])

        df_t = pd.DataFrame(
            table_data,
            columns=["Candidate", "Metric M1", "Metric M2", "Metric M3", "Benchmarking Role"]
        )

        col_widths_t = [0.18, 0.16, 0.16, 0.16, 0.34]
        row_height = 0.82 / 15.0
        t_height = (len(df_t) + 1) * row_height
        t_bottom = 0.88 - t_height

        table_t = ax_t.table(
            cellText=df_t.values,
            colLabels=df_t.columns,
            colWidths=col_widths_t,
            cellLoc="center",
            bbox=[0.0, t_bottom, 1.0, t_height]
        )
        table_t.auto_set_font_size(False)
        table_t.set_fontsize(8.5)

        for (row, col), cell in table_t.get_celld().items():
            cell.set_edgecolor("#CCCCCC")
            cell.set_linewidth(0.7)
            if row == 0:
                cell.set_facecolor("#2B3E50")
                cell.set_text_props(color="white", weight="bold", fontsize=9.0)
            else:
                cand = order_t[row - 1]
                if cand == "C1":
                    cell.set_facecolor("#E8F0FE")
                elif cand == "C2":
                    cell.set_facecolor("#E6F4EA")
                elif cand == "C4":
                    cell.set_facecolor("#FEF3E2")
                elif row == len(order_t):
                    cell.set_facecolor("#FCE8E6")
                elif row % 2 == 1:
                    cell.set_facecolor("#F8F9FA")
                else:
                    cell.set_facecolor("#FFFFFF")

        ax_t.set_title(
            f"Candidate Hierarchy (N = {N})",
            fontsize=11.5,
            fontweight="bold",
            pad=8
        )

    fig_tab.suptitle(
        "Ground-Truth Benchmark Calibration",
        fontsize=13.5,
        fontweight="bold",
        y=0.98
    )
    plt.tight_layout()
    fig_tab.subplots_adjust(top=0.88, bottom=0.04)
    fig_tab.savefig(table_png_path, dpi=300, bbox_inches="tight")
    plt.close(fig_tab)

    return csv_path, png_path
