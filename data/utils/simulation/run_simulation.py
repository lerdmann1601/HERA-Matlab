"""
Main Simulation Script for HERA Methodological Ground-Truth Validation.

This script performs a comprehensive Monte Carlo study to benchmark the linear HERA 
ranking algorithm against standard Multi-Criteria Decision Making (MCDM) baselines 
(TOPSIS, Borda Count, and Copeland Scoring) under controlled ground-truth conditions.

Methodological validation dimensions:
- Known ground-truth ranking hierarchy with varying:
  * Number of candidates (N in {10, 12, 14})
  * Sample size (n in {25, 50, 100})
  * Noise level (sigma in {2.0, 4.0, 6.0, 8.0, 10.0} %)
  * Effect magnitude (Delta in {2.5, 5.0, 8.0, 12.0} % / Cliff's d calibration)
  * Metric correlation (rho in {0.0, 0.5})
- Quantitative evaluation metrics:
  * Top-Choice Recovery Rate (Accuracy of identifying true top candidate C1)
  * Complete-Rank Recovery Rate (Exact match across all N positions)
  * Kendall's Tau Rank Correlation (Ordinal alignment with ground truth)
  * Spearman's Rho Rank Correlation (Monotonic rank distance)
  * False Superiority Rate (Proportion of pairwise inversions: Inversions / (N*(N-1)/2))
  * Expected Regret (Loss in true primary metric M1 vs chosen candidate #1)
  * Compensatory Error Rate (Mistaken selection of models with severe secondary deficit despite high tertiary scores)

High-Performance Architecture & macOS / Python Standards:
- Modular Subscript Architecture: Clean separation of concerns mirroring MATLAB's
  +HERA/+analysis package (+convergence/config, simulate, save_csv, calc_pooled_csv).
- Non-Blocking Sliding-Window Pipeline: Keeps worker pool 100% saturated across 
  scenario boundaries, eliminating processor idle time and tail latencies.
- Thread-Saturated Core Isolation: Enforces single-thread execution per worker across
  Apple Accelerate (VECLIB_MAXIMUM_THREADS=1), OpenMP, OpenBLAS, MKL, and MATLAB
  (-singleCompThread), preventing core contention on Apple Silicon.
- Dynamic Resource-Aware Scaling (DRAS): Provisions worker pools calibrated to host RAM
  and CPU architecture, avoiding host memory compression and swap thrashing.
- Streaming Incremental I/O: Results CSV is flushed step-by-step to disk after every 
  iteration; scenario summaries and checkpoints are committed incrementally.
- Sandboxed Workspace Isolation: Every parallel iteration executes in an ephemeral,
  sandboxed directory cleaned up immediately in try/finally blocks (zero disk leakage).
- Strict Deterministic Reproducibility: Independent PRNG streams per (scenario, iteration)
  injected into both data synthesis and HERA engine configuration.
- Empirical Effect Size Traceability: Calculates and logs nominal Delta, empirical Cliff's d,
  mean differences, and pooled scenario statistics (matching +HERA/+analysis).
- Configuration Reporting: Exports complete scenario grid and method parameter tables
  (HERA Thresholds only, TOPSIS, Borda, Copeland) as CSV and 300 DPI publication graphics.

Outputs generated in <data_folder>/Simulation_Output_<timestamp>/:
1. configuration_<timestamp>.json     : Complete metadata, seeds, parameter configs, and environment info.
2. scenarios_<timestamp>.csv          : Detailed index and definitions of all experimental conditions.
3. scenario_configuration_<ts>.png/pdf: 300 DPI publication graphic of scenario benchmark grid.
4. method_parameters_<timestamp>.csv  : Methodological parameter specifications (streamlined).
5. method_parameters_<ts>.png/pdf     : 300 DPI publication graphic of method parameter tables.
6. candidate_profiles_<timestamp>.csv : Candidate ground-truth scaling across M1, M2, M3 (Median & IQR for N=10, 12, 14).
7. candidate_profiles_<ts>.png/pdf    : 300 DPI publication graphic of candidate ground-truth ranking scaling.
8. simulation_results_<timestamp>.csv : Full raw dataset containing every Monte Carlo trial (streamed).
9. candidate_ranks_<timestamp>.csv    : Detailed candidate-level true vs predicted ranks.
10. global_summary_<timestamp>.csv    : Aggregated scenario metrics (Mean, Median, IQR, 95% CI, Min, Max).
11. pooled_results_<timestamp>.csv    : Grand pooled summary across all conditions (matching calc_pooled_csv.m).
12. simulation_log_<timestamp>.txt   : Execution diary logging runtimes and diagnostic details.
13. checkpoint_<timestamp>.json       : Live state tracker for progress monitoring and resumption.

Usage:
    python3 run_simulation.py

Author: Lukas von Erdmannsdorff
"""

import os
# Suppress noisy third-party package notifications to maintain pristine console logs
os.environ["OUTDATED_IGNORE"] = "1"

# Enforce strict single-thread execution per worker across all numerical libraries
# on macOS (Apple Silicon / Accelerate vecLib, OpenBLAS, MKL, OpenMP, NumExpr)
# to prevent cross-process CPU contention and thread oversubscription.
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import sys
import time
import json
import shutil
import atexit
import warnings
import platform
from datetime import datetime
from pathlib import Path

# Filter benign user warnings from third-party libraries during parallel bootstrap
warnings.filterwarnings("ignore", category=UserWarning)

import pandas as pd

# Path setup & Sub-package discovery
PACKAGE_DIR = Path(__file__).parent.resolve()
UTILS_DIR = PACKAGE_DIR.parent.resolve()
if str(UTILS_DIR) not in sys.path:
    sys.path.insert(0, str(UTILS_DIR))
if str(PACKAGE_DIR) not in sys.path:
    sys.path.insert(0, str(PACKAGE_DIR))

# Import modular components from simulation sub-package
try:
    from .config import (
        DATA_DIR, TEMP_DIR,
        NUM_CANDIDATES_LIST, SAMPLE_SIZES, NOISE_LEVELS,
        EFFECT_CALIBRATION, EFFECT_MAGNITUDES, DEFAULT_EFFECT, DEFAULT_NOISE,
        CORRELATIONS, DEFAULT_CORRELATION, ITERATIONS, SEED,
        HERA_CONFIG, TOPSIS_WEIGHTS,
        RESULT_COLUMNS, METRIC_COLS, STATISTICAL_EFFECT_COLS, CANDIDATE_RANK_COLUMNS,
        make_scenario, build_scenarios_grid
    )
    from .system import (
        suppress_stdout_stderr, format_time_duration,
        get_system_ram_gb, get_optimal_worker_count,
        setup_logger, print_study_header, print_study_completion
    )
    from .engine import (
        HERAExecutor, get_worker_executor
    )
    from .synthesis import (
        compute_cliffs_delta, get_ground_truth_means, generate_synthetic_data
    )
    from .mcdm import (
        run_hera_ranking, run_topsis_baseline, run_borda_baseline,
        run_wilcoxon_copeland_baseline, evaluate_ranking
    )
    from .reporting import (
        append_results_to_csv, compute_aggregate_summary, compute_pooled_summary,
        compute_benchmark_summary, compute_superordinate_summary,
        compute_detailed_pooled_statistics, export_summary_csvs, export_superordinate_csvs,
        update_global_summary_file, save_checkpoint,
        save_scenario_configuration, save_method_configuration, save_candidate_configuration
    )
    from .simulate import (
        execute_single_iteration, run_sliding_window_pipeline
    )
except ImportError:
    from simulation.config import (
        DATA_DIR, TEMP_DIR,
        NUM_CANDIDATES_LIST, SAMPLE_SIZES, NOISE_LEVELS,
        EFFECT_CALIBRATION, EFFECT_MAGNITUDES, DEFAULT_EFFECT, DEFAULT_NOISE,
        CORRELATIONS, DEFAULT_CORRELATION, ITERATIONS, SEED,
        HERA_CONFIG, TOPSIS_WEIGHTS,
        RESULT_COLUMNS, METRIC_COLS, STATISTICAL_EFFECT_COLS, CANDIDATE_RANK_COLUMNS,
        make_scenario, build_scenarios_grid
    )
    from simulation.system import (
        suppress_stdout_stderr, format_time_duration,
        get_system_ram_gb, get_optimal_worker_count,
        setup_logger, print_study_header, print_study_completion
    )
    from simulation.engine import (
        HERAExecutor, get_worker_executor
    )
    from simulation.synthesis import (
        compute_cliffs_delta, get_ground_truth_means, generate_synthetic_data
    )
    from simulation.mcdm import (
        run_hera_ranking, run_topsis_baseline, run_borda_baseline,
        run_wilcoxon_copeland_baseline, evaluate_ranking
    )
    from simulation.reporting import (
        append_results_to_csv, compute_aggregate_summary, compute_pooled_summary,
        compute_benchmark_summary, compute_superordinate_summary,
        compute_detailed_pooled_statistics, export_summary_csvs, export_superordinate_csvs,
        update_global_summary_file, save_checkpoint,
        save_scenario_configuration, save_method_configuration, save_candidate_configuration
    )
    from simulation.simulate import (
        execute_single_iteration, run_sliding_window_pipeline
    )

# Re-export all symbols for backwards-compatibility
__all__ = [
    "DATA_DIR", "TEMP_DIR",
    "NUM_CANDIDATES_LIST", "SAMPLE_SIZES", "NOISE_LEVELS",
    "EFFECT_CALIBRATION", "EFFECT_MAGNITUDES", "DEFAULT_EFFECT", "DEFAULT_NOISE",
    "CORRELATIONS", "DEFAULT_CORRELATION", "ITERATIONS", "SEED",
    "HERA_CONFIG", "TOPSIS_WEIGHTS",
    "RESULT_COLUMNS", "METRIC_COLS", "STATISTICAL_EFFECT_COLS", "CANDIDATE_RANK_COLUMNS",
    "make_scenario", "build_scenarios_grid",
    "suppress_stdout_stderr", "format_time_duration",
    "get_system_ram_gb", "get_optimal_worker_count",
    "setup_logger", "print_study_header", "print_study_completion",
    "HERAExecutor", "get_worker_executor",
    "compute_cliffs_delta", "get_ground_truth_means", "generate_synthetic_data",
    "run_hera_ranking", "run_topsis_baseline", "run_borda_baseline",
    "run_wilcoxon_copeland_baseline", "evaluate_ranking",
    "append_results_to_csv", "compute_aggregate_summary", "compute_pooled_summary",
    "compute_superordinate_summary", "compute_detailed_pooled_statistics", "export_superordinate_csvs",
    "update_global_summary_file", "save_checkpoint",
    "save_scenario_configuration", "save_method_configuration", "save_candidate_configuration",
    "execute_single_iteration", "run_sliding_window_pipeline",
    "main"
]


def main() -> None:
    """
    Execute the HERA Ground-Truth Monte Carlo Benchmark Simulation Suite.

    Syntax:
        main()

    Description:
        Orchestrates the end-to-end execution of the HERA Monte Carlo benchmark
        simulation suite. Parses CLI arguments and environment configurations,
        initializes the dedicated timestamped output folder structure (CSVs/,
        Graphics/, Reports/), validates the HERA computational backend, generates
        configuration tables and publication artifacts, runs the sliding-window
        parallel simulation pipeline, computes statistical and grand-pooled summaries,
        and generates high-resolution figures and multi-page PDF reports.

    Parameters:
        None. CLI arguments are parsed internally via argparse.

    Returns:
        None

    Author:
        Lukas von Erdmannsdorff
    """
    # 1. Parse CLI Arguments & Runtime Configurations
    import argparse
    parser = argparse.ArgumentParser(
        description="HERA Methodological Ground-Truth Benchmark Simulation Suite."
    )
    parser.add_argument(
        "-i", "--iterations", "-n",
        type=int,
        default=ITERATIONS,
        help=f"Number of Monte Carlo iterations per scenario (default: {ITERATIONS}, configurable via HERA_SIM_ITERATIONS env var)"
    )
    parser.add_argument(
        "-w", "--workers",
        type=int,
        default=None,
        help="Number of parallel worker processes (default: auto-detected based on system RAM/CPU)"
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=SEED,
        help=f"Base PRNG seed (default: {SEED}, configurable via HERA_SIM_SEED env var)"
    )
    parser.add_argument(
        "-o", "--output", "--output-dir",
        type=str,
        default=os.environ.get("HERA_SIM_OUTPUT_DIR", None),
        help="Custom output directory path (default: Simulation_Output_<timestamp> inside data directory, or via HERA_SIM_OUTPUT_DIR env var)"
    )
    args, _ = parser.parse_known_args()
    sim_iterations = args.iterations
    sim_seed = args.seed

    # 2. Initialize Dedicated Timestamped Directory Hierarchy (CSVs, Graphics, Reports)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    if args.output:
        custom_out = Path(args.output).expanduser()
        if not custom_out.is_absolute():
            custom_out = (DATA_DIR / custom_out).resolve()
        run_dir = custom_out
    else:
        run_dir = DATA_DIR / f"Simulation_Output_{timestamp}"
    run_dir.mkdir(parents=True, exist_ok=True)
    
    csv_dir = run_dir / "CSVs"
    csv_dir.mkdir(parents=True, exist_ok=True)
    graphics_dir = run_dir / "Graphics"
    graphics_dir.mkdir(parents=True, exist_ok=True)
    reports_dir = run_dir / "Reports"
    reports_dir.mkdir(parents=True, exist_ok=True)

    # 3. Setup Ephemeral Sandboxed Scratch Space & Process Cleanup Hook
    run_temp_dir = run_dir / "temp_workspace"
    run_temp_dir.mkdir(parents=True, exist_ok=True)

    def cleanup_run_temp() -> None:
        """Removes temporary scratch directories upon process termination."""
        if run_temp_dir.exists():
            shutil.rmtree(run_temp_dir, ignore_errors=True)
        if TEMP_DIR.exists():
            shutil.rmtree(TEMP_DIR, ignore_errors=True)
    atexit.register(cleanup_run_temp)

    # 4. Initialize Synchronized Execution Logger
    log_file = run_dir / f"simulation_log_{timestamp}.txt"
    logger = setup_logger(log_file)
    logger.info(f"[Output] Dedicated simulation output folder created: {run_dir}")
    logger.info(f"[Output] CSV Data Folder: {csv_dir}")
    logger.info(f"[Output] Graphics Folder: {graphics_dir}")
    logger.info(f"[Output] Reports Folder:  {reports_dir}")

    # 5. Profile System Hardware & Allocate Optimal DRAS Worker Pool
    ram_gb = get_system_ram_gb()
    num_workers = args.workers if args.workers is not None and args.workers > 0 else get_optimal_worker_count(ram_gb)

    # 6. Validate HERA Computational Backend (hera-matlab / Local MATLAB CLI)
    logger.info("[Environment] Validating HERA computational engine...")
    try:
        probe_executor = HERAExecutor(logger)
        hera_mode = probe_executor.mode
        probe_executor.terminate()
    except Exception as e:
        logger.error(f"[Error] Failed to initialize HERA backend: {e}")
        logger.error("Please verify that hera_matlab or local MATLAB is installed.")
        sys.exit(1)

    # 7. Construct Orthogonal Experimental Scenario Grid
    all_scenarios = build_scenarios_grid()
    scenarios = list(all_scenarios)

    if os.environ.get("HERA_SMOKE_TEST") == "1":
        logger.info("[Smoke Test] HERA_SMOKE_TEST detected. Restricting execution to 1 scenario for end-to-end verification.")
        scenarios = scenarios[:1]
        num_workers = min(num_workers, sim_iterations)

    # 8. Print Formatted Academic Study Header
    print_study_header(logger, timestamp, all_scenarios, sim_iterations, ram_gb, num_workers, hera_mode)

    # 9. Export Benchmark Parameter Tables & 300 DPI Publication Graphics
    scenarios_csv, scenarios_png = save_scenario_configuration(all_scenarios, csv_dir, timestamp, graphics_dir=graphics_dir)
    logger.info(f"[Configuration] Generated {len(all_scenarios)} experimental scenarios.")
    logger.info(f"[Configuration] Scenario index CSV saved: {scenarios_csv}")
    logger.info(f"[Configuration] Scenario grid Graphic saved: {scenarios_png}")

    methods_csv, methods_png = save_method_configuration(csv_dir, timestamp, graphics_dir=graphics_dir)
    logger.info(f"[Configuration] Method parameters CSV saved: {methods_csv}")
    logger.info(f"[Configuration] Method parameters Graphic saved: {methods_png}")

    cands_csv, cands_png = save_candidate_configuration(csv_dir, timestamp, graphics_dir=graphics_dir)
    logger.info(f"[Configuration] Candidate ground-truth profiles CSV saved: {cands_csv}")
    logger.info(f"[Configuration] Candidate ground-truth profiles Graphic saved: {cands_png}")

    # 10. Export FAIR Study Metadata JSON
    config_metadata = {
        "study_name": "HERA Methodological Ground-Truth Benchmark Study",
        "timestamp": timestamp,
        "base_seed": sim_seed,
        "iterations_per_scenario": sim_iterations,
        "parallel_workers": num_workers,
        "system": {
            "platform": platform.platform(),
            "python_version": sys.version,
            "machine": platform.machine(),
            "detected_ram_gb": round(ram_gb, 2)
        },
        "hera_config": HERA_CONFIG,
        "topsis_weights": TOPSIS_WEIGHTS.tolist(),
        "total_scenarios": len(all_scenarios),
        "executed_scenarios": len(scenarios),
        "methods_compared": ["HERA", "TOPSIS", "Borda", "Copeland"],
        "metrics_evaluated": METRIC_COLS,
        "statistical_effect_metrics": STATISTICAL_EFFECT_COLS
    }
    config_json_path = run_dir / f"configuration_{timestamp}.json"
    with open(config_json_path, "w", encoding="utf-8") as f:
        json.dump(config_metadata, f, indent=2)
    logger.info(f"[FAIR] Study metadata exported to: {config_json_path}")

    # 11. Initialize Streaming CSV Result Files with Headers (Early Allocation)
    results_csv = csv_dir / f"simulation_results_{timestamp}.csv"
    with open(results_csv, "w", newline="", encoding="utf-8") as f:
        pd.DataFrame(columns=RESULT_COLUMNS).to_csv(f, index=False)
    logger.info(f"[I/O] Streaming results initialized with headers: {results_csv}")

    candidate_ranks_csv = csv_dir / f"candidate_ranks_{timestamp}.csv"
    with open(candidate_ranks_csv, "w", newline="", encoding="utf-8") as f:
        pd.DataFrame(columns=CANDIDATE_RANK_COLUMNS).to_csv(f, index=False)
    logger.info(f"[I/O] Streaming candidate ranks initialized with headers: {candidate_ranks_csv}")

    summary_csv = csv_dir / f"global_summary_{timestamp}.csv"
    pooled_csv = csv_dir / f"pooled_results_{timestamp}.csv"
    checkpoint_json = run_dir / f"checkpoint_{timestamp}.json"
    repo_root = UTILS_DIR.parent.parent.parent.resolve()

    # 12. Execute Sliding-Window Parallel Simulation Pipeline
    all_results_cache, elapsed_total, interrupted = run_sliding_window_pipeline(
        scenarios=scenarios,
        iterations=sim_iterations,
        base_seed=sim_seed,
        num_workers=num_workers,
        run_temp_dir=run_temp_dir,
        repo_root=repo_root,
        logger=logger,
        results_csv=results_csv,
        candidate_ranks_csv=candidate_ranks_csv,
        summary_csv=summary_csv,
        checkpoint_json=checkpoint_json,
        get_ground_truth_means_fn=get_ground_truth_means
    )

    # 13. Finalize Statistical Summaries & Grand-Pooled Datasets
    benchmark_csv = None
    detailed_pooled_csv = None
    if Path(results_csv).exists() and os.path.getsize(results_csv) > 0:
        try:
            df_final_raw = pd.read_csv(results_csv)
            if not df_final_raw.empty:
                update_global_summary_file(summary_csv, df_final_raw)
                logger.info(f"[Results] Global statistical summary confirmed at: {summary_csv}")

                # Export Grand Pooled and Executive Benchmark Summary CSVs
                csv_map = export_summary_csvs(df_final_raw, csv_dir, timestamp, root_dir=run_dir)
                benchmark_csv = csv_map.get("benchmark_summary")
                detailed_pooled_csv = csv_map.get("pooled_detailed")
                pooled_csv = csv_map.get("pooled_summary", pooled_csv)

                if benchmark_csv and benchmark_csv.exists():
                    logger.info(f"[Results] Executive benchmark summary confirmed at: {benchmark_csv}")
                if detailed_pooled_csv and detailed_pooled_csv.exists():
                    logger.info(f"[Results] Detailed non-parametric pooled breakdown confirmed at: {detailed_pooled_csv}")
                logger.info(f"[Results] Grand pooled summary confirmed at: {pooled_csv}")
        except Exception as e:
            logger.warning(f"[Warning] Could not finalize summary CSVs: {e}")

    logger.info("---------------------------------------------------------------------------------------------------------")
    logger.info("[Environment] Parallel workers terminated cleanly. Cleaned up temporary workspaces.")

    cleanup_run_temp()

    # 14. Automated Publication Plotting & Multi-Page PDF Compilation
    master_pdf_path = run_dir / f"Global_Summary_{timestamp}.pdf"
    if summary_csv.exists() and os.path.getsize(summary_csv) > 0:
        try:
            logger.info(f"[Visualization] Launching publication figures & multi-page PDF generation...")
            try:
                from . import plot_results
            except ImportError:
                import plot_results
            plot_results.generate_all_plots(
                results_file=results_csv, 
                summary_file=summary_csv, 
                graphics_dir=graphics_dir, 
                reports_dir=reports_dir, 
                timestamp=timestamp
            )
            logger.info(f"[Visualization] All publication figures successfully generated in: {graphics_dir}")
            logger.info(f"[Visualization] PDF reports successfully generated in: {reports_dir}")
        except Exception as e:
            logger.warning(f"[Visualization] Automated plotting could not be completed: {e}")

    # 15. Print Study Completion Summary & Formal Citation
    print_study_completion(
        logger=logger,
        t_duration=elapsed_total,
        config_json=config_json_path,
        scenarios_csv=scenarios_csv,
        results_csv=results_csv,
        summary_csv=summary_csv,
        pooled_csv=pooled_csv,
        checkpoint_json=checkpoint_json,
        log_file=log_file,
        superordinate_csv=benchmark_csv,
        global_pdf=master_pdf_path,
        plots_dir=graphics_dir,
        reports_dir=reports_dir,
        csv_dir=csv_dir,
        candidate_profiles_csv=cands_csv
    )

    # 16. Update Main Output Symlink for Seamless Inspection
    latest_link = DATA_DIR / "Simulation_Output"
    try:
        if latest_link.is_symlink() or latest_link.is_file():
            latest_link.unlink()
        if not latest_link.exists():
            target_link = run_dir.name if run_dir.parent.resolve() == DATA_DIR.resolve() else run_dir.resolve()
            latest_link.symlink_to(target_link)
            logger.info(f"[Output] Created symlink for easy access: {latest_link} -> {target_link}")
    except Exception:
        pass


if __name__ == "__main__":
    main()
