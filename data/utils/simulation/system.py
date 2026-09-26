"""
Hardware Profiling, DRAS Worker Allocation, Logging & Academic Banners.

Provides resource management calibrated to Apple Silicon / macOS host memory,
suppression utilities, and console/diary formatting.

Author: Lukas von Erdmannsdorff
"""

import os
import sys
import platform
import logging
import subprocess
from pathlib import Path
from typing import List, Dict, Any, Optional
from contextlib import contextmanager

from .config import (
    SEED, NUM_CANDIDATES_LIST, SAMPLE_SIZES, NOISE_LEVELS,
    HERA_CONFIG
)


@contextmanager
def suppress_stdout_stderr():
    """
    Suppresses stdout and stderr at the OS/C level to silence native library prints.

    Syntax:
        with suppress_stdout_stderr():
            ...

    Description:
        Redirects standard output and standard error file descriptors (1 and 2) to os.devnull
        during the context block, preventing native C-level library logging from cluttering the console.

    Author:
        Lukas von Erdmannsdorff
    """
    sys.stdout.flush()
    sys.stderr.flush()
    try:
        null_fd = os.open(os.devnull, os.O_RDWR)
        saved_stdout_fd = os.dup(1)
        saved_stderr_fd = os.dup(2)
        os.dup2(null_fd, 1)
        os.dup2(null_fd, 2)
        try:
            yield
        finally:
            sys.stdout.flush()
            sys.stderr.flush()
            os.dup2(saved_stdout_fd, 1)
            os.dup2(saved_stderr_fd, 2)
            os.close(saved_stdout_fd)
            os.close(saved_stderr_fd)
            os.close(null_fd)
    except Exception:
        yield


def format_time_duration(seconds: float) -> str:
    """
    Formats seconds into human-readable duration string.

    Syntax:
        dur_str = format_time_duration(seconds)

    Description:
        Converts a raw elapsed second count into a zero-padded, human-readable string
        (e.g., '02h 15m 30s' or '01m 45s').

    Parameters:
        seconds (float): Elapsed time in seconds.

    Returns:
        str: Formatted duration string.

    Author:
        Lukas von Erdmannsdorff
    """
    s = int(round(seconds))
    hrs = s // 3600
    mins = (s % 3600) // 60
    secs = s % 60
    if hrs > 0:
        return f"{hrs:02d}h {mins:02d}m {secs:02d}s"
    elif mins > 0:
        return f"{mins:02d}m {secs:02d}s"
    else:
        return f"{secs:02d}s"


def get_system_ram_gb() -> float:
    """
    Detects physical system memory in Gigabytes.

    Syntax:
        ram_gb = get_system_ram_gb()

    Description:
        Queries host operating system telemetry (sysctl on macOS Darwin, /proc/meminfo on Linux)
        to identify total physical random-access memory.

    Returns:
        float: Total physical host memory in GB (defaults to 16.0 GB on detection failure).

    Author:
        Lukas von Erdmannsdorff
    """
    try:
        if platform.system() == "Darwin":
            out = subprocess.check_output(["sysctl", "-n", "hw.memsize"]).strip()
            return float(out) / (1024 ** 3)
        elif platform.system() == "Linux":
            with open("/proc/meminfo") as f:
                for line in f:
                    if line.startswith("MemTotal:"):
                        return float(line.split()[1]) / (1024 ** 2)
    except Exception:
        pass
    return 16.0


def get_optimal_worker_count(ram_gb: float) -> int:
    """
    Calculates safe parallel worker count via Dynamic Resource-Aware Scaling (DRAS).

    Syntax:
        n_workers = get_optimal_worker_count(ram_gb)

    Description:
        Balances CPU core count against RAM requirements per worker (~2.0 - 2.5 GB peak).
        Reserves memory and cores for orchestrator and OS stability. Allows manual
        override via HERA_SIM_WORKERS environment variable.

    Parameters:
        ram_gb (float): Detected physical system memory in GB.

    Returns:
        int: Recommended number of concurrent parallel worker processes.

    Author:
        Lukas von Erdmannsdorff
    """
    env_override = os.environ.get("HERA_SIM_WORKERS")
    if env_override:
        try:
            val = int(env_override)
            if val > 0:
                return val
        except ValueError:
            pass

    # Reserve 3.0 GB RAM for OS and main orchestrator process
    usable_ram = max(2.0, ram_gb - 3.0)
    ram_workers = max(1, int(usable_ram // 2.2))

    total_cpus = os.cpu_count() or 4
    # Reserve 2 CPU cores for system responsiveness and main orchestrator thread
    cpu_workers = max(1, total_cpus - 2 if total_cpus > 2 else 1)

    return min(ram_workers, cpu_workers)


def setup_logger(log_file: Path) -> logging.Logger:
    """
    Configures synchronized console and file logging with clean formatting.

    Syntax:
        logger = setup_logger(log_file)

    Description:
        Initializes a dual-handler logger directing identical log records to both
        standard output (console) and an execution diary text file.

    Parameters:
        log_file (Path): Destination log text file path.

    Returns:
        logging.Logger: Configured logger instance.

    Author:
        Lukas von Erdmannsdorff
    """
    logger = logging.getLogger("HERA_Simulation")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()

    fh = logging.FileHandler(log_file, encoding="utf-8")
    fh.setLevel(logging.INFO)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)

    formatter = logging.Formatter("%(message)s")
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)

    logger.addHandler(fh)
    logger.addHandler(ch)
    return logger


def print_study_header(
    logger: logging.Logger,
    timestamp: str,
    scenarios: List[Dict[str, Any]],
    iterations: int,
    ram_gb: float,
    num_workers: int,
    hera_mode: str
) -> None:
    """
    Prints a structured academic ASCII header matching HERA's +analysis publication standard.

    Syntax:
        print_study_header(logger, timestamp, scenarios, iterations, ram_gb, num_workers, hera_mode)

    Parameters:
        logger (logging.Logger): Logger instance for console and file output.
        timestamp (str): Execution timestamp string.
        scenarios (List[Dict[str, Any]]): List of scenario configuration dictionaries.
        iterations (int): Monte Carlo runs per scenario condition.
        ram_gb (float): Host physical memory in GB.
        num_workers (int): Number of allocated worker processes.
        hera_mode (str): Active computational backend string.

    Returns:
        None

    Author:
        Lukas von Erdmannsdorff
    """
    total_runs = len(scenarios) * iterations
    ram_per_worker = (ram_gb - 3.0) / max(1, num_workers)
    
    div_heavy = "=" * 70
    div_light = "-" * 70

    logger.info(div_heavy)
    logger.info("   HERA Multi-Criteria Ground-Truth Simulation Study")
    logger.info(f"   Benchmark Setup: {len(scenarios)} Scenarios | {iterations} Runs/Condition | 4 Decision Frameworks")
    logger.info(div_heavy)
    logger.info(" System & Computational Infrastructure:")
    logger.info(f"  Timestamp:             {timestamp}")
    logger.info(f"  Host Platform:         {platform.system()} {platform.release()} ({platform.machine()})")
    logger.info(f"  Physical Memory:       {ram_gb:.1f} GB RAM (~{ram_per_worker:.1f} GB/worker allocated)")
    logger.info(f"  Parallel Processing:   {num_workers} concurrent worker processes (Sliding-Window Pipeline)")
    logger.info(f"  Python Environment:    {platform.python_version()}")
    logger.info(f"  HERA Engine Backend:   {hera_mode} (MATLAB CLI -singleCompThread)")
    logger.info(f"  PRNG Random Seed:      {SEED} (Deterministic PCG64 stream per condition)")
    logger.info(f"  Total Evaluations:     {total_runs} Monte Carlo runs ({total_runs * 4} ranking outputs)")
    logger.info(div_light)
    logger.info(" Experimental Benchmarking Grid:")
    logger.info(f"  Candidate Scales (N):  {', '.join(map(str, NUM_CANDIDATES_LIST))} candidates")
    logger.info(f"  Sample Sizes (n):      {', '.join(map(str, SAMPLE_SIZES))} paired evaluation folds")
    logger.info(f"  Noise Levels (sigma):  {', '.join([f'{s}%' for s in NOISE_LEVELS])}")
    logger.info("  Effect Size Magnitudes (Delta -> Calibrated Median Cliff's d):")
    logger.info("    * Small:  Delta =  2.5%  ==>  Cliff's d ~ 0.25 (range: 0.20 - 0.30)")
    logger.info("    * Medium: Delta =  5.0%  ==>  Cliff's d ~ 0.50 (range: 0.45 - 0.60)")
    logger.info("    * Large:  Delta =  8.0%  ==>  Cliff's d ~ 0.80 (range: 0.75 - 0.90)")
    logger.info("    * Stress: Delta = 12.0%  ==>  Cliff's d ~ 0.95 (range: > 0.90)")
    logger.info(f"  Inter-Metric Coupling: rho in {{0.0 (Independent), 0.5 (Correlated)}}")
    logger.info(div_light)
    logger.info(" Synthetic Multi-Metric Ground-Truth Paradigm:")
    logger.info("  Primary Metric (M1):   Primary benchmark ranking criterion (Baseline: ~80.0%)")
    logger.info("  Secondary Metric (M2): Secondary priority criterion (Decisive C1 advantage: 87.0%)")
    logger.info("  Tertiary Metric (M3):  Tie-breaking criterion (Decisive C4 vs C5 resolution: 85.0%)")
    logger.info("  Trap Candidate (C_N):  Compensatory flaw model (Deficit M2: 60.0%, Inflated M3: 99.0%)")
    logger.info(div_light)
    logger.info(" Statistical & Algorithmic Configurations:")
    logger.info(f"  HERA Ranking Mode:     {HERA_CONFIG['ranking_mode']} (Sequential Non-Compensatory Hierarchical)")
    logger.info(f"  HERA Threshold Bootstr:Percentile Null Bootstrap (B_thr = {HERA_CONFIG['manual_B_thr']})")
    logger.info("  HERA Decision Logic:   Holm-Bonferroni FWER Control (alpha = 0.05), Cliff's d Dominance")
    logger.info("  MCDM Baseline Methods: TOPSIS (Weights: [0.4, 0.4, 0.2]), Borda Count,")
    logger.info("                         Wilcoxon-Copeland (Step-Down Holm FWER Control, alpha = 0.05)")
    logger.info(div_heavy + "\n")


def print_study_completion(
    logger: logging.Logger,
    t_duration: float,
    config_json: Path,
    scenarios_csv: Path,
    results_csv: Path,
    summary_csv: Path,
    pooled_csv: Path,
    checkpoint_json: Path,
    log_file: Path,
    superordinate_csv: Optional[Path] = None,
    global_pdf: Optional[Path] = None,
    plots_dir: Optional[Path] = None,
    reports_dir: Optional[Path] = None,
    csv_dir: Optional[Path] = None,
    candidate_profiles_csv: Optional[Path] = None
) -> None:
    """
    Prints completion summary and official academic citation string.

    Syntax:
        print_study_completion(logger, t_duration, config_json, ...)

    Description:
        Logs final execution metrics, file locations for reports, plots, and CSV tables,
        and provides formal BibTeX / DOI citation metadata.

    Parameters:
        logger (logging.Logger): Logger instance.
        t_duration (float): Total study execution duration in seconds.
        config_json (Path): Path to configuration JSON.
        scenarios_csv (Path): Path to scenario index CSV.
        results_csv (Path): Path to raw simulation results CSV.
        summary_csv (Path): Path to global summary CSV.
        pooled_csv (Path): Path to pooled summary CSV.
        checkpoint_json (Path): Path to final state checkpoint JSON.
        log_file (Path): Path to execution log text file.
        superordinate_csv (Optional[Path]): Optional path to superordinate summary CSV.
        global_pdf (Optional[Path]): Optional path to master PDF report.
        plots_dir (Optional[Path]): Optional path to generated graphics directory.
        reports_dir (Optional[Path]): Optional path to generated reports directory.
        csv_dir (Optional[Path]): Optional path to CSVs directory.
        candidate_profiles_csv (Optional[Path]): Optional path to candidate profiles CSV.

    Returns:
        None

    Author:
        Lukas von Erdmannsdorff
    """
    time_str = format_time_duration(t_duration)
    div_heavy = "=" * 70
    div_light = "-" * 70

    logger.info("\n" + div_heavy)
    logger.info("   Benchmark Study Completed Successfully")
    logger.info(div_heavy)
    logger.info(f" Total Study Duration:   {time_str}")
    if global_pdf is not None and global_pdf.exists():
        logger.info(f" Master Report (PDF):    {global_pdf}")
    if reports_dir is not None and reports_dir.exists():
        logger.info(f" Reports Folder:         {reports_dir}")
    if plots_dir is not None and plots_dir.exists():
        logger.info(f" Graphics Folder:        {plots_dir}")
    if csv_dir is not None and csv_dir.exists():
        logger.info(f" CSV Results Folder:     {csv_dir}")
    logger.info(f" Configuration (JSON):   {config_json}")
    logger.info(f" Scenario Index (CSV):   {scenarios_csv}")
    if candidate_profiles_csv is not None and candidate_profiles_csv.exists():
        logger.info(f" Candidate Profiles (CSV):{candidate_profiles_csv} [Ground-Truth Scaling]")
    logger.info(f" Raw Results (CSV):      {results_csv} [Streamed]")
    logger.info(f" Global Summary (CSV):   {summary_csv} [Aggregated]")
    if superordinate_csv is not None and superordinate_csv.exists():
        logger.info(f" Benchmark Summary (CSV):{superordinate_csv} [Suite Overview]")
    logger.info(f" Pooled Summary (CSV):   {pooled_csv} [FAIR Standard]")
    logger.info(f" State Checkpoint:       {checkpoint_json}")
    logger.info(f" Execution Log (TXT):    {log_file}")
    logger.info(div_light)
    logger.info(" If you use this software or benchmark in your research, please cite:")
    logger.info(" von Erdmannsdorff, L. (2026). HERA: Hierarchical-Compensatory, Effect-Size-Driven Ranking Algorithm")
    logger.info(" Preprint: https://doi.org/10.2139/ssrn.7227474")
    logger.info(" Software (Version 1.4.7): https://doi.org/10.5281/zenodo.18274870")
    logger.info(div_heavy + "\n")
