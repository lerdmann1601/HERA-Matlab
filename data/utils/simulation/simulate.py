"""
Parallel Simulation Core & Sandboxed Task Worker (Sliding-Window Pipeline).

Implements the isolated Monte Carlo trial execution with automatic retries and
the high-performance sliding-window parallel scheduler ensuring 100% worker saturation.

Author: Lukas von Erdmannsdorff
"""

import sys
import time
import uuid
import shutil
import signal
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Any, Optional
from collections import deque
import multiprocessing as mp
import concurrent.futures
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime
import numpy as np
import pandas as pd

from .config import RESULT_COLUMNS, CANDIDATE_RANK_COLUMNS
from .system import format_time_duration
from .engine import get_worker_executor
from .synthesis import compute_cliffs_delta, generate_synthetic_data
from .mcdm import (
    run_hera_ranking, run_topsis_baseline, run_borda_baseline,
    run_wilcoxon_copeland_baseline, evaluate_ranking
)
from .reporting import append_results_to_csv, update_global_summary_file, save_checkpoint


def clear_terminal_line() -> None:
    """
    Clears the entire current terminal line cleanly across any window width.

    Syntax:
        clear_terminal_line()

    Description:
        Emits ANSI escape sequences to clear the active line buffer on interactive TTY devices,
        preventing trailing artifacts during high-frequency progress reporting.

    Returns:
        None

    Author:
        Lukas von Erdmannsdorff
    """
    if sys.stdout.isatty():
        cols = shutil.get_terminal_size((120, 20)).columns
        sys.stdout.write("\r\033[2K" + " " * min(cols, 160) + "\r\033[2K")
        sys.stdout.flush()


def log_scenario_header(
    sc_track: Dict[str, Any], 
    total_scenarios: int, 
    logger: logging.Logger, 
    is_tty: bool
) -> None:
    """
    Logs standardized, aligned scenario announcement header.

    Syntax:
        log_scenario_header(sc_track, total_scenarios, logger, is_tty)

    Description:
        Formats and prints a detailed banner announcing the commencement of a new
        experimental benchmarking condition, displaying suite name, candidate count N,
        cohort sample size n, noise level sigma, effect magnitude Delta, and correlation rho.

    Parameters:
        sc_track (Dict[str, Any]): Scenario tracking metadata dictionary.
        total_scenarios (int): Total number of scenarios in the grid.
        logger (logging.Logger): Logger instance.
        is_tty (bool): Whether output stream is an interactive terminal.

    Returns:
        None

    Author:
        Lukas von Erdmannsdorff
    """
    if is_tty:
        clear_terminal_line()
    time_str_now = datetime.now().strftime("%H:%M:%S")
    sc_pct = (sc_track["sc_idx"] / total_scenarios) * 100
    sc_id = sc_track["sc"]["ScenarioID"]
    logger.info(
        f"\n[{time_str_now}] --- Commencing Scenario {sc_track['sc_idx']:02d}/{total_scenarios:02d} ({sc_pct:4.1f}%): {sc_id} ---"
    )
    logger.info(
        f"  Suite: {sc_track['sc']['Suite']:<28} | Candidates (N): {sc_track['sc']['Candidates']:02d}  | Sample Size (n): {sc_track['sc']['SampleSize']:03d}\n"
        f"  Noise (sigma): {sc_track['sc']['Noise']:4.1f}%          | Effect (Delta): {sc_track['sc']['Delta']:4.1f}% | Correlation (rho): {sc_track['sc']['Correlation']:3.1f}"
    )


def execute_single_iteration(task_args: Dict[str, Any]) -> Dict[str, Any]:
    """
    Executes a single Monte Carlo trial in an isolated worker workspace.

    Syntax:
        result = execute_single_iteration(task_args)

    Description:
        Runs an isolated Monte Carlo evaluation:
          1. Generates synthetic multi-metric data with seeded PRNG stream.
          2. Computes empirical effect sizes (Cliff's d and Delta on M1 and M2).
          3. Evaluates all MCDM methods (HERA, TOPSIS, Borda, Wilcoxon-Copeland).
          4. Scores ranking fidelity against ground truth (recovery, Kendall's tau, Spearman's rho, inversions, regret).
          5. Guarantees workspace teardown and automatic retries on transient errors.

    Parameters:
        task_args (Dict[str, Any]): Task configuration dictionary containing scenario parameters,
            candidate profiles, seeds, and paths.

    Returns:
        Dict[str, Any]: Trial outcome mapping including status, records, candidate_records, and duration.

    Author:
        Lukas von Erdmannsdorff
    """
    sc_id = task_args["sc_id"]
    sc_idx = task_args["sc_idx"]
    iter_idx = task_args["iter_idx"]
    suite = task_args["suite"]
    N = task_args["N"]
    n = task_args["n"]
    noise = task_args["noise"]
    delta = task_args["delta"]
    calibrated_median_d = task_args.get("calibrated_median_d", 0.50)
    eff_name = task_args["eff_name"]
    corr = task_args["corr"]
    true_means = task_args["true_means"]
    true_order = task_args["true_order"]
    seed = task_args["seed"]
    temp_root = Path(task_args["temp_root"])
    repo_root = Path(task_args["repo_root"])

    # Create dedicated, sandboxed workspace for this specific trial
    unique_id = uuid.uuid4().hex[:8]
    task_ws = temp_root / f"task_{sc_id}_it{iter_idx:03d}_{unique_id}"
    input_ws = task_ws / "input"
    output_ws = task_ws / "output"

    t_start = time.time()
    max_retries = 2
    last_err = None

    try:
        executor = get_worker_executor(repo_root)

        for attempt in range(1, max_retries + 1):
            try:
                # 1. Generate synthetic data using strictly deterministic RNG stream
                rng = np.random.default_rng(seed)
                data_dict = generate_synthetic_data(true_means, n, noise, corr, input_ws, rng)

                # 2. Compute empirical effect size statistics matching +analysis calculate_real_effects
                c1_m1 = data_dict["C1"][:, 0]
                c2_m1 = data_dict["C2"][:, 0]
                c1_m2 = data_dict["C1"][:, 1]
                c2_m2 = data_dict["C2"][:, 1]

                obs_cliffs_m1 = compute_cliffs_delta(c1_m1, c2_m1)
                obs_cliffs_m2 = compute_cliffs_delta(c1_m2, c2_m2)
                obs_delta_m1 = float(np.mean(c1_m1) - np.mean(c2_m1))
                obs_delta_m2 = float(np.mean(c1_m2) - np.mean(c2_m2))

                # Median pairwise Cliff's d across all candidate pairs on primary (M1) and secondary (M2) metrics
                cands = list(data_dict.keys())
                pairwise_d = []
                for i_idx in range(len(cands)):
                    for j_idx in range(i_idx + 1, len(cands)):
                        ci = cands[i_idx]
                        cj = cands[j_idx]
                        pairwise_d.append(abs(compute_cliffs_delta(data_dict[ci][:, 0], data_dict[cj][:, 0])))
                        pairwise_d.append(abs(compute_cliffs_delta(data_dict[ci][:, 1], data_dict[cj][:, 1])))
                obs_median_cliffs_d = float(np.median(pairwise_d)) if pairwise_d else 0.0

                # 3. Run Baselines
                topsis_ranks = run_topsis_baseline(data_dict)
                borda_ranks = run_borda_baseline(data_dict)
                wilcoxon_copeland_ranks = run_wilcoxon_copeland_baseline(data_dict)

                # 4. Run HERA Ranking with deterministic seed propagation
                hera_ranks = run_hera_ranking(executor, input_ws, output_ws, seed=seed)

                # 5. Evaluate all methods against ground truth
                records = []
                candidate_records = []
                for method_name, pred_ranks in [
                    ("HERA", hera_ranks),
                    ("TOPSIS", topsis_ranks),
                    ("Borda", borda_ranks),
                    ("Wilcoxon-Copeland", wilcoxon_copeland_ranks)
                ]:
                    metrics = evaluate_ranking(pred_ranks, true_order, true_means, suite=suite)
                    records.append({
                        "ScenarioID": sc_id,
                        "Suite": suite,
                        "Candidates": N,
                        "SampleSize": n,
                        "Noise": noise,
                        "EffectMagnitude": eff_name,
                        "Delta": delta,
                        "Calibrated_Median_Cliffs_d": calibrated_median_d,
                        "Correlation": corr,
                        "Iteration": iter_idx - 1,
                        "Method": method_name,
                        **metrics,
                        "Observed_Cliffs_d_M1": round(obs_cliffs_m1, 4),
                        "Observed_Cliffs_d_M2": round(obs_cliffs_m2, 4),
                        "Observed_Delta_M1": round(obs_delta_m1, 4),
                        "Observed_Delta_M2": round(obs_delta_m2, 4),
                        "Observed_Median_Cliffs_d": round(obs_median_cliffs_d, 4)
                    })
                    for c_idx, c in enumerate(true_order):
                        candidate_records.append({
                            "ScenarioID": sc_id,
                            "Suite": suite,
                            "Candidates": N,
                            "SampleSize": n,
                            "Noise": noise,
                            "EffectMagnitude": eff_name,
                            "Delta": delta,
                            "Calibrated_Median_Cliffs_d": calibrated_median_d,
                            "Correlation": corr,
                            "Iteration": iter_idx - 1,
                            "Method": method_name,
                            "Candidate": c,
                            "TrueRank": c_idx + 1,
                            "PredictedRank": int(pred_ranks[c])
                        })

                t_dur = time.time() - t_start
                return {
                    "status": "success",
                    "sc_id": sc_id,
                    "sc_idx": sc_idx,
                    "iter_idx": iter_idx,
                    "records": records,
                    "candidate_records": candidate_records,
                    "duration": t_dur
                }
            except Exception as e:
                last_err = e
                # Clean up workspace before retrying
                if task_ws.exists():
                    shutil.rmtree(task_ws, ignore_errors=True)
                if attempt < max_retries:
                    time.sleep(1.0)

        # If all retries failed:
        return {
            "status": "failed",
            "sc_id": sc_id,
            "sc_idx": sc_idx,
            "iter_idx": iter_idx,
            "error": str(last_err),
            "records": [],
            "candidate_records": [],
            "duration": time.time() - t_start
        }
    finally:
        # Guarantee instant cleanup of worker disk artifacts
        if task_ws.exists():
            shutil.rmtree(task_ws, ignore_errors=True)


def run_sliding_window_pipeline(
    scenarios: List[Dict[str, Any]],
    iterations: int,
    base_seed: int,
    num_workers: int,
    run_temp_dir: Path,
    repo_root: Path,
    logger: logging.Logger,
    results_csv: Path,
    candidate_ranks_csv: Path,
    summary_csv: Path,
    checkpoint_json: Path,
    get_ground_truth_means_fn: Any
) -> Tuple[List[Dict[str, Any]], float, bool]:
    """
    Executes the high-performance non-blocking sliding-window pipeline.

    Syntax:
        results, elapsed, interrupted = run_sliding_window_pipeline(
            scenarios, iterations, base_seed, num_workers, run_temp_dir,
            repo_root, logger, results_csv, candidate_ranks_csv,
            summary_csv, checkpoint_json, get_ground_truth_means_fn
        )

    Description:
        Schedules parallel Monte Carlo trials across concurrent worker processes
        using a bounded sliding-window buffer to ensure 100% core saturation without
        RAM ballooning. Continuously streams completed results to disk, records
        incremental state checkpoints, updates the global summary file, and prints
        real-time progress milestones.

    Parameters:
        scenarios (List[Dict[str, Any]]): List of scenario configuration dictionaries.
        iterations (int): Number of iterations per scenario.
        base_seed (int): Base PRNG seed.
        num_workers (int): Number of parallel worker processes.
        run_temp_dir (Path): Scratch directory for task workspaces.
        repo_root (Path): Root directory of the repository.
        logger (logging.Logger): Logger instance.
        results_csv (Path): Output path for streamed raw results CSV.
        candidate_ranks_csv (Path): Output path for candidate ranks CSV.
        summary_csv (Path): Output path for global aggregated summary CSV.
        checkpoint_json (Path): Output path for runtime state checkpoint JSON.
        get_ground_truth_means_fn (Any): Function returning ground-truth means and candidate order.

    Returns:
        Tuple[List[Dict[str, Any]], float, bool]:
            Tuple of (all_results_cache, total_elapsed_seconds, was_interrupted_flag).

    Author:
        Lukas von Erdmannsdorff
    """
    total_runs = len(scenarios) * iterations
    total_scenarios = len(scenarios)
    logger.info(f"[Execution] Commencing Monte Carlo benchmark ({total_runs} total runs across {total_scenarios} scenarios)...")
    logger.info("---------------------------------------------------------------------------------------------------------")
    t_start = time.time()

    all_results_cache: List[Dict[str, Any]] = []
    is_tty = sys.stdout.isatty()

    if iterations >= 50:
        log_step = 10
    elif iterations >= 20:
        log_step = 5
    elif iterations >= 6:
        log_step = 2
    else:
        log_step = 1

    interrupted = False

    def sig_handler(sig, frame):
        """Handles external SIGINT/SIGTERM termination signals gracefully."""
        nonlocal interrupted
        interrupted = True
        logger.warning("\n[Interrupt] Received stop signal! Terminating worker pool and securing partial data...")
    original_sigint = signal.signal(signal.SIGINT, sig_handler)
    original_sigterm = signal.signal(signal.SIGTERM, sig_handler)

    # 1. Initialize task generator across all scenarios
    def task_generator():
        """Yields parameterized Monte Carlo trial task dictionaries across all experimental scenarios."""
        for sc_idx, sc in enumerate(scenarios, 1):
            sc_id = sc["ScenarioID"]
            sc_num_id = sc["NumericID"]
            N = sc["Candidates"]
            n = sc["SampleSize"]
            noise = sc["Noise"]
            delta = sc["Delta"]
            eff_name = sc["EffectMagnitude"]
            corr = sc["Correlation"]
            suite = sc["Suite"]
            true_means, true_order = get_ground_truth_means_fn(suite, N, delta)

            for iter_idx in range(1, iterations + 1):
                iter_seed = int(base_seed + (sc_num_id * 10000) + iter_idx)
                yield {
                    "sc_id": sc_id,
                    "sc_idx": sc_idx,
                    "iter_idx": iter_idx,
                    "suite": suite,
                    "N": N,
                    "n": n,
                    "noise": noise,
                    "delta": delta,
                    "calibrated_median_d": sc.get("Calibrated_Median_Cliffs_d", 0.50),
                    "eff_name": eff_name,
                    "corr": corr,
                    "true_means": true_means,
                    "true_order": true_order,
                    "seed": iter_seed,
                    "temp_root": str(run_temp_dir.resolve()),
                    "repo_root": str(repo_root.resolve()),
                }

    # Tracking per-scenario status
    scenario_tracker = {}
    for sc_idx, sc in enumerate(scenarios, 1):
        scenario_tracker[sc["ScenarioID"]] = {
            "sc_idx": sc_idx,
            "sc": sc,
            "completed": 0,
            "total": iterations,
            "t_start": None,
            "t_end": None,
            "logged_start": False
        }

    completed_global_iters = 0
    max_in_flight = max(2, num_workers * 2)

    try:
        mp_context = mp.get_context("spawn")
        gen = task_generator()
        active_futures: Dict[concurrent.futures.Future, Dict[str, Any]] = {}
        recent_durations: deque = deque(maxlen=40)

        with ProcessPoolExecutor(max_workers=num_workers, mp_context=mp_context) as pool:
            # Initial filling of the sliding-window pipeline
            for task in gen:
                sc_id = task["sc_id"]
                track = scenario_tracker[sc_id]
                if track["t_start"] is None:
                    track["t_start"] = time.time()
                
                # Only announce Scenario 1 on initialization; subsequent scenarios are announced upon predecessor completion
                if track["sc_idx"] == 1 and not track["logged_start"]:
                    track["logged_start"] = True
                    log_scenario_header(track, total_scenarios, logger, is_tty)

                fut = pool.submit(execute_single_iteration, task)
                active_futures[fut] = task
                if len(active_futures) >= max_in_flight:
                    break

            # Process tasks continuously as they complete (Zero Worker Idle Time)
            while active_futures:
                if interrupted:
                    for f in list(active_futures.keys()):
                        f.cancel()
                    break

                done, _ = concurrent.futures.wait(
                    active_futures.keys(), 
                    return_when=concurrent.futures.FIRST_COMPLETED
                )

                for fut in done:
                    task = active_futures.pop(fut)
                    if interrupted:
                        continue

                    try:
                        res = fut.result()
                    except Exception as e:
                        res = {
                            "status": "failed",
                            "sc_id": task["sc_id"],
                            "sc_idx": task["sc_idx"],
                            "iter_idx": task["iter_idx"],
                            "error": str(e),
                            "records": [],
                            "candidate_records": [],
                            "duration": 0.0
                        }

                    sc_id = res["sc_id"]
                    track = scenario_tracker[sc_id]
                    track["completed"] += 1
                    completed_global_iters += 1

                    if res["status"] == "success":
                        append_results_to_csv(results_csv, res["records"], expected_cols=RESULT_COLUMNS)
                        if "candidate_records" in res and res["candidate_records"]:
                            append_results_to_csv(candidate_ranks_csv, res["candidate_records"], expected_cols=CANDIDATE_RANK_COLUMNS)
                        all_results_cache.extend(res["records"])
                        if res.get("duration", 0.0) > 0.0:
                            recent_durations.append(res["duration"])
                    else:
                        logger.error(f"[Worker Error] Trial failed in {res['sc_id']} (Iteration {res['iter_idx']}): {res.get('error')}")

                    # Calculate live timings
                    sc_completed = track["completed"]
                    sc_elapsed_now = time.time() - (track["t_start"] or time.time())
                    pace_now = sc_elapsed_now / max(1, sc_completed)
                    iter_pct = (sc_completed / iterations) * 100
                    sc_rem_sec = pace_now * (iterations - sc_completed)

                    total_elapsed_now = time.time() - t_start
                    total_study_pct = (completed_global_iters / total_runs) * 100
                    remaining_study_iters = total_runs - completed_global_iters

                    # Stable, non-spiking ETA accounting for parallel worker throughput
                    if len(recent_durations) < max(2, min(4, num_workers)):
                        est_total_rem_str = "Calibrating..."
                    else:
                        median_dur = float(np.median(recent_durations))
                        effective_pace_sec = median_dur / max(1, num_workers)
                        est_total_rem_sec = effective_pace_sec * remaining_study_iters
                        est_total_rem_str = format_time_duration(est_total_rem_sec)

                    if is_tty:
                        clear_terminal_line()
                        sys.stdout.write(
                            f"\r  -> [Overall: {track['sc_idx']:02d}/{total_scenarios:02d} ({total_study_pct:4.1f}%)] {sc_id} | "
                            f"[Sub: {sc_completed:03d}/{iterations:03d} ({iter_pct:4.1f}%)] "
                            f"Pace: {pace_now:4.2f}s/it | Rem: {format_time_duration(sc_rem_sec)} | "
                            f"Overall ETA: {est_total_rem_str}   "
                        )
                        sys.stdout.flush()

                    # Periodic milestone logging
                    if (sc_completed % log_step == 0) and (sc_completed < iterations):
                        if is_tty:
                            clear_terminal_line()
                        time_milestone = datetime.now().strftime("%H:%M:%S")
                        logger.info(
                            f"[{time_milestone}] [Overall: {track['sc_idx']:02d}/{total_scenarios:02d} ({total_study_pct:4.1f}%)] {sc_id} | "
                            f"[Sub: {sc_completed:03d}/{iterations:03d} ({iter_pct:4.1f}%)] | "
                            f"Pace: {pace_now:4.2f}s/it | Rem: {format_time_duration(sc_rem_sec)} | "
                            f"Overall ETA: {est_total_rem_str}"
                        )

                    # Scenario completion & incremental state flush
                    if sc_completed == iterations:
                        track["t_end"] = time.time()
                        sc_duration = track["t_end"] - (track["t_start"] or time.time())
                        sc_avg_iter = sc_duration / max(1, iterations)
                        time_str_end = datetime.now().strftime("%H:%M:%S")

                        if all_results_cache:
                            update_global_summary_file(summary_csv, pd.DataFrame(all_results_cache))

                        save_checkpoint(checkpoint_json, {
                            "last_completed_scenario": sc_id,
                            "completed_scenarios_count": track["sc_idx"],
                            "total_scenarios": total_scenarios,
                            "completed_iterations": completed_global_iters,
                            "total_iterations": total_runs,
                            "elapsed_seconds": total_elapsed_now,
                            "timestamp": datetime.now().isoformat()
                        })

                        if is_tty:
                            clear_terminal_line()
                        logger.info(
                            f"[{time_str_end}] [Overall: {track['sc_idx']:02d}/{total_scenarios:02d} ({total_study_pct:4.1f}%)] {sc_id} "
                            f"-> Completed in {format_time_duration(sc_duration)} ({sc_avg_iter:4.2f}s/it) | "
                            f"Overall ETA: {est_total_rem_str}"
                        )

                        # Announce subsequent scenario in strict sequential order
                        next_sc_idx = track["sc_idx"] + 1
                        for _, n_track in scenario_tracker.items():
                            if n_track["sc_idx"] == next_sc_idx and not n_track["logged_start"]:
                                n_track["logged_start"] = True
                                if n_track["t_start"] is None:
                                    n_track["t_start"] = time.time()
                                log_scenario_header(n_track, total_scenarios, logger, is_tty)
                                break

                    # Replenish pipeline with the next task from the generator
                    if not interrupted:
                        try:
                            next_task = next(gen)
                            next_sc_id = next_task["sc_id"]
                            next_track = scenario_tracker[next_sc_id]
                            if next_track["t_start"] is None:
                                next_track["t_start"] = time.time()
                            new_fut = pool.submit(execute_single_iteration, next_task)
                            active_futures[new_fut] = next_task
                        except StopIteration:
                            pass

    except KeyboardInterrupt:
        logger.warning("\n[Interrupt] KeyboardInterrupt detected. Securing collected data...")
        interrupted = True
    finally:
        signal.signal(signal.SIGINT, original_sigint)
        signal.signal(signal.SIGTERM, original_sigterm)

    elapsed_total = time.time() - t_start
    return all_results_cache, elapsed_total, interrupted
