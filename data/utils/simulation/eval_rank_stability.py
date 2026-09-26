"""
Multi-Scale Candidate Rank Stability Benchmark Suite.

Syntax:
    python3 eval_rank_stability.py

Description:
    Generates multi-scale candidate ranking data across N = 10, 12, 14 candidates
    pooled across noise levels to evaluate middle-field ordinal stability.

Author:
    Lukas von Erdmannsdorff
"""

import sys
import time
from pathlib import Path
from typing import Tuple, List, Dict, Any
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import pandas as pd

BASE_DIR = Path(__file__).parent.resolve()
if BASE_DIR.name == "simulation":
    UTILS_DIR = BASE_DIR.parent
    DATA_DIR = UTILS_DIR.parent
else:
    UTILS_DIR = BASE_DIR
    DATA_DIR = BASE_DIR.parent

if str(UTILS_DIR) not in sys.path:
    sys.path.insert(0, str(UTILS_DIR))
if str(BASE_DIR) not in sys.path:
    sys.path.insert(0, str(BASE_DIR))

try:
    from .config import CANDIDATE_RANK_COLUMNS
    from .engine import HERAExecutor
    from .synthesis import get_ground_truth_means, generate_synthetic_data
    from .mcdm import (
        run_hera_ranking, run_topsis_baseline, run_borda_baseline,
        run_wilcoxon_copeland_baseline
    )
    from .reporting import append_results_to_csv
except ImportError:
    from simulation.config import CANDIDATE_RANK_COLUMNS
    from simulation.engine import HERAExecutor
    from simulation.synthesis import get_ground_truth_means, generate_synthetic_data
    from simulation.mcdm import (
        run_hera_ranking, run_topsis_baseline, run_borda_baseline,
        run_wilcoxon_copeland_baseline
    )
    from simulation.reporting import append_results_to_csv


def evaluate_single_trial(args: Tuple[int, float, int, int]) -> List[Dict[str, Any]]:
    """
    Executes a single stability benchmark trial for candidate-level rank analysis.

    Syntax:
        records = evaluate_single_trial(args)

    Description:
        Runs an isolated multi-metric trial at candidate scale N, noise level, and
        iteration index, extracting predicted ranks from HERA and baseline MCDM methods.

    Parameters:
        args (Tuple[int, float, int, int]): Tuple containing (N, noise, iter_idx, seed).

    Returns:
        List[Dict[str, Any]]: List of candidate-level rank record dictionaries.

    Author:
        Lukas von Erdmannsdorff
    """
    N, noise, iter_idx, seed = args
    rng = np.random.default_rng(seed)
    
    # Isolate temporary directory per trial with unique UUID
    import uuid
    trial_id = f"trial_N{N}_s{int(noise*10)}_it{iter_idx}_{uuid.uuid4().hex[:6]}"
    temp_dir = BASE_DIR / "temp_stability_eval" / trial_id
    in_dir = temp_dir / "input"
    out_dir = temp_dir / "output"
    
    max_retries = 3
    last_err = None

    for attempt in range(1, max_retries + 1):
        try:
            if temp_dir.exists():
                import shutil
                shutil.rmtree(temp_dir, ignore_errors=True)
            in_dir.mkdir(parents=True, exist_ok=True)
            out_dir.mkdir(parents=True, exist_ok=True)

            # Stagger slightly to avoid MATLAB launch contention
            time.sleep(np.random.uniform(0.5, 2.5))

            executor = HERAExecutor()
            means, true_order = get_ground_truth_means("Core", N, delta=5.0)
            data = generate_synthetic_data(means, n_samples=50, noise_level=noise, correlation=0.0, output_dir=in_dir, rng=rng)
            
            hera_r = run_hera_ranking(executor, in_dir, out_dir)
            top_r = run_topsis_baseline(data)
            borda_r = run_borda_baseline(data)
            wilc_r = run_wilcoxon_copeland_baseline(data)
            
            methods = {
                "HERA": hera_r,
                "TOPSIS": top_r,
                "Borda": borda_r,
                "Wilcoxon-Copeland": wilc_r
            }
            
            records = []
            for m_name, pred_ranks in methods.items():
                for c_idx, c in enumerate(true_order):
                    records.append({
                        "ScenarioID": f"STAB_N{N}",
                        "Suite": "Core",
                        "Candidates": N,
                        "SampleSize": 50,
                        "Noise": noise,
                        "EffectMagnitude": "Medium",
                        "Correlation": 0.0,
                        "Iteration": iter_idx,
                        "Method": m_name,
                        "Candidate": c,
                        "TrueRank": c_idx + 1,
                        "PredictedRank": pred_ranks.get(c, N)
                    })
            return records
        except Exception as e:
            last_err = e
            time.sleep(2.0 * attempt)
        finally:
            import shutil
            shutil.rmtree(temp_dir, ignore_errors=True)
            
    print(f"    [Warning] Trial N={N}, noise={noise} failed after {max_retries} attempts: {last_err}")
    return []

def main() -> None:
    """
    Executes the multi-scale candidate rank stability benchmark suite.

    Syntax:
        main()

    Description:
        Runs candidate-level rank stability simulations across candidate scales
        N in {10, 12, 14} and noise levels sigma in {2, 4, 6, 8, 10}%, aggregating
        predicted ranks across MCDM methods and re-rendering publication graphics.

    Parameters:
        None

    Returns:
        None

    Author:
        Lukas von Erdmannsdorff
    """
    candidate_scales = [10, 12, 14]
    noise_levels = [2.0, 4.0, 6.0, 8.0, 10.0]
    iterations_per_noise = 2  # Total 2 * 5 = 10 iterations per candidate scale = 30 trials
    
    tasks = []
    base_seed = 42
    task_idx = 0
    for N in candidate_scales:
        for noise in noise_levels:
            for it in range(iterations_per_noise):
                seed = base_seed + task_idx * 17
                tasks.append((N, noise, it, seed))
                task_idx += 1
                
    total_tasks = len(tasks)
    print(f"[Stability Benchmark] Launching {total_tasks} trials across N={candidate_scales} (pooled across sigma in {noise_levels}%)...")
    
    out_csv = DATA_DIR / "Simulation_Output" / "candidate_ranks.csv"
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        pd.DataFrame(columns=CANDIDATE_RANK_COLUMNS).to_csv(f, index=False)
        
    num_workers = 3
    completed = 0
    t0 = time.time()
    
    with ProcessPoolExecutor(max_workers=num_workers) as pool:
        futures = [pool.submit(evaluate_single_trial, t) for t in tasks]
        for f in as_completed(futures):
            res = f.result()
            completed += 1
            if res:
                append_results_to_csv(out_csv, res, expected_cols=CANDIDATE_RANK_COLUMNS)
            elapsed = time.time() - t0
            pace = elapsed / completed
            rem = pace * (total_tasks - completed)
            print(f"[{completed:02d}/{total_tasks:02d}] Completed trial in {pace:.1f}s/trial | Est. remaining: {rem:.0f}s", flush=True)

    print(f"[Done] Multi-scale candidate rank data saved to: {out_csv}")
    
    # Re-render Figure 10
    try:
        from . import plot_results
    except ImportError:
        import plot_results
    plots_dir = DATA_DIR / "Simulation_Output" / "plots"
    plot_results.plot_candidate_rank_stability(plots_dir / "middle_field_rank_stability.png", out_csv)
    print(f"[Done] Figure 10 re-rendered with all 3 scales: {plots_dir / 'middle_field_rank_stability.png'}")

if __name__ == "__main__":
    main()
