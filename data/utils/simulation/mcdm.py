"""
MCDM Baseline Algorithms & Multi-Criteria Ranking Evaluation.

Implements TOPSIS, Borda Count, Wilcoxon-Copeland tournament with Step-Down Holm FWER,
HERA ranking invocation, and quantitative rank evaluation metrics.

Author: Lukas von Erdmannsdorff
"""

import json
import shutil
from pathlib import Path
from typing import Dict, List, Optional
import numpy as np
from scipy import stats

from .config import HERA_CONFIG, TOPSIS_WEIGHTS
from .engine import HERAExecutor


def run_hera_ranking(
    executor: HERAExecutor, 
    input_dir: Path, 
    output_dir: Path, 
    seed: Optional[int] = None
) -> Dict[str, int]:
    """
    Executes HERA ranking engine and extracts linear ranking from JSON output.

    Syntax:
        rank_dict = run_hera_ranking(executor, input_dir, output_dir, seed=seed)

    Description:
        Constructs a temporary HERA execution configuration JSON, runs the HERA
        computational core, and extracts the resulting ordinal ranking dictionary.
        Injects deterministic PRNG seed into HERA configuration to ensure 100% reproducible bootstrap.

    Parameters:
        executor (HERAExecutor): Active HERA computational backend wrapper.
        input_dir (Path): Directory containing M1.csv, M2.csv, M3.csv.
        output_dir (Path): Output directory for HERA reports and results.
        seed (Optional[int]): Random seed for HERA internal bootstrap reproducibility.

    Returns:
        Dict[str, int]: Mapping of candidate name to assigned integer rank (1 = top choice).

    Author:
        Lukas von Erdmannsdorff
    """
    if output_dir.exists():
        shutil.rmtree(output_dir, ignore_errors=True)
    output_dir.mkdir(parents=True, exist_ok=True)

    config_path = output_dir / "config.json"
    user_input = {
        "folderPath": str(input_dir.resolve()),
        "metric_names": ["M1", "M2", "M3"],
        "output_dir": str(output_dir.resolve()),
        "fileType": ".csv",
        "ranking_mode": HERA_CONFIG["ranking_mode"],
        "create_reports": HERA_CONFIG["create_reports"],
        "create_csvs": HERA_CONFIG["create_csvs"],
        "quiet_mode": HERA_CONFIG["quiet_mode"],
        "run_sensitivity_analysis": HERA_CONFIG["run_sensitivity_analysis"],
        "run_power_analysis": HERA_CONFIG["run_power_analysis"],
        "manual_B_thr": HERA_CONFIG["manual_B_thr"],
        "manual_B_ci": HERA_CONFIG["manual_B_ci"],
        "manual_B_rank": HERA_CONFIG["manual_B_rank"],
        "reproducible": True
    }
    if seed is not None:
        user_input["seed"] = int(seed)

    config = {"userInput": user_input}
    with open(config_path, "w", encoding="utf-8") as f:
        json.dump(config, f, indent=2)

    # Run ranking via active executor backend
    executor.run(config_path)

    # Locate generated JSON result
    ranking_dirs = sorted([d for d in output_dir.iterdir() if d.is_dir() and d.name.startswith("Ranking_")])
    if not ranking_dirs:
        ranking_dirs = sorted([d for d in output_dir.iterdir() if d.is_dir()])
    
    if not ranking_dirs:
        raise FileNotFoundError(f"HERA analysis directory not found in {output_dir}")

    latest_dir = ranking_dirs[-1]
    output_folder = latest_dir / "Output"
    json_files = list(output_folder.glob("*.json"))
    if not json_files:
        raise FileNotFoundError(f"HERA result JSON not found in {output_folder}")

    with open(json_files[0], "r", encoding="utf-8") as f:
        res = json.load(f)

    names = res["dataset_names"]
    ranks = res["results"]["final_rank"]
    return {n: int(r) for n, r in zip(names, ranks)}


try:
    from pymcdm.methods import TOPSIS as PyMCDMTopsis
    from pymcdm.normalizations import vector_normalization as pymcdm_vector_norm
    _PYMCDM_AVAILABLE = True
except ImportError:
    _PYMCDM_AVAILABLE = False


def run_topsis_baseline(data_dict: Dict[str, np.ndarray]) -> Dict[str, int]:
    """
    Executes standard TOPSIS with vector normalization and user weights.

    Syntax:
        rank_dict = run_topsis_baseline(data_dict)

    Description:
        Implements Technique for Order Preference by Similarity to Ideal Solution (TOPSIS).
        Uses peer-reviewed PyMCDM library (Kizielewicz & Salabun, 2020) when available,
        with seamless mathematical fallback to equivalent vectorized NumPy engine.

    Parameters:
        data_dict (Dict[str, np.ndarray]): Dictionary mapping candidate names to (n_samples, 3) arrays.

    Returns:
        Dict[str, int]: Mapping of candidate name to assigned integer rank (1 = top choice).

    Author:
        Lukas von Erdmannsdorff
    """
    cands = list(data_dict.keys())
    means = np.array([np.mean(data_dict[c], axis=0) for c in cands], dtype=np.float64)  # Shape: (N, 3)

    if _PYMCDM_AVAILABLE:
        try:
            # PyMCDM standard TOPSIS: criteria types = 1 (all benefit metrics)
            topsis_engine = PyMCDMTopsis(normalization_function=pymcdm_vector_norm)
            types = np.ones(means.shape[1], dtype=int)
            pref = topsis_engine(means, TOPSIS_WEIGHTS, types)
            # Higher relative closeness is superior -> rank 1 is highest score
            ranks = stats.rankdata(-pref, method="ordinal")
            return {c: int(r) for c, r in zip(cands, ranks)}
        except Exception:
            pass

    # Fallback to vectorized Hwang & Yoon (1981) formulation
    # Vector normalization
    denom = np.sqrt(np.sum(means ** 2, axis=0))
    denom[denom == 0] = 1e-9
    norm_matrix = means / denom

    # Weighted normalized decision matrix
    weighted_matrix = norm_matrix * TOPSIS_WEIGHTS

    # Ideal positive and negative solutions
    v_plus = np.max(weighted_matrix, axis=0)
    v_minus = np.min(weighted_matrix, axis=0)

    # Euclidean separation measures
    s_plus = np.sqrt(np.sum((weighted_matrix - v_plus) ** 2, axis=1))
    s_minus = np.sqrt(np.sum((weighted_matrix - v_minus) ** 2, axis=1))

    # Relative closeness to ideal solution
    closeness = s_minus / (s_plus + s_minus + 1e-9)

    # Highest closeness receives Rank 1
    ranks = stats.rankdata(-closeness, method="ordinal")
    return {c: int(r) for c, r in zip(cands, ranks)}


def run_borda_baseline(data_dict: Dict[str, np.ndarray]) -> Dict[str, int]:
    """
    Executes Borda Count consensus ranking by aggregating metric-wise ranks.

    Syntax:
        rank_dict = run_borda_baseline(data_dict)

    Description:
        Computes metric-wise ranks for each candidate across M1, M2, and M3,
        sums rank positions, and assigns final ranks (lowest sum = Rank 1).

    Parameters:
        data_dict (Dict[str, np.ndarray]): Dictionary mapping candidate names to (n_samples, 3) arrays.

    Returns:
        Dict[str, int]: Mapping of candidate name to assigned integer rank (1 = top choice).

    Author:
        Lukas von Erdmannsdorff
    """
    cands = list(data_dict.keys())
    means = np.array([np.mean(data_dict[c], axis=0) for c in cands])

    rank_sum = np.zeros(len(cands))
    for m in range(3):
        # Higher mean is better -> rank 1 for highest mean
        r = stats.rankdata(-means[:, m], method="ordinal")
        rank_sum += r

    # Lowest rank sum receives Rank 1
    final_ranks = stats.rankdata(rank_sum, method="ordinal")
    return {c: int(r) for c, r in zip(cands, final_ranks)}


def run_wilcoxon_copeland_baseline(data_dict: Dict[str, np.ndarray], alpha: float = 0.05) -> Dict[str, int]:
    """
    Executes paired Wilcoxon-Copeland tournament ranking with Step-Down Holm-Bonferroni FWER control.

    Syntax:
        rank_dict = run_wilcoxon_copeland_baseline(data_dict, alpha=0.05)

    Description:
        For each candidate pair (c1, c2) and each metric m in {M1, M2, M3}:
          - Evaluates H0: median difference == 0 using paired Wilcoxon signed-rank test.
          - P-values p_(1) <= p_(2) <= p_(3) are evaluated using the step-down Holm-Bonferroni procedure:
              Reject H_(k) if p_(k) < alpha / (m_total - k + 1), for k = 1, ..., m_total.
          - If significant:
              mean(diff) > 0 -> c1 receives +1 win, c2 receives -1 loss
              mean(diff) < 0 -> c1 receives -1 loss, c2 receives +1 win
          - At the first non-significant hypothesis, the step-down stops and subsequent hypotheses
            are retained as non-significant (conservative FWER control).
          - If neither candidate dominates significantly, the match is a neutral tie (0 points).
        
        Candidates are ordinally ranked by net tournament score (descending).

    Parameters:
        data_dict (Dict[str, np.ndarray]): Dictionary mapping candidate names to (n_samples, 3) arrays.
        alpha (float): Family-Wise Error Rate significance threshold (default: 0.05).

    Returns:
        Dict[str, int]: Mapping of candidate name to assigned integer rank (1 = top choice).

    Author:
        Lukas von Erdmannsdorff
    """
    cands = list(data_dict.keys())
    scores = np.zeros(len(cands))

    for i, c1 in enumerate(cands):
        for j in range(i + 1, len(cands)):
            c2 = cands[j]
            p_vals = []
            mean_diffs = []
            
            for m in range(3):
                diff = data_dict[c1][:, m] - data_dict[c2][:, m]
                if np.all(diff == 0):
                    p_vals.append(1.0)
                    mean_diffs.append(0.0)
                    continue
                try:
                    p = float(stats.wilcoxon(diff, alternative="two-sided").pvalue)
                except Exception:
                    p = 1.0
                p_vals.append(p)
                mean_diffs.append(float(np.mean(diff)))
            
            m_total = len(p_vals)
            if m_total == 0:
                continue

            # Step-Down Holm-Bonferroni procedure across the metrics for pair (c1, c2)
            sort_order = np.argsort(p_vals)
            for step_idx, orig_idx in enumerate(sort_order):
                threshold = alpha / (m_total - step_idx)
                if p_vals[orig_idx] < threshold:
                    # Statistically significant difference on metric orig_idx
                    if mean_diffs[orig_idx] > 0:
                        scores[i] += 1
                        scores[j] -= 1
                    elif mean_diffs[orig_idx] < 0:
                        scores[i] -= 1
                        scores[j] += 1
                else:
                    # Step-down ceases at the first non-rejection
                    break

    final_ranks = stats.rankdata(-scores, method="ordinal")
    return {c: int(r) for c, r in zip(cands, final_ranks)}


def evaluate_ranking(
    pred_rank_dict: Dict[str, int], 
    true_order: List[str], 
    true_means: Dict[str, List[float]],
    suite: str = "Core"
) -> Dict[str, float]:
    """
    Computes quantitative ranking evaluation metrics against ground truth.

    Syntax:
        metrics = evaluate_ranking(pred_rank_dict, true_order, true_means, suite='Core')

    Description:
        Quantifies ranking fidelity, top-choice accuracy, pairwise inversions,
        and resistance against unconstrained compensatory trade-offs:
          - TopChoice: 1.0 if true best candidate (C1) is ranked #1, 0.0 otherwise.
          - CompleteRank: 1.0 if predicted ranks match ground truth ranks perfectly.
          - KendallTau: Kendall's rank correlation coefficient tau.
          - SpearmanRho: Spearman's monotonic rank correlation coefficient rho.
          - FalseSuperiority: Pairwise inversion rate (proportion of pairs ranked backwards).
          - Regret: Difference between true M1 of C1 and true M1 of chosen #1 model.
          - CompensatoryError: 1.0 if the compensatory candidate (C_N) is ranked in top half (<= N/2).
          - RankDisplacement: Mean absolute rank displacement (|predicted - true|).
          - TopChoiceRankRegret: True rank position of selected top model minus 1.

    Parameters:
        pred_rank_dict (Dict[str, int]): Predicted integer ranks per candidate.
        true_order (List[str]): Ground-truth ordered candidate identifiers.
        true_means (Dict[str, List[float]]): Ground-truth mean performance vectors.
        suite (str): Active experimental benchmark suite name.

    Returns:
        Dict[str, float]: Dictionary mapping evaluation metric names to calculated scores.

    Author:
        Lukas von Erdmannsdorff
    """
    n_cands = len(true_order)
    pred_ranks = np.array([pred_rank_dict[c] for c in true_order])
    true_ranks = np.arange(1, n_cands + 1)

    # 1. Top-choice recovery
    top_choice = 1.0 if pred_rank_dict[true_order[0]] == 1 else 0.0

    # 2. Complete-rank recovery
    complete_rank = 1.0 if np.array_equal(pred_ranks, true_ranks) else 0.0

    # 3. Kendall Tau
    tau, _ = stats.kendalltau(true_ranks, pred_ranks)
    if np.isnan(tau):
        tau = 0.0

    # 4. Spearman Rho
    rho, _ = stats.spearmanr(true_ranks, pred_ranks)
    if np.isnan(rho):
        rho = 0.0

    # 5. False Superiority Rate (Pairwise Inversions)
    total_pairs = n_cands * (n_cands - 1) // 2
    inversions = 0
    for i in range(n_cands):
        for j in range(i + 1, n_cands):
            if pred_ranks[j] < pred_ranks[i]:
                inversions += 1
    false_superiority = inversions / max(1, total_pairs)

    # 6. Expected Regret on primary metric M1
    chosen_best = true_order[np.argmin(pred_ranks)]
    max_m1 = true_means[true_order[0]][0]
    chosen_m1 = true_means[chosen_best][0]
    regret = float(max_m1 - chosen_m1)

    # 7. Compensatory Error Rate (selection of compensatory candidate C_N in top half <= N/2)
    flawed_cand = f"C{n_cands}"
    flawed_rank = pred_rank_dict.get(flawed_cand, n_cands)
    comp_error = 1.0 if flawed_rank <= (n_cands // 2) else 0.0

    # 8. Mean Absolute Rank Displacement (Displacement from ground truth: average |pred_rank - true_rank|)
    rank_displacement = float(np.mean(np.abs(pred_ranks - true_ranks)))

    # 9. Top-Choice Selection Rank Regret: true rank of candidate chosen as #1 minus 1 (0 if true C1 is picked)
    chosen_best_cand = [c for c, r in pred_rank_dict.items() if r == 1][0]
    top_choice_rank_regret = float(true_order.index(chosen_best_cand))

    return {
        "TopChoice": top_choice,
        "CompleteRank": complete_rank,
        "KendallTau": tau,
        "SpearmanRho": rho,
        "FalseSuperiority": false_superiority,
        "Regret": regret,
        "CompensatoryError": comp_error,
        "RankDisplacement": rank_displacement,
        "TopChoiceRankRegret": top_choice_rank_regret
    }
