"""
Ground Truth Synthesis, Multivariate Data Generation & Non-Parametric Effect Sizes.

Constructs clinical evaluation profiles (Example 3 Cardiovascular Setting),
generates correlated synthetic data, and computes Cliff's Delta effect sizes.

Author: Lukas von Erdmannsdorff
"""

import shutil
from pathlib import Path
from typing import Dict, List, Tuple
import numpy as np
import pandas as pd
from scipy import stats


def compute_cliffs_delta(x: np.ndarray, y: np.ndarray) -> float:
    """
    Calculates Cliff's Delta non-parametric effect size: d = (2*U - n1*n2) / (n1*n2).

    Syntax:
        d = compute_cliffs_delta(x, y)

    Description:
        Cliff's Delta measures stochastic dominance and distribution separation on [-1.0, +1.0].
        Strictly matches the exact formula used in HERA's statistical engine (+HERA/+stats/cliffs_delta.m).

    Parameters:
        x (np.ndarray): Primary sample array (length n1).
        y (np.ndarray): Comparator sample array (length n2).

    Returns:
        float: Calculated Cliff's Delta non-parametric effect size.

    Author:
        Lukas von Erdmannsdorff
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    n1 = len(x)
    n2 = len(y)
    if n1 == 0 or n2 == 0:
        return 0.0
    res = stats.mannwhitneyu(x, y, alternative="two-sided")
    u1 = float(res.statistic)
    return float((2.0 * u1 - (n1 * n2)) / (n1 * n2))


def get_ground_truth_means(
    suite: str, 
    num_candidates: int, 
    delta: float
) -> Tuple[Dict[str, List[float]], List[str]]:
    """
    Generates the unified synthetic Ground Truth candidate profile for multi-criteria ranking benchmark.

    Syntax:
        means, true_order = get_ground_truth_means(suite, num_candidates, delta)

    Description:
        Constructs ground-truth performance distributions across 3 evaluation criteria:
          - M1 (Primary Metric): Primary ranking criterion (Baseline: ~80.0%)
          - M2 (Secondary Metric): Secondary priority criterion (Baseline: ~75.0%)
          - M3 (Tertiary Metric): Tertiary tie-breaking criterion (Baseline: ~80.0%)
          
        Unified Candidate Architecture:
          1. C1 (True Rank 1 - Multi-Criterion Leader):
             Possesses decisive M2 superiority (75.0 + delta + 7.0%), giving it the overall true top rank.
          2. C2 (True Rank 2 - Primary Metric Leader):
             Achieves highest primary metric M1 (80.0%), baseline M2 (75.0 + delta), baseline M3 (80.0%).
          3. C3 (True Rank 3 - Secondary Metric Step Jumper):
             Lower M1 (76.0%), but possesses strong M2 (75.0 + delta), beating both C4 and C5.
          4. C4 (True Rank 4 - Tertiary Metric Tie-Break Winner):
             Has M1=77.5%, M2=75.0%, M3=80.0 + delta. Decisively beats C5 on tertiary metric M3.
          5. C5 (True Rank 5 - Tertiary Metric Tie-Break Baseline):
             Has M1=77.5%, M2=75.0%, M3=80.0%. Neutrally tied with C4 on M1 and M2, falls behind on M3.
          6. C6..C_{N-1} (True Ranks 6..N-1 - Descending Baseline Gradient):
             Monotonically decreasing primary metric M1, baseline M2 (75.0%), baseline M3 (80.0%).
          7. C_N (True Rank N - Compensatory Shortcut Trap Candidate):
             Flawed model: mediocre M1 (75.0%), severe deficit on M2 (60.0%), hyper-inflated M3 (99.0%).
             Tests whether decision-making methods resist false promotion due to unconstrained compensatory trade-offs.

    Parameters:
        suite (str): Benchmark suite name.
        num_candidates (int): Total number of candidate AI models N.
        delta (float): Between-candidate performance delta (%).

    Returns:
        Tuple[Dict[str, List[float]], List[str]]:
            Tuple of (candidate_means_dict, true_rank_ordered_candidate_names).

    Author:
        Lukas von Erdmannsdorff
    """
    base_m1 = 80.0
    base_m2 = 75.0
    base_m3 = 80.0

    means: Dict[str, List[float]] = {
        "C1": [78.5, base_m2 + delta + 7.0, base_m3],
        "C2": [base_m1, base_m2 + delta, base_m3],
        "C3": [76.0, base_m2 + delta, base_m3],
        "C4": [77.5, base_m2, base_m3 + delta],
        "C5": [77.5, base_m2, base_m3],
    }

    # Step size for descending baseline candidates C6 .. C_{N-1}
    num_middle = max(1, num_candidates - 6)
    m1_start = 76.0
    m1_end = 68.0
    step = (m1_start - m1_end) / max(1, num_middle)

    for i in range(6, num_candidates):
        means[f"C{i}"] = [round(m1_start - (i - 6) * step, 2), base_m2, base_m3]

    # C_N is the compensatory trap candidate (strictly Rank N)
    means[f"C{num_candidates}"] = [75.0, 60.0, 99.0]

    true_order = [f"C{i}" for i in range(1, num_candidates + 1)]
    return means, true_order


def generate_synthetic_data(
    means: Dict[str, List[float]], 
    n_samples: int, 
    noise_level: float, 
    correlation: float, 
    output_dir: Path,
    rng: np.random.Generator
) -> Dict[str, np.ndarray]:
    """
    Generates synthetic multivariate data and saves metric CSVs matching HERA requirements.

    Syntax:
        data_dict = generate_synthetic_data(means, n_samples, noise_level, correlation, output_dir, rng)

    Description:
        Constructs correlated multivariate normal score arrays across candidate models and metrics.
        Randomly permutes column ordering to prevent input ordering bias and writes
        M1.csv, M2.csv, and M3.csv to the specified output folder.

    Parameters:
        means (Dict[str, List[float]]): Mean performance vector per candidate.
        n_samples (int): Evaluation cohort size n.
        noise_level (float): Gaussian noise standard deviation sigma (%).
        correlation (float): Inter-metric covariance correlation rho.
        output_dir (Path): Destination filesystem directory for metric CSV files.
        rng (np.random.Generator): Dedicated seeded NumPy Generator instance.

    Returns:
        Dict[str, np.ndarray]: Generated synthetic data array per candidate.

    Author:
        Lukas von Erdmannsdorff
    """
    if output_dir.exists():
        shutil.rmtree(output_dir, ignore_errors=True)
    output_dir.mkdir(parents=True, exist_ok=True)

    candidates = list(means.keys())
    
    # Construct 3x3 covariance matrix across the 3 metrics
    cov = np.full((3, 3), correlation * (noise_level ** 2))
    np.fill_diagonal(cov, noise_level ** 2)

    # Randomly permute candidate order using dedicated RNG to eliminate input ordering bias
    shuffled_candidates = rng.permutation(candidates).tolist()

    # Generate data per candidate: shape (n_samples, 3)
    data_dict = {}
    for cand in shuffled_candidates:
        mu = means[cand]
        data = rng.multivariate_normal(mu, cov, size=n_samples)
        data_dict[cand] = data

    # Write metric files M1.csv, M2.csv, M3.csv for HERA with randomized column order
    subject_ids = [f"Sub_{i+1}" for i in range(n_samples)]
    metric_names = ["M1", "M2", "M3"]

    for m_idx, m_name in enumerate(metric_names):
        df_metric = pd.DataFrame({"Subject_ID": subject_ids})
        for cand in shuffled_candidates:
            df_metric[cand] = data_dict[cand][:, m_idx]
        df_metric.to_csv(output_dir / f"{m_name}.csv", index=False)

    return data_dict
