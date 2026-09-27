"""
Simulation Configuration, Experimental Grid & Schema Definitions.

Defines all constants, benchmark scales, effect calibrations, and CSV schemas
matching HERA's +analysis package standards.

Author: Lukas von Erdmannsdorff
"""

import os
import sys
from pathlib import Path
from typing import Dict, List, Any
import numpy as np

# Directory Paths
BASE_DIR = Path(__file__).parent.parent.resolve()
DATA_DIR = BASE_DIR.parent
TEMP_DIR = BASE_DIR / "temp_workspace"

# Ensure local module imports resolve reliably
if str(BASE_DIR) not in sys.path:
    sys.path.insert(0, str(BASE_DIR))

# Synthetic multi-metric evaluation scale (percentages):
# M1 (Primary): 80.0%, M2 (Secondary): 75.0%, M3 (Tertiary): 80.0%
NUM_CANDIDATES_LIST = [10, 12, 14]
DEFAULT_CANDIDATES = 12
SAMPLE_SIZES = [25, 50, 100]
DEFAULT_SAMPLE_SIZE = 50
NOISE_LEVELS = [2.0, 4.0, 6.0, 8.0, 10.0]

# Standardized effect magnitude calibration (percentages & non-parametric Cliff's d):
# Small Effect:  Delta =  2.5% -> Median Cliff's d ~ 0.25 (Calibrated range: 0.20 - 0.30)
# Medium Effect: Delta =  5.0% -> Median Cliff's d ~ 0.50 (Calibrated range: 0.45 - 0.60)
# Large Effect:  Delta =  8.0% -> Median Cliff's d ~ 0.80 (Calibrated range: 0.75 - 0.90)
# Stress / Trap: Delta = 12.0% -> Median Cliff's d ~ 0.95 (Calibrated range: > 0.90)
EFFECT_CALIBRATION = {
    "Small": {
        "delta": 2.5,
        "median_cliffs_d": 0.25,
        "cliffs_d_range": "0.20 - 0.30",
        "label": "Small (d ~ 0.25)"
    },
    "Medium": {
        "delta": 5.0,
        "median_cliffs_d": 0.50,
        "cliffs_d_range": "0.45 - 0.60",
        "label": "Medium (d ~ 0.50)"
    },
    "Large": {
        "delta": 8.0,
        "median_cliffs_d": 0.80,
        "cliffs_d_range": "0.75 - 0.90",
        "label": "Large (d ~ 0.80)"
    },
    "Fixed": {
        "delta": 12.0,
        "median_cliffs_d": 0.95,
        "cliffs_d_range": "> 0.90",
        "label": "Stress / Trap (d ~ 0.95)"
    }
}

EFFECT_MAGNITUDES = {
    "Small": 2.5,
    "Medium": 5.0,
    "Large": 8.0
}
DEFAULT_EFFECT = "Medium"
DEFAULT_NOISE = 4.0

# Correlation between metrics (0.0 = independent, 0.5 = moderate positive correlation)
CORRELATIONS = [0.0, 0.5]
DEFAULT_CORRELATION = 0.0

# Number of Monte Carlo iterations per scenario (configurable via HERA_SIM_ITERATIONS, Default: 300)
ITERATIONS = int(os.environ.get("HERA_SIM_ITERATIONS", 300))

SEED = int(os.environ.get("HERA_SIM_SEED", 123))

# Fixed HERA bootstrap parameters optimized for Monte Carlo simulation efficiency
# Note: Threshold percentile bootstrap (B=2000) provides high precision for null distribution
# and Holm-Bonferroni significance testing. Setting BCa CI and Rank Stability bootstrap counts
# to B=50 saves substantial execution time without altering the linear ranking decisions.
HERA_CONFIG = {
    "manual_B_thr": 2000,
    "manual_B_ci": 50,
    "manual_B_rank": 50,
    "create_reports": False,
    "create_csvs": False,
    "quiet_mode": True,
    "ranking_mode": "M1_M2_M3",
    "run_sensitivity_analysis": False,
    "run_power_analysis": False,
    "reproducible": True
}

# Baseline weights for TOPSIS: [0.4, 0.4, 0.2] reflecting hierarchical criteria prioritization
TOPSIS_WEIGHTS = np.array([0.4, 0.4, 0.2])

# CSV Result Column Schema with Complete Statistical Traceability (Cliff's d & Delta)
RESULT_COLUMNS = [
    "ScenarioID", "Suite", "Candidates", "SampleSize", "Noise",
    "EffectMagnitude", "Delta", "Calibrated_Median_Cliffs_d", "Correlation", "Iteration", "Method",
    "TopChoice", "CompleteRank", "KendallTau", "SpearmanRho",
    "FalseSuperiority", "Regret", "CompensatoryError", "RankDisplacement", "TopChoiceRankRegret",
    "Observed_Cliffs_d_M1", "Observed_Cliffs_d_M2",
    "Observed_Delta_M1", "Observed_Delta_M2",
    "Observed_Median_Cliffs_d"
]

METRIC_COLS = [
    "TopChoice", "CompleteRank", "KendallTau", "SpearmanRho",
    "FalseSuperiority", "Regret", "CompensatoryError", "RankDisplacement", "TopChoiceRankRegret"
]

STATISTICAL_EFFECT_COLS = [
    "Observed_Cliffs_d_M1", "Observed_Cliffs_d_M2",
    "Observed_Delta_M1", "Observed_Delta_M2",
    "Observed_Median_Cliffs_d"
]

# Candidate-Level Ranks Schema for Ordinal Non-Parametric Stability Analysis
CANDIDATE_RANK_COLUMNS = [
    "ScenarioID", "Suite", "Candidates", "SampleSize", "Noise",
    "EffectMagnitude", "Delta", "Calibrated_Median_Cliffs_d", "Correlation", "Iteration", "Method",
    "Candidate", "TrueRank", "PredictedRank"
]


def make_scenario(
    s_id: int, 
    suite: str, 
    N: int, 
    n: int, 
    noise: float, 
    eff_name: str, 
    eff_val: float, 
    corr: float
) -> Dict[str, Any]:
    """
    Constructs a standardized experimental scenario dictionary.

    Syntax:
        scenario = make_scenario(s_id, suite, N, n, noise, eff_name, eff_val, corr)

    Description:
        Assembles parameter definitions and effect calibration metadata into
        a standardized configuration mapping for an individual benchmark condition.

    Parameters:
        s_id (int): Numeric scenario index identifier.
        suite (str): Benchmark suite name ('Core', 'Sample Size Sensitivity', etc.).
        N (int): Number of candidate AI models evaluated.
        n (int): Evaluation sample size.
        noise (float): Gaussian measurement noise level (sigma, %).
        eff_name (str): Effect size calibration label ('Small', 'Medium', 'Large').
        eff_val (float): Between-candidate performance delta (%).
        corr (float): Inter-metric collinearity correlation coefficient (rho).

    Returns:
        Dict[str, Any]: Structured scenario parameter dictionary.

    Author:
        Lukas von Erdmannsdorff
    """
    calib = EFFECT_CALIBRATION.get(eff_name, {"median_cliffs_d": 0.50, "cliffs_d_range": "0.45 - 0.60"})
    return {
        "ScenarioID": f"SC_{s_id:03d}",
        "NumericID": s_id,
        "Suite": suite,
        "Candidates": N,
        "SampleSize": n,
        "Noise": noise,
        "EffectMagnitude": eff_name,
        "Delta": eff_val,
        "Calibrated_Median_Cliffs_d": calib["median_cliffs_d"],
        "Calibrated_Cliffs_d_Range": calib["cliffs_d_range"],
        "Correlation": corr
    }


# Suite Name Constants
SUITE_CORE = "Core"
SUITE_SAMPLE_SIZE = "Sample Size Sensitivity"
SUITE_EFFECT = "Effect Magnitude Sensitivity"
SUITE_CORRELATION = "Correlation Sensitivity"


def build_scenarios_grid() -> List[Dict[str, Any]]:
    """
    Builds the streamlined 24-condition orthogonal benchmarking scenario grid.

    Syntax:
        scenarios = build_scenarios_grid()

    Description:
        Constructs the complete orthogonal parameter grid spanning four benchmark suites:
          1. Core Sweep (N in {10, 12, 14} x sigma in {2, 4, 6, 8, 10}% at n=50, Delta=5.0%, rho=0.0): 15 scenarios
          2. Sample Size Sensitivity (N in {10, 12, 14} x n in {25, 100} at sigma=4.0%, Delta=5.0%, rho=0.0): 6 scenarios
          3. Effect Magnitude Sensitivity (Delta in {Small (2.5%), Large (8.0%)} at N=12, n=50, sigma=4.0%, rho=0.0): 2 scenarios
          4. Correlation Sensitivity (rho = 0.5 at N=12, n=50, sigma=4.0%, Delta=5.0%): 1 scenario

    Returns:
        List[Dict[str, Any]]: List of 24 scenario configuration dictionaries.

    Author:
        Lukas von Erdmannsdorff
    """
    scenarios = []
    scen_id = 1

    # 1. Core Sweep: Noise x Candidate Scale at standard baseline sample size (n=50)
    for N in NUM_CANDIDATES_LIST:
        for noise in NOISE_LEVELS:
            scenarios.append(make_scenario(
                scen_id, SUITE_CORE, N, DEFAULT_SAMPLE_SIZE, noise, 
                DEFAULT_EFFECT, EFFECT_MAGNITUDES[DEFAULT_EFFECT], DEFAULT_CORRELATION
            ))
            scen_id += 1

    # 2. Sample Size Sensitivity Suite (n in {25, 100} across all N in {10, 12, 14} at Noise=4.0%, Corr=0.0)
    # Note: n=50 is already captured in the Core Sweep for each N at Noise=4.0%
    for N in NUM_CANDIDATES_LIST:
        for n in SAMPLE_SIZES:
            if n == DEFAULT_SAMPLE_SIZE:
                continue
            scenarios.append(make_scenario(
                scen_id, SUITE_SAMPLE_SIZE, N, n, DEFAULT_NOISE, 
                DEFAULT_EFFECT, EFFECT_MAGNITUDES[DEFAULT_EFFECT], DEFAULT_CORRELATION
            ))
            scen_id += 1

    # 3. Effect Magnitude Sensitivity Suite (N=12, n=50, Noise=4.0, Corr=0.0)
    # Note: Medium effect (5.0%) is already captured in the Core Sweep
    for eff_name, eff_val in EFFECT_MAGNITUDES.items():
        if eff_name == DEFAULT_EFFECT:
            continue
        scenarios.append(make_scenario(
            scen_id, SUITE_EFFECT, DEFAULT_CANDIDATES, DEFAULT_SAMPLE_SIZE, DEFAULT_NOISE, 
            eff_name, eff_val, DEFAULT_CORRELATION
        ))
        scen_id += 1

    # 4. Correlation Sensitivity Suite (N=12, n=50, Noise=4.0, Effect=Medium)
    # Note: rho=0.0 is already captured in the Core Sweep
    for corr in CORRELATIONS:
        if corr == DEFAULT_CORRELATION:
            continue
        scenarios.append(make_scenario(
            scen_id, SUITE_CORRELATION, DEFAULT_CANDIDATES, DEFAULT_SAMPLE_SIZE, DEFAULT_NOISE, 
            DEFAULT_EFFECT, EFFECT_MAGNITUDES[DEFAULT_EFFECT], corr
        ))
        scen_id += 1

    return scenarios


