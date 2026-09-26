"""
HERA Simulation Engine Sub-package.

Provides modular components for Monte Carlo ground-truth validation benchmarking:
- config: Experimental grid constants, parameter calibrations, and schemas
- system: Hardware profiling, DRAS worker scaling, logging, and ASCII banners
- engine: HERA computational wrapper across PyPI and local MATLAB CLI
- synthesis: Ground-truth generation, multivariate synthetic data, and Cliff's d
- mcdm: Baseline decision methods (TOPSIS, Borda, Wilcoxon-Copeland) & ranking evaluation
- reporting: Streaming CSV writers, checkpoints, FAIR metadata, and 300 DPI publication tables
- simulate: Sandboxed worker iteration logic and sliding-window parallel pipeline
- plot_results: Complete publication visualization suite

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

__version__ = "1.4.7"

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
    compute_superordinate_summary, compute_detailed_pooled_statistics, export_superordinate_csvs,
    update_global_summary_file, save_checkpoint,
    save_scenario_configuration, save_method_configuration, save_candidate_configuration
)

from .plot_results import (
    generate_all_plots, generate_global_summary_pdf,
    plot_pooled_core_distributions, plot_pooled_core_marginal_scales,
    plot_pooled_compensatory_stress, plot_pooled_sensitivity_summary
)

from .simulate import (
    execute_single_iteration, run_sliding_window_pipeline
)

from .run_simulation import main

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
    "generate_all_plots", "generate_global_summary_pdf",
    "plot_pooled_core_distributions", "plot_pooled_core_marginal_scales",
    "plot_pooled_compensatory_stress", "plot_pooled_sensitivity_summary",
    "execute_single_iteration", "run_sliding_window_pipeline",
    "main"
]
