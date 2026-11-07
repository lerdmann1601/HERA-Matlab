# Directory Structure

This directory contains synthetic datasets (`Example_1`, `Example_2`, and `Example_2_workflow`) used for testing and validating the HERA toolbox. These datasets are based on characteristics of real MRI images. The original MRI datasets are not included in this repository.

In addition, it includes clinical, weather, and Large Language Model (LLM) evaluation datasets (`Example_3`, `Example_4`, and `Example_5`) generated to demonstrate HERA's versatility across different domains, as well as the Monte Carlo ground-truth validation suite (`utils/`) and its reference outputs (`results/Benchmark Simulation/`).

For Examples 3 and 4, the datasets were created by training 10 different Machine Learning classifiers on real-world data: the clinical 'Cardiovascular Disease' dataset (accessed via Kaggle at <https://www.kaggle.com/datasets/sulianova/cardiovascular-disease-dataset>) for Example 3, and the 'Rain in Australia' (WeatherAUS) dataset—originally provided by the Australian Bureau of Meteorology (accessed via Kaggle at <https://www.kaggle.com/datasets/jsphyg/weather-dataset-rattle-package>) for Example 4. Prior to training, leaky variables were removed. The raw features were passed through an automated preprocessing pipeline: missing numerical values were imputed using the median and standardized, while categorical variables were imputed with the most frequent value and one-hot encoded.

Example 5 introduces a modern LLM stability analysis. It utilizes raw, instance-level prediction traces from the MMLU-Pro benchmark (Wang et al., 2024; arXiv:2406.01574) for 12 state-of-the-art open-weights foundation models (such as Qwen 2.5, Llama 3.3, and Mistral Large), sourced directly from the Hugging Face Open LLM Leaderboard V2 (Fourrier et al., 2024; URL: <https://huggingface.co/spaces/open-llm-leaderboard/open_llm_leaderboard>). These specific open-weights models were selected because they provide the transparent, instance-level result files (JSONL), a requirement not met by proprietary models whose individual prediction traces are kept private.

To create the Example 3 and 4 data and traditional analysis results for comparison against HERA, each dataset was partitioned into disjoint training and test sets. First, exactly 15,000 instances were sampled for the training set using stratified sampling to preserve the target class distribution. Next, exactly 35,000 completely disjoint test instances were sampled and partitioned into 50 independent test blocks, yielding exactly 700 instances per block. These blocks were generated using stratified sampling (Stratified K-Fold) to guarantee representative class distributions in every block. The resulting performance metrics for each block were exported into CSV files.

Furthermore, traditional statistical evaluations were applied to these metric blocks to establish a methodological baseline for comparison against analysis with HERA. Specifically, the non-parametric frequentist Friedman test paired with the Nemenyi post-hoc test (Demšar, 2006; URL: <https://jmlr.org/papers/v7/demsar06a.html>) was used for the clinical data. For the weather data, Bayesian Wilcoxon signed-rank tests (Benavoli et al., 2017; URL: <https://jmlr.org/papers/v18/16-305.html>) were applied, utilizing a Region of Practical Equivalence (ROPE) of 0.01 (1%). These traditional statistical analyses were complemented by a Borda Count aggregation on the raw performance metrics to generate a final multi-criteria point-estimate ranking for both datasets in comparison to the HERA-generated rankings.

For Example 5, the raw JSONL result files were fetched from the official Hugging Face evaluation repositories. This complete dataset was partitioned into 30 disjoint, independent test blocks using stratified k-fold splitting based on subject categories (e.g., mathematics, physics, law) to ensure representative task distributions in every block. The mean accuracy for each block was then exported into CSV files to evaluate with HERA against the strict, single-value official MMLU-Pro scores reported by the Open LLM Leaderboard (ignoring other benchmarks like MATH or GPQA).

*Note on Scientific Context:* In early 2025, Hugging Face officially retired and archived the centralized Open LLM Leaderboard V2—halting the continuous publication of standardized, instance-level evaluation traces—as static benchmarks increasingly suffered from rapid saturation, data contamination, and targeted over-optimization. While modern alternatives increasingly shift toward decentralized community evaluations (`OpenEvals`), dynamic benchmarks, or private/held-out test splits to prevent cheating, private test sets sacrifice the public instance-level transparency required for rigorous paired statistical testing. Thus, we explicitly utilize the archived V2 dataset, as it represents the last standardized, fully transparent snapshot of top open-weights models at peak benchmark saturation. Even if specific models benefited from contamination, this only exacerbates score clustering at the top, making V2 the ideal real-world case study. (For Hugging Face's official discussion of how benchmark saturation and data contamination drove their leaderboard overhauls, see their [Open LLM Leaderboard v2 Launch Blog Post](https://huggingface.co/spaces/open-llm-leaderboard/blog)).

```text
data/
├── examples/                  % Synthetic and empirical example datasets
│   ├── Example_1/             % Data for Example 1 (MRI synthetic metrics)
│   │   ├── CNR.csv            % Contrast-to-Noise Ratio data
│   │   ├── OC.csv             % Optical Contrast data
│   │   └── SNR.csv            % Signal-to-Noise Ratio data
│   ├── Example_2/             % Data for Example 2 (MRI synthetic metrics)
│   │   ├── CNR.csv            % Contrast-to-Noise Ratio data
│   │   ├── OC.csv             % Optical Contrast data
│   │   └── SNR.csv            % Signal-to-Noise Ratio data
│   ├── Example_2_workflow/    % Input datasets and baseline results for stability analysis
│   │   ├── data/              % Input data folders across alpha levels (Alpha_0 to Alpha_50)
│   │   └── results/           % Reference results generated by stability analysis
│   ├── Example_3/             % Clinical data generated from Cardiovascular Disease dataset
│   │   └── Clinical_Data/     % Data blocks for Example 3 (AUC, Sensitivity, Specificity)
│   ├── Example_4/             % Weather data generated from WeatherAUS dataset
│   │   └── Weather_Data/      % Data blocks for Example 4 (F1, Precision, Recall)
│   └── Example_5/             % LLM stability data based on MMLU-Pro benchmark
│       └── MMLU_Data/         % Data blocks for Example 5 (Accuracy)
├── results/                   % Analysis outputs and benchmark comparison baselines
│   ├── Benchmark Simulation/  % Reference outputs & plots of the Monte Carlo validation study
│   ├── HERA_Example_1/        % Full HERA analysis output for Example 1
│   ├── HERA_Example_2/        % Full HERA analysis output for Example 2
│   ├── Example_3/             % HERA_Analysis & Friedman/Nemenyi/Borda baselines for Example 3
│   ├── Example_4/             % HERA_Analysis & Bayesian Wilcoxon/Borda baselines for Example 4
│   └── Example_5/             % HERA_Analysis & official leaderboard comparison for Example 5
└── utils/                     % Python workflow, simulation, and plotting utilities
    ├── simulation/            % Encapsulated simulation package
    ├── plot_results.py        % CLI entry point for publication plotting and reporting
    ├── requirements.txt       % Python dependencies for workflow and simulation scripts
    ├── run_simulation.py      % CLI entry point for ground-truth validation simulations
    └── run_workflow.py        % Stability analysis workflow runner for Example 2
```

## Usage

You can use the datasets in `examples/` to run HERA. To reproduce the results, adjust the local `folderPath` and `output_dir` paths inside the `configuration.json` files located in `results/HERA_Example_1/`, `results/HERA_Example_2/`, and `results/Example_*/HERA_Analysis/`. Please use enough CPU cores to speed up the analysis.

### Utilities and Simulation Suite

The `utils/` folder provides automated Python tooling for running workflows, benchmarking simulations, and generating publication figures:

- `run_workflow.py`: Runs the automated stability analysis workflow across decreasing contrast enhancement levels (Example 2). Expects the `Alpha_*` input folders inside `utils/data/` (which can be copied from `examples/Example_2_workflow/data/`) and writes outputs to `utils/results/`.
- `run_simulation.py`: Runs Monte Carlo ground-truth validation simulations benchmarking HERA against standard MCDM baselines (TOPSIS, Borda Count, and Copeland Scoring). Reference results are stored in `results/Benchmark Simulation/`.
- `plot_results.py`: Generates publication-ready figures and report summaries from simulation outputs.
- `simulation/`: Modular package containing the simulation engine, rank stability evaluation, MCDM baselines, and data synthesis.

> [!NOTE]
> **Platform Support for Simulations:** While `run_workflow.py` runs cross-platform on any supported operating system, the Monte Carlo simulation suite (`run_simulation.py`) currently runs reliably only on **macOS** (specifically calibrated for Apple Silicon). Spawning dozens of concurrent worker processes that invoke the underlying HERA backend across thousands of trials requires strict per-worker thread isolation (`VECLIB_MAXIMUM_THREADS=1`, `-singleCompThread`) and Dynamic Resource-Aware Scaling (DRAS) to prevent core contention and memory swapping. On Linux and Windows, multi-process concurrency with the compiled engine backend is not yet fully optimized and may encounter thread oversubscription.

To use the Python utilities:

```bash
cd utils
pip install -r requirements.txt
python run_workflow.py
# or: python run_simulation.py
# or: python plot_results.py
```

---

## Example Walkthrough of HERA Results

👉 [HERA Example Analysis](https://lerdmann1601.github.io/HERA-Matlab/Example_Analysis)
