# Directory Structure

This directory contains synthetic datasets used for testing and validating the HERA toolbox. These datasets are based on characteristics of real MRI Images. The original MRI datasets are not included in this repository.

In addition, it includes clinical, weather, and Large Language Model (LLM) evaluation datasets (Examples 3, 4, and 5) generated to demonstrate HERA's versatility across different domains.

For Examples 3 and 4, the datasets were created by training 10 different Machine Learning classifiers on real-world data: the OpenML 'Cardiovascular Disease' dataset (OpenML ID: 45547) for Example 3, and the 'Rain in Australia' (WeatherAUS, OpenML ID: 40994) dataset for Example 4. Prior to training, leaky variables (like 'RISK_MM' in WeatherAUS) were removed. The raw features were passed through an automated preprocessing pipeline: missing numerical values were imputed using the median and standardized, while categorical variables were imputed with the most frequent value and one-hot encoded.

Example 5 introduces a modern LLM stability analysis. It utilizes raw, instance-level prediction traces from the MMLU-Pro benchmark (Wang et al., 2024; arXiv:2406.01574) for 12 state-of-the-art open-weights foundation models (such as Qwen 2.5, Llama 3.3, and Mistral Large), sourced directly from the Hugging Face Open LLM Leaderboard V2 (Fourrier et al., 2024; URL: <https://huggingface.co/spaces/open-llm-leaderboard/open_llm_leaderboard>). These specific open-weights models were selected because they provide the transparent, instance-level result files (JSONL) strictly required to perform HERA's statistical bootstrapping, a requirement not met by proprietary models whose individual prediction traces are kept private.

To create the Example 3 and 4 data and traditional analyses results for comparison against HERA, each dataset was rigorously partitioned into disjoint training and test sets. First, exactly 15,000 instances were sampled for the training set using stratified sampling to preserve the target class distribution. This uniform size was chosen to maintain comparability across the different datasets, keep computational times manageable (particularly for algorithmically expensive models like Support Vector Machines), and to provide a sufficiently large, reusable training corpus that allows adapting algorithms or performing further hyperparameter optimization in future research. Next, exactly 35,000 completely disjoint test instances were sampled and partitioned into 50 independent test blocks, yielding exactly 700 instances per block. These blocks were generated using stratified sampling (Stratified K-Fold) to guarantee representative class distributions in every block. The resulting performance metrics (such as AUROC, Sensitivity, and Specificity) for each block were exported into CSV files.

Furthermore, traditional statistical evaluations were applied to these metric blocks to establish a robust methodological baseline for comparison against analysis with HERA. Specifically, the non-parametric frequentist Friedman test paired with the Nemenyi post-hoc test (Demšar, 2006; URL: <https://jmlr.org/papers/v7/demsar06a.html>) was used for the clinical data. For the weather data, Bayesian Wilcoxon signed-rank tests (Benavoli et al., 2017; URL: <https://jmlr.org/papers/v18/16-305.html>) were applied, utilizing a Region of Practical Equivalence (ROPE) of 0.01 (1%). Both methods represent established standard procedures for comparing multiple classifiers. Finally, these traditional statistical analyses were complemented by a Borda Count aggregation—a well-established consensus ranking method—to generate a final multi-criteria point-estimate ranking for both datasets in comparison to the HERA-generated rankings.

For Example 5, the raw JSON result files were fetched from the official Hugging Face evaluation repositories. To ensure a fair comparison and avoid Missing Not At Random (MNAR) biases caused by unanswered or errored instances, exactly 12,032 MMLU-Pro instances that contained valid prediction traces across all 12 LLMs were intersected. This complete, dense dataset was partitioned into 30 disjoint, independent test blocks using stratified k-fold splitting based on subject categories (e.g., mathematics, physics, law) to ensure representative task distributions in every block. The mean accuracy for each block was then exported into CSV files to evaluate the statistical fragility of traditional point-estimate leaderboards with HERA, directly contrasting HERA's confidence intervals and statistical rankings exclusively against the strict, single-value official MMLU-Pro scores reported by the Open LLM Leaderboard (ignoring other benchmarks like MATH or GPQA).

```text
data/
└── examples/                  % Synthetic example datasets
    ├── Example_1/             % Data for Example 1
    │   ├── CNR.csv            % Contrast-to-Noise Ratio data
    │   ├── OC.csv             % Optical Contrast data
    │   └── SNR.csv            % Signal-to-Noise Ratio data
    ├── Example_2/             % Data for Example 2
    ├── Example_2_workflow/    % Python workflow for stability analysis
    │   ├── data/              % Input data for the workflow
    │   ├── results/           % Results generated by the workflow
    │   ├── run_workflow.py    % Python script to run the workflow
    │   └── requirements.txt   % Python dependencies
    ├── Example_3/             % Clinical data generated from OpenML Cardiovascular Disease dataset
    │   └── Clinical_Data/     % Data blocks for Example 3
    ├── Example_4/             % Weather data generated from WeatherAUS dataset
    │   └── Weather_Data/      % Data blocks for Example 4
    ├── Example_5/             % LLM stability data based on MMLU-Pro benchmark
    │   └── MMLU_Data/         % Data blocks for Example 5
    └── results/               % Results generated by HERA analysis or for comparison
        ├── HERA_Example_1/    % Results for Example 1
        ├── HERA_Example_2/    % Results for Example 2
        ├── Example_3/         % Results for Example 3
        ├── Example_4/         % Results for Example 4
        └── Example_5/         % Results for Example 5
```

## Usage

You can use the datasets in `examples/` to run HERA. To reproduce the results, you can adjust the `configuration.json` files from the HERA_Example folders in the results directory. Please use enough CPU cores to speed up the analysis.

---

## Example Walkthrough of HERA Results

👉 [HERA Example Analysis](https://lerdmann1601.github.io/HERA-Matlab/Example_Analysis)
