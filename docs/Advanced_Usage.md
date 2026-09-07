# MATLAB Integration & Advanced Usage

This guide covers everything for MATLAB users and developers: from installing the official toolbox (`.mltbx`) and running analyses to programmatic in-memory execution and building the toolbox from source.

---

## 1. Installation Guide (MATLAB Users)

### System Requirements

* **MATLAB**: Version R2020a or later (R2025b recommended).
* **Statistics and Machine Learning Toolbox** (Required for non-parametric tests and distributions).
* **Parallel Computing Toolbox** (Required for parallelized bootstrap execution).

### Option A: MATLAB Toolbox (.mltbx) — Recommended for Users

The easiest way to install HERA in MATLAB is via the packaged toolbox:

1. Download the latest `HERA_<version>.mltbx` from the [GitHub Releases](https://github.com/lerdmann1601/HERA-Matlab/releases) page or the [MATLAB File Exchange](https://de.mathworks.com/matlabcentral/fileexchange/183089-hera).
2. In MATLAB, navigate to your download folder and **double-click** the `.mltbx` file (or right-click and select *Install*).
3. MATLAB will automatically install the toolbox to your Add-Ons directory and configure your search path.
4. HERA is now globally available in MATLAB via the `HERA` package namespace.

### Option B: Git Clone & Path Setup (For Developers / Source Code)

If you wish to work with the active source code or contribute to development:

1. Clone the repository:

   ```bash
   git clone https://github.com/lerdmann1601/HERA-Matlab.git
   ```

2. Open MATLAB, navigate to the cloned repository root, and run the setup script:

   ```matlab
   cd HERA-Matlab
   setup_HERA
   ```

   This validates your MATLAB version and permanently adds the required folders to your MATLAB search path.

---

## 2. Standard Execution Modes

### 1. Interactive Mode (Recommended for Beginners)

Start the guided command-line interface:

```matlab
HERA.start_ranking()
```

The interactive prompt guides you through data loading, metric selection, ranking modes, and statistical parameters.

### 2. Config-Driven Mode (Single Run or Automated Pipeline)

For reproducible research or scripted workflows, pass a JSON configuration file:

```matlab
HERA.start_ranking('configFile', 'path/to/analysis_config.json')
```

For parameter specifications and templates, see [Configuration & Parameters](Configuration_&_Parameters.md).

### 3. Verification & Unit Testing

Execute the built-in 46-test validation suite to verify algorithmic integrity:

```matlab
HERA.start_ranking('runtest', 'true')
```

### 4. Convergence Analysis

Evaluate bootstrap stability and determine optimal iterations:

```matlab
HERA.start_ranking('convergence', 'true')
```

For more details, see [Convergence Analysis Documentation](Convergence_Analysis.md).

---

## 3. Programmatic Developer API (`HERA.run_ranking`)

Advanced users can call `HERA.run_ranking` directly with in-memory data matrices, completely bypassing file I/O. This is ideal for integration into simulation pipelines or larger data processing workflows.

### Syntax

```matlab
results = HERA.run_ranking(userInput);
```

### Input Structure (`userInput`)

Instead of providing a `folderPath` to CSV or Excel files, pass data matrices directly via `custom_data`:

```matlab
% 1. Prepare Data
% Cell array of matrices: each matrix is [n_Subjects x n_Methods]
data_m1 = randn(50, 5); % Metric 1 (e.g., Accuracy)
data_m2 = randn(50, 5); % Metric 2 (e.g., Runtime)
custom_data = {data_m1, data_m2};

% 2. Configure Parameters
userInput = struct();
userInput.custom_data   = custom_data;
userInput.metric_names  = {'Accuracy', 'Runtime'};
userInput.dataset_names = {'Method A', 'Method B', 'Method C', 'Method D', 'Method E'};
userInput.ranking_mode  = 'M1_M2';
userInput.output_dir    = pwd;

% 3. Execute Analysis
results = HERA.run_ranking(userInput);
```

### Output Structure (`results`)

The returned struct contains all computed statistics and ranking results:

* `.final_rank`: Final ranking vector across all methods.
* `.d_vals_all`: Effect sizes (Cliff's Delta).
* `.p_vals_all`: Raw p-values from Wilcoxon signed-rank tests.
* `.ci_d_all`: Bootstrap confidence intervals.

For the complete reference of all returned fields, see [Results Structure Reference](Results_Structure_Reference.md).

---

## 4. Packaging the MATLAB Toolbox from Source (For Developers)

> [!IMPORTANT]
> This section is strictly for **developers and maintainers** packaging HERA as a `.mltbx` file for distribution on GitHub Releases or the MATLAB File Exchange.

### Developer Requirements

* **MATLAB** (R2020a or later, R2025b recommended).
* A clean working directory with the appropriate release tag checked out (e.g., `v1.4.6`).

### Packaging Procedure

1. Open MATLAB and navigate to the project root directory.
2. Run the toolbox packaging script:

   ```matlab
   cd deploy
   package_HERA_toolbox
   ```

**What this script performs:**

1. **Pre-Build Cleanup**: Removes OS artifacts (`.DS_Store`) and Python cache directories (`__pycache__`, `*.pyc`) that can cause validation errors during File Exchange uploads.
2. **Metadata Configuration**: Uses `matlab.addons.toolbox.ToolboxOptions` to set the display name, synchronized version number, author, summary, and description.
3. **Selective File Bundling**: Explicitly packages `+HERA`, example datasets (`data/examples/`), documentation, licensing, and setup scripts, while excluding development scripts (`deploy/`), tests, and virtual environments to ensure a lean production file.
4. **Artifact Generation**: Produces the `.mltbx` file in `deploy/output/toolbox/` (e.g., `HERA_v1.4.6.mltbx`).

### Distribution

* **GitHub Releases**: Upload the generated `.mltbx` file from `deploy/output/toolbox/` as a release asset.
