# Repository Structure

The codebase is organized as a MATLAB package (`+HERA`) to ensure namespace isolation.

```text
HERA-Matlab/
├── +HERA/                     % Main Package Namespace
│   ├── +output/               % Report Generation (PDF, JSON, CSV)
│   ├── +plot/                 % Visualization (Sankey, Heatmaps)
│   ├── +run/                  % Execution Logic
│   ├── +start/                % CLI & Configuration Logic
│   ├── +stats/                % Statistical Core (Cliff's Delta, Convergence Check)
│   ├── +test/                 % Unit Test Suite
│   ├── language/              % Localization Files
│   ├── bootstrap_ranking.m    % Cluster Bootstrap Analysis
│   ├── borda_ranking.m        % Consensus Ranking
│   ├── calculate_bca_ci.m     % BCa Confidence Intervals
│   ├── calculate_ranking.m    % Core Ranking Logic
│   ├── calculate_thresholds.m % Threshold Calculation
│   ├── default.m              % Global Defaults
│   ├── design.m               % Style & Design Definitions
│   ├── generate_output.m      % Output Generation Controller
│   ├── generate_plots.m       % Plot Generation Controller
│   ├── get_language.m         % Language Loader
│   ├── get_version.m          % Version Retrieval
│   ├── load_data.m            % Data Import & Validation
│   ├── power_analysis.m       % Power Analysis
│   ├── run_ranking.m          % Core Function (Developer API)
│   ├── start_ranking.m        % Main Entry Point (User API)
│   └── run_unit_test.m        % Test Runner
├── assets/                    % Images & Logos
├── data/                      % Data Directory
│   ├── examples/              % Synthetic Example Datasets
│   └── README.md              % Data Documentation
├── deploy/                    % Build Scripts (Standalone App)
├── docs/                      % Documentation & Guides
├── paper/                     % Paper Resources
├── tests/                     % Unit Test & Analysis Reports
├── setup_HERA.m               % Path Setup Script
├── CITATION.cff               % Citation Metadata
├── CODE_OF_CONDUCT.md         % Community Standards
├── CONTRIBUTING.md            % Contribution Guidelines
├── LICENSE                    % License File
└── README.md                  % Global Documentation
```
