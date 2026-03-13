# Changelog

All notable changes to this project will be documented in this file.
> **Note:** Only the latest version is actively supported. Previous versions (including all versions prior to 1.4.0) are considered legacy. Corresponding GitHub releases may have been removed, but the code remains accessible via Git tags for reproducibility.

## [Unreleased]

### Added

- **8th Convergence Scenario**: Introduced a new `Sensitivity` scenario to the convergence analysis suite.
- **JSON Scenario Configuration**: All convergence analysis scenarios are now fully customizable through the JSON configuration file.
- **Automatic Update Checks**: HERA now automatically checks for newer versions on startup for both MATLAB (Toolbox) and Python (PyPI).
- **Min/Max Columns in CSV**: Added `Min` and `Max` columns for Error and Cost metrics in convergence CSV exports.

### Fixed

- **Parallel Scheduling**: Optimized worker saturation by separating strategies for simulation preparation and testing phases.
- **CSV Delimiter**: Switched remaining CSV output to semicolon delimiters for consistent spreadsheet compatibility.
- **Figure Rendering**: Forced layout computation before saving to prevent rendering issues in automated Convergence Analysis reports.

---

## [1.3.1] - 2026-03-09

### Added

- **Citation Request**: Added an automated citation request to the console output.
- **Dynamic Confidence Intervals**: CI levels are now dynamically displayed across plots, logs, and output tables.

### Changed

- **Plot Optimization**: Adjusted x-axis limits with dynamic margins and removed KDE from BCa distribution plots for improved accuracy.
- **Simplified Branding**: Simplified the software's full title across all documentation.

---

## [1.3.0] - 2026-03-03

### Added

- **Convergence Analysis v2.0**: Massive update featuring a JSON-driven configuration engine for automated studies.
- **Resource-Aware Scaling (DRAS)**: Intelligent workload management that calculates simulation density based on available system memory.
- **Bit-Perfect Reproducibility**: Multi-tier seeding architecture ensuring strictly isolated random streams across hardware.
- **Absolute Error Metrics**: Transitioned threshold convergence calculations to absolute deviation for better representation of small errors.

---

## [1.2.1] - 2026-02-01

### Added

- **Deployment Infrastructure**: Included the complete `deploy/` directory (build scripts for MATLAB, Python, and Standalone) in the distributed toolbox.

### Fixed

- **Ranking Robustness**: Added epsilon-based comparisons in `calculate_ranking.m` to prevent instabilities when effect sizes are near zero.

---

## [1.2.0] - 2026-01-30

### Added

- **Scientific Validation Mode**: Introduced a dedicated **Convergence Analysis** mode to validate the stability of robust convergence parameters.
- **Smart Type Conversion**: The Python wrapper now automatically detects and converts **NumPy arrays** and **Pandas DataFrames/Series**.
- **Python Type Stubs**: Added `.pyi` files for better autocompletion and static analysis.

### Changed

- **Documentation Overhaul**: Updated README and guides to match the first official toolbox and Python release structure.

---

## Archived Versions (Deprecated)

### [1.1.2] - 2026-01-26

- Maintenance update for example workflows and installation logic refinements.

### [1.1.1] - 2025-01-18

- MkDocs integration and DOI (Zenodo) metadata sync.

### [1.1.0] - 2026-01-16

- Initial release of the Python Interface and MATLAB Toolbox (`.mltbx`).

### [1.0.0] - 2025-12-17

- **Initial Release of the HERA MATLAB Runtime for MacOS.**

[Unreleased]: https://github.com/lerdmann1601/HERA-Matlab/compare/v1.3.1...HEAD
[1.3.1]: https://github.com/lerdmann1601/HERA-Matlab/releases/tag/v1.3.1
[1.3.0]: https://github.com/lerdmann1601/HERA-Matlab/releases/tag/v1.3.0
[1.2.1]: https://github.com/lerdmann1601/HERA-Matlab/releases/tag/v1.2.1
[1.2.0]: https://github.com/lerdmann1601/HERA-Matlab/releases/tag/v1.2.0
[1.1.2]: https://github.com/lerdmann1601/HERA-Matlab/releases/tag/v1.1.2
[1.1.1]: https://github.com/lerdmann1601/HERA-Matlab/releases/tag/v1.1.1
[1.1.0]: https://github.com/lerdmann1601/HERA-Matlab/releases/tag/v1.1.0
[1.0.0]: https://github.com/lerdmann1601/HERA-Matlab/releases/tag/v1.0.0
