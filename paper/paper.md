---
title: 'HERA: A Hierarchical-Compensatory, Effect-Size Driven and Non-parametric Ranking Algorithm using Data-Driven Thresholds and Bootstrap Validation'
tags:
  - MATLAB
  - ranking
  - statistics
  - non-parametric
  - effect size
  - bootstrapping
  - pairwise comparisons
  - multiple criteria
authors:
  - name: Lukas von Erdmannsdorff
    affiliation: 1
affiliations:
 - name: Institute of Neuroradiology, Goethe University Frankfurt, Germany
   index: 1
date: 01 December 2025
bibliography: paper.bib
---

# Summary

In many scientific disciplines, from clinical research to data science, there is a critical need to objectively compare the performance of different methods, algorithms, or treatments. When decisions rely on complex, paired measurements across multiple quality metrics, establishing a clear "winner" is often challenging. **HERA** (Hierarchical-Compensatory, Effect-Size Driven, and Non-parametric Ranking Algorithm) is a MATLAB toolbox designed to address this problem. It automates the generation of a robust ranking by integrating non-parametric significance testing, effect size quantification, and extensive bootstrap validation into a sequential, hierarchical decision process.

# Statement of Need

Objective ranking methods are crucial in scientific research to evaluate and compare different algorithms or treatments without subjective bias. Traditional ranking approaches often face significant limitations:
1.  **Metric Aggregation:** Simple averaging of metrics (e.g., mean scores) can be sensitive to outliers and ignores the shape of data distributions.
2.  **Statistical Rigor:** Relying solely on p-values can be misleading, especially with large sample sizes where trivial differences become significant. Conversely, small datasets may lack power.
3.  **Subjectivity:** Many Multi-Criteria Decision Making (MCDM) methods require the user to define arbitrary weights for each metric.

HERA addresses these issues by providing a **Hierarchical-Compensatory** system. It allows researchers to define an *a priori* hierarchy of metrics (e.g., Safety > Effectiveness > Efficiency). It then ranks datasets based on this hierarchy but allows for compensation: a lower-priority metric can influence the ranking only if the higher-priority metric is statistically "neutral" (i.e., no significant and practically relevant difference exists).

This approach ensures that the ranking is driven by the data rather than arbitrary weights, making it particularly suitable for benchmarking studies in fields like medical image analysis, engineering, or social sciences.

# Methodological Principles

The HERA algorithm operates on paired datasets (e.g., multiple methods applied to the same subjects). The core ranking logic proceeds in stages:

1.  **Significance and Relevance Check:** For every pair of datasets, HERA performs a Wilcoxon signed-rank test corrected for multiple comparisons using the Holm-Bonferroni method. Crucially, statistical significance alone is not enough to declare a "win". HERA simultaneously requires the effect size to exceed data-driven thresholds:
    * **Stochastic Dominance:** Measured by Cliff's Delta ($d$), ensuring robustness against outliers.
    * **Magnitude:** Measured by the Relative Difference (RelDiff).
    * **Adaptive Thresholds:** Instead of fixed constants, HERA calculates adaptive thresholds using a Percentile Bootstrap on the empirical distribution of effects within the study.

2.  **Sequential Ranking Logic:**
    * **Stage 1 (Initial Sort):** Datasets are ranked by the number of "wins" in the primary metric (Metric 1). Ties are broken by stochastic dominance or mean values.
    * **Stage 2 (Correction):** If a secondary metric (Metric 2) is defined, the algorithm iteratively checks if a lower-ranked dataset significantly "beats" a higher-ranked one in this metric. If so, they swap positions (Correction).
    * **Stage 3 (Tie-Breaking):** A tertiary metric can be used to resolve remaining ties or ambiguities where previous metrics were neutral.

3.  **Validation:**
    * **Bootstrapping:** To quantify ranking uncertainty, HERA employs a **Cluster Bootstrap**. It resamples subjects with replacement and recalculates the entire ranking $B$ times (where $B$ is determined via convergence analysis). This yields 95% confidence intervals for the final rank of each dataset.
    * **BCa Intervals:** Effect sizes are reported with Bias-Corrected and Accelerated (BCa) confidence intervals.
    * **Sensitivity Analysis:** To test the robustness of the chosen metric hierarchy, HERA permutes the order of metrics and calculates a consensus score using the **Borda Count** method.

# Features and Implementation

HERA is implemented as a modular MATLAB package. Key features include:

* **Interactive & Batch Modes:** Users can configure the analysis via a CLI-guided step-by-step process or provide a JSON configuration file for automated pipelines.
* **Automatic Convergence:** The toolbox dynamically determines the optimal number of bootstrap iterations ($B$) required for stable results, stopping when the standard error stabilizes.
* **Post-hoc Power Analysis:** Estimates the probability of detecting significant differences given the observed data distribution using simulation-based methods.
* **Visualization:** Automatically generates publication-ready figures, including:
    * Sankey diagrams visualizing rank shifts across metrics.
    * Win-Loss matrices for pairwise comparisons.
    * Detailed ranking plots with confidence intervals.
    * Bootstrap convergence plots.

The software is designed to be robust against missing data (using pairwise deletion where necessary) and includes a comprehensive unit test suite to ensure statistical accuracy.

# Acknowledgements

We acknowledge the open-source community for providing standard statistical algorithms. The methodology was developed at the Institute of Neuroradiology, Goethe University Frankfurt.

# References
