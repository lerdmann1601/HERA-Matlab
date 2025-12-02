---
title: 'HERA: A Hierarchical-Compensatory, Effect-Size Driven and Non-parametric Ranking Algorithm using Data-Driven Thresholds and Bootstrap Validation'
tags:
  - MATLAB
  - ranking
  - statistics
  - benchmarking
  - bootstrap
  - effect size
  - multi-criteria decision making
authors:
  - name: Lukas von Erdmannsdorff
    orcid: 0009-0009-3758-7363
    affiliation: 1

affiliations:
 - name: Institute of Neuroradiology, Goethe University Frankfurt
   index: 1
date: 02 December 2025
bibliography: paper.bib
---

## Summary

![HERA Logo](images/hera_logo.png)

In scientific disciplines ranging from clinical research to machine learning, researchers face the challenge of objectively comparing multiple algorithms, experimental conditions, or datasets across a variety of performance metrics. This process, often framed as Multi-Criteria Decision Making (MCDM), is critical for identifying state-of-the-art methods. However, traditional ranking approaches frequently suffer from limitations: they may rely on central tendencies that ignore data variability [@Demsar2006], depend solely on p-values which can be misleading in large samples [@Wasserstein2016], or require subjective weighting of conflicting metrics [@Taherdoost2023].

**HERA** (Hierarchical-Compensatory, Effect-Size Driven Ranking Algorithm) is a MATLAB toolbox designed to automate this comparison process. Unlike weighted-sum approaches, HERA implements a **hierarchical-compensatory logic** that integrates non-parametric significance testing (Wilcoxon signed-rank test), effect size estimation (Cliff's Delta, Relative Difference), and rigorous bootstrapping (Percentile, BCa, Cluster) to produce rankings that are both statistically robust and practically relevant.

## Statement of Need

The scientific community increasingly recognizes the pitfalls of relying on simple summary statistics or p-values alone [@Wasserstein2016]. In benchmarking studies, specifically, several issues persist:

1. **Ignoring Variance**: Ranking based on mean scores fails to account for the stability of performance across different subjects or folds.
2. **Statistical vs. Practical Significance**: A result can be statistically significant but practically irrelevant. Standard tests do not inherently distinguish between these cases [@Lakens2021].
3. **Subjectivity in Aggregation**: Many MCDM methods require users to assign arbitrary weights to metrics (e.g., "Accuracy is 0.7, Speed is 0.3"), introducing bias [@Taherdoost2023].
4. **Distributional Assumptions**: Parametric tests (e.g., t-test) assume normality, which is often violated in real-world benchmarks [@Romano2006].

HERA addresses these challenges by providing a standardized, data-driven framework. It ensures that a method is only ranked higher if it demonstrates a statistically significant and sufficiently large advantage, preventing "wins" based on negligible differences or noise.

## Methodological Framework

HERA operates on paired data matrices where rows represent subjects (or datasets) and columns represent the methods to be compared. The core innovation is its sequential logic, which allows for "compensation" between metrics based on strict statistical evidence.

### Hierarchical-Compensatory Logic

The ranking process is structured as a multi-stage tournament. It does not use a global score but refines the rank order iteratively:

![Hierarchical Logic](images/hierarchical_logic.png)

- **Stage 1 (Initial Sort)**: Methods are initially ranked based on the primary metric \(M_1\). Pairwise comparisons are performed using the Wilcoxon signed-rank test [@Wilcoxon1945], with p-values corrected for multiple comparisons using the step-down Holm-Bonferroni method [@Holm1979] to control the Family-Wise Error Rate (FWER).
- **Stage 2 (Compensatory Correction)**: This stage addresses the trade-off between metrics. A lower-ranked method can "swap" places with a higher-ranked method if it shows a statistically significant and relevant superiority in a secondary metric \(M_2\). This effectively implements a lexicographic ordering with a compensatory component [@Keeney1976].
- **Stage 3 (Tie-Breaking)**: This stage resolves ambiguity using a tertiary metric \(M_3\). It applies two sub-logics:
  - **Sublogic 3a**: A one-time correction if \(M_2\) is statistically neutral.
  - **Sublogic 3b**: An iterative tie-breaking loop if both \(M_1\) and \(M_2\) are statistically neutral.

### Statistical Rigor and Effect Sizes

A comparison between two methods is only counted as a "win" if it satisfies three conjunctive criteria:

- **Significance**: \(p < \alpha_{\text{Holm}}\).
- **Stochastic Dominance**: The effect size Cliff's Delta \(d\) must exceed a threshold \(|d| > \theta_d\) [@Cliff1993]. This ensures robustness against outliers.
- **Magnitude**: The Relative Mean Difference (RelDiff) must exceed a threshold \(r > \theta_r\).

To avoid arbitrary thresholds, HERA employs Percentile Bootstrapping [@Rousselet2021] to determine data-driven cutoffs \(\theta\) based on the distribution of effects within the experiment.

### Validation and Uncertainty

HERA integrates advanced resampling methods to quantify uncertainty:

- **BCa Confidence Intervals**: Bias-Corrected and Accelerated (BCa) intervals are calculated for all effect sizes [@DiCiccio1996].
- **Cluster Bootstrap**: To assess the stability of the final ranking, HERA performs a cluster bootstrap resampling subjects with replacement [@Field2007]. This yields a 95% confidence interval for the rank of each method.
- **Power Analysis**: A post-hoc simulation estimates the probability of detecting relevant effects given the data characteristics.
- **Sensitivity Analysis**: The algorithm permutes the metric hierarchy and aggregates the resulting rankings using a Borda Count [@Young1974] to evaluate the robustness of the decision against hierarchy changes.

## Software Features

HERA offers a flexible configuration to adapt to different study designs:

![Features](images/features.png)

- **Automated Reporting**: Generates PDF reports, Win-Loss Matrices, Sankey Diagrams, and machine-readable JSON/CSV exports.
- **Reproducibility**: Supports fixed-seed execution and configuration file-based workflows.
- **Convergence Analysis**: An adaptive algorithm automatically determines the optimal number of bootstrap iterations \(B\) to ensure stable estimates.

## Future Work

We plan to expand HERA's capabilities in several directions:

- **Python Port**: Developing a Python version to reach a broader data science audience.
- **Bayesian Methods**: Incorporating Bayesian approaches for ranking and probability estimation [@Benavoli2016].
- **Independent Samples**: Extending the framework to handle unpaired data (e.g., randomized clinical trials comparing independent groups).

## Acknowledgements

This software was developed at the Institute of Neuroradiology, Goethe University Frankfurt. I acknowledge the contributions of the open-source community for the development of the underlying statistical methods and MATLAB toolboxes that made this project possible.

## References
