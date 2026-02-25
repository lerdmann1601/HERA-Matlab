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

In scientific disciplines ranging from clinical research to machine learning, researchers face the challenge of objectively comparing multiple algorithms, experimental conditions, or datasets across a variety of performance metrics. This process, often framed as Multi-Criteria Decision Making (MCDM), is critical for identifying state-of-the-art methods. However, traditional ranking approaches frequently suffer from limitations: they may rely on central tendencies that ignore data variability [@Demsar2006; @Benavoli2016], depend solely on p-values which can be misleading in large samples [@Wasserstein2016], or require subjective weighting of conflicting metrics [@Taherdoost2023].

**HERA** (Hierarchical-Compensatory, Effect-Size Driven Ranking Algorithm) is a MATLAB toolbox designed to automate this comparison process, bridging the gap between elementary statistical tests and complex decision-making frameworks. Unlike weighted-sum approaches that collapse multi-dimensional performance into a single scalar, HERA implements a **hierarchical-compensatory logic**. This logic integrates non-parametric significance testing (Wilcoxon signed-rank test), robust effect size estimation (Cliff's Delta, Relative Difference), and bootstrapping (e.g. Percentile and Cluster) to produce rankings that are both statistically robust and practically relevant. HERA is designed for researchers in biomedical imaging, machine learning, and applied statistics who need to compare method performance across multiple quality metrics in a statistically rigorous manner without requiring subjective parameter tuning.

## Statement of Need

The scientific community increasingly recognizes the pitfalls of relying on simple summary statistics or p-values alone [@Wasserstein2016]. In benchmarking studies, specifically, several issues persist:

1. **Ignoring Variance**: Ranking based on mean scores fails to account for the stability of performance across different subjects or folds. A method might achieve a high average score due to exceptional performance on a few easy cases while failing catastrophically on others, yet still outrank a more consistent competitor.
2. **Statistical vs. Practical Significance**: A result can be statistically significant but practically irrelevant, especially in large datasets where even trivial differences yield $p < 0.05$. Standard tests do not inherently distinguish between these cases, potentially leading to the adoption of methods that offer no tangible benefit [@Sullivan2012].
3. **Subjectivity in Aggregation**: Many MCDM methods require users to assign subjective weights to metrics (e.g., "Accuracy is 0.7, Speed is 0.3"). These weights are often chosen post-hoc or lack empirical justification, introducing researcher bias that can be manipulated to favor a specific outcome [@Taherdoost2023].
4. **Distributional Assumptions**: Parametric tests (e.g., t-test) assume normality, which is often violated in real-world benchmarks where performance metrics may be skewed, bounded, or ordinal [@Romano2006].

HERA addresses these challenges by providing a standardized, data-driven framework. It ensures that a method is only ranked higher if it demonstrates a statistically significant and sufficiently large advantage, preventing "wins" based on negligible differences or noise. Unlike existing MCDM software packages such as the Python libraries pyDecision [@Pereira2024] and pymcdm [@Kizielewicz2023], or R's RMCDA [@Najafi2025], which often implement classical methods like TOPSIS [@Hwang1981], PROMETHEE [@Brans1985], and ELECTRE [@Roy1968] that require user-defined weights or preference functions, HERA eliminates subjective parameterization by using data-driven thresholds derived from bootstrap resampling. Furthermore, HERA integrates statistical hypothesis testing directly into the ranking process, a feature absent in standard MCDM toolboxes. While the MATLAB ecosystem offers robust statistical functions, it currently lacks a dedicated, open-source toolbox that unifies this advanced MCDM method with bootstrap validation, forcing researchers to rely on ad-hoc scripts.

## Methodological Framework

HERA operates on paired data matrices where rows represent subjects (or datasets) and columns represent the methods to be compared. The core innovation is its sequential logic, which allows for "compensation" between metrics based on strict statistical evidence.

### Statistical Rigor and Effect Sizes

HERA quantifies differences using statistical significance and effect sizes to ensure practical relevance independent of sample size [@Cohen1988; @Sullivan2012]. A "win" always requires satisfying three conjunctive criteria, if not it is considered "neutral":

- **Significance**: $p < \alpha_{\text{Holm}}$ (Holm-Bonferroni corrected). Pairwise comparisons use the Wilcoxon signed-rank test [@Wilcoxon1945], with p-values corrected using the step-down Holm-Bonferroni method [@Holm1979] to control the Family-Wise Error Rate (FWER).
- **Stochastic Dominance (Cliff's Delta)**:  Cliff's Delta ($d = P(X>Y) - P(Y>X)$) quantifies distribution overlap, is robust to outliers, and relates to common-language effect sizes [@Cliff1993; @Vargha2000]. The effect size $d$ must exceed a bootstrapped threshold $\theta_d$.
- **Magnitude (Relative Difference)**: The Relative Difference (RelDiff) quantifies effect magnitude on the original metric scale, normalized by the mean absolute value. This normalization is formally identical to the Symmetric Mean Absolute Percentage Error (SMAPE) used in forecasting [@Makridakis1993] and conceptually related to the Response Ratio, which uses logarithmic ratios to compare effects across studies [@Hedges1999]. The metric enables scale-independent comparisons and facilitates the interpretation of percentage changes [@Kampenes2007]. RelDiff must exceed a threshold $\delta_{\text{RelDiff}}$.

**Dual Criteria & SEM Lower Bound**
HERA's complementary logic requires both dominance and magnitude, preventing "wins" based on trivial consistent differences or noisy outliers [@Lakens2013]. Thresholds are determined via Percentile Bootstrapping (lower $\alpha/2$-quantile) [@Rousselet2021]. To filter noise in low-variance datasets, the RelDiff threshold enforces a lower bound based on the Standard Error of the Mean (SEM), ensuring $\theta_{r} \geq \theta_{\mathrm{SEM}}$. This approach is inspired by the concept of the "Smallest Worthwhile Change" [@Hopkins2004], but adapted for HERA to quantify the uncertainty of the group mean rather than individual measurement error.

### Hierarchical-Compensatory Logic

The ranking process is structured as a multi-stage tournament. It does not use a global score but refines the rank order iteratively (see Fig. 1):

![Hierarchical-Compensatory Ranking Logic](images/hierarchical_logic.png){height=85%}

- **Stage 1 (Initial Sort)**: Methods are initially ranked based on the win count of the primary metric $M_1$. In case of a tie in adjacent ranks, Cliff's Delta is used to break the tie. If Cliffs Delta is zero, the raw mean values are used to break the tie.
- **Stage 2 (Compensatory Correction)**: This stage addresses the trade-off between metrics. A lower-ranked method can "swap" places with a higher-ranked method if it shows a statistically significant and relevant superiority in a secondary metric $M_2$. This effectively implements a lexicographic ordering with a compensatory component [@Keeney1976], allowing a method that is slightly worse in the primary metric but vastly superior in a secondary metric to improve its standing.
- **Stage 3 (Tie-Breaking)**: This stage resolves "neutral" results using a tertiary metric $M_3$. It applies two sub-logics to ensure a total ordering:
  - **Sublogic 3a**: A one-time correction if the previous metric is "neutral" based on the HERA criteria. This handles cases where two methods are indistinguishable in the second metric while still respecting the initial ranking.
  - **Sublogic 3b**: To resolve groups of remaining undecided methods, an iterative correction loop is applied if both $M_1$ and $M_2$ are "neutral", iteratively using metric $M_3$ until a final stable ranking is found.

### Validation and Uncertainty

HERA integrates advanced resampling methods to quantify uncertainty:

- **BCa Confidence Intervals**: Bias-Corrected and Accelerated (BCa) intervals are calculated for all effect sizes [@DiCiccio1996].
- **Cluster Bootstrap**: To assess the stability of the final ranking, HERA performs a cluster bootstrap resampling subjects with replacement [@Field2007]. This yields a 95% confidence interval for the rank of each method.
- **Power Analysis**: A post-hoc simulation with bootstrap estimates the probability of detecting a "win", "loss" or "neutral" in all tested metrics given the data characteristics.
- **Sensitivity Analysis**: The algorithm permutes the metric hierarchy and aggregates the resulting rankings using a Borda Count [@Young1974] to evaluate the robustness of the decision against hierarchy changes.

## Software Features

HERA offers a flexible configuration of up to three metrics (see Fig. 2). This allows users to adapt the ranking logic to different study designs and needs. It also provides a range of reporting options, data integration, and reproducibility features.

- **Automated Reporting**: Generates PDF reports, Win-Loss Matrices, Sankey Diagrams, and machine-readable JSON/CSV exports.
- **Reproducibility**: Supports fixed-seed execution and configuration file-based workflows. The full analysis state, including random seeds and parameter settings, is saved in a JSON file, allowing other researchers to exactly replicate the ranking results.
- **Convergence Analysis**: To avoid the common pitfall of using an arbitrary number of bootstrap iterations, HERA implements an adaptive algorithm. It automatically monitors the stability of the estimated confidence intervals and effect size thresholds, continuing the resampling process until the estimates converge within a specified tolerance, thus determining the optimal number of iterations $B$ dynamically [@Pattengale2010]. If the characteristics of the data for bootstrapping are known, the number of bootstrap iterations can be set manually.
- **Data Integration**: HERA supports seamless data import from standard formats (CSV, Excel) and MATLAB tables, facilitating integration into existing research pipelines. Example datasets and workflows demonstrating practical applications are included in the repository.
- **Accessibility**: HERA can be easily installed by cloning the GitHub repository and running a setup script, or deployed as a standalone application that requires no MATLAB license. An interactive command-line interface guides users through the analysis without requiring programming expertise, while an API and JSON Configuration allow for automated batch processing.

![Flexible Configuration options for Ranking Logic](images/features.png)

## Acknowledgements

This software was developed at the Institute of Neuroradiology, Goethe University Frankfurt. I thank Prof. Dr. Dipl.-Phys. Ralf Deichmann (Cooperative Brain Imaging Center, Goethe University Frankfurt) for his support during the initial conceptualization of this project. I acknowledge Dr. med. Christophe Arendt (Institute of Neuroradiology, Goethe University Frankfurt) for his supervision and support throughout the project. I also thank Rejane Golbach PhD (Institute of Biostatistics and Mathematical Modeling, Goethe University Frankfurt) for her valuable feedback on the statistical methodology.

## References
