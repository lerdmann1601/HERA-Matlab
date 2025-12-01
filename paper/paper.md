---
title: 'HERA: A Hierarchical-Compensatory, Effect-Size Driven and Non-parametric Ranking Algorithm using Data-Driven Thresholds and Bootstrap Validation'
tags:
  - MATLAB
  - ranking
  - statistics
  - non-parametric
  - effect size
  - bootstrapping
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

In many scientific disciplines, there is a need to objectively compare the performance of different methods. HERA addresses this through a non-parametric sequential ranking based on significance tests (Wilcoxon) and effect sizes (Cliff's Delta). It utilizes bootstrapping for data-driven thresholds and stability analyses.

# Statement of Need

Objective ranking methods are crucial in scientific research to evaluate and compare different algorithms or treatments without subjective bias. Traditional ranking often relies on mean values or simple significance tests, which may not capture the practical significance of differences or robustness against outliers.

HERA fills this gap by providing a hierarchical and compensatory ranking system. It integrates:
- **Non-parametric tests**: Using Wilcoxon signed-rank tests to handle non-normal data distributions.
- **Effect size integration**: Incorporating Cliff's Delta to ensure that statistically significant differences are also practically meaningful.
- **Data-driven thresholds**: employing bootstrapping to determine adaptive thresholds for what constitutes a significant difference, rather than relying on arbitrary fixed values.
- **Robustness**: Features like Bias-Corrected and Accelerated (BCa) confidence intervals and power analysis ensure the reliability of the rankings.

This approach allows researchers to obtain a nuanced and robust ranking of methods, facilitating better decision-making and fairer comparisons in benchmarking studies.
