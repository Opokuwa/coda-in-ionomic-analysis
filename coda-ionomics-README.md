# Compositional Methods in Ionomic Data Analysis

My master's research report for the MSc Statistics program at Oklahoma State University (April 2026), supervised by Dr. Pratyaydipta Rudra.

## Overview

Ionomic data — measurements of elemental concentrations in biological samples — are inherently **compositional**: each observation is a set of relative abundances that sum to a constant. Despite this, ionomic datasets are routinely analyzed using standard statistical methods applied directly to raw proportions, which can produce misleading inference.

This study evaluates the practical consequences of ignoring compositional structure, using a dataset of **113 spring amphipod samples from Oklahoma**, each measured on **28 elemental components**.

## Research Questions

1. Do raw proportions, log-transformed proportions, and centered log-ratio (CLR) representations produce meaningfully different results in practice?
2. Which method offers the best statistical performance (Type I error control, power)?
3. When does the choice of representation matter most — and what property distinguishes CLR as the principled choice?

## Key Findings

- **Raw proportions are inadequate** for formal ionomic inference: they yield substantially lower statistical power and unstable model diagnostics.
- **Log-transformed proportions and CLR** perform similarly on main simulation metrics (Type I error, power), but are not interchangeable.
- **CLR is preferable** over raw proportions when the objective is inference that respects the compositional structure of the data, offering more stable diagnostics and greater power.
- **CLR is not subcompositionally coherent** — CLR values change when the set of components being analyzed changes, because the geometric mean shifts. **Pairwise log-ratios (PLR)** provide the strongest framework when stable inference across subsets of elements is required, as they are invariant to which other components are included.
- Analytical method choice — independent of data quality — can meaningfully alter biological conclusions.

## Methodology

### Real Data Analysis
- Element-by-element linear models (Sex, Aquifer, Region as predictors) fitted under Raw and CLR representations
- Side-by-side comparison of significance patterns, model diagnostics, and PCA structure
- Heatmaps of coefficient-level and ANOVA-level significance across all 28 elements

### Simulation Study
- Simulated Dirichlet-based compositional data with known ground-truth effects
- Evaluated Type I error and power for Raw, LogRaw, and CLR across 1,000 iterations
- Demonstrated subcompositional coherence failure via pairwise log-ratio (PLR) comparison

## Repository Structure

```
├── linreg.R               # Element-by-element linear models (Raw + CLR) — RUN FIRST
├── modelassumption.R      # Diagnostic plots and Shapiro-Wilk / DW tests (requires linreg.R)
├── allbarplots.R          # Stacked and grouped bar plots of elemental compositions
├── PCAs.R                 # PCA scores, loadings, biplots, highlighted aquifer plots
├── heatmaps.R             # Coefficient-level significance heatmaps (requires linreg.R)
├── anovaheatmaps.R        # ANOVA-level significance heatmaps (requires linreg.R)
├── Reproduce_mainSim.R    # Main simulation study — Type I error and power (Table 4.1, Figs 4.1–4.2)
├── Reproduce_subcomp.R    # Subcompositional coherence demonstration (Table 4.2, Fig 4.3)
└── README.md
```

### Run Order

```
1. linreg.R               ← start here for real-data analysis
2. modelassumption.R      ← requires linreg.R objects
3. heatmaps.R             ← requires linreg.R objects
4. anovaheatmaps.R        ← requires linreg.R objects
5. allbarplots.R          ← independent (loads data directly)
6. PCAs.R                 ← independent (loads data directly)
7. Reproduce_mainSim.R    ← independent simulation (loads data for n, D, mu_hat only)
8. Reproduce_subcomp.R    ← fully standalone, no data needed
```

> **Note:** The amphipod dataset is not included in this repository. Update the file path in each script that loads data to point to your local copy.

## Requirements

```r
library(compositions)   # CLR transformation, acomp()
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
library(broom)
library(lmtest)         # Durbin-Watson test
library(purrr)
library(grid)
library(scales)
```

## Report

The full master's report — *"The Practical Importance of Compositional Methods in Ionomic Data Analysis"* — is available in this repository.

**Committee:**  
Supervisor: Dr. Pratyaydipta Rudra  
Members: Dr. Joshua Habiger, Dr. Ye Liang

## Author

**Miriam Tenkorang**  
MSc Statistics, Oklahoma State University (2024–2026)  
Actuarial Candidate — SOA Exam P (Passed May 2025)  
[LinkedIn](https://www.linkedin.com/in/miriam-tenkorang-a15695263) | [GitHub](https://github.com/Opokuwa)
