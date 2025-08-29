![R version](https://img.shields.io/badge/R-version_4.4.3%20(2025--02--28)%20%22Trophy%20Case%22-276DC3?logo=r)

![tidyr version](https://img.shields.io/badge/tidyr-version_1.3.1-276DC3?logo=r)
![dplyr version](https://img.shields.io/badge/dplyr-version_1.1.4-276DC3?logo=r)
![mlVAR version](https://img.shields.io/badge/mlVAR-version_0.5.2-276DC3?logo=r)
![qgraph version](https://img.shields.io/badge/qgraph-version_1.9.8-276DC3?logo=r)
![ggplot2 version](https://img.shields.io/badge/ggplot2-version_3.5.2-276DC3?logo=r)
![mnet version](https://img.shields.io/badge/mnet-version_0.1.4-276DC3?logo=r)
![mice version](https://img.shields.io/badge/mice-version_3.17.0-276DC3?logo=r)
![imputeTS version](https://img.shields.io/badge/imputeTS-version_3.3-276DC3?logo=r)

---

# EMA network analysis
Code used for article
*"Ecological momentary assessment in patients with lung and renal cell cancer: feasibility and dynamic networks of fatigue, affect, and pain"*

This repository contains the code related to the **1) sensitivity analysis, 2) feasibility analysis, 3) main analysis,** and **4) network comparisons**. 

All code is written by the first author.

---
### Requirements
```r
# pkgs
library(tidyr)     # for tidying data
library(dplyr)     # for data manipulation
library(mlVAR)     # for multilevel vector autoregression
library(qgraph)    # for visualizing network structures
library(ggplot2)   # for plotting
library(mnet)      # install_github("jmbh/mnet"), for network comparison
library(mice)      # for multilevel mice
library(imputeTS)  # for time-series imputation, specifically Kalman filtering
```


## 1. `sensitivity-analysis.R`
**Summary**: This script performs a sensitivity analysis to determine how different data handling strategies impact the resulting network structures. It compares three methods for dealing with missing data (no imputation, Kalman filtering, and multilevel MICE) across four different participant compliance thresholds (0%, 25%, 50%, and 75%). For each of the resulting twelve datasets, the script fits a multilevel vector autoregression (mlVAR) model to estimate temporal and contemporaneous networks. Finally, it generates and saves plots of these networks.

### Overview
The script's main purpose is to perform a sensitivity analysis. It begins by defining a helper function `filter_compliance` to select participants based on their compliance rate. The core of the script involves generating 12 datasets by combining three imputation methods (no imputation, Kalman filtering, and multilevel MICE) with four participant compliance thresholds (0%, 25%, 50%, and 75%). For each of these datasets, an mlVAR model is fitted using a set of seven variables: Physical fatigue (`Fat`), Concentration problems (`ConR`), Problems with motivation to do fun things (`TdoR`), Positive affect (`PoAf`), Negative affect (`NeAf`), Pain (`Pain` or sometimes referred to as `PA`), and Physical activity (`PA` or sometimes referred to as `PH`). The script then extracts the contemporaneous and temporal network matrices from each model and visualizes them in two separate 3x4 grid plots using the `qgraph` package. 

---

## 2. `feasibility-analysis.R`
**Summary**: This script performs a feasibility analysis of an ecological momentary assessment (EMA) study. It processes raw EMA data to assess participant compliance over a 21-day period. The script calculates key metrics such as the daily percentage of completed prompts and each participant's average compliance. It then categorizes participants into `high` and `low` compliance groups based on a 75% threshold and reports descriptive statistics for each group. Finally, the script generates a heatmap visualization that displays each participant's daily compliance percentage, providing a clear visual overview of engagement throughout the study.

### Overview
The script begins by loading the full dataset and preprocessing it and to anonymize participant identifiers. It then focuses on compliance analysis using the Fat variable as an indicator of a completed prompt. The script calculates the number of completed prompts per day for each participant and the daily percentage based on a target of 5 prompts per day. These daily percentages are then used to compute each participant's average compliance over the study period.

The analysis further breaks down compliance into three key time intervals: Day 1, Days 2-20, and Day 21, reporting the mean and standard deviation for each period to identify potential changes in engagement over time. Participants are also assigned to either a `high` or `low` compliance group based on whether their average compliance is above or below a 75% threshold.

A major output of the script is a heatmap visualization created with `ggplot2`. This heatmap graphically represents the daily compliance percentage for each participant, with a color gradient ranging from red (low compliance) to green (high compliance). The plot is customized with anonymized participant labels and a simplified x-axis, making it easy to see compliance patterns for individual participants over the study's 21 days.

---

## 3. `main-analysis.R`
**Summary**: This script conducts the main network analysis for the study, using an imputed dataset derived from participants with at least 75% compliance. The analysis focuses on three primary goals: 1) estimating and visualizing networks, 2) calculating and plotting centrality metrics, and 3) assessing the stability of those metrics. It uses mlVAR models to create contemporaneous and temporal networks (at `lags=1` and `lags=2`). From these networks, it calculates the centrality indices (strength, in-strength, and out-strength), and then visualizes them with `ggplot2`. To ensure the robustness of the findings, it performs a case-drop bootstrap stability analysis, which involves repeatedly dropping 20% of participants, re-estimating the models, and computing the correlation of centrality metrics with the full model's results. The stability is then visualized using histograms.

### Overview
The script starts by loading a specific imputed dataset (Kalman filtering on participants with >75% compliance) and defines a variable set for the contemporaneous network and a seperate variable set for the temporal network, which excludes physical activity due to its retrospective phrasing being deemed incompatible for lagged analyses.

The first part of the analysis involves **mlVAR model estimation**. Three models are created:
- A model (`mlvar_con`) to estimate the contemporaneous network.
- A model (`mlvar_temp`) with the default lag of 1 to estimate the lag-1 temporal network.
- A model (`mlvar_temp2`) with a lag of 2 to estimate the lag-2 temporal network.

The script then visualizes these three networks using the `qgraph` package, defining custom node labels and colors for clarity, and saves the plots.

The second part focuses on **centrality analysis**. The script calculates the `strength` centrality for the contemporaneous network and `in-strength` and `out-strength` for the temporal networks. It then uses `ggplot2` to create two plots: one showing the strength of each node in the contemporaneous network and another comparing the in- and out-strength of the lag-1 and lag-2 temporal networks.

The final and most complex part is the **stability analysis**. A case-drop bootstrapping procedure (as suggested by [10.3758/s13428-017-0862-1](10.3758/s13428-017-0862-1) and as seen in [10.1016/j.schres.2019.10.055](10.1016/j.schres.2019.10.055)) is implemented to test the robustness of the centrality measures. For mlVAR models, the `bootnet` package doesn't seem to support bootstrapping directly within `mlVAR()` objects. However, the script contains a workaround pipeline, adapted to case-dropping subset bootstraps for mlVAR networks. For both the temporal and contemporaneous networks, the script runs 100 iterations. In each iteration, it randomly resamples 80% of the participants, fits a new mlVAR model, and calculates the centrality indices. The Spearman correlation between the bootstrapped centrality values and the full model's centrality values is computed. These correlations are then visualized in a series of histograms to show the stability distribution, providing a quantifiable measure of confidence in the network structures.

---

## 4. `network-comparisons.R`
**Summary**: This script performs two network comparisons: 
1. Evaluate whether a temporal network with a lag of 2 provides a significantly better fit than a network with a lag of 1.
2. Compare the network structures between two distinct patient groups (lung cancer vs. renal cell cancer). 
The lag comparison is done using the `mlVARcompare` function, which contrasts model fit statistics like AIC and BIC.
The group comparison is conducted using the `mlVAR_GC` function, which performs a permutation test to identify statistically significant differences in both temporal and contemporaneous network connections between the two patient cohorts. Note that these findings may be underpowered.

### Overview
The script is divided into two main sections for network comparison.

The first section compares the fit of a lag-1 temporal network against a lag-2 temporal network using the `mlVARcompare` function. This function compares the AIC and BIC for each model, variable by variable, to determine which one provides a better balance of complexity and goodness-of-fit. The output is a summary of which model is "best" based on these criteria.

The second section focuses on comparing network structures between two groups: patients with lung cancer and patients with renal cell cancer (RCC). The script first loads the dataset and summarizes the number of participants and observations in each group. The core of this analysis uses the `mlVAR_GC` function, which is specifically designed for group-level comparisons of mlVAR models. This function employs a permutation test to evaluate if the observed differences in network parameters (both temporal and contemporaneous) between the two groups are statistically significant. The script is configured to run 1,000 permutations to ensure stable p-values which took about 50 minutes.

The script's output from the `mlVAR_GC` function includes the raw differences in network parameters between the groups (`gc$EmpDiffs`) and the corresponding p-values (`gc$Pval`). The script specifically highlights the significant findings by filtering for p-values less than 0.05. It identifies one significant difference in the temporal network (the effect of `PoAf` on `NeAf_sqrt`) and two significant differences in the contemporaneous network (the connections between `TdoR` and `ConR`, and `NeAf_sqrt` and `ConR`). The script concludes with a cautionary note that the observed differences might be a result of chance due to potentially underpowered samples.

---

### Citation
If you use this code or analysis in your own work, please cite this repository or reference it as part of the appendix of any version of the paper, if available, the published version.

> **Van der Sanden, K. (2025), Schellekens, M., Mols, F., van Deun, K., van der Lee, M., ... de Rooij, B. (2025).** _Ecological momentary assessment in patients with lung and renal cell cancer: feasibility and dynamic networks of fatigue, affect, and pain._ Manuscript in preparation. 


