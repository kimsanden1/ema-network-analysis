# network-comparisons.R
# author: Kim P.C. van der Sanden
# date: 12/08/2025

# Load required packages
library(dplyr)
library(mlVAR)
library(mnet) # install_github("jmbh/mnet")

#-------------------------------------------------------------------------------
# Network comparisons
#-------------------------------------------------------------------------------

# --- Compare lag-1 vs lag-2 networks --- #
# Compare the fit of the temporal networks
# This test determines if setting the lags argument to 2 significantly improves 
# model fit compared to lags=1
mlVARcompare(mlvar_temp, mlvar_temp2)

# That's it. It extracts the AIC and BIC for each variable from each model
# Then compares these values and summarizes it with counts of which model was
# "best" across all variables


# --- Compare group networks --- #
# Load dataset
data <- readRDS("path/to/dataset.rds")

vars <- c("Fat", "ConR", "TdoR", "PoAf", "NeAf_sqrt", "Pain_sqrt", "PA")

# Check the number of patients and observations per group
data %>%
  filter(!is.na(groups)) %>%
  group_by(groups) %>%
  summarise(
    Npatients = n_distinct(subjno),
    Nobservations = n()
  )

# Use mlVAR_GC to compare network structures between group
# This function is a wrapper for mlVAR that performs a permutation test
# to compare temporal and contemporaneous network structures between groups
# Haslbeck (2025) https://doi.org/10.3758/s13428-024-02541-x

gc <- mlVAR_GC(
  data = data,              # full dataframe with the two groups
  vars = vars,              #
  idvar = "subjno",         # participant indicator
  dayvar = "Daynumber",     # handles day structure of EMA data
  beepvar = "Bleepnumber",  # handles multiple beeps per day structure
  groups = "groups",        # factor indicating group membership
  test = "permutation",     # default, non-parametric
  paired = FALSE,           # groups are independent
  estimator = "lmer",       # multilevel modeling per individual
  temporal = "orthogonal",  # estimates fixed + random effects for lagged edges
  contemporaneous = "orthogonal", # nodewise residual estimation, also uncorrelated
  scale = TRUE,             # standardize before estimation (mlVAR arg)
  nCores = 4,               # how many cores to use in computation (mlVAR arg)
  nP = 1000,                # 100 for piloting, 1000 for stable p-values
  pbar = TRUE               # show progress bar
)
# Duration of the computation in minutes
## > gc$Runtime
## elapsed 
## 50.178

# gc$EmpDiffs shows the raw diffs in network params between group 1 & group 2
# group 1 (PRO-Lung) and group 2 (PRO-RCC)
gc$EmpDiffs

# gc$Pval provides the corresponding p-values from the permutation test
# which tells us whether these observed diffs are statistically significant
gc$Pval

# Show only significant p-values
temporal_sig_pvals <- gc$Pval$Lagged_fixed
temporal_sig_pvals[temporal_sig_pvals >= 0.05] <- NA
print(temporal_sig_pvals)
# The only stat sign diff in the temporal networks between lung cancer and RCC
# [4,][,5] 0.049 (PoAf at time t-1 on NeAf_sqrt at time t)

contemp_sig_pvals <- gc$Pval$Contemp_fixed
contemp_sig_pvals[contemp_sig_pvals >= 0.05] <- NA
contemp_sig_pvals
# The only stat sign diffs in contemp networks between lung cancer and RCC
# [3,][,2] 0.046 (TdoR       <->   ConR)
# [5,][,2] 0.049 (NeAf_sqrt  <->   ConR)

# Observed differences could just be due to chance! Underpowered samples etc.

