# Replication Code for: "A Joint Model for Multivariate Longitudinal Biomarkers and Survival Outcomes with a Fast Kalman Filter"

**Authors:** Mingzhao Hu, Ya Luo, Yuedong Wang  
**Journal:** Journal of the American Statistical Association (JASA)  
**Data:** End-stage renal disease (dialysis) cohort, n = 5,707 subjects, 48 monthly visits

---

## Pipeline Structure at a Glance

| Script | Role | Sources | Key outputs |
|---|---|---|---|
| `save_theta_trans.R` | Step 0 — one-time utility | — | theta_trans_3_0504.RData |
| `albumin_fit.R` | Main step 1 — model fit | delta_llh_OU, delta_negllh_OU, filter_estimate_OU_withP, llh_survival_x_OU | theta_est_0613.RData, Figure 2 |
| `bootstrap_ci.R` | Main step 2 — bootstrap CIs | same four (inside parallel loop) | Table 2, Table S4, Figure 4 / S1 |
| `albumin_plot.R` | Main step 3 — figures | filter_estimate_OU_withP | Figure 3 |
| `model_comparison_synth.R` | Secondary — sensitivity analysis | data_preparation_functions, helper_functions_losses, llh_survival_x_OU, filter_estimate_OU_withP, llh_cs_cs, llh_cs_local, llh_joint_all_models | fit_ou/cc/cl_synth.RData, model_comparison_results.RData, model_comparison_plots.pdf |

Library files (sourced, not run directly): `delta_llh_OU.R`, `delta_negllh_OU.R`, `filter_estimate_OU_withP.R`, `llh_survival_x_OU.R`, `llh_cs_cs.R`, `llh_cs_local.R`, `llh_joint_all_models.R`, `data_preparation_functions.R`, `helper_functions_losses.R`.

---

## Overview

This repository contains all code needed to reproduce the tables, figures, and sensitivity analysis results in the paper. The analysis consists of two independent pipelines:

1. **Main analysis pipeline** (`albumin_fit.R` → `bootstrap_ci.R` → `albumin_plot.R`): fits the joint model to the real dialysis data and produces all manuscript figures and tables.
2. **Sensitivity analysis pipeline** (`model_comparison_synth.R`): fits three competing model specifications on synthetic data generated from the fitted model to justify the choice of OU process for individual trajectories.

---

## System Requirements

### R version
R 4.4.1 or later recommended.

### Required R packages
```r
install.packages(c(
  # Rcpp compilation
  "Rcpp", "RcppArmadillo", "inline",
  # Data handling
  "data.table", "dplyr", "reshape2",
  # Optimization and smoothing
  "GA", "assist",
  # Parallel computing (bootstrap)
  "parallel", "foreach", "doParallel",
  # Visualization
  "ggplot2", "gridExtra", "ggpubr",
  # Survival analysis
  "survival", "survminer",
  # Tables
  "xtable",
  # Model comparison
  "pROC"
))
```

### Compilation note
Five files (`llh_survival_x_OU.R`, `filter_estimate_OU_withP.R`, `delta_llh_OU.R`, `delta_negllh_OU.R`, `llh_joint_all_models.R`) use Rcpp/inline to compile C++ code at runtime via `source()`. Compilation takes approximately 30–60 seconds per file on first call. Requires a working C++ compiler:
- **Windows**: Rtools (https://cran.r-project.org/bin/windows/Rtools/)
- **Mac**: Xcode Command Line Tools (`xcode-select --install`) + GFortran
- **Linux**: `gcc` and `gfortran` (standard system packages)

---

## Directory Setup

Place all R scripts and RData files in a single working directory. At the top of each script, set:
```r
ROOT_DIR <- "."   # or the full path to your working directory
setwd(ROOT_DIR)
```
No other path changes are needed.

---

## File Descriptions

### Core likelihood and filter functions

These are sourced automatically by the analysis scripts. They do not need to be run directly.

| File | Description |
|---|---|
| `llh_survival_x_OU.R` | Joint log-likelihood for the cubic spline + OU + survival model. Returns −2 × log-likelihood to be minimized by `optim()`. See Section 3.4 of the paper. |
| `filter_estimate_OU_withP.R` | Kalman filter and smoother. Returns population trajectory `f_filter`, individual deviations `b_filter`, and posterior covariance matrices `P0`–`P3` used for confidence interval computation. See Section 3.3. |
| `delta_llh_OU.R` | Computes the Rosenberg `delta_est` (population spline coefficients) and `neg2llh` for the longitudinal-only OU model. Used during parameter initialization. |
| `delta_negllh_OU.R` | Returns `neg2llh` only. Used with `optim()` for univariate initialization. |
| `llh_cs_cs.R` | Longitudinal log-likelihood for the cubic spline population + cubic spline individual model (sensitivity analysis only). |
| `llh_cs_local.R` | Longitudinal log-likelihood for the cubic spline population + local level individual model (sensitivity analysis only). |
| `llh_joint_all_models.R` | Joint log-likelihoods for the cubic-cubic and cubic-local models (sensitivity analysis only). |
| `data_preparation_functions.R` | Helper functions for constructing data arrays from long-format synthetic data (sensitivity analysis only). |
| `helper_functions_losses.R` | Loss metric functions (MSE, MAE, NLL, Brier, AUC-ROC, C-index) and risk-score computation for the three-model comparison (sensitivity analysis only). Defines `compute_risk_scores_from_filter`, called by `model_comparison_synth.R`. |

### Analysis scripts (run in order for main pipeline)

#### Step 0 (optional utility): `save_theta_trans.R`
One-time helper that regenerates the starting parameter vector. Run only if `theta_trans_3_0504.RData` is missing; otherwise skip directly to Step 1.

**Inputs:**
- `albumin_data_0504.RData`: contains `fit_joint_albumin_1` (preliminary optim object) and `theta_trans_3`

**Outputs:**
- `theta_trans_3_0504.RData`: starting parameter vector for `albumin_fit.R`
- `theta_trans_est_0504.RData`: fitted vector from the preliminary fit (not used downstream)

#### Step 1: `albumin_fit.R`
Fits the joint model to the real dialysis data.

**Inputs:**
- `albumin_data_0509.RData`: preprocessed dialysis cohort data (see Data Objects section below)
- `theta_trans_3_0504.RData`: starting parameter vector for optimization (log-transformed scale)
- `theta_trans_1_0404.RData`: starting values for survival subvector

**Outputs:**
- `albumin_fit.RData`: full workspace including `fit_joint_albumin_1` (optim result object)
- `theta_est_0613.RData`: fitted parameter vector on natural scale
- Figure 2 in manuscript: population trends with 95% bootstrap CI

**Runtime:** approximately 3 hours on a standard desktop (single core).

#### Step 2: `bootstrap_ci.R`
Computes nonparametric bootstrap confidence intervals using 114 resamples (38 batches of 3, run in parallel).

**Inputs:**
- Same data files as `albumin_fit.R`
- `theta_est_0613.RData`: point estimates from Step 1
- `boostrap_realdata_0613.RData`: pre-computed bootstrap results (114 fitted parameter vectors)

**Outputs:**
- `theta_result_k.RData` (per batch): bootstrap parameter vectors
- Table 5 (parameter CIs) and Table 5 Sigma (covariance CIs) printed to console
- Figure 4 (correlation network): requires `boostrap_realdata_0613.RData`

**Runtime:** approximately 18 hours per batch (3 resamples in parallel). Pre-computed results in `boostrap_realdata_0613.RData` allow skipping directly to the table and figure generation block (Part 2 onward).

#### Step 3: `albumin_plot.R`
Generates all manuscript figures.

**Inputs:**
- `albumin_data_0509.RData`
- `theta_fit_theta_result_3.RData`: bootstrap parameter list (point estimate in `[[2]]`)
- `theta_fit_fit_joint_albumin_8.RData`: fitted optim object
- `restandardize_for_plot_0505.RData`: biomarker means and SDs for back-transformation to original clinical units

**Outputs:**
- `0509_6pt_3156.png`: intermediate diagnostic plot on standardized scale (not in manuscript)
- `0509_6pt_130_noCI_0904.pdf`: Figure 3 in manuscript (individual trajectories + six-month-ahead mortality probability, representative patient observed to month 35)

**Runtime:** approximately 10 minutes.

#### Sensitivity analysis: `model_comparison_synth.R`
Fits three model specifications on synthetic data to justify the OU individual process.

**Inputs:**
- `synth_data_train_cubic_OU.RData`: synthetic training data (n = 5,707)
- `synth_data_test_cubic_OU.RData`: synthetic test data (n = 5,707)
- `theta_trans_3_0504.RData`, `theta_trans_3_0503.RData`, `theta_trans_3_0502.RData`, `theta_trans_1_0404.RData`: initial values
- Sources: `data_preparation_functions.R`, `helper_functions_losses.R`, and all five likelihood/filter files

**Outputs:**
- `fit_ou_synth.RData`, `fit_cc_synth.RData`, `fit_cl_synth.RData`: fitted parameters for each model
- `model_comparison_results.RData`: test-set log-likelihood, MSE, AUC, C-index for all three models
- `master_comparison_synth.RData`: full workspace snapshot
- `model_comparison_plots.pdf`: comparison bar charts

**Runtime:** approximately 12 hours total (OU: 3 hrs, cubic-cubic: 6 hrs, cubic-local: 3 hrs).

---

## Data Objects

### Complete RData inventory

All `.RData`/`.Rdata` files in the bundle and which scripts use them. (Extensions are case-sensitive on Mac/Linux; filenames below match the bundle exactly.)

| RData file | Used by | Contains |
|---|---|---|
| `albumin_data_0504.RData` | save_theta_trans.R | `fit_joint_albumin_1` (preliminary optim), `theta_trans_3` |
| `albumin_data_0509.RData` | albumin_plot.R | `data_split_in`, `z_in`, `X00_cube_in` (all subjects) |
| `data_split_in.Rdata` | albumin_fit.R, bootstrap_ci.R | q×npt×m biomarker array |
| `z_in.Rdata` | albumin_fit.R, bootstrap_ci.R | m×npt survival indicator matrix |
| `X00_cube_in.Rdata` | albumin_fit.R, bootstrap_ci.R | m×p×npt covariate array |
| `A_list_in.Rdata` | albumin_fit.R, bootstrap_ci.R | list[npt] of alive-at-t subject indices |
| `AB_list_in.Rdata` | albumin_fit.R, bootstrap_ci.R | list[npt] of present-at-t subject indices |
| `theta_trans_3_0504.RData` | albumin_fit.R, bootstrap_ci.R, model_comparison_synth.R | starting parameter vector (final restart) |
| `theta_trans_3_0503.RData` | albumin_fit.R, bootstrap_ci.R, model_comparison_synth.R | starting parameter vector (2nd restart, blending) |
| `theta_trans_3_0502.RData` | albumin_fit.R, bootstrap_ci.R, model_comparison_synth.R | starting parameter vector (1st restart, blending) |
| `theta_trans_1_0404.RData` | albumin_fit.R, bootstrap_ci.R, model_comparison_synth.R | initial survival subvector |
| `theta_est_0613.RData` | bootstrap_ci.R | fitted parameters, natural scale, length 42 |
| `boostrap_realdata_0613.RData` | bootstrap_ci.R | 114 bootstrap parameter vectors (for Figure 4 / S1) |
| `restandardize_for_plot_0505.RData` | albumin_plot.R | biomarker means and SDs for back-transformation |
| `theta_fit_theta_result_3.RData` | albumin_plot.R | `theta_result` list; point estimate in `[[2]]` |
| `theta_fit_fit_joint_albumin_8.RData` | albumin_plot.R | `fit_joint_albumin` optim object; `$par` = transformed estimate |
| `synth_data_train_cubic_OU.RData` | model_comparison_synth.R | synthetic training cohort |
| `synth_data_test_cubic_OU.RData` | model_comparison_synth.R | held-out synthetic test cohort |
| `model_comparison_results.RData` | (pre-computed result for reviewers) | `comparison_results` list (llh, MSE, AUC, C-index) |

Detailed object descriptions for the primary data file follow.

### `albumin_data_0509.RData`
Contains the preprocessed dialysis cohort (n = 5,707, 48 monthly visits, 5 biomarkers).

| Object | Type | Dimensions | Description |
|---|---|---|---|
| `data_split_in` | array | q × npt × m = 5 × 48 × 5,707 | Biomarker observations. `data_split_in[k, t, i]` = biomarker k for subject i at month t. Biomarkers: ALB (albumin, g/dL), IDWG (interdialytic weight gain, %), SBP (systolic blood pressure, mmHg), NLR (neutrophil-to-lymphocyte ratio, square root), HGB (hemoglobin, g/dL). |
| `z_in` | matrix | m × npt = 5,707 × 48 | Survival indicators. `z_in[i, t]` = 1 if subject i died at month t, 0 if alive, NA if censored before month t. |
| `X00_cube_in` | array | m × p × npt = 5,707 × 6 × 48 | Covariate array. `X00_cube_in[i, , t]` = (age, sex, EKTV, diabetes, cardiovascular disease, vintage) for subject i at month t. |

Derived index lists (reconstructed in each script from `z_in`):

| Object | Type | Length | Description |
|---|---|---|---|
| `A_list_in` | list | npt = 48 | `A_list_in[[t]]` = integer vector of subject indices alive (not yet dead or censored) at month t. Used in Kalman filter to compute population trajectory. Notation: set $\mathcal{A}_t$ in paper. |
| `AB_list_in` | list | npt = 48 | `AB_list_in[[t]]` = integer vector of subject indices present at month t (alive or dying at t). Used in survival likelihood. Notation: set $\mathcal{B}_t$ in paper. |

### `theta_trans_3_0504.RData`
Starting parameter vector `theta_trans_3` for joint model optimization. All parameters on log-transformed scale (see Parameter Structure below).

### `restandardize_for_plot_0505.RData`
Biomarker-specific means (`mean_alb`, `mean_idwg`, etc.) and standard deviations (`sd_alb`, `sd_idwg`, etc.) computed from the training data. Used to back-transform standardized model outputs to original clinical units for plotting.

---

## Parameter Vector Structure

The full parameter vector `theta_trans` (length 42 for q = 5, p = 6) is optimized on the log-transformed scale. Back-transform the first 4q elements via `exp()` to recover natural-scale parameters.

| Index | Symbol (paper) | Description |
|---|---|---|
| 1–5 | log(ζ₁,...,ζ₅) | Population cubic spline smoothing variances, one per biomarker |
| 6–10 | log(ν₁²,...,ν₅²) | OU individual innovation variances |
| 11–15 | log(ξ₁,...,ξ₅) | OU mean-reversion rates |
| 16–20 | log(d₁,...,d₅) | Diagonal of Cholesky factor of cross-biomarker covariance Σ |
| 21–30 | l₁,...,l₁₀ | Lower-triangular elements of Cholesky factor of Σ (not log-transformed) |
| 31 | log(κ) | Survival baseline hazard scale |
| 32–36 | γ₂₁,...,γ₂₅ | Association parameters linking individual OU state to survival log-odds |
| 37–42 | γ₀,γ₁₁,...,γ₁₆ | Survival covariate coefficients (intercept + 6 covariates) |

General formula: total length = 0.5q² + 4.5q + p + 1

The Cholesky factorization gives Σ = LDLᵀ where D = diag(d₁,...,dq) and L is lower triangular with ones on the diagonal and l_vec in the lower triangle. Back-transform: `theta_est[1:(4*q)] <- exp(theta_trans[1:(4*q)])`.

---

## Reproducibility Notes

1. **Random seed:** `set.seed(240)` is set at the top of `albumin_fit.R` and `set.seed(1)` in `bootstrap_ci.R`. Bootstrap results may differ slightly across platforms due to platform-specific random number generation.

2. **Runtime:** Total runtime for the full pipeline (fit + bootstrap + plots) is approximately 3 hours for fitting and 18 hours per bootstrap batch. Pre-computed results (`theta_est_0613.RData`, `boostrap_realdata_0613.RData`) allow skipping the long-running steps.

3. **Rcpp compilation:** C++ code is compiled at runtime via `inline::cxxfunction`. Compilation output is printed to the console and can be ignored. If compilation fails, verify that Rtools (Windows) or Xcode/GFortran (Mac) is correctly installed.

4. **Data availability:** The dialysis cohort data (`albumin_data_0509.RData`) is available via the UCSB Box folder shared with reviewers. Due to patient privacy, the data cannot be deposited on a public repository. The synthetic data files (`synth_data_train/test_cubic_OU.RData`) are publicly available and sufficient to reproduce the sensitivity analysis results in Section 5.

