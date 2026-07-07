# Replication Code for: "A Joint Model for Multivariate Longitudinal Biomarkers and Survival Outcomes with a Fast Kalman Filter"

**Authors:** Mingzhao Hu, Ya Luo, Yuedong Wang  
**Data:** End-stage renal disease dialysis cohort, n = 5,707 subjects, 48 monthly visits

---

## Overview

This repository contains replication code for the real-data analysis, bootstrap confidence intervals, simulation summaries, and computational timing comparison reported in the manuscript.

The code is organized into four related pipelines:

1. **Main real-data analysis pipeline**  
   `save_theta_trans.R` → `albumin_fit.R` → `bootstrap_ci.R` → `albumin_plot.R`

2. **Sensitivity-analysis pipeline**  
   `model_comparison_synth.R`

3. **Simulation/bootstrap coverage replication pipeline**  
   `bootstrap_analysis.R`, with optional supporting scripts  
   `bootstrap_generate.R`, `bootstrap_generate_exact100_from_archive.R`, and `bootstrap_collect.R`

4. **Filtering-time comparison pipeline**  
   `filter_time_comparison.R`

For exact replication of the simulation/bootstrap coverage table, use the archived final bootstrap object:

```text
theta_mat_est_extra_final.Rdata
```

rather than rerunning the stochastic bootstrap-generation script.

---

## Pipeline Structure at a Glance

| Script | Role | Main inputs | Main outputs |
|---|---|---|---|
| `save_theta_trans.R` | Optional utility to regenerate starting values | `albumin_data_0504.RData` | `theta_trans_3_0504.RData` |
| `albumin_fit.R` | Main real-data model fit | `albumin_data_0509.RData`, `theta_trans_3_0504.RData`, `theta_trans_1_0404.RData` | `albumin_fit.RData`, `theta_est_0613.RData`, population trend output |
| `bootstrap_ci.R` | Real-data bootstrap confidence intervals | Same inputs as `albumin_fit.R`, plus `theta_est_0613.RData` and `boostrap_realdata_0613.RData` | Parameter CI table, covariance CI table, correlation-network figure |
| `albumin_plot.R` | Main real-data trajectory figure | `albumin_data_0509.RData`, `theta_fit_theta_result_3.RData`, `theta_fit_fit_joint_albumin_8.RData`, `restandardize_for_plot_0505.RData` | Figure of individual trajectories and six-month-ahead mortality probability |
| `model_comparison_synth.R` | Sensitivity analysis comparing OU, cubic-cubic, and cubic-local models | Synthetic training/test data and likelihood/filter functions | `fit_ou_synth.RData`, `fit_cc_synth.RData`, `fit_cl_synth.RData`, `model_comparison_results.RData`, `model_comparison_plots.pdf` |
| `bootstrap_analysis.R` | Exact replication of simulation MSE/bias/variance and bootstrap coverage table | `theta_est_mat.Rdata`, `true_val.Rdata`, `theta_mat_est_extra_final.Rdata` | `par_est_df.Rdata`, `par_est_boot_df.Rdata`, LaTeX tables printed to console |
| `bootstrap_generate.R` | Documents the standardized bootstrap-generation procedure | Simulation data lists and `llh_survival_x_OU.R` | New stochastic `theta_est_mat*.Rdata` files |
| `bootstrap_generate_exact100_from_archive.R` | Recreates the exact 100 bootstrap input files from the archived final object | `theta_mat_est_extra_final.Rdata` | `bootstrap_results/results_m1000/theta_est_mat1.Rdata` through `theta_est_mat100.Rdata` |
| `bootstrap_collect.R` | Collects 100 bootstrap estimate matrices into one final object | `bootstrap_results/results_m1000/theta_est_mat*.Rdata` | `theta_mat_est_extra_final.Rdata` |
| `filter_time_comparison.R` | Filtering-time comparison for the new algorithm and univariate treatment | `llh_OU.R`, `llh_OU_unitrt.R` | `filter_time.Rdata`, `filter_time_unitrt.Rdata`, timing table, timing plot |

Library files are sourced by the analysis scripts and do not need to be run directly.

---

## Manuscript Output Map

| Manuscript output | Script to run | Required data objects |
|---|---|---|
| Population trajectory figure | `albumin_fit.R`, with bootstrap output available if confidence intervals are shown | `albumin_data_0509.RData`, `theta_trans_3_0504.RData`, `theta_trans_1_0404.RData` |
| Individual trajectory and mortality-probability figure | `albumin_plot.R` | `albumin_data_0509.RData`, `theta_fit_theta_result_3.RData`, `theta_fit_fit_joint_albumin_8.RData`, `restandardize_for_plot_0505.RData` |
| Real-data parameter confidence-interval table | `bootstrap_ci.R` | `theta_est_0613.RData`, `boostrap_realdata_0613.RData`, real-data input files |
| Real-data covariance / correlation summaries | `bootstrap_ci.R` | `boostrap_realdata_0613.RData` |
| Correlation-network figure | `bootstrap_ci.R` | `boostrap_realdata_0613.RData` |
| Model-comparison / sensitivity-analysis table and plots | `model_comparison_synth.R` | `synth_data_train_cubic_OU.RData`, `synth_data_test_cubic_OU.RData`, likelihood/filter helper files |
| Simulation MSE, variance, and bias table | `bootstrap_analysis.R` | `theta_est_mat.Rdata`, `true_val.Rdata` |
| Simulation bootstrap coverage table | `bootstrap_analysis.R` | `theta_est_mat.Rdata`, `true_val.Rdata`, `theta_mat_est_extra_final.Rdata` |
| Filtering-time comparison table and figure | `filter_time_comparison.R` | `llh_OU.R`, `llh_OU_unitrt.R` |

Table and figure numbering may change across manuscript drafts. The script descriptions above identify the corresponding outputs by content.

---

## Recommended Directory Structure

All scripts assume that the working directory is the repository root.

Recommended structure:

```text
.
├── README.md
├── save_theta_trans.R
├── albumin_fit.R
├── bootstrap_ci.R
├── albumin_plot.R
├── model_comparison_synth.R
├── bootstrap_analysis.R
├── bootstrap_generate.R
├── bootstrap_generate_exact100_from_archive.R
├── bootstrap_collect.R
├── filter_time_comparison.R
├── theta_mat_est_extra_final.Rdata
├── R/
│   ├── llh_survival_x_OU.R
│   ├── filter_estimate_OU_withP.R
│   ├── delta_llh_OU.R
│   ├── delta_negllh_OU.R
│   ├── llh_cs_cs.R
│   ├── llh_cs_local.R
│   ├── llh_joint_all_models.R
│   ├── data_preparation_functions.R
│   ├── helper_functions_losses.R
│   ├── llh_OU.R
│   └── llh_OU_unitrt.R
├── data/
│   └── results_m1000/
│       ├── theta_est_mat.Rdata
│       ├── true_val.Rdata
│       ├── MSE_fix_mat.Rdata
│       ├── y_list.Rdata
│       ├── z_list.Rdata
│       ├── X00_cube_in_list.Rdata
│       ├── A_list.Rdata
│       └── AB_list.Rdata
├── bootstrap_results/
│   └── results_m1000/
├── analysis_results/
└── filter_time_results/
```

Large `.Rdata` files may be distributed separately as release assets or through a restricted data link. After downloading them, place them in the paths shown above before running the scripts.

---

## System Requirements

### R version

R 4.4.1 or later is recommended.

### Required R packages

```r
install.packages(c(
  # Rcpp compilation
  "Rcpp", "RcppArmadillo", "inline",

  # Data handling
  "data.table", "dplyr", "reshape2",

  # Optimization and smoothing
  "GA", "assist",

  # Parallel computing
  "parallel", "foreach", "doParallel",

  # Visualization
  "ggplot2", "gridExtra", "ggpubr",

  # Survival analysis
  "survival", "survminer",

  # Tables
  "xtable",

  # Model comparison
  "pROC",

  # Simulation / filtering-time comparison
  "mvtnorm", "clusterGeneration", "magic"
))
```

### Compilation note

Several likelihood and filtering files use `Rcpp` / `inline` to compile C++ code at runtime through `source()`. Compilation output may be printed to the console and can be ignored.

A working C++ compiler is required:

- **Windows:** Rtools
- **Mac:** Xcode Command Line Tools and GFortran
- **Linux:** `gcc` and `gfortran`

---

## Quick Start for Exact Replication of Simulation Tables

To reproduce the simulation MSE/bias/variance table and the bootstrap coverage table, run:

```r
source("bootstrap_analysis.R")
```

This script reads:

```text
data/results_m1000/theta_est_mat.Rdata
data/results_m1000/true_val.Rdata
theta_mat_est_extra_final.Rdata
```

and saves:

```text
analysis_results/par_est_df.Rdata
analysis_results/par_est_boot_df.Rdata
```

It also prints the LaTeX tables to the console.

The file:

```text
theta_mat_est_extra_final.Rdata
```

is the archived final bootstrap object used to compute the reported bootstrap coverage probabilities.

---

## Bootstrap Replication Details

### Exact bootstrap object used in the paper

The reported bootstrap coverage table is based on:

```text
theta_mat_est_extra_final.Rdata
```

This object is a list of length 100. Each element is a matrix of bootstrap parameter estimates for one simulation replicate.

For exact numerical replication of the reported bootstrap coverage table, users should run:

```r
source("bootstrap_analysis.R")
```

directly from `theta_mat_est_extra_final.Rdata`.

### Recreating the 100 individual bootstrap RData files

To recreate the exact 100 individual bootstrap input files from the archived final object, run:

```r
source("bootstrap_generate_exact100_from_archive.R")
```

This creates:

```text
bootstrap_results/results_m1000/theta_est_mat1.Rdata
bootstrap_results/results_m1000/theta_est_mat2.Rdata
...
bootstrap_results/results_m1000/theta_est_mat100.Rdata
```

Each file contains one object named:

```text
theta_est_mat1
```

Then run:

```r
source("bootstrap_collect.R")
```

to collect those 100 files back into:

```text
theta_mat_est_extra_final.Rdata
```

### Standardized bootstrap-generation script

The script:

```text
bootstrap_generate.R
```

documents the standardized bootstrap-generation procedure used in the simulation study. It performs subject-level resampling, reconstructs the at-risk lists, fits the joint model by likelihood optimization, transforms the first `4*q` parameters back from the log scale, and saves bootstrap estimate matrices.

A fresh run of `bootstrap_generate.R` is stochastic and should not be expected to recreate the exact historical 100 bootstrap objects bit-for-bit. The historical bootstrap computations were generated across multiple runs and machines before the final archived set of 100 bootstrap objects was saved. For exact reproduction of the reported bootstrap coverage results, use `theta_mat_est_extra_final.Rdata` or the 100 files recreated from it by `bootstrap_generate_exact100_from_archive.R`.

---

## Main Real-Data Analysis Pipeline

### Step 0: `save_theta_trans.R` optional utility

Run only if `theta_trans_3_0504.RData` is missing.

**Inputs**

```text
albumin_data_0504.RData
```

**Outputs**

```text
theta_trans_3_0504.RData
theta_trans_est_0504.RData
```

### Step 1: `albumin_fit.R`

Fits the joint model to the real dialysis data.

**Inputs**

```text
albumin_data_0509.RData
theta_trans_3_0504.RData
theta_trans_1_0404.RData
```

**Outputs**

```text
albumin_fit.RData
theta_est_0613.RData
```

This step produces fitted model quantities used for the population trajectory figure.

Approximate runtime: 3 hours on a standard desktop.

### Step 2: `bootstrap_ci.R`

Computes nonparametric bootstrap confidence intervals.

**Inputs**

```text
albumin_data_0509.RData
theta_est_0613.RData
boostrap_realdata_0613.RData
```

and the same data inputs required by `albumin_fit.R`.

**Outputs**

```text
theta_result_k.RData
```

The script prints parameter confidence-interval tables and covariance/correlation summaries to the console. It also produces the correlation-network figure from the precomputed bootstrap object.

Approximate runtime: 18 hours per bootstrap batch if rerun from scratch. The precomputed file `boostrap_realdata_0613.RData` allows users to skip the long bootstrap-fitting step and proceed directly to table and figure generation.

### Step 3: `albumin_plot.R`

Generates the individual trajectory and mortality-probability figure.

**Inputs**

```text
albumin_data_0509.RData
theta_fit_theta_result_3.RData
theta_fit_fit_joint_albumin_8.RData
restandardize_for_plot_0505.RData
```

**Outputs**

```text
0509_6pt_3156.png
0509_6pt_130_noCI_0904.pdf
```

Approximate runtime: 10 minutes.

---

## Sensitivity Analysis Pipeline

### `model_comparison_synth.R`

Fits three model specifications on synthetic data:

1. cubic spline population + OU individual process,
2. cubic spline population + cubic spline individual process,
3. cubic spline population + local-level individual process.

**Inputs**

```text
synth_data_train_cubic_OU.RData
synth_data_test_cubic_OU.RData
theta_trans_3_0504.RData
theta_trans_3_0503.RData
theta_trans_3_0502.RData
theta_trans_1_0404.RData
```

The script also sources:

```text
data_preparation_functions.R
helper_functions_losses.R
llh_survival_x_OU.R
filter_estimate_OU_withP.R
llh_cs_cs.R
llh_cs_local.R
llh_joint_all_models.R
```

**Outputs**

```text
fit_ou_synth.RData
fit_cc_synth.RData
fit_cl_synth.RData
model_comparison_results.RData
master_comparison_synth.RData
model_comparison_plots.pdf
```

Approximate runtime: 12 hours total.

---

## Filtering-Time Comparison Pipeline

### `filter_time_comparison.R`

Compares filtering time for the proposed Kalman-filter algorithm and the univariate-treatment alternative.

**Inputs**

```text
R/llh_OU.R
R/llh_OU_unitrt.R
```

or the same files in the repository root.

**Outputs**

```text
filter_time_results/filter_time.Rdata
filter_time_results/filter_time_unitrt.Rdata
```

The script also prints a LaTeX timing table and displays the filtering-time comparison plot.

---

## Core Likelihood and Filter Functions

These files are sourced automatically by the analysis scripts.

| File | Description |
|---|---|
| `llh_survival_x_OU.R` | Joint log-likelihood for the OU longitudinal process and survival model. Used by the real-data fit, bootstrap, and simulation bootstrap-generation scripts. |
| `filter_estimate_OU_withP.R` | Kalman filter and smoother for the OU model. Returns filtered/smoothed population and individual trajectories and posterior covariance quantities. |
| `delta_llh_OU.R` | Computes the Rosenberg `delta_est` and `neg2llh` for longitudinal-only OU initialization. |
| `delta_negllh_OU.R` | Returns `neg2llh` only. Used with `optim()` for univariate initialization. |
| `llh_cs_cs.R` | Longitudinal log-likelihood for the cubic spline population + cubic spline individual model. |
| `llh_cs_local.R` | Longitudinal log-likelihood for the cubic spline population + local-level individual model. |
| `llh_joint_all_models.R` | Joint likelihoods for the cubic-cubic and cubic-local sensitivity-analysis models. |
| `data_preparation_functions.R` | Helper functions for constructing arrays and lists from long-format synthetic data. |
| `helper_functions_losses.R` | Loss metrics and prediction metrics for model comparison, including MSE, MAE, NLL, Brier score, AUC, C-index, and risk-score computation. |
| `llh_OU.R` | Likelihood used in the filtering-time comparison for the proposed algorithm. |
| `llh_OU_unitrt.R` | Likelihood used in the filtering-time comparison for the univariate-treatment method. |

---

## Data Objects

All `.RData` and `.Rdata` extensions are case-sensitive on Mac/Linux.

### Main real-data objects

| RData file | Used by | Contains |
|---|---|---|
| `albumin_data_0504.RData` | `save_theta_trans.R` | Preliminary optimization object and starting vector |
| `albumin_data_0509.RData` | `albumin_fit.R`, `albumin_plot.R` | Preprocessed dialysis cohort data |
| `data_split_in.Rdata` | `albumin_fit.R`, `bootstrap_ci.R` | q × npt × m biomarker array |
| `z_in.Rdata` | `albumin_fit.R`, `bootstrap_ci.R` | m × npt survival indicator matrix |
| `X00_cube_in.Rdata` | `albumin_fit.R`, `bootstrap_ci.R` | m × p × npt covariate array |
| `A_list_in.Rdata` | `albumin_fit.R`, `bootstrap_ci.R` | List of alive-at-time subject indices |
| `AB_list_in.Rdata` | `albumin_fit.R`, `bootstrap_ci.R` | List of present-at-time subject indices |
| `theta_trans_3_0504.RData` | `albumin_fit.R`, `bootstrap_ci.R`, `model_comparison_synth.R` | Starting parameter vector |
| `theta_trans_3_0503.RData` | `albumin_fit.R`, `bootstrap_ci.R`, `model_comparison_synth.R` | Alternative starting parameter vector |
| `theta_trans_3_0502.RData` | `albumin_fit.R`, `bootstrap_ci.R`, `model_comparison_synth.R` | Alternative starting parameter vector |
| `theta_trans_1_0404.RData` | `albumin_fit.R`, `bootstrap_ci.R`, `model_comparison_synth.R` | Starting values for survival subvector |
| `theta_est_0613.RData` | `bootstrap_ci.R` | Fitted natural-scale parameter vector |
| `boostrap_realdata_0613.RData` | `bootstrap_ci.R` | Precomputed real-data bootstrap parameter vectors |
| `restandardize_for_plot_0505.RData` | `albumin_plot.R` | Biomarker means and standard deviations for back-transformation |
| `theta_fit_theta_result_3.RData` | `albumin_plot.R` | Bootstrap parameter list; point estimate in `[[2]]` |
| `theta_fit_fit_joint_albumin_8.RData` | `albumin_plot.R` | Fitted `optim()` object |
| `synth_data_train_cubic_OU.RData` | `model_comparison_synth.R` | Synthetic training data |
| `synth_data_test_cubic_OU.RData` | `model_comparison_synth.R` | Synthetic test data |
| `model_comparison_results.RData` | Review / verification | Precomputed model-comparison results |

### Simulation/bootstrap validation objects

| RData file | Used by | Contains |
|---|---|---|
| `data/results_m1000/theta_est_mat.Rdata` | `bootstrap_analysis.R`, `bootstrap_generate.R` | Parameter estimates from 100 simulation runs |
| `data/results_m1000/true_val.Rdata` | `bootstrap_analysis.R` | True parameter vector |
| `theta_mat_est_extra_final.Rdata` | `bootstrap_analysis.R`, `bootstrap_generate_exact100_from_archive.R` | Archived final list of 100 bootstrap estimate matrices used for reported coverage |
| `data/results_m1000/MSE_fix_mat.Rdata` | `bootstrap_generate.R` | MSE matrix used to select the simulation dataset for bootstrap resampling |
| `data/results_m1000/y_list.Rdata` | `bootstrap_generate.R` | Simulated longitudinal biomarker arrays |
| `data/results_m1000/z_list.Rdata` | `bootstrap_generate.R` | Simulated survival indicator matrices |
| `data/results_m1000/X00_cube_in_list.Rdata` | `bootstrap_generate.R` | Simulated covariate arrays |
| `data/results_m1000/A_list.Rdata` | `bootstrap_generate.R` | At-risk lists for simulated data |
| `data/results_m1000/AB_list.Rdata` | `bootstrap_generate.R` | Present-at-time lists for simulated data |

`X00_cube_in_list.Rdata` is required only for rerunning `bootstrap_generate.R`. It is not required for exact replication of the reported simulation/bootstrap coverage table if `theta_mat_est_extra_final.Rdata` is provided.

### Filtering-time objects

| RData file | Used by | Contains |
|---|---|---|
| `filter_time_results/filter_time.Rdata` | `filter_time_comparison.R` | Filtering times for the proposed algorithm |
| `filter_time_results/filter_time_unitrt.Rdata` | `filter_time_comparison.R` | Filtering times or extrapolated filtering times for the univariate-treatment method |

---

## Primary Real-Data Object Description

### `albumin_data_0509.RData`

Contains the preprocessed dialysis cohort with 5 biomarkers and 48 monthly visits.

| Object | Type | Dimensions | Description |
|---|---|---|---|
| `data_split_in` | array | 5 × 48 × 5,707 | Biomarker observations. Biomarkers: albumin, interdialytic weight gain, systolic blood pressure, neutrophil-to-lymphocyte ratio, and hemoglobin. |
| `z_in` | matrix | 5,707 × 48 | Survival indicators. `z_in[i, t] = 1` if subject `i` died at month `t`, `0` if alive, and `NA` if censored before month `t`. |
| `X00_cube_in` | array | 5,707 × 6 × 48 | Covariate array. Covariates: age, sex, EKTV, diabetes, cardiovascular disease, and vintage. |

Derived index lists are reconstructed from `z_in` in the scripts.

| Object | Type | Length | Description |
|---|---|---|---|
| `A_list_in` | list | 48 | Subject indices alive at each month. |
| `AB_list_in` | list | 48 | Subject indices present at each month. |

---

## Parameter Vector Structure

For q = 5 biomarkers and p = 6 covariates, the parameter vector has length 42.

The first `4*q` elements are optimized on the log scale and are transformed back using:

```r
theta_est[1:(4*q)] <- exp(theta_trans[1:(4*q)])
```

The simulation/bootstrap validation scripts use the following ordering:

| Index | Symbol | Description |
|---|---|---|
| 1–5 | ζ₁², ..., ζ₅² | Population smoothing variances |
| 6–10 | ν₁², ..., ν₅² | OU individual innovation variances |
| 11–15 | ξ₁, ..., ξ₅ | OU mean-reversion rates |
| 16–20 | d₁, ..., d₅ | Diagonal elements in the LDLᵀ covariance parameterization |
| 21–30 | l₁, ..., l₁₀ | Lower-triangular elements in the LDLᵀ covariance parameterization |
| 31–35 | γ₂₁, ..., γ₂₅ | Association parameters linking individual longitudinal state to survival |
| 36–42 | γ₁₀, γ₁₁, ..., γ₁₆ | Survival covariate coefficients, including intercept |

For reporting simulation results, the covariance parameters are transformed from the `(D, L)` parameterization to the entries of `Sigma_epsilon`.

The covariance is represented as:

```text
Sigma_epsilon = L D L^T
```

where `D = diag(d_1, ..., d_q)` and `L` is lower triangular with ones on the diagonal.

---

## Large Files and Data Availability

Some `.Rdata` files may be too large to upload through the GitHub web interface. In that case, they should be provided as GitHub Release assets, through a separate data archive, or through a restricted reviewer link.

For exact replication of the simulation/bootstrap coverage table, the essential large file is:

```text
theta_mat_est_extra_final.Rdata
```

For rerunning the stochastic bootstrap-generation procedure, the following large file may also be needed:

```text
data/results_m1000/X00_cube_in_list.Rdata
```

The dialysis cohort data are subject to data-use restrictions and may not be publicly deposited. Authorized users should obtain the restricted data bundle separately and place the files in the paths shown in this README.

---

## Reproducibility Notes

1. **Working directory:** Run scripts from the repository root.

2. **Exact bootstrap replication:** The reported simulation/bootstrap coverage table should be reproduced using `theta_mat_est_extra_final.Rdata`. A fresh run of `bootstrap_generate.R` documents the procedure but is not expected to recreate the exact archived bootstrap objects.

3. **Random seeds:** Seeds are set inside the scripts, but stochastic reruns may differ across platforms, R versions, and optimization environments.

4. **Long runtimes:** Real-data model fitting, real-data bootstrap fitting, synthetic model comparison, and stochastic bootstrap generation are computationally intensive. Precomputed RData files are included or provided separately to allow exact replication without rerunning all long jobs.

5. **Rcpp compilation:** C++ code may compile when likelihood/filter files are sourced. Compilation output printed to the console can be ignored.

6. **Case sensitivity:** `.RData` and `.Rdata` are distinct on Mac/Linux. Use the filenames exactly as shown.

---

## Minimal Files Needed by Replication Task

### Exact simulation/bootstrap table replication

```text
bootstrap_analysis.R
theta_mat_est_extra_final.Rdata
data/results_m1000/theta_est_mat.Rdata
data/results_m1000/true_val.Rdata
```

### Exact filtering-time table/figure replication

```text
filter_time_comparison.R
R/llh_OU.R
R/llh_OU_unitrt.R
```

### Real-data manuscript figures and tables

```text
albumin_fit.R
bootstrap_ci.R
albumin_plot.R
albumin_data_0509.RData
theta_trans_3_0504.RData
theta_trans_1_0404.RData
theta_est_0613.RData
boostrap_realdata_0613.RData
theta_fit_theta_result_3.RData
theta_fit_fit_joint_albumin_8.RData
restandardize_for_plot_0505.RData
```

plus the relevant likelihood/filter source files.

### Full bootstrap-generation documentation

```text
bootstrap_generate.R
R/llh_survival_x_OU.R
data/results_m1000/y_list.Rdata
data/results_m1000/z_list.Rdata
data/results_m1000/X00_cube_in_list.Rdata
data/results_m1000/A_list.Rdata
data/results_m1000/AB_list.Rdata
data/results_m1000/MSE_fix_mat.Rdata
data/results_m1000/theta_est_mat.Rdata
```

Again, these files are needed only to rerun the stochastic bootstrap-generation procedure, not to reproduce the final reported bootstrap coverage table.

---

## Recommended Replication Commands

Exact simulation/bootstrap table replication:

```r
source("bootstrap_analysis.R")
```

Recreate the 100 exact individual bootstrap input files from the archived final object:

```r
source("bootstrap_generate_exact100_from_archive.R")
source("bootstrap_collect.R")
```

Filtering-time comparison:

```r
source("filter_time_comparison.R")
```

Main real-data analysis, if authorized real-data files are available:

```r
source("save_theta_trans.R")   # optional, only if starting values are missing
source("albumin_fit.R")
source("bootstrap_ci.R")
source("albumin_plot.R")
```

Sensitivity analysis:

```r
source("model_comparison_synth.R")
```
