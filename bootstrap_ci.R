# =============================================================================
# Bootstrap Confidence Intervals for Joint Model Parameters
# =============================================================================
# Runs nonparametric bootstrap (114 resamples, 3 per batch, 38 batches) to
# compute confidence intervals for all parameters. Also produces:
#   - Table 5: parameter estimates with 95% bootstrap CIs
#   - Table 5 Sigma: cross-biomarker covariance/correlation matrix with CIs
#   - Figure 4: correlation network (random effects, significant edges only)
#
# Run order: 2 of 4
#   1. albumin_fit.R
#   2. bootstrap_ci.R          <- YOU ARE HERE
#   3. albumin_plot.R
#
# Inputs:
#   data_split_in.Rdata        # q x npt x m biomarker array
#   z_in.Rdata                 # m x npt survival indicator matrix
#   X00_cube_in.Rdata          # m x p x npt covariate array
#   theta_trans_3_0504.RData   # initial parameter vector, final restart
#   theta_trans_3_0503.RData   # initial parameter vector, 2nd restart (used in blending)
#   theta_trans_3_0502.RData   # initial parameter vector, 1st restart (used in blending)
#   theta_trans_1_0404.RData   # initial survival subvector
#   theta_est_0613.RData       # point estimates from albumin_fit.R
#   boostrap_realdata_0613.RData  # pre-computed bootstrap results (114 samples)
#
# Outputs (saved per batch k):
#   theta_result_k.Rdata       # list of 3 fitted parameter vectors for batch k
#   fit_joint_albumin_j.Rdata  # optim object for bootstrap sample j
#
# Parameter vector structure: see albumin_fit.R header for full description.
# =============================================================================

rm(list = ls())

# --- USER: set this to the directory containing all data and script files ---
ROOT_DIR <- "."
setwd(ROOT_DIR)

set.seed(1)

library(GA)
library(data.table)
library(dplyr)
library(assist)
library(parallel)
library(foreach)
library(doParallel)
library(xtable)

# =============================================================================
# PART 1: RUN BOOTSTRAP (38 batches x 3 samples = 114 total resamples)
# =============================================================================
# Each batch k fits the joint model on 3 bootstrap resamples in parallel.
# Runtime per batch: ~18 hours. Skip this block if theta_result_k.Rdata
# files already exist and proceed to Part 2.

for(k in 38){

  load("data_split_in.Rdata")   # data_split_in: q x npt x m biomarker array
  load("z_in.Rdata")            # z_in: m x npt survival indicator
  load("X00_cube_in.Rdata")     # X00_cube_in: m x p x npt covariates
  load("A_list_in.Rdata")       # (unused below; rebuilt from z_in after resampling)
  load("AB_list_in.Rdata")

  numCores <- detectCores(logical = FALSE) - 1
  registerDoParallel(numCores, 3)

  q   <- dim(data_split_in)[1]   # 5 biomarkers
  p   <- dim(X00_cube_in)[2]     # 6 covariates
  m   <- dim(data_split_in)[3]   # 5707 subjects
  npt <- dim(data_split_in)[2]   # 48 months

  n_sample <- 10   # samples per batch (only 3 used per foreach call)

  theta_est_mat1 <- matrix(0, nrow = 0.5 * q^2 + 4.5 * q + p + 1, ncol = n_sample)

  time_eclipse <- system.time({
    theta_result <- foreach(
      sample_id = c((3*(k-1)+1):(3*k)),
      .packages  = c("Rcpp", "inline"),
      .noexport  = c("<llh_survival_x_OU>", "<delta_llh_OU>",
                     "<delta_negllh_OU>", "<filter_estimate_OU>")
    ) %dopar% {

      source("delta_llh_OU.R")            # delta_est + neg2llh, longitudinal OU
      source("delta_negllh_OU.R")         # neg2llh only
      source("filter_estimate_OU_withP.R") # Kalman filter
      source("llh_survival_x_OU.R")       # joint likelihood

      # Standardize biomarkers and covariates (same as albumin_fit.R)
      for(i in 1:dim(data_split_in)[1]){
        data_split_in[i,,] <- (data_split_in[i,,] - mean(data_split_in[i,,], na.rm = TRUE)) /
                                sd(data_split_in[i,,], na.rm = TRUE)
      }
      for(i in c(1, 3)){
        X00_cube_in[,i,] <- (X00_cube_in[,i,] - mean(X00_cube_in[,i,], na.rm = TRUE)) /
                              sd(X00_cube_in[,i,], na.rm = TRUE)
      }

      # Bootstrap resample with replacement
      ids          <- sample(1:m, size = m, replace = TRUE)
      data_split_in <- data_split_in[,, ids]
      z_in          <- z_in[ids,]
      X00_cube_in   <- X00_cube_in[ids,,]

      # Rebuild A_list and AB_list from resampled z_in
      # A_list_in[[t]]: alive-at-t subject indices (used in Kalman filter)
      dummy_A  <- z_in; dummy_A[dummy_A == 1] <- NA
      A_list_in <- vector("list", npt)
      for(j in 1:npt) A_list_in[[j]] <- (1:m)[!is.na(dummy_A[, j])]

      # AB_list_in[[t]]: present-at-t subject indices (used in survival likelihood)
      dummy_AB   <- z_in
      AB_list_in <- vector("list", npt)
      for(j in 1:npt) AB_list_in[[j]] <- (1:m)[!is.na(dummy_AB[, j])]

      # Construct initial parameter vector (same blending as albumin_fit.R)
      load("theta_trans_3_0504.RData"); theta_trans_1 <- theta_trans_3
      load("theta_trans_3_0503.RData"); theta_trans_2 <- theta_trans_3
      load("theta_trans_3_0502.RData"); theta_trans_4 <- theta_trans_3
      rm(theta_trans_3); theta_trans_3 <- theta_trans_1
      ind <- 3
      theta_trans_3[c(0+ind, 5+ind, 10+ind, 15+ind)] <-
        theta_trans_2[c(0+ind, 5+ind, 10+ind, 15+ind)]
      ind <- 5
      theta_trans_3[c(0+ind, 5+ind, 10+ind, 15+ind)] <-
        theta_trans_4[c(0+ind, 5+ind, 10+ind, 15+ind)]
      load("theta_trans_1_0404.RData")
      theta_trans_3[31:42] <- theta_trans_1[31:42]

      # D_in: maps combined state [alpha_u (2q), alpha_v (q)] to observations
      D_in <- cbind(diag(q) %x% matrix(c(1, 0), nrow = 1, ncol = 2), diag(q))

      fit_joint_albumin <- optim(
        theta_trans_3, llh_survival_x_OU,
        data_split_in = data_split_in, z_in = z_in,
        D_in = D_in, X00_cube_in = X00_cube_in,
        A_list_in = A_list_in, AB_list_in = AB_list_in,
        hessian = FALSE
      )

      save(fit_joint_albumin,
           file = paste0("fit_joint_albumin_", sample_id, ".Rdata"))

      # Back-transform to natural scale
      theta_trans_est <- fit_joint_albumin$par
      theta_est        <- theta_trans_est
      theta_est[1:(4 * q)] <- exp(theta_trans_est[1:(4 * q)])

      # Extract named parameter subvectors (for downstream use)
      zeta_vec  <- theta_est[1:q]
      nu2_vec   <- theta_est[(q+1):(2*q)]
      xi_vec    <- theta_est[(2*q+1):(3*q)]
      d_vec     <- theta_est[(3*q+1):(4*q)]
      l_vec     <- theta_est[(4*q+1):(0.5*q^2 + 3.5*q)]
      gamma_vec <- theta_est[(0.5*q^2 + 3.5*q + 1):(0.5*q^2 + 4.5*q)]
      gamma_x   <- theta_est[(0.5*q^2 + 4.5*q + 1):(0.5*q^2 + 4.5*q + p + 1)]

      theta_est_mat1[, sample_id - 4*(k-1)] <- theta_est
    }

    save(theta_result, file = paste0("theta_result_", k, ".Rdata"))
  })

  save(time_eclipse, file = paste0("time_eclipse", k, ".Rdata"))
  rm(A_list_in, AB_list_in, dummy_A, dummy_AB, theta_est_mat1,
     theta_result, z_in, data_split_in, ids, j, m, n_sample,
     npt, numCores, p, q, time_eclipse, X00_cube_in)
}

# =============================================================================
# PART 2: AGGREGATE BOOTSTRAP RESULTS -> TABLE 5 (parameter CIs)
# =============================================================================
# Loads all 38 theta_result_k.Rdata files, assembles 42 x 114 matrix,
# computes bootstrap SE and 2.5/97.5 percentile CIs.

rm(list = ls())
library(xtable)

# Load all 38 batch results (3 samples each = 114 total)
for(k in 1:38){
  res_name <- paste0("theta_result_", k)
  load(paste0("theta_result_", k, ".Rdata"))
  assign(paste0("theta_result_", k), theta_result)
  rm(theta_result)
}

# Assemble full bootstrap matrix: 42 parameters x 114 samples
theta_mat_boot <- do.call(cbind, lapply(1:38, function(k){
  res <- get(paste0("theta_result_", k))
  matrix(unlist(res), ncol = length(res), nrow = 42)
}))

se_theta <- apply(theta_mat_boot, 1, sd)
lwr      <- apply(theta_mat_boot, 1, quantile, probs = 0.025)
upr      <- apply(theta_mat_boot, 1, quantile, probs = 0.975)

# Load point estimates from albumin_fit.R
load("theta_est_0613.RData")   # theta_est: 42-vector on natural scale

# Parameter names matching paper notation (Section 3, Table 5)
# zeta: population smoothing variances; nu2: OU innovation variances;
# xi: OU mean-reversion rates; d/l: Cholesky of cross-biomarker covariance Sigma;
# gamma_vec: individual state -> survival link; gamma_x: survival covariates
parameter <- c(
  "zeta_1^2", "zeta_2^2", "zeta_3^2", "zeta_4^2", "zeta_5^2",
  "nu_1^2",   "nu_2^2",   "nu_3^2",   "nu_4^2",   "nu_5^2",
  "xi_1",     "xi_2",     "xi_3",     "xi_4",     "xi_5",
  "d_1",      "d_2",      "d_3",      "d_4",      "d_5",
  "l_1",      "l_2",      "l_3",      "l_4",      "l_5",
  "l_6",      "l_7",      "l_8",      "l_9",      "l_10",
  "gamma_{21}", "gamma_{22}", "gamma_{23}", "gamma_{24}", "gamma_{25}",
  "gamma_{0}",  "gamma_{11}", "gamma_{12}", "gamma_{13}", "gamma_{14}",
  "gamma_{15}", "gamma_{16}"
)

par_est_df <- data.frame(estimate = theta_est, lwr = lwr, upr = upr)
rownames(par_est_df) <- parameter

print(xtable(par_est_df, digits = c(0, 4, 4, 4)))
round(par_est_df, 4)

# =============================================================================
# PART 3: SIGMA MATRIX BOOTSTRAP CIs (cross-biomarker covariance, Table 5 Sigma)
# =============================================================================
q <- 5
Sig_boot_cube  <- array(0, dim = c(q, q, 114))
Corr_boot_cube <- array(0, dim = c(q, q, 114))

for(i in 1:114){
  theta_est_here <- theta_mat_boot[, i]
  d_vec <- theta_est_here[(3*q+1):(4*q)]
  l_vec <- theta_est_here[(4*q+1):(0.5*q^2 + 3.5*q)]
  L <- diag(q); L[lower.tri(L)] <- l_vec
  D <- diag(d_vec, nrow = q, ncol = q)
  Sig_est <- L %*% D %*% t(L)
  Sig_boot_cube[,, i] <- Sig_est
  for(j in 1:q){
    Corr_boot_cube[j, j, i] <- 1
    for(kk in 1:q){
      if(kk != j){
        Corr_boot_cube[j, kk, i] <- Sig_boot_cube[j, kk, i] /
          (sqrt(Sig_boot_cube[j, j, i]) * sqrt(Sig_boot_cube[kk, kk, i]))
      }
    }
  }
}

# =============================================================================
# PART 4: FIGURE 4 — CORRELATION NETWORK (random effects)
# =============================================================================
# Loads pre-computed bootstrap results (boostrap_realdata_0613.RData) and
# plots the partial correlation graph among the 5 biomarker random effects.
# Only edges with bootstrap CIs excluding zero are shown (significant associations).
# Interpretation: correlations among INDIVIDUAL departures from population trends,
# after removing fixed covariate effects.

load("boostrap_realdata_0613.RData")
