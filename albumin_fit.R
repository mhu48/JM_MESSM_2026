# =============================================================================
# Joint Model Fitting: Cubic Spline Population + OU Individual Trajectories
# =============================================================================
# Fits the joint model (Section 3 of paper) to the dialysis cohort.
# Loads pre-processed data arrays, standardizes biomarkers and covariates,
# constructs initial parameter values, and runs numerical optimization.
#
# Run order: 1 of 4
#   1. albumin_fit.R          <- YOU ARE HERE
#   2. bootstrap_ci.R
#   3. albumin_plot.R
#
# Inputs:
#   data_split_in.Rdata       # q x npt x m biomarker array (q=5, npt=48, m=5707)
#   z_in.Rdata                # m x npt survival indicator matrix (1 = died at month t)
#   X00_cube_in.Rdata         # m x p x npt covariate array (p=6)
#   A_list_in.Rdata           # list of length npt; A_list_in[[t]] = subject indices alive at t
#   AB_list_in.Rdata          # list of length npt; AB_list_in[[t]] = subject indices present at t
#   theta_trans_3_0504.RData  # initial parameter vector (log-transformed), final restart
#   theta_trans_3_0503.RData  # initial parameter vector, 2nd restart (used in blending)
#   theta_trans_3_0502.RData  # initial parameter vector, 1st restart (used in blending)
#   theta_trans_1_0404.RData  # initial survival subvector from earliest preliminary fit
#
# Outputs:
#   albumin_fit.RData         # full workspace including fit_joint_albumin_1 (optim object)
#   theta_est_0613.RData      # theta_est: fitted parameter vector (natural scale)
#
# Parameter vector structure (theta_trans, log-transformed):
#   [1:q]          log(zeta_vec)   population smoothing variances, one per biomarker
#   [q+1:2q]       log(nu2_vec)    individual OU innovation variances
#   [2q+1:3q]      log(xi_vec)     OU mean-reversion rates
#   [3q+1:4q]      log(d_vec)      diagonal of Cholesky factor of Sigma (cross-biomarker covariance)
#   [4q+1:...]     l_vec           lower-triangular elements of Cholesky factor of Sigma
#   [n_long]       log(kappa)      survival baseline hazard scale
#   [n_long+1:     gamma_vec       association parameters linking individual state to survival
#    n_long+q]
#   [n_long+q+1:   gamma_x         survival covariate coefficients (intercept + p covariates)
#    n_long+q+p+1]
# =============================================================================

rm(list = ls())

# --- USER: set this to the directory containing all data and script files ---
ROOT_DIR <- "."
setwd(ROOT_DIR)

library(data.table)
library(dplyr)
library(assist)
library(GA)
library(xtable)
library(ggplot2)

# =============================================================================
# 1. LOAD DATA
# =============================================================================

load("data_split_in.Rdata")   # data_split_in: q x npt x m array of biomarker observations
load("z_in.Rdata")            # z_in: m x npt matrix; z_in[i,t] = 1 if subject i died at month t
load("X00_cube_in.Rdata")     # X00_cube_in: m x p x npt covariate array
load("A_list_in.Rdata")       # A_list_in: list of length npt; alive-at-t subject indices
load("AB_list_in.Rdata")      # AB_list_in: list of length npt; present-at-t subject indices

# =============================================================================
# 2. COMPUTE AND PRINT BIOMARKER SUMMARY STATISTICS (pre-standardization)
# =============================================================================
# These means and SDs describe the raw biomarker distributions and appear in Table 1.

mean_alb = mean_idwg = mean_sbp = mean_nlr = mean_hgb = c()
sd_alb = sd_idwg = sd_sbp = sd_nlr = sd_hgb = c()
for(i in 1:dim(data_split_in)[3]){
  mean_alb   = c(mean_alb,   mean(data_split_in[1,,i], na.rm = TRUE))
  mean_idwg  = c(mean_idwg,  mean(data_split_in[2,,i], na.rm = TRUE))
  mean_sbp   = c(mean_sbp,   mean(data_split_in[3,,i], na.rm = TRUE))
  mean_nlr   = c(mean_nlr,   mean(data_split_in[4,,i], na.rm = TRUE))
  mean_hgb   = c(mean_hgb,   mean(data_split_in[5,,i], na.rm = TRUE))
  sd_alb     = c(sd_alb,     sd(data_split_in[1,,i],   na.rm = TRUE))
  sd_idwg    = c(sd_idwg,    sd(data_split_in[2,,i],   na.rm = TRUE))
  sd_sbp     = c(sd_sbp,     sd(data_split_in[3,,i],   na.rm = TRUE))
  sd_nlr     = c(sd_nlr,     sd(data_split_in[4,,i],   na.rm = TRUE))
  sd_hgb     = c(sd_hgb,     sd(data_split_in[5,,i],   na.rm = TRUE))
}

mean_alb_avg  = mean(mean_alb);  sd_alb_avg  = mean(sd_alb)
mean_idwg_avg = mean(mean_idwg); sd_idwg_avg = mean(sd_idwg)
mean_sbp_avg  = mean(mean_sbp);  sd_sbp_avg  = mean(sd_sbp)
mean_nlr_avg  = mean(mean_nlr);  sd_nlr_avg  = mean(sd_nlr)
mean_hgb_avg  = mean(mean_hgb);  sd_hgb_avg  = mean(sd_hgb)

data_summary <- cbind(
  variables = c("ALB", "IDWG", "SBP", "NLR", "HGB"),
  mean      = round(c(mean_alb_avg, mean_idwg_avg, mean_sbp_avg, mean_nlr_avg, mean_hgb_avg), 2),
  stdvar    = round(c(sd_alb_avg,   sd_idwg_avg,   sd_sbp_avg,   sd_nlr_avg,   sd_hgb_avg),   2)
)
data_summary <- as.data.frame(data_summary)
colnames(data_summary) <- c("Clinical Variables", "Average of Means", "Average of SDs")
print(xtable(data_summary, include.rownames = FALSE))

# =============================================================================
# 3. SOURCE LIKELIHOOD AND FILTER FUNCTIONS
# =============================================================================

source("delta_llh_OU.R")           # delta_est + neg2llh for longitudinal-only OU model
source("delta_negllh_OU.R")        # neg2llh only, for optim-based univariate fitting
source("filter_estimate_OU_withP.R") # Kalman filter/smoother returning f_filter, b_filter, P matrices
source("llh_survival_x_OU.R")      # joint log-likelihood: longitudinal (OU) + survival

# =============================================================================
# 4. STANDARDIZE BIOMARKERS AND COVARIATES
# =============================================================================
# Biomarkers (all 5): standardize to mean 0, SD 1 across all subjects and times.
# Covariates: standardize age (col 1) and EKTV (col 3) only; binary covariates left as-is.
# Standardization constants are saved in restandardize_for_plot_0505.RData for back-transformation.

for(i in 1:dim(data_split_in)[1]){
  data_split_in[i,,] <- (data_split_in[i,,] - mean(data_split_in[i,,], na.rm = TRUE)) /
                          sd(data_split_in[i,,], na.rm = TRUE)
}

for(i in c(1, 3)){
  X00_cube_in[,i,] <- (X00_cube_in[,i,] - mean(X00_cube_in[,i,], na.rm = TRUE)) /
                        sd(X00_cube_in[,i,], na.rm = TRUE)
}

# =============================================================================
# 5. SET DIMENSIONS AND REBUILD INDEX LISTS FROM z_in
# =============================================================================
# Note: uses all m = 5707 subjects (idx = 1:m is not a subset selection)

m   <- dim(data_split_in)[3]   # number of subjects
q   <- dim(data_split_in)[1]   # number of biomarkers (5)
npt <- dim(data_split_in)[2]   # number of monthly time points (48)
p   <- dim(X00_cube_in)[2]     # number of survival covariates (6)

# A_list_in[[t]]: indices of subjects alive (not yet dead) at month t
# Used in Kalman filter to compute population curve f_t
dummy_A <- z_in
dummy_A[dummy_A == 1] <- NA
A_list_in <- vector("list", npt)
for(j in 1:npt){
  A_list_in[[j]] <- (1:m)[!is.na(dummy_A[, j])]
}

# AB_list_in[[t]]: indices of subjects present at month t (alive or dying at t)
# Used in survival likelihood to include the death event at t
dummy_AB <- z_in
AB_list_in <- vector("list", npt)
for(j in 1:npt){
  AB_list_in[[j]] <- (1:m)[!is.na(dummy_AB[, j])]
}

# =============================================================================
# 6. CONSTRUCT INITIAL PARAMETER VECTOR
# =============================================================================
# theta_trans_3_0504 is the result of sequential restarts (0502 -> 0503 -> 0504)
# converging toward a good starting point for the full joint optimization.
# theta_trans_1_0404 provides better initial values for the survival subvector (indices 31:42).

set.seed(240)
load("theta_trans_3_0504.RData")   # loads theta_trans_3
theta_trans_1 <- theta_trans_3
load("theta_trans_3_0503.RData")   # loads theta_trans_3
theta_trans_2 <- theta_trans_3
load("theta_trans_3_0502.RData")   # loads theta_trans_3
theta_trans_4 <- theta_trans_3
rm(theta_trans_3)
theta_trans_3 <- theta_trans_1

# Blend: take NLR (index 3) and HGB (index 5) entries from alternate restarts
# to improve starting values for those biomarkers
ind <- 3
theta_trans_3[c(0+ind, 5+ind, 10+ind, 15+ind)] <-
  theta_trans_2[c(0+ind, 5+ind, 10+ind, 15+ind)]
ind <- 5
theta_trans_3[c(0+ind, 5+ind, 10+ind, 15+ind)] <-
  theta_trans_4[c(0+ind, 5+ind, 10+ind, 15+ind)]

# Replace survival subvector (gamma_vec + gamma_x, indices 31:42) with 0404 values
load("theta_trans_1_0404.RData")   # loads theta_trans_1
theta_trans_3[31:42] <- theta_trans_1[31:42]

# Longitudinal subvector only (indices 1:30) for filter initialization checks
theta_trans_long_0 <- theta_trans_3[1:30]

# =============================================================================
# 7. FIT JOINT MODEL
# =============================================================================
# D_in maps the combined state vector [alpha_u (2q), alpha_v (q)] to the
# observation equation. For OU individual process: D_in selects position from
# alpha_u (cols 1,3,5,...) and the full alpha_v.
# See Section 3.2 of paper for state-space formulation.

D_in <- cbind(diag(q) %x% matrix(c(1, 0), nrow = 1, ncol = 2), diag(q))

# Evaluate initial log-likelihood as sanity check before optimization
llh_survival_x_OU(theta_trans_3, data_split_in, z_in = z_in,
                  D_in = D_in, X00_cube_in = X00_cube_in,
                  A_list_in = A_list_in, AB_list_in = AB_list_in)

# Run Nelder-Mead optimization (~3 hours runtime)
time_here <- system.time(
  fit_joint_albumin_1 <- optim(
    theta_trans_3, llh_survival_x_OU,
    data_split_in  = data_split_in, z_in = z_in,
    D_in           = D_in, X00_cube_in = X00_cube_in,
    A_list_in      = A_list_in, AB_list_in = AB_list_in,
    hessian = FALSE
  )
)

# =============================================================================
# 8. SAVE OUTPUTS
# =============================================================================

save.image("albumin_fit.RData")

# Back-transform fitted parameters to natural scale for interpretation
# theta_est[1:4q] = exp(theta_trans_est[1:4q]) restores zeta, nu2, xi, d to positive scale
theta_trans_est <- fit_joint_albumin_1$par
theta_est        <- theta_trans_est
theta_est[1:(4 * q)] <- exp(theta_trans_est[1:(4 * q)])

save(theta_est, file = "theta_est_0613.RData")
