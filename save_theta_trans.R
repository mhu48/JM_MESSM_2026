# =============================================================================
# Save Initial Parameter Vector: theta_trans_3_0504.RData
# =============================================================================
# Extracts and saves the starting parameter vector used by albumin_fit.R.
# This script is a one-time utility — run it only if theta_trans_3_0504.RData
# needs to be regenerated from the intermediate fit stored in albumin_data_0504.RData.
#
# Run order: 0 of 4 (prerequisite, only if theta_trans_3_0504.RData is missing)
#   0. save_theta_trans.R      <- YOU ARE HERE (run only if needed)
#   1. albumin_fit.R
#   2. bootstrap_ci.R
#   3. albumin_plot.R
#
# Inputs:
#   albumin_data_0504.RData    # contains fit_joint_albumin_1 (optim object from
#                              # a preliminary fit) and theta_trans_3 (starting vector)
#
# Outputs:
#   theta_trans_3_0504.RData   # starting parameter vector for albumin_fit.R
#   theta_trans_est_0504.RData # fitted parameter vector from preliminary fit
#
# Parameter vector structure: see albumin_fit.R header for full description.
# theta_trans_3 is on the log-transformed scale (same as used in optim).
# =============================================================================

rm(list = ls())

# --- USER: set this to the directory containing all data and script files ---
ROOT_DIR <- "."
setwd(ROOT_DIR)

load("albumin_data_0504.RData")
# Loads: fit_joint_albumin_1 (optim object), theta_trans_3, data_split_in, z_in, X00_cube_in

# Set dimensions
q   <- dim(data_split_in)[1]   # 5 biomarkers
m   <- dim(data_split_in)[3]   # number of subjects
npt <- dim(data_split_in)[2]   # 48 months

# Extract fitted parameter vector from preliminary model fit
# theta_trans_est: log-transformed scale (output of optim)
# theta_trans_3: starting vector used in the preliminary fit (also log-transformed)
theta_trans_est <- fit_joint_albumin_1$par

# Save both for use in albumin_fit.R
save(theta_trans_3,   file = "theta_trans_3_0504.RData")
save(theta_trans_est, file = "theta_trans_est_0504.RData")
