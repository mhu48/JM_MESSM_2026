# =============================================================================
# Helper Functions for Computing Loss Metrics (Sensitivity Analysis)
# =============================================================================
# Sourced by model_comparison_synth.R. Provides longitudinal prediction losses
# (MSE, MAE, NLL), survival prediction metrics (Brier score, time-dependent
# AUC-ROC, C-index), and risk-score computation for the three-model comparison.
#
# Functions:
#   compute_longitudinal_mse / mae / nll  -> prediction error metrics
#   compute_brier_score / integrated_brier -> survival calibration
#   compute_auc_roc / compute_cindex       -> survival discrimination
#   extract_predictions                    -> run filter, return prediction arrays
#   compute_risk_scores                    -> risk from a prediction-array result
#   compute_risk_scores_from_filter        -> risk directly from filter output
#                                             (f_filter, b_filter); used by
#                                             model_comparison_synth.R
#   compute_all_losses                     -> wrapper computing all metrics
#
# Parameter vector structure: see albumin_fit.R header.
# =============================================================================

library(survival)
library(pROC)
library(dplyr)

## ============================================================
## LONGITUDINAL LOSS FUNCTIONS
## ============================================================

#' Compute MSE for longitudinal predictions
#' @param pred_array q x T x m array of predictions
#' @param obs_array q x T x m array of observations (may contain NA)
#' @return named vector with overall MSE and per-variable MSE
compute_longitudinal_mse <- function(pred_array, obs_array) {
  q <- dim(pred_array)[1]
  
  # Overall MSE (across all variables and subjects)
  diff_sq <- (pred_array - obs_array)^2
  overall_mse <- mean(diff_sq, na.rm = TRUE)
  
  # Per-variable MSE
  per_var_mse <- numeric(q)
  for(k in 1:q) {
    per_var_mse[k] <- mean((pred_array[k,,] - obs_array[k,,])^2, na.rm = TRUE)
  }
  
  names(per_var_mse) <- c("ALB", "IDWG", "SBP", "NLR", "HGB")[1:q]
  
  return(c(overall = overall_mse, per_var_mse))
}

#' Compute MAE for longitudinal predictions
compute_longitudinal_mae <- function(pred_array, obs_array) {
  q <- dim(pred_array)[1]
  
  # Overall MAE
  diff_abs <- abs(pred_array - obs_array)
  overall_mae <- mean(diff_abs, na.rm = TRUE)
  
  # Per-variable MAE
  per_var_mae <- numeric(q)
  for(k in 1:q) {
    per_var_mae[k] <- mean(abs(pred_array[k,,] - obs_array[k,,]), na.rm = TRUE)
  }
  
  names(per_var_mae) <- c("ALB", "IDWG", "SBP", "NLR", "HGB")[1:q]
  
  return(c(overall = overall_mae, per_var_mae))
}

#' Compute negative log-likelihood for longitudinal predictions
#' Uses Gaussian assumption: -log N(y | mu, sigma^2)
compute_longitudinal_nll <- function(pred_array, obs_array, sigma_est) {
  # sigma_est is the estimated error covariance matrix from model fit
  
  q <- dim(pred_array)[1]
  npt <- dim(pred_array)[2]
  m <- dim(pred_array)[3]
  
  nll <- 0
  count <- 0
  
  for(i in 1:m) {
    for(t in 1:npt) {
      obs_t <- obs_array[,t,i]
      pred_t <- pred_array[,t,i]
      
      if(!any(is.na(obs_t))) {
        # Multivariate Gaussian log-likelihood
        residual <- obs_t - pred_t
        nll <- nll - (- 0.5 * q * log(2*pi) - 
                        0.5 * log(det(sigma_est)) - 
                        0.5 * t(residual) %*% solve(sigma_est) %*% residual)
        count <- count + 1
      }
    }
  }
  
  return(nll / count)  # Average NLL per observation
}

## ============================================================
## SURVIVAL LOSS FUNCTIONS
## ============================================================

#' Compute Brier Score for survival predictions
#' @param surv_probs_list List where each element is time-specific survival probabilities (vector length m)
#' @param event_indicator m-vector of event indicators (1=died, 0=censored)
#' @param event_times m-vector of event/censoring times
#' @param eval_times Vector of times at which to evaluate Brier score
#' @return Vector of Brier scores at each evaluation time
compute_brier_score <- function(surv_probs_list, event_indicator, event_times, eval_times) {
  
  m <- length(event_indicator)
  brier_scores <- numeric(length(eval_times))
  
  for(t_idx in seq_along(eval_times)) {
    t_eval <- eval_times[t_idx]
    
    # Get predicted survival probability at time t_eval
    pred_surv <- surv_probs_list[[t_idx]]
    
    # True survival indicator at time t_eval
    # 1 if survived past t_eval, 0 if event before t_eval
    true_surv <- as.numeric(event_times > t_eval)
    
    # For censored cases before t_eval, we don't know true status
    # Use IPCW (Inverse Probability of Censoring Weighting) approximation
    # For simplicity here, exclude censored observations before t_eval
    
    uncensored_mask <- (event_indicator == 1) | (event_times > t_eval)
    
    if(sum(uncensored_mask) > 0) {
      brier_scores[t_idx] <- mean((pred_surv[uncensored_mask] - true_surv[uncensored_mask])^2)
    } else {
      brier_scores[t_idx] <- NA
    }
  }
  
  return(brier_scores)
}

#' Compute time-dependent AUC-ROC
#' @param risk_scores m-vector of predicted risk scores (higher = more likely to die)
#' @param event_indicator m-vector of event indicators
#' @param event_times m-vector of event/censoring times
#' @param eval_time Single time point at which to evaluate AUC
#' @return AUC value
compute_auc_roc <- function(risk_scores, event_indicator, event_times, eval_time) {
  
  # Define cases and controls at eval_time
  # Cases: died before eval_time
  # Controls: survived past eval_time (either censored after or event after)
  
  cases <- (event_indicator == 1) & (event_times <= eval_time)
  controls <- (event_times > eval_time)
  
  # Exclude those censored before eval_time (unknown status)
  usable <- cases | controls
  
  if(sum(cases) > 0 & sum(controls) > 0) {
    roc_obj <- roc(response = cases[usable], predictor = risk_scores[usable], 
                   direction = ">", quiet = TRUE)
    return(as.numeric(auc(roc_obj)))
  } else {
    return(NA)
  }
}

#' Compute Harrell's C-index (concordance)
#' @param risk_scores m-vector of predicted risk scores
#' @param event_times m-vector of event/censoring times
#' @param event_indicator m-vector of event indicators
#' @return C-index value
compute_cindex <- function(risk_scores, event_times, event_indicator) {
  
  # Use survival package concordance
  surv_obj <- Surv(time = event_times, event = event_indicator)
  
  # Concordance expects lower score = better survival
  # Our risk_scores have higher = worse, so we negate
  cindex_obj <- concordance(surv_obj ~ I(-risk_scores))
  
  return(cindex_obj$concordance)
}

#' Compute integrated Brier score (across all time points)
compute_integrated_brier <- function(surv_probs_list, event_indicator, event_times, eval_times) {
  brier_vec <- compute_brier_score(surv_probs_list, event_indicator, event_times, eval_times)
  return(mean(brier_vec, na.rm = TRUE))
}

## ============================================================
## PREDICTION EXTRACTION FUNCTIONS
## ============================================================

#' Extract filtering estimates (predictions) from fitted model
#' @param theta_est Estimated parameters (natural scale)
#' @param data_split q x T x m array
#' @param A_list List of non-missing subject indices at each time
#' @param model_type One of "cubic_ou", "cubic_cubic", "cubic_local"
#' @return List with f_filter (population), b_filter (individual deviations), y_pred (predictions)
extract_predictions <- function(theta_est, data_split, A_list, model_type = "cubic_ou") {
  
  q <- dim(data_split)[1]
  npt <- dim(data_split)[2]
  m <- dim(data_split)[3]
  
  # Source appropriate filter function based on model type
  if(model_type == "cubic_ou") {
    # Use existing filter_estimate_OU_withP function
    if(!exists("filter_estimate_OU_withP")) {
      stop("filter_estimate_OU_withP function not found. Source it first.")
    }
    
    # Prepare theta for OU model
    theta_for_filter <- theta_est[1:(0.5 * q^2 + 3.5 * q)]
    filter_result <- filter_estimate_OU_withP(theta_for_filter, data_split, A_list)
    
  } else if(model_type == "cubic_cubic") {
    # Need filter function for cubic-cubic
    if(!exists("filter_estimate_cs_cs_withP")) {
      stop("filter_estimate_cs_cs_withP function not found. Create it first.")
    }
    
    theta_for_filter <- theta_est[1:(0.5 * q^2 + 2.5 * q + 1)]
    filter_result <- filter_estimate_cs_cs_withP(theta_for_filter, data_split, A_list)
    
  } else if(model_type == "cubic_local") {
    # Need filter function for cubic-local
    if(!exists("filter_estimate_cs_local_withP")) {
      stop("filter_estimate_cs_local_withP function not found. Create it first.")
    }
    
    theta_for_filter <- theta_est[1:(0.5 * q^2 + 2.5 * q + 1)]
    filter_result <- filter_estimate_cs_local_withP(theta_for_filter, data_split, A_list)
  }
  
  # Compute predictions: y_pred[k,t,i] = f[k,t] + b[k,t,i]
  f_filter <- filter_result$f_filter_mat  # q x T
  b_filter <- filter_result$b_filter_cube # q x T x m
  
  y_pred <- array(NA, dim = c(q, npt, m))
  for(i in 1:m) {
    y_pred[,,i] <- f_filter + b_filter[,,i]
  }
  
  return(list(
    f_filter = f_filter,
    b_filter = b_filter,
    y_pred = y_pred,
    P0_cube = filter_result$P0_cube,
    P1_cube = filter_result$P1_cube,
    P2_cube = filter_result$P2_cube,
    P3_cube = filter_result$P3_cube
  ))
}

#' Compute survival risk scores for each subject
#' @param theta_est Full parameter vector
#' @param data_split q x T x m data array
#' @param X00_cube m x p x T covariate array
#' @param pred_result Output from extract_predictions()
#' @param model_type Model type string
#' @return m-vector of risk scores (higher = higher risk)
compute_risk_scores <- function(theta_est, data_split, X00_cube, pred_result, model_type = "cubic_ou") {
  
  q <- dim(data_split)[1]
  m <- dim(data_split)[3]
  npt <- dim(data_split)[2]
  p <- dim(X00_cube)[2]
  
  # Extract survival parameters
  if(model_type == "cubic_ou") {
    gamma_start <- 0.5 * q^2 + 3.5 * q + 1
    dim_gamma <- q  # Assuming D_in is just identity for simplicity
    gamma_vec <- theta_est[gamma_start:(gamma_start + dim_gamma - 1)]
    gamma_x <- theta_est[(gamma_start + dim_gamma):(gamma_start + dim_gamma + p)]
  } else {
    # Cubic-cubic and cubic-local have different parameter lengths
    gamma_start <- 0.5 * q^2 + 2.5 * q + 2
    dim_gamma <- if(model_type == "cubic_cubic") 2*q else q
    gamma_vec <- theta_est[gamma_start:(gamma_start + dim_gamma - 1)]
    gamma_x <- theta_est[(gamma_start + dim_gamma):(gamma_start + dim_gamma + p)]
  }
  
  # Compute risk score for each subject based on their trajectory
  # Use time-averaged predictions
  risk_scores <- numeric(m)
  
  for(i in 1:m) {
    # Get predicted biomarkers for subject i (average over time)
    y_pred_i <- pred_result$y_pred[,,i]
    y_mean_i <- rowMeans(y_pred_i, na.rm = TRUE)
    
    # Get covariates (use baseline)
    X_i <- X00_cube[i,,1]
    
    # Compute logit(risk) = gamma_vec^T y + gamma_x^T [1, X]
    logit_risk <- sum(gamma_vec * y_mean_i) + gamma_x[1] + sum(gamma_x[-1] * X_i)
    
    risk_scores[i] <- plogis(logit_risk)  # Convert to probability
  }
  
  return(risk_scores)
}

## ============================================================
## WRAPPER FUNCTION: COMPUTE ALL LOSSES
## ============================================================

#' Compute all loss metrics for a given model
#' @param theta_est Estimated parameters
#' @param data_test q x T x m test data
#' @param z_test m x T event indicators
#' @param X00_cube_test m x p x T test covariates
#' @param A_list_test List of non-missing indices
#' @param model_type Model type string
#' @return List with all loss metrics
compute_all_losses <- function(theta_est, data_test, z_test, X00_cube_test, 
                                A_list_test, model_type = "cubic_ou") {
  
  q <- dim(data_test)[1]
  npt <- dim(data_test)[2]
  m <- dim(data_test)[3]
  
  ## 1. Extract predictions
  pred_result <- extract_predictions(theta_est, data_test, A_list_test, model_type)
  
  ## 2. Longitudinal losses
  mse_results <- compute_longitudinal_mse(pred_result$y_pred, data_test)
  mae_results <- compute_longitudinal_mae(pred_result$y_pred, data_test)
  
  # For NLL, extract Sigma from theta_est
  if(model_type == "cubic_ou") {
    d_vec <- theta_est[(3*q+1):(4*q)]
    if(q > 1) {
      l_vec <- theta_est[(4*q+1):(0.5*q^2 + 3.5*q)]
      L <- diag(q)
      L[lower.tri(L)] <- l_vec
    } else {
      L <- matrix(1)
    }
  } else {
    d_vec <- theta_est[(2*q+1):(3*q)]
    if(q > 1) {
      l_vec <- theta_est[(3*q+1):(0.5*q^2 + 2.5*q)]
      L <- diag(q)
      L[lower.tri(L)] <- l_vec
    } else {
      L <- matrix(1)
    }
  }
  
  D_mat <- diag(d_vec)
  Sigma_est <- L %*% D_mat %*% t(L)
  
  nll_result <- compute_longitudinal_nll(pred_result$y_pred, data_test, Sigma_est)
  
  ## 3. Survival losses
  # Extract event times and indicators
  event_times <- apply(z_test, 1, function(row) {
    death_time <- which(row == 1)
    if(length(death_time) > 0) return(death_time[1])
    else return(npt)  # Censored at end
  })
  
  event_indicator <- apply(z_test, 1, function(row) {
    as.numeric(any(row == 1, na.rm = TRUE))
  })
  
  # Compute risk scores
  risk_scores <- compute_risk_scores(theta_est, data_test, X00_cube_test, 
                                      pred_result, model_type)
  
  # Evaluate at multiple time horizons
  eval_times <- c(6, 12, 24, 36)  # 6-month, 1-year, 2-year, 3-year
  eval_times <- eval_times[eval_times <= npt]
  
  # For Brier score, need time-specific survival probabilities
  # Simplified: use risk_scores as proxy for all times
  surv_probs_list <- lapply(eval_times, function(t) 1 - risk_scores)
  
  brier_results <- compute_brier_score(surv_probs_list, event_indicator, 
                                        event_times, eval_times)
  names(brier_results) <- paste0("t", eval_times)
  
  # AUC at different times
  auc_results <- sapply(eval_times, function(t) {
    compute_auc_roc(risk_scores, event_indicator, event_times, t)
  })
  names(auc_results) <- paste0("t", eval_times)
  
  # C-index (overall)
  cindex_result <- compute_cindex(risk_scores, event_times, event_indicator)
  
  # Integrated Brier
  int_brier <- compute_integrated_brier(surv_probs_list, event_indicator, 
                                         event_times, eval_times)
  
  ## 4. Return all results
  return(list(
    # Longitudinal
    mse = mse_results,
    mae = mae_results,
    nll = nll_result,
    
    # Survival
    brier = brier_results,
    integrated_brier = int_brier,
    auc = auc_results,
    cindex = cindex_result,
    
    # For reference
    model_type = model_type
  ))
}


## ============================================================
## RISK SCORES DIRECTLY FROM FILTER OUTPUT
## ============================================================

#' Compute individual mortality risk scores directly from Kalman filter output.
#' Used by model_comparison_synth.R. Equivalent to compute_risk_scores() but
#' takes the filtered population trajectory and individual deviations directly
#' rather than a prediction-array result object.
#'
#' @param theta_est_nat  fitted parameters on natural scale
#' @param theta_est_trans fitted parameters on transformed (log) scale (unused; kept for signature compatibility)
#' @param f_filter       q x npt population trajectory (filtered)
#' @param b_filter       q x npt x m individual deviations (filtered);
#'                       note: permute with aperm() if filter returns q x m x npt
#' @param X_cube         m x p x npt covariate array
#' @param D_mat          design matrix mapping state to survival link
#' @param model_type     one of "cubic_ou", "cubic_cubic", "cubic_local"
#' @return numeric vector of length m, predicted one-step mortality probabilities
compute_risk_scores_from_filter <- function(theta_est_nat, theta_est_trans,
                                            f_filter, b_filter,
                                            X_cube, D_mat, model_type = "cubic_ou") {
  m       <- dim(b_filter)[3]
  q_loc   <- dim(f_filter)[1]
  npt_loc <- dim(f_filter)[2]
  p_loc   <- dim(X_cube)[2]

  if(model_type == "cubic_ou") {
    gamma_start <- as.integer(0.5 * q_loc^2 + 3.5 * q_loc) + 1
    dim_gamma   <- nrow(D_mat)
    gamma_vec   <- theta_est_nat[gamma_start:(gamma_start + dim_gamma - 1)]
    gamma_x     <- theta_est_nat[(gamma_start + dim_gamma):(gamma_start + dim_gamma + p_loc)]
  }

  risk_scores <- numeric(m)
  for(i in 1:m) {
    y_mean      <- rowMeans(f_filter + b_filter[, , i], na.rm = TRUE)
    X_base      <- X_cube[i, , 1]
    logit_risk  <- sum(gamma_vec * y_mean) + gamma_x[1] + sum(gamma_x[-1] * X_base)
    risk_scores[i] <- plogis(logit_risk)
  }
  return(risk_scores)
}
