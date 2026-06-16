# =============================================================================
# Data Preparation Helper Functions
# =============================================================================
# Helper functions for constructing model inputs from long-format data.
# Sourced by model_comparison_synth.R (synthetic data pipeline).
# Not needed for the real data pipeline (albumin_fit.R uses pre-built arrays).
#
# Functions:
#   prepare_data_array(data_long)   -> q x npt x m array of biomarker observations
#   prepare_z_matrix(data_long)     -> m x npt survival indicator matrix
#   prepare_X_cube(data_long)       -> m x p x npt covariate array
#   standardize_data_array(arr)     -> list(data_standardized, means, sds)
#   standardize_X_cube(X_cube)      -> list(X_standardized, means, sds, standardized_cols)
#   create_A_list(z_mat)            -> list of length npt; alive-at-t subject indices
#   create_AB_list(z_mat)           -> list of length npt; present-at-t subject indices
#                                      (CRITICAL: length = npt, not npt-1)
# =============================================================================

library(dplyr)

#' Convert long-format data to q x T x m array
#' @param data_long Long-format data frame with columns: dummyPATIENT_ID, time_on_trt, biomarker columns
#' @return q x T x m array
prepare_data_array <- function(data_long) {
  
  subject_ids <- sort(unique(data_long$dummyPATIENT_ID))
  m <- length(subject_ids)
  time_points <- sort(unique(data_long$time_on_trt))
  npt <- length(time_points)
  
  biomarker_cols <- c("ALBUMIN_avg", "idwg_percent_avg", "PRE_SBP_avg", 
                      "nlr_avg_sqrt", "hgb_avg")
  q <- length(biomarker_cols)
  
  data_array <- array(NA_real_, dim = c(q, npt, m))
  
  for(i in seq_along(subject_ids)) {
    subj_data <- data_long %>% 
      filter(dummyPATIENT_ID == subject_ids[i]) %>%
      arrange(time_on_trt)
    
    for(k in 1:q) {
      data_array[k, , i] <- subj_data[[biomarker_cols[k]]]
    }
  }
  
  return(data_array)
}

#' Standardize data array (matching albumin_fit_0613.R lines 66-70)
#' Each biomarker is standardized to mean 0, sd 1 across ALL observations
#' @param data_array q x T x m array
#' @return List with standardized array and scaling parameters
standardize_data_array <- function(data_array) {
  
  q <- dim(data_array)[1]
  
  means <- numeric(q)
  sds <- numeric(q)
  data_std <- data_array
  
  for(k in 1:q) {
    means[k] <- mean(data_array[k,,], na.rm = TRUE)
    sds[k] <- sd(as.vector(data_array[k,,]), na.rm = TRUE)
    data_std[k,,] <- (data_array[k,,] - means[k]) / sds[k]
  }
  
  return(list(
    data_standardized = data_std,
    means = means,
    sds = sds
  ))
}

#' Convert event data to z_in matrix matching original format
#' @param data_long Long-format data frame
#' @return m x npt matrix where:
#'   0 = alive at time t
#'   1 = died at time t (only appears once per subject)
#'   NA = already dead (times after death)
prepare_z_matrix <- function(data_long) {
  
  subject_ids <- sort(unique(data_long$dummyPATIENT_ID))
  m <- length(subject_ids)
  time_points <- sort(unique(data_long$time_on_trt))
  npt <- length(time_points)
  
  z_mat <- matrix(NA_real_, nrow = m, ncol = npt)
  
  for(i in seq_along(subject_ids)) {
    subj_data <- data_long %>%
      filter(dummyPATIENT_ID == subject_ids[i]) %>%
      arrange(time_on_trt)
    
    z_mat[i, ] <- subj_data$event_died
  }
  
  # Convert to proper format: 0 before death, 1 at death, NA after death
  for(i in 1:m) {
    death_idx <- which(z_mat[i, ] == 1)
    if(length(death_idx) > 0) {
      death_time <- min(death_idx)
      if(death_time > 1) {
        z_mat[i, 1:(death_time - 1)] <- 0
      }
      z_mat[i, death_time] <- 1
      if(death_time < npt) {
        z_mat[i, (death_time + 1):npt] <- NA
      }
    } else {
      # Censored - all 0s up to last observed time
      last_obs <- max(which(!is.na(z_mat[i, ])))
      z_mat[i, 1:last_obs] <- 0
      # If censored before end, set remaining to NA
      # For synthetic data all subjects have full follow-up or die,
      # so this handles both cases.
    }
  }
  
  return(z_mat)
}

#' Convert covariates to m x p x npt array
#' @param data_long Long-format data frame
#' @return m x p x npt array
prepare_X_cube <- function(data_long) {
  
  subject_ids <- sort(unique(data_long$dummyPATIENT_ID))
  m <- length(subject_ids)
  time_points <- sort(unique(data_long$time_on_trt))
  npt <- length(time_points)
  
  covar_cols <- c("age", "male", "EKTV_avg", "diabetic", "ACCESS_AVF_avg", "ACCESS_AVG_avg")
  p <- length(covar_cols)
  
  X_cube <- array(NA_real_, dim = c(m, p, npt))
  
  for(i in seq_along(subject_ids)) {
    subj_data <- data_long %>%
      filter(dummyPATIENT_ID == subject_ids[i]) %>%
      arrange(time_on_trt)
    
    for(t in 1:npt) {
      for(k in 1:p) {
        X_cube[i, k, t] <- subj_data[[covar_cols[k]]][t]
      }
    }
  }
  
  return(X_cube)
}

#' Standardize covariates (matching albumin_fit_0613.R lines 72-76)
#' Standardizes columns 1 (age) and 3 (EKTV_avg) only
#' @param X_cube m x p x npt array
#' @return List with standardized array and scaling parameters
standardize_X_cube <- function(X_cube) {
  
  X_std <- X_cube
  p <- dim(X_cube)[2]
  
  means_x <- numeric(p)
  sds_x <- numeric(p)
  standardize_cols <- c(1, 3)  # age and EKTV_avg
  
  for(k in 1:p) {
    means_x[k] <- mean(X_cube[, k, ], na.rm = TRUE)
    sds_x[k] <- sd(as.vector(X_cube[, k, ]), na.rm = TRUE)
    
    if(k %in% standardize_cols) {
      X_std[, k, ] <- (X_cube[, k, ] - means_x[k]) / sds_x[k]
    }
  }
  
  return(list(
    X_standardized = X_std,
    means = means_x,
    sds = sds_x,
    standardized_cols = standardize_cols
  ))
}

#' Create A_list: subjects alive (z=0) at each time
#' Matches albumin_fit_0613.R lines 85-91:
#'   dummy_A = z_in; dummy_A[dummy_A == 1] = NA
#'   A_list[[j]] = (1:m)[!is.na(dummy_A[, j])]
#' @param z_mat m x npt matrix
#' @return List of length npt
create_A_list <- function(z_mat) {
  
  npt <- ncol(z_mat)
  m <- nrow(z_mat)
  
  # Replicate original logic exactly:
  # dummy_A = z_in; dummy_A[dummy_A == 1] = NA
  # Then subjects with non-NA = subjects with z=0 (alive, not dying at this time)
  dummy_A <- z_mat
  dummy_A[dummy_A == 1] <- NA
  
  A_list <- vector("list", npt)
  
  for(t in 1:npt) {
    A_list[[t]] <- which(!is.na(dummy_A[, t]))
  }
  
  return(A_list)
}

#' Create AB_list: subjects alive or just died (z != NA) at each time
#' Matches albumin_fit_0613.R lines 92-96:
#'   dummy_AB = z_in
#'   AB_list[[j]] = (1:m)[!is.na(dummy_AB[, j])]
#' 
#' CRITICAL: Must have length npt (NOT npt-1) because C++ code
#' accesses AB_list[j+1] for j = 0..npt-2, i.e. elements 1..npt-1 (0-indexed).
#' 
#' @param z_mat m x npt matrix
#' @return List of length npt (NOT npt-1)
create_AB_list <- function(z_mat) {
  
  npt <- ncol(z_mat)
  m <- nrow(z_mat)
  
  # Replicate original logic exactly:
  # dummy_AB = z_in (no modification)
  # AB_list[[j]] = (1:m)[!is.na(dummy_AB[, j])]
  # This includes subjects with z=0 (alive) OR z=1 (dying at this time)
  # Excludes only NA (already dead)
  
  AB_list <- vector("list", npt)  # was npt-1, must be npt
  
  for(t in 1:npt) {
    AB_list[[t]] <- which(!is.na(z_mat[, t]))
  }
  
  return(AB_list)
}

#' Rescale predictions back to original scale
#' @param pred_array q x T x m array of standardized predictions
#' @param means Vector of length q
#' @param sds Vector of length q
#' @return q x T x m array in original scale
rescale_predictions <- function(pred_array, means, sds) {
  
  q <- dim(pred_array)[1]
  pred_rescaled <- pred_array
  
  for(k in 1:q) {
    pred_rescaled[k,,] <- pred_array[k,,] * sds[k] + means[k]
  }
  
  return(pred_rescaled)
}
