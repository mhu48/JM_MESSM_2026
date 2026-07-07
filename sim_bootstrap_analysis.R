rm(list = ls())

set.seed(123)

library(xtable)

#----------------------------------------------------------------------------------------------------------------------------------------
# Paths and settings
#----------------------------------------------------------------------------------------------------------------------------------------

project_dir <- normalizePath(".", mustWork = TRUE)

sim_results_dir <- file.path(project_dir, "data", "results_m1000")
output_dir <- file.path(project_dir, "analysis_results")

if(!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}

q <- 5
p <- 6
m <- 1000

#----------------------------------------------------------------------------------------------------------------------------------------
# MSE, bias, and variance based on transformed Sigma_epsilon parameters
#----------------------------------------------------------------------------------------------------------------------------------------

load(file.path(sim_results_dir, "theta_est_mat.Rdata"))
load(file.path(sim_results_dir, "true_val.Rdata"))

MSE_fn <- function(estimates, true_val){
  mean((estimates - true_val)^2)
}

bias_fn <- function(estimates, true_val){
  mean(estimates) - true_val
}

var_fn <- function(estimates){
  n <- length(estimates)
  var(estimates) * (n - 1) / n
}

estimates <- matrix(
  0,
  nrow = dim(theta_est_mat)[1],
  ncol = dim(theta_est_mat)[2]
)

for(i in seq_len(100)){
  
  d_vec <- theta_est_mat[(3 * q + 1):(4 * q), i]
  l_vec <- theta_est_mat[(4 * q + 1):(0.5 * q^2 + 3.5 * q), i]
  
  L <- diag(q)
  L[lower.tri(L)] <- l_vec
  
  D <- diag(d_vec, nrow = q, ncol = q)
  
  Sig_est <- L %*% D %*% t(L)
  
  index <- 3 * q + 1
  
  for(j in seq_len(q)){
    for(k in seq_len(j)){
      estimates[index, i] <- Sig_est[j, k]
      index <- index + 1
    }
  }
}

estimates[1:(3 * q), ] <- theta_est_mat[1:(3 * q), ]
estimates[((0.5 * q^2 + 3.5 * q) + 1):42, ] <-
  theta_est_mat[((0.5 * q^2 + 3.5 * q) + 1):42, ]

d_vec <- true_val[(3 * q + 1):(4 * q)]
l_vec <- true_val[(4 * q + 1):(0.5 * q^2 + 3.5 * q)]

L <- diag(q)
L[lower.tri(L)] <- l_vec

D <- diag(d_vec, nrow = q, ncol = q)

Sig_est <- L %*% D %*% t(L)

true_val_trans <- rep(0, length(true_val))

index <- 16

for(j in seq_len(q)){
  for(k in seq_len(j)){
    true_val_trans[index] <- Sig_est[j, k]
    index <- index + 1
  }
}

true_val_trans[1:(3 * q)] <- true_val[1:(3 * q)]
true_val_trans[((0.5 * q^2 + 3.5 * q) + 1):42] <-
  true_val[((0.5 * q^2 + 3.5 * q) + 1):42]

MSE <- c()
bias <- c()
variance <- c()

for(i in seq_along(true_val)){
  MSE[i] <- MSE_fn(estimates[i, ], true_val_trans[i])
  bias[i] <- bias_fn(estimates[i, ], true_val_trans[i])
  variance[i] <- var_fn(estimates[i, ])
}

par_names <- c(
  "zeta_1^2", "zeta_2^2", "zeta_3^2", "zeta_4^2", "zeta_5^2",
  "nu_1^2", "nu_2^2", "nu_3^2", "nu_4^2", "nu_5^2",
  "xi_1", "xi_2", "xi_3", "xi_4", "xi_5",
  "l_11", "l_21", "l_22", "l_31", "l_32", "l_33",
  "l_41", "l_42", "l_43", "l_44",
  "l_51", "l_52", "l_53", "l_54", "l_55",
  "gamma_{21}", "gamma_{22}", "gamma_{23}", "gamma_{24}", "gamma_{25}",
  "gamma_{10}", "gamma_{11}", "gamma_{12}", "gamma_{13}",
  "gamma_{14}", "gamma_{15}", "gamma_{16}"
)

par_est_df <- data.frame(
  Parameter = par_names,
  "True value" = true_val_trans,
  MSE = MSE,
  Variance = variance,
  Bias = bias
)

save(
  par_est_df,
  file = file.path(output_dir, "par_est_df.Rdata")
)

par_est_df_table <- data.frame(
  Parameter = par_est_df$Parameter,
  True_value = par_est_df$True.value,
  MSE = par_est_df$MSE * 100,
  Variance = par_est_df$Variance,
  Bias = par_est_df$Bias
)

print(
  xtable(par_est_df_table, digits = c(0, 2, 2, 2, 2, 2)),
  include.rownames = FALSE
)

#----------------------------------------------------------------------------------------------------------------------------------------
# Bootstrap coverage based on transformed Sigma_epsilon parameters
#----------------------------------------------------------------------------------------------------------------------------------------

load(file.path(project_dir, "theta_mat_est_extra_final.Rdata"))

for(i in seq_along(theta_mat_est_extra_final)){
  
  n_boot_samples <- ncol(theta_mat_est_extra_final[[i]])
  
  for(j in seq_len(n_boot_samples)){
    
    d_vec <- theta_mat_est_extra_final[[i]][(3 * q + 1):(4 * q), j]
    l_vec <- theta_mat_est_extra_final[[i]][(4 * q + 1):(0.5 * q^2 + 3.5 * q), j]
    
    L <- diag(q)
    L[lower.tri(L)] <- l_vec
    
    D <- diag(d_vec)
    
    Sig_est <- L %*% D %*% t(L)
    
    index <- 16
    
    for(k in seq_len(q)){
      for(l in seq_len(k)){
        theta_mat_est_extra_final[[i]][index, j] <- Sig_est[k, l]
        index <- index + 1
      }
    }
  }
}

CI_record <- array(
  0,
  dim = c(2, length(theta_mat_est_extra_final), 42)
)

for(i in seq_along(theta_mat_est_extra_final)){
  CI_record[1, i, ] <- apply(theta_mat_est_extra_final[[i]], 1, quantile, probs = 0.025)
  CI_record[2, i, ] <- apply(theta_mat_est_extra_final[[i]], 1, quantile, probs = 0.975)
}

cover_rate_100 <- numeric(42)

for(i in seq_len(42)){
  cover_rate_100[i] <- mean(
    par_est_df$True.value[i] >= CI_record[1, , i] &
      par_est_df$True.value[i] <= CI_record[2, , i]
  )
}

cover_rate_100
order(cover_rate_100)
sort(cover_rate_100)

par_est_boot_df <- data.frame(
  Parameter = par_names,
  "True value" = true_val_trans,
  MSE = MSE,
  Variance = variance,
  Bias = bias,
  Boot = round(cover_rate_100 * 100)
)

save(
  par_est_boot_df,
  file = file.path(output_dir, "par_est_boot_df.Rdata")
)

par_est_boot_df_table <- data.frame(
  Parameter = par_est_boot_df$Parameter,
  True_value = par_est_boot_df$True.value,
  MSE = 100 * par_est_boot_df$MSE,
  Variance = 100 * par_est_boot_df$Variance,
  Bias = 100 * par_est_boot_df$Bias,
  Boot = par_est_boot_df$Boot
)

print(
  xtable(par_est_boot_df_table, digits = c(0, 2, 2, 2, 2, 2, 0)),
  include.rownames = FALSE
)
