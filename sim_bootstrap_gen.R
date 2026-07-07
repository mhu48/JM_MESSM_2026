rm(list = ls())

set.seed(123)

library(GA)
library(data.table)
library(dplyr)
library(assist)
library(parallel)
library(foreach)
library(doParallel)

#----------------------------------------------------------------------------------------------------------------------------------------
# Generate bootstrap estimate matrices
#----------------------------------------------------------------------------------------------------------------------------------------

project_dir <- normalizePath(".", mustWork = TRUE)

sim_results_dir <- file.path(project_dir, "data", "results_m1000")
bootstrap_dir   <- file.path(project_dir, "bootstrap_results", "results_m1000")

if(!dir.exists(bootstrap_dir)){
  dir.create(bootstrap_dir, recursive = TRUE)
}

source_file <- file.path(project_dir, "R", "llh_survival_x_OU.R")
if(!file.exists(source_file)){
  source_file <- file.path(project_dir, "llh_survival_x_OU.R")
}
source(source_file)

load(file.path(sim_results_dir, "y_list.Rdata"))
load(file.path(sim_results_dir, "z_list.Rdata"))
load(file.path(sim_results_dir, "X00_cube_in_list.Rdata"))
load(file.path(sim_results_dir, "A_list.Rdata"))
load(file.path(sim_results_dir, "AB_list.Rdata"))
load(file.path(sim_results_dir, "theta_est_mat.Rdata"))
load(file.path(sim_results_dir, "MSE_fix_mat.Rdata"))

q <- 5
p <- 6
m <- 1000

n_boot_files   <- 100
n_boot_samples <- 100

#----------------------------------------------------------------------------------------------------------------------------------------
# Select the simulation dataset used for bootstrap resampling
#----------------------------------------------------------------------------------------------------------------------------------------

npt <- dim(y_list[[1]])[2]

n_time_points <- numeric(length(y_list))

for(sim_id in seq_along(y_list)){
  k <- 1
  while(k + 1 <= npt && length(A_list[[sim_id]][[k + 1]]) > 0){
    k <- k + 1
  }
  n_time_points[sim_id] <- k
}

sim_idx <- which(n_time_points == 50)

MSE_fix_mat_full <- MSE_fix_mat[, sim_idx]
MSE_fix_mat_sum <- colSums(MSE_fix_mat)
MSE_fix_mat_full_sum <- colSums(MSE_fix_mat_full)

j_min <- seq_len(100)[MSE_fix_mat_sum == min(MSE_fix_mat_full_sum)]
j_min <- j_min[1]

data_split_in <- y_list[[j_min]]
z_in <- z_list[[j_min]]
X00_cube_in <- X00_cube_in_list[[j_min]]

idx <- seq_len(m)

data_split_in_selected <- data_split_in[, , idx]
z_in_selected <- z_in[idx, ]
X00_cube_in_selected <- X00_cube_in[idx, , ]

m <- dim(data_split_in_selected)[3]
npt <- dim(data_split_in_selected)[2]

D_in <- cbind(
  diag(q) %x% matrix(c(1, 0), nrow = 1, ncol = 2),
  diag(q)
)

theta_length <- as.integer(0.5 * q^2 + 4.5 * q + p + 1)

theta_trans <- c(
  log(theta_est_mat[, j_min][1:(4 * q)]),
  theta_est_mat[, j_min][(4 * q + 1):theta_length]
)

time_elapse <- c()

#----------------------------------------------------------------------------------------------------------------------------------------
# Bootstrap model fitting
#----------------------------------------------------------------------------------------------------------------------------------------

for(r in seq_len(n_boot_files)){
  
  theta_est_mat1 <- matrix(
    0,
    nrow = theta_length,
    ncol = n_boot_samples
  )
  
  for(sample_id in seq_len(n_boot_samples)){
    
    ids <- sample(seq_len(m), size = m, replace = TRUE)
    
    data_split_boot <- data_split_in_selected[, , ids]
    z_boot <- z_in_selected[ids, ]
    X00_cube_boot <- X00_cube_in_selected[ids, , ]
    
    dummy_A <- z_boot
    dummy_A[dummy_A == 1] <- NA
    
    A_list_boot <- vector("list", npt)
    
    for(j in seq_len(npt)){
      A_list_boot[[j]] <- seq_len(m)[!is.na(dummy_A[, j])]
    }
    
    dummy_AB <- z_boot
    
    AB_list_boot <- vector("list", npt)
    
    for(j in seq_len(npt)){
      AB_list_boot[[j]] <- seq_len(m)[!is.na(dummy_AB[, j])]
    }
    
    llh_survival_x_OU(
      theta_trans,
      data_split_in = data_split_boot,
      z_in = z_boot,
      D_in = D_in,
      X00_cube_in = X00_cube_boot,
      A_list_in = A_list_boot,
      AB_list_in = AB_list_boot
    )
    
    time_here <- system.time(
      fit_joint_albumin <- optim(
        theta_trans,
        llh_survival_x_OU,
        data_split_in = data_split_boot,
        z_in = z_boot,
        D_in = D_in,
        X00_cube_in = X00_cube_boot,
        A_list_in = A_list_boot,
        AB_list_in = AB_list_boot,
        hessian = FALSE
      )
    )
    
    theta_trans_est <- fit_joint_albumin$par
    
    theta_est <- theta_trans_est
    theta_est[1:(4 * q)] <- exp(theta_trans_est[1:(4 * q)])
    
    theta_est_mat1[, sample_id] <- theta_est
    
    time_elapse <- c(time_elapse, time_here[1])
  }
  
  save(
    theta_est_mat1,
    file = file.path(bootstrap_dir, paste0("theta_est_mat", r, ".Rdata"))
  )
}

save(
  time_elapse,
  file = file.path(project_dir, "bootstrap_results", "time_elapse.Rdata")
)