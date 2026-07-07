rm(list = ls())

library(mvtnorm)
library(clusterGeneration)
library(magic)
library(assist)
library(xtable)
library(ggplot2)

#----------------------------------------------------------------------------------------------------------------------------------------
# Compare filtering time for the new algorithm and the univariate treatment
#----------------------------------------------------------------------------------------------------------------------------------------

project_dir <- normalizePath(".", mustWork = TRUE)

source_file_1 <- file.path(project_dir, "R", "llh_OU.R")
source_file_2 <- file.path(project_dir, "R", "llh_OU_unitrt.R")

if(!file.exists(source_file_1)){
  source_file_1 <- file.path(project_dir, "llh_OU.R")
}
if(!file.exists(source_file_2)){
  source_file_2 <- file.path(project_dir, "llh_OU_unitrt.R")
}

source(source_file_1)
source(source_file_2)

output_dir <- file.path(project_dir, "filter_time_results")

if(!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}

set.seed(100)

npt <- 50
q <- 5

d_u <- 2 * q
d_v <- q

tvec <- seq(0, 1, length = npt)

Sig <- genPositiveDefMat("unifcorrmat", dim = q)$Sigma / 100

ch <- chol(Sig)
dd <- diag(ch)

L <- t(ch / dd)

l_vec <- L[lower.tri(L)]
d_vec <- dd^2

D <- diag(x = d_vec, nrow = q, ncol = q)

zeta_vec <- runif(q) * 10
nu2_vec <- rchisq(q, df = 1)
xi_vec <- rchisq(q, df = 1)

dt <- 1

T_u <- diag(q) %x% matrix(c(1, 0, dt, 1), 2, 2)
T_v <- diag(x = exp(-xi_vec * dt), nrow = q, ncol = q)

Z_u <- diag(q) %x% matrix(c(1, 0), nrow = 1, ncol = 2)
Z_v <- diag(q)

Lambda <- matrix(c(dt^3 / 3, dt^2 / 2, dt^2 / 2, dt), 2, 2)

Sig_u <- diag(x = zeta_vec, nrow = q, ncol = q) %x% Lambda

Sig_v <- diag(
  x = 0.5 * nu2_vec * (1 - exp(-2 * xi_vec * dt)) / xi_vec,
  nrow = q,
  ncol = q
)

filter_time <- NULL
filter_time_unitrt <- NULL

#----------------------------------------------------------------------------------------------------------------------------------------
# Direct computation for smaller sample sizes
#----------------------------------------------------------------------------------------------------------------------------------------

mgrid <- c(10, 20, 50, 100, 200)

for(m in mgrid){
  
  delta <- rnorm(d_u)
  
  Q00 <- diag(
    x = 0.5 * nu2_vec / xi_vec,
    nrow = q,
    ncol = q
  )
  
  alpha_u <- matrix(0, nrow = d_u, ncol = npt + 1)
  alpha_u[, 1] <- delta
  
  alpha_v <- array(0, dim = c(d_v, m, npt + 1))
  
  data_split_in <- array(0, dim = c(q, npt, m))
  data_in <- matrix(0, nrow = q * m, ncol = npt)
  
  nsim <- 1
  
  time_filter_vec <- c()
  time_filter_vec_unitrt <- c()
  
  for(sim_id in seq_len(nsim)){
    
    alpha_v[, , 1] <- t(rmvnorm(n = m, sigma = Q00))
    
    for(j in seq_len(npt)){
      
      data_split_in[, j, ] <-
        Z_u %*% (matrix(1, nrow = 1, ncol = m) %x% alpha_u[, j]) +
        Z_v %*% alpha_v[, , j] +
        t(rmvnorm(n = m, sigma = Sig))
      
      data_in[, j] <- as.vector(data_split_in[, j, ])
      
      alpha_u[, j + 1] <- T_u %*% alpha_u[, j] + t(rmvnorm(n = 1, sigma = Sig_u))
      alpha_v[, , j + 1] <- T_v %*% alpha_v[, , j] + t(rmvnorm(n = m, sigma = Sig_v))
    }
    
    dummy <- matrix(0, nrow = m, ncol = npt)
    
    for(i in seq_len(m)){
      for(j in seq_len(npt)){
        if(is.na(data_split_in[1, j, i])){
          dummy[i, j] <- NA
        }
      }
    }
    
    nonNAid_list_in <- vector("list", npt)
    
    for(j in seq_len(npt)){
      nonNAid_list_in[[j]] <- seq_len(m)[!is.na(dummy[, j])]
    }
    
    if(q > 1){
      theta_trans <- c(log(zeta_vec), log(nu2_vec), log(xi_vec), log(d_vec), l_vec)
    } else {
      theta_trans <- c(log(zeta_vec), log(nu2_vec), log(xi_vec), log(d_vec))
    }
    
    time_filter_here <- system.time(
      llh_OU(theta_trans, data_split_in, nonNAid_list_in)
    )
    
    time_filter_vec[sim_id] <- time_filter_here[3]
    
    time_filter_here_unitrt <- system.time(
      llh_OU_unitrt(theta_trans, data_in, q)
    )
    
    time_filter_vec_unitrt[sim_id] <- time_filter_here_unitrt[3]
  }
  
  filter_time <- c(filter_time, mean(time_filter_vec))
  filter_time_unitrt <- c(filter_time_unitrt, mean(time_filter_vec_unitrt))
}

#----------------------------------------------------------------------------------------------------------------------------------------
# Direct computation for larger sample sizes using the new algorithm
#----------------------------------------------------------------------------------------------------------------------------------------

mgrid_large <- c(500, 1000, 10000, 100000, 1000000)

for(m in mgrid_large){
  
  delta <- rnorm(d_u)
  
  Q00 <- diag(
    x = 0.5 * nu2_vec / xi_vec,
    nrow = q,
    ncol = q
  )
  
  alpha_u <- matrix(0, nrow = d_u, ncol = npt + 1)
  alpha_u[, 1] <- delta
  
  alpha_v <- array(0, dim = c(d_v, m, npt + 1))
  
  data_split_in <- array(0, dim = c(q, npt, m))
  data_in <- matrix(0, nrow = q * m, ncol = npt)
  
  nsim <- 1
  
  time_filter_vec <- c()
  
  for(sim_id in seq_len(nsim)){
    
    alpha_v[, , 1] <- t(rmvnorm(n = m, sigma = Q00))
    
    for(j in seq_len(npt)){
      
      data_split_in[, j, ] <-
        Z_u %*% (matrix(1, nrow = 1, ncol = m) %x% alpha_u[, j]) +
        Z_v %*% alpha_v[, , j] +
        t(rmvnorm(n = m, sigma = Sig))
      
      data_in[, j] <- as.vector(data_split_in[, j, ])
      
      alpha_u[, j + 1] <- T_u %*% alpha_u[, j] + t(rmvnorm(n = 1, sigma = Sig_u))
      alpha_v[, , j + 1] <- T_v %*% alpha_v[, , j] + t(rmvnorm(n = m, sigma = Sig_v))
    }
    
    dummy <- matrix(0, nrow = m, ncol = npt)
    
    for(i in seq_len(m)){
      for(j in seq_len(npt)){
        if(is.na(data_split_in[1, j, i])){
          dummy[i, j] <- NA
        }
      }
    }
    
    nonNAid_list_in <- vector("list", npt)
    
    for(j in seq_len(npt)){
      nonNAid_list_in[[j]] <- seq_len(m)[!is.na(dummy[, j])]
    }
    
    if(q > 1){
      theta_trans <- c(log(zeta_vec), log(nu2_vec), log(xi_vec), log(d_vec), l_vec)
    } else {
      theta_trans <- c(log(zeta_vec), log(nu2_vec), log(xi_vec), log(d_vec))
    }
    
    time_filter_here <- system.time(
      llh_OU(theta_trans, data_split_in, nonNAid_list_in)
    )
    
    time_filter_vec[sim_id] <- time_filter_here[3]
  }
  
  filter_time <- c(filter_time, mean(time_filter_vec))
}

#----------------------------------------------------------------------------------------------------------------------------------------
# Extrapolate univariate-treatment time for larger sample sizes
#----------------------------------------------------------------------------------------------------------------------------------------

mgrid_small <- c(10, 20, 50, 100, 200)

mgrid2 <- mgrid_small^2
mgrid3 <- mgrid_small^3

fit <- lm(filter_time_unitrt ~ mgrid_small + mgrid2 + mgrid3)

new_vec <- c(500, 1000, 10000, 100000, 1000000)

for(i in seq_along(new_vec)){
  
  new <- new_vec[i]
  new2 <- new^2
  new3 <- new^3
  
  new_data <- data.frame(
    mgrid_small = new,
    mgrid2 = new2,
    mgrid3 = new3
  )
  
  filter_time_unitrt <- c(
    filter_time_unitrt,
    predict(fit, newdata = new_data)
  )
}

save(
  filter_time,
  file = file.path(output_dir, "filter_time.Rdata")
)

save(
  filter_time_unitrt,
  file = file.path(output_dir, "filter_time_unitrt.Rdata")
)

#----------------------------------------------------------------------------------------------------------------------------------------
# Plot and table
#----------------------------------------------------------------------------------------------------------------------------------------

load(file.path(output_dir, "filter_time.Rdata"))
load(file.path(output_dir, "filter_time_unitrt.Rdata"))

mgrid <- c(10, 20, 50, 100, 200, 500, 1000, 10000, 100000, 1000000)

filter_time_data <- data.frame(
  time = c(filter_time_unitrt[1:6], filter_time[1:6]),
  m = c(mgrid[1:6], mgrid[1:6]),
  method = c(
    rep("univariate treatment", length(mgrid[1:6])),
    rep("new algorithm", length(mgrid[1:6]))
  )
)

ggplot(filter_time_data, aes(x = m, y = time, col = method)) +
  ggtitle("Filtering time comparison") +
  geom_line() +
  geom_point() +
  xlab("number of subjects") +
  ylab("time (seconds)")

filter_time_table <- data.frame(
  m = mgrid,
  univariate_treatment = filter_time_unitrt,
  new_algorithm = filter_time
)

print(
  xtable(filter_time_table, digits = c(0, 0, 4, 4)),
  include.rownames = FALSE
)