rm(list = ls())

#----------------------------------------------------------------------------------------------------------------------------------------
# Collect bootstrap estimate matrices
#----------------------------------------------------------------------------------------------------------------------------------------

project_dir <- normalizePath(".", mustWork = TRUE)

bootstrap_dir <- file.path(project_dir, "bootstrap_results", "results_m1000")

n_boot_files <- 100

theta_mat_est_extra_final <- vector("list", n_boot_files)

for(i in seq_len(n_boot_files)){
  
  input_file <- file.path(bootstrap_dir, paste0("theta_est_mat", i, ".Rdata"))
  
  if(!file.exists(input_file)){
    stop(paste("Missing bootstrap file:", input_file))
  }
  
  load(input_file)
  
  theta_mat_est_extra_final[[i]] <- theta_est_mat1
  
  rm(theta_est_mat1)
}

save(
  theta_mat_est_extra_final,
  file = file.path(project_dir, "theta_mat_est_extra_final.Rdata")
)

length(theta_mat_est_extra_final)
dim(theta_mat_est_extra_final[[1]])