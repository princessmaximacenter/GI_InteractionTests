### ============================================================================
### PERMUTATION TEST - STEP 4C: Generate p-value null distribution
# Merge all data and write p-values to file
# Created by Saman Amini & Josephine Daub
### ============================================================================

# Collect arguments
args <- commandArgs(TRUE)

# Default setting when no argument passed
if (length(args) < 1) {
    args <- c("--help")
}

# Help section
if ("--help" %in% args) {
    cat("
        FDR correction of p-values based on a generated empirical p-values.
        Step C: merge all data and write p-values to file
        Arguments:
        --arg1=cancer type                   -name of cancer type
        --arg2=source of dataset             -project, e.g. dkfz
        --arg3=number of RDS files           -default = 5
        --arg4=number of sampled mtxs        -default = 25
        --arg5=number of permuted mtxs       -default = 250
        --arg6=startnr of sampled mtxs       -default = 1
        --arg7=endnr of sampled mtxs         -default = 25
        \n")
    q(save = "no")
}

cancer_type  <- args[1]
d_source     <- args[2]
n_RDS_files  <- ifelse(length(args) < 3, 5, as.numeric(args[3]))
n_smpl_mtx   <- ifelse(length(args) < 4, 25, as.numeric(args[4]))
n_perm_mtx   <- ifelse(length(args) < 5, 250, as.numeric(args[5]))
n_start_mtx  <- ifelse(length(args) < 6, 1, as.numeric(args[6]))
n_stop_mtx   <- ifelse(length(args) < 7, 25,  as.numeric(args[7]))

cat("ct: ", cancer_type, ", project: ", d_source, 
    ", n_RDS_files: ", n_RDS_files, ", n_smpl_mtx: ", n_smpl_mtx, 
    ", n_perm_mtx:", n_perm_mtx,  ", n_start_mtx: ", n_start_mtx, 
    ", n_stop_mtx: ", n_stop_mtx, 
    "\n", sep="")
                       
options(scipen = 100)

# import libraries
library(gtools)
source("./ped_scripts/functions.R")

mtx_path <- file.path("./ped_results", d_source, "P",cancer_type)

null_dist_out_path <- file.path("./ped_results/",
                                d_source, "/emp_p_results/null_dist/")

ct_out_path <- paste(null_dist_out_path, cancer_type, sep = "")
dir.create(file.path(ct_out_path), showWarnings = FALSE, recursive = TRUE)


for (seed in 1:n_RDS_files){
  load(file=file.path(mtx_path, paste0("perm_co_mtx_",seed,".R")))
  if (seed==1){
    co_p_DF <- co_real_perm_P
  } else {
    co_p_DF<-Map('+',co_p_DF[n_start_mtx:n_stop_mtx],
                 co_real_perm_P[n_start_mtx:n_stop_mtx])
  }
  cat("finished processing RDS file ", seed, "\n", sep="")
}

for (k in n_start_mtx:n_stop_mtx){
  # It used to be (in accordance with the calculation of the observed p-values):
  # v_p_values_CO <- (co_p_DF[[k]][,"co_perm"] + 1) / (n_perm_mtx + 1)
  # However, as we go through all permuted matrices, the random sampled matrix
  # is one of them and so is already counted
  # We keep the division by (n_perm_mtx+1) to have the same divider as the 
  # observed p-values
  v_p_values_CO <- (co_p_DF[[k]][,"co_perm"]) / (n_perm_mtx + 1)
  v_p_values_ME <- (co_p_DF[[k]][,"me_perm"]) / (n_perm_mtx + 1)

  file_ct_out_path_CO <- paste(ct_out_path, "/co_seed_", k , ".txt", sep = "")
  file_ct_out_path_ME <- paste(ct_out_path, "/me_seed_", k , ".txt", sep = "")

  write.table(v_p_values_CO, file_ct_out_path_CO, quote = FALSE,
              row.names = FALSE, col.names = FALSE)
  write.table(v_p_values_ME, file_ct_out_path_ME, quote = FALSE,
              row.names = FALSE, col.names = FALSE)
}
