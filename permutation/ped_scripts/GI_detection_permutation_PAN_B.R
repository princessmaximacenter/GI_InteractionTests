### ============================================================================
### PERMUTATION TEST - STEP 2B: Generate permuted matrices
### Saman Amini
### Changed on 7 May 2018 by Josephine Daub to allow other sources
### Change on 19 June 2018 by Josephine Daub to apply to PAN cancer test
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
        Genetic interaction detection R script
        script to perform permutation test for a matrix of individual and alteration data.

        Arguments:
        --arg2=seed number                -the seed number of the batch; 
                                           if the total number of permutations 
                                           is set to 1M, you can run the script 
                                           in parallel with seed=1 to seed=200
                                           and each batch having 5000 permutations
        --arg3=source of dataset          -project name, e.g. dkfz
        --arg4=number of permutations     -default: 50, set to low number in
                                           case of memory issues (and increase
                                           number of batches)

        \n")
         q(save = "no")
}
cancer_type   <- "PAN"
seed          <- as.numeric(args[1])
source        <- args[2]
n_perm        <- ifelse(length(args) < 3, 50, as.numeric(args[3]))

print(cancer_type)
print(seed)
# import libraries
source("./ped_scripts/functions.R")

#####################################################
## PART B: perform permutation.
#####################################################
out_path <- file.path("./ped_results/", source)
load(file=file.path(out_path, "PAN_obj.R"))
permuted_alts <- perform_perm(alteration_perm_M, ct_strata,
                                n_perm = n_perm, seed, strata_smpl=T)

dir.create(out_path, showWarnings = FALSE, recursive = T)
out_path_CT <- file.path(out_path, "PAN")
dir.create(out_path_CT, showWarnings = FALSE, recursive = T)
file_name <- paste0(paste("/perm_seed", seed, "PAN", sep = "_"), ".rds")
saveRDS(permuted_alts, file = paste0(out_path_CT, file_name))
