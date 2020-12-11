### ============================================================================
### PERMUTATION TEST - STEP 3B: Calculate empirical P
# For each RDS file, get cooccurrence counts for each of its permuted matrices 
# and compare with observed counts
# Created by Saman Amini
# Changed by Josephine Daub on May 7, 2015 to allow other sources
### ============================================================================

# Collect arguments
args <- commandArgs(TRUE)

# Default setting when less than 2 argument passed
if (length(args) < 2) {
    args <- c("--help")
}

# Help section
if ("--help" %in% args) {
    cat("
        Genetic interaction detection: empirical p_value calculation
        R script to calculate permuted p values

        Arguments:
        --arg1=cancer type                   -name of cancer type
        --arg2=source of dataset             -project, e.g. dkfz
        --arg3=number of RDS file
        \n")
         q(save = "no")
}

cancer_type      <- args[1]
d_source         <- args[2]
n_rds            <- args[3]
print(cancer_type)

# import libraries
source("./ped_scripts/functions.R")
library(gtools)
options(scipen = 100)

#############################
# read permuted data
perm_path <- file.path("./ped_results", d_source)
perm_path_CT <- file.path(perm_path, cancer_type)

# path to store intermediate data
perm_path_tmp <- file.path("./ped_results/", d_source, "/T")
perm_path_CT_tmp <- file.path(perm_path_tmp, cancer_type)
dir.create(perm_path_CT_tmp, showWarnings = F, recursive = T)

RDSfiles <- list.files(perm_path_CT, full.names = TRUE)
RDSfile <- file.path(perm_path_CT, paste0("perm_seed_", n_rds, "_", 
                                          cancer_type, ".rds"))

#############################
# get real observations
in_path <- file.path("./ped_results",d_source, "emp_p_results/")

# load co_real (created in step A)
load(file = file.path(in_path, paste0(cancer_type, ".R")))

# Run for one RDS file and save this file
co_real_perm_P <- emp_P_process_rds(co_real, RDSfile)
save(co_real_perm_P, file = file.path(perm_path_CT_tmp, paste0(n_rds,".R")))
