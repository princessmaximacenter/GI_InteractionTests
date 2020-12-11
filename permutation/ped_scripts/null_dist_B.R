### ============================================================================
### PERMUTATION TEST - STEP 4B: Generate p-value null distribution
# Compare sampled matrices with matrices in RDS file
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
        Step B: (parallel) compare sampled mtxs with mtxs in RDS file

        Arguments:
        --arg1=cancer type                   -name of cancer type
        --arg2=seed                          -provided seed
        --arg3=source of dataset             -project, e.g. dkfz
        \n")
    q(save = "no")
}

cancer_type      <- args[1]
seed             <- as.numeric(args[2])
d_source         <- args[3]

cat("ct: ", cancer_type, ", project: ", d_source, ", seed: ", 
    seed, "\n", sep="")

options(scipen = 100)

# import libraries
library(gtools)
source("./ped_scripts/functions.R")

mtx_path <- file.path("./ped_results",  d_source, "P",cancer_type)
dir.create(mtx_path,showWarnings = F, recursive = T)

perm_path <- file.path("./ped_results/", d_source)

path_data <- file.path(perm_path, cancer_type)
RDSfiles <- list.files(path_data, full.names = TRUE)
RDSfiles <- RDSfiles[!file.info(RDSfiles)$isdir] # only files.
RDSfiles <- RDSfiles[mixedorder(RDSfiles)]

load(file = file.path(mtx_path,"perm_mtx.R"))
RDSfile<-RDSfiles[seed]
co_real_perm_P <- empirical_P_dist(crossprod_perm_matrix, RDSfile)
save(co_real_perm_P, file=file.path(mtx_path, paste0("perm_co_mtx_",seed,".R")))
