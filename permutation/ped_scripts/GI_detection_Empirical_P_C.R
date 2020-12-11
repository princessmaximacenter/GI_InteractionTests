### ============================================================================
### PERMUTATION TEST - STEP 3C: Calculate empirical P
# Calculate empirical p-value by adding up all comparisons from step B
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
        --arg3=number of permutations        -default: 1000
        \n")
         q(save = "no")
}

cancer_type   <- args[1]
d_source      <- args[2]
n_perm        <- ifelse(length(args) < 3, 1000, as.numeric(args[3]))
print(cancer_type)

# import libraries
source("./ped_scripts/functions.R")
library(gtools)
options(scipen = 100)

#############################
# read permuted data
perm_path <- file.path("./ped_results/", d_source, "/T")
perm_path_CT <- file.path(perm_path, cancer_type)
Rfiles <- list.files(perm_path_CT, full.names = TRUE)
Rfiles <- Rfiles[mixedorder(Rfiles)]

#############################
# load co_real (created in step A)
in_path <- file.path("./ped_results", d_source, "emp_p_results/")

load(file = file.path(in_path, paste0(cancer_type, ".R")))
# calculate empirical p_values
co_real_perm_P <- empirical_P(co_real, Rfiles, n = n_perm) 

# write results
out_path <- file.path("./ped_results", d_source, "emp_p_results/")
dir.create(out_path, showWarnings = F, recursive = T)
file_name <- file.path(out_path, paste0(cancer_type,".txt"))
write.table(co_real_perm_P, file_name, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)
