### ============================================================================
### PERMUTATION TEST - STEP 5: FDR estimation
# Calculates the empirical FDR given the observed p-values and a 
# null distribution of p-values
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
        FDR correction of p-values based on a generated empirical p-values 
        (null distribution).
        Arguments:
        --arg1=cancer type                   -name of cancer type
        --arg2=source of dataset             -project, e.g. dkfz
        \n")
    q(save = "no")
}

cancer_type      <- args[1]
d_source         <- args[2]

source("./ped_scripts/functions.R")

null_path <- file.path("./ped_results", d_source, "emp_p_results/null_dist")
null_path_ct <- file.path(null_path, cancer_type)
co_files <- list.files(null_path_ct, full.names = TRUE, pattern = "\\co*")
me_files <- list.files(null_path_ct, full.names = TRUE, pattern = "\\me*")

co_null_dist <- unlist(lapply(co_files, read.delim, header = FALSE))
me_null_dist <- unlist(lapply(me_files, read.delim, header = FALSE))

# import results with uncorrected p values
p_results_path <- file.path("./ped_results/", d_source, "/emp_p_results")
p_results_path_ct <- file.path(p_results_path, paste0(cancer_type, ".txt"))
p_results_ct <- read.delim(p_results_path_ct)

p_results_ct_emp_FDR <- cbind(p_results_ct,
                              "emp_FDR_CO" = GetEmpericalFDR(p_results_ct$co_P,
                                                             co_null_dist)[,2])
p_results_ct_emp_FDR <- cbind(p_results_ct_emp_FDR,
                              "emp_FDR_ME" = GetEmpericalFDR(p_results_ct_emp_FDR$me_P,
                                                             me_null_dist)[,2])
file_name <- file.path(p_results_path, paste0(cancer_type, "_emp_FDR.txt"))
write.table(p_results_ct_emp_FDR, file_name, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)
