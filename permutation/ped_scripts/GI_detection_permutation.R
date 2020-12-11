### ============================================================================
### PERMUTATION TEST - STEP 2: Generate permuted matrices
### Saman Amini
### Changed on 7 May 2018 by Josephine Daub to allow other sources
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
        Script to create permuted matrices for permutation test

        Arguments:
        --arg1=cancer type                -name of cancer type
        --arg2=seed number                -the seed number of the batch; 
                                           if the total number of permutations 
                                           is set to 1M, you can run the script 
                                           in parallel with seed=1 to seed=200
                                           and each batch having 5000 permutations
        --arg3=source of dataset          -project name, e.g. dkfz
        --arg4=number of permutations     -default: 50, set to lower in
                                           case of memory issues (and increase
                                           number of batches)

        \n")
         q(save = "no")
}

cancer_type   <- args[1]
seed          <- as.numeric(args[2])
source        <- args[3]
n_perm        <- ifelse(length(args) < 4, 50, as.numeric(args[4]))

cat("ct: ", cancer_type, ", seed: ", seed, ", project: ", source, ", n_perm: ", 
    n_perm, "\n", sep="")

# import libraries
source("./ped_scripts/functions.R")

f<-file.path("./ped_results", paste0(source,"_CTs_list.txt"))
ct.tab<-read.delim(f,header = F)
n_cancer_type <- nrow(ct.tab)

path_dataset <- file.path("./ped_data/", source, "matrix_gene_sample.txt")

# read and clean data
alteration <- read.delim(path_dataset, row.names = 1, check.names = FALSE)
alteration_perm <- alteration[alteration[,cancer_type] == "YES",]
# exclude last few columns since they include cancer types
alteration_perm <- alteration_perm[,1:(ncol(alteration_perm) - n_cancer_type)]

# exclude alterations 0,1 and 100%, exclude samples with 0 alterations
alteration_perm <- threshold_alt_freq(alteration_perm, 1)
message(paste(nrow(alteration_perm), " samples and ", ncol(alteration_perm),
              " alterations are included", sep = ""))
try(if ((nrow(alteration_perm) < 2) & (ncol(alteration_perm) < 2)) {
    stop("too few samples and alterations")
    })
alteration_perm_M <- as.matrix(alteration_perm)

# define strata, a vector that contains the same integer number per alteration
vector_alts <- unlist(replicate(ncol(alteration_perm_M), 1))

#####################################################
## perform permutation.
#####################################################
permuted_alts <- perform_perm(alteration_perm_M, vector_alts, 
                              n_perm = n_perm, seed)

out_path <- file.path("./ped_results/", source)
dir.create(out_path, showWarnings = FALSE, recursive = T)
out_path_CT <- file.path(out_path, cancer_type)
dir.create(out_path_CT, showWarnings = FALSE, recursive = T)
file_name <- paste0(paste("/perm_seed", seed, cancer_type, sep = "_"), ".rds")
saveRDS(permuted_alts, file = paste0(out_path_CT, file_name))
