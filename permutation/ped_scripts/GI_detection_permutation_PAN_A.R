### ============================================================================
### PERMUTATION TEST - STEP 2A: Generate permuted matrices
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
        Script to create permuted matrices for PAN cancer permutation test

        Arguments:
        --arg1=source of dataset               -project name, e.g. dkfz
        \n")
         q(save = "no")
}
cancer_type      <- "PAN"
source           <- args[1]

print(cancer_type)
# import libraries
source("./ped_scripts/functions.R")

f<-file.path("./ped_results", paste0(source,"_CTs_list.txt"))
ct.tab<-read.delim(f,header = F)
n_cancer_type <- nrow(ct.tab)

path_dataset <- file.path("./ped_data/", source, "matrix_gene_sample.txt")

# read and clean data
alteration <- read.delim(path_dataset, row.names = 1, check.names = FALSE)
# exclude last few columns since they include cancer types
alteration_ct <- alteration[,1:(ncol(alteration) - n_cancer_type)]

# exclude alterations 0,1 and 100%, exclude samples with 0 alterations
alteration_perm <- threshold_alt_freq(alteration_ct, 1)
message(paste(nrow(alteration_perm), " samples and ", ncol(alteration_perm),
              " alterations are included", sep = ""))
try(if ((nrow(alteration_perm) < 2) & (ncol(alteration_perm) < 2)) {
  stop("too few samples and alterations")
})

alteration_perm_M <- as.matrix(alteration_perm)

ct_strata_raw<-rep(0, nrow(alteration))
for (i in 1:n_cancer_type){
  ct_strata_raw[alteration[, ncol(alteration)-n_cancer_type+i]=="YES"]<-i
}
m<-match(dimnames(alteration_perm_M)[[1]],
         rownames(alteration))
ct_strata<-ct_strata_raw[m]

# filter out cancer types with only one patient
t<-table(ct_strata)
ix<-names(t[t>1])
ix<-which(ct_strata %in% ix)
alteration_perm_M<-alteration_perm_M[ix,]
ct_strata<-as.numeric(as.factor(ct_strata[ix]))

# Save matrix and strata to use in step B (parallelized)
out_path <- file.path("./ped_results/", source)
dir.create(out_path,showWarnings = F, recursive = T)
save(alteration_perm_M, ct_strata, file=file.path(out_path, "PAN_obj.R"))
