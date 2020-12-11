### ============================================================================
### PERMUTATION TEST - STEP 3A: Calculate empirical P
# Get observed co-occurrence counts
# Created by Saman Amini
# Changed by Josephine Daub on May 7, 2015 to allow other sources
## ============================================================================

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
        R script to calculate emprical p values

        Arguments:
        --arg1=cancer type                   -name of cancer type
        --arg2=source of dataset             -project, e.g. dkfz
        --arg3=number of RDS files           -default = 20
        \n")
         q(save = "no")
}

cancer_type  <- args[1]
d_source     <- args[2]
n_rds        <- ifelse(length(args) < 3, 20, as.numeric(args[3]))
print(cancer_type)

# import libraries
source("./ped_scripts/functions.R")
library(gtools)
options(scipen = 100)

# global var
resdir<-"./ped_results/"
datadir<-"./ped_data/"
cancers<-read.delim(file.path(resdir,paste0(d_source,"_CTs_list.txt")),
                    header=F, stringsAsFactors = F)[,1]
n_cancer_type<-length(cancers)

#############################
# read permuted data
perm_path <- file.path(resdir, d_source)
perm_path_CT <- file.path(perm_path, cancer_type)
RDSfiles <- list.files(perm_path_CT, full.names = TRUE)
try(if (length(RDSfiles) < n_rds) stop(paste0("lower than ",
    n_rds, " rds files/folder is not allowed")))
RDSfiles <- RDSfiles[mixedorder(RDSfiles)]

# read real data
path_dataset <- file.path(datadir,d_source,"matrix_gene_sample.txt")

alteration <- read.delim(path_dataset, row.names = 1, check.names = FALSE)
if (cancer_type=="PAN"){
  alt_ct <- alteration
} else {
  alt_ct <- alteration[alteration[,cancer_type] == "YES",]
}

# exclude last few columns since they are not needed for now.
alt_ct_clean <- alt_ct[,1:(ncol(alt_ct) - n_cancer_type)]

# exclude alterations that happen only in one sample;
# exclude alterations that occure in 100 percent of the sample
# exclude samples that don't have any mutations
alt_ct_clean <- threshold_alt_freq(alt_ct_clean, 1)
alt_ct_clean_M <- as.matrix(alt_ct_clean)

# if PAN cancer: exclude cancer types with only one patient
# you can use alt_ct, because in PAN cancer test all samples have
# at least one alteration, so none of them are removed from previous step
# and nrow hasn't changed
if (cancer_type=="PAN"){
  ct_strata_raw<-rep(0, nrow(alt_ct))
  for (i in 1:n_cancer_type){
    ct_strata_raw[alt_ct[, ncol(alt_ct)-n_cancer_type+i]=="YES"]<-i
  }
  m<-match(dimnames(alt_ct_clean_M)[[1]], rownames(alt_ct))
  ct_strata<-ct_strata_raw[m]

  # filter out cancer types with only one patient
  t<-table(ct_strata)
  ix<-names(t[t>1])
  ix<-which(ct_strata %in% ix)
  alt_ct_clean_M<-alt_ct_clean_M[ix,]
}

# generate coocurence for real data.
crossprod_matrix_real <- crossprod_coocc(t(alt_ct_clean_M))

total_per_alteration <- colSums(alt_ct_clean)
n_samples <- nrow(alt_ct_clean)

# using "conting_table" to generate contingency table
conting_matrix <- t(apply(crossprod_matrix_real, 1,
                          conting_table, total_mut = total_per_alteration,
                          n_samples = n_samples))

colnames(conting_matrix) <- c("event1", "event2", "cooccurrence", "not_event2",
                              "event1_not", "neither", "n_samples")

co_real <- as.data.frame(conting_matrix, stringsAsFactors = FALSE)

# write results
out_path <- file.path(resdir, d_source,"emp_p_results")

dir.create(out_path, showWarnings = FALSE, recursive = T)

# save co_real to be used in step B
save(co_real, file = file.path(out_path, paste0(cancer_type, ".R")))
