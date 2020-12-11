### ============================================================================
### PERMUTATION TEST - STEP 4A: Generate p-value null distribution
# Randomly sample N permuted matrices
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
        Step A: randomly pick a matrix, repeat N times
        Arguments:
        --arg1=cancer type                   -name of cancer type
        --arg2=source of dataset             -project, e.g. dkfz
        --arg3=nb of random mtxs             -default=25
        \n")
    q(save = "no")
}

cancer_type  <- args[1]
d_source     <- args[2]
nperm        <- ifelse(length(args) < 3, 25, as.numeric(args[3]))

options(scipen = 100)

cat("ct: ", cancer_type, ", project: ", d_source, ", n_perm: ", 
    n_perm, "\n", sep="")

# import libraries
library(gtools)
source("./ped_scripts/functions.R")

mtx_path <- file.path("./ped_results", d_source, "P", cancer_type)
dir.create(mtx_path,showWarnings = F, recursive = T)

perm_path <- file.path("./ped_results", d_source)

path_data <- file.path(perm_path, cancer_type)
RDSfiles <- list.files(path_data, full.names = TRUE)
RDSfiles <- RDSfiles[!file.info(RDSfiles)$isdir] # only files.
RDSfiles <- RDSfiles[mixedorder(RDSfiles)]

n.RDSfiles<-length(RDSfiles)
RDSfile <- readRDS(RDSfiles[1])
n.mtx_per_RDS <- length(RDSfile$perm)

# Randomly select nperm matrices
# (To be able to compare with former versions of the script,
# the seed is set after each sampling)
df_r<-data.frame(r_file=rep(0,nperm), r_matrix=rep(0,nperm))
for (i in 1:nperm){
  set.seed(i)
  df_r$r_file[i]<-sample(1:n.RDSfiles, 1)
  df_r$r_matrix[i]<-sample(1:n.mtx_per_RDS,1)
}
df_r<-df_r[order(df_r$r_file, df_r$r_matrix),]

f_old <- 1
for (i in 1:nperm){
  f<-df_r$r_file[i]
  m<-df_r$r_matrix[i]
  cat("i: ", i , " f:", f , "m: ", m, "\n")
  if (f!=f_old){
    # read in the randomly selected matrix
    RDSfile <- readRDS(RDSfiles[f])
  }
  selected_matrix <- RDSfile$perm[[m]]

  if (i==1){
    n<-nrow(selected_matrix)
    crossprod_perm_matrix<-matrix(0,nrow = n*(n-1)/2,ncol = nperm)
  }
  crossprod_perm_matrix[,i]<-crossprod_coocc(selected_matrix)$cooc
  f_old <- f
}
save(crossprod_perm_matrix, file = file.path(mtx_path,"perm_mtx.R"))
