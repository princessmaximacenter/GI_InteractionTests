### ============================================================================
### PERMUTATION TEST - STEP 1: Data preparation
### this file includes scripts for pre processing the pediatric data set.
### Author: Saman Amini
### Edited by Josephine Daub on 17 April 2018 to allow other projects 
### ============================================================================

args <- commandArgs(TRUE)

# Default setting when no argument passed
argsLen <- length(args);
if (argsLen < 1) {
    args <- c("--help")
}

# Help section
if ("--help" %in% args) {
    cat("
        data prep R script
        script to generate gene sample matrix 

        Arguments:
        --arg1=input file       -input file with four columns
                                 caseID, ct, ct_sub, gene    
        --arg2=project name     -e.g. dkfz
        --arg3=use sub          -use cancer subtype? T/F(default)
        --arg4=has header       -file has header? T(default)/F
        \n")
         q(save = "no")
}

file <- args[1]
if (argsLen<2){
  cat("project is missing\n")
  q(save = "no")
} else {
  project <- args[2]
}
use_sub <- ifelse(argsLen<3, FALSE, as.logical(args[3]))
has_header <- ifelse(argsLen<4, TRUE, as.logical(args[4]))

datadir <- "./ped_data"
resdir <- "./ped_results"
dir.create(file.path(resdir,project), showWarnings = F, recursive = T)
dir.create(file.path(datadir,project), showWarnings = F, recursive = T)

cat("Directories ", datadir, " and ", resdir, " created.\n")

# read in patient data
patient_alt <- read.delim(file, header = has_header, stringsAsFactors = FALSE)
colnames(patient_alt)<-c("caseID", "ct", "ct_sub", "gene")

#===============================================================================
# generate a matrix with all patients (rows) and alterations (columns)
#===============================================================================
m <- matrix(nrow = length(unique(patient_alt$caseID)),
            ncol = length(unique(patient_alt$gene)))
rownames(m) <- unique(patient_alt$caseID)
colnames(m) <- unique(patient_alt$gene)
patients <- rownames(m)

# fill in matrix for all patients
for (i in patients) {
    vars <- patient_alt$gene[patient_alt$caseID == i] # altered genes per patient
    m[i, colnames(m) %in% vars] <- 1
    m[i, is.na(m[i, ])] <- 0
}

cat("Mutation matrix created with ", nrow(m) ," samples and ",
    ncol(m), " genes.\n", sep="")

if (use_sub) patent_alt$ct<-patient_alt$ct_sub
cancers <- unique(patient_alt$ct)

write.table(as.data.frame(cancers),
            file.path(resdir,paste0(project, "_CTs_list.txt")),
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

cat(length(cancers), " Cancer types written to ", 
    file.path(resdir,paste0(project, "_CTs_list.txt")),".\n", sep="")

# add cancer types to the matrix
m <- as.data.frame(m)
ct_DF <- as.data.frame(matrix("NO", nrow = nrow(m), ncol = length(cancers)),
                       stringsAsFactors = FALSE)
colnames(ct_DF) <- cancers
m <- cbind(m, ct_DF)

for (i in cancers) {
    ind <- unique(patient_alt$caseID[patient_alt$ct == i])
    m[rownames(m) %in% ind, i] <- "YES"
}

of <- file.path(datadir,project,"matrix_gene_sample.txt")
write.table(m, of, quote = FALSE, row.names = TRUE, col.names = TRUE,
            sep = "\t")

cat("Mutation matrix written to ", of, ".\n", sep="")