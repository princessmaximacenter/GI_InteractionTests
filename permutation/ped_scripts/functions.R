### ============================================================================#
### PERMUTATION TEST: General functions to perform permutation test
### Created by Saman Amini & Josephine Daub
### ============================================================================#

## import libraries
library(data.table)
library(reshape2)
library(vegan)

# ######################
crossprod_coocc <- function(perm_M) {
    # generate crossprod of the matrix to get concurrence.
    crossprod_M <- tcrossprod(perm_M)
    # remove the upper triangle of the matrix and diagonal
    crossprod_M[upper.tri(crossprod_M, diag = TRUE)] <- NA
    crossprod_M_low_melt <- melt(crossprod_M, na.rm = TRUE) # melt the matrix ;)
    colnames(crossprod_M_low_melt) <- c("event1", "event2", "cooc")
    return(crossprod_M_low_melt)
}

# ######################
conting_table <- function(x, total_mut, n_samples) {
    event1 <- x[1]
    event2 <- x[2]
    cooc <- as.numeric(x[3])
    event1_not <- total_mut[event1] - cooc
    not_event2 <- total_mut[event2] - cooc
    neither <- n_samples - sum(event1_not, not_event2, cooc)

    row <- c(event1, event2, cooc, not_event2, event1_not, neither,
             sum(cooc, event1_not, not_event2, neither))
}

########################
# function to perform permutation for a given matrix
perform_perm <- function(alteration_M, strata_vector, n_perm = 100, seed, 
                         strata_smpl=F) {
    set.seed(seed)
    print("start permutation")
    start_time <- Sys.time()
    if (strata_smpl){
      alteration_M<-t(alteration_M)
    }
    if (length(strata_vector) != ncol(alteration_M)) {
      stop("ncol of alteration matrix and strata should have the same length")
    }
    # perform permutation
    perm <- permatswap(t(alteration_M), strata = strata_vector, 
                       mtype = "prab", times = n_perm)
    # set rownames of permuted matrices
    if (strata_smpl){
      perm$perm <- lapply(perm$perm, function(x) {
        x<-t(x)
        rownames(x) <- rownames(alteration_M)
        return(x)
      })
    } else {
      perm$perm <- lapply(perm$perm, function(x) {
        rownames(x) <- colnames(alteration_M); x
        })
    }
   
    end_time <- Sys.time()
    calc_time <- end_time - start_time
    print(calc_time)
    return(perm)
}


# ########################
# function to process rds file to calculate empirical P values.
emp_P_process_rds <- function(co_real, RDSfile) {
  real_DF <- cbind(co_real, "co_perm" = 0, "me_perm" = 0, "co_P" = NA,
                   "me_P" = NA, "co_P_adj" = NA, "me_P_adj" = NA)

    perm_M_list <- readRDS(RDSfile)
    n_per_file <- length(perm_M_list$perm)
    print(paste("working on file: ", RDSfile, sep = ""))

    for (i in 1:n_per_file) {
      perm_M <- perm_M_list$perm[[i]]
      perm_co_M <- crossprod_coocc(perm_M)
      co_idx <- as.numeric(perm_co_M$cooc) >= as.numeric(real_DF$cooc)
      me_idx <- as.numeric(perm_co_M$cooc) <= as.numeric(real_DF$cooc)

      real_DF$co_perm[co_idx] <- real_DF$co_perm[co_idx] + 1
      real_DF$me_perm[me_idx] <- real_DF$me_perm[me_idx] + 1
    }
    print(paste("file: ", RDSfile, " DONE!",sep = ""))

  return(real_DF)
}

# ########################
# new function to calculate empirical P values.
empirical_P <- function(co_real, Rfiles, n = 1000000) {
  real_DF <- cbind(co_real, "co_perm" = 0, "me_perm" = 0, "co_P" = NA,
                   "me_P" = NA, "co_P_adj" = NA, "me_P_adj" = NA)

  for (file in Rfiles) {
    print(paste("working on file: ", file, sep = ""))
    load(file)
    real_DF$co_perm <- real_DF$co_perm + co_real_perm_P$co_perm
    real_DF$me_perm <- real_DF$me_perm + co_real_perm_P$me_perm
    print(paste("file: ", file, " DONE!",sep = ""))
  }
  # use n + 1 for preventing 0 as p_values:
  real_DF$co_P <- (real_DF$co_perm + 1) / (n + 1)
  real_DF$me_P <- (real_DF$me_perm + 1) / (n + 1)
  real_DF$co_P_adj <- p.adjust(real_DF$co_P, method = "BH")
  real_DF$me_P_adj <- p.adjust(real_DF$me_P, method = "BH")
  return(real_DF)
}


# function to exclude alterations with very low or very high frequency
threshold_alt_freq <- function(data_frame, threshold = 1) {
    data_frame_min <- data_frame[,colSums(data_frame) > threshold]
    data_frame_min_max <-
      data_frame_min[,colSums(data_frame_min) < nrow(data_frame) - threshold]
    data_frame_min_max <- data_frame_min_max[rowSums(data_frame_min_max) >= 1,]
    return(data_frame_min_max)
}

# Empirical p distribution
empirical_P_dist <- function(co.mtx, RDSfile) {

  cat("Reading RDS file ", RDSfile, "\n", sep="")

  perm_M_list <- readRDS(RDSfile)
  co_p_DF<-list()

  for (j in 1:length(perm_M_list$perm)) {
    perm_M <- perm_M_list$perm[[j]]
    perm_co_M <- crossprod_coocc(perm_M)
    cooc_num_perm<-as.numeric(perm_co_M$cooc)

    for (co in 1:ncol(co.mtx)){
      if (j==1) {
        co_p_DF[[co]] <- cbind("cooc"=co.mtx[,co], "co_perm" = 0,
                               "me_perm" = 0)
      }
      cooc_num<-as.numeric(co_p_DF[[co]][,"cooc"])

      co_idx <- cooc_num_perm >= cooc_num
      me_idx <- cooc_num_perm <= cooc_num

      co_p_DF[[co]][co_idx,"co_perm"] <- co_p_DF[[co]][co_idx,"co_perm"] + 1
      co_p_DF[[co]][me_idx,"me_perm"] <- co_p_DF[[co]][me_idx,"me_perm"] + 1
    }
    cat("matrix ", j, " done!\n", sep="")
  }
  return(co_p_DF)
}

################################################################################
# FUNCTION GetEmpericalFDR (p.obs, p.null)
# Calculate the empirical fdr
# p.obs: (named) vector with observed p-values
# p.null: vector with (pooled) p-values from randomizations OR
#         list with N vectors with p-values from N randomizations
################################################################################

GetEmpericalFDR <- function(p.obs, p.null, use.uniq=T) {

    # sort p.obs and save the original order
    ix.sort <- order(p.obs)
    p.obs <- p.obs[ix.sort]

    # remove p-values lower than 0 and larger than 1
    ix.p <- which(p.obs >= 0 & p.obs <= 1)
    p.obs.flt <- p.obs[ix.p]

    # take unique list of p-values to increase speed
    if (use.uniq){
      p.obs.uniq <- unique(p.obs.flt)
      n.p.obs.uniq<-length(p.obs.uniq)
      t<-table(p.obs.flt)
      t.cs.obs<-cumsum(t)

      # get table and running sum of table of p.null
      p.null.flt <- p.null[p.null >= 0 & p.null <= 1]
      n.p.null.flt<-length(p.null.flt)
      t<-table(p.null.flt)
      t.cs<-cumsum(t)
      t.cs<-c(0,t.cs)
      p.null.uniq<-sort(unique(p.null.flt))
      p.null.uniq<-c(0,p.null.uniq)
      n.p.null.uniq<-length(p.null.uniq)

      # create table of p.null.uniq indices linked to p.obs.uniq
      # each entry i gives the index of the largest p.null.uniq element
      # that is equal or smaller than p.obs.uniq[i]
      p.obs.null.ix<-rep(0,n.p.obs.uniq)
      i<-1
      j<-1
      while (i <= n.p.obs.uniq & j <= n.p.null.uniq){
        # for each p.obs, get the index of the highest p.null <= p.obs
        if (p.null.uniq[j]<=p.obs.uniq[i]) {
          # ok, increase j
          j <- j+1
        } else {
          # assign j-1 to index p obs, and go to next p obs
          p.obs.null.ix[i]<-(j-1)
          i <- i+1
        }
      }
      p.obs.null.ix[p.obs.null.ix==0]<-n.p.null.uniq

    } else {
      p.obs.uniq <- p.obs.flt
    }

    # for now, set the proportion of true H0's = 1 (conservative)
    pi.0 <- 1
    nsc <- length(p.obs.flt)

    # get ratio #p-values randomized tests/#p-values in real tests
    fdr.est1 <- numeric()
    for (p in 1:length(p.obs.uniq)) {
        pval <- p.obs.uniq[p]

        if (is.list(p.null)) {
            # calculate running mean of N0
            # (the proportion of random p-values <= observed pval)
            nruns <- (length(p.null))
            n <- 0
            mean.N0 <- 0
            for (run in seq(nruns)) {
                p.null.flt <- p.null[[run]]
                p.null.flt <- p.null.flt[p.null.flt >= 0 & p.null.flt <= 1]
                nsh <- length(p.null.flt)
                N0 <- sum(p.null.flt <= pval) / nsh
                n <- n + 1
                delta <- N0 - mean.N0
                mean.N0 <- mean.N0 + delta / n
            }
        } else {
            # calculate mean of N0
            # (the proportion of random p-values <= observed pval)
            if (!use.uniq) {
              mean.N0 <- sum(p.null.flt <= pval) / n.p.null.flt
            } else {
              mean.N0 <- t.cs[p.obs.null.ix[p]] / n.p.null.flt
            }
        }

        # calculate the number of rejected hypothesis if pval is threshold
        if (!use.uniq) {
          N.t <- sum(p.obs.flt <= pval)
        } else {
          N.t <- t.cs.obs[p]
        }
        # calculate first estimate of FDR
        fdr.est1[p] <- pi.0*nsc*mean.N0/N.t
    }

    # get the smallest ratio per p-value
    fdr.est2 <- numeric()
    for (i in 1:n.p.obs.uniq) {
        fdr.est2[i] <- min(fdr.est1[i:n.p.obs.uniq])
    }

    # combine results with p-values
    fdr.est3 <- data.frame(p = p.obs, fdr.est = rep(-1, length(p.obs)))
    m<-match(p.obs.flt, p.obs.uniq)
    fdr.est3[ix.p,2]<-fdr.est2[m]

    # put back in original order
    fdr.est3 <- fdr.est3[match(1:length(p.obs),ix.sort),]

    return(fdr.est3)
}


