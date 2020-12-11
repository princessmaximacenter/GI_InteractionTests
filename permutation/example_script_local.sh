#! /bin/bash

## PERMUTATION TEST EXAMPLES - LOCAL MACHINE

# Example scripts to test permutation test pipeline on a local machine
# We assume a null distribution of 250 permuted matrices, generated and
# processed in 5 batches of 50 matrices. The FDR estimation is then done
# on 25 randomly sampled permuted matrices.
# We test 3 dkfz cancer types: HGG_K27M, HGGother and MB_SHH. We also run
# a PAN cancer test including all three cancer types.

# Step 1: Generates mutation-sample matrix
Rscript ./ped_scripts/data_prep_ped.R ./ped_data/dkfz/dkfzGeneSample.txt dkfz

# Step 2: Generates 2500 permuted matrices (5 batches of 50 matrices)
for i in {1..5}
do
  Rscript ./ped_scripts/GI_detection_permutation.R HGG_K27M $i dkfz 50
  Rscript ./ped_scripts/GI_detection_permutation.R HGGother $i dkfz 50
  Rscript ./ped_scripts/GI_detection_permutation.R MB_SHH $i dkfz 50
done

# Step 2A PAN cancer test
Rscript ./ped_scripts/GI_detection_permutation_PAN_A.R dkfz
# Step 2B Pan cancer test
for i in {1..5}
do
  Rscript ./ped_scripts/GI_detection_permutation_PAN_B.R $i dkfz 50
done

# Step 3A: Get observed cooccurrence counts for each gene pair
Rscript ./ped_scripts/GI_detection_Empirical_P_A.R HGG_K27M dkfz 5
Rscript ./ped_scripts/GI_detection_Empirical_P_A.R HGGother dkfz 5
Rscript ./ped_scripts/GI_detection_Empirical_P_A.R MB_SHH dkfz 5
Rscript ./ped_scripts/GI_detection_Empirical_P_A.R PAN dkfz 5

# Step 3B: For each RDS file, get cooccurrence counts for each of its
# permuted matrices and compare with observed counts
for i in {1..5}
do
  Rscript ./ped_scripts/GI_detection_Empirical_P_B.R HGG_K27M dkfz $i
  Rscript ./ped_scripts/GI_detection_Empirical_P_B.R HGGother dkfz $i
  Rscript ./ped_scripts/GI_detection_Empirical_P_B.R MB_SHH dkfz $i
  Rscript ./ped_scripts/GI_detection_Empirical_P_B.R PAN dkfz $i
done

# Step 3C: Calculate empirical p-value by adding up all comparisons from step B
Rscript ./ped_scripts/GI_detection_Empirical_P_C.R HGG_K27M dkfz 2500
Rscript ./ped_scripts/GI_detection_Empirical_P_C.R HGGother dkfz 2500
Rscript ./ped_scripts/GI_detection_Empirical_P_C.R MB_SHH dkfz 2500
Rscript ./ped_scripts/GI_detection_Empirical_P_C.R PAN dkfz 2500

# Step 4A: Randomly pick a matrix, repeat 25 times
Rscript ./ped_scripts/null_dist_A.R HGG_K27M dkfz 25
Rscript ./ped_scripts/null_dist_A.R HGGother dkfz 25
Rscript ./ped_scripts/null_dist_A.R MB_SHH dkfz 25
Rscript ./ped_scripts/null_dist_A.R PAN dkfz 25

# Step 4B: Compare sampled mtxs with mtxs in RDS file
for i in {1..5}
do
Rscript ./ped_scripts/null_dist_B.R HGG_K27M $i dkfz
Rscript ./ped_scripts/null_dist_B.R HGGother $i dkfz
Rscript ./ped_scripts/null_dist_B.R MB_SHH $i dkfz
Rscript ./ped_scripts/null_dist_B.R PAN $i dkfz
done

# Step 4C: Merge all data and write p-values to file
Rscript ./ped_scripts/null_dist_C.R HGG_K27M dkfz 5 25 2500 1 25
Rscript ./ped_scripts/null_dist_C.R HGGother dkfz 5 25 2500 1 25
Rscript ./ped_scripts/null_dist_C.R MB_SHH dkfz 5 25 2500 1 25
Rscript ./ped_scripts/null_dist_C.R PAN dkfz 5 25 2500 1 25

# Step 5 Calculate the empirical FDR given the observed p-values
# and a null distribution of p-values
Rscript ./ped_scripts/FDR_correction.R HGG_K27M dkfz
Rscript ./ped_scripts/FDR_correction.R HGGother dkfz
Rscript ./ped_scripts/FDR_correction.R MB_SHH dkfz
Rscript ./ped_scripts/FDR_correction.R PAN dkfz
