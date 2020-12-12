# Run Permutation test, local and on hpc


The permutation folder contains 4 subfolders:
- ped_data
- ped_log (only used with cluster scripts)
- ped_results
- ped_scripts


The permutation test is written in R. To run the test, several additional R packages need to be installed:

- data.table
- reshape2
- vegan
- gtools


We have supplied **example bash scripts** to run the permutation test on a **hpc cluster** with a *SLURM* job submission framework and the *Lmod* environment module system. These scripts need to be adapted to your own hpc set up.

Below we describe the pipeline, which consists of step 0 to step 5. In our project, the permutation test was run with a null distribution of **1M permuted matrices**, generated and processed in parallel in **200 batches** of **5,000 matrices**. The FDR estimation is done on **100 randomly sampled permuted matrices**.

In the hpc examples we use the same numbers, but the **local machine examples** -which directly call the R scripts and are meant for getting familiar with the pipeline- assume a null distribution of **250 permuted matrices**, generated and processed in **5 batches** of **50 matrices**. The FDR estimation is then done on **25 randomly sampled permuted matrices**. These values are also the default settings when calling the R scripts directly. See example_script_local.sh for a full example how to run the Permutation test with 3 cancer types and how to run a PAN cancer test.

The permutation test consists of several steps:

## Step 0 - Prepare files and directories

- Create a gene sample file (e.g.: `dkfzGeneSample.txt`) with fields: `sample`, `ct`, `ct.sub`, `gene` and put it in the project directory. De `ct` and `ct.sub` fields contain the cancer type and cancer subtype. The permutation test will in principle run the test per cancer type, but you can opt in step 1 for grouping on subtype.

- Create a project directory in `ped_data` (e.g. *dkfz*) and put the gene sample file in this directory.


 
## Step 1 - Generate mutation-sample matrix


From the `permutation` folder run:

**local:**

`Rscript ./ped_scripts/data_prep_ped.R <genesamplefile> <project> <use_sub=F> <has_header=T>`

example:

`Rscript ./ped_scripts/data_prep_ped.R ./ped_data/dkfz/dkfzGeneSample.txt dkfz`

**on hpc:**

`run_script_data_prep.sh -p project -f inputfile -l logdir -s use_sub -r has_header`

``` none
parameters:
  -p : name of project folder
  -f : path to inputfile
  -l : log dir (default: ./ped_log)
  -s : use cancer subtype (TRUE/FALSE, default: FALSE)
  -r : file has header (TRUE/FALSE, default: TRUE)
```

example:

`./ped_scripts/run_script_data_prep.sh -p dkfz -f ./ped_data/dkfz/dkfzGeneSample.txt`


This script will create the *<project>* folder in `ped_data` (if not already existing) and in `ped_results`, it will create the file `matrix_gene_sample_somatic.txt` in the *<project>* folder in `ped_data`, and the file `<project>_CTs_list.txt` in the `ped_results` folder.


## Step 2 - Generate permuted matrices


**local:**

`Rscript ./ped_scripts/GI_detection_permutation.R <ct> <seed> <project> <n_perm>`

for example, to create batch 1-5 with each 50 permuted matrices for HGG_K27M cancer type:

``` bash
for i in {1..5}
do
  Rscript ./ped_scripts/GI_detection_permutation.R HGG_K27M $i dkfz 50
done
```

To run a PAN cancer test, repeat the above script for the two other cancer types in the test set: *HGGother* and *MB_SHH*. 
 
**on hpc:**

```none
./ped_scripts/batch_run_script.sh -p project -m vmem -c cantype -t hrs -s startseed \
-o stopseed -n nperm -l logdir`

Parameters:
    -p : name of project folder
    -m : memory requested (default: 30G)
    -c : cancer type (default: all in <projectname>_CTs_list.txt file)
    -t : time requested in hours (default: 10)
    -s : start seed (default: 1)
    -o : stop seed (default: 200)
    -n : number of permutations (default: 5000)
    -l : logdir
```

for example, to create 1M permuted matrices in 200 batches of 5000 matrices for all cancer types:

`./ped_scripts/batch_run_script.sh -p dkfz -s 1 -o 200 -n 5000`

The **PAN cancer test** runs in 2 steps:

**local:**

```none
Rscript ./ped_scripts/GI_detection_permutation_PAN_A.R <project>
Rscript ./ped_scripts/GI_detection_permutation_PAN_B.R <seed> <project> <n_perm>
```

for example:

```bash
Rscript ./ped_scripts/GI_detection_permutation_PAN_A.R dkfz
for i in {1..5}
do
  Rscript ./ped_scripts/GI_detection_permutation_PAN_B.R $i dkfz 50
done
```

**on hpc:**

```none
./ped_scripts/run_script_PAN_A.sh -p dkfz
./ped_scripts/batch_run_script_PAN_B.sh -p dkfz -s 1 -o 200 -n 5000
```

For the PAN cancer analysis, run all following steps, but choose *PAN* as cancer type
 
## Step 3 - Calculate p-values 


Empirical p-values are calculated in 3 steps:

### Step A - Get observed cooccurrence counts for each gene pair

**local:**

`Rscript ./ped_scripts/GI_detection_Empirical_P_A.R <ct> <project> <n_RDS=5>`

example:

`Rscript ./ped_scripts/GI_detection_Empirical_P_A.R HGG_K27M dkfz`

**on hpc:**

```none
run_script_emp_P_A.sh -p project -m vmem -c cantype -t hrs -n nb_RDS -l logdir

parameters:
  -p : name of project folder
  -m : memory requested (default: 20G)
  -c : cancer type (default: all in <projectname>_CTs_list.txt file)
  -t : time requested in hours (default: 10)
  -n : number of RDS files (default: 200)
  -l : log dir
```

example:

`./ped_scripts/run_script_emp_P_A.sh -p dkfz -n 200`

### Step B - For each RDS file, get cooccurrence counts for each of its permuted matrices and compare with observed counts

**local:**

`Rscript ./ped_scripts/GI_detection_Empirical_P_B.R <ct> <project> <n_rds>`

example:

```bash
for i in {1..5}
do
  Rscript ./ped_scripts/GI_detection_Empirical_P_B.R HGG_K27M dkfz $i
done
```

**on hpc:**

```none
batch_run_script_emp_P_B.sh -p project -m vmem -c cantype -t hrs -s startseed \
-o stopseed -l logdir

parameters:
  -p : name of project folder
  -m : memory requested (default: 30G)
  -c : cancer type (default: all in <projectname>_CTs_list.txt file)
  -t : time requested in hours (default: 10)
  -s : start seed (default: 1)
  -o : stop seed (default: 200)
  -l : log dir
```

example:

`./ped_scripts/batch_run_script_emp_P_B.sh -p dkfz -s 1 -o 200`


### Step C - Calculate empirical p-value by adding up all comparisons from step B

**local:**

`Rscript ./ped_scripts/GI_detection_Empirical_P_C.R <ct> <project> <n_perm=250>`

example:

`Rscript ./ped_scripts/GI_detection_Empirical_P_C.R HGG_K27M dkfz 250`

**on hpc:**

```none
run_script_emp_P_C.sh -p project -m vmem -c cantype -t hrs -n nperm -l logdir

parameters:
  -p : name of project folder
  -m : memory requested (default: 20G)
  -c : cancer type (default: all in <projectname>_CTs_list.txt file)
  -t : time requested in hours (default: 4)
  -n : number of permutations (default: 1e+06)
  -l : log dir
```

example:

`./ped_scripts/run_script_emp_P_C.sh -p dkfz`
 

## Step 4 – Generate a p-value null distribution for the FDR calculation

The p-value null distribution is generated in 3 steps:

### Step A: Randomly sample N permuted matrices

**local:**

`Rscript ./ped_scripts/null_dist_A.R <ct> <project> <nperm=25>`


example:

`Rscript ./ped_scripts/null_dist_A.R HGG_K27M dkfz 25`

**on hpc:**

```none
run_script_null_dist_A.sh -p project -m vmem -c cantype -t hrs -n n_mtx -l logdir

parameters:
  -p : name of project folder
  -m : memory requested (default: 30G)
  -c : cancer type (default: all in <projectname>_CTs_list.txt file)
  -t : time requested in hours (default: 8)
  -n : number of random sampled matrices (default: 100)
  -l : log dir
```

example:

`run_script_null_dist_A.sh -p dkfz`


### Step B: Compare sampled matrices with matrices in RDS file

**local:**

`Rscript ./ped_scripts/null_dist_B.R <ct> <seed> <project>`

example:

```bash
for i in {1..5}
do
  Rscript ./ped_scripts/null_dist_B.R HGG_K27M $i dkfz
done
```

**on hpc:**

```none
batch_run_script_null_dist_B.sh -p project -m vmem -c cantype -t hrs \
-s startseed -o stopseed -l logdir

parameters:
  -p : name of project folder
  -m : memory requested (default: 30G)
  -c : cancer type (default: all in <projectname>_CTs_list.txt file)
  -t : time requested in hours (default: 10)
  -s : start seed (default: 1)
  -o : stop seed (default: 200)
  -l : log dir
```


example:

`./ped_scripts/batch_run_script_null_dist_B.sh -p dkfz -s 1 -o 200`


### Step C: Merge all data and write p-values to file

**local:**

```none
Rscript ./ped_scripts/null_dist_C.R <ct> <project> <n_RDS> <n_sampled_mtxs> \
<n_perm_mtxs> <start_perm_mtxs> <end_perm_mtxs>
```

example:

`Rscript ./ped_scripts/null_dist_C.R HGG_K27M dkfz 5 25 250 1 25`

**on hpc:**

```bash
run_script_null_dist_C.sh -p project -m vmem -c cantype -t hrs -r n_rds \
                          -s n_smpl -n n_perm -l logdir -f firstmtx -e endmtx

parameters:
  -p : name of project folder
  -m : memory requested (default: 12G)
  -c : cancer type (default: all in <projectname>_CTs_list.txt file)
  -t : time requested in hours (default: 4)
  -r : number of RDS files (default: 200)
  -s : number of sampled matrices (default: 100)
  -n : number of permuted matrices (default: 1e+06)
  -l : log dir
  -f : number of first random mtx to process (default: 1)
  -e : number of last random mtx to process (default: s)
```

example:

`run_script_null_dist_C.sh -p dkfz`


## Step 5 – Estimate FDR 

**local:**

`Rscript ./ped_scripts/FDR_correction.R <ct> <project>`

example:

`Rscript ./ped_scripts/FDR_correction.R HGG_K27M dkfz`

**on hpc:**

```none
run_script_FDR_correction.sh -p project -m vmem -c cantype -t hrs -l logdir

parameters:
  -p : name of project folder
  -m : memory requested (default: 8G)
  -c : cancer type (default: all in <projectname>_CTs_list.txt file)
  -t : time requested in hours (default: 1)
  -l : log dir
```

example:

`./ped_scripts/run_script_FDR_correction.sh -p dkfz`


