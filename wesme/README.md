# Run Permutation test, local and on hpc


Copy the **wesme** folder to your local machine or high performance cluster (hpc). The WeSME test is written in Python. To run the WeSME test, Python 2.7 and additional modules need to be installed. The file `requirements.txt` lists the modules and their versions:

```
decorator==4.1.2
networkx==1.11
numpy==1.16.0
pandas==0.21.0
python-dateutil==2.6.1
pytz==2017.3
scipy==0.17.1
six==1.11.0
```

You can use a virtual Python environment to make sure the correct versions are installed and used for this project without affecting other projects.


## Set up virtual Python environment
For more info, see: http://docs.python-guide.org/en/latest/dev/virtualenvs/

#### Install virtualenv via pip

`pip install virtualenv`

#### Test your installation

`virtualenv --version`

#### Create a virtual environment for the WeSME project

```bash
cd wesme
virtualenv -p /usr/bin/python2.7 wesme_venv
````

To begin using the virtual environment, it needs to be activated

`source wesme_venv/bin/activate`

#### Install the required packages

`pip install -r requirements.txt`

if you get an error message, try:

`python -m pip install -r requirements.txt`

If you are done working in the virtual environment, you can deactivate it:

`deactivate`

Now that you have created the virtual environment 'wesme_venv' you can use it with:

`source wesme_venv/bin/activate`


## Set up a project

Add a project folder `<projectname>` to the **wesme** directory. Create a file `cantypes.txt` in this folder which contains all the cancer types of your project. In the project directory, create a `data` folder which contains a directory per cancer type. Each cancer type folder contains the file `<cancer type>_smut_list.txt`. This file contains the *gene sample matrix* in a specific format, as is shown in the example below. **Note that samples are 0-based**.

Example:

**Gene-sample matrix**

|        |Smpl1|Smpl2|Smpl3|Smpl4|
|:-------|:----|:----|:----|:----|
| GeneA  | 0   | 1   | 1   | 0   |
| GeneB  | 0   | 0   | 0   | 1   |
| GeneC  | 0   | 1   | 0   | 0   |
| GeneD  | 1   | 0   | 1   | 1   |
| GeneE  | 1   | 0   | 0   | 0   |


**`<ct>_smut_list.txt`**

```none
samples -> Smpl1,Smpl2,Smpl3,Smpl4
GeneA -> 1,2
GeneB -> 3
GeneC -> 1
GeneD -> 0,2,3
GeneE -> 0	
```
(-> stands for tab)

## Run the WeSME test

The WeSME test runs in 3 steps:
- **Step 1** calculates the p-values of each gene pair. 
- **Step 2** will create and test 300 random matrices to create a p-value null distribution. This will be done with a job array; you submit one job and this job will start 300 subjobs. 
- In **step 3** the p-value null distribution is used to calculate an empirical FDR.  See example below for project dkfz. For each cancer type, the previous step needs to be completed before you can run the next step.

Below, we will give examples of each steps for the dkfz test project. Run the code from the **wesme** folder.

### STEP 1 - Calculate the p-values for each gene pair 

**local**:

`./run_step1_can.sh <project> <ct> <n_sampl> <come=come> <seed=100> <pval_thr=1.1>`

example:

`./run_step1_can.sh dkfz HGG_K27M 10000`


**on hpc:**

```none
batch_run_step1_slurm.sh -p project -s nsampl -m vmem -c cantype
                         -t hrs -v pval -o come -e seed

parameters:
  -p : name of project folder
  -s : number of resamplings
  -m : memory requested (default: 10G)
  -c : cancer type (default: all in cantypes.txt file)
  -t : time requested in hours (default: 3)
  -v : p value threshold (default: 1.1)
  -o : come setting (default come, alternatives: co or me)
  -e : seed for weighted sampling
```

example:

`./batch_run_step1_slurm.sh -p dkfz -s 10000`


You can check the status of the job in *SLURM* with `squeue -u <username>`. After the job is finished, you can find a log file with information about the job in the **log** directory.
 
### STEP 2 - Create and test random matrices to create a p-value null distribution

**local**:

`./run_step2_can_i_slurm.sh <project> <ct> <n_sampl> <come=come> <nmut=1> <pval_thr=1.1>`

example:

```bash
for i in {1..5}
do
  `./run_step2_can_i_slurm.sh dkfz HGG_K27M 10000 come $i`
done
```


**on hpc:**

```none
batch_run_step2_slurm.sh -p project -n npermut -s nsampl -m vmem -c cantype
                         -t hrs -v pval -o come -f first -l last

parameters:
  -p : name of project folder
  -n : number of permutations for FDR calculation
  -s : number of resamplings
  -m : memory requested (default: 10G)
  -c : cancer type (default: all in cantypes.txt file)
  -t : time requested in hours (default: 8)
  -v : p value threshold (default: 1.1)
  -o : come setting (default come, alternatives: co or me)
  -f : first job nb
  -l : last job nb
```

example:

`./batch_run_step2_slurm.sh -p dkfz -n 300 -s 10000`


### STEP 3 Calculate an empirical FDR

**local:**

`./run_step3_can.sh <project> <ct> <nmut> <come=come> <pth=1.1> `

example:

`./run_step3_can.sh dkfz HGG_K27M 5`

**on hpc:**

```none
batch_run_step3_slurm.sh -p project -n npermut -m vmem -c cantype
                         -t hrs -v pval -o come

parameters:
  -p : name of project folder
  -n : number of permutations
  -m : memory requested (default: 10G)
  -c : cancer type (default: all in cantypes.txt file)
  -t : time requested in hours (default: 12)
  -v : p value threshold (default: 1.1)
  -o : come setting (default come, alternatives: co or me)
```

example:

`./batch_run_step3_slurm.sh -p dkfz -n 300`


### Run the PAN cancer test

Next, you can run the PAN cancer test. To reduce memory and time requirements, run the co-occurrence and mutual exclusivity test separately using the parameter -o.

**on hpc:**

```none
# PAN cancer step 1
./batch_run_step1_PAN_slurm.sh -p dkfz -s 10000 -o co
./batch_run_step1_PAN_slurm.sh -p dkfz -s 10000 -o me

# PAN cancer step 2. Ask for at least 10G memory, to avoid memory errors
./batch_run_step2_PAN_slurm.sh -p dkfz -n 300 -s 10000 -m 10G -o co 
./batch_run_step2_PAN_slurm.sh -p dkfz -n 300 -s 10000 -m 10G -o me 

# PAN cancer step 3
./batch_run_step3_slurm.sh -p dkfz -n 300 -m 50G -c PAN -o co
./batch_run_step3_slurm.sh -p dkfz -n 300 -m 50G -c PAN -o me
```


### Clean up big intermediate files

After step 3 is successfully finished, remove the (big) `permuted_pv` directory:

`rm -rf ./<projectname>/preproc/permuted_pv`

Use a table such as the example below to note which steps are running/finished and how much memory was needed (you can find it back in the log directory)

#### Schedule WeSME test


| cancer type | st1   |time|vmem| st2   |time|vmem| st3   |time|vmem|clean pv 
|:------------|:------|:---|:---|:------|:---|:---|:------|:---|:---|:-------
| HGG_K27M    |       |    |    |       |    |    |       |    |    |
| HGG_other   |       |    |    |       |    |    |       |    |    |
| MB_SHH      |       |    |    |       |    |    |       |    |    |














