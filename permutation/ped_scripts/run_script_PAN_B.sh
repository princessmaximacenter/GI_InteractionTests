#! /bin/bash
set -eo pipefail

project=$1
nperm=$2

if [ "${SLURM_ARRAY_TASK_ID}" = "undefined" ]; then
  seed=0
else
  seed=$((${SLURM_ARRAY_TASK_ID}))
fi

module load R/3.4.1
Rscript ./ped_scripts/GI_detection_permutation_PAN_B.R ${seed} ${project} ${nperm}
