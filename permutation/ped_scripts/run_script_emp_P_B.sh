#! /bin/bash
set -eo pipefail

ct=$1
project=$2

if [ "${SLURM_ARRAY_TASK_ID}" = "undefined" ]; then
  seed=0
else
  seed=$((${SLURM_ARRAY_TASK_ID}))
fi

module load R/3.4.1
Rscript ./ped_scripts/GI_detection_Empirical_P_B.R ${ct} ${project} ${seed}
