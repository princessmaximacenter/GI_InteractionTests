#! /bin/bash

set -eo pipefail

if [ -d "./wesme_venv" ]
then
  source wesme_venv/bin/activate
else
  echo "Error: Cannot load virtual environment. Directory wesme_venv does not exists."
  echo "Running without virtual python environment"
fi

project=$1
ns=$2
come=${3:come}
p=1.1

if [ "$SLURM_ARRAY_TASK_ID" = "undefined" ] || [ -z "$SLURM_ARRAY_TASK_ID" ]; then
  taskid=${4:-1}
  nmut=$(($taskid - 1))
else
  nmut=$(($SLURM_ARRAY_TASK_ID - 1))
fi

if [[ "$come" =~ "me" ]]; then
  # compute ME for permuted data
  python comp_come_for_all_pairs_PAN.py ${ns} smut ${p} ${project} me \
    -s ${project}/preproc/wrs/ \
    -m ${project}/preproc/permuted_mut/PAN/smut/PAN_smut_permuted_cover_${nmut}.txt \
    -o ${project}/preproc/permuted_pv/PAN/perm_PAN_smut_me_pvs_${p}_${nmut}.txt
fi

if [[ "$come" =~ "co" ]]; then
  # compute CO for permuted data
  python comp_come_for_all_pairs_PAN.py ${ns} smut ${p} ${project} co \
    -s ${project}/preproc/wrs/ \
    -m ${project}/preproc/permuted_mut/PAN/smut/PAN_smut_permuted_cover_${nmut}.txt \
    -o ${project}/preproc/permuted_pv/PAN/perm_PAN_smut_co_pvs_${p}_${nmut}.txt
fi

if [ -d "./wesme_venv" ]
then
  deactivate
fi
