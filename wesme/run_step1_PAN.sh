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
come=$3
p=1.1

if [[ "$come" =~ "me" ]]; then
  python comp_come_for_all_pairs_PAN.py ${ns} smut ${p} ${project} me \
  -m ${project}/data/PAN/PAN_smut_list.txt \
  -s ${project}/preproc/wrs/ \
  -o ${project}/results/PAN_smut_me_pvs_${p}.txt
fi

if [[ "$come" =~ "co" ]]; then
  python comp_come_for_all_pairs_PAN.py ${ns} smut ${p} ${project} co \
  -m ${project}/data/PAN/PAN_smut_list.txt \
  -s ${project}/preproc/wrs/ \
  -o ${project}/results/PAN_smut_co_pvs_${p}.txt
fi

if [ -d "./wesme_venv" ]
then
  deactivate
fi
