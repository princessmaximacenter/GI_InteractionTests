#! /bin/bash

set -eo pipefail

if [ -d "./wesme_venv" ]
then
  source wesme_venv/bin/activate
else
  echo "Error: Cannot load virtual environment. Directory wesme_venv does not exists."
fi

project=$1
ct=$2
ns=$3
come=${4:-come}
seed=${5:-100}
pths_def=(1.1)
pths=${6:-${pths_def[@]}}

if [[ "$come" == "come" ]]; then
  python comp_sample_weights.py ${ct} smut \
  -i ${project}/data/${ct}/${ct}_smut_list.txt \
  -o ${project}/data/${ct}/${ct}_smut_freq.txt

  python run_weighted_sampling.py ${ct} smut ${ns} \
  -m ${project}/data/${ct}/${ct}_smut_list.txt \
  -f ${project}/data/${ct}/${ct}_smut_freq.txt \
  -o ${project}/preproc/wrs/$ct/smut/$ns/ \
  -s ${seed}
fi

for p in ${pths[@]}
do
  if [[ "$come" =~ "me" ]]; then
    python comp_me_for_all_pairs.py ${ct} smut ${p} \
    -m ${project}/data/${ct}/${ct}_smut_list.txt \
    -s ${project}/preproc/wrs/$ct/smut/$ns/ \
    -o ${project}/results/${ct}_smut_me_pvs_${p}.txt
  fi

  if [[ "$come" =~ "co" ]]; then
    python comp_co_for_all_pairs.py ${ct} smut ${p} \
    -m ${project}/data/${ct}/${ct}_smut_list.txt \
    -s ${project}/preproc/wrs/$ct/smut/$ns/ \
    -o ${project}/results/${ct}_smut_co_pvs_${p}.txt
  fi
done

if [ -d "./wesme_venv" ]
then
  deactivate
fi
