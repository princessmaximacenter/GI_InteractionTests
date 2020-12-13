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
ct=$2
nmut=$3
come=${4:-come}
pths_def=(1.1)
pths=${5:-${pths_def[@]}}

for p in ${pths[@]}
do
  if [[ "$come" =~ "me" ]]; then
    python comp_fdr.py ${ct} smut me ${nmut} ${p} 2 \
    -m ${project}/data/${ct}/${ct}_smut_list.txt \
    -p ${project}/results/${ct}_smut_me_pvs_${p}.txt \
    -n ${project}/preproc/permuted_pv/${ct}/perm_${ct}_smut_me_pvs_${p}_ \
    -f ${project}/results/fdr/${ct}_smut_me_${p}_${nmut}_fdr_2bins_2.txt \
    -pf ${project}/results/fdr/${ct}_smut_me_${p}_${nmut}_pv_fdr_2bins_2.txt
  fi

  if [[ "$come" =~ "co" ]]; then
    python comp_fdr.py ${ct} smut co ${nmut} ${p} 2 \
    -m ${project}/data/${ct}/${ct}_smut_list.txt \
    -p ${project}/results/${ct}_smut_co_pvs_${p}.txt \
    -n ${project}/preproc/permuted_pv/${ct}/perm_${ct}_smut_co_pvs_${p}_ \
    -f ${project}/results/fdr/${ct}_smut_co_${p}_${nmut}_fdr_2bins_2.txt \
    -pf ${project}/results/fdr/${ct}_smut_co_${p}_${nmut}_pv_fdr_2bins_2.txt
  fi
done

if [ -d "./wesme_venv" ]
then
  deactivate
fi
