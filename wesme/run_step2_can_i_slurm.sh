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
ns=$3
come=${4:-come}

if [ "$SLURM_ARRAY_TASK_ID" = "undefined" ] || [ -z "$SLURM_ARRAY_TASK_ID" ]; then
  taskid=${5:-1}
  nmut=$(($taskid - 1))
else
  nmut=$(($SLURM_ARRAY_TASK_ID - 1))
fi

pths_def=(1.1)
pths=${6:-${pths_def[@]}}

# permute data matrices
if [[ "$come" == "come" ]]; then
  python run_permute_data.py ${ct} smut ${nmut} 100 \
    -i ${project}/data/${ct}/${ct}_smut_list.txt \
    -p ${project}/preproc/permuted_mut/${ct}/smut/${ct}_smut_permuted_cover
fi

for p in ${pths[@]}
do
  if [[ "$come" =~ "me" ]]; then
    # compute ME for permuted data
    python comp_me_for_all_pairs.py ${ct} smut ${p} \
      -s ${project}/preproc/wrs/$ct/smut/$ns/ \
      -m ${project}/preproc/permuted_mut/${ct}/smut/${ct}_smut_permuted_cover_${nmut}.txt \
      -o ${project}/preproc/permuted_pv/${ct}/perm_${ct}_smut_me_pvs_${p}_${nmut}.txt
  fi

  if [[ "$come" =~ "co" ]]; then
    # compute CO for permuted data
    python comp_co_for_all_pairs.py ${ct} smut ${p} \
      -s ${project}/preproc/wrs/$ct/smut/$ns/ \
      -m ${project}/preproc/permuted_mut/${ct}/smut/${ct}_smut_permuted_cover_${nmut}.txt \
      -o ${project}/preproc/permuted_pv/${ct}/perm_${ct}_smut_co_pvs_${p}_${nmut}.txt
  fi
done

if [ -d "./wesme_venv" ]
then
  deactivate
fi
