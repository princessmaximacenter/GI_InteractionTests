#! /bin/bash

#input arguments: project, nr of repeats (10)
project=$1
start_repeats=$2
end_repeats=$3

# get current directory
wsdir=$(echo $PWD)

# use canctypes file in project folder
# all cancer types
declare -a cantypes
readarray -t cantypes < ${project}/cantypes.txt

for i in $(seq ${start_repeats} ${end_repeats})
do
  # create dir project_conf_i/data
  mkdir -p ${wsdir}/${project}_conf_${i}/data/
  # create dir project_conf_i/preproc/permuted_mut/
  mkdir -p ${wsdir}/${project}_conf_${i}/preproc/permuted_mut/
  # create dir project_conf_i/preproc/wrs/
  mkdir -p ${wsdir}/${project}_conf_${i}/preproc/wrs/
  for ct in ${cantypes[@]}
  do
    ln -s ${wsdir}/${project}_conf/data/${ct}_${i} ${wsdir}/${project}_conf_${i}/data/
    ln -s ${wsdir}/${project}_conf/preproc/permuted_mut/${ct}_${i} ${wsdir}/${project}_conf_${i}/preproc/permuted_mut/
    ln -s ${wsdir}/${project}_conf/preproc/wrs/${ct}_${i} ${wsdir}/${project}_conf_${i}/preproc/wrs/
    echo ${ct}_${i} >> ${project}_conf_${i}/cantypes.txt
  done
  mkdir -p ${wsdir}/${project}_conf_${i}/data/PAN
  cp ${project}/data/PAN/PAN_smut_list.txt ${project}_conf_${i}/data/PAN/PAN_smut_list.txt
done
