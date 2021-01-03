#! /bin/bash

#input arguments: project, nr of repeats (10)
project=$1
start_repeats=$2
end_repeats=$3

# use canctypes file in project folder
# all cancer types
declare -a cantypes
readarray -t cantypes < ${project}/cantypes.txt

rm -f ${project}_conf/cantypes.txt
for ct in ${cantypes[@]}
do
  for i in $(seq ${start_repeats} ${end_repeats})
  do
    # remove subfolder
    rm -rf ${project}_conf/${ct}_$i
    # create data dir and subfolder and copy smut file from original dir
    mkdir -p ${project}_conf/data/${ct}_$i
    cp ${project}/data/${ct}/${ct}_smut_list.txt ${project}_conf/data/${ct}_$i/${ct}_${i}_smut_list.txt
    # add to cantypes file containing BALL_1, BALL_2, etc
    echo ${ct}_$i >> ${project}_conf/cantypes.txt
  done
done

# do the same for PAN
for i in $(seq ${start_repeats} ${end_repeats})
do
  # remove subfolder
  rm -rf ${project}_conf/PAN_$i
  # create data dir and subfolder and copy smut file from original dir
  mkdir -p ${project}_conf/data/PAN_$i
  cp ${project}/data/PAN/PAN_smut_list.txt ${project}_conf/data/PAN_$i/PAN_${i}_smut_list.txt
  # add to cantypes file containing BALL_1, BALL_2, etc
done

# copy cantypes.txt to cantypes_ct.txt to have an unchanged file with all
# cancer types sorted on cancer type.

cp ./${project}_conf/cantypes.txt ./ ${project}_conf/cantypes_ct.txt

# use canctypes file in project folder
# to create file with all cancer types
declare -a cantypes
readarray -t cantypes < ${project}/cantypes.txt

rm -f ${project}_conf/cantypes_all.txt
for i in $(seq ${start_repeats} ${end_repeats})
do
  for ct in ${cantypes[@]}
  do
    echo ${ct}_${i} >> ${project}_conf/cantypes_all.txt
  done
  echo PAN_${i} >> ${project}_conf/cantypes_all.txt
done
