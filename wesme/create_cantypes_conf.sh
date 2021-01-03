#! /bin/bash

#input arguments: project, nr of repeats (10)
project=$1
start_repeats=$2
end_repeats=$3

# use canctypes file in project folder
# all cancer types
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
