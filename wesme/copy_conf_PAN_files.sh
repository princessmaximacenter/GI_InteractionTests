#! /bin/bash

#input arguments: project, nr of repeats (usually 9)
project=$1
start_repeats=$2
end_repeats=$3

for i in $(seq ${start_repeats} ${end_repeats})
do
  mv -f ./${project}_conf_${i}/results/fdr/PAN_smut_me_1.1_300_pv_fdr_2bins_2.txt \
  ./${project}_conf/results/fdr/PAN_${i}_smut_me_1.1_300_pv_fdr_2bins_2.txt
  mv -f ./${project}_conf_${i}/results/fdr/PAN_smut_co_1.1_300_pv_fdr_2bins_2.txt \
  ./${project}_conf/results/fdr/PAN_${i}_smut_co_1.1_300_pv_fdr_2bins_2.txt

  mv -f ./${project}_conf_${i}/results/fdr/PAN_smut_me_1.1_300_fdr_2bins_2.txt \
  ./${project}_conf/results/fdr/PAN_${i}_smut_me_1.1_300_fdr_2bins_2.txt
  mv -f ./${project}_conf_${i}/results/fdr/PAN_smut_co_1.1_300_fdr_2bins_2.txt \
  ./${project}_conf/results/fdr/PAN_${i}_smut_co_1.1_300_fdr_2bins_2.txt

  mv -f ./${project}_conf_${i}/results/PAN_smut_me_pvs_1.1.txt \
  ./${project}_conf/results//PAN_${i}_smut_me_pvs_1.1.txt
  mv -f ./${project}_conf_${i}/results/PAN_smut_co_pvs_1.1.txt \
  ./${project}_conf/results//PAN_${i}_smut_co_pvs_1.1.txt

  #rm -r ./${project}_conf_${i}/preproc/permuted_pv/PAN
  # remove <project>_conf_xx directories
  rm -r ./${project}_conf_${i}
done
