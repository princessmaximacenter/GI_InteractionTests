#!/bin/bash -l
#$ vmem="10M"
#$ hr:="3"

#set -eo pipefail

function seconds2time ()
{
   T=$1
   D=$((T/60/60/24))
   H=$((T/60/60%24))
   M=$((T/60%60))
   S=$((T%60))

   if [[ ${D} != 0 ]]
   then
      printf '%d days %02d:%02d:%02d' $D $H $M $S
   else
      printf '%02d:%02d:%02d' $H $M $S
   fi
}

function Kb2Gb ()
{
  x=$1
  k=${x::-1}
  b=$((k * 1024))
  d=''; s=0; S=(Bytes {K,M,G,T,E,P,Y,Z}B)
  while ((b > 1024)); do
    d="$(printf ".%02d" $((b % 1024 * 100 / 1024)))"
    b=$((b / 1024))
    let s++
  done
  echo "$b$d ${S[$s]}"
}

project=$1
ct=$2
ns=$3
pths_def=(1.1)
pths=${4:-${pths_def[@]}}
vmem=$5
hrs=$6
come=${7:come}
jobid=$8
firsttask=$9
lasttask=${10}
logdir=${11}
step=${12}

# wait to make sure all info of last job is there
sleep 10

startdate=$(sacct -j ${jobid} --format=start | awk 'NR==3')

echo Started WeSME step ${step} with jobid $jobid
echo on ${startdate}
echo with the following parameters:
echo project: ${project}
echo cancertype: ${ct}
echo nb of resamplings/permutations: ${ns}
echo pvalue thresholds: ${pths}
echo vmem: ${vmem}
echo hrs: ${hrs}
echo come: ${come}
echo first task: ${firsttask}
echo last task: ${lasttask}
echo

errorcount=$(grep -o -i "error" ${logdir}/*${project}_${ct}_${come}*.e | wc -l)
failedstatus=$(sacct -j ${jobid} --format=State | awk '{gsub ("State|-|COMPLETED", ""); print $1}' |  awk '/./') | wc -l
exitstatus=$(sacct -j ${jobid}  --format=ExitCode | awk '{gsub ("ExitCode|-|0:0", ""); print $1}' |  awk '/./') | uniq
if [ ${errorcount} -ge 1  ] || [ ${failedstatus} -ge 1 ] || [ ${exitstatus} -ge 1 ]
then
  echo job finished with ERRORS "(error count: $errorcount, failed: $failedstatus, exit>0: $exitstatus)"
else
  echo job finished SUCCESSFULLY
fi

#runtime=$(qacct -j ${jobid} | grep -E -A 50 "taskid[[:space:]]+${lasttask}" | grep -o -P "(?<=ru_wallclock)[[:space:]]*\d+")
#runtime_readable=$(seconds2time ${runtime})
#echo "(last sub-)job finished with runtime:" ${runtime_readable}
#maxvmem=$(qacct -j ${jobid} | grep -E -A 50 "taskid[[:space:]]+${lasttask}" | grep -o -P "(?<=maxvmem).+")
#echo "(last sub-)job needed:" ${maxvmem}

runtimes=($(sacct -j ${jobid} --format=elapsed | awk '{gsub ("Elapsed|-", ""); print $1}' |  awk '/./'))
IFS=$'\n'
maxruntime=$(echo "${runtimes[*]}" | sort -nr | head -n1)
echo "max runtime:" ${maxruntime}

IFS=" "
maxvmems=($(sacct -j ${jobid} --format=MaxRSS | awk '{gsub ("MaxRSS|-", "", $1); print $1}' |  awk '/./'))
IFS=$'\n'
maxvmem=$(echo "${maxvmems[*]}" | sort -nr | head -n1)
maxvmem_readable=$(Kb2Gb ${maxvmem})
echo "max vmem:" ${maxvmem_readable}
