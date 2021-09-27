#!/bin/bash -l

set -eo pipefail

showUsage ()
{
  echo "usage:
  batch_run_step1_conf_slurm.sh -p project -s nsampl -m vmem -c cantype
                      -t hrs -v pval -o come -f firstseed -l lastseed

  parameters:
    -p : name of project folder
    -s : number of resamplings
    -m : memory requested (default: 10G)
    -c : cancer type (default: all in cantypes.txt file)
    -t : time requested in hours (default: 3)
    -v : p value threshold (default: 1.1)
    -o : come setting (default come, alternatives: co or me)
    -f : first seed for weighted sampling
    -l : last seed for weighted sampling
  "
}

#--------------------------------------------------------------------------
# PARSE THE arguments

COUNTER=0
ARGS=("$@")
while [ $COUNTER -lt $# ]
do
  arg=${ARGS[$COUNTER]}
  let COUNTER=COUNTER+1
  nextArg=${ARGS[$COUNTER]}

  if [[ $skipNext -eq 1 ]]; then
    #echo "Skipping"
    skipNext=0
    continue
  fi

  argKey=""
  argVal=""
  if [[ "$arg" =~ ^\- ]]; then
    # if the format is: -key=value
    if [[ "$arg" =~ \= ]]; then
      argVal=$(echo "$arg" | cut -d'=' -f2)
      argKey=$(echo "$arg" | cut -d'=' -f1)
      skipNext=0

      # if the format is: -key value
      elif [[ ! "$nextArg" =~ ^\- ]]; then
        argKey="$arg"
        argVal="$nextArg"
        skipNext=1

      # if the format is: -key (a boolean flag)
      elif [[ "$nextArg" =~ ^\- ]] || [[ -z "$nextArg" ]]; then
        argKey="$arg"
        argVal=""
        skipNext=0
      fi
  # if the format has not flag, just a value.
  else
    argKey=""
    argVal="$arg"
    skipNext=0
  fi

  case "$argKey" in
    -p)
        project="$argVal"
    ;;
    -s)
        nsampl="$argVal"
    ;;
    -m)
        vmem="$argVal"
    ;;
    -c)
        cantypes="$argVal"
    ;;
    -t)
        hr="$argVal"
    ;;
    -v)
        pth="$argVal"
    ;;
    -o)
        come="$argVal"
    ;;
    -f)
        firstseed="$argVal"
    ;;
    -l)
        lastseed="$argVal"
    ;;
    -h|--help|-help|--h)
        showUsage
        exit
    ;;
  esac
done

# check required parameters
if [[ -z "$project" ]];
then
  echo "ERROR: Missing parameter -p project folder"
  showUsage
  exit
fi

if [[ -z "$nsampl" ]];
then
  echo "ERROR: Missing parameter -s number of resamples"
  showUsage
  exit
fi

if [[ -z "$cantypes" ]];
then
  # check if cantypes.txt file exists in project folder
  # if not return error
  if [ -e ${project}/cantypes.txt ];
  then
    # else read all cancer types
    declare -a cantypes_def
    readarray -t cantypes_def < ${project}/cantypes.txt
    : ${cantypes:=${cantypes_def[@]}}
  else
    echo "ERROR: Missing both parameter -c cancertype and cantypes.txt file"
    showUsage
    exit
  fi
else
  ct=${cantypes}
  unset cantypes
  for i in $(seq ${firstseed} 1 ${lastseed}); do cantypes+=(${ct}_$i); done
fi

# set optional values
: ${vmem:="10G"}
: ${hr:="3"}
: ${pth:="1.1"}
: ${come:="come"}
: ${firstseed:="100"}
: ${lastseed:="100"}

# FINISHED PARSING ARGUMENTS
#--------------------------------------------------------------------------

# START OF SCRIPT

# create name for log file folder
d=`date +%Y%m%d_%H%M%S`

# create log directories
logdir=./log/${project}_st1_${d}
mkdir ${logdir}

for ct in ${cantypes[@]}
do
  IFS='_' read -r -a array <<< $ct
  seed=${array[-1]}
  hold_jid=$(sbatch --time=${hr}:00:00 --mem=${vmem} --parsable \
  -e ${logdir}/st1_${project}_${ct}_${come}.e \
  -o /dev/null \
  -J ${ct}_${project}_${come}_st1_ws \
  ./run_step1_can.sh ${project} ${ct} ${nsampl} ${pth} ${vmem} ${hr} ${come} ${seed})

  sbatch -d $hold_jid -J ${ct}_${project}_${come}_st1_log \
  -o ${logdir}/st1_${project}_${ct}_${come}_${hold_jid}.o \
  -e /dev/null \
  log_job.sh ${project} ${ct} ${nsampl} ${pth} ${vmem} ${hr} ${come} ${hold_jid} \
  undefined undefined ${logdir} 1
done
