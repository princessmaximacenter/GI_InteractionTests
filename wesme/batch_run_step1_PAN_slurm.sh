#! /bin/bash

set -eo pipefail

showUsage ()
{
  echo "usage:
  batch_run_step1_PAN_slurm.sh -p project -s nsampl -m vmem -t hrs -o come

  parameters:
    -p : name of project folder
    -s : number of resamplings
    -m : memory requested (default: 24G)
    -t : time requested in hours (default: 4)
    -o : come setting (default come, alternatives: co or me)
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
    -t)
        hr="$argVal"
    ;;
    -o)
        come="$argVal"
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

# set optional values
: ${vmem:="24G"}
: ${hr:="4"}
: ${come:="come"}

# FINISHED PARSING ARGUMENTS
#--------------------------------------------------------------------------

# START OF SCRIPT

# create name for log file folder
d=`date +%Y%m%d_%H%M%S`

# create log directories
logdir=./log/${project}_st1_PAN_${d}
mkdir -p ${logdir}

proj_name=${project//[\/]/_}

hold_jid=$(sbatch --time=${hr}:00:00 --mem=${vmem} --parsable \
-e ${logdir}/st1_${proj_name}_PAN_${come}.e \
-o /dev/null \
-J PAN_${proj_name}_${come}_st1_ws \
./run_step1_PAN.sh ${project} ${nsampl} ${come})

sbatch -d $hold_jid -J PAN_${proj_name}_${come}_st1_log \
-o ${logdir}/st1_${proj_name}_PAN_${come}_${hold_jid}.o \
-e /dev/null \
log_job_slurm.sh ${proj_name} PAN ${nsampl} 1.1 ${vmem} ${hr} ${come} \
${hold_jid} undefined undefined ${logdir} 1
