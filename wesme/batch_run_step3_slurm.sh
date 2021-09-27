#!/bin/bash -l

set -eo pipefail

showUsage ()
{
  echo "usage:
  batch_run_step3_slurm.sh -p project -n npermut -m vmem -c cantype
                           -t hrs -v pval -o come

  parameters:
    -p : name of project folder
    -n : number of permutations
    -m : memory requested (default: 10G)
    -c : cancer type (default: all in cantypes.txt file)
    -t : time requested in hours (default: 12)
    -v : p value threshold (default: 1.1)
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
    -n)
        nmut="$argVal"
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

if [[ -z "$nmut" ]];
then
  echo "ERROR: Missing parameter -n number of permutations for fdr computation"
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
fi

# set optional values
: ${vmem:="50G"}
: ${hr:="12"}
: ${pth:="1.1"}
: ${come:="come"}

# FINISHED PARSING ARGUMENTS
#--------------------------------------------------------------------------

# START OF SCRIPT

# create name for log file folder
d=`date +%Y%m%d_%H%M%S`

# create log directories
logdir=./log/${project}_st3_${d}
mkdir -p ${logdir}

proj_name=${project//[\/]/_}

for ct in ${cantypes[@]}
do
  hold_jid=$(sbatch --time=${hr}:00:00 --mem=${vmem} --parsable \
  -e ${logdir}/st3_${proj_name}_${ct}_${come}.e \
  -o /dev/null \
  -J ${ct}_${proj_name}_${come}_st3_ws \
  ./run_step3_can.sh ${project} ${ct} ${nmut} ${come} ${pth})

  sbatch -d $hold_jid -J ${ct}_${proj_name}_${come}_st3_log \
  -o ${logdir}/st3_${proj_name}_${ct}_${come}_${hold_jid}.o \
  -e /dev/null \
  log_job_slurm.sh ${proj_name} ${ct} ${nmut} ${pth} ${vmem} ${hr} ${come} \
  ${hold_jid} undefined undefined ${logdir} 3
done
