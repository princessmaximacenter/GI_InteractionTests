#! /bin/bash
set -eo pipefail

showUsage ()
{
  echo "usage:
  run_script_PAN_A.sh -p project -m vmem -t hrs -l logdir

  parameters:
    -p : name of project folder
    -m : memory requested (default: 5G)
    -t : time requested in hours (default: 2)
    -l : log dir
  "
}

#--------------------------------------------------------------------------
# PARSE THE arguments

ARGS=("$@")
skipNext=0
COUNTER=0

while [ ${COUNTER} -lt $# ]
do
  arg=${ARGS[${COUNTER}]}
  let COUNTER=COUNTER+1
  nextArg=${ARGS[${COUNTER}]}

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
    -m)
        vmem="$argVal"
    ;;
    -t)
        hr="$argVal"
    ;;
    -l)
        logdir="$argVal"
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

# set optional values
: ${vmem:="5G"}
: ${hr:="2"}
: ${logdir:="./ped_log"}

# FINISHED PARSING ARGUMENTS
#--------------------------------------------------------------------------

# START OF SCRIPT

# create name for log file folder
d=`date +%Y%m%d_%H%M%S`

# create log directories
proj_name=${project//[\/]/_}
logdirstep=${logdir}/step2A_PAN/${d}_${proj_name}_PAN
mkdir -p ${logdirstep}

echo "working on PAN"
echo "#!/bin/bash
module load R/3.4.1
Rscript ./ped_scripts/GI_detection_permutation_PAN_A.R ${project}" | \
sbatch --time=${hr}:00:00 --mem=${vmem} -o ${logdirstep}/%j.o \
-e ${logdirstep}/%j.e -J "J2A_"${proj_name}"_PAN"
