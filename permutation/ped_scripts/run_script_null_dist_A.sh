#! /bin/bash
set -eo pipefail

showUsage ()
{
  echo "usage:
  run_script_null_dist_A.sh -p project -m vmem -c cantype -t hrs -n n_mtx -l logdir

  parameters:
    -p : name of project folder
    -m : memory requested (default: 30G)
    -c : cancer type (default: all in <projectname>_CTs_list.txt file)
    -t : time requested in hours (default: 8)
    -n : number of random sampled matrices (default: 100)
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
    -c)
        cantypes="$argVal"
    ;;
    -t)
        hr="$argVal"
    ;;
    -n)
        nperm="$argVal"
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

if [[ -z "$cantypes" ]];
then
  # check if <projectname>_CTs_list.txt file exists in project folder
  # if not return error
  if [ -e ./ped_results/${project}_CTs_list.txt ];
  then
    # else read all cancer types
    declare -a cantypes_def
    readarray -t cantypes_def < ./ped_results/${project}_CTs_list.txt
    : ${cantypes:=${cantypes_def[@]}}
  else
    echo "ERROR: Missing both parameter -c cancertype and ${project}_CTs_list.txt file"
    showUsage
    exit
  fi
fi

# set optional values
: ${vmem:="30G"}
: ${hr:="8"}
: ${nperm:="100"}
: ${logdir:="./ped_log"}

# FINISHED PARSING ARGUMENTS
#--------------------------------------------------------------------------

# create name for log file folder
d=`date +%Y%m%d_%H%M%S`

for ct in ${cantypes[@]}
do
    # create log directories
    proj_name=${project//[\/]/_}
    logdirstep=${logdir}/step4A_null_dist/${d}_${proj_name}_${ct}
    mkdir -p ${logdirstep}

    echo "working on": "${ct}"
  
    echo "#!/bin/bash
    module load R/3.4.1
    Rscript ./ped_scripts/null_dist_A.R ${ct} ${project} ${nperm}" | \
    sbatch --time=${hr}:00:00 --mem=${vmem} -o ${logdirstep}/%j.o \
    -e ${logdirstep}/%j.e -J "J4A_"${proj_name}_${ct}

    sleep 1
done
