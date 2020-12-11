#! /bin/bash
set -eo pipefail

showUsage ()
{
  echo "usage:
  run_script_null_dist_C.sh -p project -m vmem -c cantype -t hrs -r n_rds
                            -s n_smpl -n n_perm -l logdir -f firstmtx -e endmtx

  parameters:
    -p : name of project folder
    -m : memory requested (default: 12G)
    -c : cancer type (default: all in <projectname>_CTs_list.txt file)
    -t : time requested in hours (default: 4)
    -r : number of RDS files (default: 200)
    -s : number of sampled matrices (default: 100)
    -n : number of permuted matrices (default: 1e+06)
    -l : log dir
    -f : number of first random mtx to process (default: 1)
    -e : number of last random mtx to process (default: s)

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
    -r)
        n_rds="$argVal"
    ;;
    -s)
        n_smpl="$argVal"
    ;;
    -n)
        n_perm="$argVal"
    ;;
    -l)
        logdir="$argVal"
    ;;
    -f)
        firstmtx="$argVal"
    ;;
    -e)
        endmtx="$argVal"
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
: ${vmem:="12G"}
: ${hr:="4"}
: ${n_rds:="200"}
: ${n_smpl:="100"}
: ${n_perm:="1000000"}
: ${logdir:="./ped_log"}
: ${firstmtx:="1"}
: ${endmtx:=${n_smpl}}


# FINISHED PARSING ARGUMENTS
#--------------------------------------------------------------------------

# create name for log file folder
d=`date +%Y%m%d_%H%M%S`

for ct in ${cantypes[@]}
do
    # create log directories
    proj_name=${project//[\/]/_}
    logdirstep=${logdir}/step4C_null_dist/${d}_${proj_name}_${ct}
    mkdir -p ${logdirstep}

    echo "working on": "${ct}"

    echo "#!/bin/bash
    module load R/3.4.1
    Rscript ./ped_scripts/null_dist_C.R ${ct} ${project} ${n_rds} ${n_smpl} ${n_perm} ${firstmtx} ${endmtx}" | \
    sbatch --time=${hr}:00:00 --mem=${vmem} -o ${logdirstep}/%j.o \
    -e ${logdirstep}/%j.e -J "J4C_"${proj_name}_${ct}

    sleep 1
done
