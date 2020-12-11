#! /bin/bash
set -eo pipefail

showUsage ()
{
  echo "usage:
  batch_run_script_PAN_B.sh -p project -m vmem -t hrs -s startseed -o stopseed -n nperm -l logdir

  parameters:
    -p : name of project folder
    -m : memory requested (default: 300G)
    -t : time requested in hours (default: 24)
    -s : start seed (default: 1)
    -o : stop seed (default: 200)
    -n : number of permutations (default: 5000)
    -l : logdir
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
    -s)
        startseed="$argVal"
    ;;
    -o)
        stopseed="$argVal"
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

# set optional values
: ${vmem:="300G"}
: ${hr:="24"}
: ${startseed:="1"}
: ${stopseed:="200"}
: ${nperm:="5000"}
: ${logdir:="./ped_log"}

# FINISHED PARSING ARGUMENTS
#--------------------------------------------------------------------------

# create name for log file folder
d=`date +%Y%m%d_%H%M%S`

echo "working on PAN"
# create log directories
proj_name=${project//[\/]/_}
logdirstep=${logdir}/step2_null_dist/${d}_${proj_name}_PAN/
mkdir -p ${logdirstep}

sbatch --array=${startseed}-${stopseed}%100 --time=${hr}:00:00 --mem=${vmem} \
-e ${logdirstep}/"J2_PAN"_\%A.e\%a \
-o ${logdirstep}/"J2_PAN"_\%A.o\%a \
-J "J2_"${proj_name}"_PAN" \
./ped_scripts/run_script_PAN_B.sh ${project} ${nperm}
