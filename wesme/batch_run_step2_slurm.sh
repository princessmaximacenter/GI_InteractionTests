#!/bin/bash -l

set -eo pipefail

showUsage ()
{
  echo "usage:
  batch_run_step2_slurm.sh -p project -n npermut -s nsampl -m vmem -c cantype
                           -t hrs -v pval -o come -f first -l last

  parameters:
    -p : name of project folder
    -n : number of permutations for FDR calculation
    -s : number of resamplings
    -m : memory requested (default: 10G)
    -c : cancer type (default: all in cantypes.txt file)
    -t : time requested in hours (default: 8)
    -v : p value threshold (default: 1.1)
    -o : come setting (default come, alternatives: co or me)
    -f : first job nb
    -l : last job nb
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
        first="$argVal"
    ;;
    -l)
        last="$argVal"
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
fi

# all cancer types
declare -a cantypes_def
readarray -t cantypes_def < ${project}/cantypes.txt

# set optional values
: ${vmem:="10"}
: ${hr:="8"}
: ${pth:="1.1"}
: ${come:="come"}
: ${first:="1"}
: ${last:=${nmut}}

# FINISHED PARSING ARGUMENTS
#--------------------------------------------------------------------------

# START OF SCRIPT

# create name for log file folder
d=`date +%Y%m%d_%H%M%S`

proj_name=${project//[\/]/_}

for ct in ${cantypes[@]}
do
  # create log directories
  logdir=./log/${project}_st2_${d}/${ct}
  mkdir -p ${logdir}/subjobs

  hold_jid=$(sbatch --array ${first}-${last}%50 --parsable \
  --time=${hr}:00:00 --mem=${vmem} \
  -e ${logdir}/subjobs/st2_${proj_name}_${ct}_${come}_\$TASK_ID.e \
  -o /dev/null \
  -J ${ct}_${proj_name}_${come}_st2_ws \
  ./run_step2_can_i_slurm.sh ${project} ${ct} ${nsampl} ${come} 0 ${pth} | awk -F. '{print $1}')

  sbatch -d $hold_jid -J ${ct}_${proj_name}_${come}_st2_log \
  -o ${logdir}/st2_${proj_name}_${ct}_${come}_${hold_jid}.o \
  -e /dev/null \
  log_job_slurm.sh ${proj_name} ${ct} ${nsampl} ${pth} ${vmem} ${hr} ${come} \
  ${hold_jid} ${first} ${last} ${logdir}/subjobs 2
done
