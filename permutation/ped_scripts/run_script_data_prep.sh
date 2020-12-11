#!/bin/bash -l

set -eo pipefail

showUsage ()
{
  echo "usage:
  run_script_data_prep.sh -p project -f inputfile -l logdir -s use_sub -r has_header
  parameters:
    -p : name of project folder
    -f : path to inputfile
    -l : log dir (default: ./ped_log)
    -s : use cancer subtype (TRUE/FALSE, default: FALSE)
    -r : file has header (TRUE/FALSE, default: TRUE)
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
    -f)
        inputfile="$argVal"
    ;;
    -l)
        logdir="$argVal"
    ;;
    -s)
        use_sub="$argVal"
    ;;
    -r)
        has_header="$argVal"
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

if [[ -z "$inputfile" ]];
then
  echo "ERROR: Missing parameter -i input file"
  showUsage
  exit
fi

# set optional values
: ${logdir:="./ped_log"}
: ${use_sub:="FALSE"}
: ${has_header:="TRUE"}

# FINISHED PARSING ARGUMENTS
#--------------------------------------------------------------------------

# create name for log file folder
d=`date +%Y%m%d_%H%M%S`

# create log directories
proj_name=${project//[\/]/_}
logdirstep=${logdir}/step1_data_prep/${d}_${proj_name}
mkdir -p ${logdirstep}

echo "#!/bin/bash
module load R/3.4.1
Rscript ./ped_scripts/data_prep_ped.R ${inputfile} ${project} ${use_sub} ${has_header}" | \
sbatch --time=1:00:00 --mem=1G -o ${logdirstep}/%j.o -e ${logdirstep}/%j.e -J "J1_"${proj_name}
