#!/bin/bash

# author: Fabian Buske
# date: August 2013

echo ">>>>> Create wig files with wiggler"
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"

function usage {
echo -e "usage: $(basename $0) -k NGSANE -f bam -o OUTDIR [OPTIONS]

Script running wiggler on bam files

required:
  -k | --toolkit <path>     location of the NGSANE repository 
  -f | --bam <file>         bam file
  -o | --outdir <path>      output dir
"
exit
}
# QCVARIABLES,Resource temporarily unavailable

if [ ! $# -gt 3 ]; then usage ; fi

#INPUTS                                                                                                           
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the NGSANE repository                       
        -f | --bam )            shift; f=$1 ;; # bam file                                                       
        -o | --outdir )         shift; MYOUT=$1 ;; # output dir                                                     
        --recover-from )        shift; RECOVERFROM=$1 ;; # attempt to recover from log file
        -h | --help )           usage ;;
        * )                     echo "don't understand "$1
    esac
    shift
done

#PROGRAMS
. $CONFIG
. ${NGSANE_BASE}/conf/header.sh
. $CONFIG

################################################################################
CHECKPOINT="programs"

for MODULE in $MODULE_WIGGLER; do module load $MODULE; done  # save way to load modules that itself load other modules

export PATH=$PATH_WIGGLER:$PATH
module list
echo "PATH=$PATH"
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)
#echo -e "--wiggler  --\n "$(align2rawsignal 2>&1 | head -n 3 | tail -n 1)
#[ -z "$(which align2rawsignal)" ] && echo "[ERROR] wiggler not detected (align2rawsignal)" && exit 1

echo -e "\n********* $CHECKPOINT"
################################################################################
CHECKPOINT="parameters"

# get basename of f
n=${f##*/}

# check UMAP folder exist
if [ -z "$WIGGLER_UMAPDIR" ] || [ ! -d ${WIGGLER_UMAPDIR} ]; then
    echo "[ERROR] umap dir not specified not non-existant (WIGGLER_UMAPDIR)"
    exit 1
fi

# check output format
if [ -z "$WIGGLER_OUTPUTFORMAT" ]; then
    echo "[ERROR] wiggler output format not set" && exit 1
elif [ "$WIGGLER_OUTPUTFORMAT" != "bg" ] && [ "$WIGGLER_OUTPUTFORMAT" != "wig" ] && [ "$WIGGLER_OUTPUTFORMAT" != "mat" ]; then
    echo "[ERROR] wiggler output format not known" && exit 1
fi

mkdir -p ${MYOUT}

echo -e "\n********* $CHECKPOINT"
################################################################################
CHECKPOINT="recall files from tape"

for d in ${DIR[@]}; do
    FILES=$FILES" "$( ls $OUT/$d/$TASKBWA/*$ASD.bam )
done
echo $FILES

if [ -n "$DMGET" ]; then
	dmget -a $FILES
	dmls -l $FILES
fi

echo -e "\n********* $CHECKPOINT"
################################################################################
CHECKPOINT="run align2rawsignal"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep "********* $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
    INPUTS=""
    for FILE in $FILES; do
        INPUTS="${INPUTS} -i=$FILE"
    done
    
    RUN_COMMAND="align2rawsignal $WIGGLERADDPARAMS -of=$WIGGLER_OUTPUTFORMAT ${INPUTS} -s=${FASTA_CHROMDIR} -u=${WIGGLER_UMAPDIR} -v=${MYOUT}/wiggler-${n}.log -o=${MYOUT}/${n}.$WIGGLER_OUTPUTFORMAT -mm=$MEMORY_WIGGLER"
    echo $RUN_COMMAND && eval $RUN_COMMAND

    # mark checkpoint
    [ -f ${MYOUT}/${n}.$WIGGLER_OUTPUTFORMAT ] && echo -e "\n********* $CHECKPOINT" && unset RECOVERFROM
fi 

################################################################################
echo ">>>>> wiggler - FINISHED"
echo ">>>>> enddate "`date`

