#!/bin/bash -e

# Script running hicup including reference genome digestion, read mapping for single 
# and paired DNA reads with bowtie from fastq files
# It expects a fastq file, pairdend, reference genome and digest pattern  as input.
# author: Fabian Buske
# date: Apr 2013

# messages to look out for -- relevant for the QC.sh script:
# QCVARIABLES,Resource temporarily unavailable
# RESULTFILENAME <DIR>/<TASK>/Digest_${REFERENCE_NAME}_${HICUP_RENZYME1}_${HICUP_RENZYME2}.txt

echo ">>>>> Digest genome with HiCUP "
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"

function usage {
echo -e "usage: $(basename $0) -k NGSANE -o OUTDIR [OPTIONS]"
exit
}

if [ ! $# -gt 2 ]; then usage ; fi

#INPUTS
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the NGSANE repository
        -o | --outdir )         shift; OUTDIR=$1 ;; # output dir
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

for MODULE in $MODULE_HICUP; do module load $MODULE; done  # save way to load modules that itself load other modules
export PATH=$PATH_HICUP:$PATH
module list
echo "PATH=$PATH"
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--perl        --\n "$(perl -v | grep "This is perl" )
[ -z "$(which perl)" ] && echo "[ERROR] no perl detected" && exit 1
echo -e "--HiCUP       --\n "$(hicup --version )
[ -z "$(which hicup)" ] && echo "[ERROR] no hicup detected" && exit 1

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="parameters"

if [ -z "$FASTA" ]; then
    echo "[ERROR] no reference provided (FASTA)"
    exit 1
fi

if [ -z "$HICUP_RENZYME1" ] || [ "${HICUP_RENZYME1,,}" == "none" ]; then
   echo "[ERROR] No restriction enzyme given!" && exit 1
elif [ -z "$HICUP_RCUTSITE1" ]; then
   echo "[ERROR] Restriction enzyme 1 lacks cutsite pattern!" && exit 1
fi
if [ -n "$HICUP_RENZYME2" ] && [ "${HICUP_RENZYME2,,}" != "none" ] && [ -z "$HICUP_RCUTSITE2" ]; then
   echo "[ERROR] Restriction enzyme 2 lacks cutsite pattern!" && exit 1
else
    HICUP_RENZYME2="none"
fi

if [ -z "$REFERENCE_NAME" ]; then
    echo "[ERROR] Reference assembly name not detected" && exit 1
fi

DIGESTGENOME="Digest_${REFERENCE_NAME}_${HICUP_RENZYME1}_${HICUP_RENZYME2}.txt"

# delete old bam files unless attempting to recover
if [ -z "$RECOVERFROM" ]; then
    [ -e $OUTDIR/$DIGESTGENOME ] && rm $OUTDIR/$DIGESTGENOME
fi

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="recall files from tape"

if [ -n "$DMGET" ]; then
	dmget -a $(dirname $FASTA)/*
	dmget -a $OUTDIR/*
fi

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="digest reference"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

    rm -f $OUTDIR/Digest_*
    if [ -z "$HICUP_RENZYME2" ] || [ "${HICUP_RENZYME2,,}" == "none" ]; then
       echo "Restriction Enzyme 1: $HICUP_RENZYME1:$HICUP_RCUTSITE1"
       RUN_COMMAND="$(which perl) $(which hicup_digester) --outdir $OUTDIR --genome $REFERENCE_NAME -1 $HICUP_RCUTSITE1,$HICUP_RENZYME1 $FASTA"
       echo $RUN_COMMAND && eval $RUN_COMMAND
       mv $OUTDIR/Digest_${REFERENCE_NAME}_${HICUP_RENZYME1}_*.txt $OUTDIR/${DIGESTGENOME}
    
    else
       echo "Restriction Enzyme 1: $HICUP_RENZYME1:$HICUP_RCUTSITE1 "
       echo "Restriction Enzyme 2: $HICUP_RENZYME2:$HICUP_RCUTSITE2 "
       RUN_COMMAND="$(which perl) $(which hicup_digester) --outdiroutdir $OUTDIR --genome $REFERENCE_NAME -1 $HICUP_RCUTSITE1,$HICUP_RENZYME1 -2 $HICUP_RCUTSITE2,$HICUP_RENZYME2 $FASTA"
       echo $RUN_COMMAND && eval $RUN_COMMAND
       mv $OUTDIR/Digest_${REFERENCE_NAME}_${HICUP_RENZYME1}_*.txt $OUTDIR/${DIGESTGENOME}
    fi
    
    # mark checkpoint
    if [ -f $OUTDIR/$DIGESTGENOME ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi
################################################################################
[ -e $OUTDIR/$DIGESTGENOME.dummy ] && rm $OUTDIR/$DIGESTGENOME.dummy
echo ">>>>> readmapping with hicup (bowtie) - FINISHED"
echo ">>>>> enddate "`date`

