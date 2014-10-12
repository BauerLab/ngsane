#!/bin/bash -e

# BWA calling script
# author: Fabian Buske
# date: Nov.2013

# messages to look out for -- relevant for the QC.sh script:
# QCVARIABLES,We are loosing reads,MAPQ should be 0 for unmapped read,no such file,file not found,bwa.sh: line,Resource temporarily unavailable
# RESULTFILENAME $FASTA.bwt

echo ">>>>> index generation for BWA "
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"

function usage {
echo -e "usage: $(basename $0) -k NGSANE"
exit
}

if [ ! $# -gt 1 ]; then usage ; fi

#INPUTS
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the NGSANE repository
        -o | --output )         shift; LOGFILEOUTPUT=$1 ;; # Not used but passed on by prepareJobsubmission for log name (commontask)
        --recover-from )        shift; NGSANE_RECOVERFROM=$1 ;; # attempt to recover from log file
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
NGSANE_CHECKPOINT_INIT "programs"

# save way to load modules that itself loads other modules
hash module 2>/dev/null && for MODULE in $MODULE_BWA; do module load $MODULE; done && module list 

export PATH=$PATH_BWA:$PATH
echo "PATH=$PATH"

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--bwa         --\n "$(bwa 2>&1 | head -n 3 | tail -n-2)
[ -z "$(which bwa)" ] && echo "[ERROR] no bwa detected" && exit 1

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "parameters"

if [ -z "$FASTA" ]; then
    echo "[ERROR] no reference provided (FASTA)"
    exit 1
fi

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "recall files from tape"

if [ -n "$DMGET" ]; then
	dmget -a $FASTA*
fi
    
NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "generating the index files"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then
    
    # generating the index files
    if [ ! -e $FASTA.bwt ]; then 
        echo "[NOTE] make .bwt"; 
        command="bwa index -a ${BWA_INDEXAS} $FASTA"
	echo $command && eval $command
    fi
    if [ ! -e $FASTA.fai ]; then 
        echo "[NOTE] make .fai"; 
        command="samtools faidx $FASTA"
	echo $command && eval $command
    fi

    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $FASTA.bwt
fi 

################################################################################
[ -e $FASTA.bwt.dummy ] && rm $FASTA.bwt.dummy
echo ">>>>> index generation for - FINISHED"
echo ">>>>> enddate "`date`

