#!/bin/bash -e

# Script to run yara indexer.
# author: Fabian Buske
# date: July 2014

# QCVARIABLES,Resource temporarily unavailable
# RESULTFILENAME $${FASTA%.*}.sa.val


echo ">>>>> index generation for yara"
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"

function usage {
echo -e "usage: $(basename $0) -k [OPTIONS]"
exit
}


if [ ! $# -gt 1 ]; then usage ; fi

#INPUTS                                                                                                           
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the NGSANE repository                       
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
hash module 2>/dev/null && for MODULE in $MODULE_YARA; do module load $MODULE; done && module list 

export PATH=$PATH_YARA:$PATH
echo "PATH=$PATH"

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--yara       --\n "$(yara_indexer --version)
[ -z "$(which yara_indexer)" ] && echo "[ERROR] no yara detected" && exit 1

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
	dmget -a $(dirname $FASTA)/*
fi
    
NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "generating the index files"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    if [ ! -e ${FASTA%.*}.sa.val ]; then 
        echo "[NOTE] make index"; 
        yara_indexer --tmp-folder $TMP $FASTA; 
    fi
    if [ ! -e $FASTA.fai ]; then 
        echo "[NOTE] make .fai"; 
        samtools faidx $FASTA; 
    fi

    # mark checkpoint
    [ -f ${FASTA%.*}.sa.val ] && NGSANE_CHECKPOINT_CHECK 
fi 

NGSANE_CHECKPOINT_CHECK
################################################################################
[ -e ${FASTA%.*}.sa.val.dummy ] && rm ${FASTA%.*}.sa.val.dummy
echo ">>>>> index generation for yara - FINISHED"
echo ">>>>> enddate "`date`
