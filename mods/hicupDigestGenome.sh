#!/bin/bash -e

# Script running hicup including reference genome digestion, read mapping for single 
# and paired DNA reads with bowtie from fastq files
# It expects a fastq file, pairdend, reference genome and digest pattern  as input.
# author: Fabian Buske
# date: Apr 2013

# messages to look out for -- relevant for the QC.sh script:
# QCVARIABLES,Resource temporarily unavailable
# RESULTFILENAME <DIR>/<TASK>/digested_genome.txt

echo ">>>>> Digest genome with HiCUP "
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"

function usage {
echo -e "usage: $(basename $0) -k NGSANE -f FASTQ -r REFERENCE -o OUTDIR [OPTIONS]"
exit
}

if [ ! $# -gt 3 ]; then usage ; fi

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

DIGESTGENOME="digested_genome.txt"

if [ -z "$FASTA" ]; then
    echo "[ERROR] no reference provided (FASTA)"
    exit 1
fi

# delete old bam files unless attempting to recover
if [ -z "$RECOVERFROM" ]; then
    [ -e $OUTDIR/$DIGESTGENOME ] && rm $OUTDIR/$DIGESTGENOME
fi

if [ -z "$HICUP_RENZYMES" ]; then
   echo "[ERROR] No restriction enzyme given!" && exit 1
fi
ENZYMES=(${HICUP_RENZYMES//;/ })
ENZYME1=(${ENZYMES[0]//,/ })
ENZYME2=(${ENZYMES[1]//,/ })

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

    FASTABASE=${FASTA##*/}

    rm -f $OUTDIR/Digest_*
    cd $OUTDIR
    if [ ${#ENZYMES[@]} = 1 ]; then
       echo "Restriction Enzyme 1: ${ENZYME1[1]}:${ENZYME1[0]} "
       hicup_digester -g "${FASTABASE%.*}" -1 ${ENZYME1[0]} $FASTA
       mv Digest_* ${DIGESTGENOME}
    
    elif [ ${#ENZYMES[@]} = 2 ] && [ ! -e $OUTDIR/${FASTABASE%.*}_${ENZYME1[1]}_${ENZYME2[2]}.txt ]; then
       echo "Restriction Enzyme 1: ${ENZYME1[1]}:${ENZYME1[0]} "
       echo "Restriction Enzyme 2: ${ENZYME2[1]}:${ENZYME2[0]} "
       hicup_digester -g "${FASTABASE%.*}" -1 ${ENZYME1[0]} -2 ${ENZYME2[0]} $FASTA
       mv Digest_* ${DIGESTGENOME}
    else
       echo "[ERROR] Invalid number or pattern of enzyme digest patterns."
       exit 1
    fi
    cd $SOURCE
    
    # mark checkpoint
    if [ -f $OUTDIR/$DIGESTGENOME ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi
################################################################################
[ -e $OUTDIR/$DIGESTGENOME.dummy ] && rm $OUTDIR/$DIGESTGENOME.dummy
echo ">>>>> readmapping with hicup (bowtie) - FINISHED"
echo ">>>>> enddate "`date`

