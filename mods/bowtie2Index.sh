#!/bin/bash -e

# Script to run bowtie program.
# It takes comma-seprated list of files containing short sequence reads in fasta or fastq format and bowtie index files as input.
# It produces output files: read alignments in .bam format and other files.
# author: Denis Bauer
# date: June 2012
# modified: August 2013 Fabian Buske

# QCVARIABLES,Resource temporarily unavailable
# RESULTFILENAME ${FASTA%.*}.1.bt2

echo ">>>>> index generation for Bowtie2 "
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

# save way to load modules that itself loads other modules
hash module 2>/dev/null && for MODULE in $MODULE_BOWTIE2; do module load $MODULE; done && module list 

export PATH=$PATH_BOWTIE2:$PATH
echo "PATH=$PATH"

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--bowtie2     --\n "$(bowtie2 --version)
[ -z "$(which bowtie2)" ] && echo "[ERROR] no bowtie2 detected" && exit 1

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="parameters"

if [ -z "$FASTA" ]; then
    echo "[ERROR] no reference provided (FASTA)"
    exit 1
fi

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="recall files from tape"
	
if [ -n "$DMGET" ]; then
    dmget -a $FASTA*
fi
    
echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="generating the index files"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
    
    if [[ ! -e ${FASTA%.*}.1.bt2 ]]; then 
        echo "[NOTE] make .bt2"; 
        bowtie2-build $FASTA ${FASTA%.*}; 
    fi
    if [ ! -e $FASTA.fai ]; then 
        echo "[NOTE] make .fai"; 
        samtools faidx $FASTA; 
    fi

    # mark checkpoint
    if [ -f ${FASTA%.*}.1.bt2 ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
fi 
################################################################################
[ -e ${FASTA%.*}.1.bt2.dummy ] && rm ${FASTA%.*}.1.bt2.dummy
echo ">>>>> index generation for Bowtie2 - FINISHED"
echo ">>>>> enddate "`date`
