#!/bin/bash -e

# Script for ChIP-seq QC using Chance.
# It takes read alignments in .bam format.
# It produces output files: peak regions in bed format
# author: Fabian Buske
# date: August 2013

echo ">>>>> ChIPseq QC with Chance"
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"

function usage {
echo -e "usage: $(basename $0) -k NGSANE -f FASTQ -o OUTDIR [OPTIONS]"
exit
}

# QCVARIABLES,Resource temporarily unavailable
# RESULTFILENAME <SAMPLE>_refinepeak.bed

if [ ! $# -gt 3 ]; then usage ; fi

#INPUTS                                                                                                           
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the NGSANE repository
        -f | --bam )            shift; f=$1 ;; # bam file
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

for MODULE in $MODULE_CHANCE; do module load $MODULE; done  # save way to load modules that itself load other modules

export PATH=$PATH_CHANCE:$PATH
module list
echo "PATH=$PATH"
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--chance      --\n "$(chance --version 2>&1)
[ -z "$(which run_chance_com.sh)" ] && echo "[ERROR] run_chance_com.sh not detected" && exit 1
echo -e "--Matlab (MCR)--\n "$(echo $MCRROOT)
if [ -z "$(echo $MCRROOT)" ] && echo "[ERROR] no matlab runtime environment detected" && exit 1
echo -e "--R           --\n "$(R --version | head -n 3)
[ -z "$(which R)" ] && echo "[ERROR] no R detected" && exit 1

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="parameters"

# get basename of f
n=${f##*/}
c=${CHIPINPUT##*/}

if [ -z "$CHIPINPUT" ] || [ ! -f $CHIPINPUT ]; then
    echo "[ERROR] input control not provided or invalid (CHIPINPUT)"
    exit 1
fi

if [ -z "$GENOME_ASSEMBLY" ]; then 
    echo "[ERROR] GENOME_ASSEMBLY not specified."
    exit 1
fi 

CHANCE="`which run_chance_com.sh` `echo $MCRROOT`"
echo "[NOTE] Chance Environment: $CHANCE"

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="recall files from tape"

if [ -n "$DMGET" ]; then
	dmget -a ${f}
	[ -n $CHIPINPUT ] && dmget -a $CHIPINPUT
fi

echo -e "\n********* $CHECKPOINT\n"

################################################################################
CHECKPOINT="compute IPstrength"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else

    RUN_COMMAND="${CHANCE} IPStrength -b ${GENOME_ASSEMBLY} -t bam -o ${OUTDIR}/${n}-${c}.IPstrength --ipfile $f --ipsample ${n} --inputfile ${CHIPINPUT} --inputsample ${c}"
    echo $RUN_COMMAND && eval $RUN_COMMAND
    
    # mark checkpoint
    if [ -f ${OUTDIR}/${n}-${c}.IPstrength ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi

################################################################################
CHECKPOINT="compare with ENCODE"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else

    RUN_COMMAND="${CHANCE} compENCODE -b ${GENOME_ASSEMBLY} -t bam -o ${OUTDIR}/${n}-${c}.compENCODE -e ${EXPERIMENTID} --ipfile ${f} --ipsample ${n} --inputfile ${CHIPINPUT} --inputsample ${c}"
    echo $RUN_COMMAND && eval $RUN_COMMAND
    
    # mark checkpoint
    if [ -f ${OUTDIR}/${n}-${c}.compENCODE ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi

################################################################################
CHECKPOINT="summary"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else
    
    Rscript --vanilla --quiet ${NGSANE_BASE}/tools/makeChancePlots.R $f $CHIPINPUT ${OUTDIR}

    # mark checkpoint
    if [ -f ${n/.$ASD.bam/_treat_pileup.bb} ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi

################################################################################
#[ -e $OUTDIR/${n/.$ASD.bam/_refinepeak.bed}.dummy ] && rm $OUTDIR/${n/.$ASD.bam/_refinepeak.bed}.dummy
echo ">>>>> ChIPseq QC with Chance - FINISHED"
echo ">>>>> enddate "`date`

