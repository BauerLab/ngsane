#!/bin/bash -e

# Script for ChIP-seq peak calling using MACS v2.
# It takes read alignments in .bam format.
# It produces output files: peak regions in bed format
# author: Fabian Buske
# date: August 2013

echo ">>>>> ChIPseq analysis with MACS2"
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"

function usage {
echo -e "usage: $(basename $0) -k NGSANE -f FASTQ -r REFERENCE -o OUTDIR [OPTIONS]"
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

for MODULE in $MODULE_MACS2; do module load $MODULE; done  # save way to load modules that itself load other modules

export PATH=$PATH_MACS2:$PATH
module list
echo "PATH=$PATH"
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)

echo -e "--macs2      --\n "$(macs2 --version)
[ -z "$(which macs2)" ] && echo "[ERROR] macs2 not detected" && exit 1
echo -e "--R          --\n "$(R --version | head -n 3)
[ -z "$(which R)" ] && echo "[ERROR] no R detected" && exit 1

echo -e "\n********* $CHECKPOINT"
################################################################################
CHECKPOINT="parameters"

# get basename of f
n=${f##*/}
c=${CHIPINPUT##*/}

if [ -z "$CHIPINPUT" ] || [ ! -f $CHIPINPUT ]; then
    echo "[WARN] input control not provided or invalid (CHIPINPUT)"
    unset CHIPINPUT
else
    CHIPINPUT="--control CHIPINPUT"
fi

echo -e "\n********* $CHECKPOINT"
################################################################################
CHECKPOINT="recall files from tape"

if [ -n "$DMGET" ]; then
	dmget -a ${f}
	dmls -l ${f}
fi

echo -e "\n********* $CHECKPOINT"
################################################################################
CHECKPOINT="macs 2"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep "********* $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else
    
    cd $MYOUT
    
    echo "[NOTE] data quality"
    macs2 predictd $MACS2_PREDICTD_ADDPARAM --ifile $f --gsize $MACS2_GENOMESIZE --rfile ${n/.$ASD.bam/}_predicted
    Rscript ${n}_predicted_model.R
    
    echo "[NOTE] call peaks"
    macs2 callpeak $MACS2_CALLPEAK_ADDPARAM --treatment $f $CHIPINPUT --gsize $MACS2_GENOMESIZE --name ${n/.$ASD.bam/}
    Rscript ${n}_model.R
    
    echo "[NOTE] refine peaks"
    macs2 refinepeak $MACS2_REFINEPEAK_ADDPARAM -b ${n/.$ASD.bam/}_peaks.bed -i $f $CHIPINPUT --gsize $GENOMESIZE --o-prefix ${n/.$ASD.bam/}_refined
    echo $RUN_COMMAND && eval $RUN_COMMAND
       
    echo "Peaks: `tail -n+2 ${n/.$ASD.bam/}_peaks.bed | wc -l | awk '{print $1}'`" >> ${n/.$ASD.bam/}.summary.txt
    echo "Summits: `tail -n+2 ${n/.$ASD.bam/}_summits.bed | wc -l | awk '{print $1}'`" >> ${n/.$ASD.bam/}.summary.txt

    cd $SOURCE
    
    # mark checkpoint
    [ -f ${n/.$ASD.bam/}_peaks.bed ] && echo -e "\n********* $CHECKPOINT" && unset RECOVERFROM
fi

################################################################################
CHECKPOINT="zip"

# $GZIP $MYOUT/${n/.$ASD.bam/}-${c/.$ASD.bam/}_details

echo -e "\n********* $CHECKPOINT"
################################################################################
echo ">>>>> ChIPseq analysis with MACS2 - FINISHED"
echo ">>>>> enddate "`date`

