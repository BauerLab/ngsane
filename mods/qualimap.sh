#!/bin/bash -e

# Script for QC of bam files using qualimap v2.
# It takes read alignments in .bam format.
# It produces html report
# author: Fabian Buske
# date: July 2014

# QCVARIABLES,Resource temporarily unavailable
# RESULTFILENAME <DIR>/<TASK>/<SAMPLE>-bamQC/qualimapReport.html


echo ">>>>> BAM QC with qualimap"
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"

function usage {
echo -e "usage: $(basename $0) -k NGSANE -f FASTQ -o OUTDIR [OPTIONS]"
exit
}


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

# save way to load modules that itself loads other modules
hash module 2>/dev/null && for MODULE in $MODULE_QUALIMAP; do module load $MODULE; done && module list 

export PATH=$PATH_QUALIMAP:$PATH
echo "PATH=$PATH"
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)
[ -z "$PATH_QUALIMAP" ] && PATH_QUALIMAP=$(dirname $(which qualimap.jar))

echo "[NOTE] set java parameters"
JAVAPARAMS="-Xmx"$(python -c "print int($MEMORY_QUALIMAP*0.8)")"g -Djava.io.tmpdir="$TMP" -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1" 
unset _JAVA_OPTIONS
unset DISPLAY
echo "JAVAPARAMS "$JAVAPARAMS
export JAVA_OPTS="$JAVAPARAMS"

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--qualimap    --\n "$(qualimap --version 2>&1 | tee | grep -v -P "QualiMap v")
[ -z "$(which qualimap)" ] && echo "[ERROR] qualimap not detected" && exit 1
echo -e "--R           --\n "$(R --version | head -n 3)
[ -z "$(which R)" ] && echo "[ERROR] no R detected" && exit 1

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="parameters"

# get basename of f
f=${f/%.dummy/} #if input came from pip
n=${f##*/}
SAMPLE=${n/%.$ASD.bam/}

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="recall files from tape"

if [ -n "$DMGET" ]; then
	dmget -a ${f}
	dmget -a $OUTDIR/*
fi

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="qualimap bamqc"

cd $OUTDIR
if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else
    qualimap bamqc -bam $f $QUALIMAP_BAMQC_ADDPARAM -nt $CPU_QUALIMAP -outformat HTML -outdir $OUTDIR/$SAMPLE-bamQC
        
    # mark checkpoint
    if [ -f $OUTDIR/$SAMPLE-bamQC/qualimapReport.html ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi
################################################################################
CHECKPOINT="qualimap rnaseq-qc"

cd $OUTDIR
if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else
        
    if [[ -n "$QUALIMAP_RNASEQ" ]]; then

        qualimap rnaseq -bam $f $QUALIMAP_RNASEQ_ADDPARAM -nt $CPU_QUALIMAP -outformat HTML -counts $OUTDIR/$SAMPLE"_"counts.txt -outdir $OUTDIR/$SAMPLE-rnaseqQC
           
        # mark checkpoint
        if [ -f $OUTDIR/$SAMPLE-rnaseqQC/report.html ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
        
    else
        echo "[NOTE] rnaseq qc skipped"
        echo -e "\n********* $CHECKPOINT\n"
    fi
    
fi
################################################################################
[ -e $OUTDIR/$SAMPLE-bamQC/qualimapReport.html.dummy ] && rm $OUTDIR/$SAMPLE-bamQC/qualimapReport.html.dummy
echo ">>>>> BAM QC with qualimap - FINISHED"
echo ">>>>> enddate "`date`

