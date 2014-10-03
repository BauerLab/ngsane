#!/bin/bash -e

# Script for QC of bam files using qualimap v2.
# It takes read alignments in .bam format.
# It produces html report
# author: Fabian Buske
# date: July 2014

# QCVARIABLES,Resource temporarily unavailable

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
hash module 2>/dev/null && for MODULE in $MODULE_QUALIMAP; do module load $MODULE; done && module list 

export PATH=$PATH_QUALIMAP:$PATH
echo "PATH=$PATH"
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)
[ -z "$PATH_QUALIMAP" ] && PATH_QUALIMAP=$(dirname $(which qualimap.jar))

echo "[NOTE] set java parameters"
JAVAPARAMS="-Xmx"$(python -c "print int($MEMORY_QUALIMAP*0.75)")"g -Djava.io.tmpdir="$TMP" -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1" 
unset _JAVA_OPTIONS
unset DISPLAY
echo "JAVAPARAMS "$JAVAPARAMS
export JAVA_OPTS="$JAVAPARAMS"

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--qualimap    --\n "$(qualimap --version 2>&1 | tee | grep -v -P "QualiMap v")
[ -z "$(which qualimap)" ] && echo "[ERROR] qualimap not detected" && exit 1
echo -e "--R           --\n "$(R --version | head -n 3)
[ -z "$(which R)" ] && echo "[ERROR] no R detected" && exit 1

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "parameters"

# get basename of f
f=${f/%.dummy/} #if input came from pip
n=${f##*/}
SAMPLE=${n/%$ASD.bam/}

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "recall files from tape"

if [ -n "$DMGET" ]; then
	dmget -a ${f}
	dmget -a $OUTDIR/*
fi

NGSANE_CHECKPOINT_CHECK
###############################################################################
NGSANE_CHECKPOINT_INIT "qualimap bamqc"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    cd $OUTDIR
    RUNCOMMAND="qualimap bamqc -bam $f $QUALIMAP_BAMQC_ADDPARAM -nt $CPU_QUALIMAP -outformat HTML -outdir $OUTDIR/$SAMPLE-bamQC"
    echo $RUNCOMMAND && eval $RUNCOMMAND        
    
    mv $OUTDIR/$SAMPLE-bamQC/genome_results.txt $OUTDIR/${SAMPLE}_bamQC.txt
    
    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/${SAMPLE}_bamQC.txt

fi
################################################################################
NGSANE_CHECKPOINT_INIT "qualimap rnaseq-qc"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    cd $OUTDIR    
    if [[ -n "$QUALIMAP_RNASEQ" ]]; then

        RUNCOMMAND="qualimap rnaseq -bam $f $QUALIMAP_RNASEQ_ADDPARAM -outformat HTML -counts $OUTDIR/$SAMPLE'_'counts.txt -outdir $OUTDIR/$SAMPLE-rnaseqQC"
        echo $RUNCOMMAND && eval $RUNCOMMAND
        
        echo "Aligned to genes = " $(cat $OUTDIR/$SAMPLE-rnaseqQC/qualimapReport.html | fgrep "Aligned to genes" -A 1 | tail -n 1 | sed 's/.*>\(.*\)<.*/\1/' ) > $OUTDIR/${SAMPLE}_rnaseqQC.txt
        echo "Non-unique = " $(cat $OUTDIR/$SAMPLE-rnaseqQC/qualimapReport.html | fgrep "Non-unique alignment" -A 1 | tail -n 1 | sed 's/.*>\(.*\)<.*/\1/' ) >> $OUTDIR/${SAMPLE}_rnaseqQC.txt
        echo "Ambiguous = " $(cat $OUTDIR/$SAMPLE-rnaseqQC/qualimapReport.html | fgrep "Ambiguous alignment" -A 1 | tail -n 1 | sed 's/.*>\(.*\)<.*/\1/' ) >> $OUTDIR/${SAMPLE}_rnaseqQC.txt
        echo "No feature assigned = " $(cat $OUTDIR/$SAMPLE-rnaseqQC/qualimapReport.html | fgrep "No feature assigned" -A 1 | tail -n 1 | sed 's/.*>\(.*\)<.*/\1/' ) >> $OUTDIR/${SAMPLE}_rnaseqQC.txt
        echo "Not aligned = " $(cat $OUTDIR/$SAMPLE-rnaseqQC/qualimapReport.html | fgrep "Not aligned" -A 1 | tail -n 1 | sed 's/.*>\(.*\)<.*/\1/' ) >> $OUTDIR/${SAMPLE}_rnaseqQC.txt

        echo "5 prime bias = " $(cat $OUTDIR/$SAMPLE-rnaseqQC/qualimapReport.html | fgrep "5' bias" -A 1 | tail -n 1 | sed 's/.*>\(.*\)<.*/\1/' ) >> $OUTDIR/${SAMPLE}_rnaseqQC.txt
        echo "3 prime bias = " $(cat $OUTDIR/$SAMPLE-rnaseqQC/qualimapReport.html | fgrep "3' bias" -A 1 | tail -n 1 | sed 's/.*>\(.*\)<.*/\1/' ) >> $OUTDIR/${SAMPLE}_rnaseqQC.txt
        echo "5-3 bias = " $(cat $OUTDIR/$SAMPLE-rnaseqQC/qualimapReport.html | fgrep "5'-3' bias" -A 1 | tail -n 1 | sed 's/.*>\(.*\)<.*/\1/' ) >> $OUTDIR/${SAMPLE}_rnaseqQC.txt
        echo "Reads with junctions = " $(cat $OUTDIR/$SAMPLE-rnaseqQC/qualimapReport.html | fgrep "Reads with junctions" -A 1 | tail -n 1 | sed 's/.*>\(.*\)<.*/\1/' ) >> $OUTDIR/${SAMPLE}_rnaseqQC.txt

        # fix image names
        for i in $OUTDIR/${SAMPLE}-rnaseqQC/images_qualimapReport/*.png; do
            mv "$i" "${i// /_}"
            sed -i "s|$i|${i// /_}|g" $OUTDIR/$SAMPLE-rnaseqQC/qualimapReport.html
        done
            

        # mark checkpoint
        NGSANE_CHECKPOINT_CHECK $OUTDIR/${SAMPLE}_rnaseqQC.txt
        
    else
        echo "[NOTE] rnaseq qc skipped"
        NGSANE_CHECKPOINT_CHECK
    fi
    
fi
################################################################################
echo ">>>>> BAM QC with qualimap - FINISHED"
echo ">>>>> enddate "`date`

