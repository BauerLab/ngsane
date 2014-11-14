#!/bin/bash -e

# Script running fit-hi-c to call significant chromatin interactions form HiC 
# experiments. Expects bam files as input.
# author: Fabian Buske
# date: Oct 2014

# messages to look out for -- relevant for the QC.sh script:
# QCVARIABLES,Resource temporarily unavailable
# RESULTFILENAME <DIR>/<TASK>/<SAMPLE>.txt.gz

echo ">>>>> Chromatin organization with fit-hi-c "
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
        -f | --file )           shift; f=$1 ;; # input file
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
hash module 2>/dev/null && for MODULE in $MODULE_FITHIC; do module load $MODULE; done && module list 

export PATH=$PATH_FITHIC:$PATH
echo "PATH=$PATH"
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--Python      --\n" $(python --version)
[ -z "$(which python)" ] && echo "[ERROR] no python detected" && exit 1
hash module 2>/dev/null && echo -e "--Python libs --\n "$(yolk -l)
echo -e "--fit-hi-c    --\n "$(python $(which fit-hi-c.py) --version | head -n 1)
[ -z "$(which fit-hi-c.py)" ] && echo "[ERROR] no fit-hi-c detected" && exit 1

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "parameters"

# get basename of f
n=${f##*/}
SAMPLE=${n/%$ASD.bam/}

# delete old bam files unless attempting to recover
if [ -z "$NGSANE_RECOVERFROM" ]; then
    [ -d $OUTDIR/$SAMPLE ] && rm -r $OUTDIR/$SAMPLE
    [ -f $OUTDIR/$SAMPLE.fragmentLists.gz ] && rm $OUTDIR/$SAMPLE.fragmentLists.gz
    [ -f $OUTDIR/$SAMPLE.contactCounts.gz ] && rm $OUTDIR/$SAMPLE.contactCounts.gz
    [ -f $OUTDIR/$SAMPLE.ice.txt.gz ] && rm $OUTDIR/$SAMPLE.ice.txt.gz
    [ -f $OUTDIR/$SAMPLE.log ] && rm $OUTDIR/$SAMPLE.log
fi

GENOME_CHROMSIZES=${FASTA%.*}.chrom.sizes
if [ ! -f $GENOME_CHROMSIZES ]; then
    echo "[WARN] GENOME_CHROMSIZES not found. No bigbeds will be generated"
else
    echo "[NOTE] Chromosome size: $GENOME_CHROMSIZES"
fi

if [ -z "$MAPPABILITY" ]; then
    echo "[ERROR] Mappability not specified"
    exit 1
fi

if [ -z "$HIC_RESOLUTION" ]; then
    echo "[ERROR] HiC resolution not specified"
    exit 1
fi

if [[ -z "$FITHIC_QVALUETHRESHOLD" ]]; then
    FITHIC_QVALUETHRESHOLD=0.01
fi
echo "[NOTE] Q-value threshold: $FITHIC_QVALUETHRESHOLD"

if [[ -z "$FITHIC_MAPPABILITYTHRESHOLD" ]];then
    echo "[ERROR] FITHIC_MAPPABILITYTHRESHOLD not set"
    exit 1
fi

if [[ -n "$FITHIC_CHROMOSOMES" ]]; then
    FITHIC_CHROMOSOMES="--chrompattern '$FITHIC_CHROMOSOMES'"
fi

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "recall files from tape"

if [ -n "$DMGET" ]; then
	dmget -a $f
	dmget -a $OUTDIR/*
fi

NGSANE_CHECKPOINT_CHECK

################################################################################
NGSANE_CHECKPOINT_INIT "count Interactions"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    RUN_COMMAND="python ${NGSANE_BASE}/tools/fithicCountInteractions.py --mappability=$MAPPABILITY --resolution=$HIC_RESOLUTION --chromsizes=$GENOME_CHROMSIZES $FITHIC_CHROMOSOMES --outputDir=$OUTDIR $f >> $OUTDIR/$SAMPLE.log"
    echo $RUN_COMMAND && eval $RUN_COMMAND

    [ -e $OUTDIR/${SAMPLE}$ASD.bam.fragmentLists ] && mv $OUTDIR/${SAMPLE}$ASD.bam.fragmentLists $OUTDIR/$SAMPLE.fragmentLists
    [ -e $OUTDIR/${SAMPLE}$ASD.bam.contactCounts ] && mv $OUTDIR/${SAMPLE}$ASD.bam.contactCounts $OUTDIR/$SAMPLE.contactCounts
    
    $GZIP $OUTDIR/$SAMPLE.fragmentLists $OUTDIR/$SAMPLE.contactCounts

    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE.fragmentLists.gz $OUTDIR/$SAMPLE.contactCounts.gz

fi

################################################################################
NGSANE_CHECKPOINT_INIT "ICE correction"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    RUN_COMMAND="python ${NGSANE_BASE}/tools/fithic-fixedBins/ICE-with-sparseMatrix.py $OUTDIR/$SAMPLE.contactCounts.gz $OUTDIR/$SAMPLE.fragmentLists.gz l1 $OUTDIR/$SAMPLE.ice.txt.gz 0.5 >> $OUTDIR/$SAMPLE.log"
    echo $RUN_COMMAND && eval $RUN_COMMAND
   
    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE.ice.txt.gz

fi

################################################################################
NGSANE_CHECKPOINT_INIT "fit-hi-c"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then
    cd $OUTDIR/$RESOLUTION
    RUN_COMMAND="python ${NGSANE_BASE}/tools/fithic-fixedBins/fit-hi-c-fixedSize-withBiases.py $FITHICADDPARAM --lib=${SAMPLE} --biases=$OUTDIR/$SAMPLE.ice.txt.biases.gz --fragments=$OUTDIR/$SAMPLE.fragmentLists.gz --interactions=$OUTDIR/$SAMPLE.contactCounts.gz --resolution $HIC_RESOLUTION >> $OUTDIR/$SAMPLE.log"
    echo $RUN_COMMAND && eval $RUN_COMMAND
    
    zcat $OUTDIR/$SAMPLE.spline_pass1.res$HIC_RESOLUTION.significances.txt.gz | awk -v q=$FITHIC_QVALUETHRESHOLD '$7<=q' | sort -k7g | gzip > $OUTDIR/$SAMPLE.txt.gz

    SIGCISINTERACTIONS=$(zcat $OUTDIR/$SAMPLE.txt.gz |  awk '$1==$3' | wc -l | cut -d' ' -f 2)
    SIGTRANSINTERACTIONS=$(zcat $OUTDIR/$SAMPLE.txt.gz |  awk '$1!=$3' | wc -l | cut -d' ' -f 2)
    echo "Significant cis interactions: $SIGCISINTERACTIONS" >> $OUTDIR/$SAMPLE.log
    echo "Significant trans interactions: $SIGTRANSINTERACTIONS" >> $OUTDIR/$SAMPLE.log
    
    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE.txt.gz
fi

################################################################################
[ -e $OUTDIR/$SAMPLE.txt.gz.dummy ] && rm $OUTDIR/$SAMPLE.txt.gz.dummy
echo ">>>>> Chromatin organization with fit-hi-c - FINISHED"
echo ">>>>> enddate "`date`

