#!/bin/bash -e

# Script running fit-hi-c to call significant chromatin interactions form HiC 
# experiments. Expected fragmentLists and contactCounts as input
# author: Fabian Buske
# date: Jan 2014

# messages to look out for -- relevant for the QC.sh script:
# QCVARIABLES,Resource temporarily unavailable
# RESULTFILENAME <DIR>/<TASK>/<SAMPLE>.spline_pass1.q05.txt

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

for MODULE in $MODULE_FITHIC; do module load $MODULE; done  # save way to load modules that itself load other modules
export PATH=$PATH_FITHIC:$PATH
module list
echo "PATH=$PATH"
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--Python      --\n" $(python --version)
[ -z "$(which python)" ] && echo "[ERROR] no python detected" && exit 1
echo -e "--Python libs --\n "$(yolk -l)
echo -e "--fit-hi-c    --\n "$(python $(which fit-hi-c.py) --version | head -n 1)
[ -z "$(which fit-hi-c.py)" ] && echo "[ERROR] no fit-hi-c detected" && exit 1

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="parameters"

# get basename of f
n=${f##*/}
SAMPLE=${n/%$FRAGMENTLIST/}

# delete old bam files unless attempting to recover
if [ -z "$RECOVERFROM" ]; then
    [ -d $OUTDIR/$SAMPLE ] && rm -r $OUTDIR/$SAMPLE
fi

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="recall files from tape"

if [ -n "$DMGET" ]; then
	dmget -a $OUTDIR/*
fi

echo -e "\n********* $CHECKPOINT\n"

################################################################################
CHECKPOINT="run fit-hi-c"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

    RUN_COMMAND="python $(which fit-hi-c.py) $FITHICADDPARAM --outdir $OUTDIR --fragments=$f --interactions=${f/%$FRAGMENTLIST/$CONTACTCOUNTS} --lib=${SAMPLE} &> $OUTDIR/$SAMPLE.log"
    echo $RUN_COMMAND && eval $RUN_COMMAND
    cat $OUTDIR/$SAMPLE.log # put into qout log too
    
    zcat $OUTDIR/${SAMPLE}.spline_pass1.pvals.txt.gz | awk '$7<=0.05' | sort -k7g | gzip > $OUTDIR/${SAMPLE}.spline_pass1.q05.txt.gz

    # mark checkpoint
    if [ -f $OUTDIR/${SAMPLE}.spline_pass1.q05.txt.gz ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
fi

################################################################################
[ -e $OUTDIR/${SAMPLE}.spline_pass1.q05.txt.gz.dummy ] && rm $OUTDIR/${SAMPLE}.spline_pass1.q05.txt.gz.dummy
echo ">>>>> Chromatin organization with fit-hi-c - FINISHED"
echo ">>>>> enddate "`date`

