#!/bin/bash -e

# author: Hugh French and Fabian Buske
# date: April 2014 echo ">>>>> [cuffquant]"
# RESULTFILENAME <DIR>/<TASK>/<SAMPLE>.cxb
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"

function usage {
echo -e "usage: $(basename $0) -k NGSANE -f INPUTFILE -o OUTDIR [OPTIONS]"
exit
}


if [ ! $# -gt 3 ]; then usage ; fi

#INPUTS                                                                                                           
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;;     # location of the NGSANE repository                       
        -f | --file )           shift; INPUTFILE=$1 ;;  # input file                                                       
        -o | --outdir )         shift; OUTDIR=$1 ;;     # output dir                                                     
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

for MODULE in $MODULE_CUFFLINKS; do module load $MODULE; done  # save way to load modules that itself load other modules
export PATH=$PATH_CUFFLINKS:$PATH

module list
echo "PATH=$PATH"

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--cufflinks   --\n "$(cufflinks 2>&1 | tee | head -n 2 )
[ -z "$(which cufflinks)" ] && echo "[ERROR] no cufflinks detected" && exit 1

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "parameters"

[ ! -f $INPUTFILE ] && echo "[ERROR] input file not found: $INPUTFILE" && exit 1

# get basename of f (samplename)
n=${INPUTFILE##*/}
SAMPLE=${n/%$ASD.bam/}

if [ -z "$MERGED_GTF_NAME" ]; then
    echo "[ERROR] MERGED_GTF_NAME not specified"
    exit 1
elif [ ! -f $OUT/expression/$TASK_CUFFLINKS/$MERGED_GTF_NAME.gtf ]; then 
    echo "[ERROR] merged gtf file not found in $OUT/expression/$TASK_CUFFLINKS/$MERGED_GTF_NAME.gtf "
    exit 1
fi

# unique temp folder that should be used to store temporary files
THISTMP=$TMP"/"$(whoami)"/"$(echo $OUTDIR/$SAMPLE | md5sum | cut -d' ' -f1)
[ -d $THISTMP ] && rm -r $THISTMP
mkdir -p $THISTMP

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "recall files from tape"

if [ -n "$DMGET" ]; then
    dmget -a $INPUTFILE
    dmget -a $OUTDIR/*
fi
    
NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "Run cuffquant"  
if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then
    
    mkdir -p $OUTDIR/$SAMPLE
    
    RUNCOMMAND="cuffquant --no-update-check --quiet --output-dir $OUTDIR/$SAMPLE -p $CPU_CUFFLINKS $OUT/expression/$TASK_CUFFLINKS/$MERGED_GTF_NAME.gtf $INPUTFILE"
    echo $RUNCOMMAND && eval $RUNCOMMAND

    mv $OUTDIR/$SAMPLE/abundances.cxb $OUTDIR/$SAMPLE.cxb

    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE.cxb
fi

################################################################################
NGSANE_CHECKPOINT_INIT "cleanup"

[ -f ${THISTMP}/files.txt ] && rm ${THISTMP}/files.txt
[ -d $OUTDIR/$SAMPLE ] && rm -r $OUTDIR/$SAMPLE

NGSANE_CHECKPOINT_CHECK
################################################################################
[ -e $OUTDIR/$SAMPLE.cxb.dummy ] && rm $OUTDIR/$SAMPLE.cxb.dummy
echo ">>>>> Experiment merged transcripts (cuffquant) - FINISHED"
echo ">>>>> enddate "`date`




