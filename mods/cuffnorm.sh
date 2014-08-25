#!/bin/bash -e

# author: Hugh French and Fabian Buske
# date: April 2014 echo ">>>>> [cuffnorm]"
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
        -f | --file )           shift; FILES=$1 ;;  # input file                                                       
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

hash module 2>/dev/null && for MODULE in $MODULE_CUFFLINKS; do module load $MODULE; done  # save way to load modules that itself load other modules
export PATH=$PATH_CUFFLINKS:$PATH

module list
echo "PATH=$PATH"

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--cufflinks   --\n "$(cufflinks 2>&1 | tee | head -n 2 )
[ -z "$(which cufflinks)" ] && echo "[ERROR] no cufflinks detected" && exit 1

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "parameters"

echo "[NOTE] Files: $FILES"
DATASETS=""
for f in $FILES; do
    n=$(dirname ${f})
#    n=$(dirname ${f/%$ASD.bam/.cxb}
    FILE=${n/$TASK_TOPHAT/$TASK_CUFFLINKS}/abundances.cxb

#     get directory
#    d=$(dirname $f)
#    d=${d##*/}    # add to dataset
    if [ -n "$FILE" ]; then 
        DATASETS="${DATASETS[@]} ${FILE[@]}"
    fi
done

echo "[NOTE] ${FILE[@]}"

echo "[NOTE] datasets"
echo "[NOTE] $DATASETS"

echo "[NOTE] $OUTDIR"

#mkdir -p "$OUTDIR"
#if [ -z "$NGSANE_RECOVERFROM" ]; then
#    ## TODO remove primary result files from pervious runs
#    rm ${OUTDIR}/*
#fi

# unique temp folder that should be used to store temporary files
THISTMP=$TMP"/"$(whoami)"/"$(echo $OUTDIR | md5sum | cut -d' ' -f1)
mkdir -p "$THISTMP"

#echo "[NOTE] echo $THISTMP"

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "recall files from tape"

if [ -n "$DMGET" ]; then
    dmget -a $INPUTFILE
    dmget -a $OUTDIR/*
fi
    
NGSANE_CHECKPOINT_CHECK
#################################################################################
NGSANE_CHECKPOINT_INIT "Run cuffnorm"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then
    
    RUNCOMMAND="cuffnorm --no-update-check --quiet --output-dir ${OUTDIR}/$MERGED_GTF_NAME -p $CPU_CUFFLINKS $OUT/expression/$TASK_CUFFLINKS/$MERGED_GTF_NAME.gtf $(echo $FILES | tr ',' ' ')"
    echo $RUNCOMMAND && eval $RUNCOMMAND

    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK && $OUTDIR/genes.count_table $OUTDIR/genes.fpkm_table

fi
#################################################################################
NGSANE_CHECKPOINT_INIT "cleanup"  

[ -f ${THISTMP}/files.txt ] && rm ${THISTMP}/files.txt
     
NGSANE_CHECKPOINT_CHECK
#################################################################################
echo ">>>>> Experiment merged transcripts (cuffnorm) - FINISHED"
echo ">>>>> enddate "`date`
