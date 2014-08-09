#!/bin/bash -e

echo ">>>>> HiC correlation analysis with hiclib "
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"

function usage {
echo -e "usage: $(basename $0) -k CONFIG -f FASTQ -r REFERENCE -o OUTDIR [OPTIONS]"
exit
}

# The script correlates all hiclib library in a pairwise manner.
# author: Fabian Buske
# date: August 2013

# QCVARIABLES,Resource temporarily unavailable

if [ ! $# -gt 3 ]; then usage ; fi

#INPUTS                                                                                                           
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the NGSANE repository
        -f | --files ) 		    shift; FILES=$1 ;; # input files
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
hash module 2>/dev/null && for MODULE in $MODULE_HICLIB; do module load $MODULE; done && module list 

export PATH=$PATH_HICLIB:$PATH
echo "PATH=$PATH"

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--Python      --\n" $(python --version)
hash module 2>/dev/null && echo -e "--Python libs --\n "$(yolk -l)

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "parameters"

PARAMS="--gapFile=$HICLIB_GAPFILE --referenceGenome=$FASTA_CHROMDIR"

echo "[NOTE] Files: $FILES"
OLDFS=$IFS
IFS=","
DATASETS=""
for f in $FILES; do
    # get basename of f
    n=${f##*/}
    $SAMPLE=${n/%$READONE.$FASTQ/}
    # get directory
    d=$(dirname $f)
    d=${d##*/}
    # get hdf5 file
    FILE=$(ls $SOURCE/$d/$TASK_HICLIB/$SAMPLE-fragment_dataset.hdf5)
    # add to dataset
    if [ -n "$FILE" ]; then 
    	DATASETS="${DATASETS[@]} ${FILE[@]}"
    fi
done
IFS=$OLDFS

THISTMP=$TMP"/"$(whoami)"/"$(echo $OUTDIR/$SAMPLE | md5sum | cut -d' ' -f1)
[ -d $THISTMP ] && rm -r $THISTMP
mkdir -p $THISTMP

echo "[NOTE] Datasets: $DATASETS"

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "recall files from tape"

if [ -n "$DMGET" ]; then
	dmget -a $OUTDIR/*
fi
    
NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "call hiclib"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then
    
    python ${NGSANE_BASE}/tools/hiclibCorrelate.py ${PARAMS} --outputDir=$OUT/runStats/$TASK_HICLIB --tmpDir=$THISTMP --verbose $DATASETS
    
    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK
fi
################################################################################
NGSANE_CHECKPOINT_INIT "cleanup"

[ -d $THISTMP ] && rm -r $THISTMP

fi
################################################################################
echo ">>>>> HiC correlation analysis with hiclib - FINISHED"
echo ">>>>> enddate "`date`
