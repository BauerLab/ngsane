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

for MODULE in $MODULE_HICLIB; do module load $MODULE; done  # save way to load modules that itself load other modules

export PATH=$PATH_HICLIB:$PATH
module list
echo "PATH=$PATH"

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--Python      --\n" $(python --version)
echo -e "--Python libs --\n "$(yolk -l)

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="parameters"

PARAMS="--gapFile=$HICLIB_GAPFILE --referenceGenome=$FASTA_CHROMDIR"

echo "[NOTE] Files: $FILES"
OLDFS=$IFS
IFS=","
DATASETS=""
for f in $FILES; do
    # get basename of f
    n=${f##*/}
    n=${n/%$READONE.$FASTQ/}
    # get directory
    d=$(dirname $f)
    d=${d##*/}
    # get hdf5 file
    FILE=$(ls $SOURCE/$d/hiclib/$n-fragment_dataset.hdf5)
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

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="recall files from tape"

if [ -n "$DMGET" ]; then
	dmget -a $OUTDIR/*
fi
    
echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="call hiclib"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
    
    python ${NGSANE_BASE}/tools/hiclibCorrelate.py ${PARAMS} --outputDir=$OUT/runStats/$TASK_HICLIB --tmpDir=$THISTMP --verbose $DATASETS
    
    # mark checkpoint
    echo -e "\n********* $CHECKPOINT\n"
fi
################################################################################
CHECKPOINT="cleanup"

[ -d $THISTMP ] && rm -r $THISTMP

fi
################################################################################
echo ">>>>> HiC correlation analysis with hiclib - FINISHED"
echo ">>>>> enddate "`date`
