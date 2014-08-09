#!/bin/bash -e

# author: Hugh French and Fabian Buske
# date: March 2014 echo ">>>>> [cuffnorm]"
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

# save way to load modules that itself loads other modules
hash module 2>/dev/null && for MODULE in $MODULE_CUFFLINKS; do module load $MODULE; done && module list 

export PATH=$PATH_CUFFLINKS:$PATH
echo "PATH=$PATH"

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--cufflinks   --\n "$(cufflinks 2>&1 | tee | head -n 2 )
[ -z "$(which cufflinks)" ] && echo "[ERROR] no cufflinks detected" && exit 1

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "parameters"

if [ -z "$MERGED_GTF_NAME" ]; then
    echo "[ERROR] MERGED_GTF_NAME not specified"
    exit 1
fi

if [ -z "$FASTA" ]; then
    echo "[ERROR] no reference provided (FASTA)"
    exit 1
fi

echo "[NOTE] Files: $FILES"
OLDFS=$IFS
IFS=","
DATASETS=""
for f in $FILES; do
    # get basename of f

    n=${f/%$ASD.bam/}
    FILE=${n/$TASK_TOPHAT/$TASK_CUFFLINKS}
    # get directory
    d=$(dirname $f)
    d=${d##*/}    # add to dataset
    if [ -n "$FILE" ]; then 
        DATASETS="${DATASETS[@]} ${FILE[@]}"
    fi
done
IFS=" "

echo "[NOTE] datasets: $DATASETS"

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
    # TODO add additional resources that are required and may need recovery from tape
fi
    
NGSANE_CHECKPOINT_CHECK   
################################################################################
NGSANE_CHECKPOINT_INIT "Run cuffmerge"  

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    array=(${DATASETS[@]})
    
    
    [ -f ${THISTMP}/files.txt ] &&  rm ${THISTMP}/files.txt
    touch ${THISTMP}/files.txt
    
    array=( "${array[@]/%//transcripts.gtf}" )                  
    
    for THIS_FILE in "${array[@]}"; do
        [ -f $THIS_FILE ] && echo $THIS_FILE >> ${THISTMP}/files.txt 
    done
    
    cat ${THISTMP}/files.txt
    
    RUNCOMMAND="cuffmerge -p $CPU_CUFFLINKS -o $THISTMP --ref-sequence $FASTA --ref-gtf $GTF ${THISTMP}/files.txt"
    echo $RUNCOMMAND && eval $RUNCOMMAND
    
    mv $THISTMP/merged.gtf $OUTDIR/$MERGED_GTF_NAME.gtf
    
    # mark checkpoint
    [[ -s $OUTDIR/$MERGED_GTF_NAME.gtf ]] && NGSANE_CHECKPOINT_CHECK

fi
################################################################################
NGSANE_CHECKPOINT_INIT "cleanup."  
  
[ -f ${THISTMP}/files.txt ] && rm ${THISTMP}/files.txt
  
NGSANE_CHECKPOINT_CHECK
################################################################################
echo ">>>>> Experiment merged transcripts (cuffmerge) - FINISHED"
echo ">>>>> enddate "`date`




