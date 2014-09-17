#!/bin/bash

# cufflinks calling script
# author: Denis C. Bauer
# date: March.2011

# messages to look out for -- relevant for the QC.sh script:
# QCVARIABLES,

echo ">>>>> differential expression with cuffdiff"
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"


function usage {
echo -e "usage: $(basename $0) -k NGSANE -f -o OUTDIR [OPTIONS]
"
exit
}


if [ ! $# -gt 3 ]; then usage ; fi

#DEFAULTS
THREADS=1

#INPUTS
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the NGSANE repository
        -b | --basename )       shift; fs=$1 ;; # basename
        -o | --outdir )         shift; OUTDIR=$1 ;; # output dir
        --recover-from )        shift; NGSANE_RECOVERFROM=$1 ;; # attempt to recover from log file
        -h | --help )           usage ;;
        * )                     usage
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
hash module 2>/dev/null && for MODULE in $MODULE_CUFFDIFF; do module load $MODULE; done && module list 

export PATH=$PATH_CUFFDIFF:$PATH
echo "PATH=$PATH"

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--cuffdiff    --\n "$(cuffdiff 2>&1 | head -n 1 )
[ -z "$(which cuffdiff)" ] && echo "[ERROR] no cuffdiff detected" && exit 1
echo -e "--cuffcompare --\n "$(cuffcompare 2>&1 | head -n 1 )
[ -z "$(which cuffcompare)" ] && echo "[ERROR] no cuffcompare detected" && exit 1

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "parameters"

# delete old bam files unless attempting to recover
if [ -z "$NGSANE_RECOVERFROM" ]; then
    if [ -d $OUTDIR ]; then rm -r $OUTDIR; fi
    mkdir -p $OUTDIR
fi

FILES=${FILES//,/ }

echo "[NOTE] Files: $FILES"
DATASETS=""
egrep '^REPLICATE ' $CONFIG | cut -d' ' -f 2- > $COMMAND.tmp

while read -r -a REPLICATE; do
    for f in $FILES; do
        echo $f
    #     get directory
    #    d=$(dirname $f)
    #    d=${d##*/}    # add to dataset
        if [ -n "$FILE" ]; then 
            DATASETS="${DATASETS[@]} ${FILE[@]}"
        fi
    done
done < $COMMAND.tmp
exit 1
cd $OUTDIR/
GTF=$REFSEQGTF
    
################################################################################
NGSANE_CHECKPOINT_INIT "recall files from tape"

if [ -n "$DMGET" ]; then
    dmget -a $(dirname $TOPHATOUT)/*
    dmget -a ${OUTDIR}/*
fi

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "check reference GTF"
    
if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    # Which transcript reference?
    if [ -n $GTF ]; then
        echo "[NOTE] compare to oneanother"
        
        cuffcompare -o comp.txt $CUFGTFS
        GTF=$OUTDIR/comp.combined.gtf
        
        # mark checkpoint
        NGSANE_CHECKPOINT_CHECK $OUTDIR/comp.combined.gtf
    else
    
        echo "[NOTE] skip checking reference GTF"
        # mark checkpoint
        NGSANE_CHECKPOINT_CHECK
    fi    
fi

################################################################################
NGSANE_CHECKPOINT_INIT "run cuff diff"
    
if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    cuffdiff -r $FASTA -p $CPU_CUFFDIFF -o $OUTDIR $GTF $TOPHATBAM

    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK
fi 

################################################################################
echo ">>>>> differential expression with cuffdiff - FINISHED"
echo ">>>>> enddate "`date`
