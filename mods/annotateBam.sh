#!/bin/bash -e

# Annotation of bam files
# author: Denis C. Bauer
# date: Jan.2013

# messages to look out for -- relevant for the QC.sh script:
# QCVARIABLES,Resource temporarily unavailable
# RESULTFILENAME <DIR>/$TASK_BAMANN/<SAMPLE>.anno.bed

echo ">>>>> Annotate BAM file "
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"


function usage {
echo -e "usage: $(basename $0) -k CONFIG -f BAM -o OUTDIR [OPTIONS]

Annotating BAM file with annotations in a folder specified in the config
file

required:
  -k | --toolkit <path>     config file
  -f <file>                 bam file

options:
"
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
hash module 2>/dev/null && for MODULE in $MODULE_BAMANN; do module load $MODULE; done && module list 

export PATH=$PATH_BAMANN:$PATH
echo "PATH=$PATH"

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--bedtools    --\n "$(bedtools -version)
[ -z "$(which bedtools)" ] && echo "[WARN] bedtools not detected" && exit 1

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "parameters"

# get basename of f
n=${f##*/}
SAMPLE=${n/%$ASD.bam/}

# check library variables are set
if [[ -z "$BAMANNLIB" ]]; then
    echo "[ERROR] library info not set (BAMANNLIB)"
    exit 1;
else
    echo "[NOTE] using libraries in $BAMANNLIB"
fi

# delete old bam files unless attempting to recover
if [ -z "$NGSANE_RECOVERFROM" ]; then
    [ -e $OUTDIR/$SAMPLE.anno.bed ] && rm $OUTDIR/$SAMPLE.anno.bed
    [ -e $OUTDIR/$SAMPLE.anno.stats ] && rm $OUTDIR/$SAMPLE.anno.stats
fi

NAMES="MergedReads"
NUMBER=$( ls $BAMANNLIB | grep -P "(.gtf|.bed)$" | wc -l )
ANNFILES=""

for i in $(ls $BAMANNLIB | grep -P "(.gtf|.bed)$" | tr '\n' ' '); do
	NAME=$(basename $i)
	NAMES=$NAMES" "${NAME%*.*}
	ANNFILES=$ANNFILES" "$BAMANNLIB/$i
done

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "recall files from tape"
	
if [ -n "$DMGET" ]; then
    dmget -a $BAMANNLIB/*
	dmget -a ${f}
fi
    
NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "bedmerge"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

	bedtools bamtobed -i $f | bedtools merge -i - -n  > $OUTDIR/$SAMPLE.merg.bed

	# mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE.merg.bed 

fi
 
################################################################################
NGSANE_CHECKPOINT_INIT "run annotation"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

	echo "[NOTE] annotate with $NAMES $ANNFILES"

   	# annotated counts and add column indicated non-annotated read counts
   	bedtools annotate -counts -i $OUTDIR/$SAMPLE.merg.bed -files $ANNFILES -names $NAMES | sed 's/[ \t]*$//g' | awk '{FS=OFS="\t";if (NR==1){print $0,"unannotated"} else{ for(i=5; i<=NF;i++) j+=$i; if (j>0){print $0,0}else{print $0,1}; j=0}}' > $OUTDIR/$SAMPLE.anno.bed
	
	# mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE.anno.bed 
    
    [ -e f.merg.bed  ] && rm $OUTDIR/$SAMPLE.merg.bed
fi
################################################################################
NGSANE_CHECKPOINT_INIT "summarize"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

   	COMMAND="python ${NGSANE_BASE}/tools/sumRows.py -i $OUTDIR/$SAMPLE.anno.bed -l 3 -s 3 -e $(expr $NUMBER + 5) -n $(echo $NAMES | sed 's/[\t ]/,/g'),unannotated > $OUTDIR/$SAMPLE.anno.stats"
	echo $COMMAND && eval $COMMAND

#    head -n 1 $OUTDIR/$SAMPLE.anno.bed >> $OUTDIR/$SAMPLE.merg.anno.stats
    cat $OUTDIR/$SAMPLE.anno.bed | sort -k4gr | head -n 20 >> $OUTDIR/$SAMPLE.anno.stats  
    
	# mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE.anno.stats 

fi
################################################################################
[ -e $OUTDIR/$SAMPLE.anno.bed.dummy ] && rm $OUTDIR/$SAMPLE.anno.bed.dummy
echo ">>>>> Annotate BAM file - FINISHED"
echo ">>>>> enddate "`date`

