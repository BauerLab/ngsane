#!/bin/bash -e

# Script marking duplicates 
# author: Tim Kahlke
# date: July 2016

# QCVARIABLES,Resource temporarily unavailable
# RESULTFILENAME <DIR>/<TASK>/<SAMPLE>.marked.bam

echo ">>>>> Mark duplicates"
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
        -f | --file )          shift; INPUTFILE=$1 ;;   # bam file 
        -o | --outdir )         shift; OUTDIR=$1 ;;     # output dir                                                     
        --recover-from )        shift; NGSANE_RECOVERFROM=$1 ;; # attempt to recover from log file
        -h | --help )           usage ;;
        * )                     echo "don't understand "$1
    esac
    shift
done


# Set defaults
MARKTOOL=MarkDuplicates
OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500
CREATE_INDEX=true


#PROGRAMS
. $CONFIG
. ${NGSANE_BASE}/conf/header.sh
. $CONFIG

################################################################################
NGSANE_CHECKPOINT_INIT "programs"

# save way to load modules that itself load other modules
hash module 2>/dev/null && for MODULE in $MODULE_GATK_MARKDUPS; do module load $MODULE; done && module list

echo "[NOTE] set java parameters"
JAVAPARAMS="-Xmx"$(python -c "print int($MEMORY_GATK_MARKDUPS*0.75)")"g -Djava.io.tmpdir="$TMP" -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1"
unset _JAVA_OPTIONS
echo "JAVAPARAMS "$JAVAPARAMS

PATH_PICARD_JAR=$(which picard.jar)
echo "PATH=$PATH"
echo -e "--JAVA        --\n" $(java -Xmx200m -version 2>&1)
[ -z "$(which java)" ] && echo "[ERROR] no java detected" && exit 1
[ ! -f $PATH_PICARD_JAR ] && echo "[ERROR] no picard detected" && exit 1

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)



NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "parameters"

# get basename of input file f
IFILE=${INPUTFILE##*/}
FILENAME=${IFILE%.*}

# delete old bam files unless attempting to recover
if [ -z "$NGSANE_RECOVERFROM" ]; then
    if [ -e $OUTDIR/$FILENAME.marked.bam ]; then
      rm $OUTDIR/$FILENAME.marked.bam
  fi  
fi

if [ -z "$METRICSNAME" ]; then
    echo "No metrics name given"
    exit 1
fi

# unique temp folder that should be used to store temporary files
THISTMP=$TMP"/"$(whoami)"/"$(echo $OUTDIR | md5sum | cut -d' ' -f1)
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
NGSANE_CHECKPOINT_INIT "marking duplicates"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    java $JAVAPARAMS -jar $PATH_PICARD_JAR $MARKTOOL \
    INPUT=$INPUTFILE \
    OUTPUT=$OUTDIR/$FILENAME.marked.bam \
    METRICS_FILE=$OUTDIR/$METRICSNAME \
    OPTICAL_DUPLICATE_PIXEL_DISTANCE=$OPTICAL_DUPLICATE_PIXEL_DISTANCE \
    CREATE_INDEX=$CREATE_INDEX \
    TMP_DIR=$THISTMP

    NGSANE_CHECKPOINT_CHECK $OUTDIR/$FILENAME.marked.bam 
fi

################################################################################
[ -e $OUTDIR/$FILENAME.marked.bam ] && rm $OUTDIR/*.marked.bam.dummy
echo ">>>>> creating unmapped bam - FINISHED"
echo ">>>>> enddate "`date`
