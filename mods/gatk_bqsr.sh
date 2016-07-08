#!/bin/bash -e

# Script for base qulaity score recalibration (GATK_BQSR
# author: Tim Kahlke
# date: July 2016

# QCVARIABLES,Resource temporarily unavailable
# RESULTFILENAME <DIR>/<TASK>/<SAMPLE>.bqs_recalibrated.bam

echo ">>>>> fastq to ubam conversion"
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
        -f | --file )          shift; INPUTFILE=$1 ;; # fastq file 
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

# save way to load modules that itself load other modules
hash module 2>/dev/null && for MODULE in $MODULE_GATK_BQSR; do module load $MODULE; done && module list

echo "[NOTE] set java parameters"
JAVAPARAMS="-Xmx"$(python -c "print int($MEMORY_GATK_BQSR*0.75)")"g -Djava.io.tmpdir="$TMP" -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1"
unset _JAVA_OPTIONS
echo "JAVAPARAMS "$JAVAPARAMS

PATH_GATK_JAR=$(which GenomeAnalysisTK.jar)
echo "PATH=$PATH"
echo -e "--JAVA        --\n" $(java -Xmx200m -version 2>&1)
[ -z "$(which java)" ] && echo "[ERROR] no java detected" && exit 1
[ ! -f $PATH_GATK_JAR ] && echo "[ERROR] no gatk detected" && exit 1


echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)



NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "parameters"

# get basename of input file f
IFILE=${INPUTFILE##*/}
FILENAME=${IFILE%.*}



# delete old bam files unless attempting to recover
if [ -z "$NGSANE_RECOVERFROM" ]; then
    if [ -e $OUTDIR/$FILENAME.bqs_recalibrated.bam ]; then
      rm $OUTDIR/$FILENAME.bqs_recalibrated.bam
  fi  
fi

if [ -z "$REFERENCE" ]; then
    echo "No reference file given"
    exit 1
fi

RNAME=${REFERENCE%.*}
if [ ! -f $RNAME.dict ]; then
    echo "[ERRPR] No dictionary for reference found. Create a dicitonary first"
    exit 1
fi


if [ -z "$KNOWN" ]; then
    echo "No vcf files of known sites given"
    exit 1
fi

KSTRING=""
IFS="," read -ra KNOWNS <<<$KNOWN
for k in $KNOWNS; do
    KSTRING="$KSTRING -knownSites $k" 
done

echo "knowns = $KSTRING"

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
NGSANE_CHECKPOINT_INIT "creating recal table"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    java $JAVAPARAMS -jar $PATH_GATK_JAR -T BaseRecalibrator \
    -R $REFERENCE \
    -nt 20 \
    -I $INPUTFILE \
    -o $OUTDIR/$FILENAME.recal_data.table \
    $KSTRING

    NGSANE_CHECKPOINT_CHECK $OUTDIR/$FILENAME.recal_data.table 
fi

################################################################################
NGSANE_CHECKPOINT_INIT "creating create covariation table"


if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    java $JAVAPARAMS -jar $PATH_GATK_JAR -T BaseRecalibrator \
    -R $REFERENCE \
    -I $INPUTFILE \
    -nt 20 \
    -BQSR $OUTDIR/$FILENAME.recal_data.table \
    -o $OUTDIR/$FILENAME.post_recal_data.table \
    $KSTRING

    NGSANE_CHECKPOINT_CHECK $OUTDIR/$FILENAME.post_recal_data.table    
fi

################################################################################
NGSANE_CHECKPOINT_INIT "create before/after plots"

# Taken out due to R error


if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    java $JAVAPARAMS -jar $PATH_GATK_JAR -T AnalyzeCovariates \
    -R $REFERENCE \
    -before $OUTDIR/$FILENAME.recal_data.table \
    -after $OUTDIR/$FILENAME.post_recal_data.table \
    -plots $OUTDIR/$FILENAME.recalibration_plots.pdf

    NGSANE_CHECKPOINT_CHECK $OUTDIR/$FILENAME.reclibration_plots.pdf    
fi
################################################################################
NGSANE_CHECKPOINT_INIT "apply recalibration "

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    java $JAVAPARAMS -jar $PATH_GATK_JAR -T PrintReads \
    -R $REFERENCE \
    -I $INPUTFILE \
    -nt 20 \
    -o $OUTDIR/$FILENAME.bqs_recalibrated.bam \
    -BQSR $OUTDIR/$FILENAME.recal_data.table

    NGSANE_CHECKPOINT_CHECK $OUTDIR/$FILENAME.bqs_recalibrated.bam    
fi

################################################################################
[ -e $OUTDIR/$FILENAME.bqs_recalibrated.bam ] && rm $OUTDIR/*.bqs_recalibrated.bam.dummy
echo ">>>>> GATK_BQSR - FINISHED"
echo ">>>>> enddate "`date`
