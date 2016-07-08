#!/bin/bash -e

# Script for indel realignment of bam files
# author: Tim Kahlke
# date: July 2016

# QCVARIABLES,Resource temporarily unavailable
# RESULTFILENAME <DIR>/<TASK>/<SAMPLE>.ind_real.bam

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
        -f | --file )           shift; INPUTFILE=$1 ;; # fastq file 
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
hash module 2>/dev/null && for MODULE in $MODULE_GATK_INDREAL; do module load $MODULE; done && module list

echo "[NOTE] set java parameters"
JAVAPARAMS="-Xmx"$(python -c "print int($MEMORY_GATK_INDREAL*0.75)")"g -Djava.io.tmpdir="$TMP" -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1"
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



if [ -z $REFERENCE ]; then
    echo "[ERRPR] No reference file given"
    exit 1
fi

# get basename of reference
RNAME=${REFERENCE%.*}

if [ ! -f $RNAME.dict ]; then
    echo "[ERRPR] No dictionary for reference found. Create a dicitonary first"
    exit 1
fi

if [ -z $KNOWN ]; then
    echo "[ERRPR] No known indel file given"
    exit 1
fi


# delete old bam files unless attempting to recover
if [ -z "$NGSANE_RECOVERFROM" ]; then
    if [ -e $OUTDIR/$FILENAME.ind_real.bam ]; then
      rm $OUTDIR/$FILENAME.ind_real.bam
  fi  
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
NGSANE_CHECKPOINT_INIT "create target intervals"




if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    java $JAVAPARAMS -jar $PATH_GATK_JAR -T RealignerTargetCreator \
    -R $REFERENCE \
    -nt 20 \
    -known $KNOWN \
    -I $INPUTFILE \
    -o $OUTDIR/$FILENAME.realignertargets.intervals

    NGSANE_CHECKPOINT_CHECK $OUTDIR/$FILENAME.realignertargets.intervals 
fi

################################################################################
NGSANE_CHECKPOINT_INIT "realign indels"




if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    java $JAVAPARAMS -jar $PATH_GATK_JAR -T IndelRealigner \
    -R $REFERENCE \
    -nt 20 \
    -targetIntervals $OUTDIR/$FILENAME.realignertargets.intervals \
    -known $KNOWN \
    -I $INPUTFILE \
    -o $OUTDIR/$FILENAME.ind_real.bam

    NGSANE_CHECKPOINT_CHECK $OUTDIR/$FILENAME.ind_real.bam
fi





################################################################################
[ -e $OUTDIR/$FILENAME.ind_real.bam ] && rm $OUTDIR/*.ind_real.bam.dummy
echo ">>>>> creating unmapped bam - FINISHED"
echo ">>>>> enddate "`date`
