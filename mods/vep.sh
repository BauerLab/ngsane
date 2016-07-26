#!/bin/bash -e

# Script for annotating variants using vep 
# author: Tim Kahlke
# date: July 2016

# QCVARIABLES,Resource temporarily unavailable
# RESULTFILENAME <DIR>/<TASK>/<SAMPLE>.vcf

echo ">>>>> Variant annotation using VEP"
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
        -f | --fastq )          shift; INPUTFILE=$1 ;; # fastq file 
        -o | --outdir )         shift; OUTDIR=$1 ;;     # output dir                                                     
        --recover-from )        shift; NGSANE_RECOVERFROM=$1 ;; # attempt to recover from log file
        -h | --help )           usage ;;
        * )                     echo "don't understand "$1
    esac
    shift
done


#DEFAULT
LOFTEE=1
PARAMS="--everything"
DBNSFP=1
SPECIES="homo_sapiens"


#PROGRAMS
. $CONFIG
. ${NGSANE_BASE}/conf/header.sh
. $CONFIG

################################################################################
NGSANE_CHECKPOINT_INIT "programs"

# save way to load modules that itself load other modules
hash module 2>/dev/null && for MODULE in $MODULE_VEP; do module load $MODULE; done && module list

echo "[NOTE] set java parameters"
JAVAPARAMS="-Xmx"$(python -c "print int($MEMORY_VEP*0.75)")"g -Djava.io.tmpdir="$TMP" -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1"
unset _JAVA_OPTIONS
echo "JAVAPARAMS "$JAVAPARAMS

PATH_VEP_PL=$(which variant_effect_predictor.pl)
echo "PATH=$PATH"
[ ! -f $PATH_VEP_PL ] && echo "[ERROR] No variant_effect_predictor.pl detected" && exit 1

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)



NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "parameters"

# get basename of input file f
IFILE=${INPUTFILE##*/}
FILENAME=${IFILE%.*}


if [ -z $CACHE_DIR ]; then
    echo "[ERRPR] No cache dir given"
    exit 1
fi

# delete old bam files unless attempting to recover
if [ -z "$NGSANE_RECOVERFROM" ]; then
    if [ -e $OUTDIR/$FILENAME.ubam.bam ]; then
      rm $OUTDIR/$FILENAME.ubam.bam
  fi  
fi

if [ -z $PLUGIN_DIR ]; then
    echo "[ERRPR] No plugin dir given"
    exit 1;
fi

if [ -z $FORMAT ]; then
    echo "[ERRPR] no format of input file given"
    exit 1
fi


if [ -n $LOFTEE && -z $LOFTEE_ANCESTOR ]; then
    echo "[ERRPR] No path to loftee human_ancestor.fa given"
    exit 1
fi

if [ -n $LOFTEE_ANCESTOR && ! -f "$LOFTEE_ANCESTOR.fai"]; then
    echo "No index for $LOFTEE_ANCESTOR found. Use samtools to create an index first"
    exit 1
fi

if [ -n $DBNSFP && -z $DBNSFP_DB ]; then
    echo "No database file given for dbNSFP"
    exit 1
fi

if [-n $DBNSFP_DP && ! -f "$DBNSFP_DP.tbi" ]; then
    echo "No index file for dbNSFP found. Please use tabix to create an index first"
    exit 1
fi

OPT_STRING=""





echo "String of optional params = $OPT_STRING"

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
NGSANE_CHECKPOINT_INIT "creating ubams"




if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    java $JAVAPARAMS -jar $PATH_PICARD_JAR FastqToSam $OPT_STRING \
    FASTQ=$f \
    FASTQ2=${f/%$READONE.$FASTQ/$READTWO.$FASTQ} \
    OUTPUT=$OUTDIR/$FILENAME.ubam.bam \
    READ_GROUP_NAME= $READGROUP \
    SAMPLE_NAME=$SAMPLENAME \
    LIBRARY_NAME=$LIBRARY \
    PLATFORM=$PLATFORM 


    NGSANE_CHECKPOINT_CHECK $OUTDIR/$FILENAME.ubam.bam 
fi

################################################################################
[ -e $OUTDIR/$FILENAME.ubam.bam ] && rm $OUTDIR/*.ubam.bam.dummy
echo ">>>>> creating unmapped bam - FINISHED"
echo ">>>>> enddate "`date`
