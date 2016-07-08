#!/bin/bash -e

# Script for creating clean mapped bam files
# author: Tim Kahlke
# date: July 2016

# QCVARIABLES,Resource temporarily unavailable
# RESULTFILENAME <DIR>/<TASK>/<SAMPLE>.bam

echo ">>>>> Creating a clean mapped bam file"
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


# define defaults
CLIPPING_ACTION=2
NON_PF=true
CLIPPING_ATTRIBUTE=XT
CREATE_INDEX=true
ADD_MATE_CIGAR=true
CLIP_ADAPTERS=false
CLIP_OVERLAPPING_READS=true
INCLUDE_SECONDARY_ALIGNMENTS=true
MAX_INSERTIONS_OR_DELETIONS=-1
PRIMARY_ALIGNMENT_STRATEGY=MostDistant
ATTRIBUTES_TO_RETAIN=XS

#Overwrite defaults if set
. $CONFIG
. ${NGSANE_BASE}/conf/header.sh
. $CONFIG

################################################################################
NGSANE_CHECKPOINT_INIT "programs"

# save way to load modules that itself load other modules
hash module 2>/dev/null && for MODULE in $MODULE_GATK_MAPCLEAN; do module load $MODULE; done && module list

echo "[NOTE] set java parameters"
JAVAPARAMS="-Xmx"$(python -c "print int($MEMORY_GATK_MAPCLEAN*0.75)")"g -Djava.io.tmpdir="$TMP" -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1"
unset _JAVA_OPTIONS
echo "JAVAPARAMS "$JAVAPARAMS

echo "PATH=$PATH"
echo -e "--JAVA        --\n" $(java -Xmx200m -version 2>&1)
[ -z "$(which java)" ] && echo "[ERROR] no java detected" && exit 1

PATH_PICARD_JAR=$(which picard.jar)
[ ! -f $PATH_PICARD_JAR ] && echo "[ERROR] no picard detected" && exit 1

PATH_BWA=$(which bwa)
[ ! -f $PATH_BWA ] && echo "[ERROR] no bwa detected" && exit 1

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)


NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "parameters"

# get basename of input file f
IFILE=${INPUTFILE##*/}
FILENAME=${IFILE%.*}


if [ -z $REFERENCE ]; then
    echo "[ERROR] no reference provided (FASTA)"
    exit 1
fi

if [[ ! -e $REFERENCE.bwt ]]; then
    echo "[ERROR] BWA index not detected. Exeute bwaIndex.sh first"
    exit 1
fi

# delete old bam files unless attempting to recover
if [ -z "$NGSANE_RECOVERFROM" ]; then
    if [ -e $OUTDIR/$FILENAME.mapped_cleaned.bam ]; then
      rm $OUTDIR/$FILENAME.mapped_cleaned.bam
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
NGSANE_CHECKPOINT_INIT "Create adapter clipped fastq"


if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    java $JAVAPARAMS -jar $PATH_PICARD_JAR SamToFastq \
    I=$INPUTFILE \
    FASTQ=$OUTDIR/$FILENAME.clipped.fastq \
    CLIPPING_ATTRIBUTE=$CLIPPING_ATTRIBUTE \
    CLIPPING_ACTION=$CLIPPING_ACTION \
    INTERLEAVE=true \
    NON_PF=$NON_PF \
    TMP_DIR=$THISTMP


    NGSANE_CHECKPOINT_CHECK $OUTDIR/$FILENAME.clipped.fastq 
fi

################################################################################
NGSANE_CHECKPOINT_INIT "Align to reference using BWA"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    $PATH_BWA mem -M -t $CPU_GATK_MAPCLEAN  -p $REFERENCE \
    $OUTDIR/$FILENAME.clipped.fastq > $OUTDIR/$FILENAME.bwa_mem.sam

    NGSANE_CHECKPOINT_CHECK $OUTDIR/$FILENAME.bwa_mem.sam
fi

################################################################################
NGSANE_CHECKPOINT_INIT "create a mapped and cleaned bam file"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then
    
    java $JAVAPARAMS -jar $PATH_PICARD_JAR MergeBamAlignment \
    R=$REFERENCE \
    UNMAPPED_BAM=$INPUTFILE \
    ALIGNED_BAM=$OUTDIR/$FILENAME.bwa_mem.sam \
    O=$OUTDIR/$FILENAME.mapped_cleaned.bam \
    CREATE_INDEX=$CREATE_INDEX \
    ADD_MATE_CIGAR=$ADD_MATE_CIGAR \
    CLIP_ADAPTERS=$CLIP_ADAPTERS \
    CLIP_OVERLAPPING_READS=$CLIP_OVERLAPPING_READS \
    INCLUDE_SECONDARY_ALIGNMENTS=$INCLUDE_SECONDARY_ALIGNMENTS \
    MAX_INSERTIONS_OR_DELETIONS=$MAX_INSERTIONS_OR_DELETIONS \
    PRIMARY_ALIGNMENT_STRATEGY=$PRIMARY_ALIGNMENT_STRATEGY \
    ATTRIBUTES_TO_RETAIN=$ATTRIBUTES_TO_RETAIN\
    TMP_DIR=$THISTMP

    NGSANE_CHECKPOINT_CHECK $OUTDIR/$FILENAME.mapped_cleaned.bam
fi

################################################################################
[ -e $OUTDIR/$FILENAME.mapped_cleaned.bam ] && rm $OUTDIR/*bam.dummy
echo ">>>>> creating mapped and cleaned bam - FINISHED"
echo ">>>>> enddate "`date`
