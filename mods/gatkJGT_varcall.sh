#!/bin/bash -e

# Script for creating g.vcf files from ubam 
# author: Tim Kahlke
# date: July 2016

# QCVARIABLES,Resource temporarily unavailable
# RESULTFILENAME <DIR>/<TASK>/<SAMPLE>.g.vcf

echo ">>>>> ubam to g.vcf using GATK"
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
        -f | --file )           shift; INPUTFILE=$1 ;;  # ubam file file 
        -o | --outdir )         shift; OUTDIR=$1 ;;     # output dir                                                     
        --recover-from )        shift; NGSANE_RECOVERFROM=$1 ;; # attempt to recover from log file
        -h | --help )           usage ;;
        * )                     echo "don't understand "$1
    esac
    shift
done



####DEFAULTS
# Mapclean
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

# markdups
MARKTOOL=MarkDuplicates
OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500
CREATE_INDEX=true

# gatkhaplo
STAND_EMIT_CONF=10
STAND_CALL_CONF=30



#PROGRAMS
. $CONFIG
. ${NGSANE_BASE}/conf/header.sh
. $CONFIG

################################################################################
NGSANE_CHECKPOINT_INIT "programs"

# save way to load modules that itself load other modules
hash module 2>/dev/null && for MODULE in $MODULE_GATKJGT_VARCALL; do module load $MODULE; done && module list

echo "[NOTE] set java parameters"
JAVAPARAMS="-Xmx"$(python -c "print int($MEMORY_GATKJGT_VARCALL*0.75)")"g -Djava.io.tmpdir="$TMP" -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1"
unset _JAVA_OPTIONS
echo "JAVAPARAMS "$JAVAPARAMS

PATH_PICARD_JAR=$(which picard.jar)
echo "PATH=$PATH"
echo -e "--JAVA        --\n" $(java -Xmx200m -version 2>&1)
[ -z "$(which java)" ] && echo "[ERROR] no java detected" && exit 1
[ ! -f $PATH_PICARD_JAR ] && echo "[ERROR] no picard detected" && exit 1


echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)

PATH_BWA=$(which bwa)
[ ! -f $PATH_BWA ] && echo "[ERROR] no bwa detected" && exit 1

PATH_GATK_JAR=$(which GenomeAnalysisTK.jar)
[ ! -f $PATH_GATK_JAR ] && echo "[ERROR] no gatk detected" && exit 1

NGSANE_CHECKPOINT_CHECK

################################################################################
NGSANE_CHECKPOINT_INIT "parameters"

# get basename of input file f
IFILE=${f##*/}
FILENAME=${IFILE%.*}

# delete old bam files unless attempting to recover
if [ -z "$NGSANE_RECOVERFROM" ]; then
    if [ -e $OUTDIR/$FILENAME.dedup.realigned.recalibrated.g.vcf ]; then
      rm $OUTDIR/$FILENAME.dedup.realigned.recalibrated.g.vcf
  fi
fi

if [ -z $REFERENCE ]; then
    echo "[ERROR] no reference provided (FASTA)"
    exit 1
fi

if [[ ! -e $REFERENCE.bwt ]]; then
    echo "[ERROR] BWA index not detected. Execute bwaIndex.sh first"
    exit 1
fi

RNAME=${REFERENCE%.*}

if [ ! -f $RNAME.dict ]; then
    echo "[ERRPR] No dictionary for reference found. Create a dicitonary first"
    exit 1
fi

if [ -z $INDREAL_KNOWN ]; then
    echo "[ERRPR] No known indel file given"
    exit 1
fi


KSTRING=""
IFS="," read -ra BQSR_KNOWNS <<<$BQSR_KNOWN
for k in $BQSR_KNOWNS; do
    KSTRING="$KSTRING -knownSites $k"
done



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
# MARKADAPT - Mark adapter sequences (soft clipping)
################################################################################
NGSANE_CHECKPOINT_INIT "marking adapters"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    java $JAVAPARAMS -jar $PATH_PICARD_JAR MarkIlluminaAdapters \
    I=$INPUTFILE \
    O=$OUTDIR/$FILENAME.clipped.bam \
    M=$THISTMP/markadapt.metrics \
    TMP_DIR=$THISTMP

    NGSANE_CHECKPOINT_CHECK $OUTDIR/$FILENAME.clipped.bam
fi

################################################################################
# MAPCLEAN - Map to reference and create "clean" bam file
################################################################################
NGSANE_CHECKPOINT_INIT "Create adapter clipped fastq"


if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    java $JAVAPARAMS -jar $PATH_PICARD_JAR SamToFastq \
    I=$OUTDIR/$FILENAME.clipped.bam \
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

    $PATH_BWA mem -M -t $CPU_MAPCLEAN  -p $REFERENCE \
    $OUTDIR/$FILENAME.clipped.fastq > $OUTDIR/$FILENAME.bwa_mem.sam

    NGSANE_CHECKPOINT_CHECK $OUTDIR/$FILENAME.bwa_mem.sam
fi

################################################################################
NGSANE_CHECKPOINT_INIT "create a mapped and cleaned bam file"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    java $JAVAPARAMS -jar $PATH_PICARD_JAR MergeBamAlignment \
    R=$REFERENCE \
    UNMAPPED_BAM=$OUTDIR/$FILENAME.clipped.bam \
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
    [ -e $OUTDIR/$FILENAME.clipped.bam ] && rm $OUTDIR/$FILENAME.clipped.bam
    [ -e $OUTDIR/$FILENAME.bwa_mem.sam ] && rm $OUTDIR/$FILENAME.bwa_mem.sam
fi
################################################################################
# MARKDUPS - Mark duplicate reads 
################################################################################
NGSANE_CHECKPOINT_INIT "marking duplicates"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    java $JAVAPARAMS -jar $PATH_PICARD_JAR $MARKTOOL \
    INPUT=$OUTDIR/$FILENAME.mapped_cleaned.bam \
    OUTPUT=$OUTDIR/$FILENAME.marked.bam \
    METRICS_FILE=$THISTMP/markadapt.metrics \
    OPTICAL_DUPLICATE_PIXEL_DISTANCE=$OPTICAL_DUPLICATE_PIXEL_DISTANCE \
    CREATE_INDEX=$CREATE_INDEX \
    TMP_DIR=$THISTMP

    NGSANE_CHECKPOINT_CHECK $OUTDIR/$FILENAME.marked.bam
    [ -e $OUTDIR/$FILENAME.mapped_cleaned.bam ] && rm $OUTDIR/$FILENAME.mapped_cleaned.bam
fi
################################################################################
# INDREAL - InDel re-alignment step
################################################################################
NGSANE_CHECKPOINT_INIT "create target intervals for indel realignment"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    java $JAVAPARAMS -jar $PATH_GATK_JAR -T RealignerTargetCreator \
    -R $REFERENCE \
    -known $INDREAL_KNOWN \
    -I $OUTDIR/$FILENAME.marked.bam \
    -o $OUTDIR/$FILENAME.realignertargets.intervals

    NGSANE_CHECKPOINT_CHECK $OUTDIR/$FILENAME.realignertargets.intervals
fi

################################################################################
NGSANE_CHECKPOINT_INIT "realigning indels"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    java $JAVAPARAMS -jar $PATH_GATK_JAR -T IndelRealigner \
    -R $REFERENCE \
    -targetIntervals $OUTDIR/$FILENAME.realignertargets.intervals \
    -known $INDREAL_KNOWN \
    -I $OUTDIR/$FILENAME.marked.bam \
    -o $OUTDIR/$FILENAME.ind_real.bam

    [ -e $OUTDIR/$FILENAME.realignertargets.intervals ] && rm $OUTDIR/$FILENAME.realignertargets.intervals
    [ -e $OUTDIR/$FILENAME.marked.bam ] && rm $OUTDIR/$FILENAME.marked.bam

    NGSANE_CHECKPOINT_CHECK $OUTDIR/$FILENAME.ind_real.bam
fi

################################################################################
# BQSR - Base Quality Score Re-alignment step
################################################################################
NGSANE_CHECKPOINT_INIT "creating recal table for base quality score recalibration"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    java $JAVAPARAMS -jar $PATH_GATK_JAR -T BaseRecalibrator \
    -R $REFERENCE \
    -I $OUTDIR/$FILENAME.ind_real.bam \
    -o $OUTDIR/$FILENAME.recal_data.table \
    $KSTRING

    NGSANE_CHECKPOINT_CHECK $OUTDIR/$FILENAME.recal_data.table
fi

################################################################################
NGSANE_CHECKPOINT_INIT "creating covariation table for base quality score recalibration"


if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    java $JAVAPARAMS -jar $PATH_GATK_JAR -T BaseRecalibrator \
    -R $REFERENCE \
    -I $OUTDIR/$FILENAME.ind_real.bam \
    -BQSR $OUTDIR/$FILENAME.recal_data.table \
    -o $OUTDIR/$FILENAME.post_recal_data.table \
    $KSTRING

    NGSANE_CHECKPOINT_CHECK $OUTDIR/$FILENAME.post_recal_data.table
fi

################################################################################
#NGSANE_CHECKPOINT_INIT "create before/after plots"
#
# taken out due to R error
#
#if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

#    java $JAVAPARAMS -jar $PATH_GATK_JAR -T AnalyzeCovariates \
#    -R $REFERENCE \
#    -before $OUTDIR/$FILENAME.recal_data.table \
#    -after $OUTDIR/$FILENAME.post_recal_data.table \
#    -plots $OUTDIR/$FILENAME.BQSR.recalibration_plots.pdf

#    NGSANE_CHECKPOINT_CHECK $OUTDIR/$FILENAME.BQSR.reclibration_plots.pdf
#    [ -e $OUTDIR/$FILENAME.post_recal_data.table ] && rm $OUTDIR/$FILENAME.post_recal_data.table
#fi

################################################################################
NGSANE_CHECKPOINT_INIT "apply recalibration "

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    java $JAVAPARAMS -jar $PATH_GATK_JAR -T PrintReads \
    -R $REFERENCE \
    -I $OUTDIR/$FILENAME.ind_real.bam \
    -o $OUTDIR/$FILENAME.bqs_recalibrated.bam \
    -BQSR $OUTDIR/$FILENAME.recal_data.table

    NGSANE_CHECKPOINT_CHECK $OUTDIR/$FILENAME.bqs_recalibrated.bam
    [ -e $OUTDIR/$FILENAME.ind_real.bam ] && rm $OUTDIR/$FILENAME.ind_real.bam
    [ -e $OUTDIR/$FILENAME.recal_data.table ] && rm $OUTDIR/$FILENAME.recal_data.table
fi
################################################################################
# GATKHAPLO - GATK's Haplotypecaller to create g.vcf 
################################################################################
NGSANE_CHECKPOINT_INIT "discover variants"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    java $JAVAPARAMS -jar $PATH_GATK_JAR -T HaplotypeCaller \
    -R $REFERENCE \
    -I $OUTDIR/$FILENAME.bqs_recalibrated.bam \
    --genotyping_mode DISCOVERY \
    -stand_emit_conf $STAND_EMIT_CONF \
    -stand_call_conf $STAND_CALL_CONF \
    -emitRefConfidence GVCF \
    -o $OUTDIR/$FILENAME.dedup.realigned.recalibrated.hc.g.vcf

    NGSANE_CHECKPOINT_CHECK $OUTDIR/$FILENAME.dedup.realigned.recalibrated.hc.g.vcf
fi
################################################################################
[ -e $OUTDIR/$FILENAME.ubam.bam ] && rm $OUTDIR/*.ubam.bam.dummy
echo ">>>>> creating unmapped bam - FINISHED"
echo ">>>>> enddate "`date`
