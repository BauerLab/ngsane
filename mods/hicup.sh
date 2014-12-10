#!/bin/bash -e

# Script running hicup including reference genome digestion, read mapping for single 
# and paired DNA reads with bowtie from fastq files
# It expects a fastq file, paired-end, reference genome and digest pattern as input.
# author: Fabian Buske
# date: Jan 2014

# messages to look out for -- relevant for the QC.sh script:
# QCVARIABLES,Resource temporarily unavailable
# RESULTFILENAME <DIR>/<TASK>/<SAMPLE>$ASD.bam

echo ">>>>> HiC readmapping with HiCUP "
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"

function usage {
echo -e "usage: $(basename $0) -k NGSANE -f FASTQ -o OUTDIR [OPTIONS]"
exit
}

if [ ! $# -gt 3 ]; then usage ; fi

#INPUTS
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the NGSANE repository
        -f | --fastq )          shift; f=$1 ;; # fastq file
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
hash module 2>/dev/null && for MODULE in $MODULE_HICUP; do module load $MODULE; done && module list 

export PATH=$PATH_HICUP:$PATH
echo "PATH=$PATH"
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)
[ -z "$PATH_PICARD" ] && PATH_PICARD=$(dirname $(which MarkDuplicates.jar))

echo "[NOTE] set java parameters"
JAVAPARAMS="-Xmx"$(python -c "print int($MEMORY_HICUP*0.75)")"g -Djava.io.tmpdir="$TMP" -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1" 
unset _JAVA_OPTIONS
echo "JAVAPARAMS "$JAVAPARAMS

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--bowtie      --\n "$(bowtie --version | head -n 1 )
[ -z "$(which bowtie)" ] && echo "[ERROR] no bowtie detected" && exit 1
echo -e "--samtools    --\n "$(samtools 2>&1 | head -n 3 | tail -n-2)
[ -z "$(which samtools)" ] && echo "[ERROR] no samtools detected" && exit 1
echo -e "--perl        --\n "$(perl -v | grep "This is perl" )
[ -z "$(which perl)" ] && echo "[ERROR] no perl detected" && exit 1
echo -e "--HiCUP       --\n "$(hicup --version )
[ -z "$(which hicup)" ] && echo "[ERROR] no hicup detected" && exit 1
echo -e "--PICARD      --\n "$(java $JAVAPARAMS -jar $PATH_PICARD/MarkDuplicates.jar --version 2>&1)
[ ! -f $PATH_PICARD/MarkDuplicates.jar ] && echo "[ERROR] no picard detected" && exit 1
echo -e "--R           --\n "$(R --version | head -n 3)
[ -z "$(which R)" ] && echo "[ERROR] no R detected" && exit 1

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "parameters"

# get basename of f
n=${f##*/}
SAMPLE=${n/%$READONE.$FASTQ/}

if [ -z "$FASTA" ]; then
    echo "[ERROR] no reference provided (FASTA)"
    exit 1
fi

if [[ ! -e ${FASTA%.*}.1.ebwt ]]; then
    echo "[ERROR] Bowtie index not detected. Exeute bowtieIndex.sh first"
    exit 1
fi

# delete old bam files unless attempting to recover
if [ -z "$NGSANE_RECOVERFROM" ]; then
    [ -d $OUTDIR/$SAMPLE ] && rm -r $OUTDIR/$SAMPLE
    [ -f $OUTDIR/${SAMPLE}_uniques.bam ] && rm $OUTDIR/${SAMPLE}_uniques.bam   
fi

#is paired ?
if [ "$f" != "${f/%$READONE.$FASTQ/$READTWO.$FASTQ}" ] && [ -e ${f/%$READONE.$FASTQ/$READTWO.$FASTQ} ]; then
    PAIRED="1"
else
    echo "[ERROR] HiCUP requires paired fastq libraries" && exit 1
fi

#is ziped ?
CAT="cat"
if [[ ${f##*.} == "gz" ]]; 
    then CAT="zcat"; 
elif [[ ${f##*.} == "bz2" ]]; 
    then CAT="bzcat"; 
fi

if [ -z "$HICUP_RENZYME1" ] || [ -z "$HICUP_RCUTSITE1" ]; then
    echo "[ERROR] Restriction enzyme 1 not defined" && exit 1
else
    ENZYME1PARAM="--re1 $HICUP_RCUTSITE1"
fi
if [ -z "$HICUP_RENZYME2" ] || [ -z "$HICUP_RCUTSITE2" ]; then
    echo "[NOTE] Restriction enzyme 2 not defined"
    HICUP_RENZYME2="None"
else
    ENZYME2PARAM="--re2 $HICUP_RCUTSITE2"
fi

DIGESTGENOME="$OUT/common/$TASK_HICUP/Digest_${REFERENCE_NAME}_${HICUP_RENZYME1}_${HICUP_RENZYME2}.txt"
if [ ! -f $DIGESTGENOME ]; then
    echo "[ERROR] digested genome not found: $DIGESTGENOME"; 
    exit 1
fi


# get encoding
if [ -z "$FASTQ_ENCODING" ]; then 
    echo "[NOTE] Detect fastq Phred encoding"
    FASTQ_ENCODING=$($CAT $f |  awk 'NR % 4 ==0' | python $NGSANE_BASE/tools/GuessFastqEncoding.py |  tail -n 1)
    echo "[NOTE] $FASTQ_ENCODING fastq format detected"
fi

THISTMP=$TMP"/"$(whoami)"/"$(echo $OUTDIR/$SAMPLE | md5sum | cut -d' ' -f1)
[ -d $THISTMP ] && rm -r $THISTMP
mkdir -p $THISTMP

mkdir -p $OUTDIR/$SAMPLE


NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "recall files from tape"

if [ -n "$DMGET" ]; then
	dmget -a $(dirname $FASTA)/*
	dmget -a ${f/$READONE/"*"}
	dmget -a $OUTDIR/*
fi

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "truncate"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    RUN_COMMAND="$(which perl) $(which hicup_truncater) --datestamp run --outdir $OUTDIR/$SAMPLE/ $ENZYME1PARAM $ENZYME2PARAM --threads $CPU_HICUP --zip $f ${f/%$READONE.$FASTQ/$READTWO.$FASTQ}"
    echo $RUN_COMMAND && eval $RUN_COMMAND
    
    mv $OUTDIR/$SAMPLE/hicup_truncater_summary_run.txt $OUTDIR/$SAMPLE"_truncater_summary.txt"
    
    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE/${SAMPLE}${READONE}_trunc.gz $OUTDIR/$SAMPLE/${SAMPLE}${READTWO}_trunc.gz

fi

################################################################################
NGSANE_CHECKPOINT_INIT "map"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    [ -f $OUTDIR/$SAMPLE/${SAMPLE}${READONE}_trunc_${SAMPLE}${READTWO}_trunc.pair.gz ] && rm $OUTDIR/$SAMPLE/${SAMPLE}${READONE}_trunc_${SAMPLE}${READTWO}_trunc.pair.gz
    
    RUN_COMMAND="$(which perl) $(which hicup_mapper) $HICUPMAPPER_ADDPARAM --datestamp run --bowtie $(which bowtie) --outdir $OUTDIR/$SAMPLE/ --index ${FASTA%.*} --threads $CPU_HICUP --zip $OUTDIR/$SAMPLE/${SAMPLE}${READONE}_trunc.gz $OUTDIR/$SAMPLE/${SAMPLE}${READTWO}_trunc.gz"
    echo $RUN_COMMAND && eval $RUN_COMMAND
    
    mv $OUTDIR/$SAMPLE/hicup_mapper_summary_run.txt $OUTDIR/$SAMPLE"_mapper_summary.txt"

    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE/${SAMPLE}${READONE}_trunc_${SAMPLE}${READTWO}_trunc.pair.gz

fi

################################################################################
NGSANE_CHECKPOINT_INIT "filter"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    [ -f $OUTDIR/$SAMPLE/${SAMPLE}${READONE}_trunc_${SAMPLE}${READTWO}_trunc.bam ] && rm $OUTDIR/$SAMPLE/${SAMPLE}${READONE}_trunc_${SAMPLE}${READTWO}_trunc.bam
    
    RUN_COMMAND="$(which perl) $(which hicup_filter) --datestamp run --digest $DIGESTGENOME --outdir $OUTDIR/$SAMPLE/ --threads $CPU_HICUP --zip $OUTDIR/$SAMPLE/${SAMPLE}${READONE}_trunc_${SAMPLE}${READTWO}_trunc.pair.gz"
    echo $RUN_COMMAND && eval $RUN_COMMAND
    
    mv $OUTDIR/$SAMPLE/hicup_filter_summary_run.txt $OUTDIR/$SAMPLE"_filter_summary.txt"
    
    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE/${SAMPLE}${READONE}_trunc_${SAMPLE}${READTWO}_trunc.pair.gz.bam

fi

################################################################################
NGSANE_CHECKPOINT_INIT "de-duplicate"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then
    
    RUN_COMMAND="$(which perl) $(which hicup_deduplicator) --datestamp run --pipeline_outdir $OUTDIR/$SAMPLE/ --outdir $OUTDIR/$SAMPLE/ --threads $CPU_HICUP --zip $OUTDIR/$SAMPLE/${SAMPLE}${READONE}_trunc_${SAMPLE}${READTWO}_trunc.pair.gz.bam"
    echo $RUN_COMMAND && eval $RUN_COMMAND

    mv $OUTDIR/$SAMPLE/hicup_deduplicator_summary_run.txt $OUTDIR/$SAMPLE"_deduplicator_summary.txt"
    
    # move charts
    mv $OUTDIR/$SAMPLE/${SAMPLE}$READONE.$FASTQ.truncation_barchart.svg $OUTDIR/${SAMPLE}$READONE.truncation_barchart.svg
    mv $OUTDIR/$SAMPLE/${SAMPLE}$READTWO.$FASTQ.truncation_barchart.svg $OUTDIR/${SAMPLE}$READTWO.truncation_barchart.svg
        
    mv $OUTDIR/$SAMPLE/${SAMPLE}${READONE}_trunc_${SAMPLE}${READTWO}_trunc.pair.gz.bam.deduplicator_cis_trans_piechart.svg $OUTDIR/${SAMPLE}.cis-trans.svg
    mv $OUTDIR/$SAMPLE/${SAMPLE}${READONE}_trunc_${SAMPLE}${READTWO}_trunc.pair.gz.bam.deduplicator_uniques_barchart.svg $OUTDIR/${SAMPLE}.uniques_barchart.svg
    mv $OUTDIR/$SAMPLE/${SAMPLE}${READONE}_trunc_${SAMPLE}${READTWO}_trunc.pair.gz.ditag_size_distribution.svg $OUTDIR/${SAMPLE}.ditag_size_distribution.svg
    mv $OUTDIR/$SAMPLE/${SAMPLE}${READONE}_trunc_${SAMPLE}${READTWO}_trunc.pair.gz.filter_piechart.svg $OUTDIR/${SAMPLE}.filter_piechart.svg
    
    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE/uniques_${SAMPLE}${READONE}_trunc_${SAMPLE}${READTWO}_trunc.pair.gz.bam

fi

################################################################################
NGSANE_CHECKPOINT_INIT "bam sorting"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    samtools sort -n -l 9 -O bam -@ $CPU_HICUP -o $OUTDIR/$SAMPLE$ASD.bam -T $THISTMP/$SAMPLE $OUTDIR/$SAMPLE/uniques_${SAMPLE}${READONE}_trunc_${SAMPLE}${READTWO}_trunc.pair.gz.bam 

    samtools index $OUTDIR/$SAMPLE$ASD.bam
    
    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE$ASD.bam 
  
fi
################################################################################
NGSANE_CHECKPOINT_INIT "cleanup"

[ -f $OUTDIR/$SAMPLE/${SAMPLE}${READONE}_trunc_${SAMPLE}${READTWO}_trunc.pair.gz.bam ] && rm $OUTDIR/$SAMPLE/${SAMPLE}${READONE}_trunc_${SAMPLE}${READTWO}_trunc.pair.gz*
[ -f $OUTDIR/$SAMPLE/${SAMPLE}${READONE}_trunc_${SAMPLE}${READTWO}_trunc.gz ] && rm $OUTDIR/$SAMPLE/${SAMPLE}${READONE}_trunc_${SAMPLE}${READTWO}_trunc.gz*
[ -f $OUTDIR/$SAMPLE/${SAMPLE}${READONE}_trunc_${SAMPLE}${READTWO}_trunc.map ] && rm $OUTDIR/$SAMPLE/${SAMPLE}${READONE}_trunc_${SAMPLE}${READTWO}_trunc.map

NGSANE_CHECKPOINT_CHECK
################################################################################
[ -e $OUTDIR/$SAMPLE$ASD.bam.dummy ] && rm $OUTDIR/$SAMPLE$ASD.bam.dummy
echo ">>>>> readmapping with hicup (bowtie) - FINISHED"
echo ">>>>> enddate "`date`

