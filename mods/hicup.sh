#!/bin/bash -e

# Script running hicup including reference genome digestion, read mapping for single 
# and paired DNA reads with bowtie from fastq files
# It expects a fastq file, paired-end, reference genome and digest pattern as input.
# author: Fabian Buske
# date: Jan 2014

# messages to look out for -- relevant for the QC.sh script:
# QCVARIABLES,Resource temporarily unavailable
# RESULTFILENAME <DIR>/<TASK>/<SAMPLE>.fragmentLists.gz

echo ">>>>> HiC readmapping with HiCUP "
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"

function usage {
echo -e "usage: $(basename $0) -k NGSANE -f FASTQ -r REFERENCE -o OUTDIR [OPTIONS]"
exit
}

if [ ! $# -gt 3 ]; then usage ; fi

#INPUTS
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the NGSANE repository
        -f | --fastq )          shift; f=$1 ;; # fastq file
        -o | --outdir )         shift; OUTDIR=$1 ;; # output dir
        --recover-from )        shift; RECOVERFROM=$1 ;; # attempt to recover from log file
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
CHECKPOINT="programs"

for MODULE in $MODULE_HICUP; do module load $MODULE; done  # save way to load modules that itself load other modules
export PATH=$PATH_HICUP:$PATH
module list
echo "PATH=$PATH"
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--bowtie      --\n "$(bowtie --version | head -n 1 )
[ -z "$(which bowtie)" ] && echo "[ERROR] no bowtie detected" && exit 1
echo -e "--samtools    --\n "$(samtools 2>&1 | head -n 3 | tail -n-2)
[ -z "$(which samtools)" ] && echo "[ERROR] no samtools detected" && exit 1
echo -e "--perl        --\n "$(perl -v | grep "This is perl" )
[ -z "$(which perl)" ] && echo "[ERROR] no perl detected" && exit 1
echo -e "--HiCUP       --\n "$(hicup --version )
[ -z "$(which hicup)" ] && echo "[ERROR] no hicup detected" && exit 1

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="parameters"

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
if [ -z "$RECOVERFROM" ]; then
    [ -d $OUTDIR/$SAMPLE ] && rm -r $OUTDIR/$SAMPLE
fi

#is paired ?
if [ "$f" != "${f/%$READONE.$FASTQ/$READTWO.$FASTQ}" ] && [ -e ${f/%$READONE.$FASTQ/$READTWO.$FASTQ} ]; then
    PAIRED="1"
else
    echo "[ERROR] HiCUP requires paired fastq libraries" && exit 1
fi

#is ziped ?
ZCAT="zcat"
if [[ $f != *.gz ]]; then ZCAT="cat"; fi

# get encoding
if [ -z "$FASTQ_ENCODING" ]; then 
    echo "[NOTE] Detect fastq Phred encoding"
    FASTQ_ENCODING=$($ZCAT $f |  awk 'NR % 4 ==0' | python $NGSANE_BASE/tools/GuessFastqEncoding.py |  tail -n 1)
    echo "[NOTE] $FASTQ_ENCODING fastq format detected"
fi

if [[ "$FASTQ_ENCODING" == *Phred33* ]]; then
    FASTQ_PHRED="phred33-quals"    
elif [[ "$FASTQ_ENCODING" == *Illumina* ]]; then
    FASTQ_PHRED="phred64-quals"
elif [[ "$FASTQ_ENCODING" == *Solexa* ]]; then
    FASTQ_PHRED="solexa1.3-quals"
else
    echo "[NOTE] cannot detect/don't understand fastq format: $FASTQ_ENCODING - using default (phred33-quals)"
    FASTQ_PHRED="phred33-quals"
fi

THISTMP=$TMP"/"$(whoami)"/"$(echo $OUTDIR/$SAMPLE | md5sum | cut -d' ' -f1)
[ -d $THISTMP ] && rm -r $THISTMP
mkdir -p $THISTMP


if [ -z "$HICUP_RENZYME1" ] || [ "${HICUP_RENZYME1,,}" == "none" ] || [ -z "$HICUP_RCUTSITE1" ]; then
    echo "[ERROR] Restriction enzyme 1 not defined" && exit 1
else
    ENZYME1PARAM="-re1 $HICUP_RCUTSITE1"
fi
if [ -z "$HICUP_RENZYME2" ] || [ "${HICUP_RENZYME2,,}" == "none" ] || [ -z "$HICUP_RCUTSITE2" ]; then
    echo "[NOTE] Restriction enzyme 2 not defined"
    HICUP_RENZYME2="none"
else
    ENZYME2PARAM="-re1 $HICUP_RCUTSITE2"
fi

DIGESTGENOME="$OUT/common/$TASK_HICUP/Digest_${REFERENCE_NAME}_${HICUP_RENZYME1}_${HICUP_RENZYME2}.txt"
if [ ! -f $DIGESTGENOME ]; then
    echo "[ERROR] digested genome not found: $DIGESTGENOME"; 
    exit 1
fi

mkdir -p $OUTDIR/$SAMPLE


echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="recall files from tape"

if [ -n "$DMGET" ]; then
	dmget -a $(dirname $FASTA)/*
	dmget -a ${f/$READONE/"*"}
	dmget -a $OUTDIR/*
fi

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="truncate"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

    [ -f $OUTDIR/$SAMPLE/${SAMPLE}${READONE}_trunc.gz ] && rm $OUTDIR/$SAMPLE/${SAMPLE}${READONE}_trunc.gz
    [ -f $OUTDIR/$SAMPLE/${SAMPLE}${READTWO}_trunc.gz ] && rm $OUTDIR/$SAMPLE/${SAMPLE}${READTWO}_trunc.gz
    
    RUN_COMMAND="$(which perl) $(which hicup_truncater) -datestamp run -outdir $OUTDIR/$SAMPLE/ $ENZYME1PARAM $ENZYME2PARAM -threads $CPU_HICUP -zip $f ${f/%$READONE.$FASTQ/$READTWO.$FASTQ}"
    echo $RUN_COMMAND && eval $RUN_COMMAND
    
    mv $OUTDIR/$SAMPLE/hicup_truncater_summary_run.txt $OUTDIR/$SAMPLE"_truncater_summary.txt"
    
    # mark checkpoint
    if [ -f $OUTDIR/$SAMPLE/${SAMPLE}${READONE}_trunc.gz ] && [ -f $OUTDIR/$SAMPLE/${SAMPLE}${READTWO}_trunc.gz ] ;then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi

################################################################################
CHECKPOINT="map"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

    [ -f $OUTDIR/$SAMPLE/${SAMPLE}${READONE}_trunc_${SAMPLE}${READTWO}_trunc.pair.gz ] && rm $OUTDIR/$SAMPLE/${SAMPLE}${READONE}_trunc_${SAMPLE}${READTWO}_trunc.pair.gz
    
    RUN_COMMAND="$(which perl) $(which hicup_mapper) -datestamp run -bowtie $(which bowtie) -format $FASTQ_PHRED -outdir $OUTDIR/$SAMPLE/ -index ${FASTA%.*}  -threads $CPU_HICUP -zip $OUTDIR/$SAMPLE/${SAMPLE}${READONE}_trunc.gz $OUTDIR/$SAMPLE/${SAMPLE}${READTWO}_trunc.gz"
    echo $RUN_COMMAND && eval $RUN_COMMAND
    
    mv $OUTDIR/$SAMPLE/hicup_mapper_summary_run.txt $OUTDIR/$SAMPLE"_mapper_summary.txt"
        
    # mark checkpoint
    if [ -f $OUTDIR/$SAMPLE/${SAMPLE}${READONE}_trunc_${SAMPLE}${READTWO}_trunc.pair.gz ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi

################################################################################
CHECKPOINT="filter"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

    [ -f $OUTDIR/$SAMPLE/${SAMPLE}${READONE}_trunc_${SAMPLE}${READTWO}_trunc.bam ] && rm $OUTDIR/$SAMPLE/${SAMPLE}${READONE}_trunc_${SAMPLE}${READTWO}_trunc.bam
    
    RUN_COMMAND="$(which perl) $(which hicup_filter) -datestamp run -digest $DIGESTGENOME  -outdir $OUTDIR/$SAMPLE/ -threads $CPU_HICUP -zip $OUTDIR/$SAMPLE/${SAMPLE}${READONE}_trunc_${SAMPLE}${READTWO}_trunc.pair.gz"
    echo $RUN_COMMAND && eval $RUN_COMMAND
    
    mv $OUTDIR/$SAMPLE/hicup_filter_summary_run.txt $OUTDIR/$SAMPLE"_filter_summary.txt"
    
    # mark checkpoint
    if [ -f $OUTDIR/$SAMPLE/${SAMPLE}${READONE}_trunc_${SAMPLE}${READTWO}_trunc.bam ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi

################################################################################
CHECKPOINT="de-duplicate"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

    [ -f $OUTDIR/${SAMPLE}_uniques.bam ] && rm $OUTDIR/${SAMPLE}_uniques.bam
    
    RUN_COMMAND="$(which perl) $(which hicup_deduplicator) -datestamp run -pipeline_outdir $OUTDIR/$SAMPLE/ -outdir $OUTDIR/$SAMPLE/ -zip $OUTDIR/$SAMPLE/${SAMPLE}${READONE}_trunc_${SAMPLE}${READTWO}_trunc.bam"
    echo $RUN_COMMAND && eval $RUN_COMMAND

    mv $OUTDIR/$SAMPLE/hicup_deduplicator_summary_run.txt $OUTDIR/$SAMPLE"_deduplicator_summary.txt"
    
    ln -s $SAMPLE/uniques_${SAMPLE}${READONE}_trunc_${SAMPLE}${READTWO}_trunc.bam $OUTDIR/${SAMPLE}_uniques.bam
    
    # copy piecharts
    RUNSTATS=$OUT/runStats/$TASK_HICUP
    mkdir -p $RUNSTATS
    cp -f $OUTDIR/$SAMPLE/uniques_${SAMPLE}${READONE}_trunc_${SAMPLE}${READTWO}_trunc.bam_cis-trans.png $RUNSTATS/${SAMPLE}_uniques_cis-trans.png
    cp -f $OUTDIR/$SAMPLE/${SAMPLE}${READONE}_trunc_${SAMPLE}${READTWO}_trunc.pair.gz_ditag_classification.png $RUNSTATS/${SAMPLE}_ditag_classification.png
    cp -f $OUTDIR/$SAMPLE/${SAMPLE}${READONE}_trunc_${SAMPLE}${READTWO}_trunc.pair.gz_ditag_size_distribution.png $RUNSTATS/${SAMPLE}_ditag_size_distribution.png

    # mark checkpoint
    if [ -f $OUTDIR/${SAMPLE}_uniques.bam ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi

################################################################################
CHECKPOINT="count Interactions"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

    [ -f $OUTDIR/${SAMPLE}.fragmentLists.gz ] && rm $OUTDIR/${SAMPLE}.fragmentLists.gz
    [ -f $OUTDIR/${SAMPLE}.contactCounts.gz ] && rm $OUTDIR/${SAMPLE}.contactCounts.gz

    RUN_COMMAND="python ${NGSANE_BASE}/tools/hicupCountInteractions.py --verbose --genomeFragmentFile=$DIGESTGENOME --outputDir=$OUTDIR/ $OUTDIR/${SAMPLE}_uniques.bam"
    echo $RUN_COMMAND && eval $RUN_COMMAND

    [ -e $OUTDIR/${SAMPLE}_uniques.bam.fragmentLists ] && mv $OUTDIR/${SAMPLE}_uniques.bam.fragmentLists $OUTDIR/${SAMPLE}.fragmentLists
    [ -e $OUTDIR/${SAMPLE}_uniques.bam.contactCounts ] && mv $OUTDIR/${SAMPLE}_uniques.bam.contactCounts $OUTDIR/${SAMPLE}.contactCounts
    
    $GZIP $OUTDIR/${SAMPLE}.fragmentLists $OUTDIR/${SAMPLE}.contactCounts
    
    # mark checkpoint
    if [ -f $OUTDIR/${SAMPLE}.fragmentLists.gz ] && [ -f $OUTDIR/${SAMPLE}.contactCounts.gz ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi

################################################################################
[ -e $OUTDIR/${SAMPLE}.fragmentLists.gz.dummy ] && rm $OUTDIR/${SAMPLE}.fragmentLists.gz.dummy
echo ">>>>> readmapping with hicup (bowtie) - FINISHED"
echo ">>>>> enddate "`date`

