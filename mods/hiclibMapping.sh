#!/bin/bash -e

# Script running hiclib pipeline tapping into bowtie2
# It expects a fastq file, paired end, reference genome and digest pattern  as input.
# author: Fabian Buske
# date: April 2013

echo ">>>>> HiC analysis with hiclib "
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"

function usage {
echo -e "usage: $(basename $0) -k NGSANE -f FASTQ -r REFERENCE -o OUTDIR [OPTIONS]"
exit
}

# QCVARIABLES,Resource temporarily unavailable
# RESULTFILENAME <DIR>/<TASK>/<SAMPLE>-mapped_reads.hdf5

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
for MODULE in $MODULE_HICLIB; do module load $MODULE; done  # save way to load modules that itself load other modules

export PATH=$PATH_HICLIB:$PATH
module list
echo "PATH=$PATH"

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--Python      --\n" $(python --version)
[ -z "$(which python)" ] && echo "[ERROR] no python detected" && exit 1
echo -e "--Python libs --\n "$(yolk -l)
echo -e "--bowtie      --\n "$(bowtie --version)
[ -z "$(which bowtie)" ] && echo "[ERROR] no bowtie detected" && exit 1
echo -e "--samtools    --\n "$(samtools 2>&1 | head -n 3 | tail -n-2)
[ -z "$(which samtools)" ] && echo "[ERROR] no samtools detected" && exit 1
echo -e "--PICARD      --\n "$(java $JAVAPARAMS -jar $PATH_PICARD/MarkDuplicates.jar --version 2>&1)
[ ! -f $PATH_PICARD/MarkDuplicates.jar ] && echo "[ERROR] no picard detected" && exit 1
echo -e "--samstat     --\n "$(samstat -h | head -n 2 | tail -n1)
[ -z "$(which samstat)" ] && echo "[ERROR] no samstat detected" && exit 1

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="parameters"

# get basename of f
n=${f##*/}
SAMPLE=${n/$READONE.$FASTQ/}

if [ -z "$FASTA" ]; then
    echo "[ERROR] no reference provided (FASTA)"
    exit 1
fi

if [[ ! -e ${FASTA%.*}.1.bt2 ]]; then
    echo "[ERROR] Bowtie2 index not detected. Exeute bowtie2Index.sh first"
    exit 1
fi

#is paired ?                                                                                                      
if [ "$f" != "${f/%$READONE.$FASTQ/$READTWO.$FASTQ}" ] && [ -e ${f/%$READONE.$FASTQ/$READTWO.$FASTQ} ]; then
    PAIRED="1"
else
    echo "[ERROR] hiclib requires paired-end fastq files. Could not find ${f/%$READONE.$FASTQ/$READTWO.$FASTQ}" && exit 1
fi

if [ -z "$HICLIB_RENZYMES" ]; then
	echo "[ERROR] restriction enzyme not specified" && exit 1
fi

if [ -z "BOWTIE2INDEX" ]; then
	echo "[ERROR] bowtie2 index not specified" && exit 1
fi

THISTMP=$TMP"/"$(whoami)"/"$(echo $OUTDIR/$SAMPLE | md5sum | cut -d' ' -f1)
[ -d $THISTMP ] && rm -r $THISTMP
mkdir -p $THISTMP

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="recall files from tape"

if [ -n "$DMGET" ]; then
	dmget -a $(dirname $FASTA_CHROMDIR)/*
	dmget -a $(dirname $FASTA)/*
	dmget -a ${f/$READONE/"*"}
	dmget -a $OUTDIR/*
fi

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="run hiclib"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
    
    READS="$f ${f/%$READONE.$FASTQ/$READTWO.$FASTQ}"

    PARAMS="--restrictionEnzyme=$HICLIB_RENZYMES \
       --experimentName=$SAMPLE \
       --gapFile=$HICLIB_GAPFILE \
       --referenceGenome=$FASTA_CHROMDIR \
       --index=${FASTA%.*}"
    
    if [ "$FASTQ" = "sra" ]; then
    	PARAMS="$PARAMS --inputFormat=sra --sra-reader=$(which fastq-dump)"
    elif [ "$FASTQ" = "bam" ]; then
    	PARAMS="$PARAMS --inputFormat=bam"
    else
    	PARAMS="$PARAMS --inputFormat=fastq"
    fi
    
    if [ -n "$HICLIB_READLENGTH" ]; then
    	PARAMS="$PARAMS --readLength $HICLIB_READLENGTH"
    fi

    
    RUN_COMMAND="python ${NGSANE_BASE}/tools/hiclibMapping.py ${PARAMS} $HICLIBADDPARAM --bowtie=$(which bowtie2) --cpus=$CPU_HICLIB --outputDir=$OUTDIR --tmpDir=$THISTMP --verbose $READS &> $OUTDIR/$SAMPLE.log"
    echo $RUN_COMMAND && eval $RUN_COMMAND

    # mark checkpoint
    if [ -e $OUTDIR/${SAMPLE}-mapped_reads.hdf5 ] ;then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi
################################################################################
CHECKPOINT="treat read1.bam"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

    # treat read one    
    # merge
    samtools merge $THISTMP/$SAMPLE$READONE.bam $OUTDIR/*$READONE.bam.[0-9]*
    # keep mapped
    samtools view -bh -F 4 $THISTMP/$SAMPLE$READONE.bam > $THISTMP/$SAMPLE$READONE.$ALN.bam
    rm $THISTMP/$SAMPLE$READONE.bam 
    # sort
    samtools sort -@ $CPU_HICLIB $THISTMP/$SAMPLE$READONE.$ALN.bam $THISTMP/$SAMPLE$READONE.ash
    rm $THISTMP/$SAMPLE$READONE.$ALN.bam

    mkdir -p $OUTDIR/metrices
    java $JAVAPARAMS -jar $PATH_PICARD/MarkDuplicates.jar \
        INPUT=$THISTMP/$SAMPLE$READONE.ash.bam \
        OUTPUT=$OUTDIR/$SAMPLE$READONE.$ASD.bam \
        METRICS_FILE=$OUTDIR/$SAMPLE$READONE.$ASD.bam.dupl \
        AS=true \
        VALIDATION_STRINGENCY=LENIENT \
        TMP_DIR=$THISTMP
    samtools index $OUTDIR/$SAMPLE$READONE.$ASD.bam
    
    # mark checkpoint
    if [ -e $OUTDIR/$SAMPLE$READONE.$ASD.bam ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
    
    # cleanup
    [ -f $THISTMP/$SAMPLE$READONE.ash.bam ] && rm $THISTMP/$SAMPLE$READONE.ash.bam
    if [ -n "$HICLIB_KEEPBAM" ]; then
        rm $OUTDIR/$SAMPLE$READONE.bam.[0-9]*
    fi

fi

################################################################################
CHECKPOINT="treat read2.bam"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

    # treat read two
    # merge
    samtools merge $THISTMP/$SAMPLE$READTWO.bam $OUTDIR/*$READTWO.bam.[0-9]*
    # keep mapped
    samtools view -bh -F 4 $THISTMP/$SAMPLE$READTWO.bam > $THISTMP/$SAMPLE$READTWO.$ALN.bam
    rm $THISTMP/$SAMPLE$READTWO.bam
    # sort
    samtools sort -@ $CPU_HICLIB $THISTMP/$SAMPLE$READTWO.$ALN.bam $THISTMP/$SAMPLE$READTWO.ash
    rm $THISTMP/$SAMPLE$READTWO.$ALN.bam

    java $JAVAPARAMS -jar $PATH_PICARD/MarkDuplicates.jar \
        INPUT=$THISTMP/$SAMPLE$READTWO.ash.bam \
        OUTPUT=$OUTDIR/$SAMPLE$READTWO.$ASD.bam \
        METRICS_FILE=$OUTDIR/$SAMPLE$READTWO.$ASD.bam.dupl \
        AS=true \
        VALIDATION_STRINGENCY=LENIENT \
        TMP_DIR=$THISTMP
    samtools index $OUTDIR/$SAMPLE$READTWO.$ASD.bam

    # mark checkpoint
    if [ -e $OUTDIR/$SAMPLE$READTWO.$ASD.bam ] ;then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
    
    # cleanup
    [ -f $THISTMP/$SAMPLE$READTWO.ash.bam ] && rm $THISTMP/$SAMPLE$READTWO.ash.bam
    if [ -n "$HICLIB_KEEPBAM" ]; then
        $OUTDIR/$SAMPLE$READTWO.bam.[0-9]*
    fi

fi
################################################################################
CHECKPOINT="statistics"                                                                                                

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
    
    READ1STATSOUT=$OUTDIR/$OUTDIR/$SAMPLE$READONE.$ASD.bam.stats
    samtools flagstat $OUTDIR/$OUTDIR/$SAMPLE$READONE.$ASD.bam > $READ1STATSOUT
    
    if [ -n "$SEQREG" ]; then
        echo "#custom region" >> $READ1STATSOUT
        echo $(samtools view -@ $CPU_BOWTIE -c -F 4 $OUTDIR/$SAMPLE$READONE.$ASD.bam $SEQREG )" total reads in region " >> $READ1STATSOUT
        echo $(samtools view -@ $CPU_BOWTIE -c -f 3 $OUTDIR/$SAMPLE$READONE.$ASD.bam $SEQREG )" properly paired reads in region " >> $READ1STATSOUT
    fi
    
    READ2STATSOUT=$OUTDIR/$OUTDIR/$SAMPLE$READTWO.$ASD.bam.stats
    samtools flagstat $OUTDIR/$OUTDIR/$SAMPLE$READTWO.$ASD.bam > $READ2STATSOUT
    
    if [ -n "$SEQREG" ]; then
        echo "#custom region" >> $READ2STATSOUT
        echo $(samtools view -@ $CPU_BOWTIE -c -F 4 $OUTDIR/$SAMPLE$READTWO.$ASD.bam $SEQREG )" total reads in region " >> $READ2STATSOUT
        echo $(samtools view -@ $CPU_BOWTIE -c -f 3 $OUTDIR/$SAMPLE$READTWO.$ASD.bam $SEQREG )" properly paired reads in region " >> $READ2STATSOUT
    fi

    # mark checkpoint
    if [ -f $READ1STATSOUT ] && [ -f $READ2STATSOUT ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
fi
################################################################################
CHECKPOINT="samstat"    

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
    
    samstat $OUTDIR/$SAMPLE$READONE.$ASD.bam 2>&1 | tee | grep -v -P "Bad x in routine betai"
    samstat $OUTDIR/$SAMPLE$READTWO.$ASD.bam 2>&1 | tee | grep -v -P "Bad x in routine betai"

    # mark checkpoint
    if [ -f $OUTDIR/$SAMPLE$READONE.$ASD.bam.stats ] && [ -f $OUTDIR/$SAMPLE$READTWO.$ASD.bam.stats ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
fi
################################################################################
#CHECKPOINT="verify"    
#    
#BAMREAD1=`head -n1 $OUTDIR/$SAMPLE$READONE.$ASD.bam.stats | cut -d " " -f 1`
#BAMREAD2=`head -n1 $OUTDIR/$SAMPLE$READTWO.$ASD.bam.stats | cut -d " " -f 1`
#if [ "$BAMREAD1" = "" ]; then let BAMREAD1="0"; fi
#if [ "$BAMREAD2" = "" ]; then let BAMREAD2="0"; fi
#FASTQREAD1=`$ZCAT $f | wc -l | gawk '{print int($1/4)}' `
#FASTQREAD2=`$ZCAT ${f/%$READONE.$FASTQ/$READTWO.$FASTQ} | wc -l | gawk '{print int($1/4)}' `
#
#if [ $BAMREAD1 -eq $FASTQREAD1 ] && [ $BAMREAD2 -eq $FASTQREAD2 ]; then
#    echo "[NOTE] PASS check mapping: $BAMREAD1 == $FASTQREAD1 and $BAMREAD2 == $FASTQREAD2"
#else
#    echo -e "[ERROR] We are loosing reads from .fastq -> .bam in $f: \nFastqs had $FASTQREAD1 and $FASTQREAD2 Bams has $BAMREAD1 and $BAMREAD2"
#    exit 1
#fi
#
#echo -e "\n********* $CHECKPOINT\n"
################################################################################
[ -e $OUTDIR/${SAMPLE}-mapped_reads.hdf5.dummy ] && rm $OUTDIR/${SAMPLE}-mapped_reads.hdf5.dummy
echo ">>>>> readmapping with hiclib (Bowtie2) - FINISHED"
echo ">>>>> enddate "`date`
