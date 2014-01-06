#!/bin/bash -e

# Script to run bowtie v1 program.
# It takes comma-seprated list of files containing short sequence reads in fasta or fastq format and bowtie index files as input.
# It produces output files: read alignments in .bam format and other files.
# author: Fabian Buske
# date: August 2013

# QCVARIABLES,Resource temporarily unavailable
# RESULTFILENAME <DIR>/<TASK>/<SAMPLE>.$ASD.bam


echo ">>>>> read mapping with bowtie 1"
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
        -s | --rgsi )           shift; SAMPLEID=$1 ;; # read group prefix
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

for MODULE in $MODULE_BOWTIE; do module load $MODULE; done  # save way to load modules that itself load other modules
export PATH=$PATH_BOWTIE:$PATH
module list
echo "PATH=$PATH"
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)
PATH_PICARD=$(dirname $(which MarkDuplicates.jar))

echo "[NOTE] set java parameters"
JAVAPARAMS="-Xmx"$(python -c "print int($MEMORY_BOWTIE*0.8)")"g -Djava.io.tmpdir="$TMP" -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1" 
unset _JAVA_OPTIONS
echo "JAVAPARAMS "$JAVAPARAMS

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--JAVA        --\n" $(java -Xmx200m -version 2>&1)
[ -z "$(which java)" ] && echo "[ERROR] no java detected" && exit 1
echo -e "--bowtie      --\n "$(bowtie --version)
[ -z "$(which bowtie)" ] && echo "[ERROR] no bowtie detected" && exit 1
echo -e "--samtools    --\n "$(samtools 2>&1 | head -n 3 | tail -n-2)
[ -z "$(which samtools)" ] && echo "[ERROR] no samtools detected" && exit 1
echo -e "--R           --\n "$(R --version | head -n 3)
[ -z "$(which R)" ] && echo "[ERROR] no R detected" && exit 1
echo -e "--PICARD      --\n "$(java $JAVAPARAMS -jar $PATH_PICARD/MarkDuplicates.jar --version 2>&1)
[ ! -f $PATH_PICARD/MarkDuplicates.jar ] && echo "[ERROR] no picard detected" && exit 1
echo -e "--bedtools --\n "$(bedtools --version)
[ -z "$(which bedtools)" ] && echo "[ERROR] no bedtools detected" && exit 1
echo -e "--samstat     --\n "$(samstat -h | head -n 2 | tail -n1)
[ -z "$(which samstat)" ] && echo "[ERROR] no samstat detected" && exit 1
echo -e "--convert     --\n "$(convert -version | head -n 1)
[ -z "$(which convert)" ] && echo "[WARN] imagemagick convert not detected" 

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="parameters"

# check library variables are set
if [[ -z "$EXPID" || -z "$LIBRARY" || -z "$PLATFORM" ]]; then
    echo "[ERROR] library info not set (EXPID, LIBRARY, and PLATFORM): free text needed"
    exit 1;
else
    echo "[NOTE] EXPID $EXPID; LIBRARY $LIBRARY; PLATFORM $PLATFORM"
fi

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
    [ -e $OUTDIR/$SAMPLE.$ASD.bam ] && rm $OUTDIR/$SAMPLE.$ASD.bam
    [ -e $OUTDIR/$SAMPLE.$ASD.bam.stats ] && rm $OUTDIR/$SAMPLE.$ASD.bam.stats
    [ -e $OUTDIR/$SAMPLE.$ASD.bam.dupl ] && rm $OUTDIR/$SAMPLE.$ASD.bam.dupl
    [ -e $OUTDIR/$SAMPLE.$ALN.bam ] && rm $OUTDIR/$SAMPLE.$ALN.*bam
fi

#is paired ?                                                                                                      
if [ "$f" != "${f/%$READONE.$FASTQ/$READTWO.$FASTQ}" ] && [ -e ${f/%$READONE.$FASTQ/$READTWO.$FASTQ} ]; then
    PAIRED="1"
else
    PAIRED="0"
fi

#is ziped ?                                                                                                       
ZCAT="zcat"
if [[ ${f##*.} != "gz" ]]; then ZCAT="cat"; fi

# get encoding
if [ -z "$FASTQ_ENCODING" ]; then 
    echo "[NOTE] Detect fastq Phred encoding"
    FASTQ_ENCODING=$($ZCAT $f |  awk 'NR % 4 ==0' | python $NGSANE_BASE/tools/GuessFastqEncoding.py |  tail -n 1)
    echo "[NOTE] $FASTQ_ENCODING fastq format detected"
fi

if [[ "$FASTQ_ENCODING" == *Phred33* ]]; then
    FASTQ_PHRED="--phred33-quals"    
    FASTQ_PHRED_TRIM=" -phred33"
elif [[ "$FASTQ_ENCODING" == *Illumina* ]]; then
    FASTQ_PHRED="--phred64-quals"
    FASTQ_PHRED_TRIM=" -phred64"
elif [[ "$FASTQ_ENCODING" == *Solexa* ]]; then
    FASTQ_PHRED="--solexa1.3-quals"
    FASTQ_PHRED_TRIM=" -phred64"
else
    echo "[NOTE] cannot detect/don't understand fastq format: $FASTQ_ENCODING - using default"
fi

THISTMP=$TMP"/"$(whoami)"/"$(echo $OUTDIR/$SAMPLE | md5sum | cut -d' ' -f1)
[ -d $THISTMP ] && rm -r $THISTMP
mkdir -p $THISTMP

#readgroup
FULLSAMPLEID=$SAMPLEID"${n/%$READONE.$FASTQ/}"
RG="--sam-RG \"ID:$EXPID\" --sam-RG \"SM:$FULLSAMPLEID\" --sam-RG \"LB:$LIBRARY\" --sam-RG \"PL:$PLATFORM\""

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
CHECKPOINT="bowtie"

if [ $PAIRED == "0" ]; then 
    READS="$f"
    let FASTQREADS=`$ZCAT $f | wc -l | gawk '{print int($1/4)}' `
else 

    READS="-1 $f -2 ${f/%$READONE.$FASTQ/$READTWO.$FASTQ}"
    READ1=`$ZCAT $f | wc -l | gawk '{print int($1/4)}' `
    READ2=`$ZCAT ${f/%$READONE.$FASTQ/$READTWO.$FASTQ} | wc -l | gawk '{print int($1/4)}' `
    let FASTQREADS=$READ1+$READ2
fi

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

	# Unpaired
	if [ $PAIRED == "0" ]; then
        echo "[NOTE] SINGLE READS"

        RUN_COMMAND="$ZCAT $f | bowtie $RG $BOWTIEADDPARAM $FASTQ_PHRED --threads $CPU_BOWTIE --sam $BOWTIE_OPTIONS ${FASTA%.*} - $THISTMP/$SAMPLE.$ALN.sam"

	#Paired
    else
        echo "[NOTE] PAIRED READS"
        # clever use of named pipes to avoid fastq unzipping
        [ -e $THISTMP/$SAMPLE${READONE}_pipe ] && rm $THISTMP/$SAMPLE${READONE}_pipe
        [ -e $THISTMP/$SAMPLE${READTWO}_pipe ] && rm $THISTMP/$SAMPLE${READTWO}_pipe
        mkfifo $THISTMP/$SAMPLE${READONE}_pipe $THISTMP/$SAMPLE${READTWO}_pipe
        
        $ZCAT $f > $THISTMP/$SAMPLE${READONE}_pipe &
        $ZCAT ${f/%$READONE.$FASTQ/$READTWO.$FASTQ} > $THISTMP/$SAMPLE${READTWO}_pipe &
        
		RUN_COMMAND="bowtie $RG $BOWTIEADDPARAM $FASTQ_PHRED --threads $CPU_BOWTIE --sam $BOWTIE_OPTIONS ${FASTA%.*} -1 $THISTMP/$SAMPLE${READONE}_pipe -2 $THISTMP/$SAMPLE${READTWO}_pipe $THISTMP/$SAMPLE.$ALN.sam"

    fi
    echo $RUN_COMMAND && eval $RUN_COMMAND
    
    # bam file conversion                                                                         
    samtools view -@ $CPU_BOWTIE -Sbt $FASTA.fai $THISTMP/$SAMPLE.$ALN.sam > $OUTDIR/$SAMPLE.$ALN.bam
    
    # cleanup
    [ -e $THISTMP/$SAMPLE${READONE}_pipe ] && rm $THISTMP/$SAMPLE${READONE}_pipe
    [ -e $THISTMP/$SAMPLE${READTWO}_pipe ] && rm $THISTMP/$SAMPLE${READTWO}_pipe
    [ -e $THISTMP/$SAMPLE.$ALN.sam ] && rm $THISTMP/$SAMPLE.$ALN.sam

    # mark checkpoint
    if [ -f $OUTDIR/$SAMPLE.$ALN.bam ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
fi

################################################################################
CHECKPOINT="bam conversion and sorting"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

    # create bam files for discarded reads and remove fastq files    
            
    if [ "$PAIRED" = "1" ]; then
        # fix mates
        samtools sort -@ $CPU_BOWTIE -n $OUTDIR/$SAMPLE.$ALN.bam $OUTDIR/$SAMPLE.tmp
        samtools fixmate $OUTDIR/$SAMPLE.tmp.bam $OUTDIR/$SAMPLE.$ALN.bam
        [ -e $OUTDIR/$SAMPLE.tmp.bam ] && rm $OUTDIR/$SAMPLE.tmp.bam
    fi
    
    samtools sort -@ $CPU_BOWTIE $OUTDIR/$SAMPLE.$ALN.bam $OUTDIR/$SAMPLE.ash

    # mark checkpoint
    if [ -f $OUTDIR/$SAMPLE.ash.bam ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
    
    # cleanup
    [ -e $OUTDIR/$SAMPLE.$ALN.bam ] && rm $OUTDIR/$SAMPLE.$ALN.bam
fi

################################################################################
CHECKPOINT="mark duplicates"
# create bam files for discarded reads and remove fastq files
if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
   
    if [ ! -e $OUTDIR/metrices ]; then mkdir -p $OUTDIR/metrices ; fi
    java $JAVAPARAMS -jar $PATH_PICARD/MarkDuplicates.jar \
        INPUT=$OUTDIR/$SAMPLE.ash.bam \
        OUTPUT=$OUTDIR/$SAMPLE.$ASD.bam \
        METRICS_FILE=$OUTDIR/metrices/$SAMPLE.$ASD.bam.dupl \
        AS=true \
        VALIDATION_STRINGENCY=LENIENT \
        TMP_DIR=$THISTMP
    samtools index $OUTDIR/$SAMPLE.$ASD.bam

    # mark checkpoint
    if [ -f $OUTDIR/$SAMPLE.$ASD.bam ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
    
    # cleanup
    [ -e $OUTDIR/$SAMPLE.ash.bam ] && rm $OUTDIR/$SAMPLE.ash.bam
fi

################################################################################
CHECKPOINT="statistics"                                                                                                

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
    
    STATSOUT=$OUTDIR/$SAMPLE.$ASD.bam.stats
    samtools flagstat $OUTDIR/$SAMPLE.$ASD.bam > $STATSOUT
    
    if [ -n "$SEQREG" ]; then
        echo "#custom region" >> $STATSOUT
        echo $(samtools view -@ $CPU_BOWTIE -c -F 4 $OUTDIR/$SAMPLE.$ASD.bam $SEQREG )" total reads in region " >> $STATSOUT
        echo $(samtools view -@ $CPU_BOWTIE -c -f 3 $OUTDIR/$SAMPLE.$ASD.bam $SEQREG )" properly paired reads in region " >> $STATSOUT
    fi

    # mark checkpoint
    if [ -f $STATSOUT ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
fi

################################################################################
CHECKPOINT="calculate inner distance"                                                                                                

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
    
    java $JAVAPARAMS -jar $PATH_PICARD/CollectMultipleMetrics.jar \
        INPUT=$OUTDIR/$SAMPLE.$ASD.bam \
        REFERENCE_SEQUENCE=$FASTA \
        OUTPUT=$OUTDIR/metrices/$SAMPLE.$ASD.bam \
        VALIDATION_STRINGENCY=LENIENT \
        PROGRAM=CollectAlignmentSummaryMetrics \
        PROGRAM=CollectInsertSizeMetrics \
        PROGRAM=QualityScoreDistribution \
        TMP_DIR=$THISTMP
      
    # create pdfs
    if [ -n "$(which convert)" ]; then 
        for im in $( ls $OUTDIR/metrices/*.pdf ); do
            convert $im ${im/pdf/jpg}
        done
    fi

    # mark checkpoint
    [ -f $OUTDIR/metrices/$SAMPLE.$ASD.bam.alignment_summary_metrics ] && echo -e "\n********* $CHECKPOINT\n" && unset RECOVERFROM
fi


################################################################################
CHECKPOINT="samstat"    

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
    
    samstat $OUTDIR/$SAMPLE.$ASD.bam 2>&1 | tee | grep -v -P "Bad x in routine betai"

    # mark checkpoint
    if [ -f $OUTDIR/$SAMPLE.$ASD.bam.stats ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
fi


################################################################################
CHECKPOINT="verify"    
    
BAMREADS=`head -n1 $OUTDIR/$SAMPLE.$ASD.bam.stats | cut -d " " -f 1`
if [ "$BAMREADS" = "" ]; then let BAMREADS="0"; fi
if [ $BAMREADS -eq $FASTQREADS ]; then
    echo "[NOTE] PASS check mapping: $BAMREADS == $FASTQREADS"
else
    echo -e "[ERROR] We are loosing reads from .fastq -> .bam in $f: \nFastq had $FASTQREADS Bam has $BAMREADS"
    exit 1
fi

echo -e "\n********* $CHECKPOINT\n"
################################################################################
[ -e $OUTDIR/$SAMPLE.$ASD.bam.dummy ] && rm $OUTDIR/$SAMPLE.$ASD.bam.dummy
echo ">>>>> read mapping with bowtie 1 - FINISHED"
echo ">>>>> enddate "`date`
