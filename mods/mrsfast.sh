#!/bin/bash -e

# Script to run mrsfast program.
# It takes comma-seprated list of files containing short sequence reads in fasta or fastq format and mrsfast index files as input.
# It produces output files: read alignments in .bam format and other files.
# author: Fabian Buske
# date: August 2014

# QCVARIABLES,Resource temporarily unavailable
# RESULTFILENAME <DIR>/<TASK>/<SAMPLE>$ASD.bam


echo ">>>>> read mapping with mrsfast 3.x"
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
hash module 2>/dev/null && for MODULE in $MODULE_MRSFAST; do module load $MODULE; done && module list 

export PATH=$PATH_MRSFAST:$PATH
echo "PATH=$PATH"
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)
[ -z "$PATH_PICARD" ] && PATH_PICARD=$(dirname $(which MarkDuplicates.jar))

echo "[NOTE] set java parameters"
JAVAPARAMS="-Xmx"$(python -c "print int($MEMORY_MRSFAST*0.75)")"g -Djava.io.tmpdir="$TMP" -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1" 
unset _JAVA_OPTIONS
echo "JAVAPARAMS "$JAVAPARAMS

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--JAVA        --\n" $(java -Xmx200m -version 2>&1)
[ -z "$(which java)" ] && echo "[ERROR] no java detected" && exit 1
echo -e "--mrsfast     --\n "$(mrsfast --version)
[ -z "$(which mrsfast)" ] && echo "[ERROR] no mrsfast detected" && exit 1
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

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "parameters"

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

if [[ ! -e $FASTA.index ]]; then
    echo "[ERROR] MRSFAST index not detected. Exeute mrsfastIndex.sh first"
    exit 1
fi

# delete old bam files unless attempting to recover
if [ -z "$NGSANE_RECOVERFROM" ]; then
    [ -e $OUTDIR/$SAMPLE$ASD.bam ] && rm $OUTDIR/$SAMPLE$ASD.bam
    [ -e $OUTDIR/$SAMPLE$ASD.bam.stats ] && rm $OUTDIR/$SAMPLE$ASD.bam.stats
    [ -e $OUTDIR/$SAMPLE$ASD.bam.dupl ] && rm $OUTDIR/$SAMPLE$ASD.bam.dupl
    [ -e $OUTDIR/$SAMPLE$ALN.bam ] && rm $OUTDIR/$SAMPLE$ALN.*bam
fi

#is ziped ?
CAT="cat"
if [[ ${f##*.} == "gz" ]]; 
    then CAT="zcat";
elif [[ ${f##*.} == "bz2" ]]; 
    then CAT="bzcat"; 
fi

#is paired ?                                                                                                      
if [ "$f" != "${f/%$READONE.$FASTQ/$READTWO.$FASTQ}" ] && [ -e ${f/%$READONE.$FASTQ/$READTWO.$FASTQ} ]; then
    PAIRED="1"
else
    PAIRED="0"
fi

THISTMP=$TMP"/"$(whoami)"/"$(echo $OUTDIR/$SAMPLE | md5sum | cut -d' ' -f1)
[ -d $THISTMP ] && rm -r $THISTMP
mkdir -p $THISTMP

#readgroup
FULLSAMPLEID=$SAMPLEID"${n/%$READONE.$FASTQ/}"
RG="--sam-RG \"ID:$EXPID\" --sam-RG \"SM:$FULLSAMPLEID\" --sam-RG \"LB:$LIBRARY\" --sam-RG \"PL:$PLATFORM\""

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
NGSANE_CHECKPOINT_INIT "map with mrsfast"

if [ $PAIRED == "0" ]; then 
    READS="$f"
    let FASTQREADS=`$CAT $f | wc -l | gawk '{print int($1/4)}' `
else 

    READS="-1 $f -2 ${f/%$READONE.$FASTQ/$READTWO.$FASTQ}"
    READ1=`$CAT $f | wc -l | gawk '{print int($1/4)}' `
    READ2=`$CAT ${f/%$READONE.$FASTQ/$READTWO.$FASTQ} | wc -l | gawk '{print int($1/4)}' `
    let FASTQREADS=$READ1+$READ2
fi

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

	# Unpaired
	if [ $PAIRED == "0" ]; then
        echo "[NOTE] SINGLE READS"
        RUN_COMMAND="mrsfast $MRSFASTADDPARAM --mem $MEMORY_MRSFAST --search $FASTA --seq <($CAT $f) --threads $CPU_MRSFAST -o $THISTMP/$SAMPLE$ALN.sam"

	#Paired
    else
        echo "[NOTE] PAIRED READS"
        RUN_COMMAND="mrsfast $MRSFASTADDPARAM --mem $MEMORY_MRSFAST --pe --search $FASTA --seq1 <($CAT $f) --seq2 <($CAT ${f/%$READONE.$FASTQ/$READTWO.$FASTQ}) --threads $CPU_MRSFAST -o $THISTMP/$SAMPLE$ALN.sam"

    fi
         
    echo $RUN_COMMAND && eval $RUN_COMMAND

    samtools view -bS -t $FASTA.fai $THISTMP/$SAMPLE$ALN.sam > $OUTDIR/$SAMPLE$ALN.bam
    rm $THISTMP/$SAMPLE$ALN.sam

    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE$ALN.bam
fi

################################################################################
NGSANE_CHECKPOINT_INIT "sort bam"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    if [ "$PAIRED" = "1" ]; then
        # fix and sort
        echo "[NOTE] fixmate"
        RUN_COMMAND="java $JAVAPARAMS -jar $PATH_PICARD/FixMateInformation.jar \
            I=$OUTDIR/$SAMPLE$ALN.bam \
            O=$OUTDIR/$SAMPLE.ash.bam \
            VALIDATION_STRINGENCY=SILENT \
            SORT_ORDER=coordinate \
            TMP_DIR=$THISTMP"
        echo $RUN_COMMAND && eval $RUN_COMMAND
    else
        # just sort
        samtools sort -@ $CPU_MRSFAST $OUTDIR/$SAMPLE$ALN.bam $OUTDIR/$SAMPLE.ash
    fi

    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE.ash.bam
    
    # cleanup
    [ -e $OUTDIR/$SAMPLE$ALN.bam ] && rm $OUTDIR/$SAMPLE$ALN.bam
fi

################################################################################
NGSANE_CHECKPOINT_INIT "clean sam"
# create bam files for discarded reads and remove fastq files
if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then
   
    java $JAVAPARAMS -jar $PATH_PICARD/CleanSam.jar \
        INPUT=$OUTDIR/$SAMPLE.ash.bam \
        OUTPUT=$OUTDIR/$SAMPLE.cleaned.bam \
        VALIDATION_STRINGENCY=SILENT \
        TMP_DIR=$THISTMP

    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE.cleaned.bam
    
    # cleanup
    [ -e $OUTDIR/$SAMPLE.ash.bam ] && rm $OUTDIR/$SAMPLE.ash.bam
fi

################################################################################
NGSANE_CHECKPOINT_INIT "mark duplicates"
# create bam files for discarded reads and remove fastq files
if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then
   
    if [ ! -e $OUTDIR/metrices ]; then mkdir -p $OUTDIR/metrices ; fi
    java $JAVAPARAMS -jar $PATH_PICARD/MarkDuplicates.jar \
        INPUT=$OUTDIR/$SAMPLE.cleaned.bam \
        OUTPUT=$OUTDIR/$SAMPLE$ASD.bam \
        METRICS_FILE=$OUTDIR/metrices/$SAMPLE$ASD.bam.dupl \
        AS=true \
        CREATE_MD5_FILE=true \
        COMPRESSION_LEVEL=9 \
        VALIDATION_STRINGENCY=LENIENT \
        TMP_DIR=$THISTMP
    samtools index $OUTDIR/$SAMPLE$ASD.bam

    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE$ASD.bam
    
    # cleanup
    [ -e $OUTDIR/$SAMPLE.cleaned.bam ] && rm $OUTDIR/$SAMPLE.cleaned.bam
fi

################################################################################
NGSANE_CHECKPOINT_INIT "statistics"                                                                                                

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then
    
    STATSOUT=$OUTDIR/$SAMPLE$ASD.bam.stats
    samtools flagstat $OUTDIR/$SAMPLE$ASD.bam > $STATSOUT
    
    if [ -n "$SEQREG" ]; then
        echo "#custom region" >> $STATSOUT
        echo $(samtools view -@ $CPU_MRSFAST -c -F 4 $OUTDIR/$SAMPLE$ASD.bam $SEQREG )" total reads in region " >> $STATSOUT
        echo $(samtools view -@ $CPU_MRSFAST -c -f 3 $OUTDIR/$SAMPLE$ASD.bam $SEQREG )" properly paired reads in region " >> $STATSOUT
    fi

    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $STATSOUT
fi

################################################################################
NGSANE_CHECKPOINT_INIT "calculate inner distance"                                                                                                

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then
    
    java $JAVAPARAMS -jar $PATH_PICARD/CollectMultipleMetrics.jar \
        INPUT=$OUTDIR/$SAMPLE$ASD.bam \
        REFERENCE_SEQUENCE=$FASTA \
        OUTPUT=$OUTDIR/metrices/$SAMPLE$ASD.bam \
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
    NGSANE_CHECKPOINT_CHECK $OUTDIR/metrices/$SAMPLE$ASD.bam.alignment_summary_metrics
fi

################################################################################
NGSANE_CHECKPOINT_INIT "samstat"    

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then
    
    samstat $OUTDIR/$SAMPLE$ASD.bam 2>&1 | tee | grep -v -P "Bad x in routine betai"

    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE$ASD.bam.stats
fi

################################################################################
NGSANE_CHECKPOINT_INIT "verify"    
    
BAMREADS=`head -n1 $OUTDIR/$SAMPLE$ASD.bam.stats | cut -d " " -f 1`
if [ "$BAMREADS" = "" ]; then let BAMREADS="0"; fi
if [ $BAMREADS -eq $FASTQREADS ]; then
    echo "[NOTE] PASS check mapping: $BAMREADS == $FASTQREADS"
else
    echo -e "[ERROR] We are loosing reads from .fastq -> .bam in $f: \nFastq had $FASTQREADS Bam has $BAMREADS"
    exit 1
fi

NGSANE_CHECKPOINT_CHECK
################################################################################
[ -e $OUTDIR/$SAMPLE$ASD.bam.dummy ] && rm $OUTDIR/$SAMPLE$ASD.bam.dummy
echo ">>>>> read mapping with MRSFAST 1 - FINISHED"
echo ">>>>> enddate "`date`
