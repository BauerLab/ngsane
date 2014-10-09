#!/bin/bash -e

# BWA calling script
# author: Denis C. Bauer
# date: Nov.2010
# modified: August 2013 - Fabian Buske

# messages to look out for -- relevant for the QC.sh script:
# QCVARIABLES,We are loosing reads,MAPQ should be 0 for unmapped read,no such file,file not found,bwa.sh: line,Resource temporarily unavailable
# RESULTFILENAME <DIR>/<TASK>/<SAMPLE>$ASD.bam

echo ">>>>> readmapping with BWA "
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"


function usage {
echo -e "usage: $(basename $0) -k NGSANE -f FASTQ -r REFERENCE -o OUTDIR [OPTIONS]

Script running read mapping for single and paired DNA reads from fastq files
It expects a fastq file, pairdend, reference genome  as input and 
It runs BWA, converts the output to .bam files, adds header information.

required:
  -k | --toolkit <path>     location of the NGSANE repository 
  -f | --fastq <file>       fastq file
  -o | --outdir <path>      output dir

options:
  -s | --rgsi <name>        read group prefix (default: )
"
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
hash module 2>/dev/null && for MODULE in $MODULE_BWA; do module load $MODULE; done && module list 

export PATH=$PATH_BWA:$PATH
echo "PATH=$PATH"
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)
[ -z "$PATH_PICARD" ] && PATH_PICARD=$(dirname $(which MarkDuplicates.jar))

echo "[NOTE] set java parameters"
JAVAPARAMS="-Xmx"$(python -c "print int($MEMORY_BWA*0.8)")"g -Djava.io.tmpdir="$TMP"  -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1" 
unset _JAVA_OPTIONS
echo "JAVAPARAMS "$JAVAPARAMS

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--JAVA        --\n" $(java -Xmx200m -version 2>&1)
[ -z "$(which java)" ] && echo "[ERROR] no java detected" && exit 1
echo -e "--bwa         --\n "$(bwa 2>&1 | head -n 3 | tail -n-2)
[ -z "$(which bwa)" ] && echo "[ERROR] no bwa detected" && exit 1
echo -e "--samtools    --\n "$(samtools 2>&1 | head -n 3 | tail -n-2)
[ -z "$(which samtools)" ] && echo "[ERROR] no samtools detected" && exit 1
echo -e "--R           --\n "$(R --version | head -n 3)
[ -z "$(which R)" ] && echo "[ERROR] no R detected" && exit 1
echo -e "--PICARD      --\n "$(java $JAVAPARAMS -jar $PATH_PICARD/MarkDuplicates.jar --version 2>&1)
[ ! -f $PATH_PICARD/MarkDuplicates.jar ] && echo "[ERROR] no picard detected" && exit 1
echo -e "--samstat     --\n "$(samstat -h | head -n 2 | tail -n 1 )
[ -z "$(which samstat)" ] && echo "[ERROR] no samstat detected" && exit 1

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "parameters"

# get basename of f
n=${f##*/}
SAMPLE=${n/$READONE.$FASTQ/}

if [ -z "$FASTA" ]; then
    echo "[ERROR] no reference provided (FASTA)"
    exit 1
fi

if [[ ! -e $FASTA.bwt ]]; then
    echo "[ERROR] BWA index not detected. Exeute bwaIndex.sh first"
    exit 1
fi

# check library variables are set
if [[ -z "$EXPID" || -z "$LIBRARY" || -z "$PLATFORM" ]]; then
    echo "[ERROR] library info not set (EXPID, LIBRARY, and PLATFORM): free text needed"
    exit 1;
else
    echo "[NOTE] EXPID $EXPID; LIBRARY $LIBRARY; PLATFORM $PLATFORM"
fi


# delete old bam files unless attempting to recover
if [ -z "$NGSANE_RECOVERFROM" ]; then
    [ -e $OUTDIR/$SAMPLE$ASD.bam ] && rm $OUTDIR/$SAMPLE$ASD.bam
    [ -e $OUTDIR/$SAMPLE$ASD.bam.stats ] && rm $OUTDIR/$SAMPLE$ASD.bam.stats
    [ -e $OUTDIR/$SAMPLE$ASD.bam.dupl ] && rm $OUTDIR/$SAMPLE$ASD.bam.dupl
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
    READ1=$($CAT $f | wc -l | gawk '{print int($1/4)}')
    READ2=$($CAT ${f/%$READONE.$FASTQ/$READTWO.$FASTQ} | wc -l | gawk '{print int($1/4)}')
    let FASTQREADS=$READ1+$READ2
else
    PAIRED="0"
    READS="$f"
    let FASTQREADS=`$CAT $f | wc -l | gawk '{print int($1/4)}' `
fi

if [[ -z "$BWA_ALGORITHM" ||  "$BWA_ALGORITHM" == sam* ]]; then
    if [ "$PAIRED" == "1" ]; then
        echo "[NOTE] BWA algorithm is sampe"
        BWA_ALGORITHM=sampe
    else
        echo "[NOTE] BWA algorithm is samse"
        BWA_ALGORITHM=samse
    fi
elif [[ "$BWA_ALGORITHM" == "mem" ]]; then
    echo "[NOTE] BWA algorithm is mem"
elif [[ "$BWA_ALGORITHM" == "bwasw" ]]; then
    echo "[NOTE] BWA algorithm is mem"
else
    echo "[ERROR] specified BWA algorithm invalid"
    exit 1
fi
    
# get encoding unless specified
if [ -z "$FASTQ_ENCODING" ]; then 
    echo "[NOTE] Detect fastq Phred encoding"
    FASTQ_ENCODING=$($CAT $f |  awk 'NR % 4 ==0' | python $NGSANE_BASE/tools/GuessFastqEncoding.py |  tail -n 1)
    echo "[NOTE] $FASTQ_ENCODING fastq format detected"
fi
if [[ "$FASTQ_ENCODING" == **Phred33** ]]; then
    FASTQ_PHRED=""
    FASTQ_PHRED_TRIM=" -phred33"
elif [[ "$FASTQ_ENCODING" == **Phred64** ]]; then
    FASTQ_PHRED="-I"
    FASTQ_PHRED_TRIM=" -phred64"
else
    echo "[NOTE] cannot detect/don't understand fastq format: $FASTQ_ENCODING - using default"
fi

FULLSAMPLEID=$SAMPLEID$SAMPLE

THISTMP=$TMP"/"$(whoami)"/"$(echo $OUTDIR/$SAMPLE | md5sum | cut -d' ' -f1)
[ -d $THISTMP ] && rm -r $THISTMP
mkdir -p $THISTMP

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "recall files from tape"

if [ -n "$DMGET" ]; then
	dmget -a $FASTA*
    dmget -a ${f/$READONE/"*"}
	dmget -a ${OUTDIR}
fi
    
NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "BWA"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    if [ "$PAIRED" = 1 ]; then

        echo "[NOTE] PAIRED READS"
        
        if [[ "$BWA_ALGORITHM" == "sampe" ]]; then
            bwa aln $QUAL $BWAALNADDPARAM $FASTQ_PHRED -t $CPU_BWA $FASTA $f > $THISTMP/$SAMPLE$READONE.sai
            bwa aln $QUAL $BWAALNADDPARAM $FASTQ_PHRED -t $CPU_BWA $FASTA ${f/%$READONE.$FASTQ/$READTWO.$FASTQ} > $THISTMP/$SAMPLE$READTWO.sai
            bwa sampe $BWASAMPLEADDPARAM -r "@RG\tID:$EXPID\tSM:$FULLSAMPLEID\tPL:$PLATFORM\tLB:$LIBRARY" \
                $FASTA $THISTMP/$SAMPLE$READONE.sai $THISTMP/$SAMPLE$READTWO.sai \
                $f ${f/%$READONE.$FASTQ/$READTWO.$FASTQ} | samtools view -bS -t $FASTA.fai - > $OUTDIR/$SAMPLE$ALN.bam

            [ -f $THISTMP/$SAMPLE$READONE.sai ] && rm $THISTMP/$SAMPLE$READONE.sai
            [ -f $THISTMP/$SAMPLE$READTWO.sai ] && rm $THISTMP/$SAMPLE$READTWO.sai
            
        elif [[ "$BWA_ALGORITHM" == "mem" ]]; then
        
            bwa mem $BWASAMPLEADDPARAM -t $CPU_BWA -R "@RG\tID:$EXPID\tSM:$FULLSAMPLEID\tPL:$PLATFORM\tLB:$LIBRARY" \
    	       $FASTA $f ${f/%$READONE.$FASTQ/$READTWO.$FASTQ} | samtools view -bS -t $FASTA.fai - > $OUTDIR/$SAMPLE$ALN.bam

        elif [[ "$BWA_ALGORITHM" == "bwasw" ]]; then

            bwa bwasw $BWASAMPLEADDPARAM -t $CPU_BWA
    	       $FASTA $f ${f/%$READONE.$FASTQ/$READTWO.$FASTQ} | samtools view -bS -t $FASTA.fai - > $THISTMP/$SAMPLE$ALN.bam

            # add readgroup
            java $JAVAPARAMS -jar $PATH_PICARD/AddOrReplaceReadGroups.jar \
                INPUT=$THISTMP/$SAMPLE$ALN.bam \
                OUTPUT=$OUTDIR/$SAMPLE$ALN.bam \
                RGID=$EXPID RGLB=$LIBRARY RGPL=$PLATFORM \
                RGSM=$FULLSAMPLEID RGPU="XXXXXX"
                METRICS_FILE=$OUTDIR/metrices/$SAMPLE$ASD.bam.dupl \
                VALIDATION_STRINGENCY=SILENT \
                TMP_DIR=$THISTMP

        fi
    else
        echo "[NOTE] SINGLE READS"
        
        if [[ "$BWA_ALGORITHM" == "samse" ]]; then

            bwa aln $QUAL $BWAALNADDPARAM -t $CPU_BWA $FASTA $f > $THISTMP/$SAMPLE$READONE.sai
            bwa samse $BWASAMPLEADDPARAM -r "@RG\tID:$EXPID\tSM:$FULLSAMPLEID\tPL:$PLATFORM\tLB:$LIBRARY" \ \
                $FASTA $THISTMP/$SAMPLE$READONE.sai \
            	$f | samtools view -bS -t $FASTA.fai - > $OUTDIR/$SAMPLE$ALN.bam
    
            [ -f $THISTMP/$SAMPLE$READONE.sai ] && rm $THISTMP/$SAMPLE$READONE.sai
            
        elif [[ "$BWA_ALGORITHM" == "mem" ]]; then

            bwa mem $BWASAMPLEADDPARAM -t $CPU_BWA -R "@RG\tID:$EXPID\tSM:$FULLSAMPLEID\tPL:$PLATFORM\tLB:$LIBRARY" \
                $FASTA $f | samtools view -bS -t $FASTA.fai - > $OUTDIR/$SAMPLE$ALN.bam
        
        elif [[ "$BWA_ALGORITHM" == "bwasw" ]]; then
        
            bwa bwasw $BWASAMPLEADDPARAM -t $CPU_BWA $FASTA $f | samtools view -bS -t $FASTA.fai - > $THISTMP/$SAMPLE$ALN.bam

            # add readgroup
            java $JAVAPARAMS -jar $PATH_PICARD/AddOrReplaceReadGroups.jar \
                INPUT=$THISTMP/$SAMPLE$ALN.bam \
                OUTPUT=$OUTDIR/$SAMPLE$ALN.bam \
                RGID=$EXPID RGLB=$LIBRARY RGPL=$PLATFORM \
                RGSM=$FULLSAMPLEID RGPU="XXXXXX"
                METRICS_FILE=$OUTDIR/metrices/$SAMPLE$ASD.bam.dupl \
                VALIDATION_STRINGENCY=SILENT \
                TMP_DIR=$THISTMP
                
            [ -f $THISTMP/$SAMPLE$ALN.bam ] && rm $THISTMP/$SAMPLE$ALN.bam
        fi
   
    fi
    
    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE$ALN.bam
fi


################################################################################
NGSANE_CHECKPOINT_INIT "bam sorting"

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
        samtools sort -@ $CPU_BOWTIE $OUTDIR/$SAMPLE$ALN.bam $OUTDIR/$SAMPLE.ash
    fi

    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE.ash.bam 
    
    # cleanup
    [ -e $OUTDIR/$SAMPLE$ALN.bam ] && rm $OUTDIR/$SAMPLE$ALN.bam
fi

################################################################################
NGSANE_CHECKPOINT_INIT "clean sam"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then
   
    java $JAVAPARAMS -jar $PATH_PICARD/CleanSam.jar \
        INPUT=$OUTDIR/$SAMPLE.ash.bam \
        OUTPUT=$OUTDIR/$SAMPLE.cleaned.bam \
        VALIDATION_STRINGENCY=LENIENT \
        TMP_DIR=$THISTMP

    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE.cleaned.bam
    
    # cleanup
    [ -e $OUTDIR/$SAMPLE.ash.bam  ] && rm $OUTDIR/$SAMPLE.ash.bam 
fi

################################################################################
NGSANE_CHECKPOINT_INIT "mark duplicates"
# create bam files for discarded reads and remove fastq files
if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then
    
    #TODO look at samtools for rmdup
    #val string had to be set to LENIENT (SIlENT) to avoid crash due to a definition dis-
    #agreement between bwa and picard
    #http://seqanswers.com/forums/showthread.php?t=4246
    if [ ! -e $OUTDIR/metrices ]; then mkdir -p $OUTDIR/metrices ; fi
    java $JAVAPARAMS -jar $PATH_PICARD/MarkDuplicates.jar \
        INPUT=$OUTDIR/$SAMPLE.cleaned.bam \
        OUTPUT=$OUTDIR/$SAMPLE$ASD.bam \
        METRICS_FILE=$OUTDIR/metrices/$SAMPLE$ASD.bam.dupl \
        AS=true \
        CREATE_MD5_FILE=true \
        COMPRESSION_LEVEL=9 \
        VALIDATION_STRINGENCY=SILENT \
        TMP_DIR=$THISTMP

    samtools index $OUTDIR/$SAMPLE$ASD.bam

    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE$ASD.bam
    
    # cleanup 
    [ -f $OUTDIR/$SAMPLE.cleaned.bam ] && rm $OUTDIR/$SAMPLE.cleaned.bam
fi

################################################################################
NGSANE_CHECKPOINT_INIT "statistics"                                                                                                

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then
        
    samtools flagstat $OUTDIR/$SAMPLE$ASD.bam > $OUTDIR/$SAMPLE$ASD.bam.stats
    if [ -n "$SEQREG" ]; then
        echo "#custom region" >> $OUTDIR/$SAMPLE$ASD.bam.stats
        echo $(samtools view -c -F 4 $OUTDIR/$SAMPLE$ASD.bam $SEQREG )" total reads in region " >> $OUTDIR/$SAMPLE$ASD.bam.stats
        echo $(samtools view -c -f 3 $OUTDIR/$SAMPLE$ASD.bam $SEQREG )" properly paired reads in region " >> $OUTDIR/$SAMPLE$ASD.bam.stats
    fi

    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE$ASD.bam.stats
fi

################################################################################
NGSANE_CHECKPOINT_INIT "calculate inner distance"                                                                                                

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then
      
    export PATH=$PATH:/usr/bin/
    command="java $JAVAPARAMS -jar $PATH_PICARD/CollectMultipleMetrics.jar \
        INPUT=$OUTDIR/$SAMPLE$ASD.bam \
        REFERENCE_SEQUENCE=$FASTA \
        OUTPUT=$OUTDIR/metrices/$SAMPLE$ASD.bam \
        VALIDATION_STRINGENCY=SILENT \
        PROGRAM=CollectAlignmentSummaryMetrics \
        PROGRAM=CollectInsertSizeMetrics \
        PROGRAM=QualityScoreDistribution \
        TMP_DIR=$THISTMP"
        
    echo $command && eval $command
        
    for im in $( ls $OUTDIR/metrices/*.pdf ); do
        command="convert $im ${im/pdf/jpg}"
        echo $command && eval $command
    done

    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/metrices/$SAMPLE$ASD.bam.alignment_summary_metrics
fi

################################################################################
NGSANE_CHECKPOINT_INIT "samstat"    

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then
    
    samstat $OUTDIR/$SAMPLE$ASD.bam 2>&1 | tee | fgrep -v "Bad x in routine betai"

    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE$ASD.bam.stats
    
fi

################################################################################
NGSANE_CHECKPOINT_INIT "verify"    
    
BAMREADS=$(head -n1 $OUTDIR/$SAMPLE$ASD.bam.stats | cut -d " " -f 1)
if [ "$BAMREADS" = "" ]; then let BAMREADS="0"; fi			
if [ $BAMREADS -eq $FASTQREADS ]; then
    echo "[NOTE] PASS check mapping: $BAMREADS == $FASTQREADS"
    [ -e $OUTDIR/$SAMPLE.ash.bam ] && rm $OUTDIR/$SAMPLE.ash.bam
else
    echo -e "[ERROR] We are loosing reads from .fastq -> .bam in $f: \nFastq had $FASTQREADS Bam has $BAMREADS"
    exit 1 
fi

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "cleanup" 

[ -d $THISTMP ] && rm -r $THISTMP

NGSANE_CHECKPOINT_CHECK
################################################################################
[ -e $OUTDIR/$SAMPLE$ASD.bam.dummy ] && rm $OUTDIR/$SAMPLE$ASD.bam.dummy
echo ">>>>> readmapping with BWA - FINISHED"
echo ">>>>> enddate "`date`

