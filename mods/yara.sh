#!/bin/bash -e

# Script to run yara.
# It takes comma-seprated list of files containing short sequence reads in fasta or fastq format and masai index files as input.
# It produces output files: read alignments in .bam format and other files.
# author: Fabian Buske
# date: July 2014

# QCVARIABLES,Resource temporarily unavailable
# RESULTFILENAME <DIR>/<TASK>/<SAMPLE>$ASD.bam

echo ">>>>> readmapping with Yara "
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

FORCESINGLE=0

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
hash module 2>/dev/null && for MODULE in $MODULE_YARA; do module load $MODULE; done && module list 

export PATH=$PATH_YARA:$PATH
echo "PATH=$PATH"
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)
[ -z "$PATH_IGVTOOLS" ] && PATH_IGVTOOLS=$(dirname $(which igvtools.jar))
[ -z "$PATH_PICARD" ] && PATH_PICARD=$(dirname $(which MarkDuplicates.jar))

echo "[NOTE] set java parameters"
JAVAPARAMS="-Xmx"$(python -c "print int($MEMORY_YARA*0.8)")"g -Djava.io.tmpdir="$TMP"  -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1" 
unset _JAVA_OPTIONS
echo "JAVAPARAMS "$JAVAPARAMS

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--JAVA        --\n" $(java -Xmx200m -version 2>&1)
[ -z "$(which java)" ] && echo "[ERROR] no java detected" && exit 1
echo -e "--yara       --\n "$(yara_mapper --version)
[ -z "$(which yara_mapper)" ] && echo "[ERROR] no yara detected" && exit 1
echo -e "--samtools    --\n "$(samtools 2>&1 | head -n 3 | tail -n-2)
[ -z "$(which samtools)" ] && echo "[ERROR] no samtools detected" && exit 1
echo -e "--R           --\n "$(R --version | head -n 3)
[ -z "$(which R)" ] && echo "[ERROR] no R detected" && exit 1
echo -e "--igvtools    --\n "$(java $JAVAPARAMS -jar $PATH_IGVTOOLS/igvtools.jar version 2>&1)
[ ! -f $PATH_IGVTOOLS/igvtools.jar ] && echo "[ERROR] no igvtools detected" && exit 1
echo -e "--PICARD      --\n "$(java $JAVAPARAMS -jar $PATH_PICARD/MarkDuplicates.jar --version 2>&1)
[ ! -f $PATH_PICARD/MarkDuplicates.jar ] && echo "[ERROR] no picard detected" && exit 1
echo -e "--samstat     --\n "$(samstat -h | head -n 2 | tail -n1)
[ -z "$(which samstat)" ] && echo "[ERROR] no samstat detected" && exit 1
echo -e "--convert     --\n "$(convert -version | head -n 1)
[ -z "$(which convert)" ] && echo "[WARN] imagemagick convert not detected" && exit 1

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "parameters"

# get basename of input file f
INPUTFILENAME=${f##*/}
# get sample prefix
SAMPLE=${INPUTFILENAME/%$READONE.$FASTQ/}

if [ -z "$FASTA" ]; then
    echo "[ERROR] no reference provided (FASTA)"
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
else
    PAIRED="0"
fi

#is paired ?                                                                                                      
if [ "$f" != "${f/%$READONE.$FASTQ/$READTWO.$FASTQ}" ] && [ -e ${f/%$READONE.$FASTQ/$READTWO.$FASTQ} ]; then
    echo "[NOTE] paired-end library detected"
    PAIRED="1"
    READ1=`$CAT $f | wc -l | gawk '{print int($1/4)}' `
    READ2=`$CAT ${f/%$READONE.$FASTQ/$READTWO.$FASTQ} | wc -l | gawk '{print int($1/4)}' `
    let FASTQREADS=$READ1+$READ2
else
    echo "[NOTE] single-end library detected"
    PAIRED="0"
    let FASTQREADS=`$CAT $f | wc -l | gawk '{print int($1/4)}' `
fi

if [[ ! -e ${FASTA%.*}.sa.val ]]; then
    echo "[ERROR] Yara index not detected. Execute yara_indexer first"
    exit 1
fi

THISTMP=$TMP"/"$(whoami)"/"$(echo $OUTDIR/$SAMPLE | md5sum | cut -d' ' -f1)
mkdir -p $THISTMP

FASTASUFFIX=${FASTA##*.}
    
#readgroup
#TODO check readgroups
FULLSAMPLEID=$SAMPLEID"${INPUTFILENAME/%$READONE.$FASTQ/}"
RG="--sam-rg \"ID:$EXPID\" --sam-rg \"SM:$FULLSAMPLEID\" --sam-rg \"LB:$LIBRARY\" --sam-rg \"PL:$PLATFORM\""

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "recall files from tape"
	
if [ -n "$DMGET" ]; then
    dmget -a $FASTA*
    dmget -a ${f/%$READONE.$FASTQ/"*"}
    dmget -a ${OUTDIR}
fi
    
NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "yara"
if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then
        
    if [ "$PAIRED" = "0" ]; then

        RUN_COMMAND="yara_mapper $YARA_MAPPERADDPARAM --threads $CPU_YARA --output-file $THISTMP/$SAMPLE$ALN.sam $FASTA $f"
        echo $RUN_COMMAND && eval $RUN_COMMAND

        
    elif [ "$PAIRED" = "1" ]; then

        RUN_COMMAND="yara_mapper $YARA_MAPPERADDPARAM --threads $CPU_YARA --output-file $THISTMP/$SAMPLE$ALN.sam $FASTA $f ${f/%$READONE.$FASTQ/$READTWO.$FASTQ}"
        echo $RUN_COMMAND && eval $RUN_COMMAND

    fi

    # bam file conversion                                                                        
    samtools view -Sb $THISTMP/$SAMPLE$ALN.sam > $THISTMP/$SAMPLE$ALN.bam 
    samtools sort -@ $CPU_YARA $THISTMP/$SAMPLE$ALN.bam $OUTDIR/$SAMPLE.ash

    # cleanup
    [ -e $OUTDIR/$SAMPLE$ALN.bam ] && rm $OUTDIR/$SAMPLE$ALN.bam
    [ -e $OUTDIR/$SAMPLE$ALN.sam ] && rm $OUTDIR/$SAMPLE$ALN.sam
       
    # mark checkpoint
    [[ -s $OUTDIR/$SAMPLE.ash.bam ]] && NGSANE_CHECKPOINT_CHECK 
    
fi

################################################################################
NGSANE_CHECKPOINT_INIT "add readgroup"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then
    
    java $JAVAPARAMS -jar $PATH_PICARD/AddOrReplaceReadGroups.jar \
        INPUT=$OUTDIR/$SAMPLE.ash.bam \
        OUTPUT=$OUTDIR/$SAMPLE.ashrg.bam \
        RGID=$EXPID RGLB=$LIBRARY RGPL=$PLATFORM \
        RGSM=$FULLSAMPLEID RGPU="XXXXXX"

    # mark checkpoint
    [[ -s $OUTDIR/$SAMPLE.ashrg.bam ]] && NGSANE_CHECKPOINT_CHECK 
    
    [ -e $OUTDIR/$SAMPLE.ash.bam ] && rm $OUTDIR/$SAMPLE.ash.bam
fi 

################################################################################
NGSANE_CHECKPOINT_INIT "mark duplicates"
# create bam files for discarded reads and remove fastq files
if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then
   
    if [ ! -e $OUTDIR/metrices ]; then mkdir -p $OUTDIR/metrices ; fi
    java $JAVAPARAMS -jar $PATH_PICARD/MarkDuplicates.jar \
        INPUT=$OUTDIR/$SAMPLE.ashrg.bam \
        OUTPUT=$OUTDIR/$SAMPLE$ASD.bam \
        METRICS_FILE=$OUTDIR/metrices/$SAMPLE$ASD.bam.dupl AS=true \
        VALIDATION_STRINGENCY=LENIENT \
        TMP_DIR=$THISTMP
        
    samtools index $OUTDIR/$SAMPLE$ASD.bam
          
    # mark checkpoint
    [[ -s $OUTDIR/$SAMPLE$ASD.bam ]] && NGSANE_CHECKPOINT_CHECK
    
    #cleanup
    [ -e $OUTDIR/$SAMPLE.ashrg.bam ] && rm $OUTDIR/$SAMPLE.ashrg.bam
    
fi

################################################################################
NGSANE_CHECKPOINT_INIT "statistics"                                                                                                

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then
    
    STATSOUT=$OUTDIR/$SAMPLE$ASD.bam.stats
    samtools flagstat $OUTDIR/$SAMPLE$ASD.bam > $STATSOUT
    if [ -n "$SEQREG" ]; then
        echo "#custom region" >> $STATSOUT
        echo $(samtools view -@ $CPU_YARA -c -F 4 $OUTDIR/$SAMPLE$ASD.bam $SEQREG )" total reads in region " >> $STATSOUT
        echo $(samtools view -@ $CPU_YARA -c -f 3 $OUTDIR/$SAMPLE$ASD.bam $SEQREG )" properly paired reads in region " >> $STATSOUT
    fi

    # mark checkpoint
    [[ -s $STATSOUT ]] && NGSANE_CHECKPOINT_CHECK
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
    for im in $( ls $OUTDIR/metrices/*.pdf ); do
        convert $im ${im/pdf/jpg}
    done

    # mark checkpoint
    [[ -s $OUTDIR/metrices/$SAMPLE$ASD.bam.alignment_summary_metrics ]] && NGSANE_CHECKPOINT_CHECK 
fi

################################################################################
NGSANE_CHECKPOINT_INIT "coverage track"    

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then
    
    java $JAVAPARAMS -jar $PATH_IGVTOOLS/igvtools.jar count $OUTDIR/$SAMPLE$ASD.bam $OUTDIR/$SAMPLE$ASD.bam.cov.tdf ${FASTA/.$FASTASUFFIX/}.genome
    # mark checkpoint
    [[ -s $OUTDIR/$SAMPLE$ASD.bam.cov.tdf ]] && NGSANE_CHECKPOINT_CHECK
fi

################################################################################
NGSANE_CHECKPOINT_INIT "samstat"    

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then
    
    samstat $OUTDIR/$SAMPLE$ASD.bam

    # mark checkpoint
    [[ -s $OUTDIR/$SAMPLE$ASD.bam.stats ]] && NGSANE_CHECKPOINT_CHECK
    
fi

###############################################################################
NGSANE_CHECKPOINT_INIT "verify"    

BAMREADS=$(head -n1 $OUTDIR/$SAMPLE$ASD.bam.stats | cut -d " " -f 1)
if [ "$BAMREADS" = "" ]; then let BAMREADS="0"; fi
if [ $BAMREADS -eq $FASTQREADS ]; then
    echo "[NOTE] PASS check mapping: $BAMREADS == $FASTQREADS"
else
    echo -e "[ERROR] We are loosing reads from .fastq -> .bam in $f: \nFastq had $FASTQREADS Bam has $BAMREADS"
    exit 1
fi

NGSANE_CHECKPOINT_CHECK

###############################################################################
NGSANE_CHECKPOINT_INIT "cleanup"    

[ -d $THISTMP ] && rm -r $THISTMP

NGSANE_CHECKPOINT_CHECK
################################################################################
[ -e $OUTDIR/$SAMPLE$ASD.bam.dummy ] && rm $OUTDIR/$SAMPLE$ASD.bam.dummy
echo ">>>>> readmapping with Yara - FINISHED"
echo ">>>>> enddate "`date`
