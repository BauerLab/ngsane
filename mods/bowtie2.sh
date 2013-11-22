#!/bin/bash -e

# Script to run bowtie program.
# It takes comma-seprated list of files containing short sequence reads in fasta or fastq format and bowtie index files as input.
# It produces output files: read alignments in .bam format and other files.
# author: Denis Bauer
# date: June 2012
# modified: August 2013 Fabian Buske

# QCVARIABLES,Resource temporarily unavailable
# RESULTFILENAME <DIR>/<TASK>/<SAMPLE>.$ASD.bam

echo ">>>>> readmapping with Bowtie2 "
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

FORCESINGLE=0

#INPUTS                                                                                                           
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the NGSANE repository                       
        -f | --fastq )          shift; f=$1 ;; # fastq file                                                       
        -o | --outdir )         shift; OUTDIR=$1 ;; # output dir                                                     
        -s | --rgsi )           shift; SAMPLEID=$1 ;; # SAMPLEID
        --forceSingle )         FORCESINGLE=1;;
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
for MODULE in $MODULE_BOWTIE2; do module load $MODULE; done  # save way to load modules that itself load other modules
export PATH=$PATH_BOWTIE2:$PATH
module list
echo "PATH=$PATH"
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)
PATH_PICARD=$(dirname $(which MarkDuplicates.jar))

echo "[NOTE] set java parameters"
JAVAPARAMS="-Xmx"$(python -c "print int($MEMORY_BOWTIE2*0.8)")"g -Djava.io.tmpdir="$TMP"  -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1" 
unset _JAVA_OPTIONS
echo "JAVAPARAMS "$JAVAPARAMS

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--JAVA        --\n" $(java -Xmx200m -version 2>&1)
[ -z "$(which java)" ] && echo "[ERROR] no java detected" && exit 1
echo -e "--bowtie2     --\n "$(bowtie2 --version)
[ -z "$(which bowtie2)" ] && echo "[ERROR] no bowtie2 detected" && exit 1
echo -e "--samtools    --\n "$(samtools 2>&1 | head -n 3 | tail -n-2)
[ -z "$(which samtools)" ] && echo "[ERROR] no samtools detected" && exit 1
echo -e "--R           --\n "$(R --version | head -n 3)
[ -z "$(which R)" ] && echo "[ERROR] no R detected" && exit 1
echo -e "--PICARD      --\n "$(java $JAVAPARAMS -jar $PATH_PICARD/MarkDuplicates.jar --version 2>&1)
[ ! -f $PATH_PICARD/MarkDuplicates.jar ] && echo "[ERROR] no picard detected" && exit 1
echo -e "--samstat     --\n "$(samstat -h | head -n 2 | tail -n1)
[ -z "$(which samstat)" ] && echo "[ERROR] no samstat detected" && exit 1
echo -e "--convert     --\n "$(convert -version | head -n 1)
[ -z "$(which convert)" ] && echo "[WARN] imagemagick convert not detected" && exit 1

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="parameters"

# get basename of f
n=${f##*/}

if [ -z "$FASTA" ]; then
    echo "[ERROR] no reference provided (FASTA)"
    exit 1
fi

if [[ ! -e ${FASTA%.*}.1.bt2 ]]; then
    echo "[ERROR] Bowtie2 index not detected. Exeute bowtie2Index.sh first"
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
if [ -z "$RECOVERFROM" ]; then
    [ -e $OUTDIR/${n/%$READONE.$FASTQ/.$ASD.bam} ] && rm $OUTDIR/${n/%$READONE.$FASTQ/.$ASD.bam}
    [ -e $OUTDIR/${n/%$READONE.$FASTQ/.$ASD.bam}.stats ] && rm $OUTDIR/${n/%$READONE.$FASTQ/.$ASD.bam}.stats
    [ -e $OUTDIR/${n/%$READONE.$FASTQ/.$ASD.bam}.dupl ] && rm $OUTDIR/${n/%$READONE.$FASTQ/.$ASD.bam}.dupl
fi

#is ziped ?
ZCAT="zcat"
if [[ ${f##*.} != "gz" ]]; then ZCAT="cat"; fi

#is paired ?                                                                                                      
if [ "$f" != "${f/$READONE/$READTWO}" ] && [ -e ${f/$READONE/$READTWO} ] && [ "$FORCESINGLE" = 0 ]; then
    PAIRED="1"
    READS="-1 $f -2 ${f/$READONE/$READTWO}"
    READ1=`$ZCAT $f | wc -l | gawk '{print int($1/4)}' `
    READ2=`$ZCAT ${f/$READONE/$READTWO} | wc -l | gawk '{print int($1/4)}' `
    let FASTQREADS=$READ1+$READ2
else
    PAIRED="0"
    READS="-U $f"
    let FASTQREADS=`$ZCAT $f | wc -l | gawk '{print int($1/4)}' `
fi

# get encoding
if [ -z "$FASTQ_PHRED" ]; then 
    FASTQ_ENCODING=$($ZCAT $f |  awk 'NR % 4 ==0' | python $NGSANE_BASE/tools/GuessFastqEncoding.py |  tail -n 1)
    if [[ "$FASTQ_ENCODING" == *Phred33* ]]; then
        FASTQ_PHRED="--phred33"    
    elif [[ "$FASTQ_ENCODING" == *Illumina* ]]; then
        FASTQ_PHRED="--phred64"
    elif [[ "$FASTQ_ENCODING" == *Solexa* ]]; then
        FASTQ_PHRED="--solexa-quals"
    else
        echo "[NOTE] cannot detect/don't understand fastq format: $FASTQ_ENCODING - using default"
    fi
    echo "[NOTE] $FASTQ_ENCODING fastq format detected"
fi

#readgroup
FULLSAMPLEID=$SAMPLEID"${n/%$READONE.$FASTQ/}"
RG="--sam-rg \"ID:$EXPID\" --sam-rg \"SM:$FULLSAMPLEID\" --sam-rg \"LB:$LIBRARY\" --sam-rg \"PL:$PLATFORM\""


echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="recall files from tape"
	
if [ -n "$DMGET" ]; then
    dmget -a $FASTA*
	dmget -a ${f/$READONE/"*"}
	dmget -a ${OUTDIR}
fi
    
echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="bowtie2"
if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
    
    RUN_COMMAND="bowtie2 $RG $BOWTIE2ADDPARAM $FASTQ_PHRED -t -x ${FASTA%.*} -p $CPU_BOWTIE2 $READS | samtools view -@ $CPU_BOWTIE2 -bS -t $FASTA.fai - > $OUTDIR/${n/%$READONE.$FASTQ/.$ALN.bam}"
    echo $RUN_COMMAND && eval $RUN_COMMAND
    
    if [ "$PAIRED" = "1" ]; then
        # fix mates
        samtools sort -@ $CPU_BOWTIE2 -n $OUTDIR/${n/%$READONE.$FASTQ/.$ALN.bam} $OUTDIR/${n/%$READONE.$FASTQ/.tmp}
        samtools fixmate $OUTDIR/${n/%$READONE.$FASTQ/.tmp.bam} $OUTDIR/${n/%$READONE.$FASTQ/.$ALN.bam}
        [ -e $OUTDIR/${n/%$READONE.$FASTQ/.tmp.bam} ] && rm $OUTDIR/${n/%$READONE.$FASTQ/.tmp.bam}
    fi

    samtools sort -@ $CPU_BOWTIE2 $OUTDIR/${n/%$READONE.$FASTQ/.$ALN.bam} $OUTDIR/${n/%$READONE.$FASTQ/.ash}

    [ -e $OUTDIR/${n/%$READONE.$FASTQ/.$ALN.bam} ] && rm $OUTDIR/${n/%$READONE.$FASTQ/.$ALN.bam}
       
    # mark checkpoint
    if [ -f $OUTDIR/${n/%$READONE.$FASTQ/.ash.bam} ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
    
fi

################################################################################
CHECKPOINT="mark duplicates"
# create bam files for discarded reads and remove fastq files
if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
   
    if [ ! -e $OUTDIR/metrices ]; then mkdir -p $OUTDIR/metrices ; fi
    THISTMP=$TMP/$n$RANDOM #mk tmp dir because picard writes none-unique files                                        
    mkdir -p $THISTMP
    java $JAVAPARAMS -jar $PATH_PICARD/MarkDuplicates.jar \
        INPUT=$OUTDIR/${n/%$READONE.$FASTQ/.ash.bam} \
        OUTPUT=$OUTDIR/${n/%$READONE.$FASTQ/.$ASD.bam} \
        METRICS_FILE=$OUTDIR/metrices/${n/%$READONE.$FASTQ/.$ASD.bam}.dupl AS=true \
        VALIDATION_STRINGENCY=LENIENT \
        TMP_DIR=$THISTMP
    [ -d $THISTMP ] && rm -r $THISTMP
    samtools index $OUTDIR/${n/%$READONE.$FASTQ/.$ASD.bam}
          
    # mark checkpoint
    if [ -f $OUTDIR/${n/%$READONE.$FASTQ/.$ASD.bam} ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
    
    #cleanup
    [ -e $OUTDIR/${n/%$READONE.$FASTQ/.ash.bam} ] && rm $OUTDIR/${n/%$READONE.$FASTQ/.ash.bam}
    
fi

################################################################################
CHECKPOINT="statistics"                                                                                                

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
    
    STATSOUT=$OUTDIR/${n/%$READONE.$FASTQ/.$ASD.bam}.stats
    samtools flagstat $OUTDIR/${n/%$READONE.$FASTQ/.$ASD.bam} > $STATSOUT
    if [ -n "$SEQREG" ]; then
        echo "#custom region" >> $STATSOUT
        echo $(samtools view -@ $CPU_BOWTIE2 -c -F 4 $OUTDIR/${n/%$READONE.$FASTQ/.$ASD.bam} $SEQREG )" total reads in region " >> $STATSOUT
        echo $(samtools view -@ $CPU_BOWTIE2 -c -f 3 $OUTDIR/${n/%$READONE.$FASTQ/.$ASD.bam} $SEQREG )" properly paired reads in region " >> $STATSOUT
    fi

    # mark checkpoint
    if [ -e $STATSOUT ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
fi

################################################################################
CHECKPOINT="calculate inner distance"                                                                                                

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

    THISTMP=$TMP/$n$RANDOM #mk tmp dir because picard writes none-unique files
    mkdir $THISTMP
    java $JAVAPARAMS -jar $PATH_PICARD/CollectMultipleMetrics.jar \
        INPUT=$OUTDIR/${n/%$READONE.$FASTQ/.$ASD.bam} \
        REFERENCE_SEQUENCE=$FASTA \
        OUTPUT=$OUTDIR/metrices/${n/%$READONE.$FASTQ/.$ASD.bam} \
        VALIDATION_STRINGENCY=LENIENT \
        PROGRAM=CollectAlignmentSummaryMetrics \
        PROGRAM=CollectInsertSizeMetrics \
        PROGRAM=QualityScoreDistribution \
        TMP_DIR=$THISTMP
    for im in $( ls $OUTDIR/metrices/*.pdf ); do
        convert $im ${im/pdf/jpg}
    done
    rm -r $THISTMP

    # mark checkpoint
    if [ -f $OUTDIR/metrices/${n/%$READONE.$FASTQ/.$ASD.bam}.alignment_summary_metrics ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
fi

################################################################################
CHECKPOINT="samstat"    

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
    
    samstat $OUTDIR/${n/%$READONE.$FASTQ/.$ASD.bam}

    # mark checkpoint
    if [ -f $OUTDIR/${n/%$READONE.$FASTQ/.$ASD.bam}.stats ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
    
fi

###############################################################################
CHECKPOINT="verify"    

BAMREADS=`head -n1 $OUTDIR/${n/%$READONE.$FASTQ/.$ASD.bam}.stats | cut -d " " -f 1`
if [ "$BAMREADS" = "" ]; then let BAMREADS="0"; fi
if [ $BAMREADS -eq $FASTQREADS ]; then
    echo "[NOTE] PASS check mapping: $BAMREADS == $FASTQREADS"
else
    echo -e "[ERROR] We are loosing reads from .fastq -> .bam in $f: \nFastq had $FASTQREADS Bam has $BAMREADS"
    exit 1
fi

echo -e "\n********* $CHECKPOINT\n"
################################################################################
[ -e $OUTDIR/${n/%$READONE.$FASTQ/.$ASD.bam}.dummy ] && rm $OUTDIR/${n/%$READONE.$FASTQ/.$ASD.bam}.dummy
echo ">>>>> readmapping with Bowtie2 - FINISHED"
echo ">>>>> enddate "`date`
