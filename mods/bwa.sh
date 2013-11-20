#!/bin/bash -e

# BWA calling script
# author: Denis C. Bauer
# date: Nov.2010
# modified: August 2013 - Fabian Buske

# messages to look out for -- relevant for the QC.sh script:
# QCVARIABLES,We are loosing reads,MAPQ should be 0 for unmapped read,no such file,file not found,bwa.sh: line,Resource temporarily unavailable
# RESULTFILENAME <DIR>/<TASK>/<SAMPLE>.$ASD.bam

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
  --forceSingle             run single end eventhough second read is present
  --noMapping
"
exit
}

if [ ! $# -gt 3 ]; then usage ; fi

#DEFAULTS
FORCESINGLE=0
NOMAPPING=0

#INPUTS
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the NGSANE repository
        -f | --fastq )          shift; f=$1 ;; # fastq file
        -o | --outdir )         shift; OUTDIR=$1 ;; # output dir
        -s | --rgsi )           shift; SAMPLEID=$1 ;; # read group prefix
        --forceSingle )         FORCESINGLE=1;;
        --noMapping )           NOMAPPING=1;;
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

for MODULE in $MODULE_BWA; do module load $MODULE; done  # save way to load modules that itself load other modules
export PATH=$PATH_BWA:$PATH
module list
echo "PATH=$PATH"
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)
PATH_PICARD=$(dirname $(which MarkDuplicates.jar))

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


echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="parameters"

if [[ -z "$EXPID" || -z "$LIBRARY" || -z "$PLATFORM" ]]; then
    echo "[ERROR] library info not set (EXPID, LIBRARY, and PLATFORM): free text needed"
    exit 1;
else
    echo "[NOTE] EXPID $EXPID; LIBRARY $LIBRARY; PLATFORM $PLATFORM"
fi

# get basename of f
n=${f##*/}

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
if [ -z "$RECOVERFROM" ]; then
    [ -e $OUTDIR/${n/%$READONE.$FASTQ/.$ASD.bam} ] && rm $OUTDIR/${n/%$READONE.$FASTQ/.$ASD.bam}
    [ -e $OUTDIR/${n/%$READONE.$FASTQ/.$ASD.bam}.stats ] && rm $OUTDIR/${n/%$READONE.$FASTQ/.$ASD.bam}.stats
    [ -e $OUTDIR/${n/%$READONE.$FASTQ/.$ASD.bam}.dupl ] && rm $OUTDIR/${n/%$READONE.$FASTQ/.$ASD.bam}.dupl
fi

#is ziped ?
ZCAT="zcat"
if [[ $f != *.gz ]]; then ZCAT="cat"; fi

#is paired ?
if [ "$f" != "${f/$READONE/$READTWO}" ] && [ -e ${f/$READONE/$READTWO} ] && [ "$FORCESINGLE" = 0 ]; then
    PAIRED="1"
    READ1=$($ZCAT $f | wc -l | gawk '{print int($1/4)}')
    READ2=$($ZCAT ${f/$READONE/$READTWO} | wc -l | gawk '{print int($1/4)}')
    let FASTQREADS=$READ1+$READ2
else
    PAIRED="0"
    READS="$f"
    let FASTQREADS=`$ZCAT $f | wc -l | gawk '{print int($1/4)}' `
fi

# get encoding
FASTQ_ENCODING=$($ZCAT $f |  awk 'NR % 4 ==0' | python $NGSANE_BASE/tools/GuessFastqEncoding.py |  tail -n 1)
if [[ "$FASTQ_ENCODING" == *Phred33* ]]; then
    FASTQ_PHRED=""    
elif [[ "$FASTQ_ENCODING" == *Illumina* ]]; then
    FASTQ_PHRED="-I"
elif [[ "$FASTQ_ENCODING" == *Solexa* ]]; then
    FASTQ_PHRED="-I"
else
    echo "[NOTE] cannot detect/don't understand fastq format: $FASTQ_ENCODING - using default"
fi
echo "[NOTE] $FASTQ_ENCODING fastq format detected"

FULLSAMPLEID=$SAMPLEID"${n/%$READONE.$FASTQ/}"

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
CHECKPOINT="BWA"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else
    
    if [ "$PAIRED" = 1 ]; then

        if [ "$NOMAPPING" = 0 ]; then
            echo "[NOTE] PAIRED READS"
            # clever use of named pipes to avoid IO
            [ -e $OUTDIR/${n/$FASTQ/sai} ] && rm $OUTDIR/${n/$FASTQ/sai}
            [ -e $OUTDIR/${n/$READONE.$FASTQ/$READTWO.sai} ] && rm $OUTDIR/${n/$READONE.$FASTQ/$READTWO.sai} 
            mkfifo $OUTDIR/${n/$FASTQ/sai} $OUTDIR/${n/$READONE.$FASTQ/$READTWO.sai}
        
            bwa aln $QUAL $BWAALNADDPARAM $FASTQ_PHRED -t $CPU_BWA $FASTA $f > $OUTDIR/${n/$FASTQ/sai} &
            bwa aln $QUAL $BWAALNADDPARAM $FASTQ_PHRED -t $CPU_BWA $FASTA ${f/$READONE/$READTWO} > $OUTDIR/${n/$READONE.$FASTQ/$READTWO.sai} &
            bwa sampe $FASTA $OUTDIR/${n/$FASTQ/sai} $OUTDIR/${n/$READONE.$FASTQ/$READTWO.sai} \
       	       $BWASAMPLEADDPARAM -r "@RG\tID:$EXPID\tSM:$FULLSAMPLEID\tPL:$PLATFORM\tLB:$LIBRARY" \
    	       $f ${f/$READONE/$READTWO} | samtools view -@ $CPU_BWA -bS -t $FASTA.fai - > $OUTDIR/${n/%$READONE.$FASTQ/.$ALN.bam}
    
            [ -e $OUTDIR/${n/$FASTQ/sai} ] && rm $OUTDIR/${n/$FASTQ/sai}
            [ -e $OUTDIR/${n/$READONE.$FASTQ/$READTWO.sai} ] && rm $OUTDIR/${n/$READONE.$FASTQ/$READTWO.sai}
        fi

    else
        echo "[NOTE] SINGLE READS"
        # clever use of named pipes to avoid IO
        [ -e $OUTDIR/${n/$FASTQ/sai} ] && rm $OUTDIR/${n/$FASTQ/sai}
        mkfifo $OUTDIR/${n/$FASTQ/sai}
        
        bwa aln $QUAL $BWAALNADDPARAM -t $CPU_BWA $FASTA $f > $OUTDIR/${n/$FASTQ/sai} &
    
        bwa samse $FASTA $OUTDIR/${n/$FASTQ/sai} $BWASAMPLEADDPARAM \
    	-r "@RG\tID:$EXPID\tSM:$FULLSAMPLEID\tPL:$PLATFORM\tLB:$LIBRARY" \
    	$f | samtools view -@ $CPU_BWA -bS -t $FASTA.fai - > $OUTDIR/${n/%$READONE.$FASTQ/.$ALN.bam}
    
        [ -e $OUTDIR/${n/$FASTQ/sai} ] && rm $OUTDIR/${n/$FASTQ/sai}
    fi
    
    # mark checkpoint
    if [ -f $OUTDIR/${n/%$READONE.$FASTQ/.$ALN.bam} ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
fi


################################################################################
CHECKPOINT="bam conversion and sorting"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

    samtools sort -@ $CPU_BWA $OUTDIR/${n/%$READONE.$FASTQ/.$ALN.bam} $OUTDIR/${n/%$READONE.$FASTQ/.ash}
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
    
    
    
    #TODO look at samtools for rmdup
    #val string had to be set to LENIENT (SIlENT) to avoid crash due to a definition dis-
    #agreement between bwa and picard
    #http://seqanswers.com/forums/showthread.php?t=4246
    if [ ! -e $OUTDIR/metrices ]; then mkdir -p $OUTDIR/metrices ; fi
    THISTMP=$TMP/$n$RANDOM #mk tmp dir because picard writes none-unique files
    echo $THISTMP
    mkdir -p $THISTMP
    java $JAVAPARAMS -jar $PATH_PICARD/MarkDuplicates.jar \
        INPUT=$OUTDIR/${n/%$READONE.$FASTQ/.ash.bam} \
        OUTPUT=$OUTDIR/${n/%$READONE.$FASTQ/.$ASD.bam} \
        METRICS_FILE=$OUTDIR/metrices/${n/%$READONE.$FASTQ/.$ASD.bam}.dupl AS=true \
        VALIDATION_STRINGENCY=SILENT \
        TMP_DIR=$THISTMP
    [ -d $THISTMP ] && rm -r $THISTMP
    samtools index $OUTDIR/${n/%$READONE.$FASTQ/.$ASD.bam}

    # mark checkpoint
    if [ -f $OUTDIR/${n/%$READONE.$FASTQ/.$ASD.bam} ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
fi


################################################################################
CHECKPOINT="statistics"                                                                                                

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
        
    STATSOUTDIR=$OUTDIR/${n/%$READONE.$FASTQ/.$ASD.bam}.stats
    samtools flagstat $OUTDIR/${n/%$READONE.$FASTQ/.$ASD.bam} > $STATSOUTDIR
    if [ -n "$SEQREG" ]; then
        echo "#custom region" >> $STATSOUTDIR
        echo $(samtools view -c -F 4 $OUTDIR/${n/%$READONE.$FASTQ/.$ASD.bam} $SEQREG )" total reads in region " >> $STATSOUTDIR
        echo $(samtools view -c -f 3 $OUTDIR/${n/%$READONE.$FASTQ/.$ASD.bam} $SEQREG )" properly paired reads in region " >> $STATSOUTDIR
    fi

    # mark checkpoint
    if [ -f $STATSOUT ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
fi


################################################################################
CHECKPOINT="calculate inner distance"                                                                                                

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
      
    export PATH=$PATH:/usr/bin/
    THISTMP=$TMP/$n$RANDOM #mk tmp dir because picard writes none-unique files
    mkdir -p $THISTMP
    java $JAVAPARAMS -jar $PATH_PICARD/CollectMultipleMetrics.jar \
        INPUT=$OUTDIR/${n/%$READONE.$FASTQ/.$ASD.bam} \
        REFERENCE_SEQUENCE=$FASTA \
        OUTPUT=$OUTDIR/metrices/${n/%$READONE.$FASTQ/.$ASD.bam} \
        VALIDATION_STRINGENCY=SILENT \
        PROGRAM=CollectAlignmentSummaryMetrics \
        PROGRAM=CollectInsertSizeMetrics \
        PROGRAM=QualityScoreDistribution \
        TMP_DIR=$THISTMP
    for im in $( ls $OUTDIR/metrices/*.pdf ); do
        convert $im ${im/pdf/jpg}
    done
    [ -d $THISTMP ] && rm -r $THISTMP

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

################################################################################
CHECKPOINT="verify"    
    
BAMREADS=$(head -n1 $OUTDIR/${n/%$READONE.$FASTQ/.$ASD.bam}.stats | cut -d " " -f 1)
if [ "$BAMREADS" = "" ]; then let BAMREADS="0"; fi			
if [ $BAMREADS -eq $FASTQREADS ]; then
    echo "[NOTE] PASS check mapping: $BAMREADS == $FASTQREADS"
    [ -e $OUTDIR/${n/%$READONE.$FASTQ/.ash.bam} ] && rm $OUTDIR/${n/%$READONE.$FASTQ/.ash.bam}
else
    echo -e "[ERROR] We are loosing reads from .fastq -> .bam in $f: \nFastq had $FASTQREADS Bam has $BAMREADS"
    exit 1 
fi

echo -e "\n********* $CHECKPOINT\n"
################################################################################
[ -e $OUTDIR/${n/%$READONE.$FASTQ/.$ASD.bam}.dummy ] && rm $OUTDIR/${n/%$READONE.$FASTQ/.$ASD.bam}.dummy
echo ">>>>> readmapping with BWA - FINISHED"
echo ">>>>> enddate "`date`

