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
if [ -n "$ITERATIVE_MAPPING" ]; then 
    PATH_TRIMMOMATIC=$(dirname $(which trimmomatic.jar))
    echo -e "--trimmomatic --\n " $(which $PATH_TRIMMOMATIC/trimmomatic.jar)
    [ ! -f $PATH_TRIMMOMATIC/trimmomatic.jar ] && echo "[ERROR] no trimmomatic detected" && exit 1
fi

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
    [ -e $OUTDIR/$SAMPLE.$ASD.bam ] && rm $OUTDIR/$SAMPLE.$ASD.bam
    [ -e $OUTDIR/$SAMPLE.$ASD.bam.stats ] && rm $OUTDIR/$SAMPLE.$ASD.bam.stats
    [ -e $OUTDIR/$SAMPLE.$ASD.bam.dupl ] && rm $OUTDIR/$SAMPLE.$ASD.bam.dupl
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

        echo "[NOTE] PAIRED READS"
        # clever use of named pipes to avoid IO
        [ -e $OUTDIR/$SAMPLE$READONE.sai ] && rm $OUTDIR/$SAMPLE$READONE.sai
        [ -e $OUTDIR/$SAMPLE$READTWO.sai ] && rm $OUTDIR/$SAMPLE$READTWO.sai
        mkfifo $OUTDIR/$SAMPLE$READONE.sai $OUTDIR/${n/%$READONE.$FASTQ/$READTWO.sai}
    
        bwa aln $QUAL $BWAALNADDPARAM $FASTQ_PHRED -t $CPU_BWA $FASTA $f > $OUTDIR/$SAMPLE$READONE.sai &
        bwa aln $QUAL $BWAALNADDPARAM $FASTQ_PHRED -t $CPU_BWA $FASTA ${f/%$READONE.$FASTQ/$READTWO.$FASTQ} > $OUTDIR/$SAMPLE$READTWO.sai &
        
        bwa sampe $FASTA $OUTDIR/$SAMPLE$READONE.sai $OUTDIR/$SAMPLE$READTWO.sai\
    	       $BWASAMPLEADDPARAM -r "@RG\tID:$EXPID\tSM:$FULLSAMPLEID\tPL:$PLATFORM\tLB:$LIBRARY" \
	       $f ${f/%$READONE.$FASTQ/$READTWO.$FASTQ} | samtools view -bS -t $FASTA.fai - > $OUTDIR/$SAMPLE.$ALN.bam

        [ -e $OUTDIR/$SAMPLE$READONE.sai ] && rm $OUTDIR/$SAMPLE$READONE.sai
        [ -e $OUTDIR/$SAMPLE$READTWO.sai ] && rm $OUTDIR/$SAMPLE$READTWO.sai

    else
        echo "[NOTE] SINGLE READS"
        # clever use of named pipes to avoid IO
        [ -e $OUTDIR/$SAMPLE$READONE.sai ] && rm $OUTDIR/$SAMPLE$READONE.sai
        mkfifo $OUTDIR/$SAMPLE$READONE.sai
        
        bwa aln $QUAL $BWAALNADDPARAM -t $CPU_BWA $FASTA $f > $OUTDIR/$SAMPLE$READONE.sai &
    
        bwa samse $FASTA $OUTDIR/$SAMPLE$READONE.sai \
            $BWASAMPLEADDPARAM -r "@RG\tID:$EXPID\tSM:$FULLSAMPLEID\tPL:$PLATFORM\tLB:$LIBRARY" \
        	$f | samtools view -bS -t $FASTA.fai - > $OUTDIR/$SAMPLE.$ALN.bam
    
        [ -e $OUTDIR/$SAMPLE$READONE.sai ] && rm $OUTDIR/$SAMPLE$READONE.sai
    fi
    
    # mark checkpoint
    if [ -f $OUTDIR/$SAMPLE.$ALN.bam ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
fi

################################################################################
CHECKPOINT="Iterative Mapping"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else
    # iterative mapping    
    if [ -z "$ITERATIVE_MAPPING" ]; then
        echo "[NOTE] Skip iterative mapping"
        echo -e "\n********* $CHECKPOINT\n"
        
    else
        ITERATIVEINPUT=$OUTDIR/$SAMPLE.$ALN.bam
        samtools view -H $ITERATIVEINPUT > $THISTMP/$SAMPLE.$ALN.fullsize.header
        
        if [ "$PAIRED" = 1 ]; then
            echo "[NOTE] iterative mapping of paired-end libary"
            
            for COUNTER in $(seq $ITERATIVE_MAPPING); do
                # skip ahead if possible
                if [ -f $OUTDIR/$SAMPLE.$ALN.$COUNTER.bam ]; then
                    echo "[NOTE] found results for iterative mapping $COUNTER"
                else    
                    echo "[NOTE] Iterative Mapping with read length $COUNTER"
                    # get exactly one end mapped
                    samtools view -bh -F 4 -f 8 $ITERATIVEINPUT > $THISTMP/$SAMPLE.$UNM.1.bam
                    samtools view -bh -F 8 -f 4 $ITERATIVEINPUT > $THISTMP/$SAMPLE.$UNM.2.bam
                    
                    # both two ends unmapped
                    samtools view -bh -f 12 $ITERATIVEINPUT > $THISTMP/$SAMPLE.$UNM.3.bam
                    
                    # merge and sort
                    samtools merge $THISTMP/$SAMPLE.$UNM.all.bam $THISTMP/$SAMPLE.$UNM.1.bam $THISTMP/$SAMPLE.$UNM.2.bam $THISTMP/$SAMPLE.$UNM.3.bam               
                    samtools sort $THISTMP/$SAMPLE.$UNM.all.bam $THISTMP/$SAMPLE.$UNM.tmp
                    
                    # check for premature exit from loop
                    if [ $(samtools view -c $THISTMP/$SAMPLE.$UNM.tmp.bam) == 0 ]; then 
                        break
                    fi               
    
                    # properly paired reads
                    samtools view -bh -f 1 -F 12 $ITERATIVEINPUT > $OUTDIR/$SAMPLE.$ALN.$COUNTER.iterative
    
                    rm $THISTMP/$SAMPLE.$UNM.1.bam $THISTMP/$SAMPLE.$UNM.2.bam $THISTMP/$SAMPLE.$UNM.3.bam $THISTMP/$SAMPLE.$UNM.all.bam
    
                    # convert
                    [ -e $THISTMP/$SAMPLE$READONE.fastq.gz ] && rm $THISTMP/$SAMPLE$READONE.fastq.gz
                    [ -e $THISTMP/$SAMPLE$READTWO.fastq.gz ] && rm $THISTMP/$SAMPLE$READTWO.fastq.gz
                    
                    RUN_COMMAND="java $JAVAPARAMS -jar $PATH_PICARD/CleanSam.jar \
                        I=$THISTMP/$SAMPLE.$UNM.tmp.bam \
                        O=$THISTMP/$SAMPLE.$UNM.bam \
                        VALIDATION_STRINGENCY=SILENT \
                        TMP_DIR=$THISTMP"
                    echo $RUN_COMMAND && eval $RUN_COMMAND
    
                    samtools index $THISTMP/$SAMPLE.$UNM.bam
    
                    java $JAVAPARAMS -jar $PATH_PICARD/SamToFastq.jar INPUT=$THISTMP/$SAMPLE.$UNM.bam FASTQ=$THISTMP/$SAMPLE$READONE.fastq SECOND_END_FASTQ=$THISTMP/$SAMPLE$READTWO.fastq
                    $GZIP $THISTMP/$SAMPLE$READONE.fastq $THISTMP/$SAMPLE$READTWO.fastq
                    
                    # crop reads using trimmomatic
                    RUN_COMMAND="java $JAVAPARAMS -jar $PATH_TRIMMOMATIC/trimmomatic.jar PE $FASTQ_PHRED_TRIM -threads $CPU_BWA -trimlog $THISTMP/$SAMPLE.trim.log $THISTMP/$SAMPLE$READONE.fastq.gz $THISTMP/$SAMPLE$READTWO.fastq.gz $THISTMP/$SAMPLE$READONE.cropped.fastq.gz $THISTMP/$SAMPLE$READONE.unpaired.fastq.gz $THISTMP/$SAMPLE$READTWO.cropped.fastq.gz $THISTMP/$SAMPLE$READTWO.unpaired.fastq.gz CROP:$COUNTER"
                    echo $RUN_COMMAND && eval $RUN_COMMAND
                    
                    [ -e $THISTMP/$SAMPLE$READONE.sai ] && rm $THISTMP/$SAMPLE$READONE.sai
                    [ -e $THISTMP/$SAMPLE$READTWO.sai ] && rm $THISTMP/$SAMPLE$READTWO.sai
                    mkfifo $THISTMP/$SAMPLE$READONE.sai $THISTMP/${n/%$READONE.$FASTQ/$READTWO.sai}
                
                    bwa aln $QUAL $BWAALNADDPARAM $FASTQ_PHRED -t $CPU_BWA $FASTA $THISTMP/$SAMPLE$READONE.cropped.fastq.gz > $THISTMP/$SAMPLE$READONE.sai &
                    bwa aln $QUAL $BWAALNADDPARAM $FASTQ_PHRED -t $CPU_BWA $FASTA $THISTMP/$SAMPLE$READTWO.cropped.fastq.gz > $THISTMP/$SAMPLE$READTWO.sai &
                    
                    bwa sampe $FASTA $THISTMP/$SAMPLE$READONE.sai $THISTMP/$SAMPLE$READTWO.sai\
            	       $BWASAMPLEADDPARAM -r "@RG\tID:$EXPID\tSM:$FULLSAMPLEID\tPL:$PLATFORM\tLB:$LIBRARY" \
            	       $THISTMP/$SAMPLE$READONE.cropped.fastq.gz $THISTMP/$SAMPLE$READTWO.cropped.fastq.gz | samtools -bS -t $FASTA.fai - > $THISTMP/$SAMPLE.$ALN.$COUNTER.bam
            
                    [ -e $THISTMP/$SAMPLE$READONE.sai ] && rm $THISTMP/$SAMPLE$READONE.sai
                    [ -e $THISTMP/$SAMPLE$READTWO.sai ] && rm $THISTMP/$SAMPLE$READTWO.sai
                        
                    rm $THISTMP/$SAMPLE$READONE.cropped.* $THISTMP/$SAMPLE$READTWO.cropped.* $THISTMP/$SAMPLE.$UNM.bam*
                    mv $THISTMP/$SAMPLE.$ALN.$COUNTER.bam $OUTDIR/$SAMPLE.$ALN.$COUNTER.bam
                fi
                
                ITERATIVEINPUT=$OUTDIR/$SAMPLE.$ALN.$COUNTER.bam
            done

        elif [ "$PAIRED" = 0 ]; then
            echo "[NOTE] iterative mapping of single-end libary"

            for COUNTER in $(seq $ITERATIVE_MAPPING); do
                # skip ahead if possible
                if [ -f $OUTDIR/$SAMPLE.$ALN.$COUNTER.bam ]; then
                    echo "[NOTE] found results for iterative mapping $COUNTER"
                else    
                    echo "[NOTE] Iterative Mapping with read length $COUNTER"
                    
                    # exactly unmapped
                    samtools view -bh -f 4 $ITERATIVEINPUT > $THISTMP/$SAMPLE.$UNM.tmp.bam
    
                    # check for premature exit from loop
                    if [ $(samtools view -c $THISTMP/$SAMPLE.$UNM.tmp.bam) == 0 ]; then 
                        break
                    fi
    
                    # properly mapped reads
                    samtools view -bh -F 4 $ITERATIVEINPUT > $OUTDIR/$SAMPLE.$ALN.$COUNTER.iterative
                    
                                    
                    RUN_COMMAND="java $JAVAPARAMS -jar $PATH_PICARD/CleanSam.jar \
                        I=$THISTMP/$SAMPLE.$UNM.tmp.bam \
                        O=$THISTMP/$SAMPLE.$UNM.bam \
                        VALIDATION_STRINGENCY=SILENT \
                        TMP_DIR=$THISTMP"
                    echo $RUN_COMMAND && eval $RUN_COMMAND
    
                    samtools index $THISTMP/$SAMPLE.$UNM.bam
    
                    # convert
                    [ -e $THISTMP/$SAMPLE$READONE.fastq.gz ] && rm $THISTMP/$SAMPLE$READONE.fastq.gz
                    java $JAVAPARAMS -jar $PATH_PICARD/SamToFastq.jar INPUT=$THISTMP/$SAMPLE.$UNM.bam FASTQ=$THISTMP/$SAMPLE$READONE.fastq
                    $GZIP $THISTMP/$SAMPLE$READONE.fastq
                    
                    # crop reads using trimmomatic
                    RUN_COMMAND="java $JAVAPARAMS -jar $PATH_TRIMMOMATIC/trimmomatic.jar SE $FASTQ_PHRED_TRIM -threads $CPU_BWA -trimlog $THISTMP/$SAMPLE.trim.log $THISTMP/$SAMPLE$READONE.fastq.gz $THISTMP/$SAMPLE$READONE.cropped.fastq.gz CROP:$COUNTER"
                    echo $RUN_COMMAND && eval $RUN_COMMAND
                    
                    [ -e $THISTMP/$SAMPLE$READONE.sai ] && rm $THISTMP/$SAMPLE$READONE.sai
                    mkfifo $THISTMP/$SAMPLE$READONE.sai
                
                    bwa aln $QUAL $BWAALNADDPARAM $FASTQ_PHRED -t $CPU_BWA $FASTA $THISTMP/$SAMPLE$READONE.cropped.fastq.gz > $THISTMP/$SAMPLE$READONE.sai &
                    
                    bwa samse $FASTA $THISTMP/$SAMPLE$READONE.sai \
            	       $BWASAMPLEADDPARAM -r "@RG\tID:$EXPID\tSM:$FULLSAMPLEID\tPL:$PLATFORM\tLB:$LIBRARY" \
            	       $THISTMP/$SAMPLE$READONE.cropped.fastq.gz | samtools view -bS -t $FASTA.fai - > $THISTMP/$SAMPLE.$ALN.$COUNTER.bam
            
                    [ -e $THISTMP/$SAMPLE$READONE.sai ] && rm $THISTMP/$SAMPLE$READONE.sai
                        
                    rm $THISTMP/$SAMPLE$READONE.cropped.* $THISTMP/$SAMPLE.$UNM.bam*
                    mv $THISTMP/$SAMPLE.$ALN.$COUNTER.bam $OUTDIR/$SAMPLE.$ALN.$COUNTER.bam
                fi
                ITERATIVEINPUT=$OUTDIR/$SAMPLE.$ALN.$COUNTER.bam

            done
        fi
        
        # merge iterative mapped reads
        samtools merge -h $THISTMP/$SAMPLE.$ALN.fullsize.header $THISTMP/$SAMPLE.$ALN.tmp $ITERATIVEINPUT $OUTDIR/$SAMPLE.$ALN.[0-9]*.iterative
        mv $THISTMP/$SAMPLE.$ALN.tmp $OUTDIR/$SAMPLE.$ALN.bam
        
        rm $THISTMP/$SAMPLE.$ALN.fullsize.header $ITERATIVEINPUT $OUTDIR/$SAMPLE.$ALN.[0-9]*.bam $OUTDIR/$SAMPLE.$ALN.[0-9]*.iterative

        # mark checkpoint
        if [ -f $OUTDIR/$SAMPLE.$ALN.bam ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
    fi
fi

################################################################################
CHECKPOINT="bam conversion and sorting"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

    samtools sort -@ $CPU_BWA $OUTDIR/$SAMPLE.$ALN.bam $OUTDIR/$SAMPLE.ash
    [ -e $OUTDIR/$SAMPLE.$ALN.bam ] && rm $OUTDIR/$SAMPLE.$ALN.bam

    # mark checkpoint
    if [ -f $OUTDIR/$SAMPLE.ash.bam ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
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
    java $JAVAPARAMS -jar $PATH_PICARD/MarkDuplicates.jar \
        INPUT=$OUTDIR/$SAMPLE.ash.bam \
        OUTPUT=$OUTDIR/$SAMPLE.$ASD.bam \
        METRICS_FILE=$OUTDIR/metrices/$SAMPLE.$ASD.bam.dupl \
        AS=true \
        CREATE_MD5_FILE=true \
        COMPRESSION_LEVEL=9 \
        VALIDATION_STRINGENCY=SILENT \
        TMP_DIR=$THISTMP
    samtools index $OUTDIR/$SAMPLE.$ASD.bam

    # mark checkpoint
    if [ -f $OUTDIR/$SAMPLE.$ASD.bam ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
fi


################################################################################
CHECKPOINT="statistics"                                                                                                

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
        
    STATSOUTDIR=$OUTDIR/$SAMPLE.$ASD.bam.stats
    samtools flagstat $OUTDIR/$SAMPLE.$ASD.bam > $STATSOUTDIR
    if [ -n "$SEQREG" ]; then
        echo "#custom region" >> $STATSOUTDIR
        echo $(samtools view -c -F 4 $OUTDIR/$SAMPLE.$ASD.bam $SEQREG )" total reads in region " >> $STATSOUTDIR
        echo $(samtools view -c -f 3 $OUTDIR/$SAMPLE.$ASD.bam $SEQREG )" properly paired reads in region " >> $STATSOUTDIR
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
    java $JAVAPARAMS -jar $PATH_PICARD/CollectMultipleMetrics.jar \
        INPUT=$OUTDIR/$SAMPLE.$ASD.bam \
        REFERENCE_SEQUENCE=$FASTA \
        OUTPUT=$OUTDIR/metrices/$SAMPLE.$ASD.bam \
        VALIDATION_STRINGENCY=SILENT \
        PROGRAM=CollectAlignmentSummaryMetrics \
        PROGRAM=CollectInsertSizeMetrics \
        PROGRAM=QualityScoreDistribution \
        TMP_DIR=$THISTMP
    for im in $( ls $OUTDIR/metrices/*.pdf ); do
        convert $im ${im/pdf/jpg}
    done

    # mark checkpoint
    if [ -f $OUTDIR/metrices/$SAMPLE.$ASD.bam.alignment_summary_metrics ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
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
    
BAMREADS=$(head -n1 $OUTDIR/$SAMPLE.$ASD.bam.stats | cut -d " " -f 1)
if [ "$BAMREADS" = "" ]; then let BAMREADS="0"; fi			
if [ $BAMREADS -eq $FASTQREADS ]; then
    echo "[NOTE] PASS check mapping: $BAMREADS == $FASTQREADS"
    [ -e $OUTDIR/$SAMPLE.ash.bam ] && rm $OUTDIR/$SAMPLE.ash.bam
else
    echo -e "[ERROR] We are loosing reads from .fastq -> .bam in $f: \nFastq had $FASTQREADS Bam has $BAMREADS"
    exit 1 
fi

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="cleanup" 

[ -d $THISTMP ] && rm -r $THISTMP

echo -e "\n********* $CHECKPOINT\n"
################################################################################
[ -e $OUTDIR/$SAMPLE.$ASD.bam.dummy ] && rm $OUTDIR/$SAMPLE.$ASD.bam.dummy
echo ">>>>> readmapping with BWA - FINISHED"
echo ">>>>> enddate "`date`

