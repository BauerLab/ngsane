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

if [ -z "$FASTA" ]; then
    echo "[ERROR] no reference provided (FASTA)"
    exit 1
fi

# delete old bam files unless attempting to recover
if [ -z "$RECOVERFROM" ]; then
    [ -e $OUTDIR/${n/%$READONE.$FASTQ/.$ASD.bam} ] && rm $OUTDIR/${n/%$READONE.$FASTQ/.$ASD.bam}
    [ -e $OUTDIR/${n/%$READONE.$FASTQ/.$ASD.bam}.stats ] && rm $OUTDIR/${n/%$READONE.$FASTQ/.$ASD.bam}.stats
    [ -e $OUTDIR/${n/%$READONE.$FASTQ/.$ASD.bam}.dupl ] && rm $OUTDIR/${n/%$READONE.$FASTQ/.$ASD.bam}.dupl
fi

#is paired ?                                                                                                      
if [ "$f" != "${f/$READONE/$READTWO}" ] && [ -e ${f/$READONE/$READTWO} ] && [ "$FORCESINGLE" = 0 ]; then
    PAIRED="1"
    echo 
else
    PAIRED="0"
fi

#is ziped ?                                                                                                       
ZCAT="zcat"
if [[ ${f##*.} != "gz" ]]; then ZCAT="cat"; fi


# get encoding
if [ -z "$FASTQ_PHRED" ]; then 
    FASTQ_ENCODING=$($ZCAT $f |  awk 'NR % 4 ==0' | python $NGSANE_BASE/tools/GuessFastqEncoding.py |  tail -n 1)
    if [[ "$FASTQ_ENCODING" == *Phred33* ]]; then
        FASTQ_PHRED="--phred33-quals"    
    elif [[ "$FASTQ_ENCODING" == *Illumina* ]]; then
        FASTQ_PHRED="--phred64-quals"
    elif [[ "$FASTQ_ENCODING" == *Solexa* ]]; then
        FASTQ_PHRED="--solexa1.3-quals"
    else
        echo "[NOTE] cannot detect/don't understand fastq format: $FASTQ_ENCODING - using default"
    fi
    echo "[NOTE] $FASTQ_ENCODING fastq format detected"
fi

FASTASUFFIX=${FASTA##*.}

FASTQTMP=$TMP/`whoami`_${n/%$READONE.$FASTQ/} # tmp dir to deposit fastq intermediate files
mkdir -p $FASTQTMP

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
CHECKPOINT="generating the index files"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

    if [ ! -e ${FASTA/.${FASTASUFFIX}/}.1.ebwt ]; then echo ">>>>> make .ebwt"; bowtie-build $FASTA ${FASTA/.${FASTASUFFIX}/}; fi
    if [ ! -e $FASTA.fai ]; then echo "[NOTE] make .fai"; samtools faidx $FASTA; fi

    # mark checkpoint
    if [ -f ${FASTA/.${FASTASUFFIX}/}.1.ebwt ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
fi 

################################################################################
CHECKPOINT="bowtie"

if [ $PAIRED == "0" ]; then 
    READS="$f"
    let FASTQREADS=`$ZCAT $f | wc -l | gawk '{print int($1/4)}' `
else 

    READS="-1 $f -2 ${f/$READONE/$READTWO}"
    READ1=`$ZCAT $f | wc -l | gawk '{print int($1/4)}' `
    READ2=`$ZCAT ${f/$READONE/$READTWO} | wc -l | gawk '{print int($1/4)}' `
    let FASTQREADS=$READ1+$READ2
fi

#readgroup
FULLSAMPLEID=$SAMPLEID"${n/%$READONE.$FASTQ/}"
RG="--sam-RG \"ID:$EXPID\" --sam-RG \"SM:$FULLSAMPLEID\" --sam-RG \"LB:$LIBRARY\" --sam-RG \"PL:$PLATFORM\""


if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

	# Unpaired
	if [ $PAIRED == "0" ]; then
        echo "[NOTE] SINGLE READS"

        RUN_COMMAND="$ZCAT $f | bowtie $RG $BOWTIEADDPARAM $FASTQ_PHRED --threads $CPU_BOWTIE --un $FASTQTMP/${n/%$READONE.$FASTQ/.$UNM.fq} --max $FASTQTMP/${n/%$READONE.$FASTQ/.$MUL.fq} --sam $BOWTIE_OPTIONS ${FASTA/.${FASTASUFFIX}/} - $FASTQTMP/${n/%$READONE.$FASTQ/.$ALN.sam}"

	#Paired
    else
        echo "[NOTE] PAIRED READS"
        # clever use of named pipes to avoid fastq unzipping
        [ -e $OUTDIR/${n}_pipe ] && rm $OUTDIR/${n}_pipe
        [ -e $OUTDIR/${n/$READONE/$READTWO}_pipe ] && rm $OUTDIR/${n/$READONE/$READTWO}_pipe
        mkfifo $OUTDIR/${n}_pipe $OUTDIR/${n/$READONE/$READTWO}_pipe
        
        $ZCAT $f > $OUTDIR/${n}_pipe &
        $ZCAT ${f/$READONE/$READTWO} > $OUTDIR/${n/$READONE/$READTWO}_pipe &
        
		RUN_COMMAND="bowtie $RG $BOWTIEADDPARAM $FASTQ_PHRED --threads $CPU_BOWTIE --un $FASTQTMP/${n/%$READONE.$FASTQ/.$UNM.fq} --max $FASTQTMP/${n/%$READONE.$FASTQ/.$MUL.fq} --sam $BOWTIE_OPTIONS ${FASTA/.${FASTASUFFIX}/} -1 $OUTDIR/${n}_pipe -2 $OUTDIR/${n/$READONE/$READTWO}_pipe $FASTQTMP/${n/%$READONE.$FASTQ/.$ALN.sam}"

    fi
    echo $RUN_COMMAND && eval $RUN_COMMAND
    
    if [ -e $FASTQTMP/${n/%$READONE.$FASTQ/.$MUL*.fq} ] || [ -e $FASTQTMP/${n/%$READONE.$FASTQ/.$UNM*.fq} ]; then
        $GZIP $FASTQTMP/${n/%$READONE.$FASTQ/.$MUL*.fq} $FASTQTMP/${n/%$READONE.$FASTQ/.$UNM*.fq}
    fi
    
    # bam file conversion                                                                         
    samtools view -@ $CPU_BOWTIE -Sbt $FASTA.fai $FASTQTMP/${n/%$READONE.$FASTQ/.$ALN.sam} > $OUTDIR/${n/%$READONE.$FASTQ/.$ALN.bam}
    
    # cleanup
    [ -e $OUTDIR/${n}_pipe ] && rm $OUTDIR/${n}_pipe
    [ -e $OUTDIR/${n/$READONE/$READTWO}_pipe ] && rm $OUTDIR/${n/$READONE/$READTWO}_pipe
    [ -e $FASTQTMP/${n/%$READONE.$FASTQ/.$ALN.sam} ] && rm $FASTQTMP/${n/%$READONE.$FASTQ/.$ALN.sam}

    # mark checkpoint
    if [ -f $OUTDIR/${n/%$READONE.$FASTQ/.$ALN.bam} ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
fi

################################################################################
CHECKPOINT="bam conversion and sorting"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

    # create bam files for discarded reads and remove fastq files    
    if [ $PAIRED == "1" ]; then
        if [ -e $FASTQTMP/${n/%$READONE.$FASTQ/.${UNM}_1.fq.gz} ]; then
            java $JAVAPARAMS -jar $PATH_PICARD/FastqToSam.jar \
                FASTQ=$FASTQTMP/${n/%$READONE.$FASTQ/.${UNM}_1.fq.gz} \
                FASTQ2=$FASTQTMP/${n/%$READONE.$FASTQ/.${UNM}_2.fq.gz} \
                OUTPUT=$OUTDIR/${n/%$READONE.$FASTQ/.$UNM.bam} \
                QUALITY_FORMAT=Standard \
                SAMPLE_NAME=${n/%$READONE.$FASTQ/} \
                READ_GROUP_NAME=null \
                QUIET=TRUE \
                VERBOSITY=ERROR
            samtools sort -@ $CPU_BOWTIE $OUTDIR/${n/%$READONE.$FASTQ/.$UNM.bam} $OUTDIR/${n/%$READONE.$FASTQ/.$UNM.tmp}
            mv $OUTDIR/${n/%$READONE.$FASTQ/.$UNM.tmp.bam} $OUTDIR/${n/%$READONE.$FASTQ/.$UNM.bam}
        fi
        
        if [ -e $FASTQTMP/${n/%$READONE.$FASTQ/.${MUL}_1.fq.gz} ]; then
            java $JAVAPARAMS -jar $PATH_PICARD/FastqToSam.jar \
                FASTQ=$FASTQTMP/${n/%$READONE.$FASTQ/.${MUL}_1.fq.gz} \
                FASTQ2=$FASTQTMP/${n/%$READONE.$FASTQ/.${MUL}_2.fq.gz} \
                OUTPUT=$OUTDIR/${n/%$READONE.$FASTQ/.${MUL}.bam} \
                QUALITY_FORMAT=Standard \
                SAMPLE_NAME=${n/%$READONE.$FASTQ/} \
                READ_GROUP_NAME=null \
                QUIET=TRUE \
                VERBOSITY=ERROR
        
            samtools sort -@ $CPU_BOWTIE $OUTDIR/${n/%$READONE.$FASTQ/.$MUL.bam} $OUTDIR/${n/%$READONE.$FASTQ/.$MUL.tmp}
            mv $OUTDIR/${n/%$READONE.$FASTQ/.$MUL.tmp.bam} $OUTDIR/${n/%$READONE.$FASTQ/.$MUL.bam}
        fi
    else
        if [ -e $FASTQTMP/${n/%$READONE.$FASTQ/.$UNM.fq.gz} ]; then
            java $JAVAPARAMS -jar $PATH_PICARD/FastqToSam.jar \
                FASTQ=$FASTQTMP/${n/%$READONE.$FASTQ/.$UNM.fq.gz} \
                OUTPUT=$OUTDIR/${n/%$READONE.$FASTQ/.$UNM.bam} \
                QUALITY_FORMAT=Standard \
                SAMPLE_NAME=${n/%$READONE.$FASTQ/} \
                READ_GROUP_NAME=null \
                QUIET=TRUE \
                VERBOSITY=ERROR
            samtools sort -@ $CPU_BOWTIE $OUTDIR/${n/%$READONE.$FASTQ/.$UNM.bam} $OUTDIR/${n/%$READONE.$FASTQ/.$UNM.tmp}
            mv $OUTDIR/${n/%$READONE.$FASTQ/.$UNM.tmp.bam} $OUTDIR/${n/%$READONE.$FASTQ/.$UNM.bam}
        else
            echo "[NOTE] no unmapped reads: $FASTQTMP/${n/%$READONE.$FASTQ/.$UNM.fq.gz}"
        fi
    
        if [ -e $FASTQTMP/${n/%$READONE.$FASTQ/.$MUL.fq.gz} ]; then
            java $JAVAPARAMS -jar $PATH_PICARD/FastqToSam.jar \
                FASTQ=$FASTQTMP/${n/%$READONE.$FASTQ/.$MUL.fq.gz} \
                OUTPUT=$OUTDIR/${n/%$READONE.$FASTQ/.$MUL.bam} \
                QUALITY_FORMAT=Standard \
                SAMPLE_NAME=${n/%$READONE.$FASTQ/} \
                READ_GROUP_NAME=null \
                QUIET=TRUE \
                VERBOSITY=ERROR
        
            samtools sort -@ $CPU_BOWTIE $OUTDIR/${n/%$READONE.$FASTQ/.$MUL.bam} $OUTDIR/${n/%$READONE.$FASTQ/.$MUL.tmp}
            mv $OUTDIR/${n/%$READONE.$FASTQ/.$MUL.tmp.bam} $OUTDIR/${n/%$READONE.$FASTQ/.$MUL.bam} 
        else
            echo "[NOTE] no multi-mapping reads: $FASTQTMP/${n/%$READONE.$FASTQ/.$MUL.fq.gz}"

        fi
    fi
            
    if [ "$PAIRED" = "1" ]; then
        # fix mates
        samtools sort -@ $CPU_BOWTIE -n $OUTDIR/${n/%$READONE.$FASTQ/.$ALN.bam} $OUTDIR/${n/%$READONE.$FASTQ/.tmp}
        samtools fixmate $OUTDIR/${n/%$READONE.$FASTQ/.tmp.bam} $OUTDIR/${n/%$READONE.$FASTQ/.$ALN.bam}
        [ -e $OUTDIR/${n/%$READONE.$FASTQ/.tmp.bam} ] && rm $OUTDIR/${n/%$READONE.$FASTQ/.tmp.bam}
    fi
    
    samtools sort -@ $CPU_BOWTIE $OUTDIR/${n/%$READONE.$FASTQ/.$ALN.bam} $OUTDIR/${n/%$READONE.$FASTQ/.ash}

    # mark checkpoint
    if [ -f $OUTDIR/${n/%$READONE.$FASTQ/.ash.bam} ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
    
    # cleanup
    [ -e $FASTQTMP ] && rm -r $FASTQTMP
    [ -e $OUTDIR/${n/%$READONE.$FASTQ/.$UNM.bam} ] && rm $OUTDIR/${n/%$READONE.$FASTQ/.$UNM.bam}
    [ -e $OUTDIR/${n/%$READONE.$FASTQ/.$ALN.bam} ] && rm $OUTDIR/${n/%$READONE.$FASTQ/.$ALN.bam}
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
        METRICS_FILE=$OUTDIR/metrices/${n/%$READONE.$FASTQ/.$ASD.bam}.dupl \
        AS=true \
        VALIDATION_STRINGENCY=LENIENT \
        TMP_DIR=$THISTMP
    [ -d $THISTMP ] && rm -r $THISTMP
    samtools index $OUTDIR/${n/%$READONE.$FASTQ/.$ASD.bam}

    # mark checkpoint
    if [ -f $OUTDIR/${n/%$READONE.$FASTQ/.$ASD.bam} ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
    
    # cleanup
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
        echo $(samtools view -@ $CPU_BOWTIE -c -F 4 $OUTDIR/${n/%$READONE.$FASTQ/.$ASD.bam} $SEQREG )" total reads in region " >> $STATSOUT
        echo $(samtools view -@ $CPU_BOWTIE -c -f 3 $OUTDIR/${n/%$READONE.$FASTQ/.$ASD.bam} $SEQREG )" properly paired reads in region " >> $STATSOUT
    fi

    # mark checkpoint
    if [ -f $STATSOUT ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
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
      
    # create pdfs
    if [ -n "$(which convert)" ]; then 
        for im in $( ls $OUTDIR/metrices/*.pdf ); do
            convert $im ${im/pdf/jpg}
        done
    fi
    [ -e $THISTMP ] && rm -r $THISTMP

    # mark checkpoint
    [ -f $OUTDIR/metrices/${n/%$READONE.$FASTQ/.$ASD.bam}.alignment_summary_metrics ] && echo -e "\n********* $CHECKPOINT\n" && unset RECOVERFROM
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
echo ">>>>> read mapping with bowtie 1 - FINISHED"
echo ">>>>> enddate "`date`
