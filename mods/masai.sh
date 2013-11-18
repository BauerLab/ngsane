#!/bin/bash -e

# Script to run masai.
# It takes comma-seprated list of files containing short sequence reads in fasta or fastq format and masai index files as input.
# It produces output files: read alignments in .bam format and other files.
# author: Fabian Buske
# date: Nov 2013

# QCVARIABLES,Resource temporarily unavailable
# RESULTFILENAME <DIR>/<TASK>/<SAMPLE>.$ASD.bam

echo ">>>>> readmapping with Masai "
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
for MODULE in $MODULE_MASAI; do module load $MODULE; done  # save way to load modules that itself load other modules
export PATH=$PATH_MASAI:$PATH
module list
echo "PATH=$PATH"
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)
PATH_IGVTOOLS=$(dirname $(which igvtools.jar))
PATH_PICARD=$(dirname $(which MarkDuplicates.jar))

echo "[NOTE] set java parameters"
JAVAPARAMS="-Xmx"$(python -c "print int($MEMORY_MASAI*0.8)")"g -Djava.io.tmpdir="$TMP"  -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1" 
unset _JAVA_OPTIONS
echo "JAVAPARAMS "$JAVAPARAMS

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--JAVA        --\n" $(java -Xmx200m -version 2>&1)
[ -z "$(which java)" ] && echo "[ERROR] no java detected" && exit 1
echo -e "--masai       --\n "$(masai_mapper 2>&1 | tee | grep version)
[ -z "$(which masai_mapper)" ] && echo "[ERROR] no masai detected" && exit 1
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

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="parameters"

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
if [ -z "$RECOVERFROM" ]; then
    [ -e $OUTDIR/$SAMPLE.$ASD.bam ] && rm $OUTDIR/$SAMPLE.$ASD.bam
    [ -e $OUTDIR/$SAMPLE.$ASD.bam.stats ] && rm $OUTDIR/$SAMPLE.$ASD.bam.stats
    [ -e $OUTDIR/$SAMPLE.$ASD.bam.dupl ] && rm $OUTDIR/$SAMPLE.$ASD.bam.dupl
fi

#is ziped ?
ZCAT="zcat"
if [[ ${f##*.} != "gz" ]]; then ZCAT="cat"; fi

#is paired ?                                                                                                      
if [ "$f" != "${f/%$READONE.$FASTQ/$READTWO.$FASTQ}" ] && [ -e ${f/%$READONE.$FASTQ/$READTWO.$FASTQ} ]; then
    echo "[NOTE] paired-end library detected"
    PAIRED="1"
    READ1=`$ZCAT $f | wc -l | gawk '{print int($1/4)}' `
    READ2=`$ZCAT ${f/%$READONE.$FASTQ/$READTWO.$FASTQ} | wc -l | gawk '{print int($1/4)}' `
    let FASTQREADS=$READ1+$READ2
else
    echo "[NOTE] single-end library detected"
    PAIRED="0"
    let FASTQREADS=`$ZCAT $f | wc -l | gawk '{print int($1/4)}' `
fi

if [ -z "$MASAI_INDEX" ]; then
    echo "[NOTE] default index for masai is sa"
    MASAI_INDEX="sa"
fi

THISTMP=$TMP"/"$(whoami)"/"$(echo $OUTDIR | md5sum | cut -d' ' -f1)
mkdir -p $THISTMP

FASTASUFFIX=${FASTA##*.}
    
#readgroup
#TODO check readgroups
FULLSAMPLEID=$SAMPLEID"${INPUTFILENAME/%$READONE.$FASTQ/}"
RG="--sam-rg \"ID:$EXPID\" --sam-rg \"SM:$FULLSAMPLEID\" --sam-rg \"LB:$LIBRARY\" --sam-rg \"PL:$PLATFORM\""

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="recall files from tape"
	
if [ -n "$DMGET" ]; then
    dmget -a $FASTA*
    dmget -a ${f/%$READONE.$FASTQ/"*"}
    dmget -a ${OUTDIR}
fi
    
echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="generating the index files"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
    
    if [[ "$MASAI_INDEX" = "sa" ]] && [[ ! -e ${FASTA/.$FASTASUFFIX/.sa} ]]; then 
        echo "[NOTE] make index: $MASAI_INDEX"; 
        # make sure parallel executions don't interfere with each other when generating the index.
        mkdir -p $THISTMP/reference/
        ln -s $FASTA $THISTMP/reference/
        FASTAREFERENCE=$THISTMP/reference/${FASTA##*/}
        masai_indexer -t $TMP -x $MASAI_INDEX $FASTAREFERENCE

    elif [[ "$MASAI_INDEX" = "fm" ]] && [[ ! -e ${FASTA/.$FASTASUFFIX/.fma} ]]; then
        echo "[NOTE] make index: $MASAI_INDEX";
        # make sure parallel executions don't interfere with each other when generating the index.
        mkdir -p $THISTMP/reference/
        ln -s $FASTA $THISTMP/reference/
        FASTAREFERENCE=$THISTMP/reference/${FASTA##*/}
        masai_indexer -t $TMP -x $MASAI_INDEX $FASTAREFERENCE

    elif [[ "$MASAI_INDEX" = "esa" ]] && [[ ! -e ${FASTA/.$FASTASUFFIX/.esa} ]]; then
        echo "[NOTE] make index: $MASAI_INDEX";
        # make sure parallel executions don't interfere with each other when generating the index.
        mkdir -p $THISTMP/reference/
        ln -s $FASTA $THISTMP/reference/
        FASTAREFERENCE=$THISTMP/reference/${FASTA##*/}
        masai_indexer -t $TMP -x $MASAI_INDEX $FASTAREFERENCE

    else
        FASTAREFERENCE=$FASTA
    fi

    if [ ! -e $FASTA.fai ]; then echo "[NOTE] make .fai"; samtools faidx $FASTA; fi

    # mark checkpoint
    if [ -f $FASTA.fai ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
fi 

################################################################################
CHECKPOINT="masai"
if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
        
    if [ "$PAIRED" = "0" ]; then
        
        $ZCAT $f > $THISTMP/${INPUTFILENAME/%.$FASTQ/_pipe.fastq}

        RUN_COMMAND="masai_mapper $MASAI_MAPPERADDPARAM --index $MASAI_INDEX --output-file $THISTMP/${INPUTFILENAME/%.$FASTQ/.raw} $FASTAREFERENCE $THISTMP/${INPUTFILENAME/%.$FASTQ/_pipe.fastq}"
        echo $RUN_COMMAND && eval $RUN_COMMAND
    
        RUN_COMMAND="masai_output_se $MASAI_OUTPUTADDPARAM --tmp-folder $TMP --output-file $THISTMP/$SAMPLE.$ALN.sam $FASTAREFERENCE $THISTMP/${INPUTFILENAME/%.$FASTQ/_pipe.fastq} $THISTMP/${INPUTFILENAME/%.$FASTQ/.raw}"
        echo $RUN_COMMAND && eval $RUN_COMMAND
        
    elif [ "$PAIRED" = "1" ]; then

        $ZCAT $f > $THISTMP/${INPUTFILENAME/%.$FASTQ/_pipe.fastq}
        $ZCAT ${f/%$READONE.$FASTQ/$READTWO.$FASTQ} > $THISTMP/${INPUTFILENAME/%$READONE.$FASTQ/${READTWO}_pipe.fastq}
      
        RUN_COMMAND="masai_mapper $MASAI_MAPPERADDPARAM --index $MASAI_INDEX --output-file $THISTMP/${INPUTFILENAME/%.$FASTQ/.raw} $FASTAREFERENCE $THISTMP/${INPUTFILENAME/%.$FASTQ/_pipe.fastq}"
        echo $RUN_COMMAND && eval $RUN_COMMAND

        RUN_COMMAND="masai_mapper $MASAI_MAPPERADDPARAM --index $MASAI_INDEX --output-file $THISTMP/${INPUTFILENAME/%$READONE.$FASTQ/$READTWO.raw} $FASTAREFERENCE $THISTMP/${INPUTFILENAME/%$READONE.$FASTQ/${READTWO}_pipe.fastq}"
        echo $RUN_COMMAND && eval $RUN_COMMAND

        $ZCAT $f > $THISTMP/${INPUTFILENAME/%.$FASTQ/_pipe.fastq}
        $ZCAT ${f/%$READONE.$FASTQ/$READTWO.$FASTQ} > $THISTMP/${INPUTFILENAME/%$READONE.$FASTQ/${READTWO}_pipe.fastq}
        
        RUN_COMMAND="masai_output_pe $MASAI_OUTPUTADDPARAM --tmp-folder $TMP --output-file $THISTMP/$SAMPLE.$ALN.sam $FASTAREFERENCE $THISTMP/${INPUTFILENAME/%.$FASTQ/_pipe.fastq} $THISTMP/${INPUTFILENAME/%$READONE.$FASTQ/${READTWO}_pipe.fastq} $THISTMP/${INPUTFILENAME/%.$FASTQ/.raw} $THISTMP/${INPUTFILENAME/%$READONE.$FASTQ/$READTWO.raw}"
        echo $RUN_COMMAND && eval $RUN_COMMAND

    fi

    # bam file conversion                                                                         
    samtools view -@ $CPU_MASAI -Sbt $FASTA.fai $THISTMP/$SAMPLE.$ALN.sam > $OUTDIR/$SAMPLE.$ALN.bam
    samtools sort -@ $CPU_MASAI $OUTDIR/$SAMPLE.$ALN.bam $OUTDIR/$SAMPLE.ash

    # cleanup
    [ -e $THISTMP/${INPUTFILENAME}_pipe ] && rm $THISTMP/${INPUTFILENAME}_pipe
    [ -e $THISTMP/${INPUTFILENAME/%$READONE.$FASTQ/$READTWO.$FASTQ}_pipe ] && rm $THISTMP/${INPUTFILENAME/%$READONE.$FASTQ/$READTWO.$FASTQ}_pipe
    [ -e $THISTMP/$SAMPLE.$ALN.sam ] && rm $THISTMP/$SAMPLE.$ALN.sam
    [ -e $THISTMP/${INPUTFILENAME/%.$FASTQ/.raw} ] && rm $THISTMP/${INPUTFILENAME/%.$FASTQ/.raw}
    [ -e $THISTMP/${INPUTFILENAME/%$READONE.$FASTQ/$READTWO.raw} ] && rm $THISTMP/${INPUTFILENAME/%$READONE.$FASTQ/$READTWO.raw}
    [ -e $OUTDIR/$SAMPLE.$ALN.bam ] && rm $OUTDIR/$SAMPLE.$ALN.bam
       
    # mark checkpoint
    if [ -f $OUTDIR/$SAMPLE.ash.bam ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
    
fi

################################################################################
CHECKPOINT="add readgroup"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
    
    java $JAVAPARAMS -jar $PATH_PICARD/AddOrReplaceReadGroups.jar \
        INPUT=$OUTDIR/$SAMPLE.ash.bam \
        OUTPUT=$OUTDIR/$SAMPLE.ashrg.bam \
        RGID=$EXPID RGLB=$LIBRARY RGPL=$PLATFORM \
        RGSM=$FULLSAMPLEID RGPU="XXXXXX"

    # mark checkpoint
    if [ -f $OUTDIR/$SAMPLE.ashrg.bam ];then echo -e "\n********* $CHECKPOINT\n" && unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
    
    [ -e $OUTDIR/$SAMPLE.ash.bam ] && rm $OUTDIR/$SAMPLE.ash.bam
fi 

################################################################################
CHECKPOINT="mark duplicates"
# create bam files for discarded reads and remove fastq files
if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
   
    if [ ! -e $OUTDIR/metrices ]; then mkdir -p $OUTDIR/metrices ; fi
    java $JAVAPARAMS -jar $PATH_PICARD/MarkDuplicates.jar \
        INPUT=$OUTDIR/$SAMPLE.ashrg.bam \
        OUTPUT=$OUTDIR/$SAMPLE.$ASD.bam \
        METRICS_FILE=$OUTDIR/metrices/$SAMPLE.$ASD.bam.dupl AS=true \
        VALIDATION_STRINGENCY=LENIENT \
        TMP_DIR=$THISTMP
        
    samtools index $OUTDIR/$SAMPLE.$ASD.bam
          
    # mark checkpoint
    if [ -f $OUTDIR/$SAMPLE.$ASD.bam ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
    
    #cleanup
    [ -e $OUTDIR/$SAMPLE.ashrg.bam ] && rm $OUTDIR/$SAMPLE.ashrg.bam
    
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
        echo $(samtools view -@ $CPU_MASAI -c -F 4 $OUTDIR/$SAMPLE.$ASD.bam $SEQREG )" total reads in region " >> $STATSOUT
        echo $(samtools view -@ $CPU_MASAI -c -f 3 $OUTDIR/$SAMPLE.$ASD.bam $SEQREG )" properly paired reads in region " >> $STATSOUT
    fi

    # mark checkpoint
    if [ -e $STATSOUT ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
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
    for im in $( ls $OUTDIR/metrices/*.pdf ); do
        convert $im ${im/pdf/jpg}
    done

    # mark checkpoint
    if [ -f $OUTDIR/metrices/$SAMPLE.$ASD.bam.alignment_summary_metrics ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
fi

################################################################################
CHECKPOINT="coverage track"    

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
    
    java $JAVAPARAMS -jar $PATH_IGVTOOLS/igvtools.jar count $OUTDIR/$SAMPLE.$ASD.bam $OUTDIR/$SAMPLE.$ASD.bam.cov.tdf ${FASTA/.$FASTASUFFIX/}.genome
    # mark checkpoint
    [ -f $OUTDIR/$SAMPLE.$ASD.bam.cov.tdf ] && echo -e "\n********* $CHECKPOINT\n" && unset RECOVERFROM
fi

################################################################################
CHECKPOINT="samstat"    

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
    
    samstat $OUTDIR/$SAMPLE.$ASD.bam

    # mark checkpoint
    if [ -f $OUTDIR/$SAMPLE.$ASD.bam.stats ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
    
fi

###############################################################################
CHECKPOINT="verify"    

BAMREADS=$(head -n1 $OUTDIR/$SAMPLE.$ASD.bam.stats | cut -d " " -f 1)
if [ "$BAMREADS" = "" ]; then let BAMREADS="0"; fi
if [ $BAMREADS -eq $FASTQREADS ]; then
    echo "[NOTE] PASS check mapping: $BAMREADS == $FASTQREADS"
else
    echo -e "[ERROR] We are loosing reads from .fastq -> .bam in $f: \nFastq had $FASTQREADS Bam has $BAMREADS"
    exit 1
fi

echo -e "\n********* $CHECKPOINT\n"

###############################################################################
CHECKPOINT="cleanup"    

[ -d $THISTMP ] && rm -r $THISTMP

echo -e "\n********* $CHECKPOINT\n"
    
################################################################################
[ -e $OUTDIR/$SAMPLE.$ASD.bam.dummy ] && rm $OUTDIR/$SAMPLE.$ASD.bam.dummy
echo ">>>>> readmapping with Masai - FINISHED"
echo ">>>>> enddate "`date`
