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

# get basename of f
# get basename of input file f
INPUTFILENAME=${INPUTFILE##*/}
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
    PAIRED="1"
    READS="-1 $f -2 ${f/%$READONE.$FASTQ/$READTWO.$FASTQ}"
    READ1=`$ZCAT $f | wc -l | gawk '{print int($1/4)}' `
    READ2=`$ZCAT ${f/%$READONE.$FASTQ/$READTWO.$FASTQ} | wc -l | gawk '{print int($1/4)}' `
    let FASTQREADS=$READ1+$READ2
else
    PAIRED="0"
    READS="-U $f"
    let FASTQREADS=`$ZCAT $f | wc -l | gawk '{print int($1/4)}' `
fi

THISTMP=$TMP"/"$(whoami)"/"$(echo $OUTDIR | md5sum | cut -d' ' -f1)
mkdir -p $THISTMP

#readgroup
#TODO check readgroups
#FULLSAMPLEID=$SAMPLEID"${n/%$READONE.$FASTQ/}"
#RG="--sam-rg \"ID:$EXPID\" --sam-rg \"SM:$FULLSAMPLEID\" --sam-rg \"LB:$LIBRARY\" --sam-rg \"PL:$PLATFORM\""

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
    
    FASTASUFFIX=${FASTA##*.}
    if [ ! -e ${FASTA/.${FASTASUFFIX}/}.1.bt2 ]; then echo "[NOTE] make index .$MASAI_INDEX"; 
        masai_indexer -t $TMP -x $MASAI_INDEX $FASTA
    fi
    if [ ! -e $FASTA.fai ]; then echo "[NOTE] make .fai"; samtools faidx $FASTA; fi

    # mark checkpoint
    if [ -f ${FASTA/.$FASTASUFFIX/.$MASAI_INDEX} ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
fi 

################################################################################
CHECKPOINT="masai"
if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
        
    if [ "$PAIRED" = "0" ]; then
        # avoid IO by using named pipes
        mkfifo $THISTMP/${n}_pipe
        $ZCAT $f > $THISTMP/${n}_pipe &

        RUN_COMMAND="masai_mapper $MASAI_MAPPERADDPARAM --index $MASAI_INDEX --output-format raw --output-file $THISTMP/${n/%.$FASTQ/.raw} $FASTA $THISTMP/${n}_pipe"
        echo $RUN_COMMAND && eval $RUN_COMMAND
    
        RUN_COMMAND="masai_output_se $MASAI_OUTPUTADDPARAM --tmp-folder $TMP --output-format sam --output-file $THISTMP/$SAMPLE.$ALN.sam $FASTA $f $THISTMP/${n/%.$FASTQ/.raw}"
        echo $RUN_COMMAND && eval $RUN_COMMAND
        
    else if [ "$PAIRED" = "1" ]; then
        # avoid IO by using named pipes
        mkfifo $THISTMP/${n}_pipe $THISTMP/${n/%$READONE.$FASTQ/$READTWO.$FASTQ}_pipe        
        $ZCAT $f > $THISTMP/${n}_pipe &
        $ZCAT ${f/%$READONE.$FASTQ/$READTWO.$FASTQ} > $THISTMP/${n/%$READONE.$FASTQ/$READTWO.$FASTQ}_pipe &
      
        RUN_COMMAND="masai_mapper $MASAI_MAPPERADDPARAM --index $MASAI_INDEX --output-format raw --output-file $THISTMP/${n/%.$FASTQ/.raw} $FASTA $THISTMP/${n}_pipe"
        echo $RUN_COMMAND && eval $RUN_COMMAND

        RUN_COMMAND="masai_mapper $MASAI_MAPPERADDPARAM --index $MASAI_INDEX --output-format raw --output-file $THISTMP/${n/%$READONE.$FASTQ/$READTWO.raw} $FASTA $THISTMP/${n/%$READONE.$FASTQ/$READTWO.$FASTQ}_pipe"
        echo $RUN_COMMAND && eval $RUN_COMMAND
        
        RUN_COMMAND="masai_output_pe $MASAI_OUTPUTADDPARAM --tmp-folder $TMP --output-format sam --output-file $THISTMP/$SAMPLE.$ALN.sam $FASTA $f ${f/%$READONE.$FASTQ/$READTWO.$FASTQ} $THISTMP/${n/%.$FASTQ/.raw} $THISTMP/${n/%$READONE.$FASTQ/$READTWO.raw}"
        echo $RUN_COMMAND && eval $RUN_COMMAND

    fi

    # bam file conversion                                                                         
    samtools view -@ $CPU_MASAI -Sbt $FASTA.fai $THISTMP/$SAMPLE.$ALN.sam > $OUTDIR/$SAMPLE.$ALN.bam
    samtools sort -@ $CPU_MASAI $OUTDIR/$SAMPLE.$ALN.bam $OUTDIR/$SAMPLE.ash

    # cleanup
    [ -e $THISTMP/${n}_pipe ] && rm $THISTMP/${n}_pipe
    [ -e $THISTMP/${n/%$READONE.$FASTQ/$READTWO.$FASTQ}_pipe ] && rm $THISTMP/${n/%$READONE.$FASTQ/$READTWO.$FASTQ}_pipe
    [ -e $THISTMP/$SAMPLE.$ALN.sam ] && rm $THISTMP/$SAMPLE.$ALN.sam
    [ -e $THISTMP/${n/%.$FASTQ/.raw} ] && rm $THISTMP/${n/%.$FASTQ/.raw}
    [ -e $THISTMP/${n/%$READONE.$FASTQ/$READTWO.raw} ] && rm $THISTMP/${n/%$READONE.$FASTQ/$READTWO.raw}
    [ -e $OUTDIR/$SAMPLE.$ALN.bam ] && rm $OUTDIR/$SAMPLE.$ALN.bam
       
    # mark checkpoint
    if [ -f $OUTDIR/$SAMPLE.$ALN.bam ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
    
fi

################################################################################
CHECKPOINT="mark duplicates"
# create bam files for discarded reads and remove fastq files
if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
   
    if [ ! -e $OUTDIR/metrices ]; then mkdir -p $OUTDIR/metrices ; fi
    java $JAVAPARAMS -jar $PATH_PICARD/MarkDuplicates.jar \
        INPUT=$OUTDIR/${n/%$READONE.$FASTQ/.ash.bam} \
        OUTPUT=$OUTDIR/$SAMPLE.$ASD.bam \
        METRICS_FILE=$OUTDIR/metrices/$SAMPLE.$ASD.bam.dupl AS=true \
        VALIDATION_STRINGENCY=LENIENT \
        TMP_DIR=$THISTMP
    samtools index $OUTDIR/$SAMPLE.$ASD.bam
          
    # mark checkpoint
    if [ -f $OUTDIR/$SAMPLE.$ASD.bam ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
    
    #cleanup
    [ -e $OUTDIR/${n/%$READONE.$FASTQ/.ash.bam} ] && rm $OUTDIR/${n/%$READONE.$FASTQ/.ash.bam}
    
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
    
    java $JAVAPARAMS -jar $PATH_IGVTOOLS/igvtools.jar count $OUTDIR/$SAMPLE.$ASD.bam $OUTDIR/${n/%$READONE.$FASTQ/.$ASD.bam.cov.tdf} ${FASTA/.$FASTASUFFIX/}.genome
    # mark checkpoint
    [ -f $OUTDIR/${n/%$READONE.$FASTQ/.$ASD.bam.cov.tdf} ] && echo -e "\n********* $CHECKPOINT\n" && unset RECOVERFROM
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
