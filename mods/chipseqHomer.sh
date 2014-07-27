#!/bin/bash -e

echo ">>>>> ChIPseq analysis with Homer"
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"

function usage {
echo -e "usage: $(basename $0) -k NGSANE -f FASTQ -o OUTDIR [OPTIONS]"
exit
}

# Script for ChIP-seq peak calling using Homer.
# It takes read alignments in .bam format.
# It produces output files: peak regions in bed format
# author: Fabian Buske
# date: August 2013

# QCVARIABLES,Resource temporarily unavailable
# RESULTFILENAME <DIR>/<TASK>/<SAMPLE>.summary.txt

if [ ! $# -gt 3 ]; then usage ; fi

#INPUTS                                                                                                           
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the NGSANE repository
        -f | --bam )            shift; f=$1 ;; # bam file
        -o | --outdir )         shift; OUTDIR=$1 ;; # output dir 
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
hash module 2>/dev/null && for MODULE in $MODULE_HOMERCHIPSEQ; do module load $MODULE; done && module list 

export PATH=$PATH_HOMERCHIPSEQ:$PATH
echo "PATH=$PATH"
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--samtools    --\n "$(samtools 2>&1 | head -n 3 | tail -n-2)
[ -z "$(which samtools)" ] && echo "[ERROR] no samtools detected" && exit 1
echo -e "--R           --\n "$(R --version | head -n 3)
[ -z "$(which R)" ] && echo "[ERROR] no R detected" && exit 1
echo -e "--homer       --\n "$(which makeTagDirectory)
[ -z "$(which makeTagDirectory)" ] && echo "[ERROR] homer not detected" && exit 1
echo -e "--bedToBigBed --\n "$(bedToBigBed 2>&1 | tee | head -n 1 )
[ -z "$(which bedToBigBed)" ] && echo "[WARN] bedToBigBed not detected, cannot compress bedgraphs"

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="parameters"

# get basename of f
n=${f##*/}
SAMPLE=${n/%$ASD.bam/}

GENOME_CHROMSIZES=${FASTA%.*}.chrom.sizes
if [ ! -f $GENOME_CHROMSIZES ]; then
    echo "[WARN] GENOME_CHROMSIZES not found. No bigbeds will be generated"
else
    echo "[NOTE] Chromosome size: $GENOME_CHROMSIZES"
fi

if [ -n "$CHIPINPUT" ];then
    c=${CHIPINPUT##*/}
    INPUT=${c/$ASD.bam/}
else
    INPUT="NONE"
fi

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="recall files from tape"

if [ -n "$DMGET" ]; then
    dmget -a $f
    dmget -a ${CHIPINPUT##*/}
    dmget -a $OUTDIR/*
fi

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="create tagdirectory"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else
    
    mkdir -p $OUTDIR/${SAMPLE}_homer
    RUN_COMMAND="makeTagDirectory $OUTDIR/${SAMPLE}_homer $f $HOMER_CHIPSEQ_TAGDIR_ADDPARAM &> $OUTDIR/${SAMPLE}.tagdirIP.log"
    echo $RUN_COMMAND && eval $RUN_COMMAND
    cat $OUTDIR/${SAMPLE}.tagdirIP.log # put into qout log too

    # mark checkpoint
    if [ -e $OUTDIR/${SAMPLE}_homer/tagLengthDistribution.txt ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi

################################################################################
CHECKPOINT="find peaks"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

    RUN_COMMAND="findPeaks $OUTDIR/${SAMPLE}_homer -style $HOMER_CHIPSEQ_STYLE  $HOMER_CHIPSEQ_FINDPEAKS_ADDPARAM -o auto "
    if [ -n "$CHIPINPUT" ];then
      RUN_COMMAND="$RUN_COMMAND -i $OUT/common/$TASK_HOMERCHIPSEQ/$INPUT/"  
    fi
    RUN_COMMAND="$RUN_COMMAND &> $OUTDIR/${SAMPLE}.findpeaks.log"
    echo $RUN_COMMAND && eval $RUN_COMMAND
    cat $OUTDIR/${SAMPLE}.findpeaks.log # put into qout log too
    
    if [ "$HOMER_CHIPSEQ_STYLE" == "factor" ]; then
        pos2bed.pl $OUTDIR/${SAMPLE}_homer/peaks.txt > $OUTDIR/$SAMPLE-$INPUT.bed
    
    elif [ "$HOMER_CHIPSEQ_STYLE" == "histone" ]; then
        pos2bed.pl $OUTDIR/${SAMPLE}_homer/regions.txt > $OUTDIR/$SAMPLE-$INPUT.bed
    fi

    # make bigbed
    if hash bedToBigBed && [ -f $GENOME_CHROMSIZES ]; then
        bedtools intersect -a <(cut -f1-3,5 $OUTDIR/$SAMPLE-$INPUT.bed | grep -v "^#" | sort -k1,1 -k2,2n) -b <( awk '{OFS="\t"; print $1,1,$2}' $GENOME_CHROMSIZES ) > $OUTDIR/$SAMPLE-$INPUT.tmp
        bedToBigBed -type=bed4 $OUTDIR/$SAMPLE-$INPUT.tmp $GENOME_CHROMSIZES $OUTDIR/$SAMPLE-$INPUT.bb
        rm $OUTDIR/$SAMPLE-$INPUT.tmp
    fi
        
    # mark checkpoint
    if [ -f $OUTDIR/${SAMPLE}.findpeaks.log ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi
################################################################################
CHECKPOINT="summarize"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

    cat $OUTDIR/${SAMPLE}.tagdirIP.log $OUTDIR/${SAMPLE}.findpeaks.log > $OUTDIR/${SAMPLE}.summary.txt
    

    if [ "$HOMER_CHIPSEQ_STYLE" == "factor" ]; then
        grep "^#" $OUTDIR/${SAMPLE}_homer/peaks.txt >> $OUTDIR/${SAMPLE}.summary.txt
    
    elif [ "$HOMER_CHIPSEQ_STYLE" == "histone" ]; then
        grep "^#" $OUTDIR/${SAMPLE}_homer/regions.txt >> $OUTDIR/${SAMPLE}.summary.txt
    fi
    # mark checkpoint
    if [ -f $OUTDIR/${SAMPLE}.summary.txt ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi
################################################################################
CHECKPOINT="cleanup"

if [ -z "$HOMER_KEEPTAGDIRECTORY" ]; then
   rm -f $OUTDIR/${SAMPLE}_homer/*.tags.tsv
fi

[ -f $OUTDIR/${SAMPLE}.tagdirIP.log ] && rm $OUTDIR/${SAMPLE}.tagdirIP.log
[ -f $OUTDIR/${SAMPLE}.findpeaks.log ] && rm $OUTDIR/${SAMPLE}.findpeaks.log

echo -e "\n********* $CHECKPOINT\n"
################################################################################
[ -e $OUTDIR/${SAMPLE}.summary.txt.dummy ] && rm $OUTDIR/${SAMPLE}.summary.txt.dummy
echo ">>>>> ChIPseq analysis with Homer - FINISHED"
echo ">>>>> enddate "`date`

