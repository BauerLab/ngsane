#!/bin/bash -e

# Script for ChIP-seq peak calling using peakranger.
# It takes read alignments in .bam format.
# It produces output files: peak regions in bed format
# author: Fabian Buske
# date: August 2013

# QCVARIABLES,Resource temporarily unavailable
# RESULTFILENAME <DIR>/<TASK>/<SAMPLE>_region.bed

echo ">>>>> ChIPseq analysis with Peakranger"
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

#INPUTS                                                                                                           
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the NGSANE repository
        -f | --bam )            shift; f=$1 ;; # bam file
        -o | --outdir )         shift; OUTDIR=$1 ;; # output dir 
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
hash module 2>/dev/null && for MODULE in $MODULE_PEAKRANGER; do module load $MODULE; done && module list 

export PATH=$PATH_PEAKRANGER:$PATH
echo "PATH=$PATH"
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)

echo -e "--NGSANE     --\n" $(trigger.sh -v 2>&1)
echo -e "--R          --\n "$(R --version | head -n 3)
[ -z "$(which R)" ] && echo "[ERROR] no R detected" && exit 1
echo -e "--peakranger --\n "$(peakranger | head -n 3 | tail -n 1)
[ -z "$(which peakranger)" ] && echo "[ERROR] peakranger not detected" && exit 1
echo -e "--bedtools --\n "$(bedtools --version)
[ -z "$(which bedtools)" ] && echo "[ERROR] no bedtools detected" && exit 1
echo -e "--bedToBigBed --\n "$(bedToBigBed 2>&1 | tee | head -n 1 )
[ -z "$(which bedToBigBed)" ] && echo "[WARN] bedToBigBed not detected, cannot compress bedgraphs"

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "parameters"

if [ -z "$CHIPINPUT" ] || [ ! -f $CHIPINPUT ]; then
    echo "[ERROR] input control not provided or invalid (CHIPINPUT=\"$CHIPINPUT\")"
    exit 1
fi

# get basename of f
f=${f/%.dummy/}
n=${f##*/}
SAMPLE=${n/%$ASD.bam/}
c=${CHIPINPUT##*/}
CONTROL=${c/%$ASD.bam/}

GENOME_CHROMSIZES=${FASTA%.*}.chrom.sizes
if [ ! -f $GENOME_CHROMSIZES ]; then
    echo "[WARN] GENOME_CHROMSIZES not found. Excepted at $GENOME_CHROMSIZES. Will not create bigBed file"
else
    echo "[NOTE] Chromosome size: $GENOME_CHROMSIZES"
fi

if [ -z "$NGSANE_RECOVERFROM" ]; then
    [ -e $OUTDIR/${SAMPLE}_region.bed ] && rm $OUTDIR/$SAMPLE*
fi

if [ "$PEAKRANGER_PEAKS" != "broad" ] && [ "$PEAKRANGER_PEAKS" != "sharp" ]; then
    echo "[ERROR] PEAKRANGER_PEAKS parameter not valid: $PEAKRANGER_PEAKS"
    exit 1
fi

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "recall files from tape"

if [ -n "$DMGET" ]; then
	dmget -a ${f}
	dmget -a $OUTDIR/*
fi

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "assess data quality"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then
    
    [ -f $OUTDIR/$SAMPLE.summary.txt ] && rm $OUTDIR/$SAMPLE.summary.txt
    
    RUN_COMMAND="peakranger nr --format bam --data $f --control $CHIPINPUT | tr ':' '\t' > $OUTDIR/$SAMPLE.summary.txt"
    echo $RUN_COMMAND && eval $RUN_COMMAND

    RUN_COMMAND="peakranger lc --data $f | tr ':' '\t' >> $OUTDIR/$SAMPLE.summary.txt"
    echo $RUN_COMMAND && eval $RUN_COMMAND
    
    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE.summary.txt

fi
################################################################################
NGSANE_CHECKPOINT_INIT "peakranger"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    if [ "$PEAKRANGER_PEAKS" == "broad" ]; then
        echo "[NOTE] calling broad peaks"
        RUN_COMMAND="peakranger ccat $PEAKRANGERADDPARAM --format bam --data $f --control $CHIPINPUT --output $OUTDIR/$SAMPLE -t $CPU_PEAKRANGER"

    elif [ "$PEAKRANGER_PEAKS" == "sharp" ]; then
        echo "[NOTE] calling tight peaks"
        RUN_COMMAND="peakranger ranger $PEAKRANGERADDPARAM --format bam --data $f --control $CHIPINPUT --output $OUTDIR/$SAMPLE -t $CPU_PEAKRANGER"
    fi
    echo $RUN_COMMAND && eval $RUN_COMMAND
    
    # keep only peaks that pass the FDR and merge overlapping peaks keeping the smallest score
    awk '{if($0~"fdrPassed"){print $0}}' $OUTDIR/${SAMPLE}_region.bed | bedtools sort | bedtools merge -scores min -i stdin > $OUTDIR/${SAMPLE}_region.tmp
    mv $OUTDIR/${SAMPLE}_region.tmp $OUTDIR/${SAMPLE}_region.bed

    awk '{if($0~"fdrPassed"){print $0}}' $OUTDIR/$SAMPLE"_summit.bed" | bedtools sort | bedtools merge -scores min -i stdin > $OUTDIR/$SAMPLE"_summit.tmp"
    mv $OUTDIR/$SAMPLE"_summit.tmp" $OUTDIR/$SAMPLE"_summit.bed"
    
    echo "Peaks: $(wc -l $OUTDIR/${SAMPLE}_region.bed | awk '{print $1}')" >> $OUTDIR/$SAMPLE.summary.txt
    echo "Nucleotides covered: $(awk '{sum+=$3-$2}END{print sum}' $OUTDIR/${SAMPLE}_region.bed)" >> $OUTDIR/$SAMPLE.summary.txt
    echo "ChIP input: $CONTROL" >> $OUTDIR/$SAMPLE.summary.txt

    # make bigbed
	if hash bedToBigBed && [ -f $GENOME_CHROMSIZES ] && [ -s $OUTDIR/${SAMPLE}_region.bed ] ; then
        bedtools intersect -a $OUTDIR/${SAMPLE}_region.bed -b <( awk '{OFS="\t"; print $1,1,$2}' $GENOME_CHROMSIZES ) > $OUTDIR/${SAMPLE}_region.tmp
        bedToBigBed -type=bed4 $OUTDIR/${SAMPLE}_region.tmp $GENOME_CHROMSIZES $OUTDIR/$SAMPLE.bb
        rm $OUTDIR/${SAMPLE}_region.tmp
    fi
    
    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/${SAMPLE}_region.bed

fi

################################################################################
NGSANE_CHECKPOINT_INIT "zip"

$GZIP -f $OUTDIR/$SAMPLE"_details"

NGSANE_CHECKPOINT_CHECK
################################################################################
[ -e $OUTDIR/$SAMPLE"_region.bed.dummy" ] && rm $OUTDIR/$SAMPLE"_region.bed.dummy"
echo ">>>>> ChIPseq analysis with Peakranger - FINISHED"
echo ">>>>> enddate "`date`