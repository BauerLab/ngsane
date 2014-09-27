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

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "parameters"

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
if [ "$HOMER_CHIPSEQ_STYLE" == "factor" ]; then
    PEAKSUFFIX=".narrowPeak"
elif [ "$HOMER_CHIPSEQ_STYLE" == "histone" ]; then
    PEAKSUFFIX=".broadPeak"
fi    

if [ -z "$NGSANE_RECOVERFROM" ]; then
    [ -e $OUTDIR/$SAMPLE$PEAKSUFFIX ] && rm -rf $OUTDIR/$SAMPLE $OUTDIR/$SAMPLE$PEAKSUFFIX
fi

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "recall files from tape"

if [ -n "$DMGET" ]; then
    dmget -a $f
    dmget -a ${CHIPINPUT##*/}
    dmget -a $OUTDIR/*
fi

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "create tagdirectory"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    
    mkdir -p $OUTDIR/$SAMPLE
    RUN_COMMAND="makeTagDirectory $OUTDIR/$SAMPLE $f $HOMER_CHIPSEQ_TAGDIR_ADDPARAM &> $OUTDIR/$SAMPLE.tagdirIP.log"
    echo $RUN_COMMAND && eval $RUN_COMMAND
    cat $OUTDIR/$SAMPLE.tagdirIP.log # put into qout log too

    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE/tagLengthDistribution.txt

fi
################################################################################
NGSANE_CHECKPOINT_INIT "find peaks"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    RUN_COMMAND="findPeaks $OUTDIR/$SAMPLE -style $HOMER_CHIPSEQ_STYLE  $HOMER_CHIPSEQ_FINDPEAKS_ADDPARAM -o auto "
    if [ -n "$CHIPINPUT" ];then
      RUN_COMMAND="$RUN_COMMAND -i $OUT/common/$TASK_HOMERCHIPSEQ/$INPUT/"  
    fi
    RUN_COMMAND="$RUN_COMMAND &> $OUTDIR/$SAMPLE.findpeaks.log"
    echo $RUN_COMMAND && eval $RUN_COMMAND
    cat $OUTDIR/$SAMPLE.findpeaks.log # put into qout log tooq
    
    if [ "$HOMER_CHIPSEQ_STYLE" == "factor" ]; then
        if [ "$INPUT" == "NONE" ]; then
            grep -v "#" $OUTDIR/$SAMPLE/peaks.txt | awk '{OFS="\t";print $2,$3,$4,$1,int($8),$5,$7,".",".","-1"}' > $OUTDIR/$SAMPLE$PEAKSUFFIX
        else
            grep -v "#" $OUTDIR/$SAMPLE/peaks.txt | awk 'function neglog(v){if(v==0){return 1000}else{return -log(v)/log(10)}};{OFS="\t";print $2,$3,$4,$1,int($8),$5,$6,neglog($12),".","-1"}' > $OUTDIR/$SAMPLE$PEAKSUFFIX
        fi
    elif [ "$HOMER_CHIPSEQ_STYLE" == "histone" ]; then
        if [ "$INPUT" == "NONE" ]; then
            grep -v "#" $OUTDIR/$SAMPLE/region.txt | awk '{OFS="\t";print $2,$3,$4,$1,int($8),$5,$7,".","."}' > $OUTDIR/$SAMPLE$PEAKSUFFIX
        else
            grep -v "#" $OUTDIR/$SAMPLE/region.txt | awk 'function neglog(v){if(v==0){return 1000}else{return -log(v)/log(10)}};{OFS="\t";print $2,$3,$4,$1,int($8),$5,$6,neglog($12),".","-1"}' > $OUTDIR/$SAMPLE$PEAKSUFFIX
        fi
    fi


    # make bigbed
    if hash bedToBigBed; then
        if [ -f $GENOME_CHROMSIZES ] && [ -s $OUTDIR/$SAMPLE$PEAKSUFFIX ]; then
            bedtools intersect -a <(cut -f1-3,5 $OUTDIR/$SAMPLE$PEAKSUFFIX | grep -v "^#" | sort -k1,1 -k2,2n) -b <( awk '{OFS="\t"; print $1,1,$2}' $GENOME_CHROMSIZES ) > $OUTDIR/$SAMPLE.tmp
            bedToBigBed -type=bed4 $OUTDIR/$SAMPLE.tmp $GENOME_CHROMSIZES $OUTDIR/$SAMPLE.bb
            rm $OUTDIR/$SAMPLE.tmp
        fi
    fi
        
    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE$PEAKSUFFIX

fi
################################################################################
NGSANE_CHECKPOINT_INIT "summarize"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    cat $OUTDIR/$SAMPLE.tagdirIP.log $OUTDIR/$SAMPLE.findpeaks.log > $OUTDIR/$SAMPLE.summary.txt
    

    if [ "$HOMER_CHIPSEQ_STYLE" == "factor" ]; then
        grep "^#" $OUTDIR/$SAMPLE/peaks.txt >> $OUTDIR/$SAMPLE.summary.txt
    
    elif [ "$HOMER_CHIPSEQ_STYLE" == "histone" ]; then
        grep "^#" $OUTDIR/$SAMPLE/regions.txt >> $OUTDIR/$SAMPLE.summary.txt
    fi
    
    cat > $OUTDIR/$SAMPLE/tagAutocorrelation.R <<EOF
        library(ggplot2)
        library(reshape2)
        tagAutocorrelation <- read.delim("$OUTDIR/$SAMPLE/tagAutocorrelation.txt")
        names(tagAutocorrelation)=c("x","watson","crick")
        t<-melt(tagAutocorrelation, c("x"), c("watson","crick"))
        pdf("$OUTDIR/$SAMPLE/tagAutocorrelation.pdf", width=7, height=4)
        ggplot(t,aes(x,value, group=variable, color=variable)) + geom_point(alpha = 0.5, size=1.5) + xlab("Position (nt)") + ylab("Tags") + ggtitle("Tag autocorrelation")  +  theme(legend.position = "none")
        dev.off()
EOF
    Rscript $OUTDIR/$SAMPLE/tagAutocorrelation.R
    
    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE.summary.txt

fi
################################################################################
NGSANE_CHECKPOINT_INIT "cleanup"

if [ -z "$HOMER_KEEPTAGDIRECTORY" ]; then
   rm -f $OUTDIR/$SAMPLE/*.tags.tsv
fi

[ -f $OUTDIR/$SAMPLE.tagdirIP.log ] && rm $OUTDIR/$SAMPLE.tagdirIP.log
[ -f $OUTDIR/$SAMPLE.findpeaks.log ] && rm $OUTDIR/$SAMPLE.findpeaks.log

NGSANE_CHECKPOINT_CHECK
################################################################################
[ -e $OUTDIR/$SAMPLE.summary.txt.dummy ] && rm $OUTDIR/$SAMPLE.summary.txt.dummy
echo ">>>>> ChIPseq analysis with Homer - FINISHED"
echo ">>>>> enddate "`date`

