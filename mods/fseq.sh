#!/bin/bash -e

# DNase-Seq processing using fseq
# author: Fabian Buske
# date: January 1914

# QCVARIABLES,Resource temporarily unavailable
# RESULTFILENAME <DIR>/<TASK>/<SAMPLE>.narrowPeak

echo ">>>>> Peak calling with fseq"
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"

function usage {
echo -e "usage: $(basename $0) -k NGSANE -f INPUTFILE -o OUTDIR [OPTIONS]"
exit
}

if [ ! $# -gt 3 ]; then usage ; fi

#INPUTS                                                                                                           
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;;     # location of the NGSANE repository                       
        -f | --file )           shift; INPUTFILE=$1 ;;  # input file                                                       
        -o | --outdir )         shift; OUTDIR=$1 ;;     # output dir                                                     
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
hash module 2>/dev/null && for MODULE in $MODULE_FSEQ; do module load $MODULE; done && module list 

export PATH=$PATH_FSEQ:$PATH
echo "PATH=$PATH"

echo "[NOTE] set java parameters"
JAVAPARAMS="-Xmx"$(python -c "print int($MEMORY_FSEQ*0.8)")"g -Djava.io.tmpdir="$TMP" -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1" 
unset _JAVA_OPTIONS
echo "JAVAPARAMS "$JAVAPARAMS

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--JAVA        --\n" $(java -Xmx200m -version 2>&1)
[ -z "$(which java)" ] && echo "[ERROR] no java detected" && exit 1
echo -e "--bedtools --\n "$(bedtools --version)
[ -z "$(which bedtools)" ] && echo "[ERROR] no bedtools detected" && exit 1
echo -e "--wigToBigWig --\n "$(wigToBigWig 2>&1 | tee | head -n 1)
[ -z "$(which wigToBigWig)" ] && echo "[ERROR] wigToBigWig not detected" && exit 1
echo -e "--fseq        --\n "$(fseq -v | head -n 1)
[ -z "$(which fseq)" ] && echo "[ERROR] no fseq detected" && exit 1
echo -e "--R           --\n "$(R --version | head -n 3)
[ -z "$(which R)" ] && echo "[ERROR] no R detected" && exit 1
echo -e "--bedToBigBed --\n "$(bedToBigBed 2>&1 | tee | head -n 1 )
[ -z "$(which bedToBigBed)" ] && echo "[WARN] bedToBigBed not detected, cannot compress bedgraphs"

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "parameters"

# get basename of input file f
INPUTFILENAME=${INPUTFILE##*/}
# get sample prefix
SAMPLE=${INPUTFILENAME/%$ASD.bam/}

# delete old bam files unless attempting to recover
if [ -z "$NGSANE_RECOVERFROM" ]; then
    ## TODO remove primary result files from pervious runs
fi

if [ -z "$FASTA" ] || [ ! -f $FASTA ]; then
    echo "[ERROR] no reference provided (FASTA)"
    exit 1
else
    echo "[NOTE] Reference: $FASTA"
fi

GENOME_CHROMSIZES=${FASTA%.*}.chrom.sizes
if [ ! -f $GENOME_CHROMSIZES ]; then
    echo "[ERROR] GENOME_CHROMSIZES not found. Excepted at $GENOME_CHROMSIZES"
    exit 1
else
    echo "[NOTE] Chromosome size: $GENOME_CHROMSIZES"
fi

if [ -z "$FSEQ_PVALUECUTOFF" ] || [ $FSEQ_PVALUECUTOFF -gt 1 ] ; then
    FSEQ_PVALUECUTOFF=0.05
fi

# unique temp folder that should be used to store temporary files
THISTMP=$TMP"/"$(whoami)"/"$(echo $OUTDIR | md5sum | cut -d' ' -f1)
mkdir -p $THISTMP

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "recall files from tape"

if [ -n "$DMGET" ]; then
	dmget -a $(dirname $FASTA)/*
	dmget -a $INPUTFILE
    dmget -a $OUTDIR/*
fi
    
NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "fseq bw"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    if [ -n "$FSEQ_MAKEBIGWIG" ]; then
        rm -f $THISTMP/*wig
    
        RUN_COMMAND="java $JAVAPARAMS -cp $CLASSPATH edu.duke.igsp.gkde.Main $FSEQADDPARAM -o $THISTMP -of wig <(bedtools bamtobed -i $INPUTFILE )"
        echo $RUN_COMMAND && eval $RUN_COMMAND
        
        cat $THISTMP/*.wig | wigToBigWig stdin ${GENOME_CHROMSIZES} $OUTDIR/$SAMPLE.bw
    
        # mark checkpoint
        NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE.bw
    else
        echo "[NOTE] skip bw generation"
        NGSANE_CHECKPOINT_CHECK
    fi 
fi
################################################################################
NGSANE_CHECKPOINT_INIT "fseq narrowpeak"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then
    RUN_COMMAND="java $JAVAPARAMS -cp $CLASSPATH edu.duke.igsp.gkde.Main $FSEQADDPARAM -o $THISTMP -of npf <(bedtools bamtobed -i $INPUTFILE )"
    echo $RUN_COMMAND && eval $RUN_COMMAND

    cat $THISTMP/*.npf > $OUTDIR/$SAMPLE.narrowPeak

    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE.narrowPeak
fi

################################################################################
NGSANE_CHECKPOINT_INIT "estimate p-values"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    # estimate p-value by fitting gamma distribution
    cat > $OUTDIR/$SAMPLE.R <<DELIM    
library(MASS)
x <- read.delim("$OUTDIR/$SAMPLE.narrowPeak", header=F)
g<-MASS::fitdistr(x\$V7,"gamma")
x\$V8<-sapply(x\$V7, function(y) ks.test(y, "pgamma",shape=g\$estimate["shape"],rate=g\$estimate["rate"])\$p.value)
write.table(x, "$OUTDIR/$SAMPLE.narrowPeak.tmp", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
DELIM
    
    Rscript --vanilla $OUTDIR/$SAMPLE.R

    #filter with p-value cutoff
    awk -v CUTOFF=$FSEQ_PVALUECUTOFF '{if ($8 <= CUTOFF){print $0}}' $OUTDIR/$SAMPLE.narrowPeak.tmp > $OUTDIR/$SAMPLE.narrowPeak
    [ -f $OUTDIR/$SAMPLE.narrowPeak.tmp ] && rm $OUTDIR/$SAMPLE.narrowPeak.tmp 
    [ -f $OUTDIR/$SAMPLE.R ] && rm $OUTDIR/$SAMPLE.R

    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE.narrowPeak
fi

################################################################################
NGSANE_CHECKPOINT_INIT "generate bigbed"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

	if hash bedToBigBed ; then 
        echo "[NOTE] create bigbed from peaks" 
        awk '{OFS="\t"; print $1,$2,$3,$7}' $OUTDIR/$SAMPLE.narrowPeak > $OUTDIR/$SAMPLE.peak.tmp
        bedToBigBed -type=bed4 $OUTDIR/$SAMPLE.peak.tmp $GENOME_CHROMSIZES $OUTDIR/$SAMPLE.bb
        rm $OUTDIR/$SAMPLE.peak.tmp
         # mark checkpoint
        NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE.bb
    else
        echo "[NOTE] bigbed not generated"
        NGSANE_CHECKPOINT_CHECK
    fi
fi   
###############################################################################
NGSANE_CHECKPOINT_INIT "cleanup"

[ -d $THISTMP ] && rm -r $THISTMP

NGSANE_CHECKPOINT_CHECK
################################################################################
[ -e $OUTDIR/$SAMPLE.narrowPeak.dummy ] && rm $OUTDIR/$SAMPLE.narrowPeak.dummy
echo ">>>>> Peak calling with fseq - FINISHED"
echo ">>>>> enddate "`date`
