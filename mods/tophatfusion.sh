#!/bin/bash -e

# Script to run TopHat --fusion search
# Fabian Buske and Hugh French
# date: March 2014

# RESULTFILENAME <DIR>/<TASK>/<SAMPLE>

echo ">>>>> readmapping with Tophat "
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"


function usage {
echo -e "usage: $(basename $0) -k NGSANE -f FASTQ -o OUTDIR [OPTIONS]

Script running read mapping for single and paired DNA reads from fastq files
It expects a fastq file, pairdend, reference genome  as input and 
It runs tophat, converts the output to .bam files, adds header information.

required:
  -k | --toolkit <path>     location of the NGSANE repository 
  -f | --fastq <file>       fastq file
  -o | --outdir <path>      output dir
"
exit
}

if [ ! $# -gt 3 ]; then usage ; fi

#DEFAULTS
FORCESINGLE=0

#INPUTS
while [ "$1" != "" ]; do
    case $1 in
    -k | toolkit )          shift; CONFIG=$1 ;; # ENSURE NO VARIABLE NAMES FROM CONFIG
    -f | --fastq )          shift; f=$1 ;; # fastq file
    -o | --outdir )         shift; OUTDIR=$1 ;; # output dir
    -R | --region )         shift; SEQREG=$1 ;; # (optional) region of specific interest, e.g. targeted reseq
    --forceSingle )         FORCESINGLE=1;;
    --recover-from )        shift; NGSANE_RECOVERFROM=$1 ;; # attempt to recover from log file
    -h | --help )           usage ;;
    * )                     echo "dont understand $1"
    esac
    shift
done

. $CONFIG
. ${NGSANE_BASE}/conf/header.sh
. $CONFIG

################################################################################
NGSANE_CHECKPOINT_INIT "programs"

# save way to load modules that itself loads other modules
hash module 2>/dev/null && for MODULE in $MODULE_FUSION; do module load $MODULE; done && module list 

export PATH=$PATH_FUSION:$PATH
echo "PATH=$PATH"

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--tophat      --\n "$(tophat --version)
[ -z "$(which tophat)" ] && echo "[ERROR] no tophat detected" && exit 1
echo -e "--bowtie     --\n "$(bowtie --version)
[ -z "$(which bowtie)" ] && echo "[ERROR] no bowtie detected" && exit 1
echo -e "--samtools    --\n "$(samtools 2>&1 | head -n 3 | tail -n-2)
[ -z "$(which samtools)" ] && echo "[ERROR] no samtools detected" && exit 1
echo -e "--circos      --\n "$(circos --version)
[ -z "$(which circos)" ] && echo "[ERROR] circos not detected" && exit 1

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "parameters"

[ ! -f $f ] && echo "[ERROR] input file not found: $f" 1>&2 && exit 1
# get basename of f (samplename)
n=${f##*/}
SAMPLE=${n/%$READONE.$FASTQ/}

if [ -z "$FASTA" ]; then
    echo "[ERROR] regerence genome not specified (FASTA)" 
    exit 1
elif [[ ! -e ${FASTA%.*}.1.ebwt ]]; then
    echo "[ERROR] Bowtie1 index not detected. Exeute bowtieIndex.sh first"
    exit 1
fi

# create symlinks to annotation data
echo "[NOTE] create symlinks"

#remove old files
if [ -z "$NGSANE_RECOVERFROM" ]; then
    if [ -d $OUTDIR ]; then rm -r $OUTDIR; fi
fi

if [ -n "$READONE" ] && [ "$READONE" == "$READTWO" ]; then
    echo "[ERROR] read1 == read2 " 1>&2 && exit 1
elif [ "$f" != "${f/%$READONE.$FASTQ/$READTWO.$FASTQ}" ] && [ -e ${f/%$READONE.$FASTQ/$READTWO.$FASTQ} ] && [ "$FORCESINGLE" = 0 ]; then
    PAIRED="1"
    f2=${f/%$READONE.$FASTQ/$READTWO.$FASTQ}
    echo "[NOTE] Paired library detected"
else
    PAIRED="0"
    echo "[NOTE] Single-Strand (unpaired) library detected"
fi

#is ziped ?
CAT="cat"
if [[ ${f##*.} == "gz" ]]; 
    then CAT="zcat"; 
elif [[ ${f##*.} == "bz2" ]]; 
    then CAT="bzcat"; 
fi

# get encoding
if [ -z "$FASTQ_PHRED" ]; then 
    FASTQ_ENCODING=$($CAT $f |  awk 'NR % 4 ==0' | python $NGSANE_BASE/tools/GuessFastqEncoding.py |  tail -n 1)
    if [[ "$FASTQ_ENCODING" == *Phred33* ]]; then
        FASTQ_PHRED="" # use default
    elif [[ "$FASTQ_ENCODING" == *Illumina* ]]; then
        FASTQ_PHRED="--phred64-quals"
    elif [[ "$FASTQ_ENCODING" == *Solexa* ]]; then
        FASTQ_PHRED="--solexa1.3-quals"
    else
        echo "[NOTE] cannot detect/don't understand fastq format: $FASTQ_ENCODING - using default"
    fi
    echo "[NOTE] $FASTQ_ENCODING fastq format detected"
fi

# check library info is set
if [ -z "$RNA_SEQ_LIBRARY_TYPE" ]; then
    echo "[ERROR] RNAseq library type not set (RNA_SEQ_LIBRARY_TYPE): either fr-unstranded or fr-firststrand"
    exit 1;
else
    echo "[NOTE] RNAseq library type: $RNA_SEQ_LIBRARY_TYPE"
fi
if [[ -z "$EXPID" || -z "$LIBRARY" || -z "$PLATFORM" ]]; then
    echo "[ERROR] library info not set (EXPID, LIBRARY, and PLATFORM): free text needed"
    exit 1;
else
    echo "[NOTE] EXPID $EXPID; LIBRARY $LIBRARY; PLATFORM $PLATFORM"
fi

[ ! -d $OUTDIR ] && mkdir -p $OUTDIR 

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "recall files from tape"

if [ -n "$DMGET" ]; then
    dmget -a $(dirname $FASTA)/*
    dmget -a ${f/$READONE/"*"}
    dmget -a $OUTDIR/*
fi

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "run tophat fusion"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    cd $OUTDIR
    
    [ -z "$ANNO_DIR" ] && echo "[ERROR] ANNO_DIR not specified" && exit 1
    [ -L ${PWD}/refGene.txt ] && rm ${PWD}/refGene.txt
    [ -L ${PWD}/ensGene.txt ] && rm ${PWD}/ensGene.txt
    [ -L ${PWD}/mcl ] && rm ${PWD}/mcl
    [ -L ${PWD}/blast ] && rm ${PWD}/blast
    
    ln -s ${ANNO_DIR}/refGene.txt $PWD/refGene.txt
    ln -s ${ANNO_DIR}/ensGene.txt $PWD/ensGene.txt
    ln -s ${ANNO_DIR}/mcl $PWD/mcl
    ln -s ${ANNO_DIR}/blast $PWD/blast
    echo "[NOTE] Symlinks to annotations created."

    echo "[NOTE] $f"
    echo "[NOTE] $f2"
    TOPHAT_FUSION_OUT="${OUTDIR}/tophat_sample"
    
    [ ! -d $TOPHAT_FUSION_OUT ] && mkdir -p $TOPHAT_FUSION_OUT     
    echo "[NOTE] tophat out $TOPHAT_FUSION_OUT "
    
    RUN_COMMAND="tophat -o $TOPHAT_FUSION_OUT -p $CPU_FUSION --fusion-search --keep-fasta-order --bowtie1 --no-coverage-search $TOPHATFUSIONADDPARAM ${FASTA%.*} $f $f2"
    echo $RUN_COMMAND && eval $RUN_COMMAND

    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK ${TOPHAT_FUSION_OUT}/accepted_hits.bam

fi 
################################################################################
NGSANE_CHECKPOINT_INIT "run fusion-post"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    cd $OUTDIR

    [ -L ${PWD}/refGene.txt ] && rm ${PWD}/refGene.txt  
    [ -L ${PWD}/ensGene.txt ] && rm ${PWD}/ensGene.txt  
    [ -L ${PWD}/mcl ] && rm ${PWD}/mcl  
    [ -L ${PWD}/blast ] && rm ${PWD}/blast
  
    ln -s ${ANNO_DIR}/refGene.txt $PWD/refGene.txt
    ln -s ${ANNO_DIR}/ensGene.txt $PWD/ensGene.txt
    ln -s ${ANNO_DIR}/mcl $PWD/mcl
    ln -s ${ANNO_DIR}/blast $PWD/blast

    echo "[NOTE] Symlinks to annotations checked (remade)."

    RUN_COMMAND="tophat-fusion-post -o $OUTDIR -p $CPU_FUSION $TOPHATFUSIONADDPARAMPOST ${FASTA%.*}"
    echo $RUN_COMMAND && eval $RUN_COMMAND

    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/result.txt

fi 
################################################################################
NGSANE_CHECKPOINT_INIT "run circos"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    cd $OUTDIR
    cat $OUTDIR/result.txt | awk '{ print $3, $4, $4,$6,$7,$7 }' | sed 's/chr/hs/g' > $OUTDIR/my_link.txt    
    cat $OUTDIR/result.txt | awk '{ print $3,$4,$4,$2 }' | sed 's/chr/hs/g' > $OUTDIR/link_names.txt
    cat $OUTDIR/result.txt | awk '{ print $6,$7,$7,$5 }' | sed 's/chr/hs/g' >> $OUTDIR/link_names.txt
    
    [ -f $PWD/ticks.conf ] && rm $PWD/ticks.conf
    [ -f $PWD/ideogram.conf ] && rm $PWD/ideogram.conf
    [ -f $PWD/base.conf ] && rm $PWD/base.conf
      
    ln -s ${NGSANE_BASE}/tools/circos_tophat_fusion/ticks.conf  $PWD/ticks.conf
    ln -s ${NGSANE_BASE}/tools/circos_tophat_fusion/ideogram.conf  $PWD/ideogram.conf
    ln -s ${NGSANE_BASE}/tools/circos_tophat_fusion/base.conf  $PWD/base.conf
    
    circos -conf $PWD/base.conf

    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/circos.png

fi 
###############################################################################
NGSANE_CHECKPOINT_INIT "cleanup"

[ -e $PWD/ticks.conf ] && rm $PWD/ticks.conf
[ -e $PWD/ideogram.conf ] && rm $PWD/ideogram.conf
[ -e $PWD/base.conf ] && rm $PWD/base.conf
[ -L $PWD/refGene.txt ] && rm $PWD/refGene.txt
[ -L $PWD/ensGene.txt ] && rm $PWD/ensGene.txt
[ -L $PWD/mcl ] && rm $PWD/mcl
[ -L $PWD/blast ] && rm $PWD/blast
[ -e $PWD/ticks.conf ] && rm $PWD/ticks.conf
  
[ -z "$KEEP_FUSION_BAM" ] && [  -d ${OUTDIR}/tophat_sample ] && rm -r ${OUTDIR}/tophat_sample"

NGSANE_CHECKPOINT_CHECK
################################################################################
echo ">>>>> fusion search with TopHat - FINISHED"
echo ">>>>> enddate "`date`
