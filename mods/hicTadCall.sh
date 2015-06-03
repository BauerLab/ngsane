#!/bin/bash -e

# Script running domaincall on HiC data
# Expects fithic contact matrix files as input.
# author: Fabian Buske
# date: May 2015

# messages to look out for -- relevant for the QC.sh script:
# QCVARIABLES,Resource temporarily unavailable
# RESULTFILENAME <DIR>/<TASK>/<SAMPLE>.txt.gz

echo ">>>>> Chromatin organization with domaincall"
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"

function usage {
echo -e "usage: $(basename $0) -k NGSANE -f INPUT -o OUTDIR [OPTIONS]"
exit
}

if [ ! $# -gt 3 ]; then usage ; fi

#INPUTS
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the NGSANE repository
        -f | --file )           shift; f=$1 ;; # input file
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
hash module 2>/dev/null && for MODULE in $MODULE_HICTADCALL; do module load $MODULE; done && module list 

export PATH=$PATH_FITHIC:$PATH
echo "PATH=$PATH"
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--Python      --\n" $(python --version)
[ -z "$(which python)" ] && echo "[ERROR] no python detected" && exit 1
hash module 2>/dev/null && echo -e "--Python libs --\n "$(yolk -l)
echo -e "--Matlab (MCR)--\n "$(echo "$MCRROOT")
[ -z "$(echo $MCRROOT)" ] && echo "[ERROR] no matlab runtime environment detected" && exit 1
echo -e "--TADbit      --\n "$(yolk -l | fgrep -w TADbit | fgrep -v -w "non-active")
if [[ "$(yolk -l | fgrep -w TADbit | fgrep -v -w "non-active" | wc -l | awk '{print $1}')" == 0 ]]; then echo "[WARN] no TADbit detected"; TADBIT=""; elif [ -n "$CALL_TAD_CHROMOSOMES" ]; then TADBIT="--matrixFormat tadbit"; fi
echo -e "--bedToBigBed --\n "$(bedToBigBed 2>&1 | tee | head -n 1 )
[ -z "$(which bedToBigBed)" ] && echo "[WARN] bedToBigBed not detected, cannot compress tad bed file"
echo -e "--tabix       --\n "$(tabix 2>&1 | tee | grep "Version")
[ -z "$(which tabix)" ] && echo "[WARN] tabix not detected, cannot index bed file"

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "parameters"

# Default to bam
[ -z "$INPUT_HICTADCALL_SUFFIX" ] && $INPUT_HICTADCALL_SUFFIX ="$.contactCounts.gz"

# get basename of f
n=${f##*/}
SAMPLE=${n/%$INPUT_HICTADCALL_SUFFIX/}

# delete old bam files unless attempting to recover
if [ -z "$NGSANE_RECOVERFROM" ]; then
    [ -d $OUTDIR/$SAMPLE ] && rm -r $OUTDIR/$SAMPLE
    [ -f $OUTDIR/$SAMPLE.log ] && rm $OUTDIR/$SAMPLE.log
fi

if [ -z "$HIC_RESOLUTION" ]; then
    echo "[ERROR] HiC resolution not specified"
    exit 1
fi

if [[ -n "$FITHIC_CHROMOSOMES" ]]; then
    FITHIC_CHROMOSOMES="--chrompattern '$FITHIC_CHROMOSOMES'"
fi

if [ -n "$FITHIC_START_FROM_FRAGMENTPAIRS" ]; then
    FITHIC_START_FROM_FRAGMENTPAIRS="--inputIsFragmentPairs"
fi

THISTMP=$TMP"/"$(whoami)"/"$(echo $OUTDIR/$SAMPLE | md5sum | cut -d' ' -f1)
[ -d $THISTMP ] && rm -r $THISTMP
mkdir -p $THISTMP

mkdir -p $OUTDIR/$SAMPLE

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "recall files from tape"

if [ -n "$DMGET" ]; then
	dmget -a $f
	dmget -a $OUTDIR/*
fi

NGSANE_CHECKPOINT_CHECK

################################################################################
NGSANE_CHECKPOINT_INIT "count Interactions"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    mkdir -p $OUTDIR/$SAMPLE

    if [ -n "$FITHIC_START_FROM_FRAGMENTPAIRS" ]; then 
        cp ${FASTA%.*}.chrom.sizes $OUTDIR/$SAMPLE/chromsizes
        RUN_COMMAND="python ${NGSANE_BASE}/tools/fithic-fixedBins/fithicCreate2DcontactMap.py $FITHIC_START_FROM_FRAGMENTPAIRS --resolution=$HIC_RESOLUTION --chromsizes=$OUTDIR/$SAMPLE/chromsizes $FITHIC_CHROMOSOMES --outputDir=$OUTDIR/$SAMPLE --outputFilename $SAMPLE $f > $OUTDIR/$SAMPLE.log"

    else
        # extract chrom sizes from Bam
        samtools view -H $f | fgrep -w '@SQ' | sed 's/:/\t/g' | awk '{OFS="\t";print $3,$5}' > $OUTDIR/$SAMPLE/chromsizes

        # ensure name sorted bam required
        RUN_COMMAND="samtools sort -n -O bam -@ $CPU_FITHIC -o $THISTMP/$SAMPLE.bam -T $THISTMP/$SAMPLE.tmp $f"
        echo $RUN_COMMAND && eval $RUN_COMMAND

        RUN_COMMAND="python ${NGSANE_BASE}/tools/fithic-fixedBins/fithicCreate2DcontactMap.py --resolution=$HIC_RESOLUTION --chromsizes=$OUTDIR/$SAMPLE/chromsizes $FITHIC_CHROMOSOMES --outputDir=$OUTDIR/$SAMPLE --outputFilename $SAMPLE $THISTMP/$SAMPLE.bam > $OUTDIR/$SAMPLE.log"
    fi

    echo $RUN_COMMAND && eval $RUN_COMMAND
    echo $OUTDIR/$SAMPLE/done.txt

    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE/done.txt
    
    rm -f $OUTDIR/$SAMPLE/done.txt

fi

################################################################################
NGSANE_CHECKPOINT_INIT "DI matrix"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    for MATRIX in $OUTDIR/$SAMPLE/*.matrix; do 
        # run DI matrix script one chromosome at a time
        perl `which DI_from_matrix.pl` $MATRIX $HIC_RESOLUTION 20000000 $OUTDIR/$SAMPLE/chromsizes > ${i/%.matrix/.di.txt}
    done 
    
    # combine into one big DI matrix
    cat $OUTDIR/$SAMPLE/*.di.txt > $OUTDIR/$SAMPLE.di.txt
            
    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE/$SAMPLE.di.txt

fi
exit 1
################################################################################
NGSANE_CHECKPOINT_INIT "call topological domains"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    if [[ -n "$CALL_TAD_CHROMOSOMES" && -n "$TADBIT" ]]; then

        RUN_COMMAND="${CHANCE} compENCODE -b ${GENOME_ASSEMBLY} -t bam -o $OUTDIR/$SAMPLE.compENCODE -e $EXPERIMENTID --ipfile ${f} --ipsample $SAMPLE --inputfile ${CHIPINPUT} --inputsample $CONTROL"
        echo $RUN_COMMAND && eval $RUN_COMMAND
        
        # mark checkpoint
        NGSANE_CHECKPOINT_CHECK

    else
        echo "[NOTE] skipping topological domain calling (TADbit)"
        NGSANE_CHECKPOINT_CHECK
    fi
fi


################################################################################
NGSANE_CHECKPOINT_INIT "create tabix files"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    if hash tabix; then 
        [ -f $OUTDIR/$SAMPLE.bed.gz ] && rm $OUTDIR/$SAMPLE.bed.gz*

        zcat $OUTDIR/$SAMPLE.txt.gz | awk -v R=$(( $HIC_RESOLUTION / 2)) '{OFS="\t";print $1,$2-R+1,$2+R-2,$3":"$4-R+1"-"$4+R-2","$10,"1","."; print $3,$4-R+1,$4+R-1,$1":"$2-R+1"-"$2+R-1","$10,"2","."}' \
            | bedtools sort | bedtools intersect -a - -b <(awk '{OFS="\t";print $1,0,$2}' $OUTDIR/$SAMPLE/chromsizes ) > $OUTDIR/$SAMPLE.bed
        
        bgzip $OUTDIR/$SAMPLE.bed
        tabix -p bed $OUTDIR/$SAMPLE.bed.gz
        
        # mark checkpoint
        NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE.bed.gz
    else
        echo "[NOTE] skipping tabix file creation"
    fi
fi
################################################################################
NGSANE_CHECKPOINT_INIT "cleanup"

if [ -z "$FITHIC_KEEPCONTACTMATRIX" ]; then
    rm -f $OUTDIR/$SAMPLE/$SAMPLE*.matrix.gz
fi
rm -f -r $OUTDIR/$SAMPLE/*.border
rm -f $OUTDIR/$SAMPLE/$SAMPLE.ice.txt
rm -f $OUTDIR/$SAMPLE/chromsizes

NGSANE_CHECKPOINT_CHECK
################################################################################

[ -e $OUTDIR/$SAMPLE.txt.gz.dummy ] && rm $OUTDIR/$SAMPLE.txt.gz.dummy
echo ">>>>> Chromatin organization with domaincall - FINISHED"
echo ">>>>> enddate "`date`

