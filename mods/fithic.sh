#!/bin/bash -e

# Script running fit-hi-c to call significant chromatin interactions form HiC 
# experiments. Expects bam files as input.
# author: Fabian Buske
# date: Oct 2014

# messages to look out for -- relevant for the QC.sh script:
# QCVARIABLES,Resource temporarily unavailable
# RESULTFILENAME <DIR>/<TASK>/<SAMPLE>.txt.gz

echo ">>>>> Chromatin organization with fit-hi-c "
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"

function usage {
echo -e "usage: $(basename $0) -k NGSANE -f FASTQ -o OUTDIR [OPTIONS]"
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
hash module 2>/dev/null && for MODULE in $MODULE_FITHIC; do module load $MODULE; done && module list 

export PATH=$PATH_FITHIC:$PATH
echo "PATH=$PATH"
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--Python      --\n" $(python --version)
[ -z "$(which python)" ] && echo "[ERROR] no python detected" && exit 1
hash module 2>/dev/null && echo -e "--Python libs --\n "$(yolk -l)
echo -e "--HiCorrector --\n "$(ic_mep 2>&1 | tee | grep Version)
[ -z "$(which ic_mep)" ] && echo "[ERROR] no HiCorrection detected" && exit 1
echo -e "--fit-hi-c    --\n "$(python $NGSANE_BASE/tools/fithic-fixedBins/fit-hi-c-fixedSize-withBiases.py --version | head -n 1)
echo -e "--TADbit      --\n "$(yolk -l | fgrep -w TADbit | fgrep -v -w "non-active")
if [[ "$(yolk -l | fgrep -w TADbit | fgrep -v -w "non-active" | wc -l | awk '{print $1}')" == 0 ]]; then echo "[WARN] no TADbit detected"; TADBIT=""; elif [ -n "$CALL_TAD_CHROMOSOMES" ]; then TADBIT="--create2DMatrixPerChr"; fi
echo -e "--bedToBigBed --\n "$(bedToBigBed 2>&1 | tee | head -n 1 )
[ -z "$(which bedToBigBed)" ] && echo "[WARN] bedToBigBed not detected, cannot compress tad bed file"
echo -e "--tabix       --\n "$(tabix 2>&1 | tee | grep "Version")
[ -z "$(which tabix)" ] && echo "[WARN] tabix not detected, cannot index bed file"

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "parameters"

# Default to bam
[ -z "$INPUT_FITHIC_SUFFIX" ] && $INPUT_FITHIC_SUFFIX="$ASD.bam"

# get basename of f
n=${f##*/}
SAMPLE=${n/%$INPUT_FITHIC_SUFFIX/}

# delete old bam files unless attempting to recover
if [ -z "$NGSANE_RECOVERFROM" ]; then
    [ -d $OUTDIR/$SAMPLE ] && rm -r $OUTDIR/$SAMPLE
    [ -f $OUTDIR/$SAMPLE.log ] && rm $OUTDIR/$SAMPLE.log
fi

if [ -z "$HICORRECTOR_MAXITER" ];then
    echo "[NOTE] HICORRECTOR_MAXITER set to 100 by default"
    HICORRECTOR_MAXITER=100
fi

if [ -z "$MAPPABILITY" ]; then
    echo "[ERROR] Mappability not specified"
    exit 1
fi

if [ -z "$HIC_RESOLUTION" ]; then
    echo "[ERROR] HiC resolution not specified"
    exit 1
fi

if [[ -z "$FITHIC_QVALUETHRESHOLD" ]]; then
    FITHIC_QVALUETHRESHOLD=0.01
fi
echo "[NOTE] Q-value threshold: $FITHIC_QVALUETHRESHOLD"

if [[ -z "$FITHIC_MAPPABILITYTHRESHOLD" ]];then
    echo "[ERROR] FITHIC_MAPPABILITYTHRESHOLD not set"
    exit 1
fi

if [[ -n "$FITHIC_CHROMOSOMES" ]]; then
    FITHIC_CHROMOSOMES="--chrompattern '$FITHIC_CHROMOSOMES'"
fi

if [ -n "$FITHIC_START_FROM_FRAGMENTPAIRS" ]; then
    FITHIC_START_FROM_FRAGMENTPAIRS="--inputIsFragmentPairs"
fi

if [ -n "$FITHIC_CISONLY" ]; then
    FITHIC_2DMATRIX="--create2DMatrixPerChr"
else
    FITHIC_2DMATRIX="--create2DMatrix"
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
        RUN_COMMAND="python ${NGSANE_BASE}/tools/fithic-fixedBins/fithicCountInteractions.py $FITHIC_START_FROM_FRAGMENTPAIRS $FITHIC_2DMATRIX $TADBIT --mappability=$MAPPABILITY --resolution=$HIC_RESOLUTION --chromsizes=$OUTDIR/$SAMPLE/chromsizes $FITHIC_CHROMOSOMES --outputDir=$OUTDIR/$SAMPLE --outputFilename $SAMPLE $f > $OUTDIR/$SAMPLE.log"

    else
        # extract chrom sizes from Bam
        samtools view -H $f | fgrep -w '@SQ' | sed 's/:/\t/g' | awk '{OFS="\t";print $3,$5}' > $OUTDIR/$SAMPLE/chromsizes

        # ensure name sorted bam required
        RUN_COMMAND="samtools sort -n -O bam -@ $CPU_FITHIC -o $THISTMP/$SAMPLE.bam -T $THISTMP/$SAMPLE.tmp $f"
        echo $RUN_COMMAND && eval $RUN_COMMAND

        RUN_COMMAND="python ${NGSANE_BASE}/tools/fithic-fixedBins/fithicCountInteractions.py $FITHIC_2DMATRIX $TADBIT --mappability=$MAPPABILITY --resolution=$HIC_RESOLUTION --chromsizes=$OUTDIR/$SAMPLE/chromsizes $FITHIC_CHROMOSOMES --outputDir=$OUTDIR/$SAMPLE --outputFilename $SAMPLE $THISTMP/$SAMPLE.bam > $OUTDIR/$SAMPLE.log"
    fi

    echo $RUN_COMMAND && eval $RUN_COMMAND

    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE/$SAMPLE.fragmentLists.gz $OUTDIR/$SAMPLE/$SAMPLE.contactCounts.gz

fi

################################################################################
NGSANE_CHECKPOINT_INIT "HiCorrector"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then


    if [ -n "$FITHIC_CISONLY" ]; then
        cat /dev/null > $OUTDIR/$SAMPLE/$SAMPLE.ice.txt.gz  

        for c in $OUTDIR/$SAMPLE/$SAMPLE.*.matrix; do
            CHR=$(eval "echo '${c##*/}' | sed -e 's|$SAMPLE\.\(.*\)\.matrix\$|\1|'")

            if [[ "$CPU_FITHIC" -gt 1 ]]; then
            
                RUN_COMMAND=$(which mpirun)" -np $CPU_FITHIC ic_mep --jobID=$SAMPLE --hasHeaderRow=0 --maxIteration=$HICORRECTOR_MAXITER --numRows="$(wc -l $c | awk '{print $1}')" --numTask=$CPU_FITHIC --memSizePerTask="$(echo "1 + $MEMORY_FITHIC * 900 / $CPU_FITHIC" | bc)" --inputFile=$c --outputFile=${c/%.matrix/.ice.txt} >> $OUTDIR/$SAMPLE/$SAMPLE.matrix_log"
            else
                RUN_COMMAND="ic_mes $c $MEMORY_FITHIC "$(wc -l $c | awk '{print $1}')" $HICORRECTOR_MAXITER 0 0 ${c/%.matrix/.ice.txt} > $OUTDIR/$SAMPLE/$SAMPLE.matrix_log"     
            fi
            echo $RUN_COMMAND && eval $RUN_COMMAND

            # combine ice files and convert to fit-hi-c expected bias format
            paste <(zcat $OUTDIR/$SAMPLE/$SAMPLE.fragmentLists.gz | cut -f1,2 | egrep -w "^$CHR") ${c/%.matrix/.ice.txt} | awk '{$3==0?$3=1:$3=$3; print $0}' | $GZIP >> $OUTDIR/$SAMPLE/$SAMPLE.ice.txt.gz  

        done

    # otherwise create one big matrix
    else
    
        if [[ "$CPU_FITHIC" -gt 1 ]]; then
            RUN_COMMAND=$(which mpirun)" -np $CPU_FITHIC ic_mep --jobID=$SAMPLE --hasHeaderRow=0 --maxIteration=$HICORRECTOR_MAXITER --numRows="$(wc -l $OUTDIR/$SAMPLE/$SAMPLE.matrix | awk '{print $1}')" --numTask=$CPU_FITHIC --memSizePerTask="$(echo "1 + $MEMORY_FITHIC * 900 / $CPU_FITHIC" | bc)" --inputFile=$OUTDIR/$SAMPLE/$SAMPLE.matrix --outputFile=$OUTDIR/$SAMPLE/$SAMPLE.ice.txt > $OUTDIR/$SAMPLE/$SAMPLE.matrix_log"
        else
            RUN_COMMAND="ic_mes $OUTDIR/$SAMPLE/$SAMPLE.matrix $MEMORY_FITHIC "$(wc -l $OUTDIR/$SAMPLE/$SAMPLE.matrix | awk '{print $1}')" $HICORRECTOR_MAXITER 0 0 $OUTDIR/$SAMPLE/$SAMPLE.ice.txt > $OUTDIR/$SAMPLE/$SAMPLE.matrix_log"
        fi
        echo $RUN_COMMAND && eval $RUN_COMMAND

        # convert to fit-hi-c expected bias format
        paste <(zcat $OUTDIR/$SAMPLE/$SAMPLE.fragmentLists.gz | cut -f1,2) $OUTDIR/$SAMPLE/$SAMPLE.ice.txt | awk '{$3==0?$3=1:$3=$3; print $0}' | $GZIP > $OUTDIR/$SAMPLE/$SAMPLE.ice.txt.gz  

    fi
    
    $GZIP $OUTDIR/$SAMPLE/$SAMPLE*.matrix
        
    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE/$SAMPLE.ice.txt.gz

fi

################################################################################
NGSANE_CHECKPOINT_INIT "call topological domains with TADbit"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    if [[ -n "$CALL_TAD_CHROMOSOMES" && -n "$TADBIT" ]]; then

        mkdir -p $OUTDIR/$SAMPLE/$SAMPLE
        MATRIXFILES=$(eval ls $OUTDIR/$SAMPLE/$SAMPLE.$CALL_TAD_CHROMOSOMES.matrix.gz)
        RUN_COMMAND="python ${NGSANE_BASE}/tools/fithic-fixedBins/callTADs.py --outputDir=$OUTDIR/$SAMPLE --outputFilename=$SAMPLE --threads=$CPU_FITHIC --resolution=$HIC_RESOLUTION $MATRIXFILES >> $OUTDIR/$SAMPLE.log"
        echo $RUN_COMMAND && eval $RUN_COMMAND
       
        cat /dev/null > $OUTDIR/$SAMPLE.tad.bed
        # convert .border file to beds
        for i in $OUTDIR/$SAMPLE/*.border; do
            CHROM=$(echo $i | tr '.' '\n' | tail -n2 | head -n 1)
            awk -v c=$CHROM -v r=$HIC_RESOLUTION 'BEGIN{OFS="\t"}{if (NR>1){print c,$2*r,$3*r,c"_"$1,int($4),"."}}' $i | bedtools intersect -a - -b <( awk '{OFS="\t"; print $1,0,$2}' $OUTDIR/$SAMPLE/chromsizes)  >> $OUTDIR/$SAMPLE.tad.bed
        done

        # create bigbed making sure score is between 0 and 1000
        if hash bedToBigBed; then 
            echo "[NOTE] create bigbed from tads" 
            bedtools sort -i $OUTDIR/$SAMPLE.tad.bed > $OUTDIR/$SAMPLE.tad.tmp
            bedToBigBed -type=bed6 $OUTDIR/$SAMPLE.tad.tmp $OUTDIR/$SAMPLE/chromsizes $OUTDIR/$SAMPLE.tad.bb
            rm $OUTDIR/$SAMPLE.tad.tmp
        fi

        # mark checkpoint
        NGSANE_CHECKPOINT_CHECK

    else
        echo "[NOTE] skipping topological domain calling (TADbit)"
        NGSANE_CHECKPOINT_CHECK
    fi
fi

################################################################################
NGSANE_CHECKPOINT_INIT "fit-hi-c"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    cd $OUTDIR/$RESOLUTION
    RUN_COMMAND="python ${NGSANE_BASE}/tools/fithic-fixedBins/fit-hi-c-fixedSize-withBiases.py $FITHICADDPARAM --lib=${SAMPLE} --biases=$OUTDIR/$SAMPLE/$SAMPLE.ice.txt.gz --fragments=$OUTDIR/$SAMPLE/$SAMPLE.fragmentLists.gz --interactions=$OUTDIR/$SAMPLE/$SAMPLE.contactCounts.gz --resolution $HIC_RESOLUTION >> $OUTDIR/$SAMPLE.log"
    echo $RUN_COMMAND && eval $RUN_COMMAND
    
    zcat $OUTDIR/$SAMPLE.spline_pass1.res$HIC_RESOLUTION.significances.txt.gz | awk -v q=$FITHIC_QVALUETHRESHOLD '$7<=q' | sort -k7g | gzip > $OUTDIR/$SAMPLE.txt.gz

    SIGCISINTERACTIONS=$(zcat $OUTDIR/$SAMPLE.txt.gz |  awk '$1==$3' | wc -l | cut -d' ' -f 2)
    SIGTRANSINTERACTIONS=$(zcat $OUTDIR/$SAMPLE.txt.gz |  awk '$1!=$3' | wc -l | cut -d' ' -f 2)
    echo "Significant cis interactions: $SIGCISINTERACTIONS" >> $OUTDIR/$SAMPLE.log
    echo "Significant trans interactions: $SIGTRANSINTERACTIONS" >> $OUTDIR/$SAMPLE.log
    
    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE.txt.gz
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
echo ">>>>> Chromatin organization with fit-hi-c - FINISHED"
echo ">>>>> enddate "`date`

