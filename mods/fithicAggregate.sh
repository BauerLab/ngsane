#!/bin/bash -e

# Script running fit-hi-c to call significant chromatin interactions form HiC
# experiments. Expects bam files as input.
# author: Fabian Buske
# date: Oct 2014

# messages to look out for -- relevant for the QC.sh script:
# QCVARIABLES,Resource temporarily unavailable

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
        -f | --file )           shift; FILES=$1 ;; # input files
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
echo -e "--samtools    --\n "$(samtools 2>&1 | head -n 3 | tail -n-2)
[ -z "$(which samtools)" ] && echo "[ERROR] no samtools detected" && exit 1
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

DATASETS="$(echo $FILES | tr ',' ' ')"
echo "[NOTE] Files: $DATASETS"

if [ -z "$FITHIC_POOLED_SAMPLE_NAME" ]; then
    echo "[ERROR] variable not set: FITHIC_POOLED_SAMPLE_NAME"
    exit 1
else
    SAMPLE=$FITHIC_POOLED_SAMPLE_NAME
    echo "Sample name: $SAMPLE"
fi

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

THISTMP=$TMP"/"$(whoami)"/"$(echo $OUTDIR/$SAMPLE | md5sum | cut -d' ' -f1)
[ -d $THISTMP ] && rm -r $THISTMP
mkdir -p $THISTMP

mkdir -p $OUTDIR/$SAMPLE/

# extract chrom sizes from Bam
if [ -s ${FASTA%.*}.chrom.sizes ]; then
    GENOME_CHROMSIZES=${FASTA%.*}.chrom.sizes
else
    samtools view -H ${DATASETS[0]} | fgrep -w '@SQ' | sed 's/:/\t/g' | awk '{OFS="\t";print $3,$5}' > $OUTDIR/$SAMPLE/chromsizes
fi

if [ ! -f $GENOME_CHROMSIZES ]; then
    echo "[ERROR] GENOME_CHROMSIZES not found. Excepted at $GENOME_CHROMSIZES"
    exit 1
else
    echo "[NOTE] Chromosome size: $GENOME_CHROMSIZES"
fi

if [ -n "$FITHIC_CISONLY" ]; then
    FITHIC_2DMATRIX="--create2DMatrixPerChr"
else
    FITHIC_2DMATRIX="--create2DMatrix"
fi

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

    if [ -n "$FITHIC_START_FROM_FRAGMENTPAIRS" ]; then

        RUN_COMMAND="python ${NGSANE_BASE}/tools/fithic-fixedBins/fithicCountInteractions.py $FITHIC_START_FROM_FRAGMENTPAIRS --CPU-processes $(($CPU_FITHIC<8?$CPU_FITHIC : 8)) --verbose $FITHIC_2DMATRIX --mappability=$MAPPABILITY --resolution=$HIC_RESOLUTION --chromsizes=$GENOME_CHROMSIZES $FITHIC_CHROMOSOMES --outputDir=$OUTDIR --outputFilename $SAMPLE $DATASETS > $OUTDIR/$SAMPLE.log"

  	else
        # ensure name sorted bam required
        SORTEDDATASET=""
        for DATA in $DATASETS; do
            D=${DATA##*/}
            D=${D/%$INPUT_FITHIC_SUFFIX/}
            RUNCOMMAND="samtools sort -n -O bam -@ $CPU_FITHIC -o $THISTMP/$D.bam -T $THISTMP/$D.tmp $DATA"
            echo $RUNCOMMAND && eval $RUNCOMMAND
            SORTEDDATASET="$SORTEDDATASET $THISTMP/$D.bam"
        done
        RUN_COMMAND="python ${NGSANE_BASE}/tools/fithic-fixedBins/fithicCountInteractions.py $FITHIC_2DMATRIX $TADBIT --mappability=$MAPPABILITY --resolution=$HIC_RESOLUTION --chromsizes=$GENOME_CHROMSIZES $FITHIC_CHROMOSOMES --outputDir=$OUTDIR --outputFilename $SAMPLE $SORTEDDATASET > $OUTDIR/$SAMPLE.log"

    fi

    echo $RUN_COMMAND && eval $RUN_COMMAND

    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE.fragmentLists.gz $OUTDIR/$SAMPLE.contactCounts.gz

fi

################################################################################
NGSANE_CHECKPOINT_INIT "contact matrices"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    if [ -n "$FITHIC_CISONLY" ]; then
      RUN_COMMAND="python ${NGSANE_BASE}/tools/fithic-fixedBins/fithicCreate2DcontactMap.py --verbose --fragmentFile=$OUTDIR/$SAMPLE.fragmentLists.gz --inputIsFragmentPairs --resolution=$HIC_RESOLUTION --chromsizes=$GENOME_CHROMSIZES $FITHIC_CHROMOSOMES --outputDir=$OUTDIR --outputFilename $SAMPLE --onlycis $OUTDIR/$SAMPLE.contactCounts.gz >> $OUTDIR/$SAMPLE.log"
    else
      RUN_COMMAND="python ${NGSANE_BASE}/tools/fithic-fixedBins/fithicCreate2DcontactMap.py --verbose --fragmentFile=$OUTDIR/$SAMPLE.fragmentLists.gz --inputIsFragmentPairs --resolution=$HIC_RESOLUTION --chromsizes=$GENOME_CHROMSIZES $FITHIC_CHROMOSOMES --outputDir=$OUTDIR --outputFilename $SAMPLE $OUTDIR/$SAMPLE.contactCounts.gz >> $OUTDIR/$SAMPLE.log"
    fi
    echo $RUN_COMMAND && eval $RUN_COMMAND
    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE*.matrix

fi

################################################################################
NGSANE_CHECKPOINT_INIT "HiCorrector"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    [ -f $OUTDIR/$SAMPLE.ice.txt.gz ] && rm $OUTDIR/$SAMPLE.ice.txt.gz

    if [ -n "$FITHIC_CISONLY" ]; then
        cat /dev/null > $OUTDIR/$SAMPLE.ice.txt.gz

        for CHR in $(cut -f1 $GENOME_CHROMSIZES); do
            c=$OUTDIR/$SAMPLE.$CHR.matrix

            if [ ! -f $c ]; then
                echo "[NOTE] Skipping $c (does not exist)"
            else

                if [[ "$CPU_FITHIC" -gt 1 ]]; then

                    RUN_COMMAND=$(which mpirun)" -np $CPU_FITHIC ic_mep --jobID=$SAMPLE --hasHeaderRow=0 --maxIteration=$HICORRECTOR_MAXITER --numRows="$(wc -l $c | awk '{print $1}')" --numTask=$CPU_FITHIC --memSizePerTask="$(echo "1 + $MEMORY_FITHIC * 900 / $CPU_FITHIC" | bc)" --inputFile=$c --outputFile=${c/%.matrix/.ice.txt} >> $OUTDIR/$SAMPLE.matrix_log"
                else
                    RUN_COMMAND="ic_mes $c $MEMORY_FITHIC "$(wc -l $c | awk '{print $1}')" $HICORRECTOR_MAXITER 0 0 ${c/%.matrix/.ice.txt} > $OUTDIR/$SAMPLE.matrix_log"
                fi
                echo $RUN_COMMAND && eval $RUN_COMMAND

                # combine ice files and convert to fit-hi-c expected bias format
                # add indices non-mappable bins
                paste ${c/%.matrix/.index} ${c/%.matrix/.ice.txt} > ${c/%.matrix/.ice.index.txt}
                zcat $OUTDIR/$SAMPLE.fragmentLists.gz | cut -f1,2 | egrep -w "^$CHR" | awk '{printf("%010d\t%s\n", NR, $0)}' > $OUTDIR/$SAMPLE.fragmentLists.index
                join -j 1 -o 2.2 -a 1 -e 0  $OUTDIR/$SAMPLE.fragmentLists.index ${c/%.matrix/.ice.index.txt} > ${c/%.matrix/.ice.txt}
                # combine ice files and convert to fit-hi-c expected bias format
                paste <(zcat $OUTDIR/$SAMPLE.fragmentLists.gz | cut -f1,2 | egrep -w "^$CHR") ${c/%.matrix/.ice.txt} | $GZIP >> $OUTDIR/$SAMPLE.ice.txt.gz
                # paste <(zcat $OUTDIR/$SAMPLE.fragmentLists.gz | cut -f1,2 | egrep -w "^$CHR") ${c/%.matrix/.ice.txt} | awk '{$3==0?$3=1:$3=$3; print $0}' | $GZIP >> $OUTDIR/$SAMPLE.ice.txt.gz
            fi
        done

    # otherwise create one big matrix
    else

        if [[ "$CPU_FITHIC" -gt 1 ]]; then

            RUN_COMMAND=$(which mpirun)" -np $CPU_FITHIC ic_mep --jobID=$SAMPLE --hasHeaderRow=0 --maxIteration=$HICORRECTOR_MAXITER --numRows="$(wc -l $OUTDIR/$SAMPLE.matrix | awk '{print $1}')" --numTask=$CPU_FITHIC --memSizePerTask="$(echo "1 + $MEMORY_FITHIC * 900 / $CPU_FITHIC" | bc)" --inputFile=$OUTDIR/$SAMPLE.matrix --outputFile=$OUTDIR/$SAMPLE.ice.txt > $OUTDIR/$SAMPLE.matrix_log"
        else
            RUN_COMMAND="ic_mes $OUTDIR/$SAMPLE.matrix $MEMORY_FITHIC "$(wc -l $OUTDIR/$SAMPLE.matrix | awk '{print $1}')" $HICORRECTOR_MAXITER 0 0 $OUTDIR/$SAMPLE.ice.txt > $OUTDIR/$SAMPLE.matrix_log"
        fi
        echo $RUN_COMMAND && eval $RUN_COMMAND

        # add indices non-mappable bins
        paste $OUTDIR/$SAMPLE.index $OUTDIR/$SAMPLE.ice.txt > $OUTDIR/$SAMPLE.ice.index.txt
        zcat $OUTDIR/$SAMPLE.fragmentLists.gz | cut -f1,2 | awk '{printf("%010d\t%s\n", NR, $0)}' > $OUTDIR/$SAMPLE.fragmentLists.index
        join -j 1 -o 2.2 -a 1 -e 0  $OUTDIR/$SAMPLE.fragmentLists.index $OUTDIR/$SAMPLE.ice.index.txt > $OUTDIR/$SAMPLE.ice.txt
        # convert to fit-hi-c expected bias format
        paste <(zcat $OUTDIR/$SAMPLE.fragmentLists.gz | cut -f1,2) $OUTDIR/$SAMPLE.ice.txt | $GZIP > $OUTDIR/$SAMPLE.ice.txt.gz
        # paste <(zcat $OUTDIR/$SAMPLE.fragmentLists.gz | cut -f1,2) $OUTDIR/$SAMPLE.ice.txt | awk '{$3==0?$3=1:$3=$3; print $0}' | $GZIP > $OUTDIR/$SAMPLE.ice.txt.gz

    fi

    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE.ice.txt.gz

    # cleanup
    rm -f $OUTDIR/$SAMPLE.ice.txt
    rm -f $OUTDIR/$SAMPLE*.matrix
    rm -f $OUTDIR/$SAMPLE.fragmentLists.index
    rm -f $OUTDIR/$SAMPLE*.index
fi

################################################################################
NGSANE_CHECKPOINT_INIT "fit-hi-c"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    cd $OUTDIR/$RESOLUTION
    RUN_COMMAND="python ${NGSANE_BASE}/tools/fithic-fixedBins/fit-hi-c-fixedSize-withBiases.py --lib=${SAMPLE} --biases=$OUTDIR/$SAMPLE.ice.txt.gz --fragments=$OUTDIR/$SAMPLE.fragmentLists.gz --interactions=$OUTDIR/$SAMPLE.contactCounts.gz --resolution $HIC_RESOLUTION >> $OUTDIR/$SAMPLE.log"
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
            | bedtools sort | bedtools intersect -a - -b <(awk '{OFS="\t";print $1,0,$2}' $GENOME_CHROMSIZES ) > $OUTDIR/$SAMPLE.bed

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

rm -f $OUTDIR/$SAMPLE.spline_pass1.res$HIC_RESOLUTION.significances.txt.gz

NGSANE_CHECKPOINT_CHECK
################################################################################
[ -e $OUTDIR/$FITHIC_POOLED_SAMPLE_NAME.txt.gz.dummy ] && rm $OUTDIR/$FITHIC_POOLED_SAMPLE_NAME.txt.gz.dummy
echo ">>>>> Chromatin organization with fit-hi-c - FINISHED"
echo ">>>>> enddate "`date`
