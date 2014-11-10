#!/bin/bash -e

# Batch consistency analysis using IDR
# author: Fabian Buske
# date:  August 2014

# QCVARIABLES,Resource temporarily unavailable

echo ">>>>> IDR analysis"
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
        -f | --file )           shift; FILES=$1 ;; # files
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

# save way to load modules that itself load other modules
hash module 2>/dev/null && for MODULE in $MODULE_IDR; do module load $MODULE; done && module list

export PATH=$PATH_IDR:$PATH
echo "PATH=$PATH"

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--R           --\n "$(R --version | head -n 3)
[ -z "$(which R)" ] && echo "[ERROR] no R detected" && exit 1
echo -e "--bedToBigBed --\n "$(bedToBigBed 2>&1 | tee | head -n 1 )
[ -z "$(which bedToBigBed)" ] && echo "[WARN] bedToBigBed not detected, cannot compress bedgraphs"

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "parameters"

FILES=${FILES//,/ }

if [[ -n "$IDR_ISBROADPEAK" && "$IDR_ISBROADPEAK" != "F" && "$IDR_ISBROADPEAK" != "T" ]]; then
    echo "[ERROR] IDR_ISBROADPEAK has invalid value ($IDR_ISBROADPEAK != {T|F})"
    exit 1
elif [[ -z "$IDR_ISBROADPEAK" ]]; then
    IDR_ISBROADPEAK="F"
fi

[[ -z "$IDR_PEAKHALFWIDTH" ]] && IDR_PEAKHALFWIDTH="-1"

if [ -z "$FASTA" ]; then
    echo "[ERROR] no reference provided (FASTA)"
    exit 1
fi
GENOMESIZE=${FASTA%.*}.chrom.sizes

# unique temp folder that should be used to store temporary files
THISTMP=$TMP"/"$(whoami)"/"$(echo $OUTDIR | md5sum | cut -d' ' -f1)
mkdir -p $THISTMP

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "recall files from tape"

if [ -n "$DMGET" ]; then
	dmget -a $INPUTFILE
    dmget -a $OUTDIR/*
fi
    
NGSANE_CHECKPOINT_CHECK

################################################################################
NGSANE_CHECKPOINT_INIT "pair replicates"

COMMAND=$THISTMP/idr_replicates`date +%Y%m%d`.commands

cat /dev/null > $COMMAND

egrep '^IDR_REPLICATE ' $CONFIG | cut -d' ' -f 2- > $COMMAND.tmp
for d in ${DIR[@]}; do
    OUTDIR="$OUT/$d/$TASK_IDR"
    mkdir -p $OUTDIR
    while read -r -a REPLICATE; do
        R1=
        R2=
        for f in ${FILES[@]}; do 
            if [[ "$f" =~ "${REPLICATE[1]}" ]]; then R1=$f; fi; 
            if [[ "$f" =~ "${REPLICATE[2]}" ]]; then R2=$f; fi; 
        done
        if [[ -z "$R1" || -z "$R2" ]]; then
            echo "[ERROR] one of the replicate files could not be detected"
            exit 1
        fi

        echo -ne "cd $NGSANE_BASE/tools/idr/; Rscript batch-consistency-analysis.r $R1 $R2 $IDR_PEAKHALFWIDTH $OUTDIR/${REPLICATE[0]} $IDR_MINOVERLAPRATIO $IDR_ISBROADPEAK $IDR_RANKINGMEASURE $GENOMESIZE; Rscript batch-consistency-plot.r 1 $OUTDIR/${REPLICATE[0]} $OUTDIR/${REPLICATE[0]}; awk 'BEGIN{OFS=\"\t\"}{if(NR>1){if(\$3<\$7){start=\$3}else{start=\$7};if(\$4<\$8){end=\$8}else{end=\$4}; print \$2,start,end,\"IDR_\"(NR-1), \$11, \".\"}}' $OUTDIR/${REPLICATE[0]}-overlapped-peaks.txt | sed 's/\"//g' > $OUTDIR/${REPLICATE[0]}-overlapped-peaks.bed; if hash bedToBigBed; then bedToBigBed -type=bed6+3 $OUTDIR/${REPLICATE[0]}-overlapped-peaks.bed $GENOMESIZE $OUTDIR/${REPLICATE[0]}.bb; fi" >> $COMMAND
        echo ";" >> $COMMAND
    
    done < $COMMAND.tmp
done
rm $COMMAND.tmp

echo "COMMAND:" $(cat $COMMAND)
echo $COMMAND

# mark checkpoint checking for existance/non-empty result file
NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "calculate idr"

if hash parallel ; then

    echo "[NOTE] parallel processing"
    cat $COMMAND | parallel --verbose --joblog $THISTMP/$TASK_IDR.log --gnu --eta -j $CPU_IDR "eval {}" #> /dev/null 2>&1

else
    # serial processing
    echo "[NOTE] serial processing"
    while read line; do
        echo $line
        eval $line
    done < $COMMAND
fi

rm $COMMAND
    
NGSANE_CHECKPOINT_CHECK
################################################################################
echo ">>>>> IDR analysis - FINISHED"
echo ">>>>> enddate "`date`
