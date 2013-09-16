#!/bin/bash -e

# Script for ChIP-seq peak calling using peakranger.
# It takes read alignments in .bam format.
# It produces output files: peak regions in bed format
# author: Fabian Buske
# date: August 2013

# QCVARIABLES,Resource temporarily unavailable
# RESULTFILENAME <DIR>/<TASK>/<SAMPLE>-${CHIPINPUT##*/} | sed "s/.$ASD.bam/_region.bed/"

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

for MODULE in $MODULE_PEAKRANGER; do module load $MODULE; done  # save way to load modules that itself load other modules

export PATH=$PATH_PEAKRANGER:$PATH
module list
echo "PATH=$PATH"
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)

echo -e "--NGSANE     --\n" $(trigger.sh -v 2>&1)
echo -e "--R          --\n "$(R --version | head -n 3)
[ -z "$(which R)" ] && echo "[ERROR] no R detected" && exit 1
echo -e "--peakranger --\n "$(peakranger | head -n 3 | tail -n 1)
[ -z "$(which peakranger)" ] && echo "[ERROR] peakranger not detected" && exit 1

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="parameters"

if [ -z "$CHIPINPUT" ] || [ ! -f $CHIPINPUT ]; then
    echo "[ERROR] input control not provided or invalid (CHIPINPUT=\"$CHIPINPUT\")"
    exit 1
fi

# get basename of f
f=${f/%.dummy/} #if input came from pip
n=${f##*/}
c=${CHIPINPUT##*/}

if [ -z "$RECOVERFROM" ]; then
    [ -e $OUTDIR/${n/.$ASD.bam/}-${c/.$ASD.bam/}_region.bed ] && rm $OUTDIR/${n/.$ASD.bam/}-${c/.$ASD.bam/}*
fi

if [ "$PEAKRANGER_PEAKS" != "broad" ] && [ "$PEAKRANGER_PEAKS" != "sharp" ]; then
    echo "[ERROR] PEAKRANGER_PEAKS parameter not valid: $PEAKRANGER_PEAKS"
    exit 1
fi

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="recall files from tape"

if [ -n "$DMGET" ]; then
	dmget -a ${f}
	dmls -l ${f}
fi

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="peakranger"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else
    
    echo "[NOTE] data quality"
    peakranger nr --format bam --data $f --control $CHIPINPUT > $OUTDIR/${n/.$ASD.bam/}-${c/.$ASD.bam/}.summary.txt
    
    echo "[NOTE] library complexity"
    peakranger lc --data $f >> $OUTDIR/${n/.$ASD.bam/}-${c/.$ASD.bam/}.summary.txt
    
    echo "[NOTE] make wigpe"
    peakranger wigpe --data $f --output $OUTDIR/${n/.$ASD.bam/}-${c/.$ASD.bam/}

    if [ "$PEAKRANGER_PEAKS" == "broad" ]; then
        echo "[NOTE] calling broad peaks"
        RUN_COMMAND="peakranger ccat $PEAKRANGERADDPARAM --format bam --data  $f --control $CHIPINPUT --output $OUTDIR/${n/.$ASD.bam/}-${c/.$ASD.bam/} -t $CPU_PEAKRANGER"

    elif [ "$PEAKRANGER_PEAKS" == "sharp" ]; then
        echo "[NOTE] calling tight peaks"
        RUN_COMMAND="peakranger ranger $PEAKRANGERADDPARAM --format bam --data $f --control $CHIPINPUT --output $OUTDIR/${n/.$ASD.bam/}-${c/.$ASD.bam/} -t $CPU_PEAKRANGER"
    fi
    echo $RUN_COMMAND && eval $RUN_COMMAND
    
    # echo remove comments from bed files
    grep -v "#" $OUTDIR/${n/.$ASD.bam/}-${c/.$ASD.bam/}_region.bed > $OUTDIR/${n/.$ASD.bam/}-${c/.$ASD.bam/}_region.tmp && mv $OUTDIR/${n/.$ASD.bam/}-${c/.$ASD.bam/}_region.tmp $OUTDIR/${n/.$ASD.bam/}-${c/.$ASD.bam/}_region.bed
    
    grep -v "#" $OUTDIR/${n/.$ASD.bam/}-${c/.$ASD.bam/}_summit.bed > $OUTDIR/${n/.$ASD.bam/}-${c/.$ASD.bam/}_summit.tmp && mv $OUTDIR/${n/.$ASD.bam/}-${c/.$ASD.bam/}_summit.tmp $OUTDIR/${n/.$ASD.bam/}-${c/.$ASD.bam/}_summit.bed
    
    echo "Peaks: `wc -l $OUTDIR/${n/.$ASD.bam/}-${c/.$ASD.bam/}_region.bed | awk '{print $1}'`" >> $OUTDIR/${n/.$ASD.bam/}-${c/.$ASD.bam/}.summary.txt
    echo "Summits: `wc -l $OUTDIR/${n/.$ASD.bam/}-${c/.$ASD.bam/}_summit.bed | awk '{print $1}'`" >> $OUTDIR/${n/.$ASD.bam/}-${c/.$ASD.bam/}.summary.txt


    # mark checkpoint
    if [ -f $OUTDIR/${n/.$ASD.bam/}-${c/.$ASD.bam/}_region.bed ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi

################################################################################
CHECKPOINT="zip"

$GZIP -f $OUTDIR/${n/.$ASD.bam/}-${c/.$ASD.bam/}_details

echo -e "\n********* $CHECKPOINT\n"
################################################################################
[ -e $OUTDIR/${n/.$ASD.bam/}-${c/.$ASD.bam/}_region.bed.dummy ] && rm $OUTDIR/${n/.$ASD.bam/}-${c/.$ASD.bam/}_region.bed.dummy
echo ">>>>> ChIPseq analysis with Peakranger - FINISHED"
echo ">>>>> enddate "`date`

