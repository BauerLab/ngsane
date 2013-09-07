#!/bin/bash -e

echo ">>>>> ChIPseq analysis with Homer"
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"

function usage {
echo -e "usage: $(basename $0) -k NGSANE -f FASTQ -r REFERENCE -o OUTDIR [OPTIONS]"
exit
}

# Script for ChIP-seq peak calling using Homer.
# It takes read alignments in .bam format.
# It produces output files: peak regions in bed format
# author: Fabian Buske
# date: August 2013

# QCVARIABLES,Resource temporarily unavailable
if [ ! $# -gt 3 ]; then usage ; fi

#INPUTS                                                                                                           
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the NGSANE repository
        -f | --bam )            shift; f=$1 ;; # bam file
        -o | --outdir )         shift; MYOUT=$1 ;; # output dir 
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

for MODULE in $MODULE_HOMERCHIPSEQ; do module load $MODULE; done  # save way to load modules that itself load other modules

export PATH=$PATH_HOMERCHIPSEQ:$PATH
module list
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

echo -e "\n********* $CHECKPOINT"
################################################################################
CHECKPOINT="parameters"

# get basename of f
n=${f##*/}
c=${CHIPINPUT##*/}

echo -e "\n********* $CHECKPOINT"
################################################################################
CHECKPOINT="recall files from tape"

if [ -n "$DMGET" ]; then
    dmget -a $f
    dmls -l $f
fi

echo -e "\n********* $CHECKPOINT"
################################################################################
CHECKPOINT="create tagdirectory"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep "********* $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else
    
    TAGDIRECTORY=$MYOUT/${n/%.$ASD.bam/_homer}
    mkdir -p $TAGDIRECTORY
    RUN_COMMAND="makeTagDirectory $TAGDIRECTORY $f $HOMER_CHIPSEQ_TAGDIR_ADDPARAM"
    echo $RUN_COMMAND && eval $RUN_COMMAND
    
    
    if [ -n "$CHIPINPUT" ];then
        TAGDIRECTORY=$MYOUT/${n/%.$ASD.bam/_homer}
        mkdir -p ${TAGDIRECTORY}_input
        # copy input to prevent interfering concurrent processing by homer
        cp $CHIPINPUT ${TAGDIRECTORY}
        RUN_COMMAND="makeTagDirectory ${TAGDIRECTORY}_input ${TAGDIRECTORY}/$c $HOMER_CHIPSEQ_TAGDIR_ADDPARAM"
        echo $RUN_COMMAND && eval $RUN_COMMAND
        rm $TAGDIRECTORY/$c
        INPUT=${c/.$ASD.bam/}
    else
        INPUT="NONE"
    fi

    # mark checkpoint
    [ -d $TAGDIRECTORY ] && echo -e "\n********* $CHECKPOINT" && unset RECOVERFROM
fi

################################################################################
CHECKPOINT="find peaks"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep "********* $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

    RUN_COMMAND="findPeaks $TAGDIRECTORY -style $HOMER_CHIPSEQ_STYLE  $HOMER_CHIPSEQ_FINDPEAKS_ADDPARAM -o auto"
    if [ -n "$CHIPINPUT" ];then
      RUN_COMMAND="$RUN_COMMAND -i ${TAGDIRECTORY}_input"  
    fi
    echo $RUN_COMMAND && eval $RUN_COMMAND
    
    if [ "$HOMER_CHIPSEQ_STYLE" == "factor" ]; then
        pos2bed.pl $MYOUT/${n/.$ASD.bam/_homer}/peaks.txt > $MYOUT/${n/.$ASD.bam/}-${INPUT}_peaks.bed
        grep "^#" $MYOUT/${n/.$ASD.bam/_homer}/peaks.txt > $MYOUT/${n/.$ASD.bam/}-${INPUT}.summary.txt
    
    elif [ "$HOMER_CHIPSEQ_STYLE" == "histone" ]; then
        pos2bed.pl $MYOUT/${n/.$ASD.bam/_homer}/regions.txt > $MYOUT/${n/.$ASD.bam/}-${INPUT}_regions.bed
        grep "^#" $MYOUT/${n/.$ASD.bam/_homer}/regions.txt > $MYOUT/${n/.$ASD.bam/}-${INPUT}.summary.txt
    fi

    # mark checkpoint
    [ -f $MYOUT/${n/.$ASD.bam/}-${INPUT}.summary.txt ] && echo -e "\n********* $CHECKPOINT" && unset RECOVERFROM
fi

################################################################################
CHECKPOINT="cleanup"

if [ -z "$HOMER_KEEPTAGDIRECTORY" ]; then
   rm $TAGDIRECTORY/*.tags.tsv
fi

if [ -n "$CHIPINPUT" ] && [ -d ${TAGDIRECTORY}_input ]; then
    rm -r ${TAGDIRECTORY}_input
fi

echo -e "\n********* $CHECKPOINT"
################################################################################
echo ">>>>> ChIPseq analysis with Homer - FINISHED"
echo ">>>>> enddate "`date`

