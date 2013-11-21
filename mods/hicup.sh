#!/bin/bash -e

# Script running hicup including reference genome digestion, read mapping for single 
# and paired DNA reads with bowtie from fastq files
# It expects a fastq file, pairdend, reference genome and digest pattern  as input.
# author: Fabian Buske
# date: Apr 2013

# messages to look out for -- relevant for the QC.sh script:
# QCVARIABLES

echo ">>>>> HiC readmapping with HiCUP "
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
        -f | --fastq )          shift; f=$1 ;; # fastq file
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

for MODULE in $MODULE_HICUP; do module load $MODULE; done  # save way to load modules that itself load other modules
export PATH=$PATH_HICUP:$PATH
module list
echo "PATH=$PATH"
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--bowtie      --\n "$(bowtie --version | head -n 1 )
[ -z "$(which bowtie)" ] && echo "[ERROR] no bowtie detected" && exit 1
echo -e "--perl        --\n "$(perl -v | grep "This is perl" )
[ -z "$(which perl)" ] && echo "[ERROR] no perl detected" && exit 1
echo -e "--HiCUP       --\n "$(hicup --version )
[ -z "$(which hicup)" ] && echo "[ERROR] no hicup detected" && exit 1
echo -e "--fit-hi-c    --\n "$(fit-hi-c.py --version)
[ -z "$(which fit-hi-c.py)" ] && echo "[ERROR] no fit-hi-c detected" && exit 1

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="parameters"

# get basename of f
n=${f##*/}

#output for this library
OUTDIR=${n/%$READONE.$FASTQ/}

# delete old bam files unless attempting to recover
if [ -z "$RECOVERFROM" ]; then
    [ -d $OUTDIR/$OUTDIR ] && rm -r $OUTDIR/$OUTDIR
    [ -e $OUTDIR/${n/%$READONE.$FASTQ/}.spline_pass1.q05.txt ] && rm $OUTDIR/${n/%$READONE.$FASTQ/}*.txt
fi

#is paired ?
if [ "$f" != "${f/$READONE/$READTWO}" ] && [ -e ${f/$READONE/$READTWO} ]; then
    PAIRED="1"
else
    echo "HiCUP requires paired fastq libraries" && exit 1
fi

#is ziped ?
ZCAT="zcat"
if [[ $f != *.gz ]]; then ZCAT="cat"; fi

if [ -z "$HICUP_RENZYMES" ]; then
   echo "[ERROR] No restriction enzyme given!" && exit 1
fi
ENZYMES=(${HICUP_RENZYMES//;/ })
ENZYME1=(${ENZYMES[0]//,/ })
ENZYME2=(${ENZYMES[1]//,/ })

DIGESTGENOME=""

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="recall files from tape"

if [ -n "$DMGET" ]; then
	dmget -a $(dirname $FASTA)/*
	dmget -a ${f/$READONE/"*"}
	dmget -a $OUTDIR/*
fi

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="digest reference"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

    FASTASUFFIX=${FASTA##*.}
    FASTABASE=${FASTA##*/}
    
    mkdir -p $OUTDIR/$OUTDIR
    cd $OUTDIR/$OUTDIR
    if [ ${#ENZYMES[@]} = 1 ]; then
       echo "Restriction Enzyme 1: ${ENZYME1[1]}:${ENZYME1[0]} "
       DIGESTGENOME=$OUTDIR/${FASTABASE/.$FASTASUFFIX/}_${ENZYME1[1]}_None.txt
       hicup_digester -g "${FASTABASE%.*}" -1 ${ENZYME1[0]} $FASTA
       mv Digest_* ${DIGESTGENOME}
    
    elif [ ${#ENZYMES[@]} = 2 ] && [ ! -e $OUTDIR/${FASTABASE/.$FASTASUFFIX/}_${ENZYME1[1]}_${ENZYME2[2]}.txt ]; then
       echo "Restriction Enzyme 1: ${ENZYME1[1]}:${ENZYME1[0]} "
       echo "Restriction Enzyme 2: ${ENZYME2[1]}:${ENZYME2[0]} "
       DIGESTGENOME=$OUTDIR/${FASTABASE/.$FASTASUFFIX/}_${ENZYME1[1]}_${ENZYME2[2]}.txt
       hicup_digester -g "${FASTABASE%.*}" -1 ${ENZYME1[0]} -2 ${ENZYME2[0]} $FASTA
       mv Digest_* ${DIGESTGENOME}
    else
       echo "[ERROR] Invalid number or pattern of enzyme digest patterns."
       exit 1
    fi
    cd $SOURCE
    
    # mark checkpoint
    if [ -f $DIGESTGENOME ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi

################################################################################
CHECKPOINT="create hicup conf script"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

    HICUP_CONF=$OUTDIR/${n/%$READONE.$FASTQ/.conf}
    
    cat /dev/null > $HICUP_CONF
    echo "#Number of threads to use" >> $HICUP_CONF
    echo "Threads: $CPU_HICUP" >> $HICUP_CONF
    echo "#Suppress progress updates | 0: off, 1: on" >> $HICUP_CONF
    echo "Quiet:0" >> $HICUP_CONF
    echo "#Retain all intermediate pipeline files | 0: off, 1: on" >> $HICUP_CONF
    echo "Keep:1" >> $HICUP_CONF
    echo "#Compress outputfiles | 0: off, 1: on" >> $HICUP_CONF
    echo "Zip:1" >> $HICUP_CONF
    echo "#Path to the alignment program Bowtie | include the executable Bowtie filename" >> $HICUP_CONF
    echo "Bowtie:$(which bowtie)" >> $HICUP_CONF
    echo "#Path to the reference genome indices" >> $HICUP_CONF
    echo "Index:${BOWTIE_INDEX}"  >> $HICUP_CONF
    echo "#Path to the genome digest file" >> $HICUP_CONF
    echo "DIGEST:$DIGESTGENOME" >> $HICUP_CONF
    echo "#FASTQ file format | phred33-quals, phred64-quals, solexa-quals or solexa1.3-quals" >> $HICUP_CONF
    echo "Format:phred33-quals" >> $HICUP_CONF
    echo "#Maximum di-tag length | optional parameter" >> $HICUP_CONF
    echo "#Longest:" >> $HICUP_CONF
    echo "#Minimum di-tag length | optional parameter" >> $HICUP_CONF
    echo "#Shortest:" >> $HICUP_CONF
    echo "#FASTQ files to be analysed, separating file pairs using the pipe '|' character" >> $HICUP_CONF
    echo "$f | ${f/$READONE/$READTWO} " >> $HICUP_CONF

    # mark checkpoint
    if [ -f $HICUP_CONF ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi

################################################################################
CHECKPOINT="execute hicup"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

    CURDIR=$(pwd)
    cd $OUTDIR/$OUTDIR
    RUN_COMMAND="$(which perl) $(which hicup) -c $HICUP_CONF"
    echo $RUN_COMMAND && eval $RUN_COMMAND
    
    cp -f $OUTDIR/$OUTDIR/hicup_deduplicater_summary_results_*.txt $OUTDIR/${n/%$READONE.$FASTQ/}_hicup_deduplicater_summary_results.txt
    cp -f $OUTDIR/$OUTDIR/hicup_filter_summary_results_*.txt $OUTDIR/${n/%$READONE.$FASTQ/}_hicup_filter_summary_results.txt
    cp -f $OUTDIR/$OUTDIR/hicup_mapper_summary_*.txt $OUTDIR/${n/%$READONE.$FASTQ/}_hicup_mapper_summary.txt
    cp -f $OUTDIR/$OUTDIR/hicup_truncater_summary_*.txt $OUTDIR/${n/%$READONE.$FASTQ/}_hicup_truncater_summary.txt
    ln -f -s $OUTDIR/uniques_${n/.$FASTQ/}_trunc_${n/%$READONE.$FASTQ/$READTWO}_trunc.bam $OUTDIR/${n/%$READONE.$FASTQ/}_uniques.bam
    
    cd $CURDIR

    # copy piecharts
    RUNSTATS=$OUT/runStats/$TASK_HICUP
    mkdir -p $RUNSTATS
    cp -f $OUTDIR/$OUTDIR/uniques_*_cis-trans.png $RUNSTATS/${n/%$READONE.$FASTQ/}_uniques_cis-trans.png
    cp -f $OUTDIR/$OUTDIR/*_ditag_classification.png $RUNSTATS/${n/%$READONE.$FASTQ/}_ditag_classification.png

    # mark checkpoint
    if [ -f $OUTDIR/uniques_${n/.$FASTQ/}_trunc_${n/%$READONE.$FASTQ/$READTWO}_trunc.bam ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi

################################################################################
CHECKPOINT="run fit-hi-c"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

    python ${NGSANE_BASE}/tools/hicupCountInteractions.py --verbose --genomeFragmentFile=${DIGESTGENOME} --outputDir=$OUTDIR/  $OUTDIR/${n/%$READONE.$FASTQ/}_uniques.bam
    cd $OUTDIR
    python $(which fit-hi-c.py) --mappabilityThres=2 --fragments=$OUTDIR/${n/%$READONE.$FASTQ/}_uniques.bam.fragmentLists --interactions=$OUTDIR/${n/%$READONE.$FASTQ/}_uniques.bam.contactCounts --lib=${n/%$READONE.$FASTQ/}
    cd $CURDIR
    
    awk '$7<=0.05' $OUTDIR/${n/%$READONE.$FASTQ/}.spline_pass1.significances.txt | sort -k7g > $OUTDIR/${n/%$READONE.$FASTQ/}.spline_pass1.q05.txt
    awk '$7<=0.05' $OUTDIR/${n/%$READONE.$FASTQ/}.spline_pass2.significances.txt | sort -k7g > $OUTDIR/${n/%$READONE.$FASTQ/}.spline_pass2.q05.txt
    
    $GZIP $OUTDIR/${n/%$READONE.$FASTQ/}*.significances.txt

    # mark checkpoint
    if [ -f $OUTDIR/${n/%$READONE.$FASTQ/}.spline_pass1.significances.txt ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi

################################################################################
echo ">>>>> readmapping with hicup (bowtie) - FINISHED"
echo ">>>>> enddate "`date`

