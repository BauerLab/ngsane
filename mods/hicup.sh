#!/bin/bash -e

# Script running hicup including reference genome digestion, read mapping for single 
# and paired DNA reads with bowtie from fastq files
# It expects a fastq file, paired-end, reference genome and digest pattern as input.
# author: Fabian Buske
# date: Jan 2014

# messages to look out for -- relevant for the QC.sh script:
# QCVARIABLES,Resource temporarily unavailable
# RESULTFILENAME <DIR>/<TASK>/<SAMPLE>.spline_pass1.q05.txt

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
echo -e "--Python      --\n" $(python --version)
[ -z "$(which python)" ] && echo "[ERROR] no python detected" && exit 1
echo -e "--Python libs --\n "$(yolk -l)
echo -e "--fit-hi-c    --\n "$(python $(which fit-hi-c.py) --version | head -n 1)
[ -z "$(which fit-hi-c.py)" ] && echo "[WARN] no fit-hi-c detected" && exit 1

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="parameters"

# get basename of f
n=${f##*/}
SAMPLE=${n/%$READONE.$FASTQ/}

if [ -z "$FASTA" ]; then
    echo "[ERROR] no reference provided (FASTA)"
    exit 1
fi

if [[ ! -e ${FASTA%.*}.1.ebwt ]]; then
    echo "[ERROR] Bowtie index not detected. Exeute bowtieIndex.sh first"
    exit 1
fi

# delete old bam files unless attempting to recover
if [ -z "$RECOVERFROM" ]; then
    [ -d $OUTDIR/$SAMPLE ] && rm -r $OUTDIR/$SAMPLE
    [ -e $OUTDIR/${SAMPLE}.spline_pass1.q05.txt ] && rm $OUTDIR/${SAMPLE}*.txt
fi

#is paired ?
if [ "$f" != "${f/%$READONE.$FASTQ/$READTWO.$FASTQ}" ] && [ -e ${f/%$READONE.$FASTQ/$READTWO.$FASTQ} ]; then
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

DIGESTGENOME=$OUTDIR/digested_genome.txt

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
CHECKPOINT="execute hicup"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

    [ -d $OUTDIR/$SAMPLE ] && rm -r $OUTDIR/$SAMPLE
    mkdir -p $OUTDIR/$SAMPLE
    
    cat /dev/null > $OUTDIR/${SAMPLE}.conf
    echo "#Number of threads to use" >> $OUTDIR/${SAMPLE}.conf
    echo "Threads: $CPU_HICUP" >> $OUTDIR/${SAMPLE}.conf
    echo "#Suppress progress updates | 0: off, 1: on" >> $OUTDIR/${SAMPLE}.conf
    echo "Quiet:0" >> $OUTDIR/${SAMPLE}.conf
    echo "#Retain all intermediate pipeline files | 0: off, 1: on" >> $OUTDIR/${SAMPLE}.conf
    echo "Keep:1" >> $OUTDIR/${SAMPLE}.conf
    echo "#Compress outputfiles | 0: off, 1: on" >> $OUTDIR/${SAMPLE}.conf
    echo "Zip:1" >> $OUTDIR/${SAMPLE}.conf
    echo "#Path to the alignment program Bowtie | include the executable Bowtie filename" >> $OUTDIR/${SAMPLE}.conf
    echo "Bowtie:$(which bowtie)" >> $OUTDIR/${SAMPLE}.conf
    echo "#Path to the reference genome indices" >> $OUTDIR/${SAMPLE}.conf
    echo "Index:${FASTA%.*}"  >> $OUTDIR/${SAMPLE}.conf
    echo "#Path to the genome digest file" >> $OUTDIR/${SAMPLE}.conf
    echo "DIGEST:$DIGESTGENOME" >> $OUTDIR/${SAMPLE}.conf
    echo "#FASTQ file format | phred33-quals, phred64-quals, solexa-quals or solexa1.3-quals" >> $OUTDIR/${SAMPLE}.conf
    echo "Format:phred33-quals" >> $OUTDIR/${SAMPLE}.conf
    echo "#Maximum di-tag length | optional parameter" >> $OUTDIR/${SAMPLE}.conf
    echo "#Longest:" >> $OUTDIR/${SAMPLE}.conf
    echo "#Minimum di-tag length | optional parameter" >> $OUTDIR/${SAMPLE}.conf
    echo "#Shortest:" >> $OUTDIR/${SAMPLE}.conf
    echo "#FASTQ files to be analysed, separating file pairs using the pipe '|' character" >> $OUTDIR/${SAMPLE}.conf
    echo "$f | ${f/%$READONE.$FASTQ/$READTWO.$FASTQ} " >> $OUTDIR/${SAMPLE}.conf

    cd $OUTDIR/$SAMPLE
    RUN_COMMAND="$(which perl) $(which hicup) --config $OUTDIR/${SAMPLE}.conf"
    echo $RUN_COMMAND && eval $RUN_COMMAND
    cd $SOURCE
    
    cp -f $OUTDIR/$SAMPLE/hicup_deduplicater_summary_results_*.txt $OUTDIR/${SAMPLE}_hicup_deduplicater_summary_results.txt
    cp -f $OUTDIR/$SAMPLE/hicup_filter_summary_results_*.txt $OUTDIR/${SAMPLE}_hicup_filter_summary_results.txt
    cp -f $OUTDIR/$SAMPLE/hicup_mapper_summary_*.txt $OUTDIR/${SAMPLE}_hicup_mapper_summary.txt
    cp -f $OUTDIR/$SAMPLE/hicup_truncater_summary_*.txt $OUTDIR/${SAMPLE}_hicup_truncater_summary.txt
    ln -f -s $SAMPLE/uniques_${n/.$FASTQ/}_trunc_${n/%$READONE.$FASTQ/$READTWO}_trunc.bam $OUTDIR/${SAMPLE}_uniques.bam

    # copy piecharts
    RUNSTATS=$OUT/runStats/$TASK_HICUP
    mkdir -p $RUNSTATS
    cp -f $OUTDIR/$SAMPLE/uniques_*_cis-trans.png $RUNSTATS/${SAMPLE}_uniques_cis-trans.png
    cp -f $OUTDIR/$SAMPLE/*_ditag_classification.png $RUNSTATS/${SAMPLE}_ditag_classification.png

    # mark checkpoint
    if [ -f $OUTDIR/${SAMPLE}_uniques.bam ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi

################################################################################
CHECKPOINT="run fit-hi-c"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

    if [ -z "$(which fit-hi-c.py)" ]; then
        echo "[NOTE] Skip significance testing (fit-hi-c not found)."
        echo -e "\n********* $CHECKPOINT\n"
    else
    
        python ${NGSANE_BASE}/tools/hicupCountInteractions.py --verbose --genomeFragmentFile=$DIGESTGENOME --outputDir=$OUTDIR/  $OUTDIR/${SAMPLE}_uniques.bam
        cd $OUTDIR
        python $(which fit-hi-c.py) $FITHICADDPARAM --fragments=$OUTDIR/${SAMPLE}_uniques.bam.fragmentLists --interactions=$OUTDIR/${SAMPLE}_uniques.bam.contactCounts --lib=${SAMPLE} 
        cd $SOURCE
        
        awk '$7<=0.05' $OUTDIR/${SAMPLE}.spline_pass1.significances.txt | sort -k7g > $OUTDIR/${SAMPLE}.spline_pass1.q05.txt
        
        $GZIP $OUTDIR/${SAMPLE}*.significances.txt
    
        # mark checkpoint
        if [ -f $OUTDIR/${SAMPLE}.spline_pass1.significances.txt ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
    fi
fi

################################################################################
[ -e $OUTDIR/${SAMPLE}.spline_pass1.q05.txt.dummy ] && rm $OUTDIR/${SAMPLE}.spline_pass1.q05.txt.dummy
echo ">>>>> readmapping with hicup (bowtie) - FINISHED"
echo ">>>>> enddate "`date`

