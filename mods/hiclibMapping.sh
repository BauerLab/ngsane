#!/bin/bash -e

# Script running hiclib pipeline tapping into bowtie2
# It expects a fastq file, paired end, reference genome and digest pattern  as input.
# author: Fabian Buske
# date: April 2013

echo ">>>>> HiC analysis with hiclib "
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"

function usage {
echo -e "usage: $(basename $0) -k NGSANE -f FASTQ -r REFERENCE -o OUTDIR [OPTIONS]"
exit
}

# QCVARIABLES,Resource temporarily unavailable
# RESULTFILENAME <DIR>/<TASK>/${HICLIB_RENZYMES}_<SAMPLE>-mapped_reads.hdf5

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
for MODULE in $MODULE_HICLIB; do module load $MODULE; done  # save way to load modules that itself load other modules

export PATH=$PATH_HICLIB:$PATH
module list
echo "PATH=$PATH"

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--Python      --\n" $(python --version)
[ -z "$(which python)" ] && echo "[ERROR] no python detected" && exit 1
echo -e "--Python libs --\n "$(yolk -l)

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="parameters"

# get basename of f
n=${f##*/}

if [ -z "$FASTA" ]; then
    echo "[ERROR] no reference provided (FASTA)"
    exit 1
fi

if [[ ! -e ${FASTA%.*}.1.bt2 ]]; then
    echo "[ERROR] Bowtie2 index not detected. Exeute bowtie2Index.sh first"
    exit 1
fi

#is paired ?                                                                                                      
if [ "$f" != "${f/%$READONE.$FASTQ/$READTWO.$FASTQ}" ] && [ -e ${f/%$READONE.$FASTQ/$READTWO.$FASTQ} ]; then
    PAIRED="1"
else
    echo "[ERROR] hiclib requires paired-end fastq files. Could not find ${f/%$READONE.$FASTQ/$READTWO.$FASTQ}" && exit 1
fi

if [ -z "$HICLIB_RENZYMES" ]; then
	echo "[ERROR] restriction enzyme not specified" && exit 1
fi

if [ -z "BOWTIE2INDEX" ]; then
	echo "[ERROR] bowtie2 index not specified" && exit 1
fi

READS="$f ${f/%$READONE.$FASTQ/$READTWO.$FASTQ}"
if [ -f "$OUTDIR/${n/.$FASTQ/_R1.bam.25}"  ]; then
    echo "[NOTE] using bam files from previous run"
    LIBONE=${n/.$FASTQ/_R1.bam}
    LIBTWO=${n/$READONE.$FASTQ/${READTWO}_R2.bam}
#    LIBTWO=${LIBTWO/.$FASTQ/_R2.bam}
    READS="$OUTDIR/$LIBONE $OUTDIR/$LIBTWO"
    FASTQ=bam
    echo $READS
else
    echo "[NOTE] mapping data from scratch"
fi

#EXPERIMENT=$(echo ${HICLIB_RENZYMES}_${n/%$READONE.$FASTQ/} | sed 's/_*$//g')
EXPERIMENT=$(echo ${HICLIB_RENZYMES}_${n/%$READONE.$FASTQ/})

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="recall files from tape"

if [ -n "$DMGET" ]; then
	dmget -a $(dirname $FASTA_CHROMDIR)/*
	dmget -a $(dirname $FASTA)/*
	dmget -a ${f/$READONE/"*"}
	dmget -a $OUTDIR/*
fi

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="run hiclib"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
    
    PARAMS="--restrictionEnzyme=$HICLIB_RENZYMES \
       --experimentName=$EXPERIMENT \
       --gapFile=$HICLIB_GAPFILE \
       --referenceGenome=$FASTA_CHROMDIR \
       --index=${FASTA%.*}"
    
    if [ "$FASTQ" = "sra" ]; then
    	PARAMS="$PARAMS --inputFormat=sra --sra-reader=$(which fastq-dump)"
    elif [ "$FASTQ" = "bam" ]; then
    	PARAMS="$PARAMS --inputFormat=bam"
    else
    	PARAMS="$PARAMS --inputFormat=fastq"
    fi
    
    if [ -n "$HICLIB_READLENGTH" ]; then
    	PARAMS="$PARAMS --readLength $HICLIB_READLENGTH"
    fi
    
    RUN_COMMAND="python ${NGSANE_BASE}/tools/hiclibMapping.py ${PARAMS} --bowtie=$(which bowtie2) --cpus=$CPU_HICLIB --outputDir=$OUTDIR --tmpDir=$TMP --verbose $READS &> $OUTDIR/$EXPERIMENT.hiclib.log"
    echo $RUN_COMMAND && eval $RUN_COMMAND

    RUNSTATS=$OUT/runStats/$TASK_HICLIB
    mkdir -p $RUNSTATS
    mv -f $OUTDIR/$EXPERIMENT*.pdf $RUNSTATS
    mv -f $OUTDIR/$EXPERIMENT.hiclib.log $RUNSTATS

    if [ -z "$HICLIB_KEEPBAM" ]; then 
        rm $OUTDIR/*$READONE.bam.*  $OUTDIR/*$READTWO.bam.*
    fi

    # mark checkpoint
    if [ -e $OUTDIR/${EXPERIMENT}-mapped_reads.hdf5 ] ;then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi

################################################################################
[ -e $OUTDIR/${EXPERIMENT}-mapped_reads.hdf5.dummy ] && rm $OUTDIR/${EXPERIMENT}-mapped_reads.hdf5.dummy
echo ">>>>> readmapping with hiclib (Bowtie2) - FINISHED"
echo ">>>>> enddate "`date`
