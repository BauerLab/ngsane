#!/bin/bash

echo ">>>>> HiC analysis with hiclib "
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> hiclibMapping.sh $*"

function usage {
echo -e "usage: $(basename $0) -k NGSANE -f FASTQ -r REFERENCE -e ENZYMES -o OUTDIR [OPTIONS]

Script running hiclib pipeline tapping into bowtie2
It expects a fastq file, paired end, reference genome and digest pattern  as input.

required:
  -k | --toolkit <path>     location of the NGSANE repository 
  -f | --fastq <file>       fastq file
  -e | --enzymes <name>     restriction enzyme (one per library) seperated by comma 
                            see http://biopython.org/DIST/docs/api/Bio.Restriction-module.html
  -o | --outdir <path>      output dir

options:
  -t | --threads <nr>       number of CPUs to use (default: 8)
  --fastqName               name of fastq file ending (fastq.gz)
  --oldIllumina
"
exit
}

# Script to run a hic analysis based on the hiclib framework.
# It takes comma-seprated list of files containing short sequence reads in fasta or fastq format and bowtie index files as input.
# It produces output files: read alignments in .bam format and other files.
# author: Fabian Buske
# date: April 2013

# QCVARIABLES,Resource temporarily unavailable

if [ ! $# -gt 3 ]; then usage ; fi

#DEFAULTS
THREADS=8
FASTQNAME=""
ENZYME=""
QUAL="" # standard Sanger

#INPUTS                                                                                                           
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the NGSANE repository                       
        -t | --threads )        shift; THREADS=$1 ;; # number of CPUs to use                                      
        -f | --fastq )          shift; f=$1 ;; # fastq file                                                       
        -e | --enzymes )        shift; ENZYME=$1 ;; # digestion patterns
        -o | --outdir )         shift; MYOUT=$1 ;; # output dir                                                     
        --fastqName )           shift; FASTQNAME=$1 ;; #(name of fastq or fastq.gz)
        -h | --help )           usage ;;
        * )                     echo "don't understand "$1
    esac
    shift
done

if [ -z "$ENZYME" ]; then
	echo "[ERROR] restriction enzyme not specified"
	exit 1
fi

if [ -z "BOWTIE2INDEX" ]; then
	echo "[ERROR] bowtie index not specified"
	exit 1
fi

#PROGRAMS
. $CONFIG
. ${NGSANE_BASE}/conf/header.sh
. $CONFIG

echo "********** programs"
for MODULE in $MODULE_HICLIB; do module load $MODULE; done  # save way to load modules that itself load other modules

export PATH=$PATH_HICLIB:$PATH
module list
echo "PATH=$PATH"
echo -e "--Python      --\n" $(python --version)
[ -z "$(which python)" ] && echo "[ERROR] no python detected" && exit 1
echo -e "--Python libs --\n "$(yolk -l)

# get basename of f
n=${f##*/}

# delete old bam file                                                                       
#if [ -e $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam} ]; then rm $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}; fi

#is paired ?                                                                                                      
if [ -e ${f/$READONE/$READTWO} ]; then
    PAIRED="1"
else
    echo "[ERROR] hiclib requires paired-end fastq files. Could not find ${f/$READONE/$READTWO}"
    exit 1
fi

if [ -n "$DMGET" ]; then
	echo "********** reacall files from tape"
	dmget -a $(dirname $FASTA)/*
	dmls -l $FASTA*
	dmget -a ${f/$READONE/"*"}
	dmls -l ${f/$READONE/"*"}
fi

echo "********* reads" 
FASTQNAME=${f##*/}
FASTQSUFFIX=${f##*.}
READS="$f ${f/$READONE/$READTWO}"

if [ -f "$MYOUT/${FASTQNAME/.$FASTQSUFFIX/_R1.bam.25}"  ]; then
    echo "[NOTE] using bam files from previous run"
    LIBONE=${FASTQNAME/.$FASTQSUFFIX/_R1.bam}
    LIBTWO=${FASTQNAME/$READONE/$READTWO}
    LIBTWO=${LIBTWO/.$FASTQSUFFIX/_R2.bam}
    READS="$MYOUT/$LIBONE $MYOUT/$LIBTWO"
    FASTQSUFFIX=bam
    echo $READS
else
    echo "[NOTE] mapping data from scratch"
fi
echo "********** hiclib call"

EXPERIMENT=$(echo ${ENZYME}_${FASTQNAME/$READONE.$FASTQ/} | sed 's/_*$//g')
PARAMS="--restrictionEnzyme=$ENZYME \
   --experimentName=$EXPERIMENT \
   --gapFile=$HICLIB_GAPFILE \
   --referenceGenome=$FASTA \
   --index=$BOWTIE2_INDEX"

if [ "$FASTQSUFFIX" = "sra" ]; then
	PARAMS="$PARAMS --inputFormat=sra --sra-reader=$(which fastq-dump)"
elif [ "$FASTQSUFFIX" = "bam" ]; then
	PARAMS="$PARAMS --inputFormat=bam"
else
	PARAMS="$PARAMS --inputFormat=fastq"
fi

if [ -n "$HICLIB_READLENGTH" ]; then
	PARAMS="$PARAMS --readLength $HICLIB_READLENGTH"
fi

# run hiclib.py
RUN_COMMAND="python ${NGSANE_BASE}/tools/hiclibMapping.py ${PARAMS} --bowtie=$(which bowtie2) --cpus=$THREADS --outputDir=$MYOUT --tmpDir=$TMP --verbose $READS &> $MYOUT/$EXPERIMENT.hiclib.log"
echo $RUN_COMMAND
eval $RUN_COMMAND

#rm -f $MYOUT/*$READONE.bam.*  $MYOUT/*$READTWO.bam.*

# copy heatmap
RUNSTATS=$OUT/runStats/$TASKHICLIB
mkdir -p $RUNSTATS
mv -f $MYOUT/$EXPERIMENT*.pdf $RUNSTATS
mv -f $MYOUT/$EXPERIMENT.hiclib.log $RUNSTATS

echo ">>>>> readmapping with hiclib (Bowtie2) - FINISHED"
echo ">>>>> enddate "`date`
