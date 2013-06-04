#!/bin/bash

echo ">>>>> HiC correlation analysis with hiclib "
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> hiclibCorrelate.sh $*"

function usage {
echo -e "usage: $(basename $0) -k CONFIG -f FASTQ -r REFERENCE -e ENZYMES -o OUTDIR [OPTIONS]

The script correlates all hiclib library in a pairwise manner.

required:
  -k | --toolkit <path>     location of the CONFIG repository 
  -f | --files <FILE>+      input files, separated by comma
  -o | --outdir <path>      output dir
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
FASTQNAME=""
ENZYME=""
QUAL="" # standard Sanger

#INPUTS                                                                                                           
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the NGSANE repository
	-f | --files ) 		shift; FILES=$1 ;; # input files
        -o | --outdir )         shift; MYOUT=$1 ;; # output dir                                                     
        --fastqName )           shift; FASTQNAME=$1 ;; #(name of fastq or fastq.gz)
        -h | --help )           usage ;;
        * )                     echo "don't understand "$1
    esac
    shift
done

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
echo -e "--Python libs --\n "$(yolk -l)

# find files for dataset
echo $FILES
OLDFS=$IFS
IFS=","
DATASETS=""
for f in $FILES; do
    # get basename of f
    n=${f##*/}
    n=${n/_$READONE.$FASTQ/}
    # get directory
    d=$(dirname $f)
    d=${d##*/}
    # get hdf5 file
    FILE=$(ls $SOURCE/$d/hiclib/*_$n-fragment_dataset.hdf5)
    # add to dataset
    if [ -n "$FILE" ]; then 
	DATASETS="${DATASETS[@]} ${FILE[@]}"
    fi
done
IFS=$OLDFS

if [ -n "$DMGET" ]; then
	echo "********** reacall files from tape"
	dmget -a $MYOUT/*
fi
mkdir -p $MYOUT

echo "********** hiclib call"
PARAMS="--gapFile=$HICLIB_GAPFILE \
   --referenceGenome=$FASTA "

python ${NGSANE_BASE}/tools/hiclibCorrelate.py ${PARAMS} --outputDir=$OUT/runStats/i$TASKHICLIB --tmpDir=$TMP --verbose $DATASETS

echo ">>>>> correlation with hiclib - FINISHED"
echo ">>>>> enddate "`date`
