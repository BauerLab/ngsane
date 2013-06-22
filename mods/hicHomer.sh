#!/bin/bash

echo ">>>>> HiC analysis with homer"
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> hicHomer.sh $*"

function usage {
echo -e "usage: $(basename $0) -k NGSANE -f bam -o OUTDIR [OPTIONS]

Script running HIC HOMER pipeline tapping into bowtie2
It expects bam files, paired end, as input.

required:
  -k | --toolkit <path>     location of the NGSANE repository 
  -f | --bam <file>         bam file
  -o | --outdir <path>      output dir

options:
  -t | --threads <nr>       number of CPUs to use (default: 8)
"
exit
}
# QCVARIABLES,Resource temporarily unavailable

if [ ! $# -gt 3 ]; then usage ; fi

#DEFAULTS
THREADS=8

#INPUTS                                                                                                           
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the NGSANE repository                       
        -t | --threads )        shift; THREADS=$1 ;; # number of CPUs to use                                      
        -f | --bam )            shift; f=$1 ;; # bam file                                                       
        -o | --outdir )         shift; MYOUT=$1 ;; # output dir                                                     
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
for MODULE in $MODULE_HOMERHIC; do module load $MODULE; done  # save way to load modules that itself load other modules

export PATH=$PATH_HOMERHIC:$PATH
module list
echo "PATH=$PATH"
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)
echo -e "--R       --\n "$(R --version | head -n 3)
[ -z "$(which R)" ] && echo "[ERROR] no R detected" && exit 1
echo -e "--homer   --\n "$(makeTagDirectory -version | head -n 1)
[ -z "$(which makeTagDirectory)" ] && echo "[ERROR] homer not detected" && exit 1

# get basename of f
n=${f##*/}

#is paired ?                                                                                                      
if [ -e ${f/$READONE/$READTWO} ]; then
    PAIRED="1"
else
    PAIRED="0"
fi

FASTASUFFIX=${FASTA##*.}

if [ -n "$DMGET" ]; then
	echo "********** reacall files from tape"
	dmget -a ${f/$READONE/"*"}
	dmls -l ${f/$READONE/"*"}
fi

if [ $PAIRED == "0" ]; then 
    echo "[ERROR] paired library required for HIC analysis"
    exit 1
fi

echo "********* makeTagDirectory" 
RUN_COMMAND="makeTagDirectory $MYOUT/${n/'_'$READONE.$ASD.bam/_tagdir_unfiltered} $f ${f/'_'$READONE/'_'$READTWO} $HOMER_HIC_TAGDIR_OPTIONS"
echo $RUN_COMMAND
eval $RUN_COMMAND

cp -r $MYOUT/${n/'_'$READONE.$ASD.bam/_tagdir_unfiltered} $MYOUT/${n/'_'$READONE.$ASD.bam/_tagdir_filtered}

RUN_COMMAND="makeTagDirectory $MYOUT/${n/'_'$READONE.$ASD.bam/_tagdir_filtered} -update $HOMER_HIC_TAGDIR_OPTIONS"

echo $RUN_COMMAND
eval $RUN_COMMAND

echo "********* analyzeHiC" 




echo ">>>>> HiC analysis with homer - FINISHED"
echo ">>>>> enddate "`date`

