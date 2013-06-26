#!/bin/bash

echo ">>>>> Create wig files with wiggler"
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> wiggler.sh $*"

function usage {
echo -e "usage: $(basename $0) -k NGSANE -f bam -o OUTDIR [OPTIONS]

Script running wiggler on bam files

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
MEMORY=2

#INPUTS                                                                                                           
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the NGSANE repository                       
        -t | --threads )        shift; THREADS=$1 ;; # number of CPUs to use                                      
        -m | --memory )         shift; MEMORY=$1 ;; # memory used 
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
for MODULE in $MODULE_WIGGLER; do module load $MODULE; done  # save way to load modules that itself load other modules

export PATH=$PATH_WIGGLER:$PATH
module list
echo "PATH=$PATH"
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)
echo -e "--wiggler  --\n "$(align2rawsignal 2>&1 | head -n 3 | tail -n 1)
[ -z "$(which align2rawsignal)" ] && echo "[ERROR] wiggler not detected (align2rawsignal)" && exit 1

# get basename of f
n=${f##*/}

#is paired ?                                                                                                      
if [ -e ${f/$READONE/$READTWO} ]; then
    PAIRED="1"
else
    PAIRED="0"
fi

# check UMAP folder exist
if [ -z "$WIGGLER_UMAPDIR" ] || [ ! -d ${WIGGLER_UMAPDIR} ]; then
    echo "[ERROR] umap dir not specified not non-existant (WIGGLER_UMAPDIR)"
    exit 1
fi

# check output format
if [ -z "$WIGGLER_OUTPUTFORMAT" ]; then
    echo "[ERROR] wiggler output format not set" && exit 1
else if [ "$WIGGLER_OUTPUTFORMAT" ne "bg"] && [ "$WIGGLER_OUTPUTFORMAT" ne "wig"] && [ "$WIGGLER_OUTPUTFORMAT" ne "mat"]; then
    echo "[ERROR] wiggler output format not known" && exit 1
fi

echo "********** get input"
for d in ${DIR[@]}; do
    FILES=$FILES" "$( ls $OUT/$d/$TASKBOWTIE2/*$ASD.bam )
done
echo $FILES

if [ -n "$DMGET" ]; then
	echo "********** reacall files from tape"
	dmget -a $FILES
	dmls -l $FILES
fi

echo "********* align2rawsignal" 
INPUTS=""
for FILE in FILES; do
    INPUTS="${INPUTS} -i=${$FILE}"
done

RUN_COMMAND="align2rawsignal $WIGGLERADDPARAMS -of=$WIGGLER_OUTPUTFORMAT ${INPUTS} -s=${FASTA} -u=${WIGGLER_UMAPDIR} -v=${MYOUT}/wiggler-${n/$READONE/}.log -o=${MYOUT}/${n/$READONE/}.$WIGGLER_OUTPUTFORMAT -mm=$MEMORY"
echo $RUN_COMMAND
eval $RUN_COMMAND


echo ">>>>> wiggler - FINISHED"
echo ">>>>> enddate "`date`

