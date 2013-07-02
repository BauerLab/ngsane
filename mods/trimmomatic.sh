#!/bin/bash -e

# Script to trim adapters using TRIMMOMATIC
# It takes a <Run>/*.$FASTQ[.gz] file and gets the file containing the contaminats
# via config and writes out <Run>_trimmomatic/*.$FASTQ[.gz]
#
# author: Fabian Buske
# date: June. 2013

# messages to look out for -- relevant for the QC.sh script:
# QCVARIABLES,

echo ">>>>> readtrimming with TRIMMOMATIC"
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> trimgalore.sh $*"

while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of NGSANE
        -f | --file )           shift; f=$1 ;; # fastq file
        -h | --help )           usage ;;
        * )                     echo "don't understand "$1
    esac
    shift
done

. $CONFIG
. ${NGSANE_BASE}/conf/header.sh
. $CONFIG

JAVAPARAMS="-Xmx"$(expr $MEMORY_TRIMMOMATIC - 1 )"G -Djava.io.tmpdir="$TMP
echo "JAVAPARAMS "$JAVAPARAMS

echo "********** programs"
for MODULE in $MODULE_TRIMMOMATIC; do module load $MODULE; done  # save way to load modules that itself load other modules
export PATH=$PATH_TRIMMOMATIC:$PATH;
module list
echo "PATH=$PATH"
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)
PATH_TRIMMOMATIC=$(dirname $(which trimmomatic.jar))
echo -e "--JAVA        --\n" $(java $JAVAPARAMS -version 2>&1)
[ -z "$(which java)" ] && echo "[ERROR] no java detected" && exit 1
echo -e "--trimmomatic --\n " $(which $PATH_TRIMMOMATIC/trimmomatic.jar)
[ ! -f $PATH_TRIMMOMATIC/trimmomatic.jar ] && echo "[ERROR] no trimmomatic detected" && exit 1

# check if trimming steps are set
if [ -z "$TRIMMOMATICSTEPS" ]; then
    echo "[ERROR] no trimming steps specified: TRIMMOMATICSTEPS"
    exit 1
fi

#SAMPLENAME
# get basename of f
n=${f##*/}

#is paired ?
if [ -e ${f/$READONE/$READTWO} ]; then
    echo "[NOTE] PAIRED library"
    PAIRED="1"
else
    echo "[NOTE] SINGLE library"
    PAIRED="0"
fi

FASTQDIR=$(basename $(dirname $f))
o=${f/$FASTQDIR/$FASTQDIR"_"$TASKTRIMMOMATIC}
FASTQDIRTRIM=$(dirname $o)

if [ -n "$DMGET" ]; then
    dmget -a ${f/$READONE/"*"}
fi

echo $FASTQDIRTRIM
if [ ! -d $FASTQDIRTRIM ]; then mkdir -p $FASTQDIRTRIM; fi
echo $f "->" $o
if [ "$PAIRED" = "1" ]; then echo ${f/$READONE/$READTWO} "->" ${o/$READONE/$READTWO} ; fi

echo "********** trim"
# Paired read
if [ "$PAIRED" = "1" ]
then
    RUN_COMMAND="java -jar $PATH_TRIMMOMATIC/trimmomatic.jar PE $TRIMMOMATICADDPARAM -threads $CPU_TRIMMOMATIC $f ${f/_$READONE/_$READTWO} $o ${o/_$READONE/_${READONE}_unpaired} ${o/_$READONE/_$READTWO} ${o/_$READONE/_${READTWO}_unpaired} $TRIMMOMATICSTEPS &> ${o/_$READONE.$FASTQ/}.log"
else
    RUN_COMMAND="java -jar $PATH_TRIMMOMATIC/trimmomatic.jar SE $TRIMMOMATICADDPARAM -threads $CPU_TRIMMOMATIC $f $o $TRIMMOMATICSTEPS &> ${o/_$READONE.$FASTQ/}.log"
fi
echo $RUN_COMMAND
eval $RUN_COMMAND

echo ">>>>> readtrimming with TRIMMOMATIC - FINISHED"
echo ">>>>> enddate "`date`

