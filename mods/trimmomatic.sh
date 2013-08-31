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
echo ">>>>> $(basename $0) $*"

while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of NGSANE
        -f | --file )           shift; f=$1 ;; # fastq file
        --recover-from )        shift; RECOVERFROM=$1 ;; # attempt to recover from log file
        -h | --help )           usage ;;
        * )                     echo "don't understand "$1
    esac
    shift
done

. $CONFIG
. ${NGSANE_BASE}/conf/header.sh
. $CONFIG

################################################################################
CHECKPOINT="programs"

for MODULE in $MODULE_TRIMMOMATIC; do module load $MODULE; done  # save way to load modules that itself load other modules
export PATH=$PATH_TRIMMOMATIC:$PATH;
module list
echo "PATH=$PATH"
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)
PATH_TRIMMOMATIC=$(dirname $(which trimmomatic.jar))
echo -e "--JAVA        --\n" $(java -version 2>&1)
[ -z "$(which java)" ] && echo "[ERROR] no java detected" && exit 1
echo -e "--trimmomatic --\n " $(which $PATH_TRIMMOMATIC/trimmomatic.jar)
[ ! -f $PATH_TRIMMOMATIC/trimmomatic.jar ] && echo "[ERROR] no trimmomatic detected" && exit 1

echo "[NOTE] set java parameters"
JAVAPARAMS="-Xmx"$(expr $MEMORY_TRIMMOMATIC - 1 )"G -Djava.io.tmpdir="$TMP
echo "JAVAPARAMS "$JAVAPARAMS

echo -e "\n********* $CHECKPOINT"
################################################################################
CHECKPOINT="parameters"

# get basename of f
n=${f##*/}

#is paired ?
if [ "$f" != "${f/$READONE/$READTWO}" ] && [ -e ${f/$READONE/$READTWO} ]; then
    echo "[NOTE] PAIRED library"
    PAIRED="1"
else
    echo "[NOTE] SINGLE library"
    PAIRED="0"
fi

# check if trimming steps are set
if [ -z "$TRIMMOMATICSTEPS" ]; then
    echo "[ERROR] no trimming steps specified: TRIMMOMATICSTEPS" && exit 1
fi

# get encoding
FASTQ_ENCODING=$(zcat $f |  awk 'NR % 4 ==0' | python $NGSANE_BASE/tools/GuessFastqEncoding.py |  tail -n 1)
if [[ "$FASTQ_ENCODING" == *Sanger* ]]; then
    TRIMMOMATICADDPARAM="$TRIMMOMATICADDPARAM -phred33"    
elif [[ "$FASTQ_ENCODING" == *Illumina* ]]; then
    TRIMMOMATICADDPARAM="$TRIMMOMATICADDPARAM -phred64"
elif [[ "$FASTQ_ENCODING" == *Solexa* ]]; then
    TRIMMOMATICADDPARAM="$TRIMMOMATICADDPARAM -phred64"
else
    echo "[ERROR] cannot detect/don't understand fastq format: $FASTQ_ENCODING" && exit 1
fi
echo "[NOTE] $FASTQ_ENCODING fastq format detected"

FASTQDIR=$(basename $(dirname $f))
o=${f/$FASTQDIR/$FASTQDIR"_"$TASKTRIMMOMATIC}
FASTQDIRTRIM=$(dirname $o)

echo $FASTQDIRTRIM
if [ ! -d $FASTQDIRTRIM ]; then mkdir -p $FASTQDIRTRIM; fi
echo $f "->" $o
if [ "$PAIRED" = "1" ]; then echo ${f/$READONE/$READTWO} "->" ${o/$READONE/$READTWO} ; fi


echo -e "\n********* $CHECKPOINT"
################################################################################
CHECKPOINT="recall files from tape"

if [ -n "$DMGET" ]; then
    dmget -a ${f/$READONE/"*"}
fi

echo -e "\n********* $CHECKPOINT"
################################################################################
CHECKPOINT="trim"    

if [[ -n "$RECOVERFROM" ]] && [[ $(grep "********* $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
    # Paired read
    if [ "$PAIRED" = "1" ]
    then
        RUN_COMMAND="java -jar $PATH_TRIMMOMATIC/trimmomatic.jar PE $TRIMMOMATICADDPARAM -threads $CPU_TRIMMOMATIC $f ${f/$READONE/$READTWO} $o ${o/$READONE/${READONE}_unpaired} ${o/$READONE/$READTWO} ${o/$READONE/${READTWO}_unpaired} $TRIMMOMATICSTEPS &> ${o/%$READONE.$FASTQ/}.log"
    else
        RUN_COMMAND="java -jar $PATH_TRIMMOMATIC/trimmomatic.jar SE $TRIMMOMATICADDPARAM -threads $CPU_TRIMMOMATIC $f $o $TRIMMOMATICSTEPS &> ${o/%$READONE.$FASTQ/}.log"
    fi
    echo $RUN_COMMAND && eval $RUN_COMMAND

    # mark checkpoint
    [ -f $o ] && echo -e "\n********* $CHECKPOINT" && unset RECOVERFROM
fi

################################################################################
echo ">>>>> readtrimming with TRIMMOMATIC - FINISHED"
echo ">>>>> enddate "`date`

