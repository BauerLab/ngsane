#!/bin/bash -e

# Script to trim adapters using TRIMMOMATIC
# It takes a <Run>/*.$FASTQ[.gz] file and gets the file containing the contaminats
# via config and writes out <Run>_trimmomatic/*.$FASTQ[.gz]
#
# author: Fabian Buske
# date: June. 2013

# messages to look out for -- relevant for the QC.sh script:
# QCVARIABLES,
# RESULTFILENAME fastq/<DIR>"_"$TASK_TRIMMOMATIC/<SAMPLE>$READONE.$FASTQ

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
        --recover-from )        shift; NGSANE_RECOVERFROM=$1 ;; # attempt to recover from log file
        -h | --help )           usage ;;
        * )                     echo "don't understand "$1
    esac
    shift
done

. $CONFIG
. ${NGSANE_BASE}/conf/header.sh
. $CONFIG

################################################################################
NGSANE_CHECKPOINT_INIT "programs"

# save way to load modules that itself loads other modules
hash module 2>/dev/null && for MODULE in $MODULE_TRIMMOMATIC; do module load $MODULE; done && module list 

export PATH=$PATH_TRIMMOMATIC:$PATH;
echo "PATH=$PATH"
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)
[ -z "$PATH_TRIMMOMATIC" ] && PATH_TRIMMOMATIC=$(dirname $(which trimmomatic.jar))

echo "[NOTE] set java parameters"
JAVAPARAMS="-Xmx"$(python -c "print int($MEMORY_TRIMMOMATIC*0.75)")"g -Djava.io.tmpdir="$TMP"  -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1" 
unset _JAVA_OPTIONS
echo "JAVAPARAMS "$JAVAPARAMS

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--JAVA        --\n" $(java -Xmx200m -version 2>&1)
[ -z "$(which java)" ] && echo "[ERROR] no java detected" && exit 1
echo -e "--trimmomatic --\n " $(which $PATH_TRIMMOMATIC/trimmomatic.jar)
[ ! -f $PATH_TRIMMOMATIC/trimmomatic.jar ] && echo "[ERROR] no trimmomatic detected" && exit 1


NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "parameters"

# get basename of f
n=${f##*/}

#is ziped ?
CAT="cat"
if [[ ${f##*.} == "gz" ]]; 
    then CAT="zcat"; 
elif [[ ${f##*.} == "bz2" ]]; 
    then CAT="bzcat"; 
fi

#is paired ?
if [ "$f" != "${f/%$READONE.$FASTQ/$READTWO.$FASTQ}" ] && [ -e ${f/%$READONE.$FASTQ/$READTWO.$FASTQ} ]; then
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
if [ -z $FASTQ_ENCODING ]; then
    echo "[NOTE] Detect fastq Phred encoding"
    FASTQ_ENCODING=$($CAT $f |  awk 'NR % 4 ==0' | python $NGSANE_BASE/tools/GuessFastqEncoding.py |  tail -n 1)
    echo "[NOTE] $FASTQ_ENCODING fastq format detected"
fi
if [[ "$FASTQ_ENCODING" == *Phred33* ]]; then
    FASTQ_PHRED=" -phred33"    
elif [[ "$FASTQ_ENCODING" == *Phred64* ]]; then
    FASTQ_PHRED=" -phred64"
else
    echo "[ERROR] cannot detect/don't understand fastq format: $FASTQ_ENCODING" && exit 1
fi
echo "[NOTE] $FASTQ_ENCODING fastq format detected"

FASTQDIR=$(basename $(dirname $f))
o=${f/$FASTQDIR/$FASTQDIR"_"$TASK_TRIMMOMATIC}
FASTQDIRTRIM=$(dirname $o)

echo $FASTQDIRTRIM
if [ ! -d $FASTQDIRTRIM ]; then mkdir -p $FASTQDIRTRIM; fi
echo $f "->" $o
if [ "$PAIRED" = "1" ]; then echo ${f/%$READONE.$FASTQ/$READTWO.$FASTQ} "->" ${o/%$READON.$FASTQE/$READTWO.$FASTQ} ; fi


NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "recall files from tape"

if [ -n "$DMGET" ]; then
    dmget -a ${f/$READONE/"*"}
    dmget -a ${o/$READONE/"*"}
fi

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "trim"    

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    if [ "$PAIRED" = "1" ]; then
        RUN_COMMAND="java -jar $PATH_TRIMMOMATIC/trimmomatic.jar PE $FASTQ_PHRED -threads $CPU_TRIMMOMATIC $f ${f/%$READONE.$FASTQ/$READTWO.$FASTQ} $o ${o}_unpaired ${o/%$READONE.$FASTQ/$READTWO.$FASTQ} ${o/%$READONE.$FASTQ/$READTWO.$FASTQ}_unpaired $TRIMMOMATICSTEPS > ${o/%$READONE.$FASTQ/}.log 2>&1"

    else
        RUN_COMMAND="java -jar $PATH_TRIMMOMATIC/trimmomatic.jar SE $FASTQ_PHRED -threads $CPU_TRIMMOMATIC $f $o $TRIMMOMATICSTEPS > ${o/%$READONE.$FASTQ/}.log 2>&1"
    fi
    echo $RUN_COMMAND && eval $RUN_COMMAND

    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $o

fi

################################################################################
[ -e $FASTQDIRTRIM/$n.dummy ] && rm $FASTQDIRTRIM/$n.dummy
echo ">>>>> readtrimming with TRIMMOMATIC - FINISHED"
echo ">>>>> enddate "`date`