#!/bin/bash -e
# author: Denis C. Bauer
# date: Feb.2011

echo ">>>>> fastQC"
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"

while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the NGSANE 
        -h | --help )           usage ;;
        --recover-from )        shift; RECOVERFROM=$1 ;; # attempt to recover from log file                                                  
        * )                     echo "don't understand "$1
    esac
    shift
done

. $CONFIG
. ${NGSANE_BASE}/conf/header.sh
. $CONFIG

################################################################################
CHECKPOINT="programs"

for MODULE in $MODULE_FASTQC; do module load $MODULE; done  # save way to load modules that itself load other modules
export PATH=$PATH_FASTQC:$PATH;
module list
echo "PATH=$PATH"
echo -e "--JAVA     --\n" $(java -version 2>&1)
[ -z "$(which java)" ] && echo "[ERROR] no java detected" && exit 1
echo -e "--FASTqc   --\n" $(fastqc -version 2>&1)
[ -z "$(which fastqc)" ] && echo "[ERROR] no fastqc detected" && exit 1

echo "[NOTE] set java parameters"
JAVAPARAMS="-Xmx"$(python -c "print int($MEMORY_FASTQC*0.8)")"g -Djava.io.tmpdir="$TMP" -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1" 
unset _JAVA_OPTIONS
echo "JAVAPARAMS "$JAVAPARAMS

echo -e "\n********* $CHECKPOINT"
################################################################################
CHECKPOINT="get input"

for d in ${DIR[@]}; do
    FILES=$FILES" "$( ls $OUT/fastq/$d/*$FASTQ )
done
echo $FILES

mkdir -p $OUT/runStats

if [ -d $OUT/runStats/$TASKFASTQC ]; then rm -r $OUT/runStats/$TASKFASTQC/; fi
mkdir -p $OUT/runStats/$TASKFASTQC

CPUS=$(echo $FILES | wc -w)
if [ "$CPUS" -gt "$CPU_FASTQC" ]; then echo "[NOTE] reduce to $CPU_FASTQC CPUs"; CPUS=$CPU_FASTQC; fi

echo -e "\n********* $CHECKPOINT"
################################################################################
CHECKPOINT="run fastqc"

RUN_COMMAND="fastqc --nogroup -t $CPUS --outdir $OUT/runStats/$TASKFASTQC $FILES"
echo $RUN_COMMAND && eval $RUN_COMMAND

chmod a+rx $OUT/runStats/
chmod a+rx $OUT/runStats/fastQC
chmod a+rx $OUT/runStats/fastQC/*fastqc/Images
chmod a+rx $OUT/runStats/fastQC/*fastqc/Icons
chmod a+rx $OUT/runStats/fastQC/*fastqc
chmod a+r $OUT/runStats/fastQC/*fastqc/*
chmod a+r $OUT/runStats/fastQC/*fastqc/*/*

echo -e "\n********* $CHECKPOINT"
################################################################################
echo ">>>>> fastQC - FINISHED"
echo ">>>>> enddate "`date`

