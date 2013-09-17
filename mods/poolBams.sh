#!/bin/bash -e
# author: Fabian Buske
# date: September 2013

echo ">>>>> pool bam datasets"
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

for MODULE in $MODULE_POOLBAMS; do module load $MODULE; done  # save way to load modules that itself load other modules
export PATH=$PATH_POOLBAMS:$PATH;
module list
echo "PATH=$PATH"
PATH_IGVTOOLS=$(dirname $(which igvtools.jar))
PATH_PICARD=$(dirname $(which MergeSamFiles.jar))

echo "[NOTE] set java parameters"
JAVAPARAMS="-Xmx"$(python -c "print int($MEMORY_BOWTIE*0.8)")"g -Djava.io.tmpdir="$TMP" -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1" 
unset _JAVA_OPTIONS
echo "JAVAPARAMS "$JAVAPARAMS

echo -e "--NGSANE       --\n" $(trigger.sh -v 2>&1)
echo -e "--JAVA        --\n" $(java -Xmx200m -version 2>&1)
[ -z "$(which java)" ] && echo "[ERROR] no java detected" && exit 1
echo -e "--PICARD      --\n "$(java $JAVAPARAMS -jar $PATH_PICARD/MergeSamFiles.jar --version 2>&1)
[ ! -f $PATH_PICARD/MergeSamFiles.jar ] && echo "[ERROR] no picard detected" && exit 1
echo -e "--igvtools    --\n "$(java $JAVAPARAMS -jar $PATH_IGVTOOLS/igvtools.jar version 2>&1)
[ ! -f $PATH_IGVTOOLS/igvtools.jar ] && echo "[ERROR] no igvtools detected" && exit 1
echo -e "--samstat     --\n "$(samstat -h | head -n 2 | tail -n1)
[ -z "$(which samstat)" ] && echo "[ERROR] no samstat detected" && exit 1
echo -e "--gnu parallel --\n "$(parallel --version 2>&1 | tee | head -n 1)
[ -z "$(which parallel 2> /dev/null)" ] && echo "[WARN] no gnu parallel detected, processing in serial"


echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="pair input"

if [ -z "$PATTERN" ]; then
    echo "[ERROR] no pattern defined"
    exit 1
fi
if [ -z "$FASTA" ]; then
    echo "[ERROR] no reference defined: FASTA"
    exit 1
fi
if [ -z "$REPLACEWITH" ]; then
    echo "[ERROR] no replacement defined, defaulting to _pooled"
    REPLACEWITH="_pooled"
fi

COMMAND=$OUT/tmp/pool_bams`date +%Y%m%d`.commands

cat /dev/null > $COMMAND

for d in ${DIR[@]}; do
    for i in $(ls $OUT/$d/$INPUT_POOLBAMS/*$ASD.bam | grep -P "$PATTERN"); do 
        echo "$i $i" |  sed -e "s|$PATTERN||"; 
    done | sort -k1,1 > $OUT/$d/$INPUT_POOLBAMS/pattern.tmp

    OLDIFS=$IFS
    for POOL in $(cat $OUT/$d/$INPUT_POOLBAMS/pattern.tmp | cut -d' ' -f 1 | sort -u); do 
        OUTBAM=$(grep "$POOL" $OUT/$d/$INPUT_POOLBAMS/pattern.tmp | cut -d' ' -f 2 | head -n 1 | sed -e "s|$PATTERN|$REPLACEWITH|" ) 
        INBAMS=$(grep "$POOL" $OUT/$d/$INPUT_POOLBAMS/pattern.tmp | awk '{print "INPUT="$2}' | tr '\n' ' ')
        [ -f $OUTBAM ] && rm $OUTBAM
        echo -ne "java $JAVAPARAMS -jar $PATH_PICARD/MergeSamFiles.jar CREATE_INDEX=true VERBOSITY=ERROR VALIDATION_STRINGENCY=LENIENT OUTPUT=$OUTBAM $INBAMS; samstat $OUTBAM; java $JAVAPARAMS -jar $PATH_IGVTOOLS/igvtools.jar count $OUTBAM $OUTBAM.cov.tdf ${FASTA%%.*}.genome " >> $COMMAND
        if [ "$DELETEORIGINALBAMS" = "true" ]; then
            for j in $(grep "$POOL" $OUT/$d/$INPUT_POOLBAMS/pattern.tmp | cut -d' ' -f 2 ); do
                echo -ne "; rm $j*" >> $COMMAND
            done
        fi
        echo ";" >> $COMMAND
    done
done

CPUS=$(wc -l $COMMAND | cut -d' ' -f 1)
echo $CPUS
if [[ "$CPUS" -gt "$CPU_POOLBAMS" ]]; then 
    echo "[NOTE] reduce to $CPU_POOLBAMS CPUs"; CPUS=$CPU_POOLBAMS; 
fi

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="pool data"

exit
if hash parallel ; then
    echo "[NOTE] parallel processing"
    
    cat $COMMAND | parallel --eta -j $CPUS "eval {}"

else
    # serial processing
    echo "[NOTE] serial processing"
    while read line; do
        echo $line
        eval $line
    done < $COMMAND
fi

rm $COMMAND

echo -e "\n********* $CHECKPOINT\n"
################################################################################
echo ">>>>> pool bam datasets - FINISHED"
echo ">>>>> enddate "`date`

