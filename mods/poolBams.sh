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

# save way to load modules that itself loads other modules
hash module 2>/dev/null && for MODULE in $MODULE_POOLBAMS; do module load $MODULE; done && module list 

export PATH=$PATH_POOLBAMS:$PATH;
echo "PATH=$PATH"
[ -z "$PATH_PICARD" ] && PATH_PICARD=$(dirname $(which MergeSamFiles.jar))
 
echo "[NOTE] set java parameters"
JAVAPARAMS="-Xmx"$(python -c "print int($MEMORY_POOLBAMS*0.8)")"g -Djava.io.tmpdir="$TMP" -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1" 
unset _JAVA_OPTIONS
echo "JAVAPARAMS "$JAVAPARAMS

echo -e "--NGSANE       --\n" $(trigger.sh -v 2>&1)
echo -e "--JAVA        --\n" $(java -Xmx200m -version 2>&1)
[ -z "$(which java)" ] && echo "[ERROR] no java detected" && exit 1
echo -e "--PICARD      --\n "$(java $JAVAPARAMS -jar $PATH_PICARD/MergeSamFiles.jar --version 2>&1)
[ ! -f $PATH_PICARD/MergeSamFiles.jar ] && echo "[ERROR] no picard detected" && exit 1
echo -e "--samtools    --\n "$(samtools 2>&1 | head -n 3 | tail -n-2)
[ -z "$(which samtools)" ] && echo "[ERROR] no samtools detected" && exit 1
echo -e "--samstat     --\n "$(samstat -h | head -n 2 | tail -n1)
[ -z "$(which samstat)" ] && echo "[ERROR] no samstat detected" && exit 1
echo -e "--gnu parallel --\n "$(parallel --gnu --version 2>&1 | tee | head -n 1)
[ -z "$(which parallel 2> /dev/null)" ] && echo "[WARN] no gnu parallel detected, processing in serial"


echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="pair input"

if [ -z "$FASTA" ]; then
    echo "[ERROR] no reference defined: FASTA"
    exit 1
fi

COMMAND=$OUT/tmp/pool_bams`date +%Y%m%d`.commands

cat /dev/null > $COMMAND
for d in ${DIR[@]}; do
    grep '^POOL_BAM ' $CONFIG | cut -d' ' -f 2- > $OUT/$d/$INPUT_POOLBAMS/pools.tmp

    while read -r -a POOL; do
        PATTERN=$(echo "${POOL[@]:1}" | sed -e 's/ /|/g')

        OUTBAM="$OUT/$d/$INPUT_POOLBAMS/$(echo $POOL | cut -d' ' -f1).$ASD.bam"
        [ -f $OUTBAM ] && rm $OUTBAM

        INBAMS=""
        COMMENT=""
        for i in $(ls $OUT/$d/$INPUT_POOLBAMS/*$ASD.bam | egrep "$PATTERN"); do 
            INBAMS="$INBAMS INPUT=$i"
            COMMENT="$COMMENT ${i##*/}"
        done        
                
        if [ -n "$INBAMS" ]; then   
            echo -ne "java $JAVAPARAMS -jar $PATH_PICARD/MergeSamFiles.jar ASSUME_SORTED=true QUIET=true VERBOSITY=ERROR VALIDATION_STRINGENCY=LENIENT TMP_DIR=$TMP COMPRESSION_LEVEL=9 USE_THREADING=true OUTPUT=$OUTBAM $INBAMS COMMENT='merged:$COMMENT'; samtools index $OUTBAM; samtools flagstat $OUTBAM > $OUTBAM.stats" >> $COMMAND
            echo ";" >> $COMMAND
        fi
    done < $OUT/$d/$INPUT_POOLBAMS/pools.tmp
    
    rm $OUT/$d/$INPUT_POOLBAMS/pools.tmp
    
done

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="pool data"

if hash parallel ; then
    echo "[NOTE] parallel processing"

    cat $COMMAND | parallel --verbose --joblog $TMP/$TASK_POOLBAMS.log --gnu --eta -j $CPU_POOLBAMS "eval {}" > /dev/null 2&>1

else
    # serial processing
    echo "[NOTE] serial processing"
    while read line; do
        echo $line
        eval $line
    done < $COMMAND
fi

#rm $COMMAND

echo -e "\n********* $CHECKPOINT\n"
################################################################################
echo ">>>>> pool bam datasets - FINISHED"
echo ">>>>> enddate "`date`

