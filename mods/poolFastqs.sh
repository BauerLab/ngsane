#!/bin/bash -e
# author: Fabian Buske
# date: September 2013
# pool/merge fastq files from within a experiment folder subject to a pattern 

echo ">>>>> pool fastq datasets"
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"

while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the NGSANE 
        -h | --help )           usage ;;
        --recover-from )        shift; NGSANE_RECOVERFROM=$1 ;; # attempt to recover from log file                                                  
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
hash module 2>/dev/null && for MODULE in $MODULE_POOLFASTQS; do module load $MODULE; done && module list 

export PATH=$PATH_POOLFASTQS:$PATH;
echo "PATH=$PATH"

echo -e "--NGSANE       --\n" $(trigger.sh -v 2>&1)
echo -e "--gnu parallel --\n "$(parallel --gnu --version 2>&1 | tee | head -n 1)
[ -z "$(which parallel 2> /dev/null)" ] && echo "[WARN] no gnu parallel detected, processing in serial"


NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "pair input"

COMMAND=$OUT/tmp/pool_fastqs`date +%Y%m%d`.commands

cat /dev/null > $COMMAND
for d in ${DIR[@]}; do
    grep '^POOL_FASTQ ' $CONFIG | cut -d' ' -f 2- > $OUT/$INPUT_POOLFASTQS/$d/pools.tmp

    while read -r -a POOL; do
        PATTERN=$(echo "${POOL[@]:1}" | sed -e 's/ /|/g')

        OUTFASTQ=$OUT/$INPUT_POOLFASTQS/$d/${POOL[0]}"."$FASTQ
        [ -f $OUTFASTQ ] && rm $OUTFASTQ
        mkdir -p $OUT/$d/$TASK_POOLFASTQS

        INFASTQ=""
        for i in $(ls $OUT/$INPUT_POOLFASTQS/$d/*$FASTQ | egrep "$PATTERN"); do 
            INFASTQ="$INFASTQ $i"
        done

        if [ -n "$INFASTQ" ]; then   
            echo -ne "cat  $INFASTQ > $OUTFASTQ" >> $COMMAND
            if [[ -n "$DELETEORIGINALFASTQS" ]]; then
                echo -ne ";rm $INFASTQ" >> $COMMAND            
            fi
            echo ";" >> $COMMAND
        fi
    done < $OUT/$INPUT_POOLFASTQS/$d/pools.tmp
    rm $OUT/$INPUT_POOLFASTQS/$d/pools.tmp
    
done
echo "COMMAND:" $(cat $COMMAND)

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "pool data"

if hash parallel ; then

    echo "[NOTE] parallel processing"
    cat $COMMAND | parallel --joblog $TMP/$TASK_POOLFASTQS.log --gnu -j $CPU_POOLFASTQS "eval {}" > /dev/null 2>&1

else
    # serial processing
    echo "[NOTE] serial processing"
    while read line; do
        echo $line
        eval $line
    done < $COMMAND
fi

rm $COMMAND

NGSANE_CHECKPOINT_CHECK
################################################################################
echo ">>>>> pool fastq datasets - FINISHED"
echo ">>>>> enddate "`date`

