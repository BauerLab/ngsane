#!/bin/bash -e

# Script to trim adapters using CUTADAPT
# It takes a <Run>/*.$FASTQ[.gz] file and gets the file containing the contaminats
# via config and writes out <Run>_trim/*.$FASTQ[.gz]
# contaminants need to be specified with -a, -b or -g followd by the sequence
# -a AAGGAEE
# author: Denis Bauer
# date: April. 2013

# messages to look out for -- relevant for the QC.sh script:
# QCVARIABLES,
# RESULTFILENAME fastq/<DIR>"_"$TASK_CUTADAPT/<SAMPLE>$READONE.$FASTQ


# TODO: for paired end reads the pairs need to be cleaned up (removed)
# with PICARD...

echo ">>>>> readtrimming with CUTADAPT "
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"

while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the NGSANE
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
for MODULE in $MODULE_CUTADAPT; do module load $MODULE; done  # save way to load modules that itself load other modules
export PATH=$PATH_CUTADAPT:$PATH;
module list
echo "PATH=$PATH"

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--cutadapt    --\n" $(cutadapt --version 2>&1)
[ -z "$(which cutadapt)" ] && echo "[ERROR] no cutadapt detected" && exit 1

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="parameters"

# get basename of f (samplename)
n=${f##*/}

#is paired ?
if [ "$f" != "${f/$READONE/$READTWO}" ] && [ -e ${f/$READONE/$READTWO} ] ; then
    echo "[NOTE] PAIRED library"
    echo "[WARN] paired end library detected - consider using trimmomatic or trimgalore"
    PAIRED="1"
else
    echo "[NOTE] SINGLE library"
    PAIRED="0"
fi

FASTQDIR=$(basename $(dirname $f))
o=${f/$FASTQDIR/$FASTQDIR"_"$TASK_CUTADAPT}
FASTQDIRTRIM=$(dirname $o)

echo $FASTQDIRTRIM
if [ ! -d $FASTQDIRTRIM ]; then mkdir -p $FASTQDIRTRIM; fi
echo $f "->" $o
if [ "$PAIRED" = 1 ]; then echo ${f/$READONE/$READTWO} "->" ${o/$READONE/$READTWO} ; fi

echo "[NOTE] contaminants: "$CONTAMINANTS
[ ! -n "$CONTAMINANTS" ] && echo "[ERROR] need variable CONTAMINANTS defined in $CONFIG" && exit 1

CONTAM=$(cat $CONTAMINANTS | tr '\n' ' ')
echo "[NOTE] $CONTAM"

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="recall files from tape"

if [ -n "$DMGET" ]; then
    dmget -a ${f/$READONE/"*"}
fi

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="trim"    

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
    RUN_COMMAND="cutadapt $CUTADAPTADDPARAM $CONTAM $f -o $o > $o.stats"
    echo $RUN_COMMAND
    eval $RUN_COMMAND
    cat $o.stats
    
    if [ "$PAIRED" = 1 ]; then
        RUN_COMMAND="cutadapt $CUTADAPTADDPARAM $CONTAM ${f/$READONE/$READTWO} -o ${o/$READONE/$READTWO} > ${o/$READONE/$READTWO}.stats"
        echo $RUN_COMMAND
        eval $RUN_COMMAND
        cat ${o/$READONE/$READTWO}.stats
        #TODO: clean up unmached pairs
    fi
    
    # mark checkpoint
    if [ -f $o ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi

################################################################################
CHECKPOINT="zip"    

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

    $GZIP -t $o 2>/dev/null
    if [[ $? -ne 0 ]]; then
        $GZIP -f $o
        if [ "$PAIRED" = "1" ]; then
            $GZIP -f ${o/$READONE/$READTWO}
        fi
    fi

    # mark checkpoint
    echo -e "\n********* $CHECKPOINT\n"  
fi

################################################################################
CHECKPOINT="count remaining reads"    

echo "=== Remaining reads ===" >> $FASTQDIRTRIM/${n}.stats
echo "remaining reads "$(zcat $FASTQDIRTRIM/$n | wc -l | gawk '{print int($1/4)}') >> $FASTQDIRTRIM/${n}.stats
if [ "$PAIRED" = "1" ]; then
    echo "=== Remaining reads ===" >> $FASTQDIRTRIM/${n/$READONE/$READTWO}.stats
    echo "remaining reads "$(zcat $FASTQDIRTRIM/${n/$READONE/$READTWO} | wc -l | gawk '{print int($1/4)}') >> $FASTQDIRTRIM/${n/$READONE/$READTWO}.stats
fi

echo -e "\n********* $CHECKPOINT\n"
################################################################################
[ -e $FASTQDIRTRIM/${n}.dummy ] && rm $FASTQDIRTRIM/${n}.dummy
echo ">>>>> readtrimming with CUTADAPT - FINISHED"
echo ">>>>> enddate "`date`

