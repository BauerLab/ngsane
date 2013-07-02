#!/bin/bash

# Script to trim adapters using CUTADAPT
# It takes a <Run>/*.$FASTQ[.gz] file and gets the file containing the contaminats
# via config and writes out <Run>_trim/*.$FASTQ[.gz]
# contaminants need to be specified with -a, -b or -g followd by the sequence
# -a AAGGAEE
# author: Denis Bauer
# date: April. 2013

# messages to look out for -- relevant for the QC.sh script:
# QCVARIABLES,

# TODO: for paired end reads the pairs need to be cleaned up (removed)
# with PICARD...

echo ">>>>> readtrimming with CUTADAPT "
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> cutadapt.sh $*"

while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the NGSANE
        -f | --file )           shift; f=$1 ;; # fastq file
        -h | --help )           usage ;;
        * )                     echo "don't understand "$1
    esac
    shift
done

. $CONFIG
. ${NGSANE_BASE}/conf/header.sh
. $CONFIG

#JAVAPARAMS="-Xmx"$(expr $MEMORY_CUTADAPT - 1 )"G"
#echo "JAVAPARAMS "$JAVAPARAMS

echo "********** programs"
for MODULE in $MODULE_CUTADAPT; do module load $MODULE; done  # save way to load modules that itself load other modules
export PATH=$PATH_CUTADAPT:$PATH;
module list
echo "PATH=$PATH"
echo -e "--cutadapt  --\n" $(cutadapt --version 2>&1)
[ -z "$(which cutadapt)" ] && echo "[ERROR] no cutadapt detected" && exit 1

#SAMPLENAME
# get basename of f
n=${f##*/}

if [ -e ${f/$READONE/$READTWO} ] ; then
    echo "[NOTE] PAIRED library"
    PAIRED="1"
else
    echo "[NOTE] SINGLE library"
    PAIRED="0"
fi

FASTQDIR=$(basename $(dirname $f))
o=${f/$FASTQDIR/$FASTQDIR"_"$TASKCUTADAPT}
FASTQDIRTRIM=$(dirname $o)

if [ -n "$DMGET" ]; then
    dmget -a ${f/$READONE/"*"}
fi

echo $FASTQDIRTRIM
if [ ! -d $FASTQDIRTRIM ]; then mkdir -p $FASTQDIRTRIM; fi
echo $f "->" $o
if [ "$PAIRED" = 1 ]; then echo ${f/$READONE/$READTWO} "->" ${o/$READONE/$READTWO} ; fi

echo "********** get contaminators"
echo "[NOTE] contaminants: "$CONTAMINANTS
if [ ! -n "$CONTAMINANTS" ];then echo "[WARN] need variable CONTAMINANTS defined in $CONFIG"; fi

CONTAM=$(cat $CONTAMINANTS | tr '\n' ' ')
echo "[NOTE] $CONTAM"

echo "********** trim"
RUN_COMMAND="cutadapt $CUTADAPTADDPARAM $CONTAM $f -o $o > $o.stats"
echo $RUN_COMMAND
eval $RUN_COMMAND
cat $o.stats

if [ "$PAIRED" = 1 ]; then
    echo "********** paired end"
    RUN_COMMAND="cutadapt $CUTADAPTADDPARAM $CONTAM ${f/$READONE/$READTWO} -o ${o/$READONE/$READTWO} > ${o/$READONE/$READTWO}.stats"
    echo $RUN_COMMAND
    eval $RUN_COMMAND
    cat ${o/$READONE/$READTWO}.stats
    #TODO: clean up unmached pairs
fi

echo "********** zip"
$GZIP -t $o 2>/dev/null
if [[ $? -ne 0 ]]; then
    $GZIP -f $o
    if [ "$PAIRED" = "1" ]; then
        $GZIP -f ${o/$READONE/$READTWO}
    fi
fi

# count remaining reads
echo "=== Remaining reads ===" >> $FASTQDIRTRIM/${n}.stats
echo "remaining reads "$(zcat $FASTQDIRTRIM/$n | wc -l | gawk '{print int($1/4)}') >> $FASTQDIRTRIM/${n}.stats
if [ "$PAIRED" = "1" ]; then
    echo "=== Remaining reads ===" >> $FASTQDIRTRIM/${n/$READONE/$READTWO}.stats
    echo "remaining reads "$(zcat $FASTQDIRTRIM/${n/$READONE/$READTWO} | wc -l | gawk '{print int($1/4)}') >> $FASTQDIRTRIM/${n/$READONE/$READTWO}.stats
fi

echo ">>>>> readtrimming with CUTADAPT - FINISHED"
echo ">>>>> enddate "`date`

