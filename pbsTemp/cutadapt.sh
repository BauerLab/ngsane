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
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the HiSeqInf
        -f | --file )           shift; f=$1 ;; # fastq file
        -o | --outdir )         shift; OUTDIR=$1 ;; # output dir
        -h | --help )           usage ;;
        * )                     echo "don't understand "$1
    esac
    shift
done

. $CONFIG
. $HISEQINF/pbsTemp/header.sh
. $CONFIG

#JAVAPARAMS="-Xmx"$(expr $MEMORY_CUTADAPT - 1 )"G"
#echo "JAVAPARAMS "$JAVAPARAMS

echo "********** programs"
module load $MODULE_CUTADAPT; export PATH=$PATH_CUTADAPT:$PATH;
module list
cutadapt --version

FASTQDIR=$(basename $(dirname $f))
o=${f/$FASTQDIR/$FASTQDIR"_trim"}
FASTQDIRTRIM=$(dirname $o)


if [ -e ${f/$READONE/$READTWO} ] ; then
    echo "PAIRED"
    PAIRED="1"
    o2=${o/$READONE/$READTWO}
else
    echo "SINGLE"
    PAIRED="0"
fi

if [ -n "$DMGET" ]; then
    dmget -a ${f/$READONE/"*"}
fi

#echo $FASTQDIRTRIM
if [ ! -d $FASTQDIRTRIM ]; then mkdir $FASTQDIRTRIM; fi
echo $f "->" $o
if [ "$PAIRED" = 1 ]; then echo ${f/$READONE/$READTWO} "->" $o2 ; fi
echo "contaminants: "$CONTAMINANTS
if [ ! -n "$CONTAMINANTS" ];then echo "need variable CONTAMINANTS defined in $CONFIG"; fi

echo "********** get contaminators"
CONTAM=$(cat $CONTAMINANTS | tr '\n' ' ')
echo $CONTAM

echo "********** trim"
cutadapt $CONTAM $f -o $o > $o.stats

if [ "$PAIRED" = 1 ]; then
    echo "********** paired end"
    cutadapt $CONTAM ${f/$READONE/$READTWO} -o $o2 > $o2.stats
    #TODO: clean up unmached pairs
fi

echo "********** zip"
for i in $o $o2 ; do
    $GZIP -f $i
done

echo ">>>>> readtrimming with CUTADAPT - FINISHED"
echo ">>>>> enddate "`date`

