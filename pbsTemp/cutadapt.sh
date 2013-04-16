#!/bin/bash

# Script to trim adapters using CUTADAPT
# It takes a <Run>/*.$FASTQ[.gz] file and gets the file containing the contaminats
# via config and writes out <Run>_trim/*.$FASTQ[.gz]
# author: Denis Bauer
# date: April. 2013

# messages to look out for -- relevant for the QC.sh script:
# QCVARIABLES,


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
FASTQDIRTRIM=$(basename $o)


if [ ! -d $FASTQDIRTRIM ]; then mkdir $FASTQDIRTRIM; fi
echo $f "->" $o
echo "contaminants: "$CONTAMINANTS
if [ ! -n "$CONTAMINANTS" ];then echo "need variable CONTAMINANTS defined in $CONFIG"; fi


echo "********** get contaminators"
CONTAM=""
for i in $(cat $CONTAMINANTS ); do
	CONTAM=$CONTAM+"-a $i "
done

echo $CONTAM

echo "********** trim"
$CUTADAPT $CONTAM $f -o $o > $o.stats

echo "********** zip"
gzip ${o}

echo ">>>>> readtrimming with CUTADAPT - FINISHED"
echo ">>>>> enddate "`date`

