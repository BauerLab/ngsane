#!/bin/bash

# Script to trim adapters using TRIMGALORE (tapping into CUTADAPT)
# It takes a <Run>/*.$FASTQ[.gz] file and gets the file containing the contaminats
# via config and writes out <Run>_trim/*.$FASTQ[.gz]
#
# author: Fabian Buske
# date: April. 2013

# messages to look out for -- relevant for the QC.sh script:
# QCVARIABLES,

echo ">>>>> readtrimming with TRIMGALORE "
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> trimgalore.sh $*"

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
module load $MODULE_TRIMGALORE; export PATH=$PATH_TRIMGALORE:$PATH;
module list
trim_galore --version  | grep version  | tr -d ' '

#is paired ?
if [ -e ${f/$READONE/$READTWO} ]; then
    PAIRED="1"
else
    PAIRED="0"
fi

FASTQDIR=$(basename $(dirname $f))
o=${f/$FASTQDIR/$FASTQDIR"_trim"}
FASTQDIRTRIM=$(dirname $o)

echo $FASTQDIRTRIM
if [ ! -d $FASTQDIRTRIM ]; then mkdir -p $FASTQDIRTRIM; fi
echo $f "->" $o
if [ "$PAIRED" = 1 ]; then echo ${f/$READONE/$READTWO} "->" ${o/$READONE/$READTWO} ; fi

echo "********** get contaminators"

if [ ! -n "$TRIMGALORE_ADAPTER1" ] && [ ! -n "$TRIMGALORE_ADAPTER2" ];then echo "TRIMGALORE_ADAPTER1 and 2 not defined in $CONFIG, default to 'AGATCGGAAGAGC'"; fi
CONTAM=""
if [ -n "$TRIMGALORE_ADAPTER1" ]; then
	CONTAM="$CONTAM --adapter $TRIMGALORE_ADAPTER1"
fi
if [ "$PAIRED" = 1 ] && [ -n "$TRIMGALORE_ADAPTER2" ]; then
	CONTAM="$CONTAM --adapter2 $TRIMGALORE_ADAPTER2"
fi
echo $CONTAM

echo "********** trim"
# Paired read
if [ "$PAIRED" = 1 ]
then
	trim_galore $CONTAM --paired --output_dir $FASTQDIRTRIM $f ${f/$READONE/$READTWO}
else
	trim_galore $CONTAM --output_dir $FASTQDIRTRIM $f
fi

echo ">>>>> readtrimming with TRIMGALORE - FINISHED"
echo ">>>>> enddate "`date`

