#!/bin/bash -e

# Script takes a fastq and converts it into sanger format <Run>_sanger/*.$FASTQ[.gz]
#
# author: Fabian Buske
# date: June. 2013

# messages to look out for -- relevant for the QC.sh script:
# QCVARIABLES,

echo ">>>>> fastq conversion to sanger "
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> fastq2sanger.sh $*"

while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of NGSANE
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
module list
echo "PATH=$PATH"

#SAMPLENAME
# get basename of f
n=${f##*/}

#is paired ?
if [ -e ${f/$READONE/$READTWO} ]; then
    echo "[NOTE] PAIRED library"
    PAIRED="1"
else
    echo "[NOTE] SINGLE library"
    PAIRED="0"
fi

FASTQDIR=$(basename $(dirname $f))
o=${f/$FASTQDIR/$FASTQDIR"_"$TASKFASTQ2SANGER}
FASTQDIRSANGER=$(dirname $o)

if [ -n "$DMGET" ]; then
    dmget -a ${f/$READONE/"*"}
fi

echo $FASTQDIRSANGER
if [ ! -d $FASTQDIRSANGER ]; then mkdir -p $FASTQDIRSANGER; fi
echo $f "->" $o
if [ "$PAIRED" = "1" ]; then echo ${f/$READONE/$READTWO} "->" ${o/$READONE/$READTWO} ; fi

#is ziped ?                                                                                                       
ZCAT="zcat"
COMPRESS=$GZIP
if [[ ${f##*.} != "gz" ]]; then 
    ZCAT="cat"; 
    COMPRESS="echo"
fi

echo "********** convert"

RUN_COMMAND="$ZCAT $f | perl ${NGSANE_BASE}/tools/fq_all2std.pl $FASTQ2SANGER_SOURCEFORMAT | $COMPRESS > $o"
echo $RUN_COMMAND
eval $RUN_COMMAND

if [ "$PAIRED" = "1" ]; then
    RUN_COMMAND="$ZCAT ${f/$READONE/$READTWO} | perl ${NGSANE_BASE}/tools/fq_all2std.pl $FASTQ2SANGER_SOURCEFORMAT | $COMPRESS > ${o/$READONE/$READTWO}"
    echo $RUN_COMMAND
    eval $RUN_COMMAND
fi

echo ">>>>> fastq conversion - FINISHED"
echo ">>>>> enddate "`date`

