#!/bin/bash

# Script to trim adapters using TRIMGALORE (tapping into CUTADAPT)
# It takes a <Run>/*.$FASTQ[.gz] file and gets the file containing the contaminats
# via config and writes out <Run>_trimgalore/*.$FASTQ[.gz]
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
echo ">>>>> $(basename $0) $*"

while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of NGSANE
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

for MODULE in $MODULE_TRIMGALORE; do module load $MODULE; done  # save way to load modules that itself load other modules
export PATH=$PATH_TRIMGALORE:$PATH;
module list
echo "PATH=$PATH"

echo -e "--trim galore --\n "$(trim_galore --version  | grep version  | tr -d ' ')
[ -z "$(which trim_galore)" ] && echo "[ERROR] no trim_galore detected" && exit 1

echo -n "********* $CHECKPOINT"
################################################################################
CHECKPOINT="parameters"

# get basename of f (samplename)
n=${f##*/}

#is paired ?
if [ "$f" != "${f/$READONE/$READTWO}" ] && [ -e ${f/$READONE/$READTWO} ]; then
    echo "[NOTE] PAIRED library"
    PAIRED="1"
else
    echo "[NOTE] SINGLE library"
    PAIRED="0"
fi

FASTQDIR=$(basename $(dirname $f))
o=${f/$FASTQDIR/$FASTQDIR"_"$TASKTRIMGALORE}
FASTQDIRTRIM=$(dirname $o)

echo $FASTQDIRTRIM
if [ ! -d $FASTQDIRTRIM ]; then mkdir -p $FASTQDIRTRIM; fi
echo $f "->" $o
if [ "$PAIRED" = "1" ]; then echo ${f/$READONE/$READTWO} "->" ${o/$READONE/$READTWO} ; fi

echo "********** get contaminators"

if [ ! -n "$TRIMGALORE_ADAPTER1" ] && [ ! -n "$TRIMGALORE_ADAPTER2" ];then echo "TRIMGALORE_ADAPTER1 and 2 not defined in $CONFIG, default to 'AGATCGGAAGAGC'"; fi
CONTAM=""
if [ -n "$TRIMGALORE_ADAPTER1" ]; then
	CONTAM="$CONTAM --adapter $TRIMGALORE_ADAPTER1"
fi
if [ "$PAIRED" = "1" ] && [ -n "$TRIMGALORE_ADAPTER2" ]; then
	CONTAM="$CONTAM --adapter2 $TRIMGALORE_ADAPTER2"
fi
echo $CONTAM

echo -n "********* $CHECKPOINT"
################################################################################
CHECKPOINT="recall files from tape"

if [ -n "$DMGET" ]; then
    dmget -a ${f/$READONE/"*"}
fi

echo -n "********* $CHECKPOINT"
################################################################################
CHECKPOINT="trim"    

if [[ -n "$RECOVERFROM" ]] && [[ $(grep "********* $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo -n "::::::::: passed $CHECKPOINT"
else 
    
    # Paired read
    if [ "$PAIRED" = "1" ]
    then
        trim_galore $TRIMGALOREADDPARAM $CONTAM --paired --output_dir $FASTQDIRTRIM $f ${f/$READONE/$READTWO}
        mv $FASTQDIRTRIM/${n/$READONE.$FASTQ/$READONE"_val_1".fq.gz} $FASTQDIRTRIM/$n
        mv $FASTQDIRTRIM/${n/$READONE.$FASTQ/$READTWO"_val_2".fq.gz} $FASTQDIRTRIM/${n/$READONE/$READTWO}
    else
        trim_galore $TRIMGALOREADDPARAM $CONTAM --output_dir $FASTQDIRTRIM $f
        mv $FASTQDIRTRIM/${n/$READONE.$FASTQ/$READONE"_trimmed".fq.gz} $FASTQDIRTRIM/$n
    fi

    # mark checkpoint
    [ -f $FASTQDIRTRIM/$n ] && echo -n "********* $CHECKPOINT"
fi

################################################################################
CHECKPOINT="zip"    

if [[ -n "$RECOVERFROM" ]] && [[ $(grep "********* $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo -n "::::::::: passed $CHECKPOINT"
else 
    $GZIP -t $FASTQDIRTRIM/$n 2>/dev/null
    if [[ $? -ne 0 ]]; then
        $GZIP -f $FASTQDIRTRIM/$n
        mv $FASTQDIRTRIM/$n.gz $FASTQDIRTRIM/$n
        if [ "$PAIRED" = "1" ]; then
            $GZIP -f $FASTQDIRTRIM/${n/$READONE/$READTWO}
            mv $FASTQDIRTRIM/${n/$READONE/$READTWO}.gz $FASTQDIRTRIM/${n/$READONE/$READTWO}
        fi
    fi
    # mark checkpoint
    echo -n "********* $CHECKPOINT"
fi

################################################################################
CHECKPOINT="count remaining reads"    

echo "=== Remaining reads ===" >> $FASTQDIRTRIM/${n}_trimming_report.txt
echo "remaining reads "$(zcat $FASTQDIRTRIM/$n | wc -l | gawk '{print int($1/4)}') >> $FASTQDIRTRIM/${n}_trimming_report.txt
if [ "$PAIRED" = "1" ]; then
    echo "=== Remaining reads ===" >> $FASTQDIRTRIM/${n/$READONE/$READTWO}_trimming_report.txt
    echo "remaining reads "$(zcat $FASTQDIRTRIM/${n/$READONE/$READTWO} | wc -l | gawk '{print int($1/4)}') >> $FASTQDIRTRIM/${n/$READONE/$READTWO}_trimming_report.txt
fi

echo "********* $CHECKPOINT"
################################################################################
echo ">>>>> readtrimming with TRIMGALORE - FINISHED"
echo ">>>>> enddate "`date`

