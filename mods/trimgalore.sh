#!/bin/bash -e

# Script to trim adapters using TRIMGALORE (tapping into CUTADAPT)
# It takes a <Run>/*.$FASTQ[.gz] file and gets the file containing the contaminats
# via config and writes out <Run>_trimgalore/*.$FASTQ[.gz]
#
# author: Fabian Buske
# date: April. 2013

# messages to look out for -- relevant for the QC.sh script:
# QCVARIABLES,
# RESULTFILENAME fastq/<DIR>"_"$TASK_TRIMGALORE/<SAMPLE>$READONE.$FASTQ

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
        --recover-from )        shift; NGSANE_RECOVERFROM=$1 ;; # attempt to recover from log file
        -h | --help )           usage ;;
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
hash module 2>/dev/null && for MODULE in $MODULE_TRIMGALORE; do module load $MODULE; done && module list 

export PATH=$PATH_TRIMGALORE:$PATH;
echo "PATH=$PATH"

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--trim galore --\n "$(trim_galore --version  | grep version  | tr -d ' ')
[ -z "$(which trim_galore)" ] && echo "[ERROR] no trim_galore detected" && exit 1
echo -e "--cutadapt    --\n" $(cutadapt --version 2>&1)
[ -z "$(which cutadapt)" ] && echo "[ERROR] no cutadapt detected" && exit 1

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "parameters"

# get basename of f (samplename)
n=${f##*/}

#is paired ?
if [ "$f" != "${f/%$READONE.$FASTQ/$READTWO.$FASTQ}" ] && [ -e ${f/%$READONE.$FASTQ/$READTWO.$FASTQ} ]; then
    echo "[NOTE] PAIRED library"
    PAIRED="1"
else
    echo "[NOTE] SINGLE library"
    PAIRED="0"
fi

FASTQDIRTRIM=$(dirname $f)"_"$TASK_TRIMGALORE
o=$FASTQDIRTRIM/$n

echo $FASTQDIRTRIM
if [ ! -d $FASTQDIRTRIM ]; then mkdir -p $FASTQDIRTRIM; fi
echo $f "->" $o
if [ "$PAIRED" = "1" ]; then echo ${f/%$READONE.$FASTQ/$READTWO.$FASTQ} "->" ${o/%$READONE.$FASTQ/$READTWO.$FASTQ} ; fi

#is ziped ?
CAT="cat"
if [[ ${f##*.} == "gz" ]]; 
    then CAT="zcat"; 
elif [[ ${f##*.} == "bz2" ]]; 
    then CAT="bzcat"; 
fi

if [ ! -n "$TRIMGALORE_ADAPTER1" ] && [ ! -n "$TRIMGALORE_ADAPTER2" ];then echo "TRIMGALORE_ADAPTER1 and 2 not defined in $CONFIG, default to 'AGATCGGAAGAGC'"; fi
CONTAM=""
if [ -n "$TRIMGALORE_ADAPTER1" ]; then
	CONTAM="$CONTAM --adapter $TRIMGALORE_ADAPTER1"
fi
if [ "$PAIRED" = "1" ] && [ -n "$TRIMGALORE_ADAPTER2" ]; then
	CONTAM="$CONTAM --adapter2 $TRIMGALORE_ADAPTER2"
fi
echo $CONTAM

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "recall files from tape"

if [ -n "$DMGET" ]; then
    dmget -a ${f/$READONE/"*"}
    dmget -a ${o/$READONE/"*"}
fi

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "trim"    

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then
    
    # Paired read
    if [ "$PAIRED" = "1" ]
    then
        trim_galore $TRIMGALOREADDPARAM $CONTAM --paired --output_dir $FASTQDIRTRIM <($CAT $f) <($CAT ${f/%$READONE.$FASTQ/$READTWO.$FASTQ})
        mv $FASTQDIRTRIM/${n/$READONE.$FASTQ/$READONE"_val_1".fq.gz} $FASTQDIRTRIM/$n
        mv $FASTQDIRTRIM/${n/$READONE.$FASTQ/$READTWO"_val_2".fq.gz} $FASTQDIRTRIM/${n/%$READONE.$FASTQ/$READTWO.$FASTQ}
    else
        trim_galore $TRIMGALOREADDPARAM $CONTAM --output_dir $FASTQDIRTRIM <($CAT $f)
        mv $FASTQDIRTRIM/${n/$READONE.$FASTQ/$READONE"_trimmed".fq.gz} $FASTQDIRTRIM/$n
    fi

    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $FASTQDIRTRIM/$n

fi

################################################################################
NGSANE_CHECKPOINT_INIT "zip"    

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    $GZIP -t $FASTQDIRTRIM/$n 2>/dev/null
    if [[ $? -ne 0 ]]; then
        $GZIP -f $FASTQDIRTRIM/$n
        mv $FASTQDIRTRIM/$n.gz $FASTQDIRTRIM/$n
        if [ "$PAIRED" = "1" ]; then
            $GZIP -f $FASTQDIRTRIM/${n/%$READONE.$FASTQ/$READTWO.$FASTQ}
            mv $FASTQDIRTRIM/${n/%$READONE.$FASTQ/$READTWO.$FASTQ}.gz $FASTQDIRTRIM/${n/%$READONE.$FASTQ/$READTWO.$FASTQ}
        fi
    fi
    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK
fi

################################################################################
NGSANE_CHECKPOINT_INIT "count remaining reads"    

echo "=== Remaining reads ===" >> $FASTQDIRTRIM/${n}_trimming_report.txt
echo "remaining reads "$(zcat $FASTQDIRTRIM/$n | wc -l | gawk '{print int($1/4)}') >> $FASTQDIRTRIM/${n}_trimming_report.txt
if [ "$PAIRED" = "1" ]; then
    echo "=== Remaining reads ===" >> $FASTQDIRTRIM/${n/%$READONE.$FASTQ/$READTWO.$FASTQ}_trimming_report.txt
    echo "remaining reads "$(zcat $FASTQDIRTRIM/${n/%$READONE.$FASTQ/$READTWO.$FASTQ} | wc -l | gawk '{print int($1/4)}') >> $FASTQDIRTRIM/${n/%$READONE.$FASTQ/$READTWO.$FASTQ}_trimming_report.txt
fi

NGSANE_CHECKPOINT_CHECK
################################################################################
[ -e $FASTQDIRTRIM/${n}.dummy ] && rm $FASTQDIRTRIM/${n}.dummy
echo ">>>>> readtrimming with TRIMGALORE - FINISHED"
echo ">>>>> enddate "`date`

