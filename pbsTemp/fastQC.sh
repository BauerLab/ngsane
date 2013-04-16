#!/bin/bash

while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the HiSeqInf 
        -h | --help )           usage ;;
        * )                     echo "don't understand "$1
    esac
    shift
done



. $CONFIG
. $HISEQINF/pbsTemp/header.sh
. $CONFIG

echo "********** programs"
module load $MODULE_FASTQC; export PATH=$PATH_FASTQC:$PATH;
module list
java -Xmx200M -version
fastqc -version

echo "********** get input"
for d in ${DIR[@]}; do
    FILES=$FILES" "$( ls $OUT/fastq/$d/*$FASTQ )
done

echo $FILES
if [ ! -d $OUT/runStats ]; then mkdir -p $OUT/runStats; fi
if [ -d $OUT/runStats/$TASKFASTQC ]; then rm -rf $OUT/runStats/$TASKFASTQC/; fi
mkdir -p $OUT/runStats/$TASKFASTQC

CPUS=`echo $FILES | wc -w`
if [ "$CPUS" -gt "$CPU_FASTQC" ]; then echo "reduce to $CPU_FASTQC CPUs"; CPUS=$CPU_FASTQC; fi

echo "********** run fastqc"
fastqc --nogroup -t $CPUS --outdir $OUT/runStats/$TASKFASTQC `echo $FILES`

chmod a+rx $OUT/runStats/fastQC
chmod a+rx $OUT/runStats/fastQC/*fastqc/Images
chmod a+rx $OUT/runStats/fastQC/*fastqc/Icons
chmod a+rx $OUT/runStats/fastQC/*fastqc
chmod a+r $OUT/runStats/fastQC/*fastqc/*
chmod a+r $OUT/runStats/fastQC/*fastqc/*/*
