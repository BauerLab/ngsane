#!/bin/bash
#PBS -l walltime=10:00:00
#PBS -l vmem=20gb
#PBS -l nodes=2:ppn=8
#PBS -j oe

cd ${PBS_O_WORKDIR:-.}

module load jdk

#CONFIG=$1
MAX=16

. $CONFIG
. $DATASTORE/SeqAna/apps/prod/seqaninf/pbsTemp/header.sh
. $CONFIG

for d in ${DIR[@]}; do
    FILES=$FILES" "$( ls $OUT/fastq/$d/*$FASTQ )
done

echo $FILES

CPUS=`echo $FILES | wc -w`
if [ "$CPUS" -gt "$MAX" ]; then echo "reduce to $MAX CPUs"; CPUS=$MAX; fi


$FASTQC --nogroup -t $CPUS --outdir $OUT/runStats/$TASKFASTQC `echo $FILES`
