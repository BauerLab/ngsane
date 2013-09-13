#!/bin/bash

# Call variants with samtools
# author: Denis C. Bauer
# date: Feb.2013

# messages to look out for -- relevant for the QC.sh script:
# QCVARIABLES,Resource temporarily unavailable

echo ">>>>> Variant calling with sam "
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"

function usage {
echo -e "usage: $(basename $0) -k CONFIG -f BAM -o OUTDIR [OPTIONS]"
exit
}

if [ ! $# -gt 3 ]; then usage ; fi

#INPUTS
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the NGSANE repository
        -f | --file    )        shift; f=$1 ;; # bam file
        -o | --output  )        shift; OUTDIR=$1 ;; # output directory
        --recover-from )        shift; RECOVERFROM=$1 ;; # attempt to recover from log file
        -h | --help    )        usage ;;
        * )                     echo "don't understand "$1
    esac
    shift
done

#PROGRAMS
. $CONFIG
. ${NGSANE_BASE}/conf/header.sh
. $CONFIG

################################################################################
CHECKPOINT="programs"

for MODULE in $MODULE_SAMVAR; do module load $MODULE; done  # save way to load modules that itself load other modules
export PATH=$PATH_SAMVAR:$PATH
module list
echo "PATH=$PATH"
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)
PATH_IGVTOOLS=$(dirname $(which igvtools.jar))

echo "[NOTE] set java parameters"
JAVAPARAMS="-Xmx"$(python -c "print int($MEMORY_SAMVAR*0.8)")"g -Djava.io.tmpdir="$TMP" -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1" 
unset _JAVA_OPTIONS
echo "JAVAPARAMS "$JAVAPARAMS

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--JAVA        --\n" $(java -Xmx200m version 2>&1)
[ -z "$(which java)" ] && echo "[ERROR] no java detected" && exit 1
echo -e "--samtools    --\n "$(samtools 2>&1 | head -n 3 | tail -n-2)
[ -z "$(which samtools)" ] && echo "[ERROR] no samtools detected" && exit 1
echo -e "--igvtools    --\n "$(java $JAVAPARAMS -jar $PATH_IGVTOOLS/igvtools.jar version 2>&1)
[ ! -f $PATH_IGVTOOLS/igvtools.jar ] && echo "[ERROR] no igvtools detected" && exit 1


echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="parameters"

# get basename of f
n=${f##*/}

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="recall files from tape"

if [ -n $DMGET ]; then 
    dmget -a $f; 
fi
    
echo -e "\n********* $CHECKPOINT\n"    
################################################################################
CHECKPOINT="remove duplicate reads $(date)"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

    samtools rmdup $f $OUTDIR/${n/bam/drm.bam}
    samtools index $OUTDIR/${n/bam/drm.bam}

    # mark checkpoint
    if [ -f $OUTDIR/${n/bam/drm.bam} ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi 

################################################################################
CHECKPOINT="call variants $(date)"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

    samtools mpileup -uf $FASTA -q1 -D $OUTDIR/${n/bam/drm.bam} |  bcftools view -vcg - >$OUTDIR/${n/bam/vcf}

    # mark checkpoint
    if [ -f $OUTDIR/${n/bam/vcf} ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi 

################################################################################
CHECKPOINT="convert bcf->vcf; index $(date)"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

    vcfutils.pl varFilter -D1000 -w0 -e0 $OUTDIR/${n/bam/vcf}  > $OUTDIR/${n/bam/clean.vcf}
 
    # mark checkpoint
    if [ -f $OUTDIR/${n/bam/clean.vcf} ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi 

################################################################################
CHECKPOINT="index vcf file for viewing in IGV $(date)"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

    java $JAVAPARAMS -jar $PATH_IGVTOOLS/igvtools.jar index $OUTDIR/${n/bam/clean.vcf}
    
    # mark checkpoint
    echo -e "\n********* $CHECKPOINT\n"
fi 

################################################################################
CHECKPOINT="cleanup"    

rm $OUTDIR/${n/bam/drm.bam} $OUTDIR/${n/bam/drm.bam}.bai
rm $OUTDIR/${n/bam/vcf}

echo -e "\n********* $CHECKPOINT\n"
################################################################################
echo ">>>>> Variant calling with sam - FINISHED"
echo ">>>>> enddate "`date`

