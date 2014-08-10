#!/bin/bash -e

# Script running downsample
# QC:
# author: Denis C. Bauer
# date: Sept.2011

echo ">>>>> Downsample"
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"

function usage {
echo -e "usage: $(basename $0)

required:
  -k | --toolkit <path>     location of the NGSANE repository 
  -i | --input <file>       bam file
  -o | --outdir <path>      output dir
"
exit
}

#if [ ! $# -gt 2 ]; then usage ; fi

#INPUTS
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the NGSANE repository
        -i | --input )          shift; f=$1 ;; # bam file
        -o | --outdir )         shift; OUT=$1 ;; # output dir
        -r | --reference )      shift; FASTA=$1 ;; # reference genome
        -s | --downsample )     shift; READNUMBER=$1 ;; #readnumber
        --recover-from )        shift; NGSANE_RECOVERFROM=$1 ;; # attempt to recover from log file                                                  
        -h | --help )           usage ;;
        * )                     usage
    esac
    shift
done

#PROGRAMS
. $CONFIG
. ${NGSANE_BASE}/conf/header.sh
. $CONFIG

################################################################################
NGSANE_CHECKPOINT_INIT "programs"
# save way to load modules that itself loads other modules
hash module 2>/dev/null && for MODULE in $MODULE_DOWNSAMPLE; do module load $MODULE; done && module list 

export PATH=$PATH_DOWNSAMPLE:$PATH
echo $PATH

#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)
[ -z "$PATH_GATK" ] && PATH_GATK=$(dirname $(which GenomeAnalysisTK.jar))

echo "[NOTE] set java parameters"
JAVAPARAMS="-Xmx"$(python -c "print int($MEMORY_DOWNSAMPLE*0.8)")"g -Djava.io.tmpdir="$TMP" -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1" 
unset _JAVA_OPTIONS
echo "JAVAPARAMS "$JAVAPARAMS

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--JAVA        --\n" $(java -Xmx200m -version 2>&1)
[ -z "$(which java)" ] && echo "[ERROR] no java detected" && exit 1
echo -e "--samtools    --\n "$(samtools 2>&1 | head -n 3 | tail -n-2)
[ -z "$(which samtools)" ] && echo "[ERROR] no samtools detected" && exit 1
echo -e "--PICARD      --\n "$(java $JAVAPARAMS -jar $PATH_PICARD/MarkDuplicates.jar --version 2>&1)
[ ! -f $PATH_PICARD/MarkDuplicates.jar ] && echo "[ERROR] no picard detected" && exit 1


NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "parameters"

# get basename of f
n=${f##*/}

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "recall files from tape"

if [ -n "$DMGET" ]; then
	dmget -a ${f}
fi
    
NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "extract properly paired none duplicate"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    #-F 0x0400
    samtools view -f 0x0002 -h -b $f > $OUT/${n/bam/pn.bam}
    
    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUT/${n/bam/pn.bam}

fi 

################################################################################
NGSANE_CHECKPOINT_INIT "downsample"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    READS=`samtools view -f 0x0002 -c $f`
    echo "[NOTE] number of properly paired none duplicate: $READS"
    
    PROB=`echo "$READNUMBER/$READS" | bc -l`
    echo $PROB
    java $JAVAPARAMS -jar $PICARD/DownsampleSam.jar \
        INPUT=$OUT/${n/bam/pn.bam} \
        OUTPUT=$OUT/${n/bam/pns.bam} \
        RANDOM_SEED=1 \
        VALIDATION_STRINGENCY=LENIENT \
        PROBABILITY=$PROB
    
    samtools index $OUT/${n/bam/pns.bam}
 
    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUT/${n/bam/pns.bam}

fi 

################################################################################
NGSANE_CHECKPOINT_INIT "statistics"
   
if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then
 
    samtools flagstat $OUT/${n/bam/pns.bam} > $OUT/${n/bam/pns.bam}.stats
    
    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUT/${n/bam/pns.bam}.stats

fi

################################################################################
NGSANE_CHECKPOINT_INIT "cleanup"

[ -e $OUT/${n/bam/pn.bam} ] && rm $OUT/${n/bam/pn.bam}

NGSANE_CHECKPOINT_CHECK
################################################################################
echo ">>>>> Downsample - FINISHED"
echo ">>>>> enddate "`date`
  
