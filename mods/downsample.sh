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

#DEFAULTS

#INPUTS
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the NGSANE repository
        -i | --input )          shift; f=$1 ;; # bam file
        -o | --outdir )         shift; OUT=$1 ;; # output dir
        -r | --reference )      shift; FASTA=$1 ;; # reference genome
        -s | --downsample )     shift; READNUMBER=$1 ;; #readnumber
        --recover-from )        shift; RECOVERFROM=$1 ;; # attempt to recover from log file                                                  
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
CHECKPOINT="programs"
for MODULE in $MODULE_DOWNSAMPLE; do module load $MODULE; done  # save way to load modules that itself load other modules
export PATH=$PATH_DOWNSAMPLE:$PATH
module list
echo $PATH
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)
PATH_IGVTOOLS=$(dirname $(which igvtools.jar))
PATH_GATK=$(dirname $(which GenomeAnalysisTK.jar))

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--samtools    --\n "$(samtools 2>&1 | head -n 3 | tail -n-2)
[ -z "$(which samtools)" ] && echo "[ERROR] no samtools detected" && exit 1
echo -e "--igvtools    --\n "$(java -jar $JAVAPARAMS $PATH_IGVTOOLS/igvtools.jar version 2>&1)
[ ! -f $PATH_IGVTOOLS/igvtools.jar ] && echo "[ERROR] no igvtools detected" && exit 1
echo -e "--PICARD      --\n "$(java -jar $JAVAPARAMS $PATH_PICARD/MarkDuplicates.jar --version 2>&1)
[ ! -f $PATH_PICARD/MarkDuplicates.jar ] && echo "[ERROR] no picard detected" && exit 1

echo "[NOTE] set java parameters"
JAVAPARAMS="-Xmx"$(python -c "print int($MEMORY_DOWNSAMPLE*0.8)")"g -Djava.io.tmpdir="$TMP" -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1" 
unset _JAVA_OPTIONS
echo "JAVAPARAMS "$JAVAPARAMS

echo -e "\n********* $CHECKPOINT"
################################################################################
CHECKPOINT="parameters"

# get basename of f
n=${f##*/}

################################################################################
CHECKPOINT="recall files from tape"

if [ -n "$DMGET" ]; then
	dmget -a ${f}
fi
    
echo -e "\n********* $CHECKPOINT"    
################################################################################
CHECKPOINT="extract properly paired none duplicate"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep "********* $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

    #-F 0x0400
    samtools view -f 0x0002 -h -b $f > $OUT/${n/bam/pn.bam}
    
    # mark checkpoint
    [ -f $OUT/${n/bam/pn.bam} ] && echo -e "\n********* $CHECKPOINT" && unset RECOVERFROM
fi 

################################################################################
CHECKPOINT="downsample"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep "********* $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

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
    [ -f $OUT/${n/bam/pns.bam} ] && echo -e "\n********* $CHECKPOINT" && unset RECOVERFROM
fi 

################################################################################
CHECKPOINT="statistics"
   
if [[ -n "$RECOVERFROM" ]] && [[ $(grep "********* $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
 
    samtools flagstat $OUT/${n/bam/pns.bam} > $OUT/${n/bam/pns.bam}.stats
    
    # mark checkpoint
    [ -f $OUT/${n/bam/pns.bam}.stats ] && echo -e "\n********* $CHECKPOINT" && unset RECOVERFROM
fi

################################################################################
CHECKPOINT="coverage track"
   
if [[ -n "$RECOVERFROM" ]] && [[ $(grep "********* $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
 
    GENOME=$(echo $FASTA| sed 's/.fasta/.genome/' | sed 's/.fa/.genome/' )
    java $JAVAPARAMS -jar $IGVTOOLS count $OUT/${n/bam/pns.bam} \
        $OUT/${n/bam/pns.bam.cov.tdf} $GENOME

    # mark checkpoint
    [ -f $OUT/${n/bam/pns.bam.cov.tdf} ] && echo -e "\n********* $CHECKPOINT" && unset RECOVERFROM
fi

################################################################################
CHECKPOINT="cleanup"

[ -e $OUT/${n/bam/pn.bam} ] && rm $OUT/${n/bam/pn.bam}

echo -e "\n********* $CHECKPOINT"
################################################################################
echo ">>>>> Downsample - FINISHED"
echo ">>>>> enddate "`date`
  