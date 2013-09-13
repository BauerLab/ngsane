#!/bin/bash

# Script running coverage info per list in GATK
# author: Denis C. Bauer
# date: Jan.2011

# messages to look out for -- relevant for the QC.sh script:
# QCVARIABLES,

echo ">>>>> determine the depth of coverage with GATK"
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"


function usage {
echo -e "usage: $(basename $0) -k NGSANE -f FASTQ -r REFERENCE -o OUTDIR [OPTIONS]

Script running the recalibration and realigment step (GATK)
It expects a bam file (*.asd.bam)

required:
  -k | --toolkit <path>     location of the NGSANE repository 
  -f | --fastq <file>       fastq file
  -r | --reference <file>   reference genome
  -o | --outdir <path>      output dir

options:
  -l | --genelist <file>    gene list
  -i | --interval <file>    interval list
  -t | --threads <nr>       number of CPUs to use (default: 1)
  -R | --region <ps>        region of specific interest, e.g. targeted reseq
                             format chr:pos-pos
"
exit
}

if [ ! $# -gt 4 ]; then usage ; fi

#DEFAULTS
THREADS=1

#INPUTS
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the NGSANE repository
        -t | --threads )        shift; THREADS=$1 ;; # number of CPUs to use
        -f | --fastq )          shift; f=$1 ;; # fastq file
        -r | --reference )      shift; FASTA=$1 ;; # reference genome
        -o | --outdir )         shift; OUT=$1 ;; # output dir  +++ THIS DOES NOT USE CONFIG OVERWRITE !! +++
        -G | --gene )           shift; GENE=$1 ;; # gene list
        -L | --list )           shift; LIST=$1 ;; # (optional) region of specific interest, e.g. targeted reseq
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

for MODULE in $MODULE_GATK; do module load $MODULE; done  # save way to load modules that itself load other modules
export PATH=$PATH_GATK:$PATH
module list
echo "PATH=$PATH"
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)
PATH_GATK=$(dirname $(which GenomeAnalysisTK.jar))

echo "[NOTE] set java parameters"
JAVAPARAMS="-Xmx"$(python -c "print int($MEMORY_GATK*0.8)")"g -Djava.io.tmpdir="$TMP"  -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1" 
unset _JAVA_OPTIONS
echo "JAVAPARAMS "$JAVAPARAMS

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--JAVA        --\n" $(java -Xmx200m -version 2>&1)
[ -z "$(which java)" ] && echo "[ERROR] no java detected" && exit 1
echo -e "--samtools    --\n "$(samtools 2>&1 | head -n 3 | tail -n-2)
[ -z "$(which samtools)" ] && echo "[ERROR] no samtools detected" && exit 1
echo -e "--bedtools    --\n "$(bedtools --version)
[ -z "$(which bedtools)" ] && echo "[ERROR] no bedtools detected" && exit 1
echo -e "--GATK        --\n "$(java -jar $JAVAPARAMS $PATH_GATK/GenomeAnalysisTK.jar --version)
[ ! -f $PATH_GATK/GenomeAnalysisTK.jar ] && echo "[ERROR] no GATK detected" && exit 1

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="parameters"

# get basename of f
n=${f##*/}

TASK=""

if [ -n "$LIST" ]; then TASK="-L $LIST"; fi
if [ -n "$GENE" ]; then TASK="-geneList $GENE"; fi

# delete old bam files unless attempting to recover
if [ -z "$RECOVERFROM" ]; then
    if [ -e $QOUT/$n.doc ]; then rm $QOUT/$n.doc* ]; fi
fi

#java -jar /datastore/cmis/bau04c/SeqAna/apps/prod/Picard_svn/dist/CreateSequenceDictionary.jar R=/datastore/cmis/bau04c//SeqAna/reference/prod/GRCm38/GRCm38_chr.fasta O=/datastore/cmis/bau04c//SeqAna/reference/prod/GRCm38/GRCm38_chr.dict
# BEDtools has it's own genome index file

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="calculate depthOfCoverage"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
    
    java $JAVAPARAMS -jar $PATH_GATK/GenomeAnalysisTK.jar \
        -T DepthOfCoverage \
        --minMappingQuality 10 \
        --minBaseQuality 10 \
        -R $FASTA \
        -I $f \
        -o $OUT/$n.doc \
        $TASK --omitDepthOutputAtEachBase \
        -ct 10 -ct 20 -ct 50 -ct 100 -ct 500 \
        --start 1 \
        --stop 100 \
        --nBins 10
    
    
    #    -nt $THREADS \
    #    --calculateCoverageOverGenes $LIST \
    #    --minMappingQuality 20 \
    #    --minBaseQuality 20 \
    
    # mark checkpoint
    if [ -f $OUT/$n.doc ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi 

################################################################################
CHECKPOINT="on target"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

    #statistics for on target proportion
    if [ -n "$LIST" ]; then
        # get reads overlapping exons+100 and get statistics
        # get how many reads where there to begin with
        slopBed -i $LIST -b 100 -g ${FASTA/fasta/bedtools.genome} | intersectBed \
    	-abam $f -b stdin -u | $SAMTOOLS flagstat - > $OUT/$n.stats
    
        #windowBed -abam $f -b $LIST -w 100 -u | $SAMTOOLS flagstat - > $OUT/$n.stats
        echo "# overall" >> $OUT/$n.stats
        grep "in total (QC-passed reads + QC-failed reads)" $f.stats >> $OUT/$n.stats
        grep "properly paired (" $f.stats >> $OUT/$n.stats
    
        echo "# on target 0" >> $OUT/$n.stats
        intersectBed -abam $f -b $LIST  -u | $SAMTOOLS flagstat - >> $OUT/$n.stats
        echo "# on target 200" >> $OUT/$n.stats
        slopBed -i $LIST -b 100 -g ${FASTA/fasta/bedtools.genome} | intersectBed \
    	-abam $f -b stdin -u | $SAMTOOLS flagstat - >> $OUT/$n.stats

        # mark checkpoint
        if [ -f $OUT/$n.stats ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

    fi
fi

################################################################################
echo ">>>>> determine the depth of coverage with GATK - FINISHED"
echo ">>>>> enddate "`date`
