#!/bin/bash

# Script running coverage info per list in GATK
# author: Denis C. Bauer
# date: Jan.2011

# messages to look out for -- relevant for the QC.sh script:
# QCVARIABLES,
# RESULTFILENAME <DIR>/$TASK_GATKDOC/<SAMPLE>.doc

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
        -f | --fastq )          shift; f=$1 ;; # fastq file
        -o | --outdir )         shift; OUTDIR=$1 ;; # output dir  +++ THIS DOES NOT USE CONFIG OVERWRITE !! +++
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
hash module 2>/dev/null && for MODULE in $MODULE_GATKDOC; do module load $MODULE; done && module list 

export PATH=$PATH_GATK:$PATH
echo "PATH=$PATH"
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)
[ -z "$PATH_GATK" ] && PATH_GATK=$(dirname $(which GenomeAnalysisTK.jar))

echo "[NOTE] set java parameters"
JAVAPARAMS="-Xmx"$(python -c "print int($MEMORY_GATKDOC*0.75)")"g -Djava.io.tmpdir="$TMP"  -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1" 
unset _JAVA_OPTIONS
echo "JAVAPARAMS "$JAVAPARAMS

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--JAVA        --\n" $(java -Xmx200m -version 2>&1)
[ -z "$(which java)" ] && echo "[ERROR] no java detected" && exit 1
echo -e "--samtools    --\n "$(samtools 2>&1 | head -n 3 | tail -n-2)
[ -z "$(which samtools)" ] && echo "[ERROR] no samtools detected" && exit 1
echo -e "--bedtools    --\n "$(bedtools --version)
[ -z "$(which bedtools)" ] && echo "[ERROR] no bedtools detected" && exit 1
echo -e "--GATK        --\n "$(java $JAVAPARAMS -jar $PATH_GATK/GenomeAnalysisTK.jar --version)
[ ! -f $PATH_GATK/GenomeAnalysisTK.jar ] && echo "[ERROR] no GATK detected" && exit 1

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "parameters"

# get basename of f
n=${f##*/}
SAMPLE=${n/%$ASR.bam/}

GATK_TASK=""

if [ -n "$LIST" ]; then GATK_TASK="-L $LIST"; fi
if [ -n "$GENE" ]; then GATK_TASK="-geneList $GENE"; fi

# delete old bam files unless attempting to recover
if [ -z "$NGSANE_RECOVERFROM" ]; then
    if [ -e $QOUT/$SAMPLE.doc ]; then rm $QOUT/$SAMPLE.doc* ]; fi
fi

if [ -z "$FASTA" ] || [ ! -f $FASTA ]; then
    echo "[ERROR] no reference provided (FASTA)"
    exit 1
else
    echo "[NOTE] Reference: $FASTA"
fi

GENOME_CHROMSIZES=${FASTA%.*}.chrom.sizes

if [[ $(which GenomeAnalysisTK.jar) =~ "2.8" && -z "$FORCEINTERVALS" ]]; then 
        echo "[NOTE] new GATK parallel, this will not generate IntervalStatistics"
        PARALLELENCT="-nct $CPU_GATKDOC"
		PARALLELENT="--omitIntervalStatistics -nt $CPU_GATKDOC"
fi

#java -jar /datastore/cmis/bau04c/SeqAna/apps/prod/Picard_svn/dist/CreateSequenceDictionary.jar R=/datastore/cmis/bau04c//SeqAna/reference/prod/GRCm38/GRCm38_chr.fasta O=/datastore/cmis/bau04c//SeqAna/reference/prod/GRCm38/GRCm38_chr.dict
# BEDtools has it's own genome index file

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "calculate depthOfCoverage"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then
    
    java $JAVAPARAMS -jar $PATH_GATK/GenomeAnalysisTK.jar \
		$ADDPARAMGATKDOC \
        -T DepthOfCoverage \
        --minMappingQuality 10 \
        --minBaseQuality 10 \
        -R $FASTA \
        -I $f \
        -o $OUTDIR/$SAMPLE.doc \
        $GATK_TASK --omitDepthOutputAtEachBase \
        -ct 10 -ct 20 -ct 50 -ct 100 -ct 500 \
        --start 1 \
        --stop 100 \
		$PARALLELENT \
        --nBins 10
    
    
    #    -nt $THREADS \
    #    --calculateCoverageOverGenes $LIST \
    #    --minMappingQuality 20 \
    #    --minBaseQuality 20 \
    
    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE.doc.sample_statistics

fi 

################################################################################
NGSANE_CHECKPOINT_INIT "on target"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    #statistics for on target proportion
    if [ -n "$LIST" ]; then
        # get reads overlapping exons+100 and get statistics
        # get how many reads where there to begin with
        slopBed -i $LIST -b 100 -g $GENOME_CHROMSIZES | intersectBed \
    	-abam $f -b stdin -u | samtools flagstat - > $OUTDIR/$SAMPLE$ASR.bam.stats
    
        #windowBed -abam $f -b $LIST -w 100 -u | samtools flagstat - > $OUTDIR/$SAMPLE$ASR.bam.stats
        echo "# overall" >> $OUTDIR/$SAMPLE$ASR.bam.stats
        grep "in total (QC-passed reads + QC-failed reads)" $f.stats >> $OUTDIR/$SAMPLE$ASR.bam.stats
        grep "properly paired (" $f.stats >> $OUTDIR/$SAMPLE$ASR.bam.stats
    
        echo "# on target 0" >> $OUTDIR/$SAMPLE$ASR.bam.stats
        intersectBed -abam $f -b $LIST  -u | samtools flagstat - >> $OUTDIR/$SAMPLE$ASR.bam.stats
        echo "# on target 200" >> $OUTDIR/$SAMPLE$ASR.bam.stats
        slopBed -i $LIST -b 100 -g $GENOME_CHROMSIZES | intersectBed \
    	-abam $f -b stdin -u | samtools flagstat - >> $OUTDIR/$SAMPLE$ASR.bam.stats

        # mark checkpoint
        NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE$ASR.bam.stats

    fi
fi

################################################################################
[ -e $OUTDIR/$SAMPLE.doc.dummy ] && rm $OUTDIR/$SAMPLE.doc.dummy
echo ">>>>> determine the depth of coverage with GATK - FINISHED"
echo ">>>>> enddate "`date`
