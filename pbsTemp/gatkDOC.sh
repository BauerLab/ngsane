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
echo ">>>>> gatkDOC.sh $*"


function usage {
echo -e "usage: $(basename $0) -k HISEQINF -f FASTQ -r REFERENCE -o OUTDIR [OPTIONS]

Script running the recalibration and realigment step (GATK)
It expects a bam file (*.asd.bam)

required:
  -k | --toolkit <path>     location of the HiSeqInf repository 
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
        -k | --toolkit )        shift; HISEQINF=$1 ;; # location of the HiSeqInf repository
        -t | --threads )        shift; THREADS=$1 ;; # number of CPUs to use
        -f | --fastq )          shift; f=$1 ;; # fastq file
        -r | --reference )      shift; FASTA=$1 ;; # reference genome
        -o | --outdir )         shift; OUT=$1 ;; # output dir  +++ THIS DOES NOT USE CONFIG OVERWRITE !! +++
	-G | --gene )           shift; GENE=$1 ;; # gene list
	-L | --list )           shift; LIST=$1 ;; # (optional) region of specific interest, e.g. targeted reseq
        -h | --help )           usage ;;
        * )                     usage
    esac
    shift
done

#INPUTS
#HISEQINF=$1   # location of the HiSeqInf repository
#f=$2          # fastq file
#FASTA=$3      # reference genome
#LIST=$4       # list of locations in refseq format
#OUT=$5        # output dir

TASK=""

if [ -n "$LIST" ]; then TASK="-L $LIST"; fi
if [ -n "$GENE" ]; then TASK="-geneList $GENE"; fi


#PROGRAMS
. $HISEQINF/pbsTemp/header.sh


module load jdk
export PATH=$PATH:$(basename $SAMTOOLS)
#JAVAPARAMS="-Xmx"$MYMEMORY"g" # -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1 -XX:MaxDirectMemorySize=4G"

n=`basename $f`


if [ -e $QOUT/$n.doc ]; then rm $QOUT/$n.doc* ]; fi


#java -jar /datastore/cmis/bau04c/SeqAna/apps/prod/Picard_svn/dist/CreateSequenceDictionary.jar R=/datastore/cmis/bau04c//SeqAna/reference/prod/GRCm38/GRCm38_chr.fasta O=/datastore/cmis/bau04c//SeqAna/reference/prod/GRCm38/GRCm38_chr.dict
# BEDtools has it's own genome index file


#calculate depthOfCoverage
echo "********* depthOfCoverage"

java $JAVAPARAMS -Djava.io.tmpdir=$TMP -jar $GATKJAR/GenomeAnalysisTK.jar \
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


#statistics for on target proportion
if [ -n "$LIST" ]; then
    echo "********* on target"
    # get reads overlapping exons+100 and get statistics
    # get how many reads where there to begin with
    $BEDTOOLS/slopBed -i $LIST -b 100 -g ${FASTA/fasta/bedtools.genome} | $BEDTOOLS/intersectBed \
	-abam $f -b stdin -u | $SAMTOOLS flagstat - > $OUT/$n.stats

    #$BEDTOOLS/windowBed -abam $f -b $LIST -w 100 -u | $SAMTOOLS flagstat - > $OUT/$n.stats
    echo "# overall" >> $OUT/$n.stats
    grep "in total (QC-passed reads + QC-failed reads)" $f.stats >> $OUT/$n.stats
    grep "properly paired (" $f.stats >> $OUT/$n.stats

    echo "# on target 0" >> $OUT/$n.stats
    $BEDTOOLS/intersectBed -abam $f -b $LIST  -u | $SAMTOOLS flagstat - >> $OUT/$n.stats
    echo "# on target 200" >> $OUT/$n.stats
    $BEDTOOLS/slopBed -i $LIST -b 100 -g ${FASTA/fasta/bedtools.genome} | $BEDTOOLS/intersectBed \
	-abam $f -b stdin -u | $SAMTOOLS flagstat - >> $OUT/$n.stats


fi

echo ">>>>> determine the depth of coverage with GATK - FINISHED"
echo ">>>>> enddate "`date`
