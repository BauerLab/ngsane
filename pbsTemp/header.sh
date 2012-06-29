


BINQSUB=/data/noflush/bau04c/SeqAna/apps/prod/seqaninf/bin/myqsub
BWA=/data/noflush/bau04c/SeqAna/apps/prod/bwa_git/bwa
SAMTOOLS=/data/noflush/bau04c/SeqAna/apps/prod/samtools_svn/samtools
IGVTOOLS=/data/noflush/bau04c/SeqAna/apps/prod/IGVTools/igvtools.jar
PICARD=/data/noflush/bau04c/SeqAna/apps/prod/Picard_svn/dist/
SAMSTAT=/data/noflush/bau04c/SeqAna/apps/dev/samstat/src/samstat
GATKHOME=/data/noflush/bau04c/SeqAna/apps/prod/gatk_git
GATKJAR=$GATKHOME/dist/
RSCRIPT=/apps/R/2.14.1/bin/Rscript # module load R
FASTQC=/data/noflush/bau04c/SeqAna/apps/prod/FastQC/fastqc
FASTXTK="/clusterdata/hiseq_apps/bin/devel/fastx_toolkit/"
#TMP=$TMPDIR
TMP=/home/cmis/bau04c/TMP
VCFTOOLS="/clusterdata/hiseq_apps/bin/freeze001/VCFtools_0.1.3.2/bin"
SAMUTILS="/clusterdata/hiseq_apps/bin/freeze001/tabix-0.2.3"
BEDTOOLS="/clusterdata/hiseq_apps/bin/freeze001/BEDTools-Version-2.11.2/bin/"
IMGMAGCONVERT=/usr/bin/convert # imageMagick
ANNOVAR="/clusterdata/hiseq_apps/bin/freeze001/annovar"

TOPHAT="/clusterdata/hiseq_apps/bin/freeze001/tophat-1.3.2.Linux_x86_64/tophat"
CUFFLINKSHOME="/clusterdata/hiseq_apps/bin/freeze001/cufflinks-1.1.0.Linux_x86_64"


RRBSMAP="/clusterdata/hiseq_apps/bin/devel/rrbsmap-1.5/rrbsmap"
MACS="/clusterdata/hiseq_apps/bin/devel/MACS_git"
PEAKFINDER="/clusterdata/hiseq_apps/bin/devel/vancouvershortr_svn/"

#BOWTIE="/clusterdata/hiseq_apps/bin/devel/bowtie-0.12.7/bowtie"
BOWTIETWO="/data/noflush/bau04c/SeqAna/apps/prod/bowtie2-2.0.0-beta6/"
CUTADAPT="/clusterdata/hiseq_apps/bin/devel/cutadapt/cutadapt-0.9/cutadapt"
VIENNA="/clusterdata/hiseq_apps/bin/devel/ViennaRNA/bin/"
UNAFOLD="/clusterdata/hiseq_apps/bin/devel/unafold/bin/"

# Task names
TASKFASTQC="fastQC"
TASKBWA="bwa"
TASKBOWTIE="bowtie"
TASKRCA="reCalAln"
TASKMERGE="merged"
TASKVAR="variant"
TASKDINC="dindelC"
TASKDINS="dindelS"
TASKDIN="dindel"
TASKSNP="snp"
TASKIND="indel"
TASKDOWN="downsample"
TASKDOC="coverage"
TASKDIFFEXP="diffexp"
TASKTOPHAT="tophat"
TASKCUFF="cufflinks"
TASKCUFFDIFF="cuffdiff"
TASKRRBS="rrbs"
TASKMACS="macs"
TASKANNOVAR="annovar"

#Fileabb
READONE="read1"
READTWO="read2"
FASTQ="fastq.gz"
ALN="aln" # aligned 
ASD="asd" # aligned sorted duplicate-removed
ASR="asdrr" # aligned sorted duplicate-removed raligned reacalibrated

