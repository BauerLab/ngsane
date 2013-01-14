


BINQSUB=$DATASTORE/SeqAna/apps/prod/seqaninf/bin/myqsub
BWA=$DATASTORE/SeqAna/apps/prod/bwa_git/bwa
SAMTOOLS=$DATASTORE/SeqAna/apps/prod/samtools_svn/samtools
IGVTOOLS=$DATASTORE/SeqAna/apps/prod/IGVTools/igvtools.jar
PICARD=$DATASTORE/SeqAna/apps/prod/Picard_svn/dist/
SAMSTAT=$DATASTORE/SeqAna/apps/dev/samstat/src/samstat
GATKHOME=$DATASTORE/SeqAna/apps/prod/gatk_git
GATKJAR=$GATKHOME/dist/
RSCRIPT=/apps/R/2.14.1/bin/Rscript # module load R
FASTQC=$DATASTORE/SeqAna/apps/prod/FastQC/fastqc
FASTXTK="/clusterdata/hiseq_apps/bin/devel/fastx_toolkit/"
#TMP=$TMPDIR
TMP=$DATASTORE/TMP
VCFTOOLS="/clusterdata/hiseq_apps/bin/freeze001/VCFtools_0.1.3.2/bin"
SAMUTILS="/clusterdata/hiseq_apps/bin/freeze001/tabix-0.2.3"
BEDTOOLS=$DATASTORE/SeqAna/apps/prod/bedtools/bin/
IMGMAGCONVERT=/usr/bin/convert # imageMagick
ANNOVAR="/clusterdata/hiseq_apps/bin/freeze001/annovar"

#TOPHAT=$DATASTORE/SeqAna/apps/prod/tophat-2.0.4.Linux_x86_64/tophat2 
TOPHAT=tophat/2.0.4b
CUFFLINKS=cufflinks/2.0.2
CUFFLINKSHOME=$DATASTORE/SeqAna/apps/prod/cufflinks-2.0.2.Linux_x86_64


RRBSMAP="/clusterdata/hiseq_apps/bin/devel/rrbsmap-1.5/rrbsmap"
MACS="/clusterdata/hiseq_apps/bin/devel/MACS_git"
PEAKFINDER="/clusterdata/hiseq_apps/bin/devel/vancouvershortr_svn/"

#BOWTIETWO="$DATASTORE/SeqAna/apps/prod/bowtie2-2.0.0-beta6/" <- buggy
BOWTIETWO="/apps/bowtie/2.0.0b6/"
#BOWTIETWO="/datastore/cmis/bau04c/SeqAna/apps/prod/bowtie2-2.0.0-beta7-source"
#BOWTIETWO=bowtie/2.0.0b5

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
TASKBAMANN="bamann"

#Fileabb
READONE="read1"
READTWO="read2"
FASTQ="fastq.gz"
ALN="aln" # aligned 
ASD="asd" # aligned sorted duplicate-removed
ASR="asdrr" # aligned sorted duplicate-removed raligned reacalibrated

#############
# On Cherax
#BWA
WALLTIME_BWA=50:00:00
MEMORY_BWA=50
CPU_BWA=32
NODES_BWA="nodes=4:ppn=8"
#Botie
WALLTIME_BOWTIE=10:00:00
MEMORY_BOWTIE=60
CPU_BOWTIE=8
NODES_BOWTIE="nodes=1:ppn=8"
#TOPHAT
WALLTIME_TOPHAT=60:00:00
MEMORY_TOPHAT=50
CPU_TOPHAT=8
NODES_TOPHAT="nodes=1:ppn=8"
#Recal
WALLTIME_RECAL=60:00:00
MEMORY_RECAL=50
CPU_RECAL=8
NODES_RECAL="nodes=1:ppn=8" 
#COVERAGE
WALLTIME_GATKDOC=50:00:00
MEMORY_GATKDOC=50
CPU_GATKDOC=1
NODES_GATKDOC="nodes=1:ppn=1"
#ANNOTATING BAM
WALLTIME_TASKBAMANN=05:00:00
MEMORY_TASKBAMANN=8
CPU_TASKBAMANN=1
NODES_TASKBAMANN="nodes=1:ppn=1"