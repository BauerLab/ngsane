##############################################################
# System info
##############################################################
SUBMISSIONSYSTEM="SGE"
DMGET=""

##############################################################
# SUN GRID ENGINE specific workaround for BUG (SGE 6.2u5)
##############################################################
. /etc/profile.d/modules.sh

##############################################################
# PROGRAM PATHS
##############################################################
QSUB=pbsTemp.sh
BINQSUB=jobsubmission.sh

#Additional programs not available as module
PATH_SAMTOOLS=
PATH_IGVTOOLS=/share/ClusterShare/software/contrib/fabbus/igvtools/2.2.2/
PATH_PICARD=/share/ClusterShare/software/contrib/fabbus/picard/1.89/
PATH_SAMSTAT=/share/ClusterShare/software/contrib/fabbus/samstat/1.08/

# fastQC
MODULE_FASTQC="fabbus/fastqc/0.10.1"
PATH_FASTQC=
MODULE_LATEX=
#PATH_LATEX="/data/flush/apps/texlive/2012/bin/x86_64-linux/"
PATH_LATEX=

# BWA
MODULE_BWA="fabbus/bwa/0.7.3 fabbus/R/2.15.3 fabbus/samtools/0.1.19 fabbus/picard/1.89 fabbus/samstat/1.08 fabbus/igvtools/2.2.2"
PATH_BWA=$PATH_IGVTOOLS:$PATH_PICARD:$PATH_SAMSTAT

#Bowtie2
MODULE_BOWTIETWO="fabbus/R/2.15.3 fabbus/bowtie2/2.1.0 fabbus/samtools/0.1.19"
PATH_BOWTIETWO=$PATH_IGVTOOLS:$PATH_PICARD:$PATH_SAMSTAT

#Trim adapter with CUTADAPT
MODULE_CUTADAPT="fabbus/cutadapt/1.2.1"
PATH_CUTADAPT=""

#Trim adapter with TRIMGALORE
MODULE_TRIMGALORE="fabbus/cutadapt/1.2.1 fabbus/trimgalore/0.2.8"
PATH_TRIMGALORE=""
TRIMGALORE_ADAPTER1=""
TRIMGALORE_ADAPTER2=""

#Snp calling with GATK
PATH_GATKHOME=$DATASTORE/SeqAna/apps/prod/gatk_git
#GATKHOME=$DATASTORE/SeqAna/apps/dev/gatk_git
PATH_GATKJAR=$PATH_GATKHOME/dist/
MODULE_GATKSNP="fabbus/R/2.15.3"
PATH_GATKSNP=$PATH_GATKHOME:$PATH_GATKJAR:$PATH_IGVTOOLS

# Tophat and cufflinks
MODULE_TOPHATCUFF="fabbus/python/2.7.3 fabbus/R/2.15.3 fabbus/bowtie2/2.1.0 fabbus/tophat/2.0.8b fabbus/cufflinks/2.1.1 fabbis/samtools/0.1.19"
PATH_TOPHATCUFF=$PATH_IGVTOOLS:$PATH_PICARD:$PATH_SAMSTAT

#this gzip waits for the file to migrate completly before unzipping it
#GZIP=$DATASTORE/SeqAna/apps/prod/mygzip/
GATKHOME=$DATASTORE/SeqAna/apps/prod/gatk_git
#GATKHOME=$DATASTORE/SeqAna/apps/dev/gatk_git
GATKJAR=$GATKHOME/dist/
RSCRIPT=/apps/R/2.14.1/bin/Rscript # module load R
FASTXTK="/clusterdata/hiseq_apps/bin/devel/fastx_toolkit/"
#TMP=$TMPDIR
TMP=$(pwd)/tmp
#TMP=$DATASTORE/TMP
VCFTOOLS="/clusterdata/hiseq_apps/bin/freeze001/VCFtools_0.1.3.2/bin"
SAMUTILS="/clusterdata/hiseq_apps/bin/freeze001/tabix-0.2.3"
BEDTOOLS=$DATASTORE/SeqAna/apps/prod/bedtools/bin/
IMGMAGCONVERT=/usr/bin/convert # imageMagick
ANNOVAR="/clusterdata/hiseq_apps/bin/freeze001/annovar"

RRBSMAP="/clusterdata/hiseq_apps/bin/devel/rrbsmap-1.5/rrbsmap"
MACS="/clusterdata/hiseq_apps/bin/devel/MACS_git"
PEAKFINDER="/clusterdata/hiseq_apps/bin/devel/vancouvershortr_svn/"

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
TASKSAMVAR="samvar"
TASKCUTADAPT="cutadapt"

#Fileabb
READONE="read1"
READTWO="read2"
FASTQ="fastq.gz"
ALN="aln" # aligned 
ASD="asd" # aligned sorted duplicate-removed
ASR="asdrr" # aligned sorted duplicate-removed raligned reacalibrated

#############
# On Wolfpack
#FASTQC
WALLTIME_FASTQC=10:00:00
MEMORY_FASTQC=20
CPU_FASTQC=16
NODES_FASTQC="nodes=2:ppn=8"
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
WALLTIME_BAMANN=5:00:00
MEMORY_BAMANN=32
CPU_BAMANN=1
NODES_BAMANN="nodes=1:ppn=1"
#GATK VARCALL
WALLTIME_VAR=100:00:00
MEMORY_VAR=20
CPU_VAR=1
NODES_VAR="nodes=1:ppn=1"
#SAMTOOLS VARCALL
WALLTIME_SAMVAR=40:00:00
MEMORY_SAMVAR=40
CPU_SAMVAR=1
NODES_SAMVAR="nodes=1:ppn=1"
#CUTADAPT
WALLTIME_CUTADAPT=4:00:00
MEMORY_CUTADAPT=40
CPU_CUTADAPT=1
NODES_CUTADAPT="nodes=1:ppn=1"
#TRIMGALORE
WALLTIME_TRIMGALORE=4:00:00
MEMORY_TRIMGALORE=40
CPU_TRIMGALORE=1
NODES_TRIMGALORE="nodes=1:ppn=1"
