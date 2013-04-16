##############################################################
# System info
##############################################################
SUBMISSIONSYSTEM="SGE"
DMGET=""
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

MODULE_FASTQC="marcow/fastqc/prebuilt/0.10.1"
PATH_FASTQC=
MODULE_LATEX=
#PATH_LATEX="/data/flush/apps/texlive/2012/bin/x86_64-linux/"
PATH_LATEX=

MODULE_BWA="fabbus/bwa/0.7.3 fabbus/R/2.15.3 fabbus/samtools/0.1.19 fabbus/picard/1.89 fabbus/samstat/1.08 fabbus/igvtools/2.2.2"
PATH_BWA=$PATH_IGVTOOLS:$PATH_PICARD:$PATH_SAMSTAT

MODULE_BOWTIETWO="fabbus/R/2.15.3 fabbus/bowtie2/2.1.0 fabbus/samtools/0.1.19"
PATH_BOWTIETWO=$PATH_IGVTOOLS:$PATH_PICARD:$PATH_SAMSTAT

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

#TOPHAT=$DATASTORE/SeqAna/apps/prod/tophat-2.0.4.Linux_x86_64/tophat2 
TOPHAT=tophat/2.0.4b
CUFFLINKS=cufflinks/2.0.2
CUFFLINKSHOME=$DATASTORE/SeqAna/apps/prod/cufflinks-2.0.2.Linux_x86_64


RRBSMAP="/clusterdata/hiseq_apps/bin/devel/rrbsmap-1.5/rrbsmap"
MACS="/clusterdata/hiseq_apps/bin/devel/MACS_git"
PEAKFINDER="/clusterdata/hiseq_apps/bin/devel/vancouvershortr_svn/"
CUTADAPT="/home/cmis/bau04c/SeqAna/apps/prod/cutadapt-1.2.1/bin/cutadapt"

#
# Need to use 2.0.0b6 or higher because of --sam-RG
# only cherax has this install so far
#BOWTIETWO="$DATASTORE/SeqAna/apps/prod/bowtie2-2.0.0-beta6/" <- buggy
#BOWTIETWO="/apps/bowtie/2.0.0b6/"
#BOWTIETWO="/datastore/cmis/bau04c/SeqAna/apps/prod/bowtie2-2.0.0-beta7-source"
#BOWTIETWO=bowtie/2.0.0b5


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

#Fileabb
READONE="read1"
READTWO="read2"
FASTQ="fastq.gz"
ALN="aln" # aligned 
ASD="asd" # aligned sorted duplicate-removed
ASR="asdrr" # aligned sorted duplicate-removed raligned reacalibrated

#############
# On Cherax
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
