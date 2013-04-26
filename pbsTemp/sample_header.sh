##############################################################
# System info
##############################################################
SUBMISSIONSYSTEM=""                               # SGE or PBS
DMGET=""                    # or Yes when storing data on tape
TMP=$(pwd)/tmp                                       # TMP dir

##############################################################
# SUN GRID ENGINE specific workaround for BUG (SGE 6.2u5)
##############################################################
## uncomment if running on SGE
#. /etc/profile.d/modules.sh

##############################################################
# Task Names
##############################################################
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
TASKHICUP="hicup"
TASKHICLIB="hiclib"

##############################################################
# PROGRAM PATHS
##############################################################
QSUB=pbsTemp.sh
BINQSUB=jobsubmission.sh

#Additional programs not available as module
PATH_SAMTOOLS=
PATH_IGVTOOLS=
PATH_PICARD=
PATH_SAMSTAT=

##############################################################
# gzip alternatives, e.g. pigz
MODULE_GZIP=
GZIP=gzip			# command, e.g. gzip or pigz

##############################################################
# FASTQC
WALLTIME_FASTQC=10:00:00
MEMORY_FASTQC=20
CPU_FASTQC=16
NODES_FASTQC="nodes=2:ppn=8"

MODULE_FASTQC=
PATH_FASTQC=
MODULE_LATEX=
PATH_LATEX=

##############################################################
# SAMTOOLS 
WALLTIME_SAMVAR=40:00:00
MEMORY_SAMVAR=40
CPU_SAMVAR=1
NODES_SAMVAR="nodes=1:ppn=1"

##############################################################
# BWA
WALLTIME_BWA=50:00:00
MEMORY_BWA=50
CPU_BWA=32
NODES_BWA="nodes=4:ppn=8"

MODULE_BWA=
PATH_BWA=$PATH_IGVTOOLS:$PATH_PICARD:$PATH_SAMSTAT

##############################################################
# Bowtie2
WALLTIME_BOWTIE=10:00:00
MEMORY_BOWTIE=60
CPU_BOWTIE=8
NODES_BOWTIE="nodes=1:ppn=8"

MODULE_BOWTIETWO=
PATH_BOWTIETWO=$PATH_IGVTOOLS:$PATH_PICARD:$PATH_SAMSTAT
BOWTIE2_INDEX=

##############################################################
# Trim adapter with CUTADAPT
WALLTIME_CUTADAPT=4:00:00
MEMORY_CUTADAPT=40
CPU_CUTADAPT=1
NODES_CUTADAPT="nodes=1:ppn=1"

MODULE_CUTADAPT=
PATH_CUTADAPT=

##############################################################
# Trim adapter with TRIMGALORE
WALLTIME_TRIMGALORE=4:00:00
MEMORY_TRIMGALORE=40
CPU_TRIMGALORE=1
NODES_TRIMGALORE="nodes=1:ppn=1"

MODULE_TRIMGALORE=
PATH_TRIMGALORE=
TRIMGALORE_ADAPTER1=
TRIMGALORE_ADAPTER2=

##############################################################
# Snp calling with GATK
# COVERAGE
WALLTIME_GATKDOC=50:00:00
MEMORY_GATKDOC=50
CPU_GATKDOC=1
NODES_GATKDOC="nodes=1:ppn=1"
# GATK VARCALL
WALLTIME_VAR=100:00:00
MEMORY_VAR=20
CPU_VAR=1
NODES_VAR="nodes=1:ppn=1"

PATH_GATKHOME=
PATH_GATKJAR=
MODULE_GATKSNP=
PATH_GATKSNP=$PATH_GATKHOME:$PATH_GATKJAR:$PATH_IGVTOOLS

##############################################################
# Tophat and cufflinks
WALLTIME_TOPHAT=60:00:00
MEMORY_TOPHAT=50
CPU_TOPHAT=8
NODES_TOPHAT="nodes=1:ppn=8"

MODULE_TOPHATCUFF=
PATH_TOPHATCUFF=$PATH_IGVTOOLS:$PATH_PICARD:$PATH_SAMSTAT

##############################################################
# HICLIB
WALLTIME_HICLIB=10:00:00
MEMORY_HICLIB=60
CPU_HICLIB=8
NODES_HICLIB="nodes=1:ppn=8"

MODULE_HICLIB=
PATH_HICLIB=
HICLIB_GAPFILE=
HICLIB_RENZYMES=
HICLIB_READLENGTH=


##############################################################
# R
WALLTIME_R=1:00:00
MEMORY_R=10
CPU_R=1
NODES_R="nodes=1:ppn=1"

MODULE_R=
PATH_R=
RSCRIPT=Rscript

##############################################################
#this gzip waits for the file to migrate completly before unzipping it
#GZIP=$DATASTORE/SeqAna/apps/prod/mygzip/
#GATKHOME=$DATASTORE/SeqAna/apps/prod/gatk_git
#GATKHOME=$DATASTORE/SeqAna/apps/dev/gatk_git
#GATKJAR=$GATKHOME/dist/

FASTXTK="/clusterdata/hiseq_apps/bin/devel/fastx_toolkit/"

VCFTOOLS="/clusterdata/hiseq_apps/bin/freeze001/VCFtools_0.1.3.2/bin"
SAMUTILS="/clusterdata/hiseq_apps/bin/freeze001/tabix-0.2.3"
BEDTOOLS=$DATASTORE/SeqAna/apps/prod/bedtools/bin/
IMGMAGCONVERT=/usr/bin/convert # imageMagick
ANNOVAR="/clusterdata/hiseq_apps/bin/freeze001/annovar"

RRBSMAP="/clusterdata/hiseq_apps/bin/devel/rrbsmap-1.5/rrbsmap"
MACS="/clusterdata/hiseq_apps/bin/devel/MACS_git"
PEAKFINDER="/clusterdata/hiseq_apps/bin/devel/vancouvershortr_svn/"

VIENNA=
UNAFOLD=


#Fileabb
READONE="read1"
READTWO="read2"
FASTQ="fastq.gz"
FASTA=            # fasta file usually from the reference genome
FASTA_CHROMDIR=   # folder containing individual fasta files for each chromosome of the reference genome 
ALN="aln" # aligned 
ASD="asd" # aligned sorted duplicate-removed
ASR="asdrr" # aligned sorted duplicate-removed raligned reacalibrated

#############
# On Wolfpack
#Recal
WALLTIME_RECAL=60:00:00
MEMORY_RECAL=50
CPU_RECAL=8
NODES_RECAL="nodes=1:ppn=8" 
#ANNOTATING BAM
WALLTIME_BAMANN=5:00:00
MEMORY_BAMANN=32
CPU_BAMANN=1
NODES_BAMANN="nodes=1:ppn=1"
