##############################################################
# System info
##############################################################
SUBMISSIONSYSTEM="SGE"                                # or PBS
DMGET=""                    # or Yes when storing data on tape
TMP=/share/Temp                                      # TMP dir

##############################################################
# SUN GRID ENGINE specific workaround for BUG (SGE 6.2u5)
##############################################################
. /etc/profile.d/modules.sh

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
TASKTRIMGALORE="trimgalore"
TASKHICUP="hicup"
TASKHICLIB="hiclib"
TASKFASTQ2SANGER="sanger"
TASKHOMERHIC="homerhic"
TASKWIGGLER="wiggler"
TASKTRIMMOMATIC="trimmomatic"

TASKTRINITY="trinity"
 TASKINCHWORM="inchworm"
 TASKCHRYSALIS="chrysalis"
 TASKBUTTERFLY="butterfly"

##############################################################
# PROGRAM PATHS
##############################################################
QSUB=prepareJobSubmission.sh
BINQSUB=jobsubmission.sh
QSUBEXTRA="-m e -M `whoami`@garvan.unsw.edu.au" # email notification

#Additional programs not available as module
PATH_SAMTOOLS=
PATH_IGVTOOLS=
PATH_PICARD=
PATH_SAMSTAT=

##############################################################
# gzip alternatives, e.g. pigz
MODULE_GZIP="gi/pigz/2.3"
GZIP="pigz -11"
[ -n "$MODULE_GZIP" ] && module load $MODULE_GZIP

##############################################################
# FASTQC
WALLTIME_FASTQC=10:00:00
MEMORY_FASTQC=20
CPU_FASTQC=16
NODES_FASTQC="nodes=2:ppn=8"

MODULE_FASTQC="gi/fastqc/0.10.1"
PATH_FASTQC=
MODULE_LATEX=
PATH_LATEX=

##############################################################
# FASTQ2SANGER (requires perl)
#
WALLTIME_FASTQ2SANGER=10:00:00
MEMORY_FASTQ2SANGER=20
CPU_FASTQ2SANGER=1
NODES_FASTQ2SANGER="nodes=1:ppn=1"

MODULE_FASTQ2SANGER=
PATH_FASTQ2SANGER=
FASTQ2SANGER_SOURCEFORMAT=

##############################################################
# SAMTOOLS 
WALLTIME_SAMVAR=48:00:00
MEMORY_SAMVAR=40
CPU_SAMVAR=1
NODES_SAMVAR="nodes=1:ppn=1"
MODULE_SAMVAR="gi/samtools/0.1.19 gi/igvtools/2.3.5 gi/gatk/2.5"

##############################################################
# BWA
WALLTIME_BWA=48:00:00
MEMORY_BWA=50
CPU_BWA=16
NODES_BWA="nodes=4:ppn=8"

MODULE_BWA="gi/bwa/0.7.4 gi/R/3.0.0 gi/samtools/0.1.19 gi/picard-tools/1.91 gi/samstat/1.08 gi/igvtools/2.3.5 gi/samstat/1.08"
PATH_BWA=

##############################################################
# Bowtie2 & Bowtie1
WALLTIME_BOWTIE=72:00:00
MEMORY_BOWTIE=60
CPU_BOWTIE=16
NODES_BOWTIE="nodes=1:ppn=8"

MODULE_BOWTIETWO="gi/R/3.0.0 gi/bowtie/2.1.0 gi/samtools/0.1.19 gi/igvtools/2.3.5 gi/picard-tools/1.91 gi/samstat/1.08"
PATH_BOWTIETWO=
BOWTIE2_INDEX=

MODULE_BOWTIE="gi/R/3.0.0 gi/bowtie/1.i0.0 gi/samtools/0.1.19 gi/igvtools/2.3.5 gi/picard-tools/1.91 gi/samstat/1.08 gi/homer/4.2"
PATH_BOWTIE=
BOWTIE_INDEX=

##############################################################
# Wiggler
# https://sites.google.com/site/anshulkundaje/projects/wiggler
WALLTIME_WIGGLER=10:00:00
MEMORY_WIGGLER=60
CPU_WIGGLER=8
NODES_WIGGLER="nodes=1:ppn=8"

MODULE_WIGGLER="gi/samtools/0.1.19 fabbus/matlab/mcr2010b fabbus/wiggler/2.0"
PATH_WIGGLER=

WIGGLER_UMAPDIR=
WIGGLERADDPARAMS=

##############################################################
# HOMER HIC
# http://biowhat.ucsd.edu/homer/index.html
WALLTIME_HOMERHIC=60:00:00
MEMORY_HOMERHIC=60
CPU_HOMERHIC=8
NODES_HOMERHIC="nodes=1:ppn=8"

MODULE_HOMERHIC="gi/R/3.0.0 fabbus/perl/5.14.2 fabbus/circos/0.62.1 gi/homer/4.2"
PATH_HOMERHIC=

##############################################################
# Trim adapter with CUTADAPT
WALLTIME_CUTADAPT=16:00:00
MEMORY_CUTADAPT=40
CPU_CUTADAPT=1
NODES_CUTADAPT="nodes=1:ppn=1"

MODULE_CUTADAPT="fabbus/cutadapt/1.2.1"
PATH_CUTADAPT=""

##############################################################
# Trim adapter with TRIMGALORE
WALLTIME_TRIMGALORE=16:00:00
MEMORY_TRIMGALORE=40
CPU_TRIMGALORE=1
NODES_TRIMGALORE="nodes=1:ppn=1"

MODULE_TRIMGALORE="gi/fastx_toolkit/0.0.13.2 fabbus/cutadapt/1.2.1 fabbus/trimgalore/0.2.8"
PATH_TRIMGALORE=""
TRIMGALORE_ADAPTER1=""
TRIMGALORE_ADAPTER2=""

##############################################################
# Trimming with TRIMMOMATIC
# http://www.usadellab.org/cms/index.php?page=trimmomatic
WALLTIME_TRIMMOMATIC=4:00:00
MEMORY_TRIMMOMATIC=40
CPU_TRIMMOMATIC=1
NODES_TRIMMOMATIC="nodes=1:ppn=1"

MODULE_TRIMMOMATIC="gi/trimmomatic/0.30"
PATH_TRIMMOMATIC=
TRIMMOMATICADDPARAMS=
TRIMMOMATICSTEPS=

##############################################################
# GATK suite
#Recal
WALLTIME_RECAL=60:00:00
MEMORY_RECAL=50
CPU_RECAL=8
NODES_RECAL="nodes=1:ppn=8"
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

MODULE_GATK="gi/R/3.0.0 gi/gatk/2.5 gi/samtools/0.1.19 gi/igvtools/2.3.5"
MODULE_GATKSNP="gi/R/3.0.0 gi/gatk/2.5 gi/igvtools/2.3.5"
PATH_GATKSNP=
PATH_GATKHOME=
PATH_GATKJAR=
##############################################################
# Tophat and cufflinks
WALLTIME_TOPHAT=192:00:00
MEMORY_TOPHAT=60
CPU_TOPHAT=16
NODES_TOPHAT="nodes=2:ppn=8"

MODULE_TOPHATCUFF="fabbus/python/2.7.3 gi/R/3.0.0 gi/bowtie/2.1.0 gi/tophat/2.0.9 gi/cufflinks/2.1.1 gi/samtools/0.1.19 gi/bedtools/2.17.0 gi/igvtools/2.3.5 gi/picard-tools/1.91 gi/samstat/1.08 hugfre/HTSeq/0.5.4p3 gi/rnaseqc/1.1.7 gi/bwa/0.7.4"
PATH_TOPHATCUFF=

##############################################################
# R (3.0.0)
# http://www.r-project.org/
WALLTIME_R=12:00:00
MEMORY_R=50
CPU_R=1
NODES_R="nodes=1:ppn=1"

MODULE_R=
PATH_R=
RSCRIPT=Rscript

##############################################################
# HICLIB
WALLTIME_HICLIB=48:00:00
MEMORY_HICLIB=60
CPU_HICLIB=32
NODES_HICLIB="nodes=1:ppn=8"
CPU_HICLIB_POSTCOMMAND=1
NODES_HICLIB_POSTCOMMAND="nodes=1:ppn=1"

MODULE_HICLIB="gi/bowtie/2.1.0 gi/samtools/0.1.19 gi/hdf5/1.8.10-patch1 fabbus/hiclib/30_04_13 gi/picard-tools/1.91"
PATH_HICLIB=
HICLIB_GAPFILE= #/share/ClusterShare/software/contrib/fabbus/hiclib/19_04_13/hiclib/gap.txt
HICLIB_RENZYMES=
HICLIB_READLENGTH=

##############################################################
# HICUP
WALLTIME_HICUP=48:00:00
MEMORY_HICUP=60
CPU_HICUP=32
NODES_HICUP="nodes=1:ppn=8"

MODULE_HICUP="fabbus/hicup/0.3.0 fabbus/python/2.7.3 fabbus/fit-hi-c/28_12_2012"
PATH_HICUP=
HICUP_RENZYMES=

FASTXTK=
VCFTOOLS=
SAMUTILS=
BEDTOOLS=
IMGMAGCONVERT=
ANNOVAR=

#Fileabb
READONE="R1"
READTWO="R2"
FASTQ="fastq.gz"
FASTA=            # fasta file usually from the reference genome
FASTA_CHROMDIR=   # folder containing individual fasta files for each chromosome of the reference genome 
UNM="unm" # unmapped
ALN="aln" # aligned 
MUL="mul" # non-unique aligned
ASD="asd" # aligned sorted duplicate-removed
ASR="asdrr" # aligned sorted duplicate-removed raligned reacalibrated

#############
# On Wolfpack
#ANNOTATING BAM
WALLTIME_BAMANN=5:00:00
MEMORY_BAMANN=32
CPU_BAMANN=1
NODES_BAMANN="nodes=1:ppn=1"
