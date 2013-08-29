##############################################################
# System info
##############################################################
SUBMISSIONSYSTEM="PBS"                            # SGE or PBS
QUEUEWAIT=" -W depend=afterok:"                   # PBS
QUEUEWAITSEP=":"
#QUEUEWAIT=" -hold_jid "                          # SGE
#QUEUEWAITSEP=","       
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
TASKBOWTIE2="bowtie2"
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
TASKRRBSMAP="rrbs"
TASKMACS="macs"
TASKANNOVAR="annovar"
TASKBAMANN="bamann"
TASKSAMVAR="samvar"
TASKCUTADAPT="cutadapt"
TASKTRIMGALORE="trimgalore"
TASKHICUP="hicup"
TASKHICLIB="hiclib"
TASKFASTQ2SANGER="sanger"
TASKWIGGLER="wiggler"
TASKTRIMMOMATIC="trimmomatic"
TASKHOMERHIC="homerhic"
TASKHOMERCHIPSEQ="homerchipseq"
TASKPEAKRANGER="peakranger"
TASKMACS2="macs2"
TASKMEMECHIP="memechip"

##############################################################
# PROGRAM PATHS
##############################################################
QSUB=prepareJobSubmission.sh
BINQSUB=jobsubmission.sh
QSUBEXTRA=""            # any extra such as email notification

#Additional programs not necessarily available as module
PATH_SAMTOOLS=
PATH_IGVTOOLS=
PATH_PICARD=
PATH_SAMSTAT=
PATH_FASTXTK=

# Commonly used file abbreviations
READONE="read1"
READTWO="read2"
FASTQ="fastq.gz"
FASTA=            # fasta file usually from the reference genome
FASTA_CHROMDIR=   # folder containing individual fasta files for each chromosome of the reference genome 

# file infixes
UNM="unm"   # unmapped
ALN="aln"   # aligned 
MUL="mul"   # non-unique aligned
ASD="asd"   # aligned sorted duplicate-removed
ASR="asdrr" # aligned sorted duplicate-removed raligned recalibrated

##############################################################
# Summary specifics
# html2pdf conversion via PRINCE
MODULE_SUMMARY=
PATH_SUMMARY=
HTMLOUT="Summary"

##############################################################
# gzip alternatives, e.g.
# pigz (2.3) - http://zlib.net/pigz/
MODULE_GZIP=
GZIP="gzip -9"			# command, e.g. gzip or pigz
[ -n "$MODULE_GZIP" ] && module load $MODULE_GZIP

##############################################################
# FASTQC (0.10.1) 
# http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
WALLTIME_FASTQC=10:00:00
MEMORY_FASTQC=20
CPU_FASTQC=16
NODES_FASTQC="nodes=2:ppn=8"

MODULE_FASTQC=
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
# SAMTOOLS (0.1.19) SNP calling
# http://samtools.sourceforge.net/
WALLTIME_SAMVAR=40:00:00
MEMORY_SAMVAR=40
CPU_SAMVAR=1
NODES_SAMVAR="nodes=1:ppn=1"

MODULE_SAMVAR="jdk samtools/0.1.18 igvtools/2.3.5 gatk/2.5"

##############################################################
# BWA (0.7.4) 
# http://bio-bwa.sourceforge.net/
WALLTIME_BWA=50:00:00
MEMORY_BWA=50
CPU_BWA=32
NODES_BWA="nodes=4:ppn=8"

MODULE_BWA=
PATH_BWA=$PATH_IGVTOOLS:$PATH_PICARD:$PATH_SAMSTAT

##############################################################
# Bowtie (1.0.0)
# http://bowtie-bio.sourceforge.net/index.shtml
WALLTIME_BOWTIE=72:00:00
MEMORY_BOWTIE=60
CPU_BOWTIE=8
NODES_BOWTIE="nodes=1:ppn=8"

MODULE_BOWTIE=
PATH_BOWTIE=$PATH_IGVTOOLS:$PATH_PICARD:$PATH_SAMSTAT
BOWTIE_INDEX=


##############################################################
# Bowtie2 (2.1.0) 
# http://bowtie-bio.sourceforge.net/index.shtml

WALLTIME_BOWTIE2=72:00:00
MEMORY_BOWTIE2=60
CPU_BOWTIE2=8
NODES_BOWTIE2="nodes=1:ppn=8"

MODULE_BOWTIE2=
PATH_BOWTIE2=$PATH_IGVTOOLS:$PATH_PICARD:$PATH_SAMSTAT
BOWTIE2_INDEX=

##############################################################
# Wiggler
# https://sites.google.com/site/anshulkundaje/projects/wiggler
WALLTIME_WIGGLER=10:00:00
MEMORY_WIGGLER=60
CPU_WIGGLER=1
NODES_WIGGLER="nodes=1:ppn=1"

MODULE_WIGGLER=
PATH_WIGGLER=

WIGGLER_UMAPDIR=
WIGGLERADDPARAMS=""

##############################################################
# HOMER HIC 
# http://biowhat.ucsd.edu/homer/index.html
WALLTIME_HOMERHIC=60:00:00
MEMORY_HOMERHIC=60
CPU_HOMERHIC=8
NODES_HOMERHIC="nodes=1:ppn=8"

MODULE_HOMERHIC=
PATH_HOMERHIC=$PATH_IGVTOOLS:$PATH_PICARD:$PATH_SAMSTAT

##############################################################
# HOMER CHIPSEQ
# http://biowhat.ucsd.edu/homer/index.html
WALLTIME_HOMERCHIPSEQ=12:00:00
MEMORY_HOMERCHIPSEQ=20
CPU_HOMERCHIPSEQ=1
NODES_HOMERCHIPSEQ="nodes=1:ppn=1"

MODULE_HOMERCHIPSEQ=
PATH_HOMERCHIPSEQ=

##############################################################
# Peakranger
# http://ranger.sourceforge.net/
WALLTIME_PEAKRANGER=12:00:00
MEMORY_PEAKRANGER=20
CPU_PEAKRANGER=1
NODES_PEAKRANGER="nodes=1:ppn=1"

MODULE_PEAKRANGER=
PATH_PEAKRANGER=

##############################################################
# MACS2
# https://github.com/taoliu/MACS/
WALLTIME_MACS2=12:00:00
MEMORY_MACS2=20
CPU_MACS2=1
NODES_MACS2="nodes=1:ppn=1"

MODULE_MACS2=
PATH_MACS2=

##############################################################
# MEMECHIP
# http://meme.nbcr.net/
WALLTIME_MEMECHIP=48:00:00
MEMORY_MEMECHIP=40
CPU_MEMECHIP=8
NODES_MEMECHIP="nodes=1:ppn=8"

MODULE_MEMECHIP=
PATH_MEMECHIP=

MEMECHIPDATABASES=

##############################################################
# Trim adapter with CUTADAPT ()
# https://code.google.com/p/cutadapt/
WALLTIME_CUTADAPT=4:00:00
MEMORY_CUTADAPT=40
CPU_CUTADAPT=1
NODES_CUTADAPT="nodes=1:ppn=1"

MODULE_CUTADAPT=
PATH_CUTADAPT=
CUTADAPT_OPTIONS=

##############################################################
# Trim adapter with TRIMGALORE
# http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/
WALLTIME_TRIMGALORE=4:00:00
MEMORY_TRIMGALORE=40
CPU_TRIMGALORE=1
NODES_TRIMGALORE="nodes=1:ppn=1"

MODULE_TRIMGALORE=
PATH_TRIMGALORE=
TRIMGALORE_ADAPTER1=
TRIMGALORE_ADAPTER2=
TRIMGALORE_OPTIONS=

##############################################################
# Trimming with TRIMMOMATIC
# http://www.usadellab.org/cms/index.php?page=trimmomatic
WALLTIME_TRIMMOMATIC=8:00:00
MEMORY_TRIMMOMATIC=40
CPU_TRIMMOMATIC=1
NODES_TRIMMOMATIC="nodes=1:ppn=1"

MODULE_TRIMMOMATIC=
PATH_TRIMMOMATIC=
TRIMMOMATICADDPARAMS=
TRIMMOMATICSTEPS=

##############################################################
# Snp calling with GATK
# http://www.broadinstitute.org/gatk/
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
# Tophat (2.0.8b) and cufflinks (2.1.1)
# http://tophat.cbcb.umd.edu/
# http://cufflinks.cbcb.umd.edu/
WALLTIME_TOPHAT=60:00:00
MEMORY_TOPHAT=50
CPU_TOPHAT=8
NODES_TOPHAT="nodes=1:ppn=8"

MODULE_TOPHATCUFF=
PATH_TOPHATCUFF=$PATH_IGVTOOLS:$PATH_PICARD:$PATH_SAMSTAT:$PATH_RNASEQC

##############################################################
# HICLIB 
# https://bitbucket.org/mirnylab/hiclib
WALLTIME_HICLIB=50:00:00
MEMORY_HICLIB=60
CPU_HICLIB=16
NODES_HICLIB="nodes=1:ppn=8"
CPU_HICLIB_POSTCOMMAND=1
NODES_HICLIB_POSTCOMMAND="nodes=1:ppn=1"

MODULE_HICLIB=
PATH_HICLIB=
HICLIB_GAPFILE=
HICLIB_RENZYMES=
HICLIB_READLENGTH=

##############################################################
# HICUP + fit-hi-c
# http://www.bioinformatics.babraham.ac.uk/projects/hicup/
WALLTIME_HICUP=10:00:00
MEMORY_HICUP=60
CPU_HICUP=8
NODES_HICUP="nodes=1:ppn=8"

MODULE_HICUP=
PATH_HICUP=
HICUP_RENZYMES=

##############################################################
# R (3.0.0)
# http://www.r-project.org/
WALLTIME_R=1:00:00
MEMORY_R=10
CPU_R=1
NODES_R="nodes=1:ppn=1"

MODULE_R=
PATH_R=
RSCRIPT=Rscript

##############################################################
# Bam Annotations
# 
WALLTIME_BAMANN=5:00:00
MEMORY_BAMANN=32
CPU_BAMANN=1
NODES_BAMANN="nodes=1:ppn=1"

MODULE_BAMANN=
PATH_BAMANN=

##############################################################
# Read re-calibration
# 
WALLTIME_RECAL=60:00:00
MEMORY_RECAL=50
CPU_RECAL=8
NODES_RECAL="nodes=1:ppn=8" 

MODULE_RECAL=
PATH_RECAL=

##############################################################
# reduced representation bisulfite sequencing mapping 
# https://code.google.com/p/bsmap/
WALLTIME_RRBSMAP=60:00:00
MEMORY_RRBSMAP=50
CPU_RRBSMAP=32
NODES_RRBSMAP="nodes=4:ppn=8"

MODULE_RRBSMAP=
PATH_RRBSMAP=


##############################################################
# downsample
# 
WALLTIME_DOWNSAMPLE=5:00:00
MEMORY_DOWNSAMPLE=20
CPU_DOWNSAMPLE=1
NODES_DOWNSAMPLE="nodes=1:ppn=1"

MODULE_DOWNSAMPLE=
PATH_DOWNSAMPLE=$PATH_IGVTOOLS:$PATH_PICARD:$PATH_SAMTOOLS


##############################################################
# demultiplex with Fastxtoolkit
# http://hannonlab.cshl.edu/fastx_toolkit/
WALLTIME_DEMULTIPLEX=5:00:00
MEMORY_DEMULTIPLEX=20
CPU_DEMULTIPLEX=1
NODES_DEMULTIPLEX="nodes=1:ppn=1"

MODULE_DEMULTIPLEX=
PATH_DEMULTIPLEX=$PATH_FASTXTK

##############################################################
#VCFTOOLS="/clusterdata/hiseq_apps/bin/freeze001/VCFtools_0.1.3.2/bin"
#SAMUTILS="/clusterdata/hiseq_apps/bin/freeze001/tabix-0.2.3"
#ANNOVAR="/clusterdata/hiseq_apps/bin/freeze001/annovar"
#
#MACS="/clusterdata/hiseq_apps/bin/devel/MACS_git"
#PEAKFINDER="/clusterdata/hiseq_apps/bin/devel/vancouvershortr_svn/"
#
#VIENNA=
#UNAFOLD=

