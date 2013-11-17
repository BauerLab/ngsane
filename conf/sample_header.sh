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
## uncomment within CSIRO
#module use /apps/gi/modulefiles

##############################################################
# Software Modules
##############################################################
NG_R=
NG_PYTHON=
NG_GZIP=
NG_JAVA=
NG_FASTQC=
NG_SAMTOOLS=
NG_IGVTOOLS=
NG_GATK=
NG_BWA=
NG_IMAGEMAGIC=
NG_PICARD=
NG_SAMSTAT=
NG_UCSCTOOLS=
NG_BEDTOOLS=
NG_BOWTIE=
NG_BOWTIE2=
NG_BOOST=
NG_PEAKRANGER=
NG_MEME=
NG_TOPHAT=
NG_RNASEQC=
NG_CUFFLINKS=
NG_MONO=
NG_BLUE=
NG_PERL=
NG_PRINCE=
NG_WIGGLER=
NG_HOMER=
NG_CUTADAPT=
NG_TRIMGALORE=
NG_TRIMMOMATIC=
NG_HDF5=
NG_HICLIB=
NG_HICUP=
NG_FITHIC=
NG_FASTQSCREEN=
NG_MATLAB=
NG_CHANCE=
NG_PARALLEL=
NG_TRINITY=

##############################################################
# Task Names
##############################################################
TASKFASTQC="fastQC"
TASKBWA="bwa"
TASKBOWTIE="bowtie"
TASKBOWTIE2="bowtie2"
TASKRCA="reCalAln"
TASKVAR="variant"
TASKSNP="snp"
TASKIND="indel"
TASKDOWN="downsample"
TASKDOC="coverage"
TASKDIFFEXP="diffexp"
TASKTOPHAT="tophat"
TASKHTSEQCOUNT="htseqcount"
TASKCUFFLINKS="cufflinks"
TASKCUFFDIFF="cuffdiff"
TASKRRBSMAP="rrbs"
TASKANNOVAR="annovar"
TASKBAMANN="bamann"
TASKSAMVAR="samvar"
TASKCUTADAPT="cutadapt"
TASKTRIMGALORE="trimgalore"
TASKHICUP="hicup"
TASKHICLIB="hiclib"
TASKWIGGLER="wiggler"
TASKTRIMMOMATIC="trimmomatic"
TASKHOMERHIC="homerhic"
TASKHOMERCHIPSEQ="homerchipseq"
TASKPEAKRANGER="peakranger"
TASKMACS2="macs2"
TASKMEMECHIP="memechip"
TASKTRINITY="trinity"
TASKINCHWORM="trinity_inchworm"
TASKCHRYSALIS="trinity_chrysalis"
TASKBUTTERFLY="trinity_butterfly"
TASKFASTQSCREEN="fastqscreen"
TASKBIGWIG="bigwig"
TASKBLUE="blue"
TASKCHANCE="chance"
TASKPOOLBAMS="pooledbam"

##############################################################
# PROGRAM PATHS
##############################################################
QSUB=${NGSANE_BASE}/core/prepareJobSubmission.sh
BINQSUB=${NGSANE_BASE}/core/jobSubmission.sh
QSUBEXTRA=""            # any extra such as email notification

#Additional programs not necessarily available as module
PATH_SAMTOOLS=
PATH_IGVTOOLS=
PATH_PICARD=
PATH_SAMSTAT=
PATH_FASTXTK=

# Commonly used file abbreviations
READONE="_read1"
READTWO="_read2"
FASTQ="fastq.gz"
FASTA=            # fasta file usually from the reference genome
FASTA_CHROMDIR=   # folder containing individual fasta files for each chromosome of the reference genome 

# file infixes
UNM="unm"   # unmapped
ALN="aln"   # aligned 
MUL="mul"   # non-unique aligned
ASD="asd"   # aligned sorted duplicate-removed
ASR="asdrr" # aligned sorted duplicate-removed raligned recalibrated

MODULES_DEFAULT=
for MODULE in $MODULES_DEFAULT; do module load $MODULES_DEFAULT; done

##############################################################
# Summary specifics
# html2pdf conversion via PRINCE
MODULE_SUMMARY="${NG_R} ${NG_PYTHON} ${NG_PRINCE}"
PATH_SUMMARY=
HTMLOUT="Summary"

##############################################################
# gzip alternatives, e.g.
# pigz (2.3) - http://zlib.net/pigz/
MODULE_GZIP=${NG_GZIP}
GZIP="gzip -9"			# command, e.g. gzip or pigz
[ -n "$MODULE_GZIP" ] && module load $MODULE_GZIP

##############################################################
# FASTQC (0.10.1) 
# http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
WALLTIME_FASTQC=10:00:00
MEMORY_FASTQC=10
CPU_FASTQC=2
NODES_FASTQC="nodes=1:ppn=2"
INPUT_FASTQC="fastq"
MODULE_FASTQC="${NG_JAVA} ${NG_FASTQC}"
PATH_FASTQC=

##############################################################
# SAMTOOLS (0.1.19) SNP calling
# http://samtools.sourceforge.net/
WALLTIME_SAMVAR=40:00:00
MEMORY_SAMVAR=40
CPU_SAMVAR=1
NODES_SAMVAR="nodes=1:ppn=1"
INPUT_SAMVAR=$TASKBWA
MODULE_SAMVAR="${NG_SAMTOOLS} ${NG_JAVA} ${NG_IGVTOOLS} ${NG_GATK}"
PATH_SAMVAR=

##############################################################
# BWA (0.7.4) 
# http://bio-bwa.sourceforge.net/
WALLTIME_BWA=50:00:00
MEMORY_BWA=50
CPU_BWA=32
NODES_BWA="nodes=4:ppn=8"
INPUT_BWA="fastq"
MODULE_BWA="${NG_BWA} ${NG_JAVA} ${NG_SAMTOOLS} ${NG_IGVTOOLS} ${NG_R} ${NG_IMAGEMAGIC} ${NG_PICARD} ${NG_SAMSTAT} ${NG_UCSCTOOLS} ${NG_BEDTOOLS}"
PATH_BWA=

##############################################################
# Bowtie (1.0.0)
# http://bowtie-bio.sourceforge.net/index.shtml
WALLTIME_BOWTIE=72:00:00
MEMORY_BOWTIE=60
CPU_BOWTIE=8
NODES_BOWTIE="nodes=1:ppn=8"
INPUT_BOWTIE="fastq"
MODULE_BOWTIE="${NG_BOWTIE} ${NG_JAVA} ${NG_SAMTOOLS} ${NG_IGVTOOLS} ${NG_R} ${NG_IMAGEMAGIC} ${NG_PICARD} ${NG_SAMSTAT} ${NG_UCSCTOOLS} ${NG_BEDTOOLS}"
PATH_BOWTIE=$PATH_IGVTOOLS:$PATH_PICARD:$PATH_SAMSTAT

##############################################################
# Bowtie2 (2.1.0) 
# http://bowtie-bio.sourceforge.net/index.shtml
WALLTIME_BOWTIE2=72:00:00
MEMORY_BOWTIE2=60
CPU_BOWTIE2=8
NODES_BOWTIE2="nodes=1:ppn=8"
INPUT_BOWTIE2="fastq"
MODULE_BOWTIE2="${NG_BOWTIE2} ${NG_JAVA} ${NG_SAMTOOLS} ${NG_IGVTOOLS} ${NG_R} ${NG_IMAGEMAGIC} ${NG_PICARD} ${NG_SAMSTAT} ${NG_UCSCTOOLS} ${NG_BEDTOOLS}"
PATH_BOWTIE2=$PATH_IGVTOOLS:$PATH_PICARD:$PATH_SAMSTAT

##############################################################
# Wiggler
# https://sites.google.com/site/anshulkundaje/projects/wiggler
WALLTIME_WIGGLER=10:00:00
MEMORY_WIGGLER=60
CPU_WIGGLER=1
NODES_WIGGLER="nodes=1:ppn=1"
INPUT_WIGGLER=$TASKBWA
MODULE_WIGGLER="${NG_SAMTOOLS} ${NG_MATLAB} ${NG_WIGGLER}"
PATH_WIGGLER=

WIGGLER_UMAPDIR=

##############################################################
# HOMER HIC 
# http://biowhat.ucsd.edu/homer/index.html
WALLTIME_HOMERHIC=60:00:00
MEMORY_HOMERHIC=60
CPU_HOMERHIC=8
NODES_HOMERHIC="nodes=1:ppn=8"
INPUT_HOMERHIC=$TASKBWA
MODULE_HOMERHIC="${NG_PERL} ${NG_HOMER}"
PATH_HOMERHIC=$PATH_IGVTOOLS:$PATH_PICARD:$PATH_SAMSTAT

##############################################################
# HOMER CHIPSEQ
# http://biowhat.ucsd.edu/homer/index.html
WALLTIME_HOMERCHIPSEQ=12:00:00
MEMORY_HOMERCHIPSEQ=20
CPU_HOMERCHIPSEQ=1
NODES_HOMERCHIPSEQ="nodes=1:ppn=1"
INPUT_HOMERCHIPSEQ=$TASKBOWTIE
MODULE_HOMERCHIPSEQ="${NG_HOMER} ${NG_JAVA} ${NG_R} ${NG_SAMTOOLS} ${NG_PERL}"
PATH_HOMERCHIPSEQ=

##############################################################
# Peakranger
# http://ranger.sourceforge.net/
WALLTIME_PEAKRANGER=12:00:00
MEMORY_PEAKRANGER=20
CPU_PEAKRANGER=1
NODES_PEAKRANGER="nodes=1:ppn=1"
INPUT_PEAKRANGER=$TASKBOWTIE
MODULE_PEAKRANGER="${NG_BOOST} ${NG_R} ${NG_PEAKRANGER}"
PATH_PEAKRANGER=

##############################################################
# MACS2
# https://github.com/taoliu/MACS/
WALLTIME_MACS2=12:00:00
MEMORY_MACS2=20
CPU_MACS2=1
NODES_MACS2="nodes=1:ppn=1"
INPUT_MACS2=$TASKBOWTIE
MODULE_MACS2="${NG_PYTHON} ${NG_R} ${NG_UCSCTOOLS}"
PATH_MACS2=

##############################################################
# MEMECHIP
# http://meme.nbcr.net/
WALLTIME_MEMECHIP=48:00:00
MEMORY_MEMECHIP=40
CPU_MEMECHIP=8
NODES_MEMECHIP="nodes=1:ppn=8"
INPUT_MEMECHIP=$TASKMACS2
MODULE_MEMECHIP="${NG_MEME} ${NG_BEDTOOLS} ${NG_PERL}"
PATH_MEMECHIP=

MEMECHIPDATABASES=

##############################################################
# Trim adapter with CUTADAPT ()
# https://code.google.com/p/cutadapt/
WALLTIME_CUTADAPT=4:00:00
MEMORY_CUTADAPT=40
CPU_CUTADAPT=1
NODES_CUTADAPT="nodes=1:ppn=1"
INPUT_CUTADAPT="fastq"
MODULE_CUTADAPT="${NG_CUTADAPT}"
PATH_CUTADAPT=

##############################################################
# Trim adapter with TRIMGALORE
# http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/
WALLTIME_TRIMGALORE=4:00:00
MEMORY_TRIMGALORE=40
CPU_TRIMGALORE=1
NODES_TRIMGALORE="nodes=1:ppn=1"
INPUT_TRIMGALORE="fastq"
MODULE_TRIMGALORE="${NG_TRIMGALORE} ${NG_CUTADAPT}"
PATH_TRIMGALORE=

##############################################################
# Trimming with TRIMMOMATIC
# http://www.usadellab.org/cms/index.php?page=trimmomatic
WALLTIME_TRIMMOMATIC=8:00:00
MEMORY_TRIMMOMATIC=40
CPU_TRIMMOMATIC=1
NODES_TRIMMOMATIC="nodes=1:ppn=1"
INPUT_TRIMMOMATIC="fastq"
MODULE_TRIMMOMATIC="${NG_JAVA} ${NG_TRIMMOMATIC}"
PATH_TRIMMOMATIC=

##############################################################
# Snp calling with GATK
# http://www.broadinstitute.org/gatk/
# COVERAGE
WALLTIME_GATKDOC=50:00:00
MEMORY_GATKDOC=50
CPU_GATKDOC=1
NODES_GATKDOC="nodes=1:ppn=1"
INPUT_GATKDOC=$TASKRCA
# GATK VARCALL
WALLTIME_VAR=100:00:00
MEMORY_VAR=20
CPU_VAR=1
NODES_VAR="nodes=1:ppn=1"
INPUT_VAR=$TASKRCA

MODULE_GATK="${NG_GATK} ${NG_JAVA} ${NG_R} ${NG_SAMTOOLS} ${NG_IGVTOOLS}"
MODULE_GATKSNP="${NG_GATK} ${NG_JAVA} ${NG_R} ${NG_IGVTOOLS}"

##############################################################
# Tophat (2.0.9) 
# http://tophat.cbcb.umd.edu/
WALLTIME_TOPHAT=62:00:00
MEMORY_TOPHAT=50
CPU_TOPHAT=16
NODES_TOPHAT="nodes=2:ppn=8"
INPUT_TOPHAT="fastq"
MODULE_TOPHAT="${NG_TOPHAT} ${NG_BOOST} ${NG_JAVA} ${NG_PYTHON} ${NG_R} ${NG_BOWTIE2} ${NG_SAMTOOLS} ${NG_IMAGEMAGIC} ${NG_IGVTOOLS} ${NG_PICARD} ${NG_SAMSTAT} ${NG_BEDTOOLS} ${NG_RNASEQC}"
PATH_TOPHAT=$PATH_IGVTOOLS:$PATH_PICARD:$PATH_SAMSTAT:$PATH_RNASEQC

##############################################################
# Cufflinks (2.1.1)
# http://cufflinks.cbcb.umd.edu/
WALLTIME_CUFFLINKS=192:00:00
MEMORY_CUFFLINKS=50
CPU_CUFFLINKS=4
NODES_CUFFLINKS="nodes=1:ppn=4"
INPUT_CUFFLINKS=$TASKTOPHAT
MODULE_CUFFLINKS="${NG_CUFFLINKS}"
PATH_CUFFLINKS=

##############################################################
# HTSEQ-count (0.5.4.p3)
# http://www-huber.embl.de/users/anders/HTSeq/doc/index.html
WALLTIME_HTSEQCOUNT=24:00:00
MEMORY_HTSEQCOUNT=50
CPU_HTSEQCOUNT=1
NODES_HTSEQCOUNT="nodes=1:ppn=1"
INPUT_HTSEQCOUNT=$TASKTOPHAT
MODULE_HTSEQCOUNT="${NG_PYTHON} ${NG_R} ${NG_BEDTOOLS} ${NG_SAMTOOLS}"
PATH_HTSEQCOUNT=

##############################################################
# HICLIB 
# https://bitbucket.org/mirnylab/hiclib
WALLTIME_HICLIB=50:00:00
MEMORY_HICLIB=60
CPU_HICLIB=16
NODES_HICLIB="nodes=1:ppn=8"
CPU_HICLIB_POSTCOMMAND=1
NODES_HICLIB_POSTCOMMAND="nodes=1:ppn=1"
INPUT_HICLIB="fastq"
MODULE_HICLIB="${NG_HICLIB} ${NG_BOWTIE2} ${NG_SAMTOOLS} ${NG_HDF5} ${NG_PICARD}"
PATH_HICLIB=
HICLIB_GAPFILE=

##############################################################
# HICUP + fit-hi-c
# http://www.bioinformatics.babraham.ac.uk/projects/hicup/
WALLTIME_HICUP=10:00:00
MEMORY_HICUP=60
CPU_HICUP=8
NODES_HICUP="nodes=1:ppn=8"
INPUT_HICUP="fastq"
MODULE_HICUP="${NG_HICUP} ${NG_PYTHON} ${NG_FITHIC}"
PATH_HICUP=

##############################################################
# Bam Annotations
# 
WALLTIME_BAMANN=5:00:00
MEMORY_BAMANN=32
CPU_BAMANN=1
NODES_BAMANN="nodes=1:ppn=1"
INPUT_BAMANN=$TASKBWA
MODULE_BAMANN="${NG_BEDTOOLS}"
PATH_BAMANN=

##############################################################
# Read re-calibration
# 
WALLTIME_RECAL=60:00:00
MEMORY_RECAL=50
CPU_RECAL=8
NODES_RECAL="nodes=1:ppn=8" 
INPUT_REALRECAL=$TASKBWA
MODULE_RECAL="${NG_JAVA} ${NG_GATK} ${NG_R} ${NG_SAMTOOLS} ${NG_IGVTOOLS}"
PATH_RECAL=

##############################################################
# reduced representation bisulfite sequencing mapping 
# https://code.google.com/p/bsmap/
WALLTIME_RRBSMAP=60:00:00
MEMORY_RRBSMAP=50
CPU_RRBSMAP=32
NODES_RRBSMAP="nodes=4:ppn=8"
INPUT_RRBSMAP="fastq"
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
# Screen reads against multiple indices
# http://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/
WALLTIME_FASTQSCREEN=48:00:00
MEMORY_FASTQSCREEN=60
CPU_FASTQSCREEN=8
NODES_FASTQSCREEN="nodes=1:ppn=8"
INPUT_FASTQSCREEN="fastq"
MODULE_FASTQSCREEN="${NG_PERL} ${NG_FASTQSCREEN} ${NG_BOWTIE2}"
PATH_FASTQSCREEN=

FASTQSCREEN_DBCONF=

##############################################################
# Create bigwigs from bam files
# http://genome.ucsc.edu/util.html/
WALLTIME_BIGWIG=12:00:00
MEMORY_BIGWIG=12
CPU_BIGWIG=2
NODES_BIGWIG="nodes=1:ppn=2"
INPUT_BIGWIG=$TASKBOWTIE
MODULE_BIGWIG="${NG_UCSCTOOLS} ${NG_JAVA} ${NG_SAMTOOLS} ${NG_BEDTOOLS}"
PATH_BIGWIG=

##############################################################
# Blue Read error correction
# http://www.bioinformatics.csiro.au/blue/
WALLTIME_BLUE=10:00:00
MEMORY_BLUE=60
CPU_BLUE=4
NODES_BLUE="nodes=1:ppn=4"
INPUT_BLUE="fastq"
MODULE_BLUE="${NG_MONO} ${NG_BLUE}"
PATH_BLUE=

##############################################################
# ChIP QC with CHANCE
# https://github.com/songlab/chance/downloads
WALLTIME_CHANCE=4:00:00
MEMORY_CHANCE=20
CPU_CHANCE=1
NODES_CHANCE="nodes=1:ppn=1"
INPUT_CHANCE=$TASKBOWTIE
MODULE_CHANCE="${NG_CHANCE} ${NG_JAVA} ${NG_MATLAB} ${NG_R}"
PATH_CHANCE=

##############################################################
# Pool bam files (e.g. replicates)
# 
WALLTIME_POOLBAMS=10:00:00
MEMORY_POOLBAMS=60
CPU_POOLBAMS=16
NODES_POOLBAMS="nodes=2:ppn=8"
INPUT_POOLBAMS=$TASKBOWTIE
MODULE_POOLBAMS="${NG_PARALLEL} ${NG_PICARD} ${NG_SAMTOOLS} ${NG_IGVTOOLS} ${NG_SAMSTAT}"
PATH_POOLBAMS=$PATH_IGVTOOLS:$PATH_PICARD:$PATH_SAMSTAT

##############################################################
# RNA-Seq De novo Assembly Using Trinity
# http://trinityrnaseq.sourceforge.net/

### Stage P1: Time and resources required for Inchworm stage
### Only use at maximum, half the available CPUs on a node
# - Inchworm will not efficiently use any more than 4 CPUs and you will have to take longer for resources to be assigned
# —min_kmer_cov 2 to reduce memory requirements with large read sets.
WALLTIME_INCHWORM="4:00:00"             # optional on Wolfpack
MEMORY_INCHWORM="40"                    # will use it for --JM
NCPU_INCHWORM="4"                               # Use less than half of the CPUs on a node. This algorithm is limited by cache memory
NODES_INCHWORM="1"
NODETYPE_INCHWORM="all.q"
INPUT_INCHWORM="fastq"
#NODETYPE_INCHWORM="intel.q"    # Inchworm performs faster when Trinity was installed using the Intell compiler (Intell systems only

### Stage P2: Time and resources required for Chrysalis stage
### Starts with Bowtie alignment and post-processing of alignment file
### All CPUs presenct can be used for the Chrysalis parts.
#They may take a while to be provisioned, so the less request, possibly the faster the jobs turnaround.
# For one step (the parallel sort) it needs as much memory as specified in P1. Less memory, means more I/O for sorting
WALLTIME_CHRYSALIS="24:00:00"           # optional on Wolfpack
MEMORY_CHRYSALIS="40"                           # will use it for --JM
NCPU_CHRYSALIS="16"                             # For very large datasets, besides normalisation, maybe use 32 cores
NODES_CHRYSALIS="1"
NODETYPE_CHRYSALIS="all.q"              # dont use intel.q on Wolfpack for this
INPUT_CHRYSALIS="fastq"

# This stage is actually Chrysalis::readsToTranscript and Butterfly. Both should ideally be run through a SGE/PBS array
# The Chrysalis bit is I/O heavy, so a local memory node is used. If files take up over 500GB, this will cause problems.
# You may want to normalise your data and/or run Martin's optimised, standalone Trinity module
WALLTIME_BUTTERFLY="72:00:00"
MEMORY_BUTTERFLY="40"
NCPU_BUTTERFLY="32"
NODES_BUTTERFLY="1"
NODETYPE_BUTTERFLY="all.q"
INPUT_BUTTERFLY="fastq"

MODULES_TRINITY="${NG_TRINITY} ${NG_BOWTIE} ${NG_JAVA}"
PATH_TRINITY=

