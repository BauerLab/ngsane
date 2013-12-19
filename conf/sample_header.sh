##############################################################
# System info
##############################################################
SUBMISSIONSYSTEM="PBS"                            # SGE or PBS
QUEUEWAIT=" -W depend=afterok:"                   # PBS
QUEUEWAITSEP=":"
#QUEUEWAIT=" -hold_jid "                          # SGE
#QUEUEWAITSEP=","
#QUEUEPARENV="smp"     
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
NG_WIGGLER=
NG_HOMER=
NG_CUTADAPT=
NG_TRIMGALORE=
NG_TRIMMOMATIC=
NG_HDF5=
NG_HICLIB=${NG_PYTHON}
NG_HICUP=
NG_FITHIC=
NG_FASTQSCREEN=
NG_MATLAB=
NG_CHANCE=
NG_PARALLEL=
NG_TRINITY=
<<<<<<< HEAD
NG_PINDEL=

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
TASKPINDEL="pindel"
=======
NG_MACS2=${NG_PYTHON}
NG_HTSEQ=${NG_PYTHON}

##############################################################
# Software reference
##############################################################
NG_CITE_NGSANE="(in review); 'NGSANE: A Lightweight Production Informatics Framework for High Throughput Data Analysis; Buske FA, French HJ, Smith MA, Clark SJ, Bauer DC"
NG_CITE_R="R: A language and environment for statistical computing; R Core Team. R Foundation for Statistical Computing, Vienna, Austria, (2013)"
NG_CITE_PYTHON="Guido van Rossum, Jelke de Boer: Linking a Stub Generator (AIL) to a Prototyping Language (Python). Spring 1991 EurOpen Conference Proceedings (May 20-24, 1991) Tromso, Norway."
NG_CITE_FASTQC="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/"
NG_CITE_SAMTOOLS="Bioinformatics. 2009 Aug 15;25(16):2078-9. doi: 10.1093/bioinformatics/btp352. Epub 2009 Jun 8. 'The Sequence Alignment/Map format and SAMtools.'; Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R; 1000 Genome Project Data Processing Subgroup."
NG_CITE_IGVTOOLS="Brief Bioinform. 2013 Mar;14(2):178-92. doi: 10.1093/bib/bbs017. Epub 2012 Apr 19. 'Integrative Genomics Viewer (IGV): high-performance genomics data visualization and exploration.'; Thorvaldsdóttir H, Robinson JT, Mesirov JP."
NG_CITE_GATK="Genome Res. 2010 Sep;20(9):1297-303. doi: 10.1101/gr.107524.110. Epub 2010 Jul 19. 'The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. McKenna A, Hanna M, Banks E, Sivachenko A, Cibulskis K, Kernytsky A, Garimella K, Altshuler D, Gabriel S, Daly M, DePristo MA."
NG_CITE_BWA="Bioinformatics. 2010 Mar 1;26(5):589-95. doi: 10.1093/bioinformatics/btp698. Epub 2010 Jan 15. 'Fast and accurate long-read alignment with Burrows-Wheeler transform.'; Li H, Durbin R."
NG_CITE_PICARD="http://picard.sourceforge.net."
NG_CITE_SAMSTAT="Bioinformatics. 2011 Jan 1;27(1):130-1. doi: 10.1093/bioinformatics/btq614. Epub 2010 Nov 18. 'SAMStat: monitoring biases in next generation sequencing data.'; Lassmann T, Hayashizaki Y, Daub CO."
NG_CITE_UCSCTOOLS="Brief Bioinform. 2013 Mar;14(2):144-61. doi: 10.1093/bib/bbs038. Epub 2012 Aug 20. 'The UCSC genome browser and associated tools.'; Kuhn RM, Haussler D, Kent WJ."
NG_CITE_BEDTOOLS="Bioinformatics. 2010 Mar 15;26(6):841-2. doi: 10.1093/bioinformatics/btq033. Epub 2010 Jan 28. 'BEDTools: a flexible suite of utilities for comparing genomic features.; Quinlan AR, Hall IM."
NG_CITE_BOWTIE="Genome Biol. 2009;10(3):R25. doi: 10.1186/gb-2009-10-3-r25. Epub 2009 Mar 4. 'Ultrafast and memory-efficient alignment of short DNA sequences to the human genome.'; Langmead B, Trapnell C, Pop M, Salzberg SL."
NG_CITE_BOWTIE2="Nat Methods. 2012 Mar 4;9(4):357-9. doi: 10.1038/nmeth.1923. 'Fast gapped-read alignment with Bowtie 2.';Langmead B, Salzberg SL."
NG_CITE_PEAKRANGER="BMC Bioinformatics. 2011 May 9;12:139. doi: 10.1186/1471-2105-12-139. 'PeakRanger: a cloud-enabled peak caller for ChIP-seq data.'; Feng X, Grossman R, Stein L."
NG_CITE_MEME="Nucleic Acids Res. 2009 Jul;37(Web Server issue):W202-8. doi: 10.1093/nar/gkp335. Epub 2009 May 20. 'MEME SUITE: tools for motif discovery and searching.';Bailey TL, Boden M, Buske FA, Frith M, Grant CE, Clementi L, Ren J, Li WW, Noble WS."
NG_CITE_TOPHAT="Bioinformatics. 2009 May 1;25(9):1105-11. doi: 10.1093/bioinformatics/btp120. Epub 2009 Mar 16. 'TopHat: discovering splice junctions with RNA-Seq.'; Trapnell C, Pachter L, Salzberg SL."
NG_CITE_RNASEQC="Bioinformatics. 2012 Jun 1;28(11):1530-2. doi: 10.1093/bioinformatics/bts196. Epub 2012 Apr 25.; 'RNA-SeQC: RNA-seq metrics for quality control and process optimization.'; DeLuca DS, Levin JZ, Sivachenko A, Fennell T, Nazaire MD, Williams C, Reich M, Winckler W, Getz G."
NG_CITE_CUFFLINKS="Nat Biotechnol. 2010 May;28(5):511-5. doi: 10.1038/nbt.1621. Epub 2010 May 2. 'Transcript assembly and quantification by RNA-Seq reveals unannotated transcripts and isoform switching during cell differentiation.'; Trapnell C, Williams BA, Pertea G, Mortazavi A, Kwan G, van Baren MJ, Salzberg SL, Wold BJ, Pachter L"
NG_CITE_BLUE="(in review); 'Blue: correcting sequencing errors using consensus and context'; Greenfield P, Duesing K, Papanicolaou A, Bauer DC"
NG_CITE_WIGGLER="Nature. 2012 Sep 6;489(7414):57-74. doi: 10.1038/nature11247. 'An integrated encyclopedia of DNA elements in the human genome.'; ENCODE Project Consortium, Bernstein BE, Birney E, Dunham I, Green ED, Gunter C, Snyder M."
NG_CITE_HOMER="Mol Cell. 2010 May 28;38(4):576-89. doi: 10.1016/j.molcel.2010.05.004. 'Simple combinations of lineage-determining transcription factors prime cis-regulatory elements required for macrophage and B cell identities.'; Heinz S, Benner C, Spann N, Bertolino E, Lin YC, Laslo P, Cheng JX, Murre C, Singh H, Glass CK."
NG_CITE_CUTADAPT="MARTIN, M.; 'Cutadapt removes adapter sequences from high-throughput sequencing reads.'; EMBnet.journal, North America, 17, may. 2011. Available at: <http://journal.embnet.org/index.php/embnetjournal/article/view/200>. Date accessed: 18 Nov. 2013."
NG_CITE_TRIMGALORE="http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/"
NG_CITE_TRIMMOMATIC="Nucleic Acids Res. 2012 Jul;40(Web Server issue):W622-7. doi: 10.1093/nar/gks540. Epub 2012 Jun 8. 'RobiNA: a user-friendly, integrated software solution for RNA-Seq-based transcriptomics.'; Lohse M, Bolger AM, Nagel A, Fernie AR, Lunn JE, Stitt M, Usadel B."
NG_CITE_HICLIB="Nat Methods. 2012 Oct;9(10):999-1003. doi: 10.1038/nmeth.2148. Epub 2012 Sep 2. 'Iterative correction of Hi-C data reveals hallmarks of chromosome organization.'; Imakaev M, Fudenberg G, McCord RP, Naumova N, Goloborodko A, Lajoie BR, Dekker J, Mirny LA."
NG_CITE_HICUP="http://www.bioinformatics.babraham.ac.uk/projects/hicup/"
NG_CITE_FITHIC=
NG_CITE_FASTQSCREEN="http://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/"
NG_CITE_MATLAB="MathWorks, (2012). Bioinformatics Toolbox: User's Guide (R2012a). Retrieved July 14, 2012 from www.mathworks.com/help/pdf_doc/bioinfo/bioinfo_ug.pdf"
NG_CITE_CHANCE="Genome Biol. 2012 Oct 15;13(10):R98. [Epub ahead of print] 'CHANCE: comprehensive software for quality control and validation of ChIP-seq data.'; Diaz A, Nellore A, Song JS."
NG_CITE_PARALLEL="http://www.gnu.org/s/parallel‎"
NG_CITE_TRINITY="Nat Biotechnol. 2011 May 15;29(7):644-52. doi: 10.1038/nbt.1883. 'Full-length transcriptome assembly from RNA-Seq data without a reference genome.'; Grabherr MG, Haas BJ, Yassour M, Levin JZ, Thompson DA, Amit I, Adiconis X, Fan L, Raychowdhury R, Zeng Q, Chen Z, Mauceli E, Hacohen N, Gnirke A, Rhind N, di Palma F, Birren BW, Nusbaum C, Lindblad-Toh K, Friedman N, Regev A."
NG_CITE_MACS2="Genome Biol. 2008;9(9):R137. doi: 10.1186/gb-2008-9-9-r137. Epub 2008 Sep 17. 'Model-based analysis of ChIP-Seq (MACS).'; Zhang Y, Liu T, Meyer CA, Eeckhoute J, Johnson DS, Bernstein BE, Nusbaum C, Myers RM, Brown M, Li W, Liu XS."
NG_CITE_HTSEQ="Conf Proc IEEE Eng Med Biol Soc. 2013 Jul;2013:647-50. doi: 10.1109/EMBC.2013.6609583. 'Benchmarking RNA-Seq quantification tools.'; Chandramohan R, Wu PY, Phan JH, Wang MD."
>>>>>>> origin/master

##############################################################
# PROGRAM PATHS
##############################################################
QSUB=${NGSANE_BASE}/core/prepareJobSubmission.sh
BINQSUB=${NGSANE_BASE}/core/jobSubmission.sh
QSUBEXTRA=""            # any extra such as email notification

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
# gzip alternatives, e.g.
# pigz (2.3) - http://zlib.net/pigz/
MODULE_GZIP=${NG_GZIP}
GZIP="gzip -9"			# command, e.g. gzip or pigz
[ -n "$MODULE_GZIP" ] && module load $MODULE_GZIP

# source content in default folder twice to properly set up crosslinked variables
for NGSANE_DEFAULTS in ${NGSANE_BASE}/conf/header.d/* ; do source $NGSANE_DEFAULTS; done
for NGSANE_DEFAULTS in ${NGSANE_BASE}/conf/header.d/* ; do source $NGSANE_DEFAULTS; done

# Location of the meme motif database of choice
MEMECHIPDATABASES=
<<<<<<< HEAD

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
# RNA-Seq De novo Assembly Using Trinity
# http://trinityrnaseq.sourceforge.net/

### Stage P1: Time and resources required for Inchworm stage
### Only use at maximum, half the available CPUs on a node 
# - Inchworm will not efficiently use any more than 4 CPUs and you will have to take longer for resources to be assigned
# —min_kmer_cov 2 to reduce memory requirements with large read sets.
WALLTIME_INCHWORM="4:00:00"		# optional on Wolfpack 
MEMORY_INCHWORM="40" 			# will use it for --JM
NCPU_INCHWORM="4" 				# Use less than half of the CPUs on a node. This algorithm is limited by cache memory
NODES_INCHWORM="1"
NODETYPE_INCHWORM="all.q"  		
#NODETYPE_INCHWORM="intel.q" 	# Inchworm performs faster when Trinity was installed using the Intell compiler (Intell systems only

### Stage P2: Time and resources required for Chrysalis stage
### Starts with Bowtie alignment and post-processing of alignment file
### All CPUs presenct can be used for the Chrysalis parts. 
#They may take a while to be provisioned, so the less request, possibly the faster the jobs turnaround.
# For one step (the parallel sort) it needs as much memory as specified in P1. Less memory, means more I/O for sorting
WALLTIME_CHRYSALIS="24:00:00"		# optional on Wolfpack 
MEMORY_CHRYSALIS="40"	 			# will use it for --JM
NCPU_CHRYSALIS="16" 				# For very large datasets, besides normalisation, maybe use 32 cores
NODES_CHRYSALIS="1"
NODETYPE_CHRYSALIS="all.q"  		# dont use intel.q on Wolfpack for this

# This stage is actually Chrysalis::readsToTranscript and Butterfly. Both should ideally be run through a SGE/PBS array 
# The Chrysalis bit is I/O heavy, so a local memory node is used. If files take up over 500GB, this will cause problems. 
# You may want to normalise your data and/or run Martin's optimised, standalone Trinity module
WALLTIME_BUTTERFLY="72:00:00"		
MEMORY_BUTTERFLY="40"	 			
NCPU_BUTTERFLY="32" 				
NODES_BUTTERFLY="1"
NODETYPE_BUTTERFLY="all.q"  		

MODULES_TRINITY=
PATH_TRINITY=

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
MODULE_BLUE="${NG_MONO} ${NG_BLUE} ${NG_R} ${NG_IMAGEMAGIC}"
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

##############################################################
# pindel
# http://gmt.genome.wustl.edu/pindel/current/
WALLTIME_PINDEL=60:00:00
MEMORY_PINDEL=50
CPU_PINDEL=32
NODES_PINDEL="nodes=4:ppn=8"
INPUT_PINDEL=$TASKBWA
MODULE_PINDEL="${NG_PINDEL} ${NG_PERL}"
PATH_PINDEL=
=======
>>>>>>> origin/master
