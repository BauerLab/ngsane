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
NG_MACS2=${NG_PYTHON}
NG_HTSEQ=${NG_PYTHON}
NG_CIRCOS=
NG_BLAT=
NG_PINDEL=
NG_SEQLOGO=

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
NG_CITE_PARALLEL="http://www.gnu.org/s/parallel"
NG_CITE_TRINITY="Nat Biotechnol. 2011 May 15;29(7):644-52. doi: 10.1038/nbt.1883. 'Full-length transcriptome assembly from RNA-Seq data without a reference genome.'; Grabherr MG, Haas BJ, Yassour M, Levin JZ, Thompson DA, Amit I, Adiconis X, Fan L, Raychowdhury R, Zeng Q, Chen Z, Mauceli E, Hacohen N, Gnirke A, Rhind N, di Palma F, Birren BW, Nusbaum C, Lindblad-Toh K, Friedman N, Regev A."
NG_CITE_MACS2="Genome Biol. 2008;9(9):R137. doi: 10.1186/gb-2008-9-9-r137. Epub 2008 Sep 17. 'Model-based analysis of ChIP-Seq (MACS).'; Zhang Y, Liu T, Meyer CA, Eeckhoute J, Johnson DS, Bernstein BE, Nusbaum C, Myers RM, Brown M, Li W, Liu XS."
NG_CITE_HTSEQ="Conf Proc IEEE Eng Med Biol Soc. 2013 Jul;2013:647-50. doi: 10.1109/EMBC.2013.6609583. 'Benchmarking RNA-Seq quantification tools.'; Chandramohan R, Wu PY, Phan JH, Wang MD."
NG_CITE_CIRCOS="Genome Res. 2009 Sep;19(9):1639-45. doi: 10.1101/gr.092759.109. Epub 2009 Jun 18. 'Circos: an information aesthetic for comparative genomics.'; Krzywinski M, Schein J, Birol I, Connors J, Gascoyne R, Horsman D, Jones SJ, Marra MA."
NG_CITE_BLAT="Genome Res. 2002 Apr;12(4):656-64. 'BLAT--the BLAST-like alignment tool.'; Kent WJ."
NG_CITE_PINDEL="Ye K, Schulz MH, Long Q, Apweiler R, Ning Z. Pindel: a pattern growth approach to detect break points of large deletions and medium sized insertions from paired-end short reads. Bioinformatics. 2009 Nov 1;25(21):2865-71. Epub 2009 Jun 26."
NG_CITE_SEQLOGO="Nucleic Acids Res. 1990 Oct 25;18(20):6097-100.; 'Sequence logos: a new way to display consensus sequences.'; Schneider TD, Stephens RM."

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
