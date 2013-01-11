#!/bin/bash

# Script to run TopHat program
# It takes comma-seprated list of files containing short sequence reads in fasta or fastq format and bowtie index files as input.
# It produces output files: read alignments in .bam format and other files.
# author: Chikako Ragan, Denis Bauer
# date: Jan. 2011



# messages to look out for -- relevant for the QC.sh script:
# QCVARIABLES,truncated file


echo ">>>>> readmapping with Tophat "
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> tophatcuff.sh $*"


function usage {
echo -e "usage: $(basename $0) -k HISEQINF -f FASTA -r REFERENCE -o OUTDIR [OPTIONS]

Script running read mapping for single and paired DNA reads from fastq files
It expects a fastq file, pairdend, reference genome  as input and 
It runs BWA, converts the output to .bam files, adds header information and
writes the coverage information for IGV.

required:
  -k | --toolkit <path>     location of the HiSeqInf repository 
  -f | --fastq <file>       fastq file
  -r | --reference <file>   reference genome
  -o | --outdir <path>      output dir

options:
  -t | --threads <nr>       number of CPUs to use (default: 1)
  -i | --rgid <name>        read group identifier RD ID (default: exp)
  -l | --rglb <name>        read group library RD LB (default: qbi)
  -p | --rgpl <name>        read group platform RD PL (default: illumna)
  -s | --rgsi <name>        read group sample RG SM prefac (default: )
  -R | --region <ps>        region of specific interest, e.g. targeted reseq
                             format chr:pos-pos
  -S | --sam                do not convert to bam file (default confert); not the
                             resulting sam file is not duplicate removed
  --forceSingle             run single end eventhough second read is present
"
exit
}

if [ ! $# -gt 3 ]; then usage ; fi

#DEFAULTS
THREADS=1
EXPID="exp"           # read group identifier RD ID
LIBRARY="qbi"         # read group library RD LB
PLATFORM="illumina"   # read group platform RD PL
DOBAM=1               # do the bam file
FORCESINGLE=0
INSERT=200
MEMORY=2G

#INPUTS
while [ "$1" != "" ]; do
    case $1 in
#        -k | --toolkit )        shift; HISEQINF=$1 ;; # location of the HiSeqInf repository
	-k | toolkit )          shift; CONFIG=$1 ;; # ENSURE NO VARIABLE NAMES FROM CONFIG
        -t | --threads )        shift; THREADS=$1 ;; # number of CPUs to use
        -f | --fastq )          shift; f=$1 ;; # fastq file
        -r | --reference )      shift; FASTA=$1 ;; # reference genome
        -o | --outdir )         shift; OUTDIR=$1 ;; # output dir
	-a | --annot )          shift; REFSEQGTF=$1 ;; # refseq annotation
	-i | --insert )         shift; INSERT=$1 ;; #mate insert size
	
	-l | --rglb )           shift; LIBRARY=$1 ;; # read group library RD LB
	-p | --rgpl )           shift; PLATFORM=$1 ;; # read group platform RD PL
	-s | --rgsi )           shift; SAMPLEID=$1 ;; # read group sample RG SM (pre)
	-R | --region )         shift; SEQREG=$1 ;; # (optional) region of specific interest, e.g. targeted reseq
        -S | --sam )            DOBAM=0 ;;
	
	--forceSingle )         FORCESINGLE=1;;
        -h | --help )           usage ;;
        * )                     echo "dont understand $1"
    esac
    shift
done


#PROGRAMS (note, both configs are necessary to overwrite the default, here:e.g.  TASKTOPHAT)
. $CONFIG
. $HISEQINF/pbsTemp/header.sh
. $CONFIG

#$BOWTIEDIR=`echo $BOWTIE | sed 's/\(.*\)bowtie/\1/'`

#export bowtie directory
SAMDIR=`dirname $SAMTOOLS`
export PATH=$PATH:$BOWTIETWO:$SAMDIR
module load python
module load jdk/1.7.0_03
module load boost

#module load tophat/2.0.4b
module load $TOPHAT
module load $CUFFLINKS
#module load $BOWTIETWO
#module load samtools

echo $(module list)
echo $(which bowtie2)

JAVAPARAMS="-Xmx"$MEMORY"g -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1 -XX:MaxDirectMemorySize=4G"

echo $PATH
#module load samtools
#module load tophat


n=`basename $f`

#INPUTS
#HISEQINF=$1    # location of the HighSeqInf reposository
#BOWTIEINDEX=$2   # bowtie index files
#SHORTREAD=$3     # a seq read file in fasta or fastq format
#THREADS=$4       # number of CPUs to use
#INNERDIS=$5      # mean inner distance between mate pairs
#OUTDIR=$6           # output dir

CUFOUT=${OUTDIR/$TASKTOPHAT/$TASKCUFF}

#remove old files
if [ -d $OUTDIR ]; then rm -r $OUTDIR; fi
if [ -d $CUFOUT ]; then rm -r $CUFOUT; fi

#is paired ?
if [ -e ${f/$READONE/$READTWO} ] && [ "$FORCESINGLE" = 0 ]; then
    echo "********* PAIRED READS"
    f2=${f/$READONE/$READTWO}
    dmget -a $f
    dmget -a $f2
else
    echo "********* SINGLE READS"
    dmget -a $f
fi

#is ziped ?
ZCAT="cat" # always cat
if [[ $f = *.gz ]]; then
    ZCAT="zcat";
#    echo "unzip first"
#    f=${f/fastq.gz/fastq}
#    if [ ! -e $f ]; then gunzip -c $f.gz >$f; fi
#    if [ -n "$f2" ]; then
#	f2=${f2/fastq.gz/fastq}
#	if [ ! -e $f2 ]; then gunzip -c $f2.gz >$f2; fi
#    fi
#    n=${n/fastq.gz/fastq}
fi

# generating the index files
#if [ ! -e ${FASTA/.fasta/}.1.bt2 ]; then echo ">>>>> make .bt2"; $BOWTIETWO/bowtie2-build $FASTA ${FASTA/.fasta/}; fi
#if [ ! -e $FASTA.fai ]; then echo ">>>>> make .fai"; $SAMTOOLS faidx $FASTA; fi

echo "********* tophat"
tophat -p $THREADS -o $OUTDIR ${FASTA/.fasta/} $f $f2

echo "********* merge mapped and unmapped"
BAMFILE=$OUTDIR/../${n/_$READONE.$FASTQ/.tph.bam}
#ln -f  $OUTDIR/accepted_hits.bam $BAMFILE
$SAMTOOLS merge -f $BAMFILE.tmp.bam $OUTDIR/accepted_hits.bam $OUTDIR/unmapped.bam
$SAMTOOLS sort $BAMFILE.tmp.bam ${BAMFILE/.bam/}
rm $BAMFILE.tmp.bam

##mv $BAMFILE $OUTDIR/../${n/_$READONE.fastq/.tph.bam}.tmp
##$SAMTOOLS sort $OUTDIR/../${n/_$READONE.fastq/.tph.bam}.tmp ${BAMFILE/.bam/}


##statistics
echo "********* flagstat"
$SAMTOOLS flagstat $BAMFILE >$BAMFILE.stats
READ1=`$ZCAT $f | wc -l | gawk '{print int($1/4)}' `
if [ -n "$f2" ]; then READ2=`$ZCAT $f2 | wc -l | gawk '{print int($1/4)}' `; fi
let FASTQREADS=$READ1+$READ2
echo $FASTQREADS" fastq reads" >>$BAMFILE.stats
JUNCTION=`wc -l $OUTDIR/junctions.bed | cut -d' ' -f 1 `
echo $JUNCTION" junction reads" >> $BAMFILE.stats
# get junction genes overlapping exons +-200bp
JUNCTGENE=`$BEDTOOLS/windowBed -a $OUTDIR/junctions.bed -b $REFSEQGTF -u -w 200 | wc -l | cut -d' ' -f 1 `
echo $JUNCTGENE" junction reads NCBIM37" >> $BAMFILE.stats

##index
echo "********* index"
$SAMTOOLS index $BAMFILE

#TODO: the UCSC FASTA is not compattible with the bowtie index
# http://seqanswers.com/forums/showthread.php?t=6478
#echo "********* calculate inner distance"
#if [ ! -e $OUTDIR/metrices ]; then mkdir $OUTDIR/metrices ; fi
#java -Xmx4g -jar $PICARD/CollectMultipleMetrics.jar \
#    INPUT=$BAMFILE \
#    OUTPUT=$OUTDIR/metrices/met \
#    VALIDATION_STRINGENCY=LENIENT \
#    PROGRAM=CollectAlignmentSummaryMetrics \
#    PROGRAM=CollectInsertSizeMetrics \
#    PROGRAM=QualityScoreDistribution \
#    TMP_DIR=$TMP 
#for f in $( ls $OUTDIR/metrices/*.pdf ); do
#    $IMGMAGCONVERT $f ${f/pdf/jpg}
#done

#coverage for IGV
echo "********* coverage track"
java -Xmx1g -jar $IGVTOOLS count $BAMFILE \
    $BAMFILE.cov.tdf ${FASTA/.fasta/}.genome

#run cufflinks
echo "********* cufflinks"
echo ">>>>> from $BAMFILE to $CUFOUT"
echo "$CUFFLINKSHOME/cufflinks --quiet -r $FASTA -p $THREADS -o $CUFOUT $BAMFILE"
cufflinks --quiet --GTF-guide $REFSEQGTF -p $THREADS -o $CUFOUT \
    $BAMFILE 




echo ">>>>> allignment with TopHat - FINISHED"
echo ">>>>> enddate "`date`


