#!/bin/bash

# Modified bwa.sh - produce output in sam file
# Script running read mapping for single and paired DNA reads from fastq files
# It expects a fastq file, pairdend, reference genome  as input and 
# It runs BWA, converts the output to .bam files, adds header information and
# writes the coverage information for IGV.
# QC:
# it tests that the .bam file has the same number of reads as the fastq file.
# author: Denis C. Bauer
# date: Nov.2010

#INPUTS
HISEQINF=$1   # location of the HiSeqInf repository
THREADS=$2    # number of CPUs to use
f=$3          # fastq file
FASTA=$4      # reference genome
OUT=$5        # output dir
EXPID=$6      # read group identifier RD ID
LIBRARY=$7    # read group library RD LB
PLATFORM=$8   # read group platform RD PL
SAMPLEID=$9   # read group sample RG SM (pre)
SEQREG=$10    # (optional) region of specific interest, e.g. targeted reseq

#PROGRAMS
. $HISEQINF/pbsTemp/header.sh


n=`basename $f`
echo ">>>>> readmapping with BWA "
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> bwa.sh $HISEQINF $THREADS $f $FASTA $OUT $EXPID $LIBRARY $PLATFORM $SAMPLEID $SEQREG"


# delete old bam file
if [ -e $OUT/${n/_read1.fastq/.asd.bam} ]; then rm $OUT/${n/_read1.fastq/.asd.bam}; fi
if [ -e $OUT/${n/_read1.fastq/.asd.bam}.stats ]; then rm $OUT/${n/_read1.fastq/.asd.bam}.stats; fi
if [ -e $OUT/${n/_read1.fastq/.asd.bam}.dupl ]; then rm $OUT/${n/_read1.fastq/.asd.bam}.dupl; fi


#is paired ?
if [ -e ${f/read1/read2} ]; then
    PAIRED="1"
else
    PAIRED="0"
fi

FULLSAMPLEID=$SAMPLEID"${n/_read1.fastq/}"
echo ">>>>> full sample ID "$FULLSAMPLEID

# generating the index files
if [ ! -e $FASTA.bwt ]; then $BWA index -a bwtsw $FASTA; fi
if [ ! -e $FASTA.fai ]; then $SAMTOOLS faidx $FASTA; fi

echo "********* mapping"
# Paired read
if [ "$PAIRED" = 1 ]
then 	
    echo "********* PAIRED"
    $BWA aln -t $THREADS $FASTA $f > $OUT/${n/fastq/sai}
    $BWA aln -t $THREADS $FASTA ${f/read1/read2} > $OUT/${n/read1.fastq/read2.sai}
    $BWA sampe $FASTA $OUT/${n/fastq/sai} $OUT/${n/read1.fastq/read2.sai} \
	-i $EXPID \
	-m $FULLSAMPLEID \
	-p $PLATFORM \
	-l $LIBRARY \
	$f ${f/read1/read2} >$OUT/${n/_read1.fastq/.aln.sam}
    
    rm $OUT/${n/fastq/sai}
    rm $OUT/${n/read1.fastq/read2.sai}
    READ1=`wc -l $f | gawk '{print $i/4}' `
    READ2=`wc -l ${f/read1/read2} | gawk '{print $i/4}' `
    let FASTQREADS=$READ1+$READ2
# Single read
else
    echo "********* SINGLE"
    $BWA aln -t $THREADS $FASTA $f > $OUT/${n/fastq/sai}
    $BWA samse $FASTA $OUT/${n/fastq/sai} \
	-i $EXPID \
	-m $FULLSAMPLEID \
	-p $PLATFORM \
	-l $LIBRARY \
	$f >$OUT/${n/_read1.fastq/.aln.sam}
    rm $OUT/${n/fastq/sai}
    let FASTQREADS=`wc -l $f | gawk '{print $i/4}' `
fi

echo "********* sorting, bam-conversion"
$SAMTOOLS view -bt $FASTA.fai $OUT/${n/_read1.fastq/.aln.sam} | $SAMTOOLS sort - $OUT/${n/_read1.fastq/.ash}

#TODO look at samtools for rmdup
#val string had to be set to LENIENT to avoid crash due to a definition dis-
#agreement between bwa and picard
#http://seqanswers.com/forums/showthread.php?t=4246
echo "********* mark duplicates"
java -Xmx4g -jar $PICARD/MarkDuplicates.jar \
    INPUT=$OUT/${n/_read1.fastq/.ash.bam} \
    OUTPUT=$OUT/${n/_read1.fastq/.asd.bam} \
    METRICS_FILE=$OUT/${n/_read1.fastq/.asd.bam}.dupl AS=true \
    VALIDATION_STRINGENCY=LENIENT \
    TMP_DIR=$TMP

$SAMTOOLS index $OUT/${n/_read1.fastq/.asd.bam}


# statistics
echo "********* statistics"
STATSOUT=$OUT/${n/_read1.fastq/.asd.bam}.stats
$SAMTOOLS flagstat $OUT/${n/_read1.fastq/.asd.bam} > $STATSOUT
if [ -n $SEQREG ]; then
    echo "#custom for "$SEQREG >> $STATSOUT
    echo `$SAMTOOLS view $OUT/${n/_read1.fastq/.ash.bam} $SEQREG | wc -l`" total reads in region " >> $STATSOUT
    echo `$SAMTOOLS view -f 2 $OUT/${n/_read1.fastq/.ash.bam} $SEQREG | wc -l`" properly paired reads in region " >> $STATSOUT
fi


BAMREADS=`head -n1 $OUT/${n/_read1.fastq/.asd.bam}.stats | cut -d " " -f 1`
if [ "$BAMREADS" = "" ]; then let BAMREADS="0"; fi			
if [ $BAMREADS -eq $FASTQREADS ]; then
    echo "-----------------> PASS check mapping: $BAMREADS == $FASTQREADS"
    #rm $OUT/${n/_read1.fastq/.aln.sam}
    rm $OUT/${n/_read1.fastq/.ash.bam}
else
    echo -e "***ERROR**** We are loosing reads from .fastq -> .bam in $f: \nFastq had $FASTQREADS Bam has $BAMREADS"
    exit 1
      
fi

echo "********* coverage track"
GENOME=$(echo $FASTA| sed 's/.fasta/.genome/' | sed 's/.fa/.genome/' )
java -Xmx1g -jar $IGVTOOLS count $OUT/${n/_read1.fastq/.asd.bam} \
    $OUT/${n/_read1.fastq/.asd.bam.cov.tdf} $GENOME


echo ">>>>> readmapping with BWA - FINISHED"
echo ">>>>> enddate "`date`
