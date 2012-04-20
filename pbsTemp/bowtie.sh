#!/bin/bash

#maybe newline before bash?

# Script to run bowtie program.
# It takes comma-seprated list of files containing short sequence reads in fasta or fastq format and bowtie index files as input.
# It produces output files: read alignments in .bam format and other files.
# author: Chikako Ragan, Denis Bauer
# date: Jan. 2011

#INPUTS
HISEQINF=$1      # location of the HighSeqInf reposository
BOWTIEINDEX=$2   # bowtie index files
SHORTREAD=$3     # location of read sequences
MISMATCH=$4      # number of max mismatch allowed
TRIMEDSEQ=$5     # num of nt trimed from 3' end
OUT=$6           # output dir
THREADS=$7

#PROGRAMS
. $HISEQINF/pbsTemp/header.sh


echo ">>>>> alignment with with bowtie "
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> bowtie.sh $*"

name=`basename $SHORTREAD`
name=${name/.fastq/}".mm"

#remove old files
if [ -e $OUT/$name.bam ]; then rm -r $OUT/$name.bam; fi

#run bowtie command -v $MISMATCH -m 1
echo "********* bowtie" 
$BOWTIE -S -a -p $THREADS -3 $TRIMEDSEQ $BOWTIEINDEX $SHORTREAD $OUT/$name.sam

echo "********* sorting, bam-conversion"
$SAMTOOLS view -bt $FASTA.fai $OUT/$name.sam | $SAMTOOLS sort - $OUT/$name.ash

#index
echo "********* index"
$SAMTOOLS index $OUT/$name.ash.bam

#coverage for IGV
echo "********* coverage track"
java -Xmx1g -jar $IGVTOOLS count $OUT/$name.ash.bam \
    $OUT/$name.ash.bam.tdf $BOWTIEINDEX.genome

# calculate statistics
echo "********* statistics"
$SAMTOOLS flagstat $OUT/$name.ash.bam > $OUT/$name.ash.bam.stats

# test that all the reads in the fastq file are in the bam file
# and clean up
let FASTQREADS=`wc -l $SHORTREAD | gawk '{print $i/4}' `
BAMREADS=`head -n1 $OUT/$name.ash.bam.stats | cut -d " " -f 1`
if [ "$BAMREADS" = "" ]; then let BAMREADS="0"; fi			
if [ $BAMREADS -eq $FASTQREADS ]; then
    echo "-----------------> PASS check mapping: $BAMREADS == $FASTQREADS"
    rm $OUT/$name.sam
else
    echo -e "***ERROR**** We are loosing reads from .fastq -> .bam in $SHORTREAD: \nFastq had $FASTQREADS Bam has $BAMREADS"
    exit 1
      
fi

echo ">>>>> allignment with bowtie - FINISHED"
echo ">>>>> enddate "`date`

