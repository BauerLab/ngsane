#!/bin/bash

# Script running indel calling with GATK
# QC:
# author: Denis C. Bauer
# date: Nov.2010

#INPUTS
HISEQINF=$1   # location of the HiSeqInf repository
f=$2          # fastq file
FASTA=$3      # reference genome
DBROD=$4
REFSEQROD=$5
OUT=$6        # output dir
SEQREG=$7     # (optional) region of specific interest, e.g. targeted reseq

#PROGRAMS
. $HISEQINF/conf/header.sh

#PARAMETERS


n=`basename $f`
echo ">>>>> call INDELSs using GATK"
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> gatkIndel.sh $HISEQINF $f $FASTA $DBROD $REFSEQROD $OUT $SEQREG"

# delete old indel files
if [ -e $OUT/${n/bam/indel.vcf} ]; then
    rm $OUT/${n/bam/indel.vcf}
    rm $OUT/${n/bam/indel.bed}
    rm $OUT/${n/bam/indel.vcf}.verbose
    rm $OUT/${n/bam/indel.vcf}.idx
fi


echo "********* call indels"
java -Xmx4g -jar $GATKJAR/GenomeAnalysisTK.jar -l WARN -quiet \
    -T IndelGenotyperV2 \
    -R $FASTA \
    -I $f \
    -o $OUT/${n/bam/indel.vcf} \
    -bed $OUT/${n/bam/indel.bed} \
    -mnr 1000000 \
    -verbose $OUT/${n/bam/indel.vcf}.verbose \
    --refseq $REFSEQROD

#			echo "filter indels"
#			$GATKHOME/perl/filterSingleSampleCalls.pl \
#				--calls $dir/ind/${n/bam/indel.list.bed} \
#				--max_cons_av_mm 3.0 \
#				--max_cons_nqs_av_mm 0.5 \
#				--mode ANNOTATE > $dir/ind/${n/bam/indel.list.fi.bed}


echo "********* make index for IGV"
java -Xmx1g -jar $IGVTOOLS index $OUT/${n/bam/indel.vcf}

echo ">>>>> call INDELs using GATK - FINISHED"
echo ">>>>> enddate "`date`