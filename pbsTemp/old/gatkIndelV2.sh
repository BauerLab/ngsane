#!/bin/bash

# Script running indel calling with GATK
# QC:
# author: Denis C. Bauer
# date: Nov.2010

#INPUTS
HISEQINF=$1   # location of the HiSeqInf repository
FILES=$2          # file with paths to bamfiles
FASTA=$3      # reference genome
DBROD=$4
REFSEQROD=$5
OUT=$6        # output dir
SEQREG=$7     # (optional) region of specific interest, e.g. targeted reseq

#PROGRAMS
. $HISEQINF/pbsTemp/header.sh

#PARAMETERS


echo ">>>>> call INDELSs using GATK"
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> gatkIndelV2.sh $HISEQINF $FILES $FASTA $DBROD $REFSEQROD $OUT $SEQREG"

# delete old indel files
if [ -e $OUT/genotype.indel.vcf ]; then
    rm $OUT/genotype.indel.vcf
    rm $OUT/genotype.indel.bed
    rm $OUT/genotype.indel.vcf.verbose
    rm $OUT/genotype.indel.vcf.idx
fi

echo "********* call indels"
echo "java -Xmx16g -jar $GATKJAR/GenomeAnalysisTK.jar -l WARN -quiet \
    -T IndelGenotyperV2 \
    -R $FASTA \
    -o $OUT/genotype.indel.vcf \
    -bed $OUT/genotype.indel.bed \
    -mnr 1000000 \
    -L $SEQREG \
    -verbose $OUT/genotype.indel.vcf.verbose \\" >indel.tmp
for f in $( less $FILES ); do echo "-I $f \\" >>indel.tmp ; done
echo "--refseq $REFSEQROD" >>indel.tmp

chmod -u=rwx indel.tmp
./indel.tmp

echo "********* make index for IGV"
java -Xmx1g -jar $IGVTOOLS index $OUT/genotype.indel.vcf

rm indel.tmp
rm $FILES

echo ">>>>> call INDELs using GATK - FINISHED"
echo ">>>>> enddate "`date`