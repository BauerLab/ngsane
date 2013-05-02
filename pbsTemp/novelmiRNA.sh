#!/bin/bash

# Find novel miRNA by finding genomic regions that do not overlap with known
# annotations (exons, miRNA, ncRNA) but still have a good coverage in RNA seq
# data. 
# author: Denis C. Bauer
# date: Feb.2011

#INPUTS
HISEQINF=$1   # location of the HiSeqInf repository
f=$2          # bam file
FASTA=$3      # reference genome
MIRNA=$4      # all annotated areas in the genome
GENES=$5    # 
CONS=$6
OUT=$7        # output dir

#PROGRAMS
. $HISEQINF/conf/header.sh


n=`basename $f`
name=${n/.bam/}
GENOME=$(echo $FASTA| sed 's/.fasta/.BEDgenome/' | sed 's/.fa/.BEDgenome/' )
echo ">>>>> identify novel RNAs "
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> novelmiRNA.sh $HISEQINF $f $FASTA $MIRNA $GENES $OUT"

transcrFilter="1"
consFilter="1"
miRBaseFilter=""


# delete old files
<<EOF

echo "********* bed conversion"
# extract only mapped, (paired) and high quality reads
if [ "$PAIRED" = "1" ]; then echo "PAIRED"; SPEC="-f 0x2"; fi
$SAMTOOLS view -q 10 $SPEC -b $f > $OUT/$name.clean.bam


# get the coverage of these area
$BEDTOOLS/genomeCoverageBed -ibam $OUT/$name.clean.bam -g $GENOME -bga > $OUT/$name.clean.bedg

# get all with sufficient depth
awk '$4 >= 10' $OUT/$name.clean.bedg | $BEDTOOLS/mergeBed -i stdin > $OUT/$name.gd10.bed


echo ">>>>> Found "`wc -l $OUT/$name.gd10.bed | cut -d " " -f1 `" locations with 10x coverage ( miRBase "`$BEDTOOLS/intersectBed -a $OUT/$name.gd10.bed -b $MIRNA -u | wc -l`")"

echo "********* intersect"
# find locations that are overlapping between a and b 
# and print the original region from b ensure the right strand
# genome  ====================================================
# ANN       ==========      ========            ==========
# ANNFREE ==          ======        ============          ====     
# RNAseq           =====      ==         ====
# inters           =====                 ====

FI=""
if [ -n "$transcrFilter" ]; then
    $BEDTOOLS/subtractBed -a $OUT/$name.gd10.bed -b $GENES -s > $OUT/$name.gd10$FI"T".bed
    FI="T"
    echo ">>>>> "`wc -l $OUT/$name.gd10$FI.bed | cut -d " " -f1 `" pass exon filter ( miRBase "`$BEDTOOLS/intersectBed -a $OUT/$name.gd10$FI.bed -b $MIRNA -u | wc -l`")"
fi

if [ -n "$consFilter" ]; then
    $BEDTOOLS/intersectBed -a $OUT/$name.gd10$FI.bed -b $CONS -u -s > $OUT/$name.gd10$FI"C".bed
    FI=$FI"C"
    echo ">>>>> "`wc -l $OUT/$name.gd10$FI.bed | cut -d " " -f1 `" pass conservation filter ( miRBase "`$BEDTOOLS/intersectBed -a $OUT/$name.gd10$FI.bed -b $MIRNA -u | wc -l`")"
fi

if [ -n "$miRBaseFilter" ]; then
    # this is not correct...
    $BEDTOOLS/intersectBed -a $OUT/$name.gd10$FI.bed -b $MIRNA -v > $OUT/$name.gd10$FI"R".bed
    FI=$FI"R"
    echo ">>>>> "`wc -l $OUT/$name.gd10$FI.bed | cut -d " " -f1 `" are not in miRBase"
fi

echo "********* extract sequence"
gawk '{OFS="\t"; print $1,$2-60,$3+60,$4}' $OUT/$name.gd10$FI.bed > $OUT/$name.gd10.novel$FI.ext.bed
$BEDTOOLS/fastaFromBed -s -fi $FASTA -bed $OUT/$name.gd10.novel$FI.ext.bed -fo $OUT/$name.gd10.novel$FI.fasta

echo "********* fold miRNA"
if [ ! -e $OUT/pics ]; then mkdir $OUT/pics; else rm -r $OUT/pics/*; fi
cd $OUT/pics
$VIENNA/RNAfold < $OUT/$name.gd10.novel$FI.fasta >$OUT/$name.gd10.novel$FI.fold
for i in $( ls ); do
    mv $i ${i/:/_}
done
cd ../../

EOF

 echo ">>>>> "`wc -l $OUT/$name.gd10.novelTCS.bed | cut -d " " -f1 `" pass conservation filter ( miRBase "`$BEDTOOLS/intersectBed -a $OUT/$name.gd10.novelTCS.bed -b $MIRNA -u | wc -l`")"


echo ">>>>> identify novel RNAs - FINISHED"
echo ">>>>> enddate "`date`