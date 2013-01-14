#!/bin/bash

# Annotation of bam files
# author: Denis C. Bauer
# date: Jan.2013

# messages to look out for -- relevant for the QC.sh script:
# 

echo ">>>>> Annotate BAM file "
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> annotateBam.sh $*"


function usage {
echo -e "usage: $(basename $0) -k CONFIG -f BAM -o OUTDIR [OPTIONS]

Annotating BAM file with annotations in a folder specified in the config
file

required:
  -k | --toolkit <path>     config file
  -f <file>                 bam file

options:

"
exit
}


if [ ! $# -gt 3 ]; then usage ; fi

#DEFAULTS

#INPUTS
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the HiSeqInf repository
        -f           )          shift; f=$1 ;; # bam file
        -h | --help )           usage ;;
        * )                     echo "don't understand "$1
    esac
    shift
done

#PROGRAMS
. $CONFIG
. $HISEQINF/pbsTemp/header.sh
. $CONFIG


n=`basename $f`

# delete old bam file
#if [ -e $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam} ]; then rm $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}; fi
#if [ -e $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}.stats ]; then rm $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}.stats; fi
#if [ -e $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}.dupl ]; then rm $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}.dupl; fi


#RRNA=$DATASTORE/SeqAna/reference/prod/b37/rRNA_b37.gtf
#TRNA=$DATASTORE/SeqAna/reference/prod/b37/tRNA_b37.gtf
#LINCRNA=$DATASTORE/SeqAna/reference/prod/b37/lincRNA_b37.gtf
#MIRNA=$DATASTORE/SeqAna/reference/prod/b37/sno_miRNA_b37.gtf
#SRPRNA=$DATASTORE/SeqAna/reference/prod/b37/srpRNA_b37.gtf
#GENES=$DATASTORE/SeqAna/reference/prod/b37/UCSC_REfSeq_genes_b37.gtf
#SEGDUPS=$DATASTORE/SeqAna/reference/prod/b37/genomicSuperDups_b37.gtf

#Gencode v14
PGENES=$DATASTORE/SeqAna/reference/prod/b37/annotation/gencode.v14.protein_coding.b37.gtf
RRNA=$DATASTORE/SeqAna/reference/prod/b37/annotation/gencode.v14.rRNA.b37.gtf
TRNA=$DATASTORE/SeqAna/reference/prod/b37/annotation/gencode.v14.Mt_tRNA.b37.gtf
LINCRNA=$DATASTORE/SeqAna/reference/prod/b37/annotation/gencode.v14.long_noncoding_RNAs_b37.gtf
MIRNA=$DATASTORE/SeqAna/reference/prod/b37/annotation/gencode.v14.miRNA.b37.gtf
SNORNA=$DATASTORE/SeqAna/reference/prod/b37/annotation/gencode.v14.snoRNA.b37.gtf
SNRNA=$DATASTORE/SeqAna/reference/prod/b37/annotation/gencode.v14.snRNA.b37.gtf
TEC=$DATASTORE/SeqAna/reference/prod/b37/annotation/gencode.v14.TEC.b37.gtf
MISCRNA=$DATASTORE/SeqAna/reference/prod/b37/annotation/gencode.v14.misc_RNA.b37.gtf
POLYA=$DATASTORE/SeqAna/reference/prod/b37/gencode.v14.polyAs_b37.gtf
Other=$DATASTORE/SeqAna/reference/prod/b37/gencode.v14.annotation.b37.gtf

#o=$(head -n ${PBS_ARRAYID} tasks.txt | tail -1)
#echo $o

# excract only the duplicates
# [ if -d $MYOUT/dups ]; then mkdir 
#d=${o/bam/dups.bam}
#d=${d/bwa/"bwa/dups"}
#$SAMTOOLS view -f 0x400 -b $o > ${o/bam/dups.bam}

for o in $f; do

    dmget -a $PGENES $RRNA $TRNA $LINCRNA $MIRNA $SNORNA $SNRNA $TEC $MISCRNA $POLYA $OTHER
    dmget -a $o

    echo $o
    exit

    $BEDTOOLS/bedtools bamtobed -i $o > $o.bed
    $BEDTOOLS/bedtools merge -i $o.bed -n > $o.merg.bed
    rm $o.bed
    $BEDTOOLS/bedtools annotate -counts -i $o.merg.bed -files $PGENES $RRNA $TRNA $LINCRNA $MIRNA $SNORNA $SNRNA $TEC $MISCRNA $POLYA $OTHER -names genes rRNA tRNA lincRNA miRNA snoRNA snrna TEC miscRNA polyA other >$o.merg.anno.bed
    rm $o.merg.bed
    
    echo "total Pgenes rRNA tRNA lincRNA miRNA snoRNA snRNA TEC miscRNA PolyA other rest" >$o.merg.anno.stats
    gawk 'BEGIN{a=0; b=0; c=0; d=0; e=0; f=0; g=0; h=0; i=0; j=0; k=0}{
         a=a+$4; 
         if( $5 !=0){b=b+$4; next}; 
         if( $6 !=0){c=c+$4; next};
         if( $7 !=0){d=d+$4; next}; 
         if( $8 !=0){e=e+$4; next}; 
         if( $9 !=0){f=f+$4; next};
         if( $10!=0){g=g+$4; next};
         if( $11!=0){h=h+$4; next};
         if( $12!=0){i=i+$4; next};
         if( $13!=0){j=j+$4; next};
         if( $14!=0){k=k+$4; next};
         if( $5+$6+$7+$8+$9+$10+$11+$12+$13+$14==0){x=x+$4}}
         END{print a" "b" "c" "d" "e" "f" "g" "h" "i" "j" "k" "x}' $o.merg.anno.bed >> $o.merg.anno.stats

    #echo "total GLxx rRNA tRNA lincRNA sno_miRNA srpRNA REFSEQ segdups rest" >$o.merg.anno.stats
    #gawk 'BEGIN{r=0; t=0; l=0; m=0; s=0; g=0; a=0; u=0; d=0; x=0}{a=a+$4; 
    #     if( $5 !=0){r=r+$4; next}; 
    #     if( $6 !=0){t=t+$4; next};
    #     if( $7 !=0){l=l+$4; next}; 
    #     if( $8 !=0){m=m+$4; next}; 
    #     if( $9 !=0){s=s+$4; next};
    #     if( $10!=0){g=g+$4; next};
    #     if( $11!=0){d=d+$4; next};
    #     if( $1~GL){u=u+$4; next};
    #     if( $11+$5+$6+$7+$8+$9+$10+$11!=0){x=x+$4}}
    #     END{print a" "u" "r" "t" "l" "m" "s" "g" "d" "x}' $o.merg.anno.bed >> $o.merg.anno.stats

    head -n 1 $o.merg.anno.bed >> $o.merg.anno.stats
    cat $o.merg.anno.bed | sort -k4gr | head -n 20 >>$o.merg.anno.stats  
    
done



echo ">>>>> Annotate BAM file"
echo ">>>>> enddate "`date`

