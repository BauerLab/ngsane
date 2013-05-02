#!/bin/bash

# Call variants with samtools
# author: Denis C. Bauer
# date: Feb.2013

# messages to look out for -- relevant for the QC.sh script:
# 

echo ">>>>> Variant calling with sam "
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> samSNPs.sh $*"


function usage {
echo -e "usage: $(basename $0) -k CONFIG -f BAM -o OUTDIR [OPTIONS]

Variant calling with sam

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
        -f | --file    )        shift; f=$1 ;; # bam file
	-o | --output  )        shift; MYOUT=$1 ;; # output directory
        -h | --help    )        usage ;;
        * )                     echo "don't understand "$1
    esac
    shift
done

#PROGRAMS
. $CONFIG
. $HISEQINF/conf/header.sh
. $CONFIG


n=`basename $f`

module load samtools

# delete old bam file
#if [ -e $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam} ]; then rm $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}; fi
#if [ -e $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}.stats ]; then rm $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}.stats; fi
#if [ -e $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}.dupl ]; then rm $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}.dupl; fi

echo "********* rm duplicate reads $(date)"
samtools rmdup $f $MYOUT/${n/bam/drm.bam}
samtools index $MYOUT/${n/bam/drm.bam}

echo "********* call variants $(date)"
samtools mpileup -uf $FASTA -q1 -D $MYOUT/${n/bam/drm.bam} |  bcftools view -vcg - >$MYOUT/${n/bam/vcf}

echo "********* convert bcf->vcf; index $(date)"
vcfutils.pl varFilter -D1000 -w0 -e0 $MYOUT/${n/bam/vcf}  >$MYOUT/${n/bam/clean.vcf}

echo "********* index vcf file for viewing in IGV $(date)"
java -Xmx2g -jar $IGVTOOLS index $MYOUT/${n/bam/clean.vcf}


rm $MYOUT/${n/bam/drm.bam} $MYOUT/${n/bam/drm.bam}.bai
rm $MYOUT/${n/bam/vcf}


echo ">>>>> Variant calling with sam"
echo ">>>>> enddate "`date`

