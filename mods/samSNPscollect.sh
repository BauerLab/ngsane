#!/bin/bash

# Call variants with samtools
# author: Denis C. Bauer
# date: Feb.2013

# messages to look out for -- relevant for the QC.sh script:
# 

echo ">>>>> Collect Variants after calling with sam "
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> samSNPscollect.sh $*"


function usage {
echo -e "usage: $(basename $0) -k CONFIG -f BAM -o OUTDIR [OPTIONS]

Variant calling with sam

required:
  -k | --toolkit <path>     config file
  -f <files>                 bam files
  -o <dir>
options:

"
exit
}


if [ ! $# -gt 3 ]; then usage ; fi

#DEFAULTS

#INPUTS
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the NGSANE repository
        -f           )          shift; FILES=${1//,/ } ;; # bam files
	    -o           )          shift; MYOUT=$1 ;; # outputdir
        -h | --help )           usage ;;
        * )                     echo "don't understand "$1
    esac
    shift
done

#PROGRAMS
. $CONFIG
. ${NGSANE_BASE}/conf/header.sh
. $CONFIG


#module load samtools
module load jdk

# delete old bam file
#if [ -e $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam} ]; then rm $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}; fi
#if [ -e $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}.stats ]; then rm $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}.stats; fi
#if [ -e $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}.dupl ]; then rm $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}.dupl; fi

# ensure dir is there
if [ ! -d $MYOUT ]; then mkdir -p $MYOUT; fi


# prep for joining
VARIANTS=""
NAMES=""
for f in $FILES; do
    i=${f/$TASKBWA/$TASKBWA"-"$TASKSAMVAR} #point to var folder
    i=${i/bam/"clean.vcf"} # correct ending
    b=$(basename $i)
    arrIN=(${b//./ })
    name=${arrIN[0]}
    NAMES=$NAMES"$name,"
    VARIANTS=$VARIANTS" --variant:$name $i "
done

echo $VARIANTS

REGION=""
if [ -n "$REF" ]; then echo $REF; REGION="-L $REF"; fi


echo "********** join with GATK"
java -Xmx2g -jar $GATKJAR/GenomeAnalysisTK.jar \
   -R $FASTA \
   -T CombineVariants \
   $VARIANTS \
   -o $MYOUT/joined.vcf \
   -genotypeMergeOptions PRIORITIZE \
   $REGION \
   -priority $NAMES

echo ">>>>> Join variant after calling with sam"
echo ">>>>> enddate "`date`

