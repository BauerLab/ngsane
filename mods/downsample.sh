#!/bin/bash

# Script running downsample
# QC:
# author: Denis C. Bauer
# date: Sept.2011

echo ">>>>> Downsample"
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"

function usage {
echo -e "usage: $(basename $0)

required:
  -k | --toolkit <path>     location of the NGSANE repository 
  -i | --input <file>       bam file
  -o | --outdir <path>      output dir

"
exit
}


#if [ ! $# -gt 2 ]; then usage ; fi

#DEFAULTS

#INPUTS
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the NGSANE repository
        -i | --input )          shift; f=$1 ;; # bam file
        -o | --outdir )         shift; OUT=$1 ;; # output dir
        -r | --reference )      shift; FASTA=$1 ;; # reference genome
        -s | --downsample )     shift; READNUMBER=$1 ;; #readnumber
        -h | --help )           usage ;;
        * )                     usage
    esac
    shift
done



#PROGRAMS
. $CONFIG
. ${NGSANE_BASE}/conf/header.sh
. $CONFIG

# get basename of f
n=${f##*/}

#-F 0x0400
echo "********* extract properly paired none duplicate"
$SAMTOOLS view -f 0x0002 -h -b $f >$OUT/${n/bam/pn.bam}
echo "********* get number of properly paired none duplicate"
READS=`$SAMTOOLS view -f 0x0002 -c $f`
echo $READS

echo "********* downsample"
PROB=`echo "$READNUMBER/$READS" | bc -l`
echo $PROB
java -jar -Xmx4g $PICARD/DownsampleSam.jar \
    INPUT=$OUT/${n/bam/pn.bam} \
    OUTPUT=$OUT/${n/bam/pns.bam} \
    RANDOM_SEED=1 \
    VALIDATION_STRINGENCY=LENIENT \
    PROBABILITY=$PROB

$SAMTOOLS index $OUT/${n/bam/pns.bam}
    
# statistics
echo "********* statistics"
$SAMTOOLS flagstat $OUT/${n/bam/pns.bam} > $OUT/${n/bam/pns.bam}.stats


echo "********* coverage track"
GENOME=$(echo $FASTA| sed 's/.fasta/.genome/' | sed 's/.fa/.genome/' )
java -Xmx1g -jar $IGVTOOLS count $OUT/${n/bam/pns.bam} \
    $OUT/${n/bam/pns.bam.cov.tdf} $GENOME


rm $OUT/${n/bam/pn.bam}

echo ">>>>> readmapping with Downsample - FINISHED"
echo ">>>>> enddate "`date`
  