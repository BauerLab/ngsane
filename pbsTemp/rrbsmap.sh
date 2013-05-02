#!/bin/bash

# RRBS mapping
# author: Denis C. Bauer
# date: Sept.2011

# messages to look out for -- relevant for the QC.sh script:
# QCVARIABLES,

echo ">>>>> readmapping with RRBSmap "
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> rrbsmap.sh $*"


function usage {
echo -e "usage: $(basename $0) -k HISEQINF -f FASTQ -r REFERENCE -o OUTDIR [OPTIONS]

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
  -u | --rgpu <name>        read group platform unit RG PU (default:flowcell )
  -A | --adapter <str>      3-prime adapter sequence
"
exit
}


if [ ! $# -gt 3 ]; then usage ; fi

#DEFAULTS
THREADS=1
EXPID="exp"           # read group identifier RD ID
LIBRARY="qbi"         # read group library RD LB
PLATFORM="illumina"   # read group platform RD PL
UNIT="flowcell"       # read group platform unit RG PU
FASTQNAME=""
ADAPTER=""
RSITE="C-CGG"         # restriction enzyme cutting site

#INPUTS
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; HISEQINF=$1 ;; # location of the HiSeqInf repository
        -t | --threads )        shift; THREADS=$1 ;; # number of CPUs to use
        -f | --fastq )          shift; f=$1 ;; # fastq file
        -r | --reference )      shift; FASTA=$1 ;; # reference genome
        -o | --outdir )         shift; OUT=$1 ;; # output dir
	-i | --rgid )           shift; EXPID=$1 ;; # read group identifier RD ID
	-l | --rglb )           shift; LIBRARY=$1 ;; # read group library RD LB
	-p | --rgpl )           shift; PLATFORM=$1 ;; # read group platform RD PL
	-s | --rgsi )           shift; SAMPLEID=$1 ;; # read group sample RG SM (pre)
	-u | --rgpu )           shift; UNIT=$1 ;; # read group platform unit RG PU 
	-A | --adapter )        shift; ADAPTER="-A "$1 ;; # adapter
	-R | --restriction )    shift; RSITE=$1 ;;
	--fastqName )           shift; FASTQNAME=$1 ;; #(name of fastq or fastq.gz)
        -h | --help )           usage ;;
        * )                     usage
    esac
    shift
done


#PROGRAMS
. $HISEQINF/conf/header.sh
if [ -n "$FASTQNAME" ]; then FASTQ=$FASTQNAME ; fi

n=`basename $f`


# delete old bam file
if [ -e $OUT/${n/'_'$READONE.$FASTQ/.rrbs.bam} ]; then rm $OUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}; fi
if [ -e $OUT/${n/'_'$READONE.$FASTQ/.rrbs.bam}.stats ]; then rm $OUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}.stats; fi


#is paired ?
if [ -e ${f/$READONE/$READTWO} ] && [ "$FORCESINGLE" = 0 ]; then
    PAIRED="true"
else
    PAIRED="false"
fi

#is ziped ?
ZCAT="zcat"
if [[ $f != *.fastq.gz ]]; then ZCAT="cat"; fi


FULLSAMPLEID=$SAMPLEID"${n/'_'$READONE.$FASTQ/}"
echo ">>>>> full sample ID "$FULLSAMPLEID

echo "********* mapping"
# Paired read
if [ "$PAIRED" = "true" ]
then
    echo "********* PAIRED READS"
    $RRBSMAP -a $f -b ${f/$READONE/$READTWO} -d $FASTA -o $OUT/${n/$FASTQ/mapped.bam} \
	$ADAPTER -D $RSITE -p $THREADS -2 $OUT/${n/$FASTQ/unmapped.bam}
    READ1=`$ZCAT $f | wc -l | gawk '{print int($1/4)}' `
    READ2=`$ZCAT ${f/$READONE/$READTWO} | wc -l | gawk '{print int($1/4)}' `
    let FASTQREADS=$READ1+$READ2
# Single read
else
    echo "********* SINGLE READS"
    $RRBSMAP -a $f -d $FASTA -o $OUT/${n/$FASTQ/mapped.bam} \
	$ADAPTER -D $RSITE -p $THREADS -2 $OUT/${n/$FASTQ/unmapped.bam}
    let FASTQREADS=`$ZCAT $f | wc -l | gawk '{print int($1/4)}' `
fi


echo "********* merge to single file"
java -Xmx4g -jar $PICARD/MergeBamAlignment.jar \
    UNMAPPED_BAM=$OUT/${n/$FASTQ/unmapped.bam} \
    ALIGNED_BAM=$OUT/${n/$FASTQ/mapped.bam} \
    OUTPUT=$OUT/${n/$FASTQ/merged.bam} \
    REFERENCE_SEQUENCE=$FASTA \
    PAIRED_RUN=$PAIRED \
    IS_BISULFITE_SEQUENCE=true 
    

echo "********* add readgroup"
java -Xmx4g -jar $PICARD/AddOrReplaceReadGroups.jar \
    INPUT=$OUT/${n/$FASTQ/merged.bam} \
    OUTPUT=$OUT/${n/$FASTQ/rrbs.bam} \
    RGID=$EXPID RGLB=$LIBRARY RGPL=$PLATFORM \
    RGPU=$UNIT RGSM=$FULLSAMPLEID 

echo "********* mark duplicates"
if [ ! -e $OUT/metrices ]; then mkdir $OUT/metrices ; fi
THISTMP=$TMP/$n$RANDOM #mk tmp dir because picard writes none-unique files
mkdir $THISTMP
java -Xmx4g -jar $PICARD/MarkDuplicates.jar \
    INPUT=$OUT/${n/$FASTQ/rrbs.bam} \
    OUTPUT=$OUT/${n/$FASTQ/rrbsd.bam} \
    METRICS_FILE=$OUT/metrices/${n/$FASTQ/rrbsd.bam}.dupl AS=true \
    VALIDATION_STRINGENCY=LENIENT \
    TMP_DIR=$THISTMP
rm -r $THISTMP
$SAMTOOLS index $OUT/${n/$FASTQ/rrbsd.bam}


# statistics
echo "********* statistics"
STATSOUT=$OUT/${n/$FASTQ/rrbsd.bam}.stats
$SAMTOOLS flagstat $OUT/${n/$FASTQ/rrbsd.bam} > STATSOUT
fi

echo "********* verify"
BAMREADS=`head -n1 $STATSOUT | cut -d " " -f 1`
if [ "$BAMREADS" = "" ]; then let BAMREADS="0"; fi			
if [ $BAMREADS -eq $FASTQREADS ]; then
    echo "-----------------> PASS check mapping: $BAMREADS == $FASTQREADS"
    rm $OUT/${n/$FASTQ/mapped.bam}
    rm $OUT/${n/$FASTQ/unmapped.bam}
    rm $OUT/${n/$FASTQ/merged.bam}
    rm $OUT/${n/$FASTQ/rrbs.bam}
else
    echo -e "***ERROR**** We are loosing reads from .fastq -> .bam in $f: \nFastq had $FASTQREADS Bam has $BAMREADS"
    exit 1
      
fi

echo "********* coverage track"
GENOME=$(echo $FASTA| sed 's/.fasta/.genome/' | sed 's/.fa/.genome/' )
java -Xmx1g -jar $IGVTOOLS count $OUT/${n/$FASTQ/rrbsd.bam} \
    $OUT/${n/'_'$READONE.$FASTQ/.$ASD.bam.cov.tdf} $GENOME


echo ">>>>> readmapping with rrbsmap - FINISHED"
echo ">>>>> enddate "`date`