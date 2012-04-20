#!/bin/bash

# BWA calling script
# author: Denis C. Bauer
# date: Nov.2010

# messages to look out for -- relevant for the QC.sh script:
# QCVARIABLES,We are loosing reads,for unmapped read,no such file,file not found,bwa.sh: line

echo ">>>>> readmapping with BWA "
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> bwa.sh $*"


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
  -R | --region <ps>        region of specific interest, e.g. targeted reseq
                             format chr:pos-pos
  -S | --sam                do not convert to bam file (default confert); not the
                             resulting sam file is not duplicate removed
  --forceSingle             run single end eventhough second read is present
  --oldIllumina             old Illumina encoding 1.3+
  --noMapping
  --fastqName               name of fastq file ending (fastq.gz)
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
DOBAM=1               # do the bam file
FORCESINGLE=0
NOMAPPING=0
FASTQNAME=""
QUAL="" # standard Sanger

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
	-R | --region )         shift; SEQREG=$1 ;; # (optional) region of specific interest, e.g. targeted reseq
        -S | --sam )            DOBAM=0 ;;
	--fastqName )           shift; FASTQNAME=$1 ;; #(name of fastq or fastq.gz)
	--forceSingle )         FORCESINGLE=1;;
	--noMapping )           NOMAPPING=1;;
	--oldIllumina )         QUAL="-S";;   # old illumina encoding 1.3+
        -h | --help )           usage ;;
        * )                     usage
    esac
    shift
done


#HISEQINF=$1   # location of the HiSeqInf repository
#THREADS=$2    # number of CPUs to use
#f=$3          # fastq file
#FASTA=$4      # reference genome
#OUT=$5        # output dir
#EXPID=$6      # read group identifier RD ID
#LIBRARY=$7    # read group library RD LB
#PLATFORM=$8   # read group platform RD PL
#SAMPLEID=$9   # read group sample RG SM (pre)
#SEQREG=$10    # (optional) region of specific interest, e.g. targeted reseq

#PROGRAMS
. $HISEQINF/pbsTemp/header.sh
if [ -n "$FASTQNAME" ]; then FASTQ=$FASTQNAME ; fi

n=`basename $f`

# delete old bam file
if [ -e $OUT/${n/'_'$READONE.$FASTQ/.$ASD.bam} ]; then rm $OUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}; fi
if [ -e $OUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}.stats ]; then rm $OUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}.stats; fi
if [ -e $OUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}.dupl ]; then rm $OUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}.dupl; fi

#is paired ?
if [ -e ${f/$READONE/$READTWO} ] && [ "$FORCESINGLE" = 0 ]; then
    PAIRED="1"
else
    PAIRED="0"
fi

#is ziped ?
ZCAT="zcat"
if [[ $f != *.fastq.gz ]]; then ZCAT="cat"; fi

FULLSAMPLEID=$SAMPLEID"${n/'_'$READONE.$FASTQ/}"
echo ">>>>> full sample ID "$FULLSAMPLEID

# generating the index files
if [ ! -e $FASTA.bwt ]; then echo ">>>>> make .bwt"; $BWA index -a bwtsw $FASTA; fi
if [ ! -e $FASTA.fai ]; then echo ">>>>> make .fai"; $SAMTOOLS faidx $FASTA; fi

echo "********* mapping"
# Paired read
if [ "$PAIRED" = 1 ]
then
    if [ "$NOMAPPING" = 0 ]; then
    echo "********* PAIRED READS"
    $BWA aln $QUAL -t $THREADS $FASTA $f > $OUT/${n/$FASTQ/sai}
    $BWA aln $QUAL -t $THREADS $FASTA ${f/$READONE/$READTWO} > $OUT/${n/$READONE.$FASTQ/$READTWO.sai}

    #TODO INCLUDE
    #-r "@RG\tID:$EXPID\tSM:$FULLSAMPLEID\tPL:$PLATFORM\tLB:$LIBRARY\tPU:$UNIT" \
    $BWA sampe $FASTA $OUT/${n/$FASTQ/sai} $OUT/${n/$READONE.$FASTQ/$READTWO.sai} \
	-r "@RG\tID:$EXPID\tSM:$FULLSAMPLEID\tPL:$PLATFORM\tLB:$LIBRARY" \
	$f ${f/$READONE/$READTWO} >$OUT/${n/'_'$READONE.$FASTQ/.$ALN.sam}

#   $BWA sampe $FASTA $OUT/${n/$FASTQ/sai} $OUT/${n/$READONE.$FASTQ/$READTWO.sai} \
#	-i $EXPID \
#	-m $FULLSAMPLEID \
#	-p $PLATFORM \
#	-l $LIBRARY \
#	$f ${f/$READONE/$READTWO} >$OUT/${n/'_'$READONE.$FASTQ/.$ALN.sam}
    
    rm $OUT/${n/$FASTQ/sai}
    rm $OUT/${n/$READONE.$FASTQ/$READTWO.sai}
    fi
    READ1=`$ZCAT $f | wc -l | gawk '{print int($1/4)}' `
    READ2=`$ZCAT ${f/$READONE/$READTWO} | wc -l | gawk '{print int($1/4)}' `
    let FASTQREADS=$READ1+$READ2
# Single read
else
    echo "********* SINGLE READS"
    $BWA aln $QUAL -t $THREADS $FASTA $f > $OUT/${n/$FASTQ/sai}

    $BWA samse $FASTA $OUT/${n/$FASTQ/sai} \
	-r "@RG\tID:$EXPID\tSM:$FULLSAMPLEID\tPL:$PLATFORM\tLB:$LIBRARY" \
	$f >$OUT/${n/'_'$READONE.$FASTQ/.$ALN.sam}

#    $BWA samse $FASTA $OUT/${n/$FASTQ/sai} \
#	-i $EXPID \
#	-m $FULLSAMPLEID \
#	-p $PLATFORM \
#	-l $LIBRARY \
#	$f >$OUT/${n/'_'$READONE.$FASTQ/.$ALN.sam}

    rm $OUT/${n/$FASTQ/sai}
    let FASTQREADS=`$ZCAT $f | wc -l | gawk '{print int($1/4)}' `
fi

# exit if only the sam file is required
if [ "$DOBAM" = 0 ]; then
    SAMREADS=`grep -v @ $OUT/${n/'_'$READONE.$FASTQ/.$ALN.sam} | wc -l`
    if [ "$SAMREADS" = "" ]; then let SAMREADS="0"; fi			
    if [ $SAMREADS -eq $FASTQREADS ]; then
	echo "-----------------> PASS check mapping: $SAMREADS == $FASTQREADS"
	echo ">>>>> readmapping with BWA - FINISHED"
	echo ">>>>> enddate "`date`
	exit
    else
	echo -e "***ERROR**** We are loosing reads from .fastq -> .sam in $f: 
                 Fastq had $FASTQREADS Bam has $SAMREADS"
	exit 1
    fi
fi


# continue for normal bam file conversion
echo "********* sorting and bam-conversion"
$SAMTOOLS view -bt $FASTA.fai $OUT/${n/'_'$READONE.$FASTQ/.$ALN.sam} | $SAMTOOLS sort - $OUT/${n/'_'$READONE.$FASTQ/.ash}


#TODO look at samtools for rmdup
#val string had to be set to LENIENT to avoid crash due to a definition dis-
#agreement between bwa and picard
#http://seqanswers.com/forums/showthread.php?t=4246
echo "********* mark duplicates"
if [ ! -e $OUT/metrices ]; then mkdir $OUT/metrices ; fi
THISTMP=$TMP/$n$RANDOM #mk tmp dir because picard writes none-unique files
mkdir $THISTMP
java -Xmx4g -jar $PICARD/MarkDuplicates.jar \
    INPUT=$OUT/${n/'_'$READONE.$FASTQ/.ash.bam} \
    OUTPUT=$OUT/${n/'_'$READONE.$FASTQ/.$ASD.bam} \
    METRICS_FILE=$OUT/metrices/${n/'_'$READONE.$FASTQ/.$ASD.bam}.dupl AS=true \
    VALIDATION_STRINGENCY=LENIENT \
    TMP_DIR=$THISTMP
rm -r $THISTMP
$SAMTOOLS index $OUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}



# statistics
echo "********* statistics"
STATSOUT=$OUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}.stats
$SAMTOOLS flagstat $OUT/${n/'_'$READONE.$FASTQ/.$ASD.bam} > $STATSOUT
if [ -n $SEQREG ]; then
    echo "#custom region" >> $STATSOUT
    echo `$SAMTOOLS view $OUT/${n/'_'$READONE.$FASTQ/.ash.bam} $SEQREG | wc -l`" total reads in region " >> $STATSOUT
    echo `$SAMTOOLS view -f 2 $OUT/${n/'_'$READONE.$FASTQ/.ash.bam} $SEQREG | wc -l`" properly paired reads in region " >> $STATSOUT
fi


echo "********* calculate inner distance"
export PATH=$PATH:/usr/bin/
THISTMP=$TMP/$n$RANDOM #mk tmp dir because picard writes none-unique files
mkdir $THISTMP
java -Xmx4g -jar $PICARD/CollectMultipleMetrics.jar \
    INPUT=$OUT/${n/'_'$READONE.$FASTQ/.$ASD.bam} \
    REFERENCE_SEQUENCE=$FASTA \
    OUTPUT=$OUT/metrices/${n/'_'$READONE.$FASTQ/.$ASD.bam} \
    VALIDATION_STRINGENCY=LENIENT \
    PROGRAM=CollectAlignmentSummaryMetrics \
    PROGRAM=CollectInsertSizeMetrics \
    PROGRAM=QualityScoreDistribution \
    TMP_DIR=$THISTMP
for im in $( ls $OUT/metrices/*.pdf ); do
    $IMGMAGCONVERT $im ${im/pdf/jpg}
done
rm -r $THISTMP

echo "********* verify"
BAMREADS=`head -n1 $OUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}.stats | cut -d " " -f 1`
if [ "$BAMREADS" = "" ]; then let BAMREADS="0"; fi			
if [ $BAMREADS -eq $FASTQREADS ]; then
    echo "-----------------> PASS check mapping: $BAMREADS == $FASTQREADS"
    rm $OUT/${n/'_'$READONE.$FASTQ/.$ALN.sam}
    rm $OUT/${n/'_'$READONE.$FASTQ/.ash.bam}
else
    echo -e "***ERROR**** We are loosing reads from .fastq -> .bam in $f: \nFastq had $FASTQREADS Bam has $BAMREADS"
    exit 1
      
fi

echo "********* coverage track"
GENOME=$(echo $FASTA| sed 's/.fasta/.genome/' | sed 's/.fa/.genome/' )
java -Xmx1g -jar $IGVTOOLS count $OUT/${n/'_'$READONE.$FASTQ/.$ASD.bam} \
    $OUT/${n/'_'$READONE.$FASTQ/.$ASD.bam.cov.tdf} $GENOME


echo ">>>>> readmapping with BWA - FINISHED"
echo ">>>>> enddate "`date`

