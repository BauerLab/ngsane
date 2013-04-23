#!/bin/bash

echo ">>>>> HiC analysis with hiclib "
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> hiclib.sh $*"

function usage {
echo -e "usage: $(basename $0) -k HISEQINF -f FASTQ -r REFERENCE -o OUTDIR [OPTIONS]"
exit
}

# Script to run a hic analysis based on the hiclib framework.
# It takes comma-seprated list of files containing short sequence reads in fasta or fastq format and bowtie index files as input.
# It produces output files: read alignments in .bam format and other files.
# author: Fabian Buske
# date: April 2013

# QCVARIABLES,Resource temporarily unavailable

if [ ! $# -gt 3 ]; then usage ; fi

THREADS=1
MEMORY=2
EXPID="exp"           # read group identifier RD ID                                                               
LIBRARY="tkcc"        # read group library RD LB                                                                  
PLATFORM="illumina"   # read group platform RD PL                                                                 
UNIT="flowcell"       # read group platform unit RG PU                                                            
FASTQNAME=""
FORCESINGLE=0

#INPUTS                                                                                                           
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the HiSeqInf repository                       
        -t | --threads )        shift; THREADS=$1 ;; # number of CPUs to use                                      
        -m | --memory )         shift; MEMORY=$1 ;; # memory used 
        -f | --fastq )          shift; f=$1 ;; # fastq file                                                       
        -r | --reference )      shift; FASTA=$1 ;; # reference genome                                             
        -o | --outdir )         shift; MYOUT=$1 ;; # output dir                                                     
        -i | --rgid )           shift; EXPID=$1 ;; # read group identifier RD ID                                  
        -l | --rglb )           shift; LIBRARY=$1 ;; # read group library RD LB                                   
        -p | --rgpl )           shift; PLATFORM=$1 ;; # read group platform RD PL                                 
        -s | --rgsi )           shift; SAMPLEID=$1 ;; # read group sample RG SM (pre)                             
        -u | --rgpu )           shift; UNIT=$1 ;; # read group platform unit RG PU
	--fastqName )           shift; FASTQNAME=$1 ;; #(name of fastq or fastq.gz)
        -h | --help )           usage ;;
        * )                     echo "don't understand "$1
    esac
    shift
done

#PROGRAMS
. $CONFIG
. $HISEQINF/pbsTemp/header.sh
. $CONFIG

echo "********** programs"
module load $MODULE_HICLIB; 
export PATH=$PATH_HICLIB:$PATH
module list
echo $PATH
echo -e "--Python      --\n" $(python --version)
echo -e "--Python libs --\n "$(yolk -l)

n=`basename $f`

# delete old bam file                                                                                             
#if [ -e $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam} ]; then rm $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}; fi

#is paired ?                                                                                                      
if [ -e ${f/$READONE/$READTWO} ] && [ "$FORCESINGLE" = 0 ]; then
    PAIRED="1"
else
	echo "[ERROR] library needs to be paired"
	exit 1
fi

#is ziped ?                                                                                                       
ZCAT="zcat"
if [[ ${f##*.} != "gz" ]]; then ZCAT="cat"; fi

echo "********* generating the index files"
FASTASUFFIX=${FASTA##*.}
if [ ! -e ${FASTA/.${FASTASUFFIX}/}.1.bt2 ]; then echo ">>>>> make .bt2"; bowtie2-build $FASTA ${FASTA/.${FASTASUFFIX}/}; fi
if [ ! -e $FASTA.fai ]; then echo ">>>>> make .fai"; samtools faidx $FASTA; fi

if [ -n "$DMGET" ]; then
	echo "********** reacall files from tape"
	dmget -a $(dirname $FASTA)/*
	dmls -l $FASTA*
	dmget -a ${f/$READONE/"*"}
	dmls -l ${f/$READONE/"*"}
fi

#run bowtie command -v $MISMATCH -m 1
echo "********* bowtie" 
READS="$f ${f/$READONE/$READTWO}"
READ1=`$ZCAT $f | wc -l | gawk '{print int($1/4)}' `
READ2=`$ZCAT ${f/$READONE/$READTWO} | wc -l | gawk '{print int($1/4)}' `
let FASTQREADS=$READ1+$READ2

# run hiclib.py
PARAMS="--restrictionEnzyme $HICLIB_ENZYME \
   --experimentName $HICLIB_EXPERIMENT \
   --gapFile $HICLIB_GAPFILE \
   --referenceGenome $FASTA \
   --index ${FASTA/.${FASTASUFFIX}/} "
   
if [ $FASTQSUFFIX = "sra" ]; then
	PARAMS="$PARAMS --inputFormat sra --sra-reader $(which fastq-dump)"
elif [ $FASTQSUFFIX = "bam" ]; then
	PARAMS="$PARAMS --inputFormat bam"
else
	PARAMS="$PARAMS --inputFormat fastq"
fi

if [ -n "$HICLIB_READLENGTH" ]; then
	PARAMS="$PARAMS --readLength $HICLIB_READLENGTH"
fi

python run_hiclib.py \
   ${PARAMS} \
   --bowtie $(which bowtie2)
   --cpus $THREADS \
   --outputDir $MYOUT \
   --tmpDir $TMP \
   --verbose \
   $READS

#
## statistics                                                                                                      
#echo "********* statistics"
#STATSOUT=$MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}.stats
#samtools flagstat $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam} > $STATSOUT
#if [ -n $SEQREG ]; then
#    echo "#custom region" >> $STATSOUT
#    echo `samtools view $MYOUT/${n/'_'$READONE.$FASTQ/.ash.bam} $SEQREG | wc -l`" total reads in region " >> $STAT\
#SOUT
#    echo `samtools view -f 2 $MYOUT/${n/'_'$READONE.$FASTQ/.ash.bam} $SEQREG | wc -l`" properly paired reads in re\
#gion " >> $STATSOUT
#fi
#
#echo "********* calculate inner distance"
#export PATH=$PATH:/usr/bin/
#THISTMP=$TMP/$n$RANDOM #mk tmp dir because picard writes none-unique files
#mkdir $THISTMP
#java $JAVAPARAMS -jar $PATH_PICARD/CollectMultipleMetrics.jar \
#    INPUT=$MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam} \
#    REFERENCE_SEQUENCE=$FASTA \
#    OUTPUT=$MYOUT/metrices/${n/'_'$READONE.$FASTQ/.$ASD.bam} \
#    VALIDATION_STRINGENCY=LENIENT \
#    PROGRAM=CollectAlignmentSummaryMetrics \
#    PROGRAM=CollectInsertSizeMetrics \
#    PROGRAM=QualityScoreDistribution \
#    TMP_DIR=$THISTMP
#for im in $( ls $MYOUT/metrices/*.pdf ); do
#    convert $im ${im/pdf/jpg}
#done
#rm -r $THISTMP
#
#
#
#echo "********* verify"
#BAMREADS=`head -n1 $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}.stats | cut -d " " -f 1`
#if [ "$BAMREADS" = "" ]; then let BAMREADS="0"; fi
#if [ $BAMREADS -eq $FASTQREADS ]; then
#    echo "-----------------> PASS check mapping: $BAMREADS == $FASTQREADS"
#    rm $MYOUT/${n/'_'$READONE.$FASTQ/.$ALN.sam}
#    rm $MYOUT/${n/'_'$READONE.$FASTQ/.$ALN.un.sam}
#    rm $MYOUT/${n/'_'$READONE.$FASTQ/.ash.bam}
#    rm $MYOUT/${n/'_'$READONE.$FASTQ/.unm}.bam
#    rm $MYOUT/${n/'_'$READONE.$FASTQ/.map}.bam
#else
#    echo -e "***ERROR**** We are loosing reads from .fastq -> .bam in $f: \nFastq had $FASTQREADS Bam has $BAMREA\
#DS"
#    exit 1
#fi
#
##coverage for IGV
#echo "********* coverage track"
#java $JAVAPARAMS -jar $PATH_IGVTOOLS/igvtools.jar count $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam} \
#$MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam.cov.tdf} ${FASTA/$FASTASUFFIX/}.genome
#
#echo "********* samstat"
#samstat $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}
#
#echo ">>>>> readmapping with BWA - FINISHED"
#echo ">>>>> enddate "`date`
#
