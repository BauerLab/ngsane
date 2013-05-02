#!/bin/bash

echo ">>>>> HiC analysis with hiclib "
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> hiclibMapping.sh $*"

function usage {
echo -e "usage: $(basename $0) -k HISEQINF -f FASTQ -r REFERENCE -e ENZYMES -o OUTDIR [OPTIONS]

Script running hiclib pipeline tapping into bowtie2
It expects a fastq file, paired end, reference genome and digest pattern  as input.

required:
  -k | --toolkit <path>     location of the HiSeqInf repository 
  -f | --fastq <file>       fastq file
  -e | --enzymes <name>     restriction enzyme (one per library) seperated by comma
  -o | --outdir <path>      output dir

options:
  -t | --threads <nr>       number of CPUs to use (default: 1)
  -m | --memory <nr>        memory available (default: 2)

  --fastqName               name of fastq file ending (fastq.gz)
  --oldIllumina
"
exit
}

# Script to run a hic analysis based on the hiclib framework.
# It takes comma-seprated list of files containing short sequence reads in fasta or fastq format and bowtie index files as input.
# It produces output files: read alignments in .bam format and other files.
# author: Fabian Buske
# date: April 2013

# QCVARIABLES,Resource temporarily unavailable

if [ ! $# -gt 3 ]; then usage ; fi


#DEFAULTS
MYTHREADS=8
MYMEMORY=16
EXPID="exp"           # read group identifier RD ID
LIBRARY="tkcc"        # read group library RD LB
PLATFORM="illumina"   # read group platform RD PL
UNIT="flowcell"       # read group platform unit RG PU
DOBAM=1               # do the bam file
FORCESINGLE=0
NOMAPPING=0
FASTQNAME=""
ENZYME=""
QUAL="" # standard Sanger

#INPUTS                                                                                                           
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the HiSeqInf repository                       
        -t | --threads )        shift; THREADS=$1 ;; # number of CPUs to use                                      
        -m | --memory )         shift; MEMORY=$1 ;; # memory used 
        -f | --fastq )          shift; f=$1 ;; # fastq file                                                       
        -e | --enzymes )        shift; ENZYME=$1 ;; # digestion patterns
        -o | --outdir )         shift; MYOUT=$1 ;; # output dir                                                     
        --fastqName )           shift; FASTQNAME=$1 ;; #(name of fastq or fastq.gz)
        -h | --help )           usage ;;
        * )                     echo "don't understand "$1
    esac
    shift
done

if [ -z "$ENZYME" ]; then
	echo "[ERROR] restriction enzyme not specified"
	exit 1
fi

if [ -z "BOWTIE2INDEX" ]; then
	echo "[ERROR] bowtie index not specified"
	exit 1
fi

#PROGRAMS
. $CONFIG
. $HISEQINF/conf/header.sh
. $CONFIG

#JAVAPARAMS="-Xmx"$MYMEMORY"g -Djava.io.tmpdir="$TMP # -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1 -XX:MaxDirectMemorySize=4G"
#echo "JAVAPARAMS "$JAVAPARAMS

echo "********** programs"
for MODULE in $MODULE_HICLIB; do module load $MODULE; done  # save way to load modules that itself load other modules

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
	echo "[ERROR] hiclib requires paired-end fastq files"
	exit 1
fi

if [ -n "$DMGET" ]; then
	echo "********** reacall files from tape"
	dmget -a $(dirname $FASTA)/*
	dmls -l $FASTA*
	dmget -a ${f/$READONE/"*"}
	dmls -l ${f/$READONE/"*"}
fi

#is ziped ?                                                                                                       
ZCAT="zcat"
if [[ ${f##*.} != "gz" ]]; then ZCAT="cat"; fi

echo "********* reads" 
FASTQNAME=${f##*/}
READS="$f ${f/$READONE/$READTWO}"
#READ1=`$ZCAT $f | wc -l | gawk '{print int($1/4)}' `
#READ2=`$ZCAT ${f/$READONE/$READTWO} | wc -l | gawk '{print int($1/4)}' `
#let FASTQREADS=$READ1+$READ2

echo "********** hiclib call"
# run hiclib.py
PARAMS="--restrictionEnzyme $ENZYME \
   --experimentName $(echo ${ENZYME}_${FASTQNAME/$READONE.$FASTQ/} | sed 's/_*$//g') \
   --gapFile $HICLIB_GAPFILE \
   --referenceGenome $FASTA \
   --index $BOWTIE2_INDEX"

if [ "$FASTQSUFFIX" = "sra" ]; then
	PARAMS="$PARAMS --inputFormat sra --sra-reader $(which fastq-dump)"
elif [ "$FASTQSUFFIX" = "bam" ]; then
	PARAMS="$PARAMS --inputFormat bam"
else
	PARAMS="$PARAMS --inputFormat fastq"
fi

if [ -n "$HICLIB_READLENGTH" ]; then
	PARAMS="$PARAMS --readLength $HICLIB_READLENGTH"
fi

python $HISEQINF/bin/hiclibMapping.py ${PARAMS} --bowtie $(which bowtie2) --cpus $THREADS --outputDir $MYOUT --tmpDir $TMP --verbose $READS

echo "********* merge bam files"

for R in $READS; do
	READNAME=${R##*/}
	ALIGNMENTS=""
	for i in $(ls $MYOUT/${READQNAME/$FASTQ/}*.bam.*); do
		ALIGNMENTS=$ALIGNMENTS" INPUT=$i"
	done

	java $JAVAPARAMS -jar $PATH_PICARD/MergeSamFiles \
	    $ALIGNMENTS \
    	    OUTPUT=$MYOUT/${READNAME}.bam \ 
	    USE_THREADING=true
done

# copy heatmap
RUNSTATS=$OUT/runStats/hiclib
mkdir -p $RUNSTATS
mv $MYOUT/*.pdf $RUNSTATS

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
