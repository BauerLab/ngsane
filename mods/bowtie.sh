#!/bin/bash

echo ">>>>> alignment with bowtie v1"
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> bowtie.sh $*"

function usage {
echo -e "usage: $(basename $0) -k NGSANE -f FASTQ -r REFERENCE -o OUTDIR [OPTIONS]"
exit
}

# Script to run bowtie v1 program.
# It takes comma-seprated list of files containing short sequence reads in fasta or fastq format and bowtie index files as input.
# It produces output files: read alignments in .bam format and other files.
# author: Fabian Buske
# date: June 2013

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
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the NGSANE repository                       
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
. ${NGSANE_BASE}/conf/header.sh
. $CONFIG

JAVAPARAMS="-Xmx"$MEMORY"G -Djava.io.tmpdir="$TMP #-XX:ConcGCThreads=1 -XX:ParallelGCThreads=1 -XX:MaxDirectMemorySize=10G"
echo "JAVAPARAMS "$JAVAPARAMS

echo "********** programs"
for MODULE in $MODULE_BOWTIETWO; do module load $MODULE; done  # save way to load modules that itself load other modules
export PATH=$PATH_BOWTIETWO:$PATH
module list
echo "PATH=$PATH"
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)
PATH_IGVTOOLS=$(dirname $(which igvtools.jar))
PATH_PICARD=$(dirname $(which MarkDuplicates.jar))
echo -e "--JAVA    --\n" $(java $JAVAPARAMS -version 2>&1)
[ -z "$(which java)" ] && echo "[ERROR] no java detected" && exit 1
echo -e "--bowtie  --\n "$(bowtie --version)
[ -z "$(which bowtie)" ] && echo "[ERROR] no bowtie detected" && exit 1
echo -e "--samtools--\n "$(samtools 2>&1 | head -n 3 | tail -n-2)
[ -z "$(which samtools)" ] && echo "[ERROR] no samtools detected" && exit 1
echo -e "--R       --\n "$(R --version | head -n 3)
[ -z "$(which R)" ] && echo "[ERROR] no R detected" && exit 1
echo -e "--igvtools--\n "$(java -jar $JAVAPARAMS $PATH_IGVTOOLS/igvtools.jar version 2>&1)
[ ! -f $PATH_IGVTOOLS/igvtools.jar ] && echo "[ERROR] no igvtools detected" && exit 1
echo -e "--PICARD  --\n "$(java -jar $JAVAPARAMS $PATH_PICARD/MarkDuplicates.jar --version 2>&1)
[ ! -f $PATH_PICARD/MarkDuplicates.jar ] && echo "[ERROR] no picard detected" && exit 1
echo -e "--samstat --\n "$(samstat -h | head -n 2 | tail -n1)
[ -z "$(which samstat)" ] && echo "[ERROR] no samstat detected" && exit 1

# get basename of f
n=${f##*/}


# delete old bam file                                                                                             
if [ -e $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam} ]; then rm $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}; fi
if [ -e $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}.stats ]; then rm $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}.stats; fi
if [ -e $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}.dupl ]; then rm $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}.dupl; fi

#is paired ?                                                                                                      
if [ -e ${f/$READONE/$READTWO} ] && [ "$FORCESINGLE" = 0 ]; then
    PAIRED="1"
else
    PAIRED="0"
fi

#is ziped ?                                                                                                       
ZCAT="zcat"
if [[ ${f##*.} != "gz" ]]; then ZCAT="cat"; fi

echo "********* generating the index files"
FASTASUFFIX=${FASTA##*.}
if [ ! -e ${FASTA/.${FASTASUFFIX}/}.1.ebwt ]; then echo ">>>>> make .ebwt"; bowtie-build $FASTA ${FASTA/.${FASTASUFFIX}/}; fi
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
if [ $PAIRED == "0" ]; then 
    READS="-U $f"
    let FASTQREADS=`$ZCAT $f | wc -l | gawk '{print int($1/4)}' `
else 
    READS="-1 $f -2 ${f/$READONE/$READTWO}"
    READ1=`$ZCAT $f | wc -l | gawk '{print int($1/4)}' `
    READ2=`$ZCAT ${f/$READONE/$READTWO} | wc -l | gawk '{print int($1/4)}' `
    let FASTQREADS=$READ1+$READ2
fi

#readgroup
FULLSAMPLEID=$SAMPLEID"${n/'_'$READONE.$FASTQ/}"
RG="--sam-RG \"ID:$EXPID\" --sam-RG \"SM:$FULLSAMPLEID\" --sam-RG \"LB:$LIBRARY\" --sam-RG \"PL:$PLATFORM\""

RUN_COMMAND="bowtie $RG --tryhard --best --strata --time --threads $THREADS -v 3 -m 1 --un $MYOUT/${n/'_'$READONE.$FASTQ/.$ALN.un.fq} --max $MYOUT/${n/'_'$READONE.$FASTQ/.$ALN.mul.fq} --sam $BOWTIE_OPTIONS $BOWTIEINDEX ${FASTA/.${FASTASUFFIX}/}$READS $MYOUT/${n/'_'$READONE.$FASTQ/.$ALN.sam}"
echo $RUN_COMMAND
eval $RUN_COMMAND

# create bam files for discarded reads and remove fastq files
java $JAVAPARAMS -jar $PATH_PICARD/FastqToSam.jar \
    FASTQ=$MYOUT/${n/'_'$READONE.$FASTQ/.$ALN.un.fq} \
    OUTPUT=$MYOUT/${n/'_'$READONE.$FASTQ/.$ALN.un.bam} \
    QUALITY_FORMAT=Standard \
    SAMPLE_NAME=${n/'_'$READONE.$FASTQ/} \
    READ_GROUP_NAME=null \
    QUIET=TRUE \
    VERBOSITY=ERROR

samtools sort $MYOUT/${n/'_'$READONE.$FASTQ/.$ALN.un.bam} $MYOUT/${n/'_'$READONE.$FASTQ/.$ALN.un.bam.tmp} && mv $MYOUT/${n/'_'$READONE.$FASTQ/.$ALN.un.bam.tmp} $MYOUT/${n/'_'$READONE.$FASTQ/.$ALN.un.bam}

java $JAVAPARAMS -jar $PATH_PICARD/FastqToSam.jar \
    FASTQ=$MYOUT/${n/'_'$READONE.$FASTQ/.$ALN.mul.fq} \
    OUTPUT=$MYOUT/${n/'_'$READONE.$FASTQ/.$ALN.mul.bam} \
    QUALITY_FORMAT=Standard \
    SAMPLE_NAME=${n/'_'$READONE.$FASTQ/} \
    READ_GROUP_NAME=null \
    QUIET=TRUE \
    VERBOSITY=ERROR

samtools sort $MYOUT/${n/'_'$READONE.$FASTQ/.$ALN.mul.bam} $MYOUT/${n/'_'$READONE.$FASTQ/.$ALN.mul.bam.tmp} && mv $MYOUT/${n/'_'$READONE.$FASTQ/.$ALN.mul.bam.tmp} $MYOUT/${n/'_'$READONE.$FASTQ/.$ALN.mul.bam}

rm $MYOUT/${n/'_'$READONE.$FASTQ/.$ALN.un.fq}
rm $MYOUT/${n/'_'$READONE.$FASTQ/.$ALN.mul.fq}

# continue for normal bam file conversion                                                                         
echo "********* sorting and bam-conversion"
samtools view -F4 -S -bt $FASTA.fai $MYOUT/${n/'_'$READONE.$FASTQ/.$ALN.sam} $MYOUT/${n/'_'$READONE.$FASTQ/.$ALN.bam}
samtools sort $MYOUT/${n/'_'$READONE.$FASTQ/.$ALN.bam} $MYOUT/${n/'_'$READONE.$FASTQ/.map}

# merge mappend and unmapped
samtools merge -f $MYOUT/${n/'_'$READONE.$FASTQ/.ash}.bam $MYOUT/${n/'_'$READONE.$FASTQ/.map}.bam $MYOUT/${n/'_'$READONE.$FASTQ/.unm}.bam  $MYOUT/${n/'_'$READONE.$FASTQ/.$ALN.mul.bam}

if [ "$PAIRED" = "1" ]; then
    # fix mates
    samtools sort -n $MYOUT/${n/'_'$READONE.$FASTQ/.ash}.bam $MYOUT/${n/'_'$READONE.$FASTQ/.ash}.bam.tmp
    samtools fixmate $MYOUT/${n/'_'$READONE.$FASTQ/.ash}.bam.tmp.bam - | samtools sort - $MYOUT/${n/'_'$READONE.$FASTQ/.ash}
    rm $MYOUT/${n/'_'$READONE.$FASTQ/.ash}.bam.tmp.bam
fi

echo "********* mark duplicates"
if [ ! -e $MYOUT/metrices ]; then mkdir -p $MYOUT/metrices ; fi
THISTMP=$TMP/$n$RANDOM #mk tmp dir because picard writes none-unique files                                        
mkdir -p $THISTMP
java $JAVAPARAMS -jar $PATH_PICARD/MarkDuplicates.jar \
    INPUT=$MYOUT/${n/'_'$READONE.$FASTQ/.ash.bam} \
    OUTPUT=$MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam} \
    METRICS_FILE=$MYOUT/metrices/${n/'_'$READONE.$FASTQ/.$ASD.bam}.dupl \
    AS=true \
    VALIDATION_STRINGENCY=LENIENT \
    TMP_DIR=$THISTMP
rm -rf $THISTMP
samtools index $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}


# statistics                                                                                                      
echo "********* statistics"
STATSOUT=$MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}.stats
samtools flagstat $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam} > $STATSOUT
if [ -n $SEQREG ]; then
    echo "#custom region" >> $STATSOUT
    echo $(samtools view $MYOUT/${n/'_'$READONE.$FASTQ/.ash.bam} $SEQREG | wc -l)" total reads in region " >> $STATSOUT
    echo $(samtools view -f 2 $MYOUT/${n/'_'$READONE.$FASTQ/.ash.bam} $SEQREG | wc -l)" properly paired reads in region " >> $STATSOUT
else
    echo $(samtools view -c $MYOUT/${n/'_'$READONE.$FASTQ/.$ALN.una.bam})" unaligned_reads " >> $STATSOUT
    echo $(samtools view -c $MYOUT/${n/'_'$READONE.$FASTQ/.$ALN.mul.bam})" multiple_reads " >> $STATSOUT
    echo $(samtools view -c $MYOUT/${n/'_'$READONE.$FASTQ/.$ALN.bam})" aligned_reads" >> $STATSOUT
    echo $(samtools view -c $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam})" unique_reads " >> $STATSOUT
fi


echo "********* calculate inner distance"
export PATH=$PATH:/usr/bin/
THISTMP=$TMP/$n$RANDOM #mk tmp dir because picard writes none-unique files
mkdir $THISTMP
java $JAVAPARAMS -jar $PATH_PICARD/CollectMultipleMetrics.jar \
    INPUT=$MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam} \
    REFERENCE_SEQUENCE=$FASTA \
    OUTPUT=$MYOUT/metrices/${n/'_'$READONE.$FASTQ/.$ASD.bam} \
    VALIDATION_STRINGENCY=LENIENT \
    PROGRAM=CollectAlignmentSummaryMetrics \
    PROGRAM=CollectInsertSizeMetrics \
    PROGRAM=QualityScoreDistribution \
    TMP_DIR=$THISTMP
for im in $( ls $MYOUT/metrices/*.pdf ); do
    convert $im ${im/pdf/jpg}
done
rm -r $THISTMP

echo "********* create Homer Tag directory"
TAGDIRECTORY=$MYOUT/${n/'_'$READONE.$FASTQ/_homer}
mkdir -p $TAGDIRECTORY
makeTagDirectory $TAGDIRECTORY $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}

echo "********* create Bigwig"
fetchChromSizes $GENOME > $MYOUT/$GENOME.chrom.sizes
makeUCSCfile $TAGDIRECTORY -o $MYOUT/${n/'_'$READONE.$FASTQ/.bw} -bigWig $GENOME.chrom.sizes -norm 1e6 -fsize 1e20 -fragLength 300

echo "********* verify"
BAMREADS=`head -n1 $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}.stats | cut -d " " -f 1`
if [ "$BAMREADS" = "" ]; then let BAMREADS="0"; fi
if [ $BAMREADS -eq $FASTQREADS ]; then
    echo "-----------------> PASS check mapping: $BAMREADS == $FASTQREADS"
    rm $MYOUT/${n/'_'$READONE.$FASTQ/.$ALN.sam}
    rm $MYOUT/${n/'_'$READONE.$FASTQ/.$ALN.un.sam}
    rm $MYOUT/${n/'_'$READONE.$FASTQ/.ash.bam}
    rm $MYOUT/${n/'_'$READONE.$FASTQ/.unm}.bam
    rm $MYOUT/${n/'_'$READONE.$FASTQ/.map}.bam
else
    echo -e "[ERROR] We are loosing reads from .fastq -> .bam in $f: \nFastq had $FASTQREADS Bam has $BAMREADS"
    exit 1
fi

#coverage for IGV
echo "********* coverage track"
java $JAVAPARAMS -jar $PATH_IGVTOOLS/igvtools.jar count $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam} \
$MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam.cov.tdf} ${FASTA/$FASTASUFFIX/}.genome

echo "********* samstat"
samstat $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}

echo ">>>>> readmapping with bowtie2 - FINISHED"
echo ">>>>> enddate "`date`

