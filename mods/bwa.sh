#!/bin/bash

# BWA calling script
# author: Denis C. Bauer
# date: Nov.2010

# messages to look out for -- relevant for the QC.sh script:
# QCVARIABLES,We are loosing reads,MAPQ should be 0 for unmapped read,no such file,file not found,bwa.sh: line,Resource temporarily unavailable

#module load bwa
#BWA=bwa

echo ">>>>> readmapping with BWA "
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> bwa.sh $*"


function usage {
echo -e "usage: $(basename $0) -k NGSANE -f FASTQ -r REFERENCE -o OUTDIR [OPTIONS]

Script running read mapping for single and paired DNA reads from fastq files
It expects a fastq file, pairdend, reference genome  as input and 
It runs BWA, converts the output to .bam files, adds header information and
writes the coverage information for IGV.

required:
  -k | --toolkit <path>     location of the NGSANE repository 
  -f | --fastq <file>       fastq file
  -r | --reference <file>   reference genome
  -o | --outdir <path>      output dir

options:
  -t | --threads <nr>       number of CPUs to use (default: 1)
  -m | --memory <nr>        memory available (default: 2)
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
MYTHREADS=1
MYMEMORY=2
#EXPID="exp"           # read group identifier RD ID
#LIBRARY="qbi"         # read group library RD LB
#PLATFORM="illumina"   # read group platform RD PL
#UNIT="flowcell"       # read group platform unit RG PU
DOBAM=1               # do the bam file
FORCESINGLE=0
NOMAPPING=0
FASTQNAME=""
QUAL="" # standard Sanger


#INPUTS
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the NGSANE repository
        -t | --threads )        shift; MYTHREADS=$1 ;; # number of CPUs to use
        -m | --memory )         shift; MYMEMORY=$1 ;; # memory used
        -f | --fastq )          shift; f=$1 ;; # fastq file
        -r | --reference )      shift; FASTA=$1 ;; # reference genome
        -o | --outdir )         shift; MYOUT=$1 ;; # output dir
        -i | --rgid )           shift; EXPID=$1 ;; # read group identifier RD ID
        -l | --rglb )           shift; LIBRARY=$1 ;; # read group library RD LB
        -p | --rgpl )           shift; PLATFORM=$1 ;; # read group platform RD PL
        -s | --rgsi )           shift; SAMPLEID=$1 ;; # read group sample RG SM (pre)
        -R | --region )         shift; SEQREG=$1 ;; # (optional) region of specific interest, e.g. targeted reseq
        -S | --sam )            DOBAM=0 ;;
        --fastqName )           shift; FASTQNAME=$1 ;; #(name of fastq or fastq.gz)
        --forceSingle )         FORCESINGLE=1;;
        --noMapping )           NOMAPPING=1;;
        --oldIllumina )         QUAL="-S";;   # old illumina encoding 1.3+
        -h | --help )           usage ;;
        * )                     echo "don't understand "$1
    esac
    shift
done

#PROGRAMS
. $CONFIG
. ${NGSANE_BASE}/conf/header.sh
. $CONFIG

JAVAPARAMS="-Xmx"$MYMEMORY"g -Djava.io.tmpdir="$TMP # -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1 -XX:MaxDirectMemorySize=4G"
echo "JAVAPARAMS "$JAVAPARAMS

echo "********** programs"
for MODULE in $MODULE_BWA; do module load $MODULE; done  # save way to load modules that itself load other modules
export PATH=$PATH_BWA:$PATH
module list
echo "PATH=$PATH"
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)
PATH_IGVTOOLS=$(dirname $(which igvtools.jar))
PATH_PICARD=$(dirname $(which MarkDuplicates.jar))
echo -e "--JAVA    --\n" $(java $JAVAPARAMS -version 2>&1)
[ -z "$(which java)" ] && echo "[ERROR] no java detected" && exit 1
echo -e "--bwa     --\n "$(bwa 2>&1 | head -n 3 | tail -n-2)
[ -z "$(which bwa)" ] && echo "[ERROR] no bwa detected" && exit 1
echo -e "--samtools--\n "$(samtools 2>&1 | head -n 3 | tail -n-2)
[ -z "$(which samtools)" ] && echo "[ERROR] no samtools detected" && exit 1
echo -e "--R       --\n "$(R --version | head -n 3)
[ -z "$(which R)" ] && echo "[ERROR] no R detected" && exit 1
echo -e "--igvtools--\n "$(java -jar $JAVAPARAMS $PATH_IGVTOOLS/igvtools.jar version 2>&1)
[ ! -f $PATH_IGVTOOLS/igvtools.jar ] && echo "[ERROR] no igvtools detected" && exit 1
echo -e "--PICARD  --\n "$(java -jar $JAVAPARAMS $PATH_PICARD/MarkDuplicates.jar --version 2>&1)
[ ! -f $PATH_PICARD/MarkDuplicates.jar ] && echo "[ERROR] no picard detected" && exit 1
echo -e "--samstat --\n "$(samstat -h | head -n 2 | tail -n 1 )
[ -z "$(which samstat)" ] && echo "[ERROR] no samstat detected" && exit 1

if [[ -z "$EXPID" || -z "$LIBRARY" || -z "$PLATFORM" ]]; then
    echo "[ERROR] library info not set (EXPID, LIBRARY, and PLATFORM): free text needed"
    exit 1;
else
    echo "[NOTE] EXPID $EXPID; LIBRARY $LIBRARY; PLATFORM $PLATFORM"
fi

# get basename of f
n=${f##*/}

# check library variables are set
if [[ -z "$EXPID" || -z "$LIBRARY" || -z "$PLATFORM" ]]; then
    echo "[ERROR] library info not set (EXPID, LIBRARY, and PLATFORM): free text needed"
    exit 1;
else
    echo "[NOTE] EXPID $EXPID; LIBRARY $LIBRARY; PLATFORM $PLATFORM"
fi


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

if [ -n "$DMGET" ]; then
	echo "********** reacall files from tape"
	dmget -a $(dirname $FASTA)/*
	dmget -a $(dirname $(which samtools))/*
	dmget -a $(dirname $(which bwa))/*
	dmget -a $PATH_PICARD/*
	dmget -a ${f/$READONE/"*"}
fi

#is ziped ?
ZCAT="zcat"
if [[ $f != *.gz ]]; then ZCAT="cat"; fi

FULLSAMPLEID=$SAMPLEID"${n/'_'$READONE.$FASTQ/}"
echo ">>>>> full sample ID "$FULLSAMPLEID
FASTASUFFIX=${FASTA##*.}

# generating the index files
if [ ! -e $FASTA.bwt ]; then echo ">>>>> make .bwt"; bwa index -a bwtsw $FASTA; fi
if [ ! -e $FASTA.fai ]; then echo ">>>>> make .fai"; samtools faidx $FASTA; fi

echo "********* mapping"
# Paired read
if [ "$PAIRED" = 1 ]
then
    if [ "$NOMAPPING" = 0 ]; then
    echo "********* PAIRED READS"
    bwa aln $QUAL -t $MYTHREADS $FASTA $f > $MYOUT/${n/$FASTQ/sai}
    bwa aln $QUAL -t $MYTHREADS $FASTA ${f/$READONE/$READTWO} > $MYOUT/${n/$READONE.$FASTQ/$READTWO.sai}
    bwa sampe $FASTA $MYOUT/${n/$FASTQ/sai} $MYOUT/${n/$READONE.$FASTQ/$READTWO.sai} \
	-r "@RG\tID:$EXPID\tSM:$FULLSAMPLEID\tPL:$PLATFORM\tLB:$LIBRARY" \
	$f ${f/$READONE/$READTWO} >$MYOUT/${n/'_'$READONE.$FASTQ/.$ALN.sam}

    rm $MYOUT/${n/$FASTQ/sai}
    rm $MYOUT/${n/$READONE.$FASTQ/$READTWO.sai}
    fi
    READ1=$($ZCAT $f | wc -l | gawk '{print int($1/4)}')
    READ2=$($ZCAT ${f/$READONE/$READTWO} | wc -l | gawk '{print int($1/4)}')
    let FASTQREADS=$READ1+$READ2
# Single read
else
    echo "********* SINGLE READS"
    bwa aln $QUAL -t $MYTHREADS $FASTA $f > $MYOUT/${n/$FASTQ/sai}

    bwa samse $FASTA $MYOUT/${n/$FASTQ/sai} \
	-r "@RG\tID:$EXPID\tSM:$FULLSAMPLEID\tPL:$PLATFORM\tLB:$LIBRARY" \
	$f >$MYOUT/${n/'_'$READONE.$FASTQ/.$ALN.sam}

#    bwa samse $FASTA $MYOUT/${n/$FASTQ/sai} \
#	-i $EXPID \
#	-m $FULLSAMPLEID \
#	-p $PLATFORM \
#	-l $LIBRARY \
#	$f >$MYOUT/${n/'_'$READONE.$FASTQ/.$ALN.sam}

    rm $MYOUT/${n/$FASTQ/sai}
    let FASTQREADS=$($ZCAT $f | wc -l | gawk '{print int($1/4)}')
fi

# exit if only the sam file is required
if [ "$DOBAM" = 0 ]; then
    SAMREADS=`grep -v @ $MYOUT/${n/'_'$READONE.$FASTQ/.$ALN.sam} | wc -l`
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
samtools view -bt $FASTA.fai $MYOUT/${n/'_'$READONE.$FASTQ/.$ALN.sam} | samtools sort - $MYOUT/${n/'_'$READONE.$FASTQ/.ash}

if [ "$PAIRED" = "1" ]; then
    # fix mates
    samtools sort -n $MYOUT/${n/'_'$READONE.$FASTQ/.ash}.bam $MYOUT/${n/'_'$READONE.$FASTQ/.ash}.tmp
    samtools fixmate $MYOUT/${n/'_'$READONE.$FASTQ/.ash}.tmp.bam $MYOUT/${n/'_'$READONE.$FASTQ/.ash}
    rm $MYOUT/${n/'_'$READONE.$FASTQ/.ash}.tmp.bam
fi

#TODO look at samtools for rmdup
#val string had to be set to LENIENT (SIlENT) to avoid crash due to a definition dis-
#agreement between bwa and picard
#http://seqanswers.com/forums/showthread.php?t=4246
echo "********* mark duplicates"
if [ ! -e $MYOUT/metrices ]; then mkdir -p $MYOUT/metrices ; fi
THISTMP=$TMP/$n$RANDOM #mk tmp dir because picard writes none-unique files
echo $THISTMP
mkdir -p $THISTMP
java $JAVAPARAMS -jar $PATH_PICARD/MarkDuplicates.jar \
    INPUT=$MYOUT/${n/'_'$READONE.$FASTQ/.ash.bam} \
    OUTPUT=$MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam} \
    METRICS_FILE=$MYOUT/metrices/${n/'_'$READONE.$FASTQ/.$ASD.bam}.dupl AS=true \
    VALIDATION_STRINGENCY=SILENT \
    TMP_DIR=$THISTMP
rm -rf $THISTMP
samtools index $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}



# statistics
echo "********* statistics"
STATSMYOUT=$MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}.stats
samtools flagstat $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam} > $STATSMYOUT
if [ -n $SEQREG ]; then
    echo "#custom region" >> $STATSMYOUT
    echo $(samtools view $MYOUT/${n/'_'$READONE.$FASTQ/.ash.bam} $SEQREG | wc -l)" total reads in region " >> $STATSMYOUT
    echo $(samtools view -f 2 $MYOUT/${n/'_'$READONE.$FASTQ/.ash.bam} $SEQREG | wc -l)" properly paired reads in region " >> $STATSMYOUT
fi


echo "********* calculate inner distance"
export PATH=$PATH:/usr/bin/
THISTMP=$TMP/$n$RANDOM #mk tmp dir because picard writes none-unique files
mkdir -p $THISTMP
java $JAVAPARAMS -jar $PATH_PICARD/CollectMultipleMetrics.jar \
    INPUT=$MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam} \
    REFERENCE_SEQUENCE=$FASTA \
    OUTPUT=$MYOUT/metrices/${n/'_'$READONE.$FASTQ/.$ASD.bam} \
    VALIDATION_STRINGENCY=SILENT \
    PROGRAM=CollectAlignmentSummaryMetrics \
    PROGRAM=CollectInsertSizeMetrics \
    PROGRAM=QualityScoreDistribution \
    TMP_DIR=$THISTMP
for im in $( ls $MYOUT/metrices/*.pdf ); do
    convert $im ${im/pdf/jpg}
done
rm -rf $THISTMP

echo "********* verify"
BAMREADS=`head -n1 $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}.stats | cut -d " " -f 1`
if [ "$BAMREADS" = "" ]; then let BAMREADS="0"; fi			
if [ $BAMREADS -eq $FASTQREADS ]; then
    echo "-----------------> PASS check mapping: $BAMREADS == $FASTQREADS"
    rm $MYOUT/${n/'_'$READONE.$FASTQ/.$ALN.sam}
    rm $MYOUT/${n/'_'$READONE.$FASTQ/.ash.bam}
else
    echo -e "***ERROR**** We are loosing reads from .fastq -> .bam in $f: \nFastq had $FASTQREADS Bam has $BAMREADS"
    exit 1 
fi

echo "********* coverage track"
GENOME=$(echo $FASTA | sed 's/.${FASTASUFFIX}/.genome/' )
java $JAVAPARAMS -jar $PATH_IGVTOOLS/igvtools.jar count $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam} \
    $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam.cov.tdf} $GENOME


echo "********* samstat"
samstat $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}

echo ">>>>> readmapping with BWA - FINISHED"
echo ">>>>> enddate "`date`

