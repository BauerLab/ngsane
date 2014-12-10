#!/bin/bash -e

# STAR calling script
# author: Boris Guennewig
# date: September 2014
# modified: 

# messages to look out for -- relevant for the QC.sh script:
# QCVARIABLES,We are loosing reads,MAPQ should be 0 for unmapped read,no such file,file not found,star.sh: line,Resource temporarily unavailable
# RESULTFILENAME <DIR>/<TASK>/<SAMPLE>$ASD.bam

echo ">>>>> readmapping with STAR "
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"


function usage {
echo -e "usage: $(basename $0) -k NGSANE -f FASTQ -r REFERENCE -o OUTDIR [OPTIONS]

Script running read mapping for single and paired DNA reads from fastq files
It expects a fastq file, pairdend, STAR/GMAP index as input and 
It runs STAR, converts the output to .bam files, adds header information.

required:
  -k | --toolkit <path>     location of the NGSANE repository 
  -f | --fastq <file>       fastq file
  -o | --outdir <path>      output dir
"
exit
}

if [ ! $# -gt 3 ]; then usage ; fi

#INPUTS
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the NGSANE repository
        -f | --fastq )          shift; f=$1 ;; # fastq file
        -o | --outdir )         shift; OUTDIR=$1 ;; # output dir
        -s | --rgsi )           shift; SAMPLEID=$1 ;; # read group prefix
        --recover-from )        shift; RECOVERFROM=$1 ;; # attempt to recover from log file
        -h | --help )           usage ;;
        * )                     echo "don't understand "$1
    esac
    shift
done

#PROGRAMS
. $CONFIG
. ${NGSANE_BASE}/conf/header.sh
. $CONFIG

################################################################################
NGSANE_CHECKPOINT_INIT "programs"

for MODULE in $MODULE_STAR; do module load $MODULE; done  # save way to load modules that itself load other modules
export PATH=$PATH_STAR:$PATH
module list
echo "PATH=$PATH"
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)
PATH_PICARD=$(dirname $(which MarkDuplicates.jar))

echo "[NOTE] set java parameters"
JAVAPARAMS="-Xmx"$(python -c "print int($MEMORY_STAR*0.8)")"g -Djava.io.tmpdir="$TMP"  -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1" 
unset _JAVA_OPTIONS
echo "JAVAPARAMS "$JAVAPARAMS

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--JAVA        --\n" $(java -Xmx200m -version 2>&1)
[ -z "$(which java)" ] && echo "[ERROR] no java detected" && exit 1
echo -e "--STAR        --\n "$(STAR --version 2>&1 | head -n 3 | tail -n-2)
[ -z "$(which STAR)" ] && echo "[ERROR] no STAR detected" && exit 1
echo -e "--samtools    --\n "$(samtools 2>&1 | head -n 3 | tail -n-2)
[ -z "$(which samtools)" ] && echo "[ERROR] no samtools detected" && exit 1
echo -e "--PICARD      --\n "$(java $JAVAPARAMS -jar $PATH_PICARD/MarkDuplicates.jar --version 2>&1)
[ ! -f $PATH_PICARD/MarkDuplicates.jar ] && echo "[ERROR] no picard detected" && exit 1
echo -e "--samstat     --\n "$(samstat -h | head -n 2 | tail -n1)
[ -z "$(which samstat)" ] && echo "[ERROR] no samstat detected" && exit 1
echo -e "--bedToBigBed --\n "$(bedToBigBed 2>&1 | tee | head -n 1 )
[ -z "$(which bedToBigBed)" ] && echo "[WARN] bedToBigBed not detected, cannot compress bedgraphs"

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "parameters"

# get basename of f
n=${f##*/}
SAMPLE=${n/%$READONE.$FASTQ/}

if [ -z "$FASTA" ]; then
    echo "[ERROR] no reference provided (FASTA)"
    exit 1
fi

# if [[ ! -f $FASTA/*.maps ]]; then
#     echo "[ERROR] GMAP/STAR index not detected. Exeute: gmap_build —D <directory> -d <genome_name> [-s none] [-k <kmer size>] <file.fasta>"
#     exit 1
# fi

#remove old files
if [ -z "$RECOVERFROM" ]; then
    [ -d $OUTDIR/$SAMPLE ] && rm -r $OUTDIR/$SAMPLE
    [ -f $OUTDIR/$SAMPLE$ASD.bam ] && rm -f $OUTDIR/$SAMPLE$ASD.bam*
fi

#is paired ?
if [ "$f" != "${f/%$READONE.$FASTQ/$READTWO.$FASTQ}" ] && [ -e ${f/%$READONE.$FASTQ/$READTWO.$FASTQ} ]; then
    echo "[NOTE] PAIRED library"
    PAIRED="1"
    f2=${f/%$READONE.$FASTQ/$READTWO.$FASTQ}
    BAM2BW_OPTION_ISPAIRED="True"
else
    echo "[NOTE] SINGLE library"
    BAM2BW_OPTION_ISPAIRED="False"
    PAIRED="0"
fi

## is ziped ?
ZCAT="cat" # always cat
if [[ ${f##*.} == "gz" ]]; then # unless its zipped
    ZCAT="zcat" 
elif [[ ${f##*.} == "bz2" ]]; then 
    ZCAT="bzcat"; 
fi

if [[ -z "$EXPID" || -z "$LIBRARY" || -z "$PLATFORM" ]]; then
    echo "[ERROR] library info not set (EXPID, LIBRARY, and PLATFORM): free text needed"
    exit 1;
else
    echo "[NOTE] EXPID $EXPID; LIBRARY $LIBRARY; PLATFORM $PLATFORM"
fi
FULLSAMPLEID=$SAMPLEID"$SAMPLE"
RG="--outSAMattrRGline \"ID:$EXPID\" \"SM:$FULLSAMPLEID\" \"LB:$LIBRARY\" \"PL:$PLATFORM\""

THISTMP=$TMP"/"$(whoami)"/"$(echo $OUTDIR/$SAMPLE | md5sum | cut -d' ' -f1)
[ -d $THISTMP ] && rm -r $THISTMP


NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "recall files from tape"

if [ -n "$DMGET" ]; then
	dmget -a $(dirname $FASTA)/*
    dmget -a $(dirname $INDEX)/*
    dmget -a $GMAP_index
    dmget -a ${f/$READONE/"*"}
	dmget -a ${OUTDIR}/$SAMPLE
fi
    
NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "STAR "

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then
    
    mkdir -p $OUTDIR/$SAMPLE
    if [ "$PAIRED" = 1 ]; then
        time STAR \
            --runMode alignReads --readFilesIn $f $f2 --readFilesCommand $ZCAT \
            --genomeDir $INDEX \
            --outFileNamePrefix $OUTDIR/$SAMPLE/ \
            --outSAMtype BAM SortedByCoordinate \
            --outTmpDir $THISTMP \
            --runThreadN $CPU_STAR $STARADDPARAM $RG
    
     else
        time STAR \
            --runMode alignReads --readFilesIn $f --readFilesCommand $ZCAT \
            --genomeDir $INDEX \
            --outFileNamePrefix $OUTDIR/$SAMPLE/ \
            --outSAMtype BAM SortedByCoordinate \
            --outTmpDir $THISTMP \
            --runThreadN $CPU_STAR $STARADDPARAM $RG
    fi
    
    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE/Aligned.sortedByCoord.out.bam
fi
################################################################################
NGSANE_CHECKPOINT_INIT "create bigwig "

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    time STAR \
        --runMode inputAlignmentsFromBAM --inputBAMfile $OUTDIR/$SAMPLE/Aligned.sortedByCoord.out.bam \
        --outWigType bedGraph --outWigNorm RPM \
        --outTmpDir $THISTMP \
        --outWigReferencesPrefix $OUTDIR/$SAMPLE/ \
        --runThreadN $CPU_STAR
        
    # extract chrom sizes from Bam
    samtools view -H $OUTDIR/$SAMPLE/Aligned.sortedByCoord.out.bam | fgrep -w '@SQ' | sed 's/:/\t/g' | awk '{OFS="\t";print $3,$5}' > $OUTDIR/$SAMPLE/Aligned.sortedByCoord.out.chromsizes

    bedGraphToBigWig $OUTDIR/$SAMPLE/Aligned.sortedByCoord.out.bg $OUTDIR/$SAMPLE/Aligned.sortedByCoord.out.chromsizes $OUTDIR/$SAMPLE.bw

    #grep -E "NH:i:1|^@" "$OUTDIR"Aligned.out.sam| samtools view -@ $CPU_STAR -bS - -o -| samtools sort -@ $CPU_STAR - "$OUTDIR"Aligned.out
    #samtools view -@ $CPU_STAR -h -S "$OUTDIR"Aligned.out.sam -b -o "$OUTDIR"Aligned.out
    #maybe mark duplicates for future
   
    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE.bw
fi
###############################################################################
NGSANE_CHECKPOINT_INIT "clean sam "

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    mkdir -p $THISTMP
   
    java $JAVAPARAMS -jar $PATH_PICARD/CleanSam.jar \
        INPUT=$OUTDIR/$SAMPLE/Aligned.sortedByCoord.out.bam \
        OUTPUT=$OUTDIR/$SAMPLE.cleaned.bam \
        VALIDATION_STRINGENCY=LENIENT \
        TMP_DIR=$THISTMP
        
    samtools sort -@ $CPU_STAR -o $OUTDIR/$SAMPLE$ASD.bam -O bam -T $SAMPLE -l 9 $THISTMP/$SAMPLE$ALN.bam

    samtools index $OUTDIR/$SAMPLE$ASD.bam
  
    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE.cleaned.bam

fi
################################################################################
NGSANE_CHECKPOINT_INIT "mark duplicates "

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then
   
    if [ ! -e $OUTDIR/metrices ]; then mkdir -p $OUTDIR/metrices ; fi

    java $JAVAPARAMS -jar $PATH_PICARD/MarkDuplicates.jar \
        INPUT=$OUTDIR/$SAMPLE.cleaned.bam \
        OUTPUT=$OUTDIR/$SAMPLE$ASD.bam \
        METRICS_FILE=$OUTDIR/metrices/$SAMPLE$ASD.bam.dupl \
        AS=true \
        CREATE_MD5_FILE=true \
        COMPRESSION_LEVEL=9 \
        VALIDATION_STRINGENCY=LENIENT \
        TMP_DIR=$THISTMP

    samtools index $OUTDIR/$SAMPLE$ASD.bam

    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE$ASD.bam 
    
    # cleanup
    [ -e $OUTDIR/$SAMPLE.cleaned.bam ] && rm $OUTDIR/$SAMPLE.cleaned.bam
fi
###############################################################################
NGSANE_CHECKPOINT_INIT "flagstat"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    echo "[NOTE] samtools flagstat"
    samtools flagstat $OUTDIR/$SAMPLE$ASD.bam > $OUTDIR/$SAMPLE$ASD.bam.stats
    READ1=$($ZCAT $f | wc -l | gawk '{print int($1/4)}' )
    FASTQREADS=$READ1
    if [ -n "$f2" ]; then 
        READ2=$($ZCAT $f2 | wc -l | gawk '{print int($1/4)}' );
        let FASTQREADS=$READ1+$READ2
    fi
    echo $FASTQREADS" fastq reads" >> $OUTDIR/$SAMPLE$ASD.bam.stats

    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE$ASD.bam.stats
fi 
################################################################################
NGSANE_CHECKPOINT_INIT "samstat"    

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    echo "[NOTE] samstat"
    samstat $OUTDIR/$SAMPLE$ASD.bam 2>&1 | tee | grep -v -P "Bad x in routine betai"
  
    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK

fi
###############################################################################
NGSANE_CHECKPOINT_INIT "cleanup"

    rm $OUTDIR/$SAMPLE/Aligned.out.sam $OUTDIR/$SAMPLE$ASD.bam 

NGSANE_CHECKPOINT_CHECK
################################################################################
[ -e $OUTDIR/$SAMPLE$ASD.bam.dummy ] && rm $OUTDIR/$SAMPLE$ASD.bam.dummy
echo ">>>>> readmapping with STAR - FINISHED"
echo ">>>>> enddate "`date`