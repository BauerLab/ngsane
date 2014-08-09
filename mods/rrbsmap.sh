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
echo ">>>>> $(basename $0) $*"


function usage {
echo -e "usage: $(basename $0) -k NGSANE -f FASTQ -r REFERENCE -o OUTDIR [OPTIONS]

Script running read mapping for single and paired DNA reads from fastq files
It expects a fastq file, paired-end, reference genome  as input and 
It runs RRBSmap, converts the output to .bam files, adds header information and
writes the coverage information for IGV.

required:
  -k | --toolkit <path>     location of the NGSANE repository 
  -f | --fastq <file>       fastq file
  -r | --reference <file>   reference genome
  -o | --outdir <path>      output dir

options:
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
ADAPTER=""
RSITE="C-CGG"         # restriction enzyme cutting site

#INPUTS
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the NGSANE repository
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
        --recover-from )        shift; NGSANE_RECOVERFROM=$1 ;; # attempt to recover from log file
        -h | --help )           usage ;;
        * )                     usage
    esac
    shift
done


#PROGRAMS
. $CONFIG
. ${NGSANE_BASE}/conf/header.sh
. $CONFIG

################################################################################
NGSANE_CHECKPOINT_INIT "programs"

# save way to load modules that itself loads other modules
hash module 2>/dev/null && for MODULE in $MODULE_RRBSMAP; do module load $MODULE; done && module list 

export PATH=$PATH_RRBS:$PATH
echo "PATH=$PATH"
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)
[ -z "$PATH_IGVTOOLS" ] && PATH_IGVTOOLS=$(dirname $(which igvtools.jar))
[ -z "$PATH_PICARD" ] && PATH_PICARD=$(dirname $(which MarkDuplicates.jar))

echo "[NOTE] set java parameters"
JAVAPARAMS="-Xmx"$(python -c "print int($MEMORY_RRBS*0.8)")"g -Djava.io.tmpdir="$TMP" -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1" 
unset _JAVA_OPTIONS
echo "JAVAPARAMS "$JAVAPARAMS

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--JAVA        --\n" $(java -Xmx200m -version 2>&1)
[ -z "$(which java)" ] && echo "[ERROR] no java detected" && exit 1
echo -e "--samtools    --\n "$(samtools 2>&1 | head -n 3 | tail -n-2)
[ -z "$(which samtools)" ] && echo "[ERROR] no samtools detected" && exit 1
echo -e "--igvtools    --\n "$(java $JAVAPARAMS -jar $PATH_IGVTOOLS/igvtools.jar version 2>&1)
[ ! -f $PATH_IGVTOOLS/igvtools.jar ] && echo "[ERROR] no igvtools detected" && exit 1
echo -e "--PICARD      --\n "$(java $JAVAPARAMS -jar $PATH_PICARD/MarkDuplicates.jar --version 2>&1)
[ ! -f $PATH_PICARD/MarkDuplicates.jar ] && echo "[ERROR] no picard detected" && exit 1
echo -e "--samstat     --\n "$(samstat -h | head -n 2 | tail -n1)
[ -z "$(which samstat)" ] && echo "[ERROR] no samstat detected" && exit 1
echo -e "--convert     --\n "$(convert -version | head -n 1)
[ -z "$(which convert)" ] && echo "[WARN] imagemagick convert not detected" 


NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "parameters"

# check library variables are set
if [[ -z "$EXPID" || -z "$LIBRARY" || -z "$PLATFORM" ]]; then
    echo "[ERROR] library info not set (EXPID, LIBRARY, and PLATFORM): free text needed"
    exit 1;
else
    echo "[NOTE] EXPID $EXPID; LIBRARY $LIBRARY; PLATFORM $PLATFORM"
fi

# get basename of f
n=${f##*/}

#is ziped ?
CAT="cat"
if [[ ${f##*.} == "gz" ]]; 
    then CAT="zcat"; 
elif [[ ${f##*.} == "bz2" ]]; 
    then CAT="bzcat"; 
fi

FASTASUFFIX=${FASTA##*.}

# delete old bam files unless attempting to recover
if [ -z "$NGSANE_RECOVERFROM" ]; then
    if [ -e $OUT/${n/%$READONE.$FASTQ/$ASD.bam} ]; then rm $OUT/${n/%$READONE.$FASTQ/$ASD.bam}; fi
    if [ -e $OUT/${n/%$READONE.$FASTQ/$ASD.bam}.stats ]; then rm $OUT/${n/%$READONE.$FASTQ/$ASD.bam}.stats; fi
fi

#is paired ?
if [ "$f" != "${f/%$READONE.$FASTQ/$READTWO.$FASTQ}" ] && [ -e ${f/%$READONE.$FASTQ/$READTWO.$FASTQ} ] && [ "$FORCESINGLE" = 0 ]; then
    PAIRED="1"
    READ1=`$CAT $f | wc -l | gawk '{print int($1/4)}' `
    READ2=`$CAT ${f/%$READONE.$FASTQ/$READTWO.$FASTQ} | wc -l | gawk '{print int($1/4)}' `
    let FASTQREADS=$READ1+$READ2

else
    PAIRED="0"
    let FASTQREADS=`$CAT $f | wc -l | gawk '{print int($1/4)}' `
fi

FULLSAMPLEID=$SAMPLEID"${n/%$READONE.$FASTQ/}"
echo "[NOTE] full sample ID "$FULLSAMPLEID

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "recall files from tape"

if [ -n "$DMGET" ]; then
	dmget -a $(dirname $FASTA)/*
	dmget -a ${f/$READONE/"*"}
	dmget -a $OUTDIR/*
fi
    
NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "mapping"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    if [ "$PAIRED" = "1" ]; then
        echo "[NOTE] PAIRED READS"
        $RRBSMAP -a $f -b ${f/%$READONE.$FASTQ/$READTWO.$FASTQ} -d $FASTA -o $OUT/${n/%$READONE.$FASTQ$ALN.bam} \
    	$ADAPTER -D $RSITE -p $CPU_RRBSMAP -2 $OUT/${n/%$READONE.$FASTQ$UNM.bam}
    else
    
        echo "[NOTE] SINGLE READS"
        $RRBSMAP -a $f -d $FASTA -o $OUT/${n/%$READONE.$FASTQ$ALN.bam} \
    	$ADAPTER -D $RSITE -p $CPU_RRBSMAP -2 $OUT/${n/%$READONE.$FASTQ$UNM.bam}
    
    fi
    
    # mark checkpoint
    [[ -s $OUT/${n/%$READONE.$FASTQ$ALN.bam} ]] && NGSANE_CHECKPOINT_CHECK 
fi 

################################################################################
NGSANE_CHECKPOINT_INIT "merge to single file"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    java $JAVAPARAMS -jar $PATH_PICARD/MergeBamAlignment.jar \
        UNMAPPED_BAM=$OUT/${n/%$READONE.$FASTQ$UNM.bam} \
        ALIGNED_BAM=$OUT/${n/%$READONE.$FASTQ$ALN.bam} \
        OUTPUT=$OUT/${n/%$READONE.$FASTQ/.ash.bam} \
        REFERENCE_SEQUENCE=$FASTA \
        PAIRED_RUN=$PAIRED \
        IS_BISULFITE_SEQUENCE=true 

    # mark checkpoint
    [[ -s $OUT/${n/%$READONE.$FASTQ/.ash.bam} ]] && NGSANE_CHECKPOINT_CHECK
fi 

################################################################################
NGSANE_CHECKPOINT_INIT "add readgroup"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then
    
    java $JAVAPARAMS -jar $PATH_PICARD/AddOrReplaceReadGroups.jar \
        INPUT=$OUT/${n/%$READONE.$FASTQ/.ash.bam} \
        OUTPUT=$OUT/${n/%$READONE.$FASTQ/.ashrg.bam} \
        RGID=$EXPID RGLB=$LIBRARY RGPL=$PLATFORM \
        RGPU=$UNIT RGSM=$FULLSAMPLEID 

    # mark checkpoint
    [[ -s $OUT/${n/%$READONE.$FASTQ/.ashrg.bam} ]] && NGSANE_CHECKPOINT_CHECK
fi 

################################################################################
NGSANE_CHECKPOINT_INIT "mark duplicates"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    if [ ! -e $OUT/metrices ]; then mkdir $OUT/metrices ; fi
    THISTMP=$TMP/$n$RANDOM #mk tmp dir because picard writes none-unique files
    mkdir $THISTMP
    java $JAVAPARAMS -jar $PATH_PICARD/MarkDuplicates.jar \
        INPUT=$OUT/${n/%$READONE.$FASTQ/.ashrg.bam} \
        OUTPUT=$OUT/${n/%$READONE.$FASTQ/$ASD.bam} \
        METRICS_FILE=$OUT/metrices/${n/%$READONE.$FASTQ/$ASD.bam}.dupl AS=true \
        VALIDATION_STRINGENCY=LENIENT \
        TMP_DIR=$THISTMP
    rm -r $THISTMP
    $SAMTOOLS index $OUT/${n/%$READONE.$FASTQ/$ASD.bam}

    # mark checkpoint
    [[ -s $OUT/${n/%$READONE.$FASTQ/$ASD.bam} ]] && NGSANE_CHECKPOINT_CHECK
fi 

################################################################################
NGSANE_CHECKPOINT_INIT "statistics"                                                                                                

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    STATSOUT=$OUT/${n/%$READONE.$FASTQ/$ASD.bam}.stats
    $SAMTOOLS flagstat $OUT/${n/%$READONE.$FASTQ/$ASD.bam} > STATSOUT
    
    if [ -n "$SEQREG" ]; then
        echo "#custom region" >> $STATSOUT
        echo $(samtools view -c -F 4 $MYOUT/${n/%$READONE.$FASTQ/.ash.bam} $SEQREG )" total reads in region " >> $STATSOUT
        echo $(samtools view -c -f 3 $MYOUT/${n/%$READONE.$FASTQ/.ash.bam} $SEQREG )" properly paired reads in region " >> $STATSOUT
    fi

    # mark checkpoint
    [[ -s $STATSOUT ]] && NGSANE_CHECKPOINT_CHECK
fi


################################################################################
NGSANE_CHECKPOINT_INIT "coverage track"    

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    java $JAVAPARAMS -jar $IGVTOOLS count $OUT/${n/%$READONE.$FASTQ/$ASD.bam} \
        $OUT/${n/%$READONE.$FASTQ/$ASD.bam.cov.tdf} ${FASTA/.$FASTASUFFIX/.genome}

    # mark checkpoint
    [[ -s $OUT/${n/%$READONE.$FASTQ/$ASD.bam.cov.tdf} ]] && NGSANE_CHECKPOINT_CHECK
fi

################################################################################
NGSANE_CHECKPOINT_INIT "verify"    
    

BAMREADS=`head -n1 $STATSOUT | cut -d " " -f 1`
if [ "$BAMREADS" = "" ]; then let BAMREADS="0"; fi			
if [ $BAMREADS -eq $FASTQREADS ]; then
    echo "[NOTE] PASS check mapping: $BAMREADS == $FASTQREADS"
    rm $OUT/${n/%$READONE$FASTQ$ALN.bam}
    rm $OUT/${n/%$READONE$FASTQ$UNM.bam}
    rm $OUT/${n/%$READONE$FASTQ/.ash.bam}
else
    echo -e "[ERROR] We are loosing reads from .fastq -> .bam in $f: \nFastq had $FASTQREADS Bam has $BAMREADS"
    exit 1
      
fi

NGSANE_CHECKPOINT_CHECK
################################################################################
echo ">>>>> readmapping with rrbsmap - FINISHED"
echo ">>>>> enddate "`date`
