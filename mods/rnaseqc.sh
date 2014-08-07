#!/bin/bash -e

# Script to run TopHat program
# It takes comma-seprated list of files containing short sequence reads in fasta or fastq format and bowtie index files as input.
# It produces output files: read alignments in .bam format and other files.
# author: Chikako Ragan, Denis Bauer
# date: Jan. 2011
# modified by Fabian Buske and Hugh French
# date: 2013-


# messages to look out for -- relevant for the QC.sh script:
# QCVARIABLES,truncated file

echo ">>>>> RNA QC with RNAseq-QC"
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"


function usage {
echo -e "usage: $(basename $0) -k NGSANE -f FASTQ -o OUTDIR [OPTIONS]

Script running RNAseq-QC on a bam file

required:
  -k | --toolkit <path>     location of the NGSANE repository 
  -f | --fastq <file>       fastq file
  -o | --outdir <path>      output dir
"
exit
}

if [ ! $# -gt 3 ]; then usage ; fi

#DEFAULTS
FORCESINGLE=0

#INPUTS
while [ "$1" != "" ]; do
	case $1 in
	-k | toolkit )          shift; CONFIG=$1 ;; # ENSURE NO VARIABLE NAMES FROM CONFIG
	-f | --file )           shift; f=$1 ;; # input file
	-o | --outdir )         shift; OUTDIR=$1 ;; # output dir
    --recover-from )        shift; RECOVERFROM=$1 ;; # attempt to recover from log file
	-h | --help )           usage ;;
	* )                     echo "dont understand $1"
	esac
	shift
done

. $CONFIG
. ${NGSANE_BASE}/conf/header.sh
. $CONFIG

################################################################################
CHECKPOINT="programs"

# save way to load modules that itself loads other modules
hash module 2>/dev/null && for MODULE in $MODULE_RNASEQC; do module load $MODULE; done && module list 

export PATH=$PATH_RNASEQC:$PATH
echo "PATH=$PATH"
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)
[ -z "$PATH_PICARD" ] && PATH_PICARD=$(dirname $(which MarkDuplicates.jar))
[ -z "$PATH_RNASEQC" ] && PATH_RNASEQC=$(dirname $(which RNA-SeQC.jar))

echo "[NOTE] set java parameters"
JAVAPARAMS="-Xmx"$(python -c "print int($MEMORY_RNASEQC*0.8)")"g -Djava.io.tmpdir="$TMP"  -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1" 
unset _JAVA_OPTIONS
echo "JAVAPARAMS "$JAVAPARAMS

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--JAVA        --\n" $(java -Xmx200m -version 2>&1)
[ -z "$(which java)" ] && echo "[ERROR] no java detected" && exit 1
echo -e "--samtools    --\n "$(samtools 2>&1 | head -n 3 | tail -n-2)
[ -z "$(which samtools)" ] && echo "[ERROR] no samtools detected" && exit 1
echo -e "--picard      --\n "$(java $JAVAPARAMS -jar $PATH_PICARD/MarkDuplicates.jar --version 2>&1)
[ ! -f $PATH_PICARD/MarkDuplicates.jar ] && echo "[ERROR] no picard detected" && exit 1
echo -e "--RNA-SeQC    --\n "$(java $JAVAPARAMS -jar ${PATH_RNASEQC}/RNA-SeQC.jar --version  2>&1 | head -n 1 )
[ -z "$(which RNA-SeQC.jar)" ] && echo "[ERROR] no RNA_SeQC.jar detected" && exit 1

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="parameters"

[ ! -f $f ] && echo "[ERROR] input file not found: $f" 1>&2 && exit 1

# get basename of f (samplename)
n=${f##*/}
SAMPLE=${n/$ASD.bam/}

if [ -z "$FASTA" ]; then
    echo "[ERROR] no reference provided (FASTA)"
    exit 1
fi

#remove old files
if [ -z "$RECOVERFROM" ]; then
    if [ -d $OUTDIR ]; then rm -r $OUTDIR; fi
fi

## GTF provided?
if [ -n "$GTF" ]; then
    echo "[NOTE] GTF: $GTF"
    if [ ! -f $GTF ]; then
        echo "[ERROR] GTF specified but not found!"
        exit 1
    fi 
    if [ ! -z "$DOCTOREDGTFSUFFIX" ]; then
        if [ ! -f ${GTF/%.gtf/$DOCTOREDGTFSUFFIX} ] ; then
            echo "[ERROR] Doctored GTF suffix specified but gtf not found: ${GTF/%.gtf/$DOCTOREDGTFSUFFIX}"
            exit 1
        else 
            echo "[NOTE] Doctored GTF: ${GTF/%.gtf/$DOCTOREDGTFSUFFIX}"
        fi
    fi
else
    echo "[NOTE] no GTF specified!"
fi

THISTMP=$TMP"/"$(whoami)"/"$(echo $OUTDIR/$SAMPLE | md5sum | cut -d' ' -f1)
[ -d $THISTMP ] && rm -r $THISTMP
mkdir -p $THISTMP

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="recall files from tape"

if [ -n "$DMGET" ]; then
    dmget -a $(dirname $FASTA)/*
    dmget -a ${f}
    dmget -a $OUTDIR/*
fi

echo -e "\n********* $CHECKPOINT\n"

################################################################################
CHECKPOINT="RNA-SeQC"    

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
	## ensure bam is properly ordered for GATK
	#reheader bam
	RUN_COMMAND="java -jar $JAVAPARAMS $PATH_PICARD/ReorderSam.jar INPUT=$f OUTPUT=$THISTMP/${SAMPLE}_unsorted.bam REFERENCE=$FASTA VALIDATION_STRINGENCY=SILENT"
    echo $RUN_COMMAND && eval $RUN_COMMAND	

	#sort
	samtools sort $THISTMP/${SAMPLE}_unsorted.bam $THISTMP/${SAMPLE}_sorted
	
	#index
	samtools index $THISTMP/${SAMPLE}_sorted.bam
	
	rm $THISTMP/${SAMPLE}_unsorted.bam
    
    # take doctored GTF if available
    if [ -n "$DOCTOREDGTFSUFFIX" ]; then 
        RNASEQC_GTF=${GTF/%.gtf/$DOCTOREDGTFSUFFIX}; 
    else
        RNASEQC_GTF=$GTF
    fi
    # run GC stratification if available
    if [ -f ${RNASEQC_GTF}.gc ]; then RNASEQC_CG="-strat gc -gc ${RNASEQC_GTF}.gc"; fi
    
    
    RUN_COMMAND="java $JAVAPARAMS -jar ${PATH_RNASEQC}/RNA-SeQC.jar $RNASEQCADDPARAM -n 1000 -s '${SAMPLE}|$THISTMP/${SAMPLE}_sorted.bam|${SAMPLE}' -t ${RNASEQC_GTF}  -r ${FASTA} -o $OUTDIR $RNASEQC_CG  2>&1 | tee | grep -v 'Ignoring SAM validation error: ERROR:'"
    # TODO: find a way to fix the bamfile reg. the SAM error
    # http://seqanswers.com/forums/showthread.php?t=28155
    echo $RUN_COMMAND && eval $RUN_COMMAND

	rm $THISTMP/${SAMPLE}_sorted.bam
	rm $THISTMP/${SAMPLE}_sorted.bam.bai

    # mark checkpoint
    if [ -d $OUTDIR/ ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi
################################################################################
echo ">>>>> RNA QC with RNAseq-QC - FINISHED"
echo ">>>>> enddate "`date`
