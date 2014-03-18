#!/bin/bash -e

# Script to run TopHat --fusion search
# Fabian Buske and Hugh French
# date: March 2014


# messages to look out for -- relevant for the QC.sh script:
# QCVARIABLES,truncated file
# RESULTFILENAME <DIR>/<TASK>/<SAMPLE>.$ASD.bam

echo ">>>>> readmapping with Tophat "
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"


function usage {
echo -e "usage: $(basename $0) -k NGSANE -f FASTQ -o OUTDIR [OPTIONS]

Script running read mapping for single and paired DNA reads from fastq files
It expects a fastq file, pairdend, reference genome  as input and 
It runs tophat, converts the output to .bam files, adds header information.

required:
  -k | --toolkit <path>     location of the NGSANE repository 
  -f | --fastq <file>       fastq file
  -o | --outdir <path>      output dir

options:
  -R | --region <ps>        region of specific interest, e.g. targeted reseq
                             format chr:pos-pos
  --forceSingle             run single end eventhough second read is present
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
	-f | --fastq )          shift; f=$1 ;; # fastq file
	-o | --outdir )         shift; OUTDIR=$1 ;; # output dir
	-R | --region )         shift; SEQREG=$1 ;; # (optional) region of specific interest, e.g. targeted reseq
	--forceSingle )         FORCESINGLE=1;;
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

for MODULE in $MODULE_TOPHAT; do module load $MODULE; done  # save way to load modules that itself load other modules
export PATH=$PATH_TOPHAT:$PATH
module list
echo "PATH=$PATH"

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)

echo -e "--tophat      --\n "$(tophat --version)
[ -z "$(which tophat)" ] && echo "[ERROR] no tophat detected" && exit 1
echo -e "--bowtie1     --\n "$(bowtie1 --version)
[ -z "$(which bowtie1)" ] && echo "[ERROR] no bowtie1 detected" && exit 1
echo -e "--samtools    --\n "$(samtools 2>&1 | head -n 3 | tail -n-2)
[ -z "$(which samtools)" ] && echo "[ERROR] no samtools detected" && exit 1



echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="parameters"

[ ! -f $f ] && echo "[ERROR] input file not found: $f" 1>&2 && exit 1

# get basename of f (samplename)
n=${f##*/}

if [[ ! -e ${FASTA%.*}.1.bt2 ]]; then
    echo "[ERROR] Bowtie2 index not detected. Exeute bowtieIndex.sh first"
    exit 1
fi

# get info about input file
BAMFILE=$OUTDIR/../${n/%$READONE.$FASTQ/.$ASD.bam}

#remove old files
if [ -z "$RECOVERFROM" ]; then
    if [ -d $OUTDIR ]; then rm -r $OUTDIR; fi
fi


if [ -n "$READONE" ] && [ "$READONE" == "$READTWO" ]; then
	echo "[ERROR] read1 == read2 " 1>&2 && exit 1
elif [ "$f" != "${f/%$READONE.$FASTQ/$READTWO.$FASTQ}" ] && [ -e ${f/%$READONE.$FASTQ/$READTWO.$FASTQ} ] && [ "$FORCESINGLE" = 0 ]; then
    PAIRED="1"
    f2=${f/%$READONE.$FASTQ/$READTWO.$FASTQ}
    echo "[NOTE] Paired library detected"
else
    PAIRED="0"
    echo "[NOTE] Single-Strand (unpaired) library detected"
fi


## is ziped ?
ZCAT="cat" # always cat
if [[ $f = *.gz ]]; then # unless its zipped
    ZCAT="zcat";
fi

# get encoding
if [ -z "$FASTQ_PHRED" ]; then 
    FASTQ_ENCODING=$($ZCAT $f |  awk 'NR % 4 ==0' | python $NGSANE_BASE/tools/GuessFastqEncoding.py |  tail -n 1)
    if [[ "$FASTQ_ENCODING" == *Phred33* ]]; then
        FASTQ_PHRED="" # use default
    elif [[ "$FASTQ_ENCODING" == *Illumina* ]]; then
        FASTQ_PHRED="--phred64-quals"
    elif [[ "$FASTQ_ENCODING" == *Solexa* ]]; then
        FASTQ_PHRED="--solexa1.3-quals"
    else
        echo "[NOTE] cannot detect/don't understand fastq format: $FASTQ_ENCODING - using default"
    fi
    echo "[NOTE] $FASTQ_ENCODING fastq format detected"
fi


#ln -s /share/ClusterShare/biodata/contrib/fusion/refGene.txt $PWD/refGene.txt
#ln -s /share/ClusterShare/biodata/contrib/fusion/ensGene.txt $PWD/ensGene.txt
#ln -s /share/ClusterShare/biodata/contrib/fusion/mcl $PWD/mcl

#cp -rs /share/ClusterShare/biodata/contrib/fusion/blast/ $PWD/blast/


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

# check library info is set
if [ -z "$RNA_SEQ_LIBRARY_TYPE" ]; then
    echo "[ERROR] RNAseq library type not set (RNA_SEQ_LIBRARY_TYPE): either fr-unstranded or fr-firststrand"
    exit 1;
else
    echo "[NOTE] RNAseq library type: $RNA_SEQ_LIBRARY_TYPE"
fi
if [[ -z "$EXPID" || -z "$LIBRARY" || -z "$PLATFORM" ]]; then
    echo "[ERROR] library info not set (EXPID, LIBRARY, and PLATFORM): free text needed"
    exit 1;
else
    echo "[NOTE] EXPID $EXPID; LIBRARY $LIBRARY; PLATFORM $PLATFORM"
fi

if [ -n "$TOPHATTRANSCRIPTOMEINDEX" ]; then
    echo "[NOTE] RNAseq --transcriptome-index specified: ${TOPHATTRANSCRIPTOMEINDEX}"
    TOPHAT_TRANSCRIPTOME_PARAM="--transcriptome-index=${TOPHATTRANSCRIPTOMEINDEX}"

else
    echo "[NOTE] no --transcriptome-index specified."
    TOPHAT_TRANSCRIPTOME_PARAM=
fi

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="recall files from tape"

if [ -n "$DMGET" ]; then
    dmget -a $(dirname $FASTA)/*
    dmget -a ${f/$READONE/"*"}
    dmget -a $OUTDIR/*
fi

echo -e "\n********* $CHECKPOINT\n"


################################################################################
CHECKPOINT="run tophat"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

    echo "[NOTE] tophat $(date)"


   
  RUN_COMMAND="tophat $FASTQ_PHRED --fusion-search --keep-fasta-order --bowtie1 --no-coverage-search --max-intron-length 100000 --fusion-min-dist 100000 --fusion-anchor-length 13 --fusion-ignore-chromosomes chrM --num-threads $CPU_TOPHAT --library-type $RNA_SEQ_LIBRARY_TYPE --rg-id $EXPID --rg-sample $PLATFORM --rg-library $LIBRARY --output-dir $OUTDIR ${FASTA%.*} $f $f2"
    echo $RUN_COMMAND && eval $RUN_COMMAND
    echo "[NOTE] tophat end $(date)"

  RUN_COMMAND="tophat-fusion-post --num-fusion-reads 1 --num-fusion-pairs 2 --num-fusion-both 5 --num-threads $CPU_TOPHAT --output-dir $OUTDIR ${FASTA%.*} $f $f2"
    echo $RUN_COMMAND && eval $RUN_COMMAND
    echo "[NOTE] tophat-fusion-post end $(date)"


    # mark checkpoint
    if [ -e $OUTDIR/accepted_hits.bam ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi 



###############################################################################
CHECKPOINT="cleanup"


echo -e "\n********* $CHECKPOINT\n"
################################################################################
[ -e ${BAMFILE}.dummy ] && rm ${BAMFILE}.dummy
echo ">>>>> alignment with TopHat - FINISHED"
echo ">>>>> enddate "`date`
