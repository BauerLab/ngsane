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
# RESULTFILENAME <DIR>/<TASK>/<SAMPLE>$ASD.bam

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

# save way to load modules that itself loads other modules
hash module 2>/dev/null && for MODULE in $MODULE_TOPHAT; do module load $MODULE; done && module list 

export PATH=$PATH_TOPHAT:$PATH
echo "PATH=$PATH"
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)
[ -z "$PATH_PICARD" ] && PATH_PICARD=$(dirname $(which CleanSam.jar))

echo "[NOTE] set java parameters"
JAVAPARAMS="-Xmx"$(python -c "print int($MEMORY_TOPHAT*0.8)")"g -Djava.io.tmpdir="$TMP"  -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1" 
unset _JAVA_OPTIONS
echo "JAVAPARAMS "$JAVAPARAMS

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--JAVA        --\n" $(java -Xmx200m -version 2>&1)
[ -z "$(which java)" ] && echo "[ERROR] no java detected" && exit 1
echo -e "--tophat      --\n "$(tophat --version)
[ -z "$(which tophat)" ] && echo "[ERROR] no tophat detected" && exit 1
echo -e "--bowtie2     --\n "$(bowtie2 --version)
[ -z "$(which bowtie2)" ] && echo "[ERROR] no bowtie2 detected" && exit 1
echo -e "--samtools    --\n "$(samtools 2>&1 | head -n 3 | tail -n-2)
[ -z "$(which samtools)" ] && echo "[ERROR] no samtools detected" && exit 1
echo -e "--R           --\n "$(R --version | head -n 3)
[ -z "$(which R)" ] && echo "[ERROR] no R detected" && exit 1
echo -e "--picard      --\n "$(java $JAVAPARAMS -jar $PATH_PICARD/CleanSam.jar --version 2>&1)
[ ! -f $PATH_PICARD/CleanSam.jar ] && echo "[ERROR] no picard detected" && exit 1
echo -e "--samstat     --\n "$(samstat -h | head -n 2 | tail -n1)
[ -z "$(which samstat)" ] && echo "[ERROR] no samstat detected" && exit 1
echo -e "--bedtools    --\n "$(bedtools --version)
[ -z "$(which bedtools)" ] && echo "[ERROR] no bedtools detected" && exit 1


echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="parameters"

[ ! -f $f ] && echo "[ERROR] input file not found: $f" 1>&2 && exit 1

# get basename of f (samplename)
n=${f##*/}
SAMPLE=${n/%$READONE.$FASTQ/}

if [[ ! -e ${FASTA%.*}.1.bt2 ]]; then
    echo "[ERROR] Bowtie2 index not detected. Exeute bowtieIndex.sh first"
    exit 1
fi

# get info about input file
BAMFILE=$OUTDIR/../$SAMPLE$ASD.bam

#remove old files
if [ -z "$RECOVERFROM" ]; then
    if [ -d $OUTDIR ]; then rm -r $OUTDIR; fi
fi


if [ -n "$READONE" ] && [ "$READONE" == "$READTWO" ]; then
	echo "[ERROR] read1 == read2 " 1>&2 && exit 1
elif [ "$f" != "${f/%$READONE.$FASTQ/$READTWO.$FASTQ}" ] && [ -e ${f/%$READONE.$FASTQ/$READTWO.$FASTQ} ] && [ "$FORCESINGLE" = 0 ]; then
    PAIRED="1"
    f2=${f/%$READONE.$FASTQ/$READTWO.$FASTQ}
    BAM2BW_OPTION_ISPAIRED="True"
    echo "[NOTE] Paired library detected"
else
    PAIRED="0"
    BAM2BW_OPTION_ISPAIRED="False"
    echo "[NOTE] Single-Strand (unpaired) library detected"
fi

#is ziped ?
CAT="cat"
if [[ ${f##*.} == "gz" ]]; 
    then CAT="zcat"; 
elif [[ ${f##*.} == "bz2" ]]; 
    then CAT="bzcat"; 
fi

# get encoding
if [ -z "$FASTQ_PHRED" ]; then 
    FASTQ_ENCODING=$($CAT $f |  awk 'NR % 4 ==0' | python $NGSANE_BASE/tools/GuessFastqEncoding.py |  tail -n 1)
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


if [ $RNA_SEQ_LIBRARY_TYPE = "fr-unstranded" ]; then
    echo "[NOTE] make bigwigs; library is fr-unstranded "
    BAM2BW_OPTION_1="FALSE"
    BAM2BW_OPTION_2="FALSE"
elif [ $RNA_SEQ_LIBRARY_TYPE = "fr-firststrand" ]; then
    echo "[NOTE] make bigwigs; library is fr-firststrand "
    BAM2BW_OPTION_1="TRUE"
    BAM2BW_OPTION_2="TRUE"
elif [ $RNA_SEQ_LIBRARY_TYPE = "fr-secondstrand" ]; then
    echo "[NOTE] make bigwigs; library is fr-secondstrand "
    BAM2BW_OPTION_1="TRUE"
    BAM2BW_OPTION_2="FALSE"	    
fi

BIGWIGSDIR=$OUTDIR/../

mkdir -p $OUTDIR

if [ -n "$TOPHATTRANSCRIPTOMEINDEX" ]; then
    echo "[NOTE] RNAseq --transcriptome-index specified: ${TOPHATTRANSCRIPTOMEINDEX}"
    TOPHAT_TRANSCRIPTOME_PARAM="--transcriptome-index=${TOPHATTRANSCRIPTOMEINDEX}"
    PICARD_REFERENCE=${TOPHATTRANSCRIPTOMEINDEX}.fa
else
    echo "[NOTE] no --transcriptome-index specified."
    TOPHAT_TRANSCRIPTOME_PARAM=
    PICARD_REFERENCE=$FASTA
fi

THISTMP=$TMP"/"$(whoami)"/"$(echo $OUTDIR/$SAMPLE | md5sum | cut -d' ' -f1)
[ -d $THISTMP ] && rm -r $THISTMP
mkdir -p $THISTMP

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
    
    RUN_COMMAND="tophat $TOPHATADDPARAM $TOPHAT_TRANSCRIPTOME_PARAM $FASTQ_PHRED --keep-fasta-order --num-threads $CPU_TOPHAT --library-type $RNA_SEQ_LIBRARY_TYPE --rg-id $EXPID --rg-sample $PLATFORM --rg-library $LIBRARY --output-dir $OUTDIR ${FASTA%.*} $f $f2"
    echo $RUN_COMMAND && eval $RUN_COMMAND
    echo "[NOTE] tophat end $(date)"

    # mark checkpoint
    if [ -e $OUTDIR/accepted_hits.bam ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi 

################################################################################
CHECKPOINT="merge mapped and unmapped"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

    echo "[NOTE] add unmapped reads"
    if [ -f $OUTDIR/unmapped.bam ]; then
        samtools merge -@ $CPU_TOPHAT -f $THISTMP/$SAMPLE.tmp.bam $OUTDIR/accepted_hits.bam $OUTDIR/unmapped.bam
    else
        ln -s $OUTDIR/accepted_hits.bam $THISTMP/$SAMPLE.tmp.bam
    fi
    
    if [ "$PAIRED" = "1" ]; then
        # fix and sort
        echo "[NOTE] fixmate"
        RUN_COMMAND="java $JAVAPARAMS -jar $PATH_PICARD/FixMateInformation.jar \
            I=$THISTMP/$SAMPLE.tmp.bam \
            O=$THISTMP/$SAMPLE.sorted.bam \
            VALIDATION_STRINGENCY=SILENT \
            SORT_ORDER=coordinate \
            TMP_DIR=$THISTMP"
        echo $RUN_COMMAND && eval $RUN_COMMAND
    else
        # just sort
        samtools sort -@ $CPU_TOPHAT $THISTMP/$SAMPLE.tmp.bam $THISTMP/$SAMPLE.sorted
    fi
    rm $THISTMP/$SAMPLE.tmp.bam   
    
    echo "[NOTE] add read group"
    RUN_COMMAND="java $JAVAPARAMS -jar $PATH_PICARD/AddOrReplaceReadGroups.jar \
        I=$THISTMP/$SAMPLE.sorted.bam \
        O=$THISTMP/$SAMPLE.rg.bam \
        LB=$EXPID PL=Illumina PU=XXXXXX SM=$EXPID \
        VALIDATION_STRINGENCY=SILENT \
        TMP_DIR=$THISTMP"
    echo $RUN_COMMAND && eval $RUN_COMMAND
    
    rm $THISTMP/$SAMPLE.sorted.bam
        
    RUN_COMMAND="java $JAVAPARAMS -jar $PATH_PICARD/CleanSam.jar \
        I=$THISTMP/$SAMPLE.rg.bam \
        O=$BAMFILE \
        CREATE_MD5_FILE=true \
        COMPRESSION_LEVEL=9 \
        VALIDATION_STRINGENCY=SILENT \
        TMP_DIR=$THISTMP"
    echo $RUN_COMMAND && eval $RUN_COMMAND

    rm $THISTMP/$SAMPLE.rg.bam

    # mark checkpoint
    if [ -f $BAMFILE ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
fi 

################################################################################
CHECKPOINT="flagstat"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
    echo "[NOTE] samtools flagstat"
    samtools flagstat $BAMFILE > $BAMFILE.stats
    READ1=$($CAT $f | wc -l | gawk '{print int($1/4)}' )
    FASTQREADS=$READ1
    if [ -n "$f2" ]; then 
        READ2=$($CAT $f2 | wc -l | gawk '{print int($1/4)}' );
        let FASTQREADS=$READ1+$READ2
    fi
    echo $FASTQREADS" fastq reads" >> $BAMFILE.stats
    JUNCTION=$(wc -l $OUTDIR/junctions.bed | cut -d' ' -f 1)
    echo $JUNCTION" junction reads" >> $BAMFILE.stats
    ## get junction genes overlapping exons +-200bp
    
    if [ -n "$GTF" ]; then
        JUNCTGENE=$(windowBed -a $OUTDIR/junctions.bed -b $GTF -u -w 200 | wc -l | cut -d' ' -f 1)
        echo $JUNCTGENE" junction reads GTF" >> $BAMFILE.stats
    else 
        echo "0 junction reads (no gtf given)" >> $BAMFILE.stats
    fi

    # mark checkpoint
    if [ -f $BAMFILE.stats ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi 

################################################################################
CHECKPOINT="index and calculate inner distance"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
    echo "[NOTE] samtools index"
    samtools index $BAMFILE

    echo "********* calculate inner distance"
    echo "[NOTE] picard CollectMultipleMetrics"
    if [ ! -e $OUTDIR/../metrices ]; then mkdir -p $OUTDIR/../metrices ; fi
    RUN_COMMAND="java $JAVAPARAMS -jar $PATH_PICARD/CollectMultipleMetrics.jar \
        INPUT=$BAMFILE \
        REFERENCE_SEQUENCE=$PICARD_REFERENCE \
        OUTPUT=$OUTDIR/../metrices/$(basename $BAMFILE) \
        VALIDATION_STRINGENCY=SILENT \
        PROGRAM=CollectAlignmentSummaryMetrics \
        PROGRAM=CollectInsertSizeMetrics \
        PROGRAM=QualityScoreDistribution \
        TMP_DIR=$THISTMP"
    echo $RUN_COMMAND && eval $RUN_COMMAND
    
    for im in $( ls $OUTDIR/../metrices/$(basename $BAMFILE)*.pdf ); do
        convert $im ${im/pdf/jpg}
    done
   
    # mark checkpoint
    if [ -f $OUTDIR/../metrices/${BAMFILE##*/}.alignment_summary_metrics ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi 

################################################################################
CHECKPOINT="samstat"    

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

    echo "[NOTE] samstat"
    samstat $BAMFILE 2>&1 | tee | grep -v -P "Bad x in routine betai"
  
    # mark checkpoint
    if [ -f $BAMFILE.stats ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi

################################################################################
CHECKPOINT="extract mapped reads"    

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

    echo "[NOTE] extract mapped reads"
    if [ "$PAIRED" = "1" ]; then
        samtools view -@ $CPU_TOPHAT -f 3 -h -b $BAMFILE > ${BAMFILE/$ASD/$ALN}
    else
        samtools view -@ $CPU_TOPHAT -F 4 -h -b $BAMFILE > ${BAMFILE/$ASD/$ALN}
    fi
    samtools index ${BAMFILE/$ASD/$ALN}

    # mark checkpoint
    if [ -f ${BAMFILE/$ASD/$ALN} ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi

###############################################################################
CHECKPOINT="create bigwigs"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
    
    #file_arg sample_arg stranded_arg firststrand_arg paired_arg
	RUN_COMMAND="Rscript --vanilla ${NGSANE_BASE}/tools/BamToBw.R ${BAMFILE/$ASD/$ALN} ${n/%$READONE.$FASTQ/} $BAM2BW_OPTION_1 $OUTDIR/../ $BAM2BW_OPTION_2 $BAM2BW_OPTION_ISPAIRED"
	echo $RUN_COMMAND && eval $RUN_COMMAND

    # mark checkpoint
    if [ -f $OUTDIR/../${n/%$READONE.$FASTQ/.bw} ] || [ -f $OUTDIR/../${n/%$READONE.$FASTQ/_+.bw} ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi

###############################################################################
CHECKPOINT="cleanup"

[ -e ${BAMFILE/$ASD/$ALN} ] && rm ${BAMFILE/$ASD/$ALN} 
[ -e ${BAMFILE/$ASD/$ALN}.bai ] && rm ${BAMFILE/$ASD/$ALN}.bai
[ -d $THISTMP ] && rm -r $THISTMP

echo -e "\n********* $CHECKPOINT\n"
################################################################################
[ -e ${BAMFILE}.dummy ] && rm ${BAMFILE}.dummy
echo ">>>>> alignment with TopHat - FINISHED"
echo ">>>>> enddate "`date`
