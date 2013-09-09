#!/bin/bash -e

# Script to run HTseq-count
# It takes tophat bam files as input.
# It produces feature count files.
# author: Hugh French
# date: 2013

# messages to look out for -- relevant for the QC.sh script:
# QCVARIABLES,truncated file
# RESULTFILENAME <SAMPLE>_masked

echo ">>>>> feature counting with HTSEQ-COUNT "
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"


function usage {
echo -e "usage: $(basename $0) -k NGSANE -o OUTDIR [OPTIONS]"
exit
}

if [ ! $# -gt 3 ]; then usage ; fi

#INPUTS
while [ "$1" != "" ]; do
	case $1 in
	-k | toolkit )          shift; CONFIG=$1 ;; # ENSURE NO VARIABLE NAMES FROM CONFIG
	-f | --bam )            shift; f=$1 ;; # fastq file
	-o | --outdir )         shift; OUTDIR=$1 ;; # output dir
    --recover-from )        shift; RECOVERFROM=$1 ;; # attempt to recover from log file
	-h | --help )           usage ;;
	* )                     echo "dont understand $1"
	esac
	shift
done


#PROGRAMS (note, both configs are necessary to overwrite the default, here:e.g.  TASKTOPHAT)
. $CONFIG
. ${NGSANE_BASE}/conf/header.sh
. $CONFIG

################################################################################
CHECKPOINT="programs"

for MODULE in $MODULE_HTSEQCOUNT; do module load $MODULE; done  # save way to load modules that itself load other modules
export PATH=$PATH_HTSEQCOUNT:$PATH
module list
echo "PATH=$PATH"
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--samtools    --\n "$(samtools 2>&1 | head -n 3 | tail -n-2)
[ -z "$(which samtools)" ] && echo "[ERROR] no samtools detected" && exit 1
echo -e "--R           --\n "$(R --version | head -n 3)
[ -z "$(which R)" ] && echo "[ERROR] no R detected" && exit 1
echo -e "--bedtools    --\n "$(bedtools --version)
[ -z "$(which bedtools)" ] && echo "[ERROR] no bedtools detected" && exit 1
echo -e "--htSeq       --\n "$(htseq-count | tail -n 1)
[ -z "$(which htseq-count)" ] && [ -n "$GENCODEGTF" ] && echo "[ERROR] no htseq-count or GENCODEGTF detected" && exit 1
echo -e "--Python      --\n" $(python --version 2>&1 | tee | head -n 1 )
[ -z "$(which python)" ] && echo "[ERROR] no python detected" && exit 1
[ hash yolk ] && echo -e "--Python libs --\n "$(yolk -l)

echo "[NOTE] set java parameters"
JAVAPARAMS="-Xmx"$(python -c "print int($MEMORY_TOPHAT*0.8)")"g -Djava.io.tmpdir="$TMP"  -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1" 
unset _JAVA_OPTIONS
echo "JAVAPARAMS "$JAVAPARAMS

echo -e "\n********* $CHECKPOINT"
################################################################################
CHECKPOINT="recall files from tape"

if [ -n "$DMGET" ]; then
    dmget -a $(dirname $FASTA)/*
    dmget -a ${f}*
fi

echo -e "\n********* $CHECKPOINT"
################################################################################
CHECKPOINT="parameters"

[ ! -f $f ] && echo "[ERROR] input file not found: $f" && exit 1

# get basename of f (samplename)
n=${f##*/}

# get info about input file
FASTASUFFIX=${FASTA##*.}

#remove old files
if [ -z "$RECOVERFROM" ]; then
    if [ -d $OUTDIR ]; then rm -r $OUTDIR; fi
fi

## GTF provided?
if [ -n "$GENCODEGTF" ]; then
    echo "[NOTE] Gencode GTF: $GENCODEGTF"
    if [ ! -f $GENCODEGTF ]; then
        echo "[ERROR] GENCODE GTF specified but not found!"
        exit 1
    fi 
elif [ -n "$REFSEQGTF" ]; then
    echo "[NOTE] Refseq GTF: $REFSEQGTF"
    if [ ! -f $REFSEQGTF ]; then
        echo "[ERROR] REFSEQ GTF specified but not found!"
        exit 1
    fi
fi

if [ -n "$REFSEQGTF" ] && [ -n "$GENCODEGTF" ]; then
    echo "[WARN] GENCODE and REFSEQ GTF found. GENCODE takes preference."
fi
if [ ! -z "$DOCTOREDGTFSUFFIX" ]; then
    if [ ! -f ${GENCODEGTF/%.gtf/$DOCTOREDGTFSUFFIX} ] ; then
        echo "[ERROR] Doctored GTF suffix specified but gtf not found: ${GENCODEGTF/%.gtf/$DOCTOREDGTFSUFFIX}"
        exit 1
    else 
        echo "[NOTE] Doctored GTF: ${GENCODEGTF/%.gtf/$DOCTOREDGTFSUFFIX}"
    fi
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

annoF=${GENCODEGTF##*/}
anno_version=${annoF%.*}

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
RPKMSSDIR=$OUTDIR/../


# run flagstat if no stats available for bam file
[ ! -e $f.stats ] && samtools flagstat > $f.stats
# check "paired in sequencing" entry to detect library
if [[ $(cat $f.stats | head -n 4 | tail -n 1 | cut -d' ' -f 1) -gt 0 ]]; then
    PAIRED=1
    echo "[NOTE] paired library detected"
else 
    PAIRED=0
    echo "[NOTE] single-end library detected"
fi

if [ "$RNA_SEQ_LIBRARY_TYPE" = "fr-unstranded" ]; then
       echo "[NOTE] library is fr-unstranded; do not run htseq-count stranded"
       HT_SEQ_OPTIONS="--stranded=no"
elif [ "$RNA_SEQ_LIBRARY_TYPE" = "fr-firststrand" ]; then
       echo "[NOTE] library is fr-firststrand; run htseq-count stranded"
       HT_SEQ_OPTIONS="--stranded=reverse"
elif [ "$RNA_SEQ_LIBRARY_TYPE" = "fr-secondstrand" ]; then
       echo "[NOTE] library is fr-secondstrand; run htseq-count stranded"
       HT_SEQ_OPTIONS="--stranded=yes"
fi


mkdir -p $OUTDIR

echo -e "\n********* $CHECKPOINT"
################################################################################
CHECKPOINT="run htseq-count"    
if [[ -n "$RECOVERFROM" ]] && [[ $(grep "********* $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
	## htseq-count 

	samtools sort -n $f $OUTDIR/${n/%.$ASD.bam/.tmp}
	samtools fixmate $OUTDIR/${n/%.$ASD.bam/.tmp.bam} $OUTDIR/${n}

	samtools view $OUTDIR/${n} -f 3 | htseq-count --quiet $HT_SEQ_OPTIONS - $GENCODEGTF | grep ENSG > $OUTDIR/${anno_version}.gene
	samtools view $OUTDIR/${n} -f 3 | htseq-count --quiet --mode=intersection-strict $HT_SEQ_OPTIONS - $GENCODEGTF | grep ENSG > $OUTDIR/${anno_version}_strict.gene
	samtools view $OUTDIR/${n} -f 3 | htseq-count --quiet --mode=intersection-nonempty $HT_SEQ_OPTIONS - $GENCODEGTF | grep ENSG > $OUTDIR/${anno_version}_nonempty.gene
#	samtools view $OUTDIR/${n}  | htseq-count --quiet --idattr="transcript_id" $HT_SEQ_OPTIONS - $GENCODEGTF | grep ENST > $OUTDIR/${anno_version}.transcript
    
    # mark checkpoint
    [ -f $OUTDIR/${anno_version}.gene ] && echo -e "\n********* $CHECKPOINT" && unset RECOVERFROM    
fi

################################################################################
CHECKPOINT="Create bigwigs"    
if [[ -n "$RECOVERFROM" ]] && [[ $(grep "********* $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

	#make a paired only (f -3 ) bam so bigwigs are comparable to counts.
	samtools view -f 3 -h -b $f > $OUTDIR/accepted_hits_f3.bam
	
    #file_arg sample_arg stranded_arg firststrand_arg paired_arg
    Rscript --vanilla ${NGSANE_BASE}/tools/BamToBw.R $OUTDIR/accepted_hits_f3.bam ${n/%.$ASD.bam/} $BAM2BW_OPTION_1 $BIGWIGSDIR $BAM2BW_OPTION_2
	
    # mark checkpoint
    [ -f ${BIGWIGSDIR}/${n/%.$ASD.bam/.bw} ] && echo -e "\n********* $CHECKPOINT" && unset RECOVERFROM    
fi
	
################################################################################
CHECKPOINT="calculate RPKMs per Gencode Gene"    
if [[ -n "$RECOVERFROM" ]] && [[ $(grep "********* $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
    echo "[NOTE] Gencode RPKM calculation"
	
    Rscript --vanilla ${NGSANE_BASE}/tools/CalcGencodeGeneRPKM.R $GENCODEGTF $OUTDIR/${anno_version}.gene $RPKMSSDIR/${n/%.$ASD.bam/_gene} ${anno_version}

    echo "[NOTE] Create filtered bamfile"
	
    ##remove r_RNA and create counts.
	python ${NGSANE_BASE}/tools/extractFeature.py -f $GENCODEGTF --keep rRNA Mt_tRNA Mt_rRNA tRNA rRNA_pseudogene tRNA_pseudogene Mt_tRNA_pseudogene Mt_rRNA_pseudogene > $OUTDIR/mask.gff
	python ${NGSANE_BASE}/tools/extractFeature.py -f $GENCODEGTF --keep RNA18S5 RNA28S5 -l 17 >> $OUTDIR/mask.gff
	        
    intersectBed -v -abam $f -b $OUTDIR/mask.gff > $OUTDIR/tophat_aligned_reads_masked.bam    
	    
    samtools index $OUTDIR/tophat_aligned_reads_masked.bam

    rm $OUTDIR/mask.gff
    
    samtools sort -n $OUTDIR/tophat_aligned_reads_masked.bam $OUTDIR/tophat_aligned_reads_masked_sorted.tmp
    samtools fixmate $OUTDIR/tophat_aligned_reads_masked_sorted.tmp.bam $OUTDIR/tophat_aligned_reads_masked_sorted.bam
    rm $OUTDIR/tophat_aligned_reads_masked_sorted.tmp.bam
	
    samtools view $OUTDIR/tophat_aligned_reads_masked_sorted.bam -f 3  | htseq-count --quiet $HT_SEQ_OPTIONS - $GENCODEGTF | grep ENSG > $OUTDIR/${anno_version}_masked.gene
    # add intersect
    samtools view $OUTDIR/tophat_aligned_reads_masked_sorted.bam -f 3  | htseq-count --quiet --mode intersection-strict $HT_SEQ_OPTIONS - $GENCODEGTF | grep ENSG > $OUTDIR/${anno_version}_masked_strict.gene
    samtools view $OUTDIR/tophat_aligned_reads_masked_sorted.bam -f 3  | htseq-count --quiet --mode intersection-nonempty $HT_SEQ_OPTIONS - $GENCODEGTF | grep ENSG > $OUTDIR/${anno_version}_masked_nonempty.gene
    
    #  samtools view $OUTDIR/tophat_aligned_reads_masked_sorted.bam  | htseq-count --quiet --idattr="transcript_id" $HT_SEQ_OPTIONS - $GENCODEGTF | grep ENST > $OUTDIR/${anno_version}_masked.transcript

    echo "[NOTE] calculate RPKMs per Gencode Gene masked"

    Rscript --vanilla ${NGSANE_BASE}/tools/CalcGencodeGeneRPKM.R $GENCODEGTF $OUTDIR/${anno_version}_masked.gene $RPKMSSDIR/${n/%.$ASD.bam/_gene_masked} ${anno_version}

    rm $OUTDIR/tophat_aligned_reads_masked_sorted.bam

    rm $OUTDIR/${n}

    #file_arg sample_arg stranded_arg firststrand_arg paired_arg
    Rscript --vanilla ${NGSANE_BASE}/tools/BamToBw.R $OUTDIR/accepted_hits_f3.bam ${n/%.$ASD.bam/}_masked $BAM2BW_OPTION_1 $BIGWIGSDIR $BAM2BW_OPTION_2

    # mark checkpoint
    [ -f $OUTDIR/${n/%.$ASD.bam/_masked} ] && echo -e "\n********* $CHECKPOINT" && unset RECOVERFROM    
fi

################################################################################
[ -e $OUTDIR/../${n/%.$ASD.bam/_masked}.dummy ] && rm $OUTDIR/../${n/%.$ASD.bam/_masked}.dummy
echo ">>>>> feature counting with HTSEQ-COUNT - FINISHED"
echo ">>>>> enddate "`date`
