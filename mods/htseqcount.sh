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
[ -z "$(which htseq-count)" ] && [ -n "$GTF" ] && echo "[ERROR] no htseq-count or GTF detected" && exit 1
echo -e "--Python      --\n" $(python --version 2>&1 | tee | head -n 1 )
[ -z "$(which python)" ] && echo "[ERROR] no python detected" && exit 1
[ hash yolk ] && echo -e "--Python libs --\n "$(yolk -l)

echo "[NOTE] set java parameters"
JAVAPARAMS="-Xmx"$(python -c "print int($MEMORY_TOPHAT*0.8)")"g -Djava.io.tmpdir="$TMP"  -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1" 
unset _JAVA_OPTIONS
echo "JAVAPARAMS "$JAVAPARAMS

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="recall files from tape"

if [ -n "$DMGET" ]; then
    dmget -a ${f}*
fi

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="parameters"

[ ! -f $f ] && echo "[ERROR] input file not found: $f" && exit 1

# get basename of f (samplename)
n=${f##*/}

#remove old files
if [ -z "$RECOVERFROM" ]; then
    if [ -d $OUTDIR ]; then rm -r $OUTDIR; fi
fi

## GTF provided?
if [ -n "$GTF" ]; then
    echo "[NOTE] Gencode GTF: $GTF"
    if [ ! -f $GTF ]; then
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

if [ -n "$REFSEQGTF" ] && [ -n "$GTF" ]; then
    echo "[WARN] GENCODE and REFSEQ GTF found. GENCODE takes preference."
fi
if [ ! -z "$DOCTOREDGTFSUFFIX" ]; then
    if [ ! -f ${GTF/%.gtf/$DOCTOREDGTFSUFFIX} ] ; then
        echo "[ERROR] Doctored GTF suffix specified but gtf not found: ${GTF/%.gtf/$DOCTOREDGTFSUFFIX}"
        exit 1
    else 
        echo "[NOTE] Doctored GTF: ${GTF/%.gtf/$DOCTOREDGTFSUFFIX}"
    fi
fi

# check library info is set
if [ -z "$RNA_SEQ_LIBRARY_TYPE" ]; then
    echo "[ERROR] RNAseq library type not set (RNA_SEQ_LIBRARY_TYPE): either fr-unstranded or fr-firststrand"
    exit 1;
else
    echo "[NOTE] RNAseq library type: $RNA_SEQ_LIBRARY_TYPE"
fi

annoF=${GTF##*/}
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

if [ -z "$HTSEQCOUNT_MODES" ]; then
    echo "[ERROR] HTSEQCOUNT_MODES not defined" && exit 1
fi

if [ -z "$HTSEQCOUNT_ATTRIBUTES" ]; then
    echo "[ERROR] HTSEQCOUNT_ATTRIBUTES not defined" && exit 1
fi

mkdir -p $OUTDIR

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="run htseq-count"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
    cat /dev/null > $OUTDIR/${anno_version}.summary.txt

	samtools sort -n $f $OUTDIR/${n/%.$ASD.bam/.tmp}
	samtools fixmate $OUTDIR/${n/%.$ASD.bam/.tmp.bam} $OUTDIR/${n}
    
    for ATTR in $HTSEQCOUNT_ATTRIBUTES; do 
        for MODE in $HTSEQCOUNT_MODES; do 
            echo $ATTR $MODE
        	samtools view $OUTDIR/${n} -f 3 | htseq-count --quiet --idattr=$ATTR --mode=$MODE $HT_SEQ_OPTIONS - $GTF > $OUTDIR/${anno_version}.$MODE.$ATTR.tmp
            head -n-5 $OUTDIR/${anno_version}.$MODE.$ATTR.tmp > $OUTDIR/${anno_version}.$MODE.$ATTR
            echo "${ATTR} ${MODE} "$(tail -n 5 $OUTDIR/${anno_version}.$MODE.$ATTR.tmp | sed 's/\s\+/ /g' | tr '\n' ' ') >> $OUTDIR/${anno_version}.summary.txt
            rm $OUTDIR/${anno_version}.$MODE.$ATTR.tmp
        done
    done
    
    Rscript --vanilla ${NGSANE_BASE}/tools/CalcGencodeGeneRPKM.R $GTF $OUTDIR/${anno_version}.union.gene_id $RPKMSSDIR/${n/%.$ASD.bam/_gene_} ${anno_version}
#    Rscript --vanilla ${NGSANE_BASE}/tools/CalcGencodeGeneRPKM.R $GTF $OUTDIR/${anno_version}.union.transcript_id $RPKMSSDIR/${n/%.$ASD.bam/_transcript_} ${anno_version}

    # mark checkpoint
    if [ -f $OUTDIR/${anno_version}.union.gene_id ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
fi
################################################################################
CHECKPOINT="mask GTF"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
    echo "[NOTE] Create filtered bamfile"
	
    ##remove r_RNA and create counts.
	python ${NGSANE_BASE}/tools/extractFeature.py -f $GTF --keep rRNA Mt_tRNA Mt_rRNA tRNA rRNA_pseudogene tRNA_pseudogene Mt_tRNA_pseudogene Mt_rRNA_pseudogene > $OUTDIR/mask.gff
	python ${NGSANE_BASE}/tools/extractFeature.py -f $GTF --keep RNA18S5 RNA28S5 -l 17 >> $OUTDIR/mask.gff
	        
    intersectBed -v -abam $f -b $OUTDIR/mask.gff > $OUTDIR/${n/%.$ASD.bam/.$ALN.masked.bam}    
	    
    samtools index $OUTDIR/${n/%.$ASD.bam/.$ALN.masked.bam}

    [ -e $OUTDIR/mask.gff ] && rm $OUTDIR/mask.gff
    
    samtools sort -n $OUTDIR/${n/%.$ASD.bam/.$ALN.masked.bam} $OUTDIR/${n/%.$ASD.bam/.$ASD.masked.tmp}
    samtools fixmate $OUTDIR/${n/%.$ASD.bam/.$ASD.masked.tmp}.bam $OUTDIR/${n/%.$ASD.bam/.$ASD.masked.bam}

    [ -e $OUTDIR/${n/%.$ASD.bam/.$ASD.masked.tmp}.bam ] && rm $OUTDIR/${n/%.$ASD.bam/.$ASD.masked.tmp}.bam
	
    # mark checkpoint
    if [ -f $OUTDIR/${n/%.$ASD.bam/.$ASD.masked.bam} ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
   
fi	
################################################################################
CHECKPOINT="calculate RPKMs per Gencode Gene"    

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
    echo "[NOTE] Gencode RPKM calculation"
	
	cat /dev/null > $OUTDIR/${anno_version}_masked.summary.txt
	
    for ATTR in $HTSEQCOUNT_ATTRIBUTES; do 
        for MODE in $HTSEQCOUNT_MODES; do 
            echo $ATTR $MODE
            samtools view $OUTDIR/${n/%.$ASD.bam/.$ASD.masked.bam} -f 3 | htseq-count --quiet --idattr=$ATTR --mode=$MODE $HT_SEQ_OPTIONS - $GTF > $OUTDIR/${anno_version}_masked.$MODE.$ATTR.tmp
            head -n-5 $OUTDIR/${anno_version}_masked.$MODE.$ATTR.tmp > $OUTDIR/${anno_version}_masked.$MODE.$ATTR
            echo "${ATTR} ${MODE} "$(tail -n 5 $OUTDIR/${anno_version}_masked.$MODE.$ATTR.tmp | sed 's/\s\+/ /g' | tr '\n' ' ') >> $OUTDIR/${anno_version}_masked.summary.txt
            rm $OUTDIR/${anno_version}_masked.$MODE.$ATTR.tmp
        done
    done
    
    echo "[NOTE] calculate RPKMs per Gencode Gene masked"

    Rscript --vanilla ${NGSANE_BASE}/tools/CalcGencodeGeneRPKM.R $GTF $OUTDIR/${anno_version}_masked.union.gene_id $RPKMSSDIR/${n/%.$ASD.bam/_gene_masked_} ${anno_version} 
#    Rscript --vanilla ${NGSANE_BASE}/tools/CalcGencodeGeneRPKM.R $GTF $OUTDIR/${anno_version}_masked.union.transcript_id $RPKMSSDIR/${n/%.$ASD.bam/_transcript_masked_} ${anno_version} 

    # mark checkpoint
    if [ -f $RPKMSSDIR/${n/%.$ASD.bam/_gene_RPKM}${anno_version}.csv ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
   
fi

################################################################################
CHECKPOINT="create bigwigs"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"

	#make a paired only (f -3 ) bam so bigwigs are comparable to counts.
	samtools view -f 3 -h -b $f > $OUTDIR/${n/%.$ASD.bam/.$ALN.bam}
	
    #file_arg sample_arg stranded_arg firststrand_arg paired_arg
    Rscript --vanilla ${NGSANE_BASE}/tools/BamToBw.R $OUTDIR/${n/%.$ASD.bam/.$ALN.bam} ${n/%.$ASD.bam/} $BAM2BW_OPTION_1 $BIGWIGSDIR $BAM2BW_OPTION_2

    #file_arg sample_arg stranded_arg firststrand_arg paired_arg
    Rscript --vanilla ${NGSANE_BASE}/tools/BamToBw.R $OUTDIR/${n/%.$ASD.bam/.$ALN.bam} ${n/%.$ASD.bam/}_masked $BAM2BW_OPTION_1 $BIGWIGSDIR $BAM2BW_OPTION_2

    [ -e $OUTDIR/${n/%.$ASD.bam/.$ALN.bam} ] && rm $OUTDIR/${n/%.$ASD.bam/.$ALN.bam}
    
    # mark checkpoint
    if [ -f ${BIGWIGSDIR}/${n/%.$ASD.bam/.bw} ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi

################################################################################
CHECKPOINT="summarize"

cat $OUTDIR/${anno_version}.summary.txt | awk '{print all,$0}' > $OUTDIR/../${n}_${anno_version}.summary.txt
cat $OUTDIR/${anno_version}_masked.summary.txt | awk '{print masked,$0}' >> $OUTDIR/../${n}_${anno_version}.summary.txt
   
echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="cleanup"    

[ -e $OUTDIR/${n} ] && rm $OUTDIR/${n}
[ -e $OUTDIR/${n/%.$ASD.bam/.$ASD.masked.bam} ] && rm $OUTDIR/${n/%.$ASD.bam/.$ASD.masked.bam}
[ -e $OUTDIR/${anno_version}.summary.txt ] && rm $OUTDIR/${anno_version}.summary.txt
[ -e $OUTDIR/${anno_version}_masked.summary.txt ] && rm $OUTDIR/${anno_version}_masked.summary.txt

echo -e "\n********* $CHECKPOINT\n"
################################################################################
[ -e $OUTDIR/../${n/%.$ASD.bam/_masked}.dummy ] && rm $OUTDIR/../${n/%.$ASD.bam/_masked}.dummy
echo ">>>>> feature counting with HTSEQ-COUNT - FINISHED"
echo ">>>>> enddate "`date`
