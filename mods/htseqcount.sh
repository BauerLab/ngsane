#!/bin/bash -e

# Script to run HTseq-count
# It takes tophat bam files as input.
# It produces feature count files.
# author: Hugh French
# date: 2013

# messages to look out for -- relevant for the QC.sh script:
# QCVARIABLES,truncated file
# RESULTFILENAME <DIR>/<TASK>/<SAMPLE>.summary.txt

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

. $CONFIG
. ${NGSANE_BASE}/conf/header.sh
. $CONFIG

################################################################################
CHECKPOINT="programs"

# save way to load modules that itself loads other modules
hash module 2>/dev/null && for MODULE in $MODULE_HTSEQCOUNT; do module load $MODULE; done && module list 

export PATH=$PATH_HTSEQCOUNT:$PATH
echo "PATH=$PATH"
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)
[ -z "$PATH_PICARD" ] && PATH_PICARD=$(dirname $(which MarkDuplicates.jar))

echo "[NOTE] set java parameters"
JAVAPARAMS="-Xmx"$(python -c "print int($MEMORY_HTSEQCOUNT*0.8)")"g -Djava.io.tmpdir="$TMP" -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1" 
unset _JAVA_OPTIONS
echo "JAVAPARAMS "$JAVAPARAMS

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--samtools    --\n "$(samtools 2>&1 | head -n 3 | tail -n-2)
[ -z "$(which samtools)" ] && echo "[ERROR] no samtools detected" && exit 1
echo -e "--R           --\n "$(R --version | head -n 3)
[ -z "$(which R)" ] && echo "[ERROR] no R detected" && exit 1
echo -e "--bedtools    --\n "$(bedtools --version)
[ -z "$(which bedtools)" ] && echo "[ERROR] no bedtools detected" && exit 1
echo -e "--PICARD      --\n "$(java $JAVAPARAMS -jar $PATH_PICARD/FixMateInformation.jar --version 2>&1)
[ ! -f $PATH_PICARD/FixMateInformation.jar ] && echo "[ERROR] no picard detected" && exit 1
echo -e "--htSeq       --\n "$(htseq-count | tail -n 1)
[ -z "$(which htseq-count)" ] && [ -n "$GTF" ] && echo "[ERROR] no htseq-count or GTF detected" && exit 1
echo -e "--Python      --\n" $(python --version 2>&1 | tee | head -n 1 )
[ -z "$(which python)" ] && echo "[ERROR] no python detected" && exit 1
hash module 2>/dev/null && echo -e "--Python libs --\n "$(yolk -l)

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="recall files from tape"

if [ -n "$DMGET" ]; then
    dmget -a ${f}*
	dmget -a $OUTDIR/*
fi

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="parameters"

[ ! -f $f ] && echo "[ERROR] input file not found: $f" && exit 1

# get basename of f (samplename)
n=${f##*/}
SAMPLE=${n/%$ASD.bam/}

#remove old files
if [ -z "$RECOVERFROM" ]; then
    if [ -d $OUTDIR ]; then rm -r $OUTDIR; fi
fi

## GTF provided?
if [ -z "$GTF" ] || [ ! -f $GTF ]; then
    echo "[ERROR] GTF not specified or not found!"
    exit 1
else
    echo "[NOTE] GTF: $GTF"
fi

if [ ! -z "$DOCTOREDGTFSUFFIX" ]; then
    if [ ! -f ${GTF/%.gtf/$DOCTOREDGTFSUFFIX} ] ; then
        echo "[ERROR] Doctored GTF suffix specified but gtf not found: ${GTF/%.gtf/$DOCTOREDGTFSUFFIX}"
        exit 1
    else 
        echo "[NOTE] Using detected doctored GTF: ${GTF/%.gtf/$DOCTOREDGTFSUFFIX}"
        GTF=${GTF/%.gtf/$DOCTOREDGTFSUFFIX}
    fi
fi

annoF=${GTF##*/}
anno_version=${annoF%.*}


echo "[NOTE] GTF: $anno_version"


# check library info is set
if [ -z "$RNA_SEQ_LIBRARY_TYPE" ]; then
    echo "[ERROR] RNAseq library type not set (RNA_SEQ_LIBRARY_TYPE): either fr-unstranded or fr-firststrand"
    exit 1;
else
    echo "[NOTE] RNAseq library type: $RNA_SEQ_LIBRARY_TYPE"
fi


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
       HTSEQCOUNT_ADDPARAMS="--stranded=no"
elif [ "$RNA_SEQ_LIBRARY_TYPE" = "fr-firststrand" ]; then
       echo "[NOTE] library is fr-firststrand; run htseq-count stranded"
       HTSEQCOUNT_ADDPARAMS="--stranded=reverse"
elif [ "$RNA_SEQ_LIBRARY_TYPE" = "fr-secondstrand" ]; then
       echo "[NOTE] library is fr-secondstrand; run htseq-count stranded"
       HTSEQCOUNT_ADDPARAMS="--stranded=yes"
fi

if [ -z "$HTSEQCOUNT_MODES" ]; then
    echo "[ERROR] HTSEQCOUNT_MODES not defined" && exit 1
fi

if [ -z "$HTSEQCOUNT_ATTRIBUTES" ]; then
    echo "[ERROR] HTSEQCOUNT_ATTRIBUTES not defined" && exit 1
fi

mkdir -p $OUTDIR

THISTMP=$TMP"/"$(whoami)"/"$(echo $OUTDIR/$SAMPLE | md5sum | cut -d' ' -f1)
[ -d $THISTMP ] && rm -r $THISTMP
mkdir -p $THISTMP

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="fix mates"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

	if	[[ -n "$HTSEQCOUNT_UNIQUE" ]] ; then

		echo "[NOTE] Filter for uniquely mapped reads"
   		samtools view -h $f | grep -E 'NH:i:1|^@' | samtools view -b -S - > $THISTMP/$SAMPLE.tmp
		samtools sort -@ $CPU_HTSEQCOUNT -n $THISTMP/$SAMPLE.tmp $THISTMP/$SAMPLE.tmp
		rm $THISTMP/$SAMPLE.tmp
    else 
    	samtools sort -@ $CPU_HTSEQCOUNT -n $f $THISTMP/$SAMPLE.tmp
	fi	

    # fix mate pairs
    echo "[NOTE] fixmate"
    RUN_COMMAND="java $JAVAPARAMS -jar $PATH_PICARD/FixMateInformation.jar \
        I=$THISTMP/$SAMPLE.tmp.bam \
        O=$OUTDIR/$SAMPLE.fixed.bam \
        VALIDATION_STRINGENCY=SILENT \
        SORT_ORDER=coordinate \
        TMP_DIR=$THISTMP"
        echo $RUN_COMMAND && eval $RUN_COMMAND
        
	samtools index $OUTDIR/$SAMPLE.fixed.bam
    
    # cleanup
    
    [ -e $THISTMP/$SAMPLE.tmp.bam ] && rm $THISTMP/$SAMPLE.tmp.bam
    # mark checkpoint
    if [ -f $OUTDIR/$SAMPLE.fixed.bam ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
   
fi
################################################################################
CHECKPOINT="Calculate Counts"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
    cat /dev/null > $OUTDIR/GTF.summary.txt
    echo "[NOTE] Summary file - $OUTDIR/${anno_version}.summary.txt"    
    for ATTR in $HTSEQCOUNT_ATTRIBUTES; do 
        for MODE in $HTSEQCOUNT_MODES; do 
            echo "[NOTE] processing $ATTR $MODE"
            if [ "$PAIRED" = 1 ]; then 
            	samtools view -f 3 $OUTDIR/$SAMPLE.fixed.bam | htseq-count  --order=pos --idattr=$ATTR --mode=$MODE $HTSEQCOUNT_ADDPARAMS - $GTF > $THISTMP/GTF.$MODE.$ATTR.tmp
            else
            	samtools view -F 4 $OUTDIR/$SAMPLE.fixed.bam | htseq-count  --order=pos --idattr=$ATTR --mode=$MODE $HTSEQCOUNT_ADDPARAMS - $GTF > $THISTMP/GTF.$MODE.$ATTR.tmp
        	fi
            head -n-5 $THISTMP/GTF.$MODE.$ATTR.tmp > $OUTDIR/$anno_version.$MODE.$ATTR
            echo "${ATTR} ${MODE} "$(tail -n 5 $THISTMP/GTF.$MODE.$ATTR.tmp | sed 's/\s\+/ /g' | tr '\n' ' ') >> $OUTDIR/GTF.summary.txt
            rm $THISTMP/GTF.$MODE.$ATTR.tmp
            
        done
    done
    
    # mark checkpoint
    if [ -f $OUTDIR/$anno_version.$MODE.$ATTR ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
fi
	

################################################################################
CHECKPOINT="summarize"

cat $OUTDIR/GTF.summary.txt | awk '{print "all",$0}' > ${OUTDIR}/${SAMPLE}.${anno_version}.summary.txt
   
echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="cleanup"    

if [ -z "$HTSEQCOUNT_KEEPBAMS"]; then
    [ -e $OUTDIR/$SAMPLE.fixed.bam ] && rm $OUTDIR/$SAMPLE.fixed.bam
    [ -e $OUTDIR/$SAMPLE.fixed.bam.bai ] && rm $OUTDIR/$SAMPLE.fixed.bam.bai
fi
[ -e $OUTDIR/GTF.summary.txt ] && rm $OUTDIR/GTF.summary.txt
[ -e $OUTDIR/..//${SAMPLE}.summary.txt.dummy ] && rm $OUTDIR/..//${SAMPLE}.summary.txt.dummy


echo -e "\n********* $CHECKPOINT\n"
################################################################################

echo ">>>>> feature counting with HTSEQ-COUNT - FINISHED"
echo ">>>>> enddate "`date`
