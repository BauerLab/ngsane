#!/bin/bash -e

# Script to run CUFFLINKS program
# It takes tophat bam files as input.
# author: Chikako Ragan, Denis Bauer
# date: Jan. 2011
# modified by Fabian Buske and Hugh French
# date: 2013-

# messages to look out for -- relevant for the QC.sh script:
# QCVARIABLES,truncated file
# RESULTFILENAME <DIR>/<TASK>/<SAMPLE>_transcripts.gtf

echo ">>>>> transcript assembly with cufflinks "
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"


function usage {
echo -e "usage: $(basename $0) -k NGSANE -f BAM -o OUTDIR [OPTIONS]"
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
hash module 2>/dev/null && for MODULE in $MODULE_CUFFLINKS; do module load $MODULE; done && module list 

export PATH=$PATH_CUFFLINKS:$PATH
echo "PATH=$PATH"
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--cufflinks   --\n "$(cufflinks 2>&1 | tee | head -n 2 )
[ -z "$(which cufflinks)" ] && echo "[ERROR] no cufflinks detected" && exit 1

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="parameters"

[ ! -f $f ] && echo "[ERROR] input file not found: $f" && exit 1

# get basename of f (samplename)
n=${f##*/}

#remove old files unless recovering
if [ -z "$RECOVERFROM" ]; then
    if [ -e $OUTDIR/../${n/%.$ASD.bam/_transcripts.gtf} ]; then rm $OUTDIR/../${n/%.$ASD.bam/_transcripts.gtf}; fi
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

# check library info is set
if [ -z "$RNA_SEQ_LIBRARY_TYPE" ]; then
    echo "[ERROR] RNAseq library type not set (RNA_SEQ_LIBRARY_TYPE): either fr-unstranded or fr-firststrand"
    exit 1;
else
    echo "[NOTE] RNAseq library type: $RNA_SEQ_LIBRARY_TYPE"
fi

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="recall files from tape"

if [ -n "$DMGET" ]; then
    dmget -a $(dirname $FASTA)/*
    dmget -a ${f}*
    dmget -a ${OUTDIR}/*
fi

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="run cufflinks"    

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
#    echo ">>>>> from $f to $OUTPUT"
    echo "[NOTE] cufflinks $(date)"

    ## add GTF file if present
    if [ -n "$GTF" ]; then 
        echo "[NOTE] execute cufflinks with guide (GTF provided)"
        RUN_COMMAND="cufflinks --quiet $CUFFLINKSADDPARAM --GTF-guide $GTF -p $CPU_CUFFLINKS --library-type $RNA_SEQ_LIBRARY_TYPE -o $OUTDIR $f"

    else
        # non reference guided
        echo "[NOTE] de novo run (no GTF provided)"
        RUN_COMMAND="cufflinks --quiet $CUFFLINKSADDPARAM --frag-bias-correct $FASTA -p $CPU_CUFFLINKS --library-type $RNA_SEQ_LIBRARY_TYPE -o $OUTDIR $f"
    fi
    echo $RUN_COMMAND && eval $RUN_COMMAND

    CURDIR=$(pwd) 
    cd $OUTDIR/.. 
    ln -f -s ${n/%.$ASD.bam/}/transcripts.gtf ${n/%.$ASD.bam/_transcripts.gtf} 
    cd $CURDIR
    
    echo "[NOTE] cufflinks end $(date)"

    # mark checkpoint
    if [ -e $OUTDIR/../${n/%.$ASD.bam/_transcripts.gtf} ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi

################################################################################
CHECKPOINT="statistics"    

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

    SUMMARYFILE=$OUTDIR/../${n/%.$ASD.bam/.summary.txt}
    cat /dev/null > $SUMMARYFILE
    
    for i in $(ls $OUTDIR/*.fpkm_tracking); do
        echo "${i##*/} $(cat $i | tail -n+2 | wc -l ); $(cat $i | cut -f 13 | tail -n+2 | sort | uniq -c | tr '\n' ';')" >> $SUMMARYFILE
    done
    
    for i in $(ls $OUTDIR/*.gtf); do
        if [ -e $i ]; then
            echo "${i##*/} $(cat $i | wc -l )" >> $SUMMARYFILE
        fi
    done
    
    # mark checkpoint
    if [ -e $SUMMARYFILE ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi


fi
################################################################################
[ -e $OUTDIR/../${n/%.$ASD.bam/_transcripts.gtf}.dummy ] && rm $OUTDIR/../${n/%.$ASD.bam/_transcripts.gtf}.dummy
echo ">>>>> transcript assembly with cufflinks - FINISHED"
echo ">>>>> enddate "`date`
