#!/bin/bash

# cufflinks calling script
# author: Denis C. Bauer
# date: March.2011

# messages to look out for -- relevant for the QC.sh script:
# QCVARIABLES,

echo ">>>>> differential expression with cuffdiff"
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"


function usage {
echo -e "usage: $(basename $0) -k NGSANE -f  -r REFERENCE -o OUTDIR [OPTIONS]

Script running read mapping for single and paired DNA reads from fastq files
It expects a fastq file, pairdend, reference genome  as input and 
It runs BWA, converts the output to .bam files, adds header information and
writes the coverage information for IGV.

required:
  -k | --toolkit <path>     location of the NGSANE repository 
  -b | --basename <b1[,b2]> basename comma separated
  -r | --reference <file>   reference genome
  -o | --outdir <path>      output dir

options:
  -t | --threads <nr>       number of CPUs to use (default: 1)
  -i | --rgid <name>        read group identifier RD ID (default: exp)
  -l | --rglb <name>        read group library RD LB (default: qbi)
  -p | --rgpl <name>        read group platform RD PL (default: illumna)
  -s | --rgsi <name>        read group sample RG SM prefac (default: )
  -R | --region <ps>        region of specific interest, e.g. targeted reseq
                             format chr:pos-pos
"
exit
}


if [ ! $# -gt 3 ]; then usage ; fi

#DEFAULTS
THREADS=1

#INPUTS
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the NGSANE repository
        -t | --threads )        shift; THREADS=$1 ;; # number of CPUs to use
        -b | --basename )       shift; fs=$1 ;; # basename
        -r | --reference )      shift; FASTA=$1 ;; # reference genome
        -o | --outdir )         shift; OUTDIR=$1 ;; # output dir
        -a | --annot )          shift; REFSEQGTF=$1 ;; # refseq annotation
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
hash module 2>/dev/null && for MODULE in $MODULE_CUFFDIFF; do module load $MODULE; done && module list 

export PATH=$PATH_CUFFDIFF:$PATH
echo "PATH=$PATH"

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--cuffdiff    --\n "$(cuffdiff 2>&1 | head -n 1 )
[ -z "$(which cuffdiff)" ] && echo "[ERROR] no cuffdiff detected" && exit 1
echo -e "--cuffcompare --\n "$(cuffcompare 2>&1 | head -n 1 )
[ -z "$(which cuffcompare)" ] && echo "[ERROR] no cuffcompare detected" && exit 1

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "parameters"

# delete old bam files unless attempting to recover
if [ -z "$NGSANE_RECOVERFROM" ]; then
    if [ -d $OUTDIR ]; then rm -r $OUTDIR; fi
    mkdir -p $OUTDIR
fi

# get basename of f (samplename)
n=${fs/,/:}
O=${OUT/$n/}
CUFOUT=${O/$TASK_CUFFDIFF/Run\/$TASK_CUFF/}
TOPHATOUT=${O/$TASK_CUFFDIFF/Run\/$TASK_TOPHAT/}

CUFGTFS=""
TOPHATBAM=""
for v in ${fs//,/ }; do
    f=$(basename $v)
    CUFGTFS=$CUFGTFS" "$CUFOUT/$f/transcripts.gtf
    TOPHATBAM=$TOPHATBAM" "$TOPHATOUT/$f/accepted_hits.bam
done

cd $OUTDIR/
GTF=$REFSEQGTF
    
################################################################################
NGSANE_CHECKPOINT_INIT "recall files from tape"

if [ -n "$DMGET" ]; then
    dmget -a $(dirname $TOPHATOUT)/*
    dmget -a ${OUTDIR}/*
fi

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "check reference GTF"
    
if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    # Which transcript reference?
    if [ -n $GTF ]; then
        echo "[NOTE] compare to oneanother"
        
        cuffcompare -o comp.txt $CUFGTFS
        GTF=$OUTDIR/comp.combined.gtf
        
        # mark checkpoint
        [[ -s $OUTDIR/comp.combined.gtf ]] && NGSANE_CHECKPOINT_CHECK
    else
    
        echo "[NOTE] skip checking reference GTF"
        # mark checkpoint
        NGSANE_CHECKPOINT_CHECK
    fi    
fi

################################################################################
NGSANE_CHECKPOINT_INIT "run cuff diff"
    
if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    cuffdiff -r $FASTA -p $CPU_CUFFDIFF -o $OUTDIR $GTF $TOPHATBAM

    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK
fi 

################################################################################
echo ">>>>> differential expression with cuffdiff - FINISHED"
echo ">>>>> enddate "`date`
