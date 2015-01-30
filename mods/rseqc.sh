#!/bin/bash -e

# author: Fabian Buske
# date: November 2014

# QCVARIABLES,Resource temporarily unavailable

echo ">>>>> RSeQC"
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"

function usage {
echo -e "usage: $(basename $0) -k NGSANE -f INPUTFILE -o OUTDIR [OPTIONS]"
exit
}


if [ ! $# -gt 3 ]; then usage ; fi

FORCESINGLE=0

#INPUTS                                                                                                           
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;;     # location of the NGSANE repository                       
        -f | --bam )            shift; f=$1 ;;  # input file                                                       
        -o | --outdir )         shift; OUTDIR=$1 ;;     # output dir                                                     
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

for MODULE in $MODULE_RSEQC; do module load $MODULE; done  # save way to load modules that itself load other modules
export PATH=$PATH_RSEQC:$PATH
module list
echo "PATH=$PATH"

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--R           --\n "$(R --version | head -n 3)
[ -z "$(which R)" ] && echo "[ERROR] no R detected" && exit 1
echo -e "--Python      --\n" $(python --version)
[ -z "$(which python)" ] && echo "[ERROR] no python detected" && exit 1

NGSANE_CHECKPOINT_CHECK 
################################################################################
NGSANE_CHECKPOINT_INIT "parameters"

# get basename of input file f
n=${f##*/}
SAMPLE=${n/%$ASD.bam/}

# delete old bam files unless attempting to recover
if [ -z "$RECOVERFROM" ]; then
    [ -d $OUTDIR/$SAMPLE ] && rm -r $OUTDIR/$SAMPLE
fi

if [ -z "$GENEMODEL" ] || [ ! -f $GENEMODEL ]; then
    echo "[ERROR] no reference provided (GENEMODEL)"
    exit 1
else
    echo "[NOTE] Reference: $GENEMODEL"
fi


# check library info is set
if [ -z "$RNA_SEQ_LIBRARY_TYPE" ]; then
    echo "[ERROR] RNAseq library type not set (RNA_SEQ_LIBRARY_TYPE): either fr-unstranded or fr-firststrand"
    exit 1;
else
    echo "[NOTE] RNAseq library type: $RNA_SEQ_LIBRARY_TYPE"
fi


# unique temp folder that should be used to store temporary files
THISTMP=$TMP"/"$(whoami)"/"$(echo $OUTDIR | md5sum | cut -d' ' -f1)
mkdir -p $THISTMP

NGSANE_CHECKPOINT_CHECK 
################################################################################
NGSANE_CHECKPOINT_INIT "recall files from tape"

if [ -n "$DMGET" ]; then
    dmget -a $(dirname $GENEMODEL)/*
    dmget -a ${f}*
    dmget -a ${OUTDIR}/*
fi
    
NGSANE_CHECKPOINT_CHECK 
###############################################################################
NGSANE_CHECKPOINT_INIT "Run RSeQC - bam_stat"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

	bam_stat.py --input-file $f > $OUTDIR/$SAMPLE.bam_stat.txt 2>&1

  NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE.bam_stat.txt
fi

#################################################################################

NGSANE_CHECKPOINT_INIT "Run RSeQC - clipping profile"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

  clipping_profile.py --input-file $f --out-prefix $OUTDIR/$SAMPLE


  NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE.clipping_profile.pdf
fi

#################################################################################

NGSANE_CHECKPOINT_INIT "Run RSeQC - geneBody coverage"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

  geneBody_coverage.py --input-file $f --refgene $GENEMODEL --out-prefix $OUTDIR/$SAMPLE

  NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE.geneBodyCoverage.pdf
fi

#################################################################################
NGSANE_CHECKPOINT_INIT "Run RSeQC - infer-experiment"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

  infer_experiment.py --input-file $f --refgene $GENEMODEL > $OUTDIR/$SAMPLE.infer_experiment.txt

  NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE.infer_experiment.txt
fi

#################################################################################  
NGSANE_CHECKPOINT_INIT "Run RSeQC - innerDistancee"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

  inner_distance.py --input-file $f --refgene $GENEMODEL --out-prefix $OUTDIR/$SAMPLE 

  NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE.inner_distance_plot.pdf
fi

#################################################################################
NGSANE_CHECKPOINT_INIT "Run RSeQC - junction annotation"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

  junction_annotation.py --input-file $f --refgene $GENEMODEL --out-prefix $OUTDIR/$SAMPLE > $OUTDIR/$SAMPLE.junction.annot.txt 2>&1

  NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE.junction.annot.txt
fi

#################################################################################    
NGSANE_CHECKPOINT_INIT "Run RSeQC - junction_saturation"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

  junction_saturation.py --input-file $f --refgene $GENEMODEL --out-prefix $OUTDIR/$SAMPLE

  NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE.junctionSaturation_plot.pdf
fi

################################################################################# 
NGSANE_CHECKPOINT_INIT "Run RSeQC - read distribution"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

  read_distribution.py --input-file $f --refgene $GENEMODEL > $OUTDIR/$SAMPLE.read.distribution.txt 2>&1

  NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE.read.distribution.txt
fi

#################################################################################     
NGSANE_CHECKPOINT_INIT "Run RSeQC - read_duplication"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

  read_duplication.py --input-file $f --out-prefix $OUTDIR/$SAMPLE

  NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE.DupRate_plot.pdf
fi

################################################################################# 
NGSANE_CHECKPOINT_INIT "Run RSeQC - read_GC"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

  read_GC.py --input-file $f --out-prefix $OUTDIR/$SAMPLE

  NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE.GC_plot.pdf
fi

################################################################################# 
NGSANE_CHECKPOINT_INIT "Run RSeQC - read_NVC"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

  read_NVC.py --input-file $f --out-prefix $OUTDIR/$SAMPLE
  
  NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE.NVC_plot.pdf
fi

#################################################################################    
NGSANE_CHECKPOINT_INIT "Run RSeQC - read_quality"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

  read_quality.py --input-file $f --out-prefix $OUTDIR/$SAMPLE

  NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE.qual.heatmap.pdf
fi

#################################################################################   
NGSANE_CHECKPOINT_INIT "Run RSeQC - RPKM saturation"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

  RPKM_saturation.py --input-file $f --out-prefix $OUTDIR/$SAMPLE --refgene $GENEMODEL 
  
  NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE.saturation.pdf
fi

#################################################################################   
NGSANE_CHECKPOINT_INIT "Run RSeQC - RPKM count"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

  RPKM_count.py --input-file $f --out-prefix $OUTDIR/$SAMPLE --refgene $GENEMODEL --skip-multi-hits

  NGSANE_CHECKPOINT_CHECK
fi

################################################################################
echo ">>>>> RSeQC - FINISHED"
echo ">>>>> enddate "`date`
