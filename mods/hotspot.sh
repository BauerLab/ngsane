#!/bin/bash -e

# DNase-Seq processing using hotspot
# author: Fabian Buske
# date: January 1914

# QCVARIABLES,Resource temporarily unavailable
# RESULTFILENAME <DIR>/<TASK>/<SAMPLE>.narrowPeak

echo ">>>>> Peak calling with hotspot"
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

#INPUTS                                                                                                           
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;;     # location of the NGSANE repository                       
        -f | --file )           shift; INPUTFILE=$1 ;;  # input file                                                       
        -o | --outdir )         shift; OUTDIR=$1 ;;     # output dir                                                     
        --recover-from )        shift; NGSANE_RECOVERFROM=$1 ;; # attempt to recover from log file
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

# save way to load modules that itself loads other modules
hash module 2>/dev/null && for MODULE in $MODULE_HOTSPOT; do module load $MODULE; done && module list 

export PATH=$PATH_HOTSPOT:$PATH
echo "PATH=$PATH"

echo "[NOTE] set java parameters"
JAVAPARAMS="-Xmx"$(python -c "print int($MEMORY_HOTSPOT*0.75)")"g -Djava.io.tmpdir="$TMP" -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1" 
unset _JAVA_OPTIONS
echo "JAVAPARAMS "$JAVAPARAMS

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--JAVA        --\n" $(java -Xmx200m -version 2>&1)
[ -z "$(which java)" ] && echo "[ERROR] no java detected" && exit 1
echo -e "--Python      --\n" $(python --version 2>&1 | tee | head -n 1 )
[ -z "$(which python)" ] && echo "[ERROR] no python detected" && exit 1
hash module 2>/dev/null && echo -e "--Python libs --\n "$(yolk -l)
echo -e "--bedtools --\n "$(bedtools --version)
[ -z "$(which bedtools)" ] && echo "[ERROR] no bedtools detected" && exit 1
echo -e "--bedops --\n "$(bedops --version 2>&1 | tee | head -n 3 | tail -n 1)
[ -z "$(which bedops)" ] && echo "[ERROR] no bedtools detected" && exit 1
echo -e "--hotspot        --\n "$(hotspot -v | head -n 1)
[ -z "$(which hotspot)" ] && echo "[ERROR] no hotspot detected" && exit 1
echo -e "--R           --\n "$(R --version | head -n 3)
[ -z "$(which R)" ] && echo "[ERROR] no R detected" && exit 1
echo -e "--bedToBigBed --\n "$(bedToBigBed 2>&1 | tee | head -n 1 )
[ -z "$(which bedToBigBed)" ] && echo "[WARN] bedToBigBed not detected, cannot compress bedgraphs"

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "parameters"

# get basename of input file f
INPUTFILENAME=${INPUTFILE##*/}
# get sample prefix
SAMPLE=${INPUTFILENAME/%$ASD.bam/}

# delete old bam files unless attempting to recover
if [ -z "$NGSANE_RECOVERFROM" ]; then
    [ -f $OUTDIR/$SAMPLE.bw ] && rm $OUTDIR/$SAMPLE.bw
fi

if [ -z "$FASTA" ] || [ ! -f $FASTA ]; then
    echo "[ERROR] no reference provided (FASTA)"
    exit 1
else
    echo "[NOTE] Reference: $FASTA"
fi

if [ -z "$UCSC_GENOMEASSEMBLY" ]; then
    echo "[ERROR] Genome assembly not set"
    exit 1
fi

GENOME_CHROMSIZES=${FASTA%.*}.chrom.sizes
if [ ! -f $GENOME_CHROMSIZES ]; then
    echo "[ERROR] GENOME_CHROMSIZES not found. Excepted at $GENOME_CHROMSIZES"
    exit 1
else
    echo "[NOTE] Chromosome size: $GENOME_CHROMSIZES"
fi

# unique temp folder that should be used to store temporary files
THISTMP=$TMP"/"$(whoami)"/"$(echo $OUTDIR/$SAMPLE | md5sum | cut -d' ' -f1)
mkdir -p $THISTMP

GENOME="hg19"

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "recall files from tape"

if [ -n "$DMGET" ]; then
	dmget -a $(dirname $FASTA)/*
	dmget -a $INPUTFILE
    dmget -a $OUTDIR/*
fi
    
NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "hotspot config"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    echo "[NOTE] create tokenizer"
    echo "[script-tokenizer]" > $OUTDIR/$SAMPLE.token
    echo "_TAGS_ = $INPUTFILE" >> $OUTDIR/$SAMPLE.token
    echo "_OUTDIR_ = $OUTDIR/$SAMPLE" >> $OUTDIR/$SAMPLE.token
    echo "_RANDIR_ = $THISTMP" >> $OUTDIR/$SAMPLE.token
    echo "$HOTSPOTCONFIG" | sed -e "s|HOTSPOT_HOME|$HOTSPOT_HOME|g" >> $OUTDIR/$SAMPLE.token

    echo "[NOTE] create config"
    cat <<EOF > $OUTDIR/$SAMPLE.sh
#!/bin/bash
## Do everything, including badspots and final cleanup
scripts="$HOTSPOT_HOME/pipeline-scripts/run_badspot \
    $HOTSPOT_HOME/pipeline-scripts/run_make_lib \
    $HOTSPOT_HOME/pipeline-scripts/run_wavelet_peak_finding \
    $HOTSPOT_HOME/pipeline-scripts/run_10kb_counts \
    $HOTSPOT_HOME/pipeline-scripts/run_generate_random_lib \
    $HOTSPOT_HOME/pipeline-scripts/run_pass1_hotspot \
    $HOTSPOT_HOME/pipeline-scripts/run_pass1_merge_and_thresh_hotspots \
    $HOTSPOT_HOME/pipeline-scripts/run_pass2_hotspot \
    $HOTSPOT_HOME/pipeline-scripts/run_rescore_hotspot_passes \
    $HOTSPOT_HOME/pipeline-scripts/run_spot \
    $HOTSPOT_HOME/pipeline-scripts/run_thresh_hot.R \
    $HOTSPOT_HOME/pipeline-scripts/run_both-passes_merge_and_thresh_hotspots \
    $HOTSPOT_HOME/pipeline-scripts/run_add_peaks_per_hotspot \
    $HOTSPOT_HOME/pipeline-scripts/run_final"

python $HOTSPOT_HOME/ScriptTokenizer/src/script-tokenizer.py \
    --clobber \
    --execute-scripts \
    --output-dir=$OUTDIR \
    $OUTDIR/$SAMPLE.token \
    \$scripts
EOF

    chmod 777 $OUTDIR/$SAMPLE.sh

    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE.sh
fi

################################################################################
NGSANE_CHECKPOINT_INIT "run hotspot"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    cd $THISTMP
    $OUTDIR/$SAMPLE.sh
    cd $SOURCE
    
    # mark checkpoint
    if [ -f $OUTDIR/$SAMPLE.narrowPeak ];then echo -e "\n********* $CHECKPOINT\n"; unset NGSANE_RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
fi

################################################################################
#NGSANE_CHECKPOINT_INIT "generate bigbed"
#
#if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then
#
#    # convert to bigbed
#	if hash bedToBigBed ; then 
#        echo "[NOTE] create bigbed from peaks" 
#        awk '{OFS="\t"; print $1,$2,$3,$7}' $OUTDIR/$SAMPLE.narrowPeak > $OUTDIR/$SAMPLE.peak.tmp
#        bedToBigBed -type=bed4 $OUTDIR/$SAMPLE.peak.tmp $GENOME_CHROMSIZES $OUTDIR/$SAMPLE.bb
#        rm $OUTDIR/$SAMPLE.peak.tmp
#         # mark checkpoint
#        NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE.bb
#    else
#        echo "[NOTE] bigbed not generated"
#        NGSANE_CHECKPOINT_CHECK
#    fi
#fi   
###############################################################################
NGSANE_CHECKPOINT_INIT "cleanup"

[ -d $THISTMP ] && rm -r $THISTMP

NGSANE_CHECKPOINT_CHECK
################################################################################
[ -e $OUTDIR/$SAMPLE.narrowPeak.dummy ] && rm $OUTDIR/$SAMPLE.narrowPeak.dummy
echo ">>>>> Peak calling with hotspot - FINISHED"
echo ">>>>> enddate "`date`
