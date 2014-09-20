#!/bin/bash -e 

# author: Fabian Buske
# date: November 2013

echo ">>>>> Create wig files with wiggler"
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"

function usage {
echo -e "usage: $(basename $0) -k NGSANE -f bam -o OUTDIR [OPTIONS]

Script running wiggler on bam files

required:
  -k | --toolkit <path>     location of the NGSANE repository 
  -f | --bam <file>         bam file
  -o | --outdir <path>      output dir
"
exit
}
# QCVARIABLES,Resource temporarily unavailable
# RESULTFILENAME <DIR>/<TASK>/<SAMPLE>.$WIGGLER_OUTPUTFORMAT.gz

if [ ! $# -gt 3 ]; then usage ; fi

#INPUTS                                                                                                           
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the NGSANE repository                       
        -f | --bam )            shift; f=$1 ;; # bam file                                                       
        -o | --outdir )         shift; OUTDIR=$1 ;; # output dir                                                     
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
hash module 2>/dev/null && for MODULE in $MODULE_WIGGLER; do module load $MODULE; done && module list

export PATH=$PATH_WIGGLER:$PATH
echo "PATH=$PATH"
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--wiggler     --\n "$(align2rawsignal 2>&1 | head -n 3 | tail -n 1)
[ -z "$(which align2rawsignal)" ] && echo "[ERROR] wiggler not detected (align2rawsignal)" && exit 1

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "parameters"

# get basename of f
n=${f##*/}
SAMPLE=${n/%$ASD.bam/}

# check UMAP folder exist
if [ -z "$WIGGLER_UMAPDIR" ] || [ ! -d ${WIGGLER_UMAPDIR} ]; then
    echo "[ERROR] umap dir not specified not non-existant (WIGGLER_UMAPDIR)"
    exit 1
fi

# check reference folder exist
if [ -z "$FASTA_CHROMDIR" ] || [ ! -d ${FASTA_CHROMDIR} ]; then
    echo "[ERROR] Chromosome directory not specified (FASTA_CHROMDIR)"
    exit 1
fi

# check output format
if [ -z "$WIGGLER_OUTPUTFORMAT" ]; then
    echo "[ERROR] wiggler output format not set" && exit 1
elif [ "$WIGGLER_OUTPUTFORMAT" != "bg" ] && [ "$WIGGLER_OUTPUTFORMAT" != "wig" ]; then
    echo "[ERROR] wiggler output format not known" && exit 1
fi

THISTMP=$TMP"/"$(whoami)"/"$(echo $OUTDIR/$SAMPLE | md5sum | cut -d' ' -f1)
mkdir -p $THISTMP

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "recall files from tape"

if [ -n "$DMGET" ]; then
	dmget -a ${f}
	dmget -a $OUTDIR/*
fi

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "calculate fragment size"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then
    
    if [ -z "$WIGGLER_FRAGMENTSIZE" ]; then
        macs2 predictd --rfile $OUTDIR/$SAMPLE -i $f > $OUTDIR/$SAMPLE.fragment_size.txt 2>&1
    fi
    
    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE.fragment_size.txt 
    
fi 
################################################################################
NGSANE_CHECKPOINT_INIT "run align2rawsignal"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then
    
    if [ -z "$WIGGLER_FRAGMENTSIZE" ]; then
        WIGGLER_FRAGMENTSIZE=$(grep 'alternative fragment length(s) may be' $OUTDIR/$SAMPLE.fragment_size.txt | sed 's/.* be //' | cut -d' ' -f 1 | tr ',' '\n' | egrep -v "^-" | head -n 1)
    fi
    
    [ -f $OUTDIR/$SAMPLE.log ] && rm $OUTDIR/$SAMPLE.log
    
    RUN_COMMAND="align2rawsignal $WIGGLERADDPARAMS -f=$WIGGLER_FRAGMENTSIZE -of=$WIGGLER_OUTPUTFORMAT -i=$f -s=${FASTA_CHROMDIR} -u=${WIGGLER_UMAPDIR} -v=$OUTDIR/$SAMPLE.log -o=$THISTMP/$SAMPLE.$WIGGLER_OUTPUTFORMAT -mm=$MEMORY_WIGGLER"
    echo $RUN_COMMAND && eval $RUN_COMMAND

    if hash bedToBigBed 2>&- && [[ "$WIGGLER_OUTPUTFORMAT" == "bg" ]] && [[ -f $GENOME_CHROMSIZES ]]; then
        bedGraphToBigWig $THISTMP/$SAMPLE.$WIGGLER_OUTPUTFORMAT $GENOME_CHROMSIZES $OUTDIR/$SAMPLE.bw
        CHECKPOINTFILE=$OUTDIR/$SAMPLE.bw
    elif hash wigToBigWig 2>&- && [[ "$WIGGLER_OUTPUTFORMAT" == "wig" ]] && [[ -f $GENOME_CHROMSIZES ]]; then 
        wigToBigWig $THISTMP/$SAMPLE.$WIGGLER_OUTPUTFORMAT $GENOME_CHROMSIZES $OUTDIR/$SAMPLE.bw
        CHECKPOINTFILE=$OUTDIR/$SAMPLE.bw
    else
        $GZIP -c $THISTMP/$SAMPLE.$WIGGLER_OUTPUTFORMAT > $OUTDIR/$SAMPLE.$WIGGLER_OUTPUTFORMAT.gz
        CHECKPOINTFILE=$OUTDIR/$SAMPLE.$WIGGLER_OUTPUTFORMAT.gz
    fi
        
    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $CHECKPOINTFILE
    
fi 
################################################################################
NGSANE_CHECKPOINT_INIT "cleanup"

if [ -d $THISTMP ]; then rm -r $THISTMP; fi

NGSANE_CHECKPOINT_CHECK
################################################################################
[ -e $OUTDIR/$SAMPLE.$WIGGLER_OUTPUTFORMAT.gz.dummy ] && rm $OUTDIR/$SAMPLE.$WIGGLER_OUTPUTFORMAT.gz.dummy
echo ">>>>> wiggler - FINISHED"
echo ">>>>> enddate "`date`
