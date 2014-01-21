#!/bin/bash -e

# DNase-Seq processing using fseq
# author: Fabian Buske
# date: January 1914

# QCVARIABLES,Resource temporarily unavailable
# RESULTFILENAME <DIR>/<TASK>/<SAMPLE>.bw

echo ">>>>> DNase-Seq with fseq"
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
CHECKPOINT="programs"

for MODULE in $MODULE_FSEQ; do module load $MODULE; done  # save way to load modules that itself load other modules
export PATH=$PATH_FSEQ:$PATH
module list
echo "PATH=$PATH"

echo "[NOTE] set java parameters"
JAVAPARAMS="-Xmx"$(python -c "print int($MEMORY_FSEQ*0.8)")"g -Djava.io.tmpdir="$TMP" -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1" 
unset _JAVA_OPTIONS
echo "JAVAPARAMS "$JAVAPARAMS

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--JAVA        --\n" $(java -Xmx200m -version 2>&1)
[ -z "$(which java)" ] && echo "[ERROR] no java detected" && exit 1
echo -e "--bedtools --\n "$(bedtools --version)
[ -z "$(which bedtools)" ] && echo "[ERROR] no bedtools detected" && exit 1
echo -e "--wigToBigWig --\n "$(wigToBigWig 2>&1 | tee | head -n 1)
[ -z "$(which wigToBigWig)" ] && echo "[ERROR] wigToBigWig not detected" && exit 1
echo -e "--fseq        --\n "$(fseq -v | head -n 1)
[ -z "$(which | head -n 1)" ] && echo "[ERROR] no | head -n 1 detected" && exit 1

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="parameters"

# get basename of input file f
INPUTFILENAME=${INPUTFILE##*/}
# get sample prefix
SAMPLE=${INPUTFILENAME/%$bam/}

# delete old bam files unless attempting to recover
if [ -z "$RECOVERFROM" ]; then
    ## TODO remove primary result files from pervious runs
fi

if [ -z "$FASTA" ] || [ ! -f $FASTA ]; then
    echo "[ERROR] no reference provided (FASTA)"
    exit 1
else
    echo "[NOTE] Reference: $FASTA"
fi

GENOME_CHROMSIZES=${FASTA%.*}.chrom.sizes
if [ ! -f $GENOME_CHROMSIZES ]; then
    echo "[ERROR] GENOME_CHROMSIZES not found. Excepted at $GENOME_CHROMSIZES"
    exit 1
else
    echo "[NOTE] Chromosome size: $GENOME_CHROMSIZES"
fi


if [ -n "$FSEQ_PLOIDY" ]; then
    FSEQ_PLOIDY_PARAM="-p $FSEQ_PLOIDY"
fi

# unique temp folder that should be used to store temporary files
THISTMP=$TMP"/"$(whoami)"/"$(echo $OUTDIR | md5sum | cut -d' ' -f1)
mkdir -p $THISTMP

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="recall files from tape"

if [ -n "$DMGET" ]; then
	dmget -a $(dirname $FASTA)/*
	dmget -a $INPUTFILE
    dmget -a $OUTDIR/*
fi
    
echo -e "\n********* $CHECKPOINT\n"

################################################################################
CHECKPOINT="fseq"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else

    rm -f $THISTMP/*wig

#    RUN_COMMAND="bedtools bamtobed -i $f | fseq $FSEQADDPARAM $FSEQ_PLOIDY_PARAM -b $FSEQ_BACKGROUNDDIR -o $THISTMP -of wig "
    RUN_COMMAND="bedtools bamtobed -i $f | java $JAVAOPTS -cp $CLASSPATH edu.duke.igsp.gkde.Main $FSEQADDPARAM $FSEQ_PLOIDY_PARAM -b $FSEQ_BACKGROUNDDIR -o $THISTMP -of wig "

    echo $RUN_COMMAND && eval $RUN_COMMAND
    
    RUN_COMMAND="cat $THISTMP/*.wig | wigToBigWig stdin ${GENOME_CHROMSIZES} $OUTDIR/SAMPLE.bw"
    echo $RUN_COMMAND && eval $RUN_COMMAND

    # mark checkpoint
    if [ -f $OUTDIR/$SAMPLE.bw ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
fi
###############################################################################
CHECKPOINT="cleanup"

[ -d $THISTMP ] && rm -r $THISTMP

echo -e "\n********* $CHECKPOINT\n"
################################################################################
[ -e $OUTDIR/$SAMPLE.bw.dummy ] && rm $OUTDIR/$SAMPLE.bw.dummy
echo ">>>>> DNase-Seq with fseq - FINISHED"
echo ">>>>> enddate "`date`
