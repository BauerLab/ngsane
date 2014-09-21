#!/bin/bash -e

# author: Denis Bauer
# date: April 2014

# QCVARIABLES,Resource temporarily unavailable
# RESULTFILENAME <DIR>/<TASK>/<SAMPLE>$READONE.$FASTQ

echo ">>>>> SRA converter"
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
NGSANE_CHECKPOINT_INIT "programs"

for MODULE in $MODULE_SRACONV; do module load $MODULE; done  # save way to load modules that itself load other modules
export PATH=$PATH_SRACONV:$PATH
module list
echo "PATH=$PATH"

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
## [TODO] test and output versions of software utilized in this mod 

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "parameters"

# get basename of input file f
INPUTFILENAME=${INPUTFILE##*/}
# get sample prefix
SAMPLE=${INPUTFILENAME/%.sra/}

# delete old bam files unless attempting to recover
if [ -z "$RECOVERFROM" ]; then
    [ -e $OUTDIR/$SAMPLE$READONE.$FASTQ ] && rm $OUTDIR/$SAMPLE$READONE.$FASTQ
fi

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "recall files from tape"

if [ -n "$DMGET" ]; then
	dmget -a $INPUTFILE
fi
    
NGSANE_CHECKPOINT_CHECK

################################################################################
NGSANE_CHECKPOINT_INIT "SRA conversion"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    command="fastq-dump $FASTQDUMADDPARAM $INPUTFILE -O $OUTDIR/ --gzip $SPLIT"
    echo $command && eval $command

    [ -e $OUTDIR/$SAMPLE"_1".$FASTQ ] && mv $OUTDIR/$SAMPLE"_1".$FASTQ $OUTDIR/$SAMPLE$READONE.$FASTQ
    [ -e $OUTDIR/$SAMPLE"_2".$FASTQ ] && mv $OUTDIR/$SAMPLE"_2".$FASTQ $OUTDIR/$SAMPLE$READTWO.$FASTQ
 
    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE$READONE.$FASTQ
fi

################################################################################
## TODO: specify primary output file suffix (same as at the top RESULTFILENAME section)
[ -e $OUTDIR/$SAMPLE$READONE.$FASTQ.dummy ] && rm $OUTDIR/$SAMPLE$READONE.$FASTQ.dummy
echo ">>>>> SRA converter - FINISHED"
echo ">>>>> enddate "`date`
