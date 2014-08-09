#!/bin/bash -e

# Template for new mods
# Work through the TODOs and replace all [TEMPLATE...] patterns
# author: Fabian Buske
# date: November 2013, revised August 2014

# QCVARIABLES,Resource temporarily unavailable
## [TODO] specify primary output file suffix 
# RESULTFILENAME <DIR>/<TASK>/<SAMPLE>.[TEMPLATE_RESULT_FILE_SUFFIX]

## [TODO] change [TEMPLATE] to program
echo ">>>>> [TEMPLATE purpose]"
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

## [TODO] change [TEMPLATE] to program
# save way to load modules that itself load other modules
hash module 2>/dev/null && for MODULE in $MODULE_[TEMPLATE]; do module load $MODULE; done && module list

export PATH=$PATH_[TEMPLATE]:$PATH
echo "PATH=$PATH"

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
## [TODO] test and output versions of software utilized in this mod 

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "parameters"

# get basename of input file f
INPUTFILENAME=${INPUTFILE##*/}
# get sample prefix
SAMPLE=${INPUTFILENAME/%$READONE.$FASTQ/}

# delete old bam files unless attempting to recover
if [ -z "$NGSANE_RECOVERFROM" ]; then
    ## TODO remove primary result files from pervious runs
fi

## TODO remove comments if paired/single library preps should be detected based on the READ identifier patterns
#if [ "$INPUTFILE" != "${INPUTFILE/$READONE/$READTWO}" ] && [ -e ${INPUTFILE/$READONE/$READTWO} ]; then
#    PAIRED="1"
#else
#    PAIRED="0"
#fi

## TODO remove comments if compressed status of input files should be detected
#CAT="cat"
#if [[ ${f##*.} == "gz" ]]; 
#    then CAT="zcat"; 
#elif [[ ${f##*.} == "bz2" ]]; 
#    then CAT="bzcat"; 
#fi

# unique temp folder that should be used to store temporary files
THISTMP=$TMP"/"$(whoami)"/"$(echo $OUTDIR | md5sum | cut -d' ' -f1)
mkdir -p $THISTMP

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "recall files from tape"

if [ -n "$DMGET" ]; then
	dmget -a $INPUTFILE
    dmget -a $OUTDIR/*
    # TODO add additional resources that are required and may need recovery from tape
fi
    
NGSANE_CHECKPOINT_CHECK

## TODO choose a checkpoint section (with or without result file check) from the 
##      two blocks below, add as many such sections as needed.
################################################################################
# TODO add template checkpoint name
NGSANE_CHECKPOINT_INIT "[TEMPLATE CHECKPOINT WITH RESULT FILE CHECK]"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    ## TODO: add program calls that constitute a self-contained block
    [TEMPLATE PROGRAM CALL 1] [TEMPLATE ADDITIONAL PROGRAM PARAMETERS] [TEMPLATE INPUTFILENAME] [TEMPLATE INTERMEDIATE FILE 1]
    [TEMPLATE PROGRAM CALL 2] [TEMPLATE ADDITIONAL PROGRAM PARAMETERS] [TEMPLATE INTERMEDIATE FILE 1] [TEMPLATE INTERMEDIATE FILE 2]
    

    ## TODO: specify intermediate result file pattern for testing 
    # mark checkpoint checking for existance/non-empty result file
    [[ -s $OUTDIR/$SAMPLE.[TEMPLATE_INTERMEDIATE_FILE_SUFFIX] ]] && NGSANE_CHECKPOINT_CHECK 
fi

################################################################################
# TODO add template checkpoint name
NGSANE_CHECKPOINT_INIT "[TEMPLATE CHECKPOINT WITHOUT RESULT FILE CHECK]"
    
## TODO: add program calls that constitute a self-contained block 
[TEMPLATE PROGRAM CALL 1] [TEMPLATE ADDITIONAL PROGRAM PARAMETERS] [TEMPLATE INPUTFILENAME] [TEMPLATE INTERMEDIATE FILE 1]
[TEMPLATE PROGRAM CALL 2] [TEMPLATE ADDITIONAL PROGRAM PARAMETERS] [TEMPLATE INTERMEDIATE FILE 1] [TEMPLATE INTERMEDIATE FILE 2]

# mark checkpoint without checking for result
NGSANE_CHECKPOINT_CHECK
################################################################################
## TODO: specify primary output file suffix (same as at the top RESULTFILENAME section)
[ -e $OUTDIR/$SAMPLE.[TEMPLATE_RESULT_FILE_SUFFIX].dummy ] && rm $OUTDIR/$SAMPLE.[TEMPLATE_RESULT_FILE_SUFFIX].dummy
echo ">>>>> [TEMPLATE purpose] - FINISHED"
echo ">>>>> enddate "`date`
