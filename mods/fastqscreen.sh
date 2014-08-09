#!/bin/bash -e

# Script to screen a fastq library for contamination
# It takes a <Run>/*.$FASTQ[.gz] file and gets the file containing the contaminats
#
# author: 
# date: 

# messages to look out for -- relevant for the QC.sh script:
# QCVARIABLES,Resource temporarily unavailable
# RESULTFILENAME <DIR>/<TASK>/<SAMPLE>_screen.txt

echo ">>>>> read screening with FASTQSCREEN"
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"

while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of NGSANE
        -f | --file )           shift; f=$1 ;; # fastq file
        -o | --outdir )         shift; OUTDIR=$1 ;; # output dir
        --recover-from )        shift; NGSANE_RECOVERFROM=$1 ;; # attempt to recover from log file
        -h | --help )           usage ;;
        * )                     echo "don't understand "$1
    esac
    shift
done

. $CONFIG
. ${NGSANE_BASE}/conf/header.sh
. $CONFIG

################################################################################
NGSANE_CHECKPOINT_INIT "programs"

# save way to load modules that itself loads other modules
hash module 2>/dev/null && for MODULE in $MODULE_FASTQSCREEN; do module load $MODULE; done && module list 

export PATH=$PATH_FASTQSCREEN:$PATH;
echo "PATH=$PATH"
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)

echo -e "--NGSANE       --\n" $(trigger.sh -v 2>&1)
echo -e "--perl         --\n "$(perl -v | grep "This is perl" )
[ -z "$(which perl)" ] && echo "[ERROR] no perl detected" && exit 1
echo -e "--fastq_screen --\n "  $(`which perl` `which fastq_screen` --version)
[ ! -f $(which fastq_screen) ] && echo "[ERROR] no fastq_screen detected" && exit 1

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "parameters"

# get basename of f
n=${f##*/}
SAMPLE=${n/%$READONE.$FASTQ/}

#is ziped ?
CAT="cat"
if [[ ${f##*.} == "gz" ]]; 
    then CAT="zcat"; 
elif [[ ${f##*.} == "bz2" ]]; 
    then CAT="bzcat"; 
fi

#is paired ?
if [ "$f" != "${f/%$READONE.$FASTQ/$READTWO.$FASTQ}" ] && [ -e ${f/%$READONE.$FASTQ/$READTWO.$FASTQ} ]; then
    echo "[NOTE] PAIRED library"
    PAIRED="1"
else
    echo "[NOTE] SINGLE library"
    PAIRED="0"
fi

if [ -z "$FASTQSCREEN_DBCONF" ] || [ ! -f $FASTQSCREEN_DBCONF ]; then
    echo "[NOTE] FASTQSCREEN_DBCONF file not specified or found" && exit 1
fi
 
mkdir -p $OUTDIR

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "recall files from tape"

if [ -n "$DMGET" ]; then
    dmget -a ${f/$READONE/"*"}
    dmget -a ${OUTDIR}/*
fi

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "fastq screening"    

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then


    # Paired read
    if [ "$PAIRED" = "1" ]; then
        RUN_COMMAND="`which perl` `which fastq_screen` $FASTQSCREENADDPARAM --outdir $OUTDIR --conf $FASTQSCREEN_DBCONF --paired --threads $CPU_FASTQSCREEN <($CAT $f) <($CAT ${f/%$READONE.$FASTQ/$READTWO.$FASTQ})"
    else
        RUN_COMMAND="`which perl` `which fastq_screen` $FASTQSCREENADDPARAM --outdir $OUTDIR --conf $FASTQSCREEN_DBCONF --threads $CPU_FASTQSCREEN <($CAT $f)"
    fi
    echo $RUN_COMMAND && eval $RUN_COMMAND

    mv $OUTDIR/${n}_screen.txt $OUTDIR/$SAMPLE"_"screen.txt
    mv $OUTDIR/${n}_screen.png $OUTDIR/$SAMPLE"_"screen.png

    # mark checkpoint
    [[ -s $OUTDIR/$SAMPLE"_"screen.txt ]] && NGSANE_CHECKPOINT_CHECK

fi

################################################################################
[ -e $OUTDIR/$SAMPLE"_"screen.txt.dummy ] && rm $OUTDIR/$SAMPLE"_"screen.txt.dummy
echo ">>>>> read screening with  FASTQSCREEN - FINISHED"
echo ">>>>> enddate "`date`

