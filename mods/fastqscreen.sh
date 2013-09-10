#!/bin/bash -e

# Script to ... 
# It takes a <Run>/*.$FASTQ[.gz] file and gets the file containing the contaminats
#
# author: 
# date: 

# messages to look out for -- relevant for the QC.sh script:
# QCVARIABLES,

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
        -o | --outdir )         shift; MYOUT=$1 ;; # output dir
        --recover-from )        shift; RECOVERFROM=$1 ;; # attempt to recover from log file
        -h | --help )           usage ;;
        * )                     echo "don't understand "$1
    esac
    shift
done

. $CONFIG
. ${NGSANE_BASE}/conf/header.sh
. $CONFIG

################################################################################
CHECKPOINT="programs"

for MODULE in $MODULE_FASTQSCREEN; do module load $MODULE; done  # save way to load modules that itself load other modules
export PATH=$PATH_FASTQSCREEN:$PATH;
module list
echo "PATH=$PATH"
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--JAVA        --\n" $(java -version 2>&1)
[ -z "$(which java)" ] && echo "[ERROR] no java detected" && exit 1
echo -e "--fastq_screen --\n "  $(`which perl` `which fastq_screen` --version)
[ ! -f $(which fastq_screen) ] && echo "[ERROR] no fastq_screen detected" && exit 1

echo "[NOTE] set java parameters"
JAVAPARAMS="-Xmx"$(expr $MEMORY_FASTQSCREEN - 1 )"G -Djava.io.tmpdir="$TMP
echo "JAVAPARAMS "$JAVAPARAMS

echo -e "\n[CHECKPOINT] $CHECKPOINT\n"
################################################################################
CHECKPOINT="parameters"

# get basename of f
n=${f##*/}

#is paired ?
if [ "$f" != "${f/$READONE/$READTWO}" ] && [ -e ${f/$READONE/$READTWO} ]; then
    echo "[NOTE] PAIRED library"
    PAIRED="1"
else
    echo "[NOTE] SINGLE library"
    PAIRED="0"
fi

if [ -z "$FASTQSCREEN_DBCONF" ] || [ ! -f $FASTQSCREEN_DBCONF ]; then
    echo "[NOTE] FASTQSCREEN_DBCONF file not specified or found" && exit 1
fi
 
mkdir -p $MYOUT

echo -e "\n[CHECKPOINT] $CHECKPOINT\n"
################################################################################
CHECKPOINT="recall files from tape"

if [ -n "$DMGET" ]; then
    dmget -a ${f/$READONE/"*"}
fi

echo -e "\n[CHECKPOINT] $CHECKPOINT\n"
################################################################################
CHECKPOINT="fastq screening"    

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\[CHECKPOINT\] $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 


    # Paired read
    if [ "$PAIRED" = "1" ]; then
        RUN_COMMAND="`which perl` `which fastq_screen` $FASTQSCREENADDPARAM --outdir $MYOUT --conf $FASTQSCREEN_DBCONF --paired --threads $CPU_FASTQSCREEN $f ${f/$READONE/$READTWO}"
    else
        RUN_COMMAND="`which perl` `which fastq_screen` $FASTQSCREENADDPARAM --outdir $MYOUT --conf $FASTQSCREEN_DBCONF --threads $CPU_FASTQSCREEN $f"
    fi
    echo $RUN_COMMAND && eval $RUN_COMMAND

    mv $MYOUT/${n}_screen.txt $MYOUT/${n/$READONE.$FASTQ/}_screen.txt
    mv $MYOUT/${n}_screen.png $MYOUT/${n/$READONE.$FASTQ/}_screen.png

    # mark checkpoint
    if [ -f $MYOUT/${n/$READONE.$FASTQ/}_screen.txt ];then echo -e "\n[CHECKPOINT] $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi

################################################################################
echo ">>>>> read screening with  FASTQSCREEN - FINISHED"
echo ">>>>> enddate "`date`

