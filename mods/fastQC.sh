#!/bin/bash -e
# author: Denis C. Bauer
# date: Feb.2011

# QCVARIABLES,Resource temporarily unavailable
# RESULTFILENAME <DIR>/<TASK>/<SAMPLE>$READONE.html

echo ">>>>> fastQC"
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"

while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the NGSANE 
        -f | --file )           shift; INPUTFILE=$1 ;;  # input file                                                       
        -o | --outdir )         shift; OUTDIR=$1 ;;     # output dir                                                     
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
hash module 2>/dev/null && for MODULE in $MODULE_FASTQC; do module load $MODULE; done && module list 

export PATH=$PATH_FASTQC:$PATH;
echo "PATH=$PATH"
echo "[NOTE] set java parameters"
JAVAPARAMS="-Xmx"$(python -c "print int($MEMORY_FASTQC*0.8)")"g -Djava.io.tmpdir="$TMP" -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1" 
unset _JAVA_OPTIONS
echo "JAVAPARAMS "$JAVAPARAMS
export _JAVA_OPTIONS=$JAVAPARAMS

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--JAVA        --\n" $(java -Xmx200m -version 2>&1)
[ -z "$(which java)" ] && echo "[ERROR] no java detected" && exit 1
echo -e "--FASTqc      --\n" $(fastqc -version 2>&1)
[ -z "$(which fastqc)" ] && echo "[ERROR] no fastqc detected" && exit 1

NGSANE_CHECKPOINT_CHECK

################################################################################
NGSANE_CHECKPOINT_INIT "parameters"

# get basename of input file f
INPUTFILENAME=${INPUTFILE##*/}
# get sample prefix
SAMPLE=${INPUTFILENAME/%$READONE.$FASTQ/}

#is paired ?                                                                                                      
if [ "$INPUTFILE" != "${INPUTFILE/%$READONE.$FASTQ/$READTWO.$FASTQ}" ] && [ -e ${INPUTFILE/%$READONE.$FASTQ/$READTWO.$FASTQ} ]; then
    PAIRED="1"
else
    PAIRED="0"
fi


    [ -d $OUTDIR/$SAMPLE$READONE ] && rm -r $OUTDIR/$SAMPLE$READONE
    [ -e $OUTDIR/$SAMPLE$READONE.html ] && rm $OUTDIR/$SAMPLE$READONE.html
    [ -d $OUTDIR/$SAMPLE$READTWO ] && rm -r $OUTDIR/$SAMPLE$READTWO
    [ -e $OUTDIR/$SAMPLE$READTWO.html ] && rm $OUTDIR/$SAMPLE$READTWO.html


NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "recall files from tape"

if [ -n "$DMGET" ]; then
	dmget -a $INPUTFILE
    dmget -a $OUTDIR/*
fi
    
NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "fastqc read 1"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    RUN_COMMAND="fastqc $FASTQCADDPARAM --threads $CPU_FASTQC --outdir $OUTDIR $INPUTFILE"
    echo $RUN_COMMAND && eval $RUN_COMMAND
    # check for ".fastq.gz" suffix as FASTQC removes both suffixes then
    if [ "$FASTQ" != "fastq.gz" ];then 

        unzip -q -d $OUTDIR -o -u $OUTDIR/${INPUTFILENAME%.*}"_"fastqc.zip
        mv $OUTDIR/${INPUTFILENAME%.*}"_"fastqc $OUTDIR/$SAMPLE$READONE
        mv $OUTDIR/${INPUTFILENAME%.*}"_"fastqc.html $OUTDIR/$SAMPLE$READONE.html
    else
        unzip -q -d $OUTDIR -o -u $OUTDIR/$SAMPLE$READONE"_"fastqc.zip
        mv $OUTDIR/$SAMPLE$READONE"_"fastqc $OUTDIR/$SAMPLE$READONE
        mv $OUTDIR/$SAMPLE$READONE"_"fastqc.html $OUTDIR/$SAMPLE$READONE.html
    fi
    cp $OUTDIR/$SAMPLE$READONE/fastqc_data.txt $OUTDIR/$SAMPLE$READONE.txt

    chmod -R a+rx $OUTDIR/
    
    NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE$READONE.txt
fi
################################################################################
NGSANE_CHECKPOINT_INIT "fastqc read 2"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    if [ "$PAIRED" = "1" ];then
        RUN_COMMAND="fastqc $FASTQCADDPARAM --threads $CPU_FASTQC --outdir $OUTDIR ${INPUTFILE/%$READONE.$FASTQ/$READTWO.$FASTQ}"
        echo $RUN_COMMAND && eval $RUN_COMMAND
        R2=${INPUTFILENAME/%$READONE.$FASTQ/$READTWO.$FASTQ}
        if [ "$FASTQ" != "fastq.gz" ];then 

            unzip -q -d $OUTDIR -o -u $OUTDIR/${R2%.*}"_"fastqc.zip
            mv $OUTDIR/${R2%.*}"_"fastqc $OUTDIR/$SAMPLE$READTWO
            mv $OUTDIR/${R2%.*}"_"fastqc.html $OUTDIR/$SAMPLE$READTWO.html
        else
            unzip -q -d $OUTDIR -o -u $OUTDIR/$SAMPLE$READTWO"_"fastqc.zip
            mv $OUTDIR/$SAMPLE$READTWO"_"fastqc $OUTDIR/$SAMPLE$READTWO
            mv $OUTDIR/$SAMPLE$READTWO"_"fastqc.html $OUTDIR/$SAMPLE$READTWO.html
        fi
        cp $OUTDIR/$SAMPLE$READTWO/fastqc_data.txt $OUTDIR/$SAMPLE$READTWO.txt

        chmod -R a+rx $OUTDIR/
        
        NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE$READTWO.txt

    else
        echo "[NOTE] Single-end pair library detected. No second read to process."
        NGSANE_CHECKPOINT_CHECK
    fi
    
fi
################################################################################
[ -e $OUTDIR/$SAMPLE$READONE.html.dummy ] && rm $OUTDIR/$SAMPLE$READONE.html.dummy
echo ">>>>> fastQC - FINISHED"
echo ">>>>> enddate "`date`

