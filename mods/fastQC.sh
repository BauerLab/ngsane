#!/bin/bash -e
# author: Denis C. Bauer
# date: Feb.2011

# QCVARIABLES,Resource temporarily unavailable
# RESULTFILENAME <DIR>/<TASK>/<SAMPLE>$READONE"_"fastqc.zip

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

for MODULE in $MODULE_FASTQC; do module load $MODULE; done  # save way to load modules that itself load other modules
export PATH=$PATH_FASTQC:$PATH;
module list
echo "PATH=$PATH"
echo "[NOTE] set java parameters"
JAVAPARAMS="-Xmx"$(python -c "print int($MEMORY_BOWTIE*0.8)")"g -Djava.io.tmpdir="$TMP" -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1" 
unset _JAVA_OPTIONS
echo "JAVAPARAMS "$JAVAPARAMS

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--JAVA        --\n" $(java -Xmx200m -version 2>&1)
[ -z "$(which java)" ] && echo "[ERROR] no java detected" && exit 1
echo -e "--FASTqc      --\n" $(fastqc -version 2>&1)
[ -z "$(which fastqc)" ] && echo "[ERROR] no fastqc detected" && exit 1

echo -e "\n********* $CHECKPOINT\n"

################################################################################
CHECKPOINT="parameters"

# get basename of input file f
INPUTFILENAME=${INPUTFILE##*/}
# get sample prefix
SAMPLE=${INPUTFILENAME/%$READONE.$INPUTFILEASTQ/}

#is paired ?                                                                                                      
if [ "$INPUTFILE" != "${INPUTFILE/%$READONE.$FASTQ/$READTWO.$FASTQ}" ] && [ -e ${INPUTFILE/%$READONE.$FASTQ/$READTWO.$FASTQ} ]; then
    PAIRED="1"
else
    PAIRED="0"
fi

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="recall files from tape"

if [ -n "$DMGET" ]; then
	dmget -a $INPUTFILE
    dmget -a $OUTDIR/*
fi
    
echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="fastqc read 1"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

    RUN_COMMAND="fastqc --nogroup -t $CPU_FASTQC --outdir $OUTDIR $INPUTFILE"
    echo $RUN_COMMAND && eval $RUN_COMMAND
    # check for ".fastq.gz" suffix as FASTQC removes both suffixes then
    if [ "$FASTQ" != "fastq.gz" ];then 
        mv $OUTDIR/${INPUTFILENAME%.*}"_"fastqc.zip $OUTDIR/${INPUTFILENAME/%.$FASTQ/"_"fastqc.zip}; 
        mv $OUTDIR/${INPUTFILENAME%.*}"_"fastqc $OUTDIR/${INPUTFILENAME/%.$FASTQ/"_"fastqc}
    fi

    chmod -R a+rx $OUTDIR/
    
    if [ -f $OUTDIR/${INPUTFILENAME/%.$FASTQ/"_"fastqc.zip} ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
fi
################################################################################
CHECKPOINT="fastqc read 2"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

    if [ "$PAIRED" = "1" ];then
        RUN_COMMAND="fastqc --nogroup -t $CPU_FASTQC --outdir $OUTDIR ${INPUTFILE/%$READONE.$FASTQ/$READTWO.$FASTQ}"
        echo $RUN_COMMAND && eval $RUN_COMMAND
        R2=${INPUTFILENAME/%$READONE.$FASTQ/$READTWO.$FASTQ}
        if [ "$FASTQ" != "fastq.gz" ];then 
            mv $OUTDIR/${R2%.*}"_"fastqc.zip $OUTDIR/${INPUTFILENAME/%$READONE.$FASTQ/$READTWO"_"fastqc.zip}; 
            mv $OUTDIR/${R2%.*}"_"fastqc $OUTDIR/${INPUTFILENAME/%$READONE.$FASTQ/$READTWO"_"fastqc}; 
        fi
    else
        echo "[NOTE] Single-end pair library detected. No second read to process."
    fi
    
    chmod -R a+rx $OUTDIR/
    
    if [ -f $OUTDIR/${INPUTFILENAME/%.$FASTQ/"_"fastqc.zip} ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
fi
################################################################################
[ -e $OUTDIR/${INPUTFILENAME/%.$FASTQ/"_"fastqc.zip}.dummy ] && rm $OUTDIR/${INPUTFILENAME/%.$FASTQ/"_"fastqc.zip}.dummy
echo ">>>>> fastQC - FINISHED"
echo ">>>>> enddate "`date`

