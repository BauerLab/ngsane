#!/bin/bash -e
# Martin Smith, August 2013
# tested with trinityrnaseq_r2013-02-25

# QCVARIABLES,Resource temporarily unavailable
# RESULTFILENAME <DIR>/<TASK>/<SAMPLE>.$ASD.bam

echo ">>>>> transcriptome assembly with trinity chrysalis"
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"

function usage {
echo -e "usage: $(basename $0) -k NGSANE -f FASTQ -r REFERENCE -o OUTDIR [OPTIONS]"
exit
}

if [ ! $# -gt 3 ]; then usage ; fi

while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the configuration file
        -f | --fastq )          shift; f=$1 ;; # fastq file
        -o | --outdir )         shift; OUTDIR=$1 ;; # output dir
        --recover-from )        shift; RECOVERFROM=$1 ;; # attempt to recover from log file
        -h | --help )           usage ;;
        * )                     echo "don't understand "$1
    esac
    shift
done

. $CONFIG
. ${NGSANE_BASE}/conf/header.sh
. $CONFIG

if [ -e $NODE_TMP ]; then TMP=$NODE_TMP; fi     # overwrites temp folder variable to node-local temp folder, if specified

################################################################################
CHECKPOINT="programs"
for MODS in $MODULES_TRINITY; do module load $MODS; done 
export PATH=$PATH_TRINITY:$PATH
module list
echo "PATH=$PATH"

echo "[NOTE] set java parameters"
JAVAPARAMS="-Xmx"$(python -c "print int($MEMORY_CHRYSALIS*0.8)")"g -Djava.io.tmpdir="$TMP" -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1" 
unset _JAVA_OPTIONS
echo "JAVAPARAMS "$JAVAPARAMS

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--JAVA        --\n" $(java -Xmx200m -version 2>&1)
[ -z "$(which java)" ] && echo "[ERROR] no java detected" && exit 1
echo -e "--bowtie      --\n "$(bowtie --version)
[ -z "$(which bowtie)" ] && echo "[ERROR] no bowtie detected" && exit 1
echo -e "--perl        --\n "$(perl -v | grep "version" )
[ -z "$(which perl)" ] && echo "[ERROR] no perl detected" && exit 1
echo -e "--trinity     --\n "$(trinity --version)
[ -z "$(which trinity)" ] && echo "[ERROR] no trinity detected" && exit 1

ulimit -s unlimited

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="recall files from tape"

if [ -n "$DMGET" ]; then
	dmget -a $(dirname $FASTA)/*
	dmget -a ${f/$READONE/"*"}
fi
    
echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="detect library"

# get basename of f
n=${f##*/}

## is paired ?                                                                                                      
if [ -e ${f/$READONE/$READTWO} ]; then
    PAIRED="1"
    f2=${f/$READONE/$READTWO}
    echo "[NOTE] Paired library detected: "$f" + "$f2
else
    PAIRED="0"
    echo "[ERROR] Single-Strand (unpaired) library detected !!!! currently unsupported !!!!"
    exit 1
fi

#TODO check tmp dir
mkdir $TMP/$JOB_ID
TMP_LOC=${f#*$SOURCE/fastq/}  # project name, i.e. sample1
cp $SOURCE/${TMP_LOC%/*}/trinity/ $TMP/$JOB_ID
cd $TMP/$JOB_ID

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="Chrysalis"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

    echo "[NOTE] --max_reads_per_graph set to 1 million because very high I/O is needed otherwise. It is unlikely that a transcript needs more than 1 million reads to be assembled"

    if [ -z "$SS_LIBRARY_TYPE" ]; then
        echo "[WARN] strand-specific RNAseq library type not set ($SS_LIBRARY_TYPE): treating input as non-stranded"
        RUN_COMMAND="$(which perl) Trinity.pl --seqType fq --left $f --right $f2 --output /tmp/$JOB_ID --JM $MEMORY_CHRYSALIS"G" --CPU $NCPU_CHRYSALIS --no_run_quantifygraph"
    else
        echo "[NOTE] RNAseq library type: $SS_LIBRARY_TYPE"
        RUN_COMMAND="$(which perl) Trinity.pl --seqType fq --left $f --right $f2 --SS_lib_type $SS_LIBRARY_TYPE --output /tmp/$JOB_ID --JM $MEMORY_CHRYSALIS"G" --CPU $NCPU_CHRYSALIS --no_run_quantifygraph"
    fi
    
    mv $TMP/$JOB_ID/* $SOURCE/${TMP_LOC%/*}/trinity/
    
    echo $RUN_COMMAND && eval $RUN_COMMAND

    echo "[NOTE] Chrysalis pt.1 has completed properly! thank Martin by buying him another beer"

    # mark checkpoint
    #TODO: result file?
    if [ -f  ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
fi 

################################################################################
#TODO correct dummy
[ -e $OUTDIR/${n/%$READONE.$FASTQ/.$ASD.bam}.dummy ] && rm $OUTDIR/${n/%$READONE.$FASTQ/.$ASD.bam}.dummy
echo ">>>>> transcriptome assembly with trinity chrysalis - FINISHED"
echo ">>>>> enddate "`date`