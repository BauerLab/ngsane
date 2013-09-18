#!/bin/bash -e
# Martin Smith, August 2013
# tested with trinityrnaseq_r2013-02-25

# QCVARIABLES,Resource temporarily unavailable
# RESULTFILENAME <DIR>/<TASK>/<SAMPLE>

echo ">>>>> transcriptome assembly with trinity inchworm"
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

# $GE_TASK_ID

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

################################################################################
CHECKPOINT="programs"
for MODS in $MODULES_TRINITY; do module load $MODS; done 
export PATH=$PATH_TRINITY:$PATH
module list
echo "PATH=$PATH"

echo "[NOTE] set java parameters"
JAVAPARAMS="-Xmx"$(python -c "print int($MEMORY_INCHWORM*0.8)")"g -Djava.io.tmpdir="$TMP" -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1" 
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


#if [ $NODETYPE_INCHWORM -eq "intel.q" ] ; then
	#some fancy optimisation stuff for intel nodes. Consider the intel compiler! 
#	echo "[NOTE] forcing copact thread placement on intel nodes"
#	export OMP_NUM_THREADS=$NCPU_INCHWORM
#	export KMP_AFFINITY=compact
#fi

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

# cd /tmp ##??????

# this runs Inchworm only

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="Inchworm"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

    echo "[NOTE] --max_reads_per_graph set to 1 million because very high I/O is needed otherwise. It is unlikely that a transcript needs more than 1 million reads to be assembled"

    if [ -z "$SS_LIBRARY_TYPE" ]; then
        echo "[WARNING] strand-specific RNAseq library type not set ($SS_LIBRARY_TYPE): treating input as non-stranded"
        RUN_COMMAND="Trinity.pl --seqType fq --left $f --right $f2 --max_reads_per_graph 1000000 --output $OUTDIR --JM $MEMORY_INCHWORM"G" --CPU $NCPU_INCHWORM --no_run_chrysalis"
    else
        echo "[NOTE] RNAseq library type: $SS_LIBRARY_TYPE"
        RUN_COMMAND="Trinity.pl --seqType fq --left $f --right $f2 --SS_lib_type $SS_LIBRARY_TYPE --max_reads_per_graph 1000000 --output $OUTDIR --JM $MEMORY_INCHWORM"G" --CPU $NCPU_INCHWORM --no_run_chrysalis"
    fi
    echo $RUN_COMMAND && eval $RUN_COMMAND
    echo "[NOTE] inchworm has completed properly! thank Martin by buying him a beer"

    # mark checkpoint
    #TODO: result file?
    if [ -f  ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
fi 

################################################################################
CHECKPOINT="cleanup"

if [ -n $CLEANUP ]; then 
#    rm $OUTDIR/*ebwt $OUTDIR/*ebwt  #bowtie files required for chrysalis, and this is the WRONG file location (log path)
    TMP_LOC=${f#*$SOURCE/fastq/}
    echo "[NOTE] cleaning up: rm "$SOURCE/${TMP_LOC%/*}"/trinity/jellyfish.kmers.fa"
    rm $SOURCE/${TMP_LOC%/*}/trinity/jellyfish.kmers.fa
else
    echo "[WARN] no cleaning up: I hope you plan on running some more analyses on the chrysalis/inchworm stuff. \
            Otherwise, free up some space you dirty, dirty hacker.... "
fi

echo -e "\n********* $CHECKPOINT\n"
################################################################################
#TODO correct dummy
#[ -e $OUTDIR/${n/%$READONE.$FASTQ/.$ASD.bam}.dummy ] && rm $OUTDIR/${n/%$READONE.$FASTQ/.$ASD.bam}.dummy
echo ">>>>> transcriptome assembly with trinity chrysalis - FINISHED"
echo ">>>>> enddate "`date`
