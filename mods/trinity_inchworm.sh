#!/bin/bash -e
# Martin Smith, August 2013
# tested with trinityrnaseq_r2013-02-25

# QCVARIABLES,Resource temporarily unavailable
# RESULTFILENAME <DIR>/<TASK>/<SAMPLE>.inchworm.timing

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
hash module 2>/dev/null && for MODULE in $MODULES_TRINITY; do module load $MODULE; done  && module list 

export PATH=$PATH_TRINITY:$PATH
echo "PATH=$PATH"

echo "[NOTE] set java parameters"
JAVAPARAMS="-Xmx"$(python -c "print int($MEMORY_INCHWORM*0.75)")"g -Djava.io.tmpdir="$TMP" -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1" 
unset _JAVA_OPTIONS
echo "JAVAPARAMS "$JAVAPARAMS

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--JAVA        --\n" $(java -Xmx200m -version 2>&1)
[ -z "$(which java)" ] && echo "[ERROR] no java detected" && exit 1
echo -e "--bowtie      --\n "$(bowtie --version)
[ -z "$(which bowtie)" ] && echo "[ERROR] no bowtie detected" && exit 1
echo -e "--perl        --\n "$(perl -v | grep "version" )
[ -z "$(which perl)" ] && echo "[ERROR] no perl detected" && exit 1
echo -e "--trinity     --\n "$(Trinity.pl --version)
[ -z "$(which Trinity.pl)" ] && echo "[ERROR] no trinity detected" && exit 1


#if [ $NODETYPE_INCHWORM -eq "intel.q" ] ; then
	#some fancy optimisation stuff for intel nodes. Consider the intel compiler! 
#	echo "[NOTE] forcing copact thread placement on intel nodes"
#	export OMP_NUM_THREADS=$NCPU_INCHWORM
#	export KMP_AFFINITY=compact
#fi

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "detect library"

# get basename of f
n=${f##*/}
SAMPLE=${n/%$READONE.$FASTQ/}

## is paired ?                                                                                                      
if [ -e ${f/%$READONE.$FASTQ/$READTWO.$FASTQ} ]; then
    PAIRED="1"
    f2=${f/%$READONE.$FASTQ/$READTWO.$FASTQ}
    echo "[NOTE] Paired library detected: "$f" + "$f2
else
    PAIRED="0"
    echo "[ERROR] Single-Strand (unpaired) library detected !!!! currently unsupported !!!!"
    exit 1
fi

# fancy symbolic link generation to work in common trinity folder
mkdir -p $OUTDIR/../$TASK_TRINITY/$SAMPLE
ln -f -s ../$TASK_TRINITY/$SAMPLE $OUTDIR/$SAMPLE

# make sure we use the same tmp folder name all time so stick with one jobid
if [ ! -f $OUTDIR/$SAMPLE/persistent_id.tmp ]; then 
    PERSISTENT_ID="trinity_$JOB_ID"
    echo "$PERSISTENT_ID" >  $OUTDIR/$SAMPLE/persistent_id.tmp;
    echo "[NOTE] Persistent ID set to: $PERSISTENT_ID"
else
    PERSISTENT_ID="$(head -n 1 $OUTDIR/$SAMPLE/persistent_id.tmp)"
    echo "[NOTE] Persistent ID found: $PERSISTENT_ID"
fi

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "recall files from tape"

if [ -n "$DMGET" ]; then
	dmget -a ${f/$READONE/"*"}
	dmget -a $OUTDIR/$SAMPLE/*
fi
    
NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "Inchworm"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    echo "[NOTE] --max_reads_per_graph set to 1 million because very high I/O is needed otherwise. It is unlikely that a transcript needs more than 1 million reads to be assembled"

    if [ -z "$SS_LIBRARY_TYPE" ]; then
        echo "[WARNING] strand-specific RNAseq library type not set ($SS_LIBRARY_TYPE): treating input as non-stranded"
        RUN_COMMAND="$(which perl) $(which Trinity.pl) --seqType fq --left $f --right $f2 --max_reads_per_graph 1000000 --output $OUTDIR/$SAMPLE --JM $MEMORY_INCHWORM"G" --CPU $NCPU_INCHWORM --no_run_chrysalis"
    else
        echo "[NOTE] RNAseq library type: $SS_LIBRARY_TYPE"
        RUN_COMMAND="$(which perl) $(which Trinity.pl) --seqType fq --left $f --right $f2 --SS_lib_type $SS_LIBRARY_TYPE --max_reads_per_graph 1000000 --output $OUTDIR/$SAMPLE --JM $MEMORY_INCHWORM"G" --CPU $NCPU_INCHWORM --no_run_chrysalis"
    fi
    echo $RUN_COMMAND && eval $RUN_COMMAND
    echo "[NOTE] inchworm has completed properly! thank Martin by buying him a beer"

    cp $OUTDIR/$SAMPLE/Trinity.timing $OUTDIR/$SAMPLE.inchworm.timing

    echo "jellyfish kmers: "$(grep ">" $OUTDIR/$SAMPLE/jellyfish.kmers.fa |  wc -l) > $OUTDIR/$SAMPLE.summary.txt
    echo "inchworm fasta: "$(grep ">" <$(ls $OUTDIR/$SAMPLE/inchworm.*.fa | head -n 1) | wc -l) >> $OUTDIR/$SAMPLE.summary.txt
    echo "both fasta: "$(grep ">" $OUTDIR/$SAMPLE/both.fa |  wc -l) >> $OUTDIR/$SAMPLE.summary.txt
    
    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE.inchworm.timing $OUTDIR/$SAMPLE.summary.txt
fi 

################################################################################
NGSANE_CHECKPOINT_INIT "cleanup"

if [ -n $CLEANUP ]; then 
#    rm $OUTDIR/*ebwt $OUTDIR/*ebwt  #bowtie files required for chrysalis, and this is the WRONG file location (log path)
    [ -e $OUTDIR/$SAMPLE/jellyfish.kmers.fa ] && rm $OUTDIR/$SAMPLE/jellyfish.kmers.fa
else
    echo "[WARN] no cleaning up: I hope you plan on running some more analyses on the chrysalis/inchworm stuff. \
            Otherwise, free up some space you dirty, dirty hacker.... "
fi

NGSANE_CHECKPOINT_CHECK
################################################################################
#TODO correct dummy
[ -e $OUTDIR/$SAMPLE.inchworm.timing.dummy ] && rm $OUTDIR/$SAMPLE.inchworm.timing.dummy
echo ">>>>> transcriptome assembly with trinity inchworm - FINISHED"
echo ">>>>> enddate "`date`
