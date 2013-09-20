#!/bin/bash -e
# Martin Smith, August 2013
# tested with trinityrnaseq_r2013-02-25

# QCVARIABLES,Resource temporarily unavailable
# RESULTFILENAME <DIR>/<TASK>/<SAMPLE>.fasta.gz

echo ">>>>> transcriptome assembly with trinity butterfly"
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

################################################################################
CHECKPOINT="programs"
for MODS in $MODULES_TRINITY; do module load $MODS; done 
export PATH=$PATH_TRINITY:$PATH
module list
echo "PATH=$PATH"

echo "[NOTE] set java parameters"
JAVAPARAMS="-Xmx"$(python -c "print int($MEMORY_BUTTERFLY*0.8)")"g -Djava.io.tmpdir="$TMP" -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1" 
unset _JAVA_OPTIONS
echo "JAVAPARAMS "$JAVAPARAMS

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--JAVA        --\n" $(java -Xmx200m -version 2>&1)
[ -z "$(which java)" ] && echo "[ERROR] no java detected" && exit 1
echo -e "--bowtie      --\n "$(bowtie --version)
[ -z "$(which bowtie)" ] && echo "[ERROR] no bowtie detected" && exit 1
echo -e "--perl        --\n "$(perl -v | grep "This is perl" )
[ -z "$(which perl)" ] && echo "[ERROR] no perl detected" && exit 1
echo -e "--trinity     --\n "$(Trinity.pl --version)
[ -z "$(which Trinity.pl)" ] && echo "[ERROR] no trinity detected" && exit 1

ulimit -s unlimited

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="detect library"

# get basename of f
n=${f##*/}
SAMPLE=${n/%$READONE.$FASTQ/}

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

# fancy symbolic link generation to work in common trinity folder
mkdir -p $OUTDIR/../$TASKTRINITY/$SAMPLE
ln -f -s ../$TASKTRINITY/$SAMPLE $OUTDIR/$SAMPLE

# make sure we use the same tmp folder name all time so stick with one jobid
if [ ! -f $OUTDIR/$SAMPLE/persistent_id.tmp ]; then 
    echo "[ERROR] persistent ID not detected"
    exit 1
else
    PERSISTENT_ID="$(head -n 1 $OUTDIR/$SAMPLE/persistent_id.tmp)"
    echo "[NOTE] Persistent ID found: $PERSISTENT_ID"
fi

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="recall files from tape"

if [ -n "$DMGET" ]; then
	dmget -a ${f/$READONE/"*"}
	dmget -a $OUTDIR/$SAMPLE/*
fi
    
echo -e "\n********* $CHECKPOINT\n"
###############################################################################
CHECKPOINT="Butterfly"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

#   edit files to replace temporary file paths to local file paths
    sed -i "s?${TMP}\/${PERSISTENT_ID}?${OUTDIR}\/${SAMPLE}?g"  $OUTDIR/../$TASKTRINITY/$SAMPLE/chrysalis/quantifyGraph_commands $OUTDIR/../$TASKTRINITY/$SAMPLE/chrysalis/butterfly_commands $OUTDIR/../$TASKTRINITY/$SAMPLE/chrysalis/component_file_listing.txt

    # this runs buterfly only
    echo "[NOTE] --max_reads_per_graph set to 1 million because very high I/O is needed otherwise. It is unlikely that a transcript needs more than 1 million reads to be assembled"
    
    if [ -z "$SS_LIBRARY_TYPE" ]; then
        echo "[WARNING] strand-specific RNAseq library type not set ($SS_LIBRARY_TYPE): treating input as non-stranded"
        RUN_COMMAND="$(which perl) $(which Trinity.pl) --seqType fq --left $f --right $f2 --max_reads_per_graph 1000000 --output $OUTDIR/../$TASKTRINITY/$SAMPLE  \
                        --JM $MEMORY_BUTTERFLY"G" --CPU $NCPU_BUTTERFLY "
    else
        echo "[NOTE] RNAseq library type: $SS_LIBRARY_TYPE"
        RUN_COMMAND="$(which perl) $(which Trinity.pl) --seqType fq --left $f --right $f2 --SS_lib_type $SS_LIBRARY_TYPE --max_reads_per_graph 1000000  \
                        --output $OUTDIR/../$TASKTRINITY/$SAMPLE --JM $MEMORY_BUTTERFLY"G" --CPU $NCPU_BUTTERFLY "
    fi
    echo $RUN_COMMAND && eval $RUN_COMMAND

    echo "[NOTE] butterly has completed properly! thank Martin by buying him yet another beer (hic!)"

    # mark checkpoint
    if [ -f $OUTDIR/$SAMPLE/Trinity.fasta ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
fi 

################################################################################
#CHECKPOINT="statistics"
#if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
#    echo "::::::::: passed $CHECKPOINT"
#else 
#    RUN_COMMAND="$(which perl) $(which Trinity.pl) --seqType fq --left $f --right $f2 --SS_lib_type $SS_LIBRARY_TYPE --max_reads_per_graph 1000000 \
#                        --output $OUTDIR/../$TASKTRINITY/$SAMPLE --JM $MEMORY_BUTTERFLY"G" --CPU $NCPU_BUTTERFLY "
#    echo $RUN_COMMAND && eval $RUN_COMMAND
#fi

################################################################################
CHECKPOINT="zip"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

    echo "Transcripts: "$(grep ">" $OUTDIR/$SAMPLE/Trinity.fasta | wc -l) > $OUTDIR/${n/%$READONE.$FASTQ/.summary.txt}
    $GZIP $OUTDIR/$SAMPLE/Trinity.fasta 
    # put into trinity  <<<-vvv--- this may cause problems with downstream components
    mv $OUTDIR/$SAMPLE/Trinity.fasta.gz $OUTDIR/../$TASKTRINITY/$SAMPLE.fasta.gz

    # mark checkpoint
    if [ -f $OUTDIR/../$TASKTRINITY/$SAMPLE.fasta.gz ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
fi 
################################################################################
CHECKPOINT="cleanup"

if [ -n $CLEANUP ]; then 
    echo "[WARNING] cleaning up, all intermediary files will disappear: rm $OUTDIR/$SAMPLE/... \"POOF\""
    cd $OUTDIR/$SAMPLE /
    ls ./ | grep -v "Trinity.*" | xargs -d"\n" rm -r 
#    mv $SOURCE/${TMP_LOC%/*}.fasta.gz $SOURCE/${TMP_LOC%/*}/trinity/ 

else
    echo "[NOTE] no cleaning up: I hope you plan on running some more analyses on the chrysalis/inchworm stuff. \
            Otherwise, free up some space you dirty, dirty hacker.... "
    # mv /tmp/$JOB_ID/* $SOURCE/${TMP_LOC%/*}/trinity/
fi

echo -e "\n********* $CHECKPOINT\n"
################################################################################
[ -e $OUTDIR/$SAMPLE.fasta.gz.dummy ] && rm $OUTDIR/$SAMPLE.fasta.gz.dummy
echo ">>>>> transcriptome assembly with trinity butterfly - FINISHED"
echo ">>>>> enddate "`date`
