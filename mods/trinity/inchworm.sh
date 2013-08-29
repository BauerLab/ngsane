#!/bin/bash 
# inchworm.sh
# tested with trinityrnaseq_r2013-02-25

# $GE_TASK_ID

FORCESINGLE=0
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the configuration file
        -f | --fastq )          shift; f=$1 ;; # fastq file
        -o | --outdir )         shift; MYOUT=$1 ;; # output dir
        --forceSingle )         FORCESINGLE=1;;
        * )                     echo "don't understand "$1
    esac
    shift
done

. $CONFIG
. ${NGSANE_BASE}/conf/header.sh
. $CONFIG

## TO DO : append $PATH for systems without modules for bowtie, java, etc...
for MODS in $MODULES_TRINITY; do module load $MODS; done 

#if [ $NODETYPE_INCHWORM -eq "intel.q" ] ; then
	#some fancy optimisation stuff for intel nodes. Consider the intel compiler! 
#	echo "[NOTE] forcing copact thread placement on intel nodes"
#	export OMP_NUM_THREADS=$NCPU_INCHWORM
#	export KMP_AFFINITY=compact
#fi

echo "********* detect library"
## is paired ?                                                                                                      
if [ -e ${f/$READONE/$READTWO} ] && [ "$FORCESINGLE" = 0 ]; then
    PAIRED="1"
    f2=${f/$READONE/$READTWO}
    echo "[NOTE] Paired library detected: "$f" + "$f2
else
    PAIRED="0"
    echo "[NOTE] Single-Strand (unpaired) library detected !!!! currently unsupported !!!!"
    exit 1
fi

# cd /tmp ##??????

# this runs Inchworm only
# we set --max_reads_per_graph to 1million because very high I/O is needed otherwise. It is unlikely that a transcript needs more than 1 million reads to be assembled....
if [ -z "$SS_LIBRARY_TYPE" ]; then
    echo "[WARNING] strand-specific RNAseq library type not set ($SS_LIBRARY_TYPE): treating input as non-stranded"
    RUN_COMMAND="Trinity.pl --seqType fq --left $f --right $f2 --max_reads_per_graph 1000000 --output $MYOUT --JM $MEMORY_INCHWORM"G" --CPU $NCPU_INCHWORM --no_run_chrysalis"
else
    echo "[NOTE] RNAseq library type: $SS_LIBRARY_TYPE"
    RUN_COMMAND="Trinity.pl --seqType fq --left $f --right $f2 --SS_lib_type $SS_LIBRARY_TYPE --max_reads_per_graph 1000000 --output $MYOUT --JM $MEMORY_INCHWORM"G" --CPU $NCPU_INCHWORM --no_run_chrysalis"
fi
echo $RUN_COMMAND && eval $RUN_COMMAND
if [ $? -ne 0 ]; then
    echo "[ERROR] inchworm.sh did not complete properly"
else
    echo "[NOTE] inchworm has completed properly! thank Martin by buying him a beer"
fi

if [ -n $CLEANUP ]; then 
#    rm $MYOUT/*ebwt $MYOUT/*ebwt  #bowtie files required for chrysalis, and this is the WRONG file location (log path)
    TMP_LOC=${f#*$SOURCE/fastq/}
    echo "[NOTE] cleaning up: rm "$SOURCE/${TMP_LOC%/*}"/trinity/jellyfish.kmers.fa"
    rm $SOURCE/${TMP_LOC%/*}/trinity/jellyfish.kmers.fa
fi


