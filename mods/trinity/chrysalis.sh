#!/bin/bash 
# chrysalis_1.sh
# tested with trinityrnaseq_r2013-02-25


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

ulimit -s unlimited

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

mkdir $TMP/$JOB_ID
TMP_LOC=${f#*$SOURCE/fastq/}  # project name, i.e. sample1
cp $SOURCE/${TMP_LOC%/*}/trinity/ $TMP/$JOB_ID
cd $TMP/$JOB_ID

# we set --max_reads_per_graph to 1million because very high I/O is needed otherwise. It is unlikely that a transcript needs more than 1 million reads to be assembled....
if [ -z "$SS_LIBRARY_TYPE" ]; then
    echo "[WARNING] strand-specific RNAseq library type not set ($SS_LIBRARY_TYPE): treating input as non-stranded"
    RUN_COMMAND="Trinity.pl --seqType fq --left $f --right $f2 --output /tmp/$JOB_ID --JM $MEMORY_CHRYSALIS"G" --CPU $NCPU_CHRYSALIS --no_run_quantifygraph"
else
    echo "[NOTE] RNAseq library type: $SS_LIBRARY_TYPE"
    RUN_COMMAND="Trinity.pl --seqType fq --left $f --right $f2 --SS_lib_type $SS_LIBRARY_TYPE --output /tmp/$JOB_ID --JM $MEMORY_CHRYSALIS"G" --CPU $NCPU_CHRYSALIS --no_run_quantifygraph"
fi

mv $TMP/$JOB_ID/* $SOURCE/${TMP_LOC%/*}/trinity/

echo $RUN_COMMAND && eval $RUN_COMMAND
if [ $? -ne 0 ]; then
    echo "[ERROR] chrysalis.sh did not complete properly"
else
    echo "[NOTE] Chrysalis pt.1 has completed properly! thank Martin by buying him another beer"
fi



