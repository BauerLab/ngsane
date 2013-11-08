#!/bin/bash -e

# Pindel script
# author: Denis C. Bauer
# date: Nov.2013

# messages to look out for -- relevant for the QC.sh script:
# QCVARIABLES,
# RESULTFILENAME <DIR>/<TASK>/<SAMPLE>

echo ">>>>> Structural Variant Calling with Pindel"
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"


function usage {
echo -e "usage: $(basename $0) -k NGSANE -f bam -o OUTDIR [OPTIONS]

Script running pindel
It expects a bam file (*.$ASD.bam)

required:
  -k | --toolkit <path>     location of the NGSANE repository 
  -f | --fastq <file>       bam file
  -o | --outdir <path>      output dir

options:
"
exit
}

if [ ! $# -gt 4 ]; then usage ; fi

#INPUTS
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the NGSANE repository
        -f | --bam )            shift; f=$1 ;; # bam file
        -o | --outdir )         shift; OUTDIR=$1 ;; # output dir
        -L | --region )         shift; SEQREG=$1 ;; # (optional) region of specific interest, e.g. targeted reseq
        --recover-from )        shift; RECOVERFROM=$1 ;; # attempt to recover from log file
        -h | --help )           usage ;;
        * )                     usage
    esac
    shift
done


#PROGRAMS
. $CONFIG
. ${NGSANE_BASE}/conf/header.sh
. $CONFIG

################################################################################
CHECKPOINT="programs"

for MODULE in $MODULE_PINDEL; do module load $MODULE; done  # save way to load modules that itself load other modules
export PATH=$PATH_RECAL:$PATH
module list
echo $PATH
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)

#echo "[NOTE] set java parameters"
#JAVAPARAMS="-Xmx"$(python -c "print int($MEMORY_RECAL*0.8)")"g -Djava.io.tmpdir="$TMP" -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1" 
#unset _JAVA_OPTIONS
#echo "JAVAPARAMS "$JAVAPARAMS

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
#echo -e "--JAVA        --\n" $(java -Xmx200m -version 2>&1)
#[ -z "$(which java)" ] && echo "[ERROR] no java detected" && exit 1
#echo -e "--samtools    --\n "$(samtools 2>&1 | head -n 3 | tail -n-2)
#[ -z "$(which samtools)" ] && echo "[ERROR] no samtools detected" && exit 1
echo -e "--Pindel       --\n "$(pindel --v | grep version)
[ -z "$(pindel)" ] && echo "[ERROR] no pindel detected" && exit 1
#echo -e "--igvtools    --\n "$(java $JAVAPARAMS -jar $PATH_IGVTOOLS/igvtools.jar version 2>&1)
#[ ! -f $PATH_IGVTOOLS/igvtools.jar ] && echo "[ERROR] no igvtools detected" && exit 1
#echo -e "--GATK        --\n "$(java $JAVAPARAMS -jar $PATH_GATK/GenomeAnalysisTK.jar --version)
#[ ! -f $PATH_GATK/GenomeAnalysisTK.jar ] && echo "[ERROR] no GATK detected" && exit 1

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="parameters"

# get basename of f
f=${f/%.dummy/} #if input came from pipe
n=${f##*/}

#BAMREADS=`head -n1 $f.stats | cut -d " " -f 1`

#if [ -n "$SEQREG" ]; then REGION="-L $SEQREG"; fi

#is paired ?
p=`grep "paired in " $f.stats | cut -d " " -f 1`
if [ ! "$p" -eq "0" ]; then
    PAIRED="1"
    echo "[NOTE] PAIRED"
else
    echo "[NOTE] SINGLE"
    PAIRED="0"
fi

NAME=$(basename $f)
NAME=${NAME/.$ASD.bam/}


# delete old bam files unless attempting to recover
if [ -z "$RECOVERFROM" ]; then
    [ -e $OUTDIR/$NAME.summary.txt ] && rm $OUTDIR/$NAME*
fi


echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="recall files from tape"

if [ -n "$DMGET" ]; then
	dmget -a $f
fi
    
echo -e "\n********* $CHECKPOINT\n"    

################################################################################
CHECKPOINT="create bam config file"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
    echo "[NOTE] $CHECKPOINT"

	# grep the median insert size for this file from the previously run Picard metrices
	command="grep MEDIAN_INSERT_SIZE  ${f/$INPUT_PINDEL/$INPUT_PINDEL\/metrices}.insert_size_metrics -A 1 | tail -n 1 | cut -f 1"
	echo $command 
	INSERT=$(eval $command)
	
	echo -e "$f\t$INSERT\t$NAME" > $OUTDIR/$NAME.conf.txt
        
    # mark checkpoint
    if [ -f $OUTDIR/$NAME.conf.txt ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi 

################################################################################
CHECKPOINT="pindel"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
    echo "[NOTE] $CHECKPOINT"

	command="pindel $PINDELADDPARAM -f $FASTA -i $OUTDIR/$NAME.conf.txt -c ALL -o $OUTDIR/$NAME -T $CPU_PINDEL"
	echo $command && eval $command


	for i in D SI LI INV TD BP ; do 
		grep "#" $OUTDIR/$NAME"_"$i -A 1 | grep -v "#" | grep -v "\-\-" ;
	done | sort -k8,8 -k 10,10n > $OUTDIR/$NAME.summary.txt
        
	rm $OUTDIR/$NAME.conf.txt

    # mark checkpoint
    if [ -f $OUTDIR/$NAME.summary.txt ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi 



################################################################################
[ -e $OUTDIR/$NAME.dummy ] && rm $OUTDIR/$NAME.dummy
echo ">>>>> Structural Variant Calling with Pindel - FINISHED"
echo ">>>>> enddate "`date`
