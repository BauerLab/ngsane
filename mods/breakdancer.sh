#!/bin/bash -e

# BREAKDANCER script
# author: Denis C. Bauer
# date: Jan.2014

# messages to look out for -- relevant for the QC.sh script:
# QCVARIABLES,config not found
# RESULTFILENAME <DIR>/<TASK>/<SAMPLE>

echo ">>>>> Structural Variant Calling with Breakdancer"
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"


function usage {
echo -e "usage: $(basename $0) -k NGSANE -f bam -o OUTDIR [OPTIONS]

Script running breakdancer
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

for MODULE in $MODULE_BREAKDANCER; do module load $MODULE; done  # save way to load modules that itself load other modules
export PATH=$PATH_RECAL:$PATH
module list
echo $PATH
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)
PATH_PICARD=$(dirname $(which MarkDuplicates.jar))
PATH_GATK=$(dirname $(which GenomeAnalysisTK.jar))
PATH_BREAKDANCER=$(dirname $(which breakdancer-max))


echo "[NOTE] set java parameters"
JAVAPARAMS="-Xmx"$(python -c "print int($MEMORY_BREAKDANCER*0.8)")"g -Djava.io.tmpdir="$TMP" -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1" 
unset _JAVA_OPTIONS
echo "JAVAPARAMS "$JAVAPARAMS

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--JAVA        --\n" $(java -Xmx200m -version 2>&1)
[ -z "$(which java)" ] && echo "[ERROR] no java detected" && exit 1
#echo -e "--samtools    --\n "$(samtools 2>&1 | head -n 3 | tail -n-2)
#[ -z "$(which samtools)" ] && echo "[ERROR] no samtools detected" && exit 1
echo -e "--Breakdancer       --\n "$(breakdancer-max 2>&1 | grep version)
[ -z "$(which breakdancer-max)" ] && echo "[ERROR] no BREAKDANCER detected" && exit 1
echo -e "--PICARD      --\n "$(java $JAVAPARAMS -jar $PATH_PICARD/MarkDuplicates.jar --version 2>&1)
[ ! -f $PATH_PICARD/MarkDuplicates.jar ] && echo "[ERROR] no picard detected" && exit 1
echo -e "--GATK        --\n "$(java $JAVAPARAMS -jar $PATH_GATK/GenomeAnalysisTK.jar --version)
[ ! -f $PATH_GATK/GenomeAnalysisTK.jar ] && echo "[ERROR] no GATK detected" && exit 1

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="parameters"

# get basename of f
f=${f/%.dummy/} #if input came from pipe
n=${f##*/}

NAME=$(basename $f)
NAME=${NAME/.bam/}

REFERENCE_NAME=$(basename $FASTA)
REFERENCE_NAME=${REFERENCE_NAME/%.*/}
REFERENCE_ENDING=${FASTA/*./}

echo $REFERENCE_ENDING" "$REFERENCE_DATE" "$REFERENCE_NAME


# delete old bam files unless attempting to recover
if [ -z "$RECOVERFROM" ]; then
    [ -e $OUTDIR/$NAME.SV ] && rm $OUTDIR/$NAME.*
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

#	if [ ! -e $f.bai ]; then
#		echo "[NOTE] create index for bam file"
#		samtools index $f
#	fi

	cd $OUTDIR
	command="perl ${PATH_BREAKDANCER}/bam2cfg.pl $BREAKDANCERCONFIGADDPARAM -g -h ${f} > $OUTDIR/$NAME.cfg"
	echo $command && eval $command 
	cd ../
	
    # mark checkpoint
    if [ -f $OUTDIR/$NAME.cfg ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi 


################################################################################
CHECKPOINT="Run breakdancer"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
    echo "[NOTE] $CHECKPOINT"

	# -d $OUTDIR/$NAME
	command="breakdancer-max $BREAKDANCERADDPARAM -g $OUTDIR/$NAME.bed $OUTDIR/$NAME.cfg > $OUTDIR/$NAME.SV"
	echo $command && eval $command


    # mark checkpoint
    if [ -f $OUTDIR/$NAME.bed ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi 


#############################################################################
#CHECKPOINT="convert to vcf"
#
#if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
#    echo "::::::::: passed $CHECKPOINT"
#else 
#    echo "[NOTE] $CHECKPOINT"
#
#	python ${NGSANE_BASE}/tools/breakdancer2vcf.py -i $OUTDIR/$NAME.SV -o $OUTDIR/$NAME.vcf  
#
#   #mark checkpoint
#	if [ -f $OUTDIR/$NAME.vcf ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
#
#fi 



################################################################################
[ -e $OUTDIR/$NAME.dummy ] && rm $OUTDIR/$NAME.dummy
echo ">>>>> Structural Variant Calling with BREAKDANCER - FINISHED"
echo ">>>>> enddate "`date`
