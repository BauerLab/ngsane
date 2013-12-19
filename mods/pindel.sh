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
PATH_PICARD=$(dirname $(which MarkDuplicates.jar))
PATH_GATK=$(dirname $(which GenomeAnalysisTK.jar))

echo "[NOTE] set java parameters"
JAVAPARAMS="-Xmx"$(python -c "print int($MEMORY_PINDEL*0.8)")"g -Djava.io.tmpdir="$TMP" -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1" 
unset _JAVA_OPTIONS
echo "JAVAPARAMS "$JAVAPARAMS

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--JAVA        --\n" $(java -Xmx200m -version 2>&1)
[ -z "$(which java)" ] && echo "[ERROR] no java detected" && exit 1
#echo -e "--samtools    --\n "$(samtools 2>&1 | head -n 3 | tail -n-2)
#[ -z "$(which samtools)" ] && echo "[ERROR] no samtools detected" && exit 1
echo -e "--Pindel       --\n "$(pindel --v | grep version)
[ -z "$(pindel)" ] && echo "[ERROR] no pindel detected" && exit 1
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
NAME=${NAME/.$ASD.bam/}

REFERENCE_NAME=$(basename $FASTA)
REFERENCE_NAME=${REFERENCE_NAME/%.*/}
REFERENCE_ENDING=${FASTA/*./}

echo $REFERENCE_ENDING" "$REFERENCE_DATE" "$REFERENCE_NAME


# delete old bam files unless attempting to recover
if [ -z "$RECOVERFROM" ]; then
    [ -e $OUTDIR/$NAME"_"BD ] && rm $OUTDIR/$NAME*
fi


echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="recall files from tape"

if [ -n "$DMGET" ]; then
	dmget -a $f
fi
    
echo -e "\n********* $CHECKPOINT\n"    


################################################################################
CHECKPOINT="calculate inner distance"                                                                                                

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
    

	if [ ! -e ${f/$INPUT_PINDEL/$INPUT_PINDEL/metrices}.insert_size_metrics ]; then

		echo "[NOTE] ${f/$INPUT_PINDEL/$INPUT_PINDEL/metrices}.insert_size_metrics does not exist: start insert-size calculation"

		INDIR=${OUTDIR/$TASKPINDEL/$INPUT_PINDEL}
		#mkdir $INDIR/metrices


		echo "[NOTE] reorder to be conform with reference"
	    java $JAVAPARAMS -jar $PATH_PICARD/ReorderSam.jar \
	        INPUT=$f \
	        REFERENCE=$FASTA \
			VALIDATION_STRINGENCY=SILENT \
	        OUTPUT=$TMP/$NAME.bam


		echo "[NOTE] run picard metrices to get insert size"
	    export PATH=$PATH:/usr/bin/
	    THISTMP=$TMP/$NAME$RANDOM #mk tmp dir because picard writes none-unique files
	    mkdir -p $THISTMP
	    java $JAVAPARAMS -jar $PATH_PICARD/CollectMultipleMetrics.jar \
	        INPUT=$TMP/$NAME.bam  \
	        REFERENCE_SEQUENCE=$FASTA \
	        OUTPUT=${f/$INPUT_PINDEL/$INPUT_PINDEL\/metrices} \
	        VALIDATION_STRINGENCY=SILENT \
	        PROGRAM=CollectAlignmentSummaryMetrics \
	        PROGRAM=CollectInsertSizeMetrics \
	        PROGRAM=QualityScoreDistribution \
	        TMP_DIR=$THISTMP
	    for im in $( ls $INDIR/metrices/*.pdf ); do
	        convert $im ${im/pdf/jpg}
	    done
	    [ -d $THISTMP ] && rm -r $THISTMP && rm $TMP/$NAME.bam

		
	fi

    # mark checkpoint
    if [ -f ${f/$INPUT_PINDEL/$INPUT_PINDEL/metrices}.insert_size_metrics ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
fi




################################################################################
CHECKPOINT="create bam config file"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
    echo "[NOTE] $CHECKPOINT"

	if [ ! -e $f.bai ]; then
		echo "[NOTE] create index for bam file"
		samtools index $f
	fi

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


	if [ "$INPUT_IS_BWA"=="false" ]; then
		echo "[NOTE] need to convert non-bwa input"
		command="samtools view $f | sam2pindel - $OUTDIR/$NAME.s2b.txt $INSERT $NAME 0 $INPUT_TYPE"
		echo $command && eval $command

		echo "[NOTE] run Pindel witn non-bwa alignment"
		command="pindel $PINDELADDPARAM -f $FASTA -p $OUTDIR/$NAME.s2b.txt -c ALL -o $OUTDIR/$NAME -T $CPU_PINDEL"
		echo $command && eval $command

	else
		echo "[NOTE] run Pindel with bwa alignment"
		command="pindel $PINDELADDPARAM -f $FASTA -i $OUTDIR/$NAME.conf.txt -c ALL -o $OUTDIR/$NAME -T $CPU_PINDEL"
		echo $command && eval $command

	fi

    # mark checkpoint
    if [ -f $OUTDIR/$NAME"_"BP ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi 


################################################################################
CHECKPOINT="convert to vcf and join"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
    echo "[NOTE] $CHECKPOINT"

	if [ ! -e ${FASTA/REFERENCE_ENDING/fa} ]; then echo "[ERROR] pindel2vcf needs fasta file with .fa ending - please create a symbolic link"; exit -1; fi 

	echo "[NOTE] convert to vcf"
	VAR=""
	for i in D SI LI INV TD BP ; do 
		pindel2vcf -G -p $OUTDIR/$NAME"_"$i -r ${FASTA/fasta/fa} -d $REFERENCE_DATE -R $REFERENCE_NAME -v $OUTDIR/$NAME"_"$i.vcf 
		VAR=$VAR" -V "$OUTDIR"/"$NAME"_"$i".vcf"
	done


	echo "[NOTE] join to single vcf"
   	command="java $JAVAPARAMS -cp $PATH_GATK/GenomeAnalysisTK.jar org.broadinstitute.sting.tools.CatVariants \
		-R $FASTA $VAR -out $OUTDIR/$NAME.summary.unsorted.vcf"
    echo $command && eval $command

	echo "[NOTE] order to reference"
    perl ${NGSANE_BASE}/tools/vcfsorter.pl ${FASTA/$REFERENCE_ENDING/dict} $OUTDIR/$NAME.summary.unsorted.vcf > $OUTDIR/$NAME.summary.sorteduncorr.vcf

	echo "[NOTE] correct END values from pindel2vcf errors (http://www.biostars.org/p/88908/)"
   	python ${NGSANE_BASE}/tools/checkpindelENDs.py $OUTDIR/$NAME.summary.sorteduncorr.vcf $OUTDIR/$NAME.summary.vcf
        
    # mark checkpoint
    if [ -f $OUTDIR/$NAME.summary.vcf ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

	rm $OUTDIR/$NAME.conf.txt $OUTDIR/$NAME.summary.sorteduncorr.vcf $OUTDIR/$NAME.summary.unsorted.vcf
	for i in D SI LI INV TD BP ; do rm $OUTDIR/$NAME"_"$i.vcf; done

fi 



################################################################################
[ -e $OUTDIR/$NAME.dummy ] && rm $OUTDIR/$NAME.dummy
echo ">>>>> Structural Variant Calling with Pindel - FINISHED"
echo ">>>>> enddate "`date`
