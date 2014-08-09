#!/bin/bash -e

# Pindel script
# author: Denis C. Bauer
# date: Nov.2013

# messages to look out for -- relevant for the QC.sh script:
# QCVARIABLES,config not found
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
It expects a bam file (*$ASD.bam)

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
        --recover-from )        shift; NGSANE_RECOVERFROM=$1 ;; # attempt to recover from log file
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
NGSANE_CHECKPOINT_INIT "programs"

 # save way to load modules that itself load other modules
hash module 2>/dev/null && for MODULE in $MODULE_PINDEL; do module load $MODULE; done && module list

export PATH=$PATH_RECAL:$PATH
echo $PATH
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)
[ -z "$PATH_PICARD" ] && PATH_PICARD=$(dirname $(which ReorderSam.jar))
[ -z "$PATH_GATK" ] && PATH_GATK=$(dirname $(which GenomeAnalysisTK.jar))

echo "[NOTE] set java parameters"
JAVAPARAMS="-Xmx"$(python -c "print int($MEMORY_PINDEL*0.8)")"g -Djava.io.tmpdir="$TMP" -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1" 
unset _JAVA_OPTIONS
echo "JAVAPARAMS "$JAVAPARAMS

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--JAVA        --\n" $(java -Xmx200m -version 2>&1)
[ -z "$(which java)" ] && echo "[ERROR] no java detected" && exit 1
echo -e "--samtools    --\n "$(samtools 2>&1 | head -n 3 | tail -n-2)
[ -z "$(which samtools)" ] && echo "[ERROR] no samtools detected" && exit 1
echo -e "--Pindel       --\n "$(pindel --v | grep version)
[ -z "$(pindel)" ] && echo "[ERROR] no pindel detected" && exit 1
echo -e "--PICARD      --\n "$(java $JAVAPARAMS -jar $PATH_PICARD/ReorderSam.jar --version 2>&1)
[ ! -f $PATH_PICARD/ReorderSam.jar ] && echo "[ERROR] no picard detected" && exit 1
echo -e "--GATK        --\n "$(java $JAVAPARAMS -jar $PATH_GATK/GenomeAnalysisTK.jar --version)
[ ! -f $PATH_GATK/GenomeAnalysisTK.jar ] && echo "[ERROR] no GATK detected" && exit 1

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "parameters"

# get basename of f
f=${f/%.dummy/} #if input came from pipe
n=${f##*/}
SAMPLE=${n/%$ASD.bam/}

REFERENCE_NAME=$(basename $FASTA)
REFERENCE_NAME=${REFERENCE_NAME/%.*/}
REFERENCE_ENDING=${FASTA/*./}

echo $REFERENCE_ENDING" "$REFERENCE_DATE" "$REFERENCE_NAME


# delete old bam files unless attempting to recover
if [ -z "$NGSANE_RECOVERFROM" ]; then
    [ -e $OUTDIR/$SAMPLE"_"BD ] && rm $OUTDIR/$SAMPLE*
fi

if [ ! -e ${FASTA/$REFERENCE_ENDING/fa} ]; then echo "[ERROR] pindel2vcf needs fasta file with .fa ending - please create a symbolic link"; exit -1; fi 
if [ ! -e ${FASTA/$REFERENCE_ENDING/fasta} ]; then echo "[ERROR] GenomeAnalysisTK.jar needs fasta file with .fasta ending - please create a symbolic link"; exit -1; fi 


NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "recall files from tape"

if [ -n "$DMGET" ]; then
	dmget -a $f
fi
    
NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "calculate inner distance"                                                                                                

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then
    
	INDIR=${OUTDIR/-$TASK_PINDEL/}

	if [ ! -e $INDIR/metrices/$n.insert_size_metrics ]; then

		echo "[NOTE] $INDIR/metrices/$n.insert_size_metrics does not exist: start insert-size calculation"
		if [ ! -d $INDIR/metrices ]; then mkdir -p $INDIR/metrices; fi


	    export PATH=$PATH:/usr/bin/
	    THISTMP=$TMP/$SAMPLE$RANDOM #mk tmp dir because picard writes none-unique files
	    mkdir -p $THISTMP

		echo "[NOTE] reorder to be conform with reference"
	    java $JAVAPARAMS -jar $PATH_PICARD/ReorderSam.jar \
	        INPUT=$f \
	        REFERENCE=$FASTA \
			VALIDATION_STRINGENCY=SILENT \
	        OUTPUT=$THISTMP/$SAMPLE.bam


		echo "[NOTE] run picard metrices to get insert size"
	    java $JAVAPARAMS -jar $PATH_PICARD/CollectMultipleMetrics.jar \
	        INPUT=$THISTMP/$SAMPLE.bam  \
	        REFERENCE_SEQUENCE=$FASTA \
	        OUTPUT=$INDIR/metrices/$n \
	        VALIDATION_STRINGENCY=SILENT \
	        PROGRAM=CollectAlignmentSummaryMetrics \
	        PROGRAM=CollectInsertSizeMetrics \
	        PROGRAM=QualityScoreDistribution \
	        TMP_DIR=$THISTMP
	    for im in $( ls $INDIR/metrices/*.pdf ); do
	        convert $im ${im/pdf/jpg}
	    done
	    [ -d $THISTMP ] && rm -r $THISTMP

		
	fi

    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK ${f/$INPUT_PINDEL/$INPUT_PINDEL/metrices}.insert_size_metrics
fi

################################################################################
NGSANE_CHECKPOINT_INIT "create bam config file"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

	if [ ! -e $f.bai ]; then
		echo "[NOTE] create index for bam file"
		samtools index $f
	fi

	# grep the median insert size for this file from the previously run Picard metrices
	if [ -z "$INSERT" ]; then # sometimes Picard insert size calculation is wrong...
		command="grep MEDIAN_INSERT_SIZE  ${f/$INPUT_PINDEL/$INPUT_PINDEL/metrices}.insert_size_metrics -A 1 | tail -n 1 | cut -f 1"
		echo $command 
		INSERT=$(eval $command)
	fi
	
	echo -e "$f\t$INSERT\t$SAMPLE" > $OUTDIR/$SAMPLE.conf.txt
        
    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE.conf.txt

fi 

################################################################################
NGSANE_CHECKPOINT_INIT "pindel"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

	if [[ ${INPUT_IS_BWA}="false" ]]; then

		if [ -z "$INSERT" ]; then # sometimes Picard insert size calculation is wrong...
			INSERT=$(grep MEDIAN_INSERT_SIZE  ${f/$INPUT_PINDEL/$INPUT_PINDEL/metrices}.insert_size_metrics -A 1 | tail -n 1 | cut -f 1)
		fi

		echo "[NOTE] need to convert non-bwa input"
		command="samtools view $f | sam2pindel - $OUTDIR/$SAMPLE.s2b.txt $INSERT $SAMPLE 0 $INPUT_TYPE"
		echo $command && eval $command

		echo "[NOTE] run Pindel witn non-bwa alignment"
		command="pindel $PINDELADDPARAM -f $FASTA -p $OUTDIR/$SAMPLE.s2b.txt -c ALL -o $OUTDIR/$SAMPLE -T $CPU_PINDEL"
		echo $command && eval $command

	else
		echo "[NOTE] run Pindel with bwa alignment"
		command="pindel $PINDELADDPARAM -f $FASTA -i $OUTDIR/$SAMPLE.conf.txt -c ALL -o $OUTDIR/$SAMPLE -T $CPU_PINDEL"
		echo $command && eval $command

	fi

    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE"_"BP

fi 

################################################################################
NGSANE_CHECKPOINT_INIT "convert to vcf"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then
    echo "[NOTE] $CHECKPOINT"

	for i in D SI LI INV TD BP ; do 
		pindel2vcf -G -p $OUTDIR/$SAMPLE"_"$i -r ${FASTA/$REFERENCE_ENDING/fa} -d $REFERENCE_DATE -R $REFERENCE_NAME -v $OUTDIR/$SAMPLE"_"$i.vcf 
	done

    # mark checkpoint
	if [ -f $OUTDIR/$SAMPLE"_BP".vcf ];then echo -e "\n********* $CHECKPOINT\n"; unset NGSANE_RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi 


################################################################################
NGSANE_CHECKPOINT_INIT "join vcf files"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

	VAR=""
	for i in D SI LI INV TD BP ; do  VAR=$VAR" -V "$OUTDIR"/"$SAMPLE"_"$i".vcf"; done

	echo "[NOTE] join to single vcf"
   	command="java $JAVAPARAMS -cp $PATH_GATK/GenomeAnalysisTK.jar org.broadinstitute.sting.tools.CatVariants \
		-R ${FASTA/$REFERENCE_ENDING/fasta} $VAR -out $OUTDIR/$SAMPLE.summary.unsorted.vcf"
    echo $command && eval $command

	echo "[NOTE] order to reference"
    perl ${NGSANE_BASE}/tools/vcfsorter.pl ${FASTA/$REFERENCE_ENDING/dict} $OUTDIR/$SAMPLE.summary.unsorted.vcf > $OUTDIR/$SAMPLE.summary.sorteduncorr.vcf

	echo "[NOTE] correct END values from pindel2vcf errors (http://www.biostars.org/p/88908/)"
   	python ${NGSANE_BASE}/tools/checkpindelENDs.py $OUTDIR/$SAMPLE.summary.sorteduncorr.vcf $OUTDIR/$SAMPLE$ASD.vcf
        
    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE$ASD.vcf

	rm $OUTDIR/$SAMPLE.conf.txt $OUTDIR/$SAMPLE.summary.sorteduncorr.vcf $OUTDIR/$SAMPLE.summary.unsorted.vcf*
	for i in D SI LI INV TD BP ; do rm $OUTDIR/$SAMPLE"_"$i.vcf; done

fi 

################################################################################
[ -e $OUTDIR/$SAMPLE.dummy ] && rm $OUTDIR/$SAMPLE.dummy
echo ">>>>> Structural Variant Calling with Pindel - FINISHED"
echo ">>>>> enddate "`date`
