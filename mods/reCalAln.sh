#!/bin/bash -e

# Recalibation script; this is for gatk v2.5
# author: Denis C. Bauer
# date: May.2013

# messages to look out for -- relevant for the QC.sh script:
# QCVARIABLES,We are loosing reads,for unmapped read,no such file,file not found,reCalAln.sh: line
# RESULTFILENAME <DIR>/<TASK>/<SAMPLE>$ASR.bam

echo ">>>>> recalibration and realignment using GATK"
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"


function usage {
echo -e "usage: $(basename $0) -k NGSANE -f bam -r REFERENCE -o OUTDIR -d SNPDB [OPTIONS]

Script running the recalibration and realigment step (GATK)
It expects a bam file (*$ASD.bam)

required:
  -k | --toolkit <path>     location of the NGSANE repository 
  -f | --fastq <file>       bam file
  -r | --reference <file>   reference genome
  -o | --outdir <path>      output dir
  -d | --snpdb <path>       path to a SNPdb instance (rod)

options:
  -L | --region <ps>        region of specific interest, e.g. targeted reseq
                             format chr:pos-pos
"
exit
}

if [ ! $# -gt 4 ]; then usage ; fi

#INPUTS
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the NGSANE repository
        -f | --bam )            shift; f=$1 ;; # bam file
        -r | --reference )      shift; FASTA=$1 ;; # reference genome
        -d | --snpdb )          shift; DBSNPVCF=$1 ;; # snpdb
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

# save way to load modules that itself loads other modules
hash module 2>/dev/null && for MODULE in $MODULE_RECAL; do module load $MODULE; done && module list 

export PATH=$PATH_RECAL:$PATH
echo $PATH
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)
[ -z "$PATH_GATK" ] && PATH_GATK=$(dirname $(which GenomeAnalysisTK.jar))
[ -z "$PATH_PICARD" ] && PATH_PICARD=$(dirname $(which FixMateInformation.jar))

echo "[NOTE] set java parameters"
JAVAPARAMS="-Xmx"$(python -c "print int($MEMORY_RECAL*0.8)")"g -Djava.io.tmpdir="$TMP" -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1" 
unset _JAVA_OPTIONS
echo "JAVAPARAMS "$JAVAPARAMS

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--JAVA        --\n" $(java -Xmx200m -version 2>&1)
[ -z "$(which java)" ] && echo "[ERROR] no java detected" && exit 1
echo -e "--samtools    --\n "$(samtools 2>&1 | head -n 3 | tail -n-2)
[ -z "$(which samtools)" ] && echo "[ERROR] no samtools detected" && exit 1
echo -e "--PICARD      --\n "$(java $JAVAPARAMS -jar $PATH_PICARD/FixMateInformation.jar --version 2>&1)
[ ! -f $PATH_PICARD/FixMateInformation.jar ] && echo "[ERROR] no picard detected" && exit 1
echo -e "--R           --\n "$(R --version | head -n 3)
[ -z "$(which R)" ] && echo "[ERROR] no R detected" && exit 1
echo -e "--GATK        --\n "$(java $JAVAPARAMS -jar $PATH_GATK/GenomeAnalysisTK.jar --version)
[ ! -f $PATH_GATK/GenomeAnalysisTK.jar ] && echo "[ERROR] no GATK detected" && exit 1

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "parameters"

# get basename of f
f=${f/%.dummy/} #if input came from pipe
n=${f##*/}
SAMPLE=${n/%$ASD.bam/}

BAMREADS=`head -n1 $f.stats | cut -d " " -f 1`

if [ -n "$SEQREG" ]; then REGION="-L $SEQREG"; fi

#is paired ?SAMPLE=${n/%$ASD.bam/}
p=`grep "paired in " $f.stats | cut -d " " -f 1`
if [ ! "$p" -eq "0" ]; then
    PAIRED="1"
    echo "[NOTE] PAIRED"
else
    echo "[NOTE] SINGLE"
    PAIRED="0"
fi

# delete old bam files unless attempting to recover
if [ -z "$NGSANE_RECOVERFROM" ]; then
    [ -e $OUTDIR/${SAMPLE}$ASR.bam ] && rm $OUTDIR/${SAMPLE}$ASR.bam
    [ -e $OUTDIR/${SAMPLE}$ASR.bam.stats ] && rm $OUTDIR/${SAMPLE}$ASR.bam.stats
fi

if [ -z "$DBSNPVCF" ] || [ ! -e $DBSNPVCF ] ; then
    echo "[ERROR] DBSNPVCF parameter not set or data not found"
    exit 1
fi

# parallelization supported (2013): http://gatkforums.broadinstitute.org/discussion/1975/recommendations-for-parallelizing-gatk-tools
if [[ $(which GenomeAnalysisTK.jar) =~ "2.8" ]]; then 
    echo "[NOTE] new GATK parallele"
    PARALLELENCT="-nct $CPU_RECAL"
	PARALLELENT="-nt $CPU_RECAL"
fi
#
## bwa/name$ASD.bam -> /reCalAln/name$ASD.bam
#f2=${f/$INPUT_REALRECAL/$TASK_RECAL}
## /reCalAln/name$ASD.bam -> /reCalAln/name$ASD.real
#f3=${f2/%bam/real.bam}

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "recall files from tape"

if [ -n "$DMGET" ]; then
	dmget -a $f
	dmget -a $OUTDIR/*
fi
    

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "find intervals to improve"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then
    echo "[NOTE] $CHECKPOINT"

    java $JAVAPARAMS -jar $PATH_GATK/GenomeAnalysisTK.jar -l WARN \
        -T RealignerTargetCreator \
        -I $f \
        -R $FASTA \
        -o $OUTDIR/${SAMPLE}.intervals \
        -known $DBSNPVCF \
        $REGION \
        -nt $CPU_RECAL
        
    # mark checkpoint
    [[ -s $OUTDIR/${SAMPLE}.intervals ]] && NGSANE_CHECKPOINT_CHECK

fi 

################################################################################
NGSANE_CHECKPOINT_INIT "realign"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then
    echo "[NOTE] $CHECKPOINT"

	# parallelization not supported (2013): http://gatkforums.broadinstitute.org/discussion/1975/recommendations-for-parallelizing-gatk-tools
    java $JAVAPARAMS -jar $PATH_GATK/GenomeAnalysisTK.jar -l WARN \
        -T IndelRealigner \
        -I $f \
        -R $FASTA \
        -targetIntervals $OUTDIR/${SAMPLE}.intervals \
        --out $OUTDIR/${SAMPLE}.real.bam \
        -known $DBSNPVCF \
        -compress 0 

    samtools index $OUTDIR/${SAMPLE}.real.bam

    # mark checkpoint
    [[ -s $OUTDIR/${SAMPLE}.real.bam ]] && NGSANE_CHECKPOINT_CHECK

fi 

################################################################################
NGSANE_CHECKPOINT_INIT "count covariantes "

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then
    echo "[NOTE] $CHECKPOINT"
    
    # there seems to be a warning reg. Rscript but it does not affect output
    # BadCigar is necessary in case of "MESSAGE: START (101) > (100) STOP -- this should never happen, please check read" (bwa)
    java $JAVAPARAMS -jar $PATH_GATK/GenomeAnalysisTK.jar -l WARN \
        -T BaseRecalibrator \
        -R $FASTA \
        -knownSites $DBSNPVCF \
        -I $OUTDIR/${SAMPLE}.real.bam \
        -dcov 1000 \
        -cov ReadGroupCovariate \
        -cov QualityScoreCovariate \
        -cov CycleCovariate \
		$PARALLELENCT \
        -rf BadCigar \
        -o $OUTDIR/${SAMPLE}.real.covar.grp
	#   --plot_pdf_file $OUTDIR/GATKorig.pdf \


    # mark checkpoint
    [[ -s $OUTDIR/${SAMPLE}.real.covar.grp ]] && NGSANE_CHECKPOINT_CHECK

fi 


################################################################################
NGSANE_CHECKPOINT_INIT "adjust score"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then
    echo "[NOTE] $CHECKPOINT"

    java $JAVAPARAMS -jar $PATH_GATK/GenomeAnalysisTK.jar -l WARN \
        -T PrintReads \
        -R $FASTA \
        -I $OUTDIR/${SAMPLE}.real.bam \
        -o $OUTDIR/${SAMPLE}.real.recal.bam \
        -BQSR $OUTDIR/${SAMPLE}.real.covar.grp \
        -nct $CPU_RECAL
       
    samtools index $OUTDIR/${SAMPLE}.real.recal.bam

    # mark checkpoint
    [[ -s $OUTDIR/${SAMPLE}.real.recal.bam ]] && NGSANE_CHECKPOINT_CHECK

fi 

################################################################################
NGSANE_CHECKPOINT_INIT "evaluate performace"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then
    echo "[NOTE] $CHECKPOINT"

    java $JAVAPARAMS -jar $PATH_GATK/GenomeAnalysisTK.jar -l WARN \
         -T RecalibrationPerformance \
         -o $OUTDIR/${SAMPLE}.real.recal.performace \
         -R $FASTA \
         --recal $OUTDIR/${SAMPLE}.real.covar.grp \
         -nct $CPU_RECAL 
    
    # mark checkpoint
    [[ -s $OUTDIR/${SAMPLE}.real.recal.performace ]] && NGSANE_CHECKPOINT_CHECK

fi 

################################################################################
NGSANE_CHECKPOINT_INIT "count covariantes after recalibration"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then
    echo "[NOTE] $CHECKPOINT"


	if [[ $(which GenomeAnalysisTK.jar) =~ "2.5" ]]; then 
	echo "[NOTE] old GATK plotting "

    # TODO remove the below once it is not experimental anymore
    echo "********* counting covariantes after recalibration"
    java $JAVAPARAMS -jar $PATH_GATK/GenomeAnalysisTK.jar  -l WARN \
        -T BaseRecalibrator \
        -R $FASTA \
        -knownSites $DBSNPVCF \
        -I $OUTDIR/${SAMPLE}.real.recal.bam  \
        -dcov 1000 \
        -cov ReadGroupCovariate \
        -cov QualityScoreCovariate \
        -cov CycleCovariate \
        --plot_pdf_file $OUTDIR/GATKrecal.pdf \
        -rf BadCigar \
        -o $OUTDIR/${SAMPLE}.real.recal.covar.grp 
    #    -nt $CPU_RECAL

	else
		touch $OUTDIR/${SAMPLE}.real.recal.covar.grp
	fi
    
    #echo " plotting both"
    #java $JAVAPARAMS -jar $PATH_GATK/AnalyzeCovariates.jar \
    #    -recalFile ${f3/%.bam/.covar.csv} \
    #    -outputDir $OUTDIR/GATKorig/$n \
    #    -ignoreQ 5
    
    #    -Rscript $RSCRIPT \
    #    -resources $GATKHOME/R/ \
    
    
    #java $JAVAPARAMS -jar $PATH_GATK/AnalyzeCovariates.jar \
    #    -recalFile ${f3/%.bam/.recal.covar.csv} \
    #    -outputDir $OUTDIR/GATKrcal/$n  \
    #    -ignoreQ 5
    
    # mark checkpoint
    [[ -s $OUTDIR/${SAMPLE}.real.recal.covar.grp ]] && NGSANE_CHECKPOINT_CHECK

fi 

################################################################################
NGSANE_CHECKPOINT_INIT "sort/index"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then
    echo "[NOTE] $CHECKPOINT"

    if [ "$PAIRED" == "1" ]; then
        echo "[NOTE] fixmate"
        RUN_COMMAND="java $JAVAPARAMS -jar $PATH_PICARD/FixMateInformation.jar \
            I=$OUTDIR/${SAMPLE}.real.recal.bam \
            O=$OUTDIR/${SAMPLE}$ASR.bam \
            VALIDATION_STRINGENCY=SILENT \
            COMPRESSION_LEVEL=0 \
            SORT_ORDER=coordinate \
            TMP_DIR=$THISTMP"
            
        echo $RUN_COMMAND && eval $RUN_COMMAND       
    else
        samtools sort -@ $CPU_RECAL $OUTDIR/${SAMPLE}.real.recal.bam $OUTDIR/${SAMPLE}$ASR 
    fi
    
    samtools index $OUTDIR/${SAMPLE}$ASR.bam

    # mark checkpoint
    [[ -s $OUTDIR/${SAMPLE}$ASR.bam ]] && NGSANE_CHECKPOINT_CHECK

fi 

################################################################################
NGSANE_CHECKPOINT_INIT "statistics"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then
    echo "[NOTE] $CHECKPOINT"

    samtools flagstat $OUTDIR/${SAMPLE}$ASR.bam >> $OUTDIR/${SAMPLE}$ASR.bam.stats
    if [ -n "$SEQREG" ]; then
        echo "#custom region " >> $OUTDIR/${SAMPLE}$ASR.bam.stats
        echo $(samtools view -@ $CPU_RECAL -c -F 4 $OUTDIR/${SAMPLE}$ASR.bam $SEQREG )" total reads in region " >> $OUTDIR/${SAMPLE}$ASR.bam.stats
        echo $(samtools view -@ $CPU_RECAL -c -f 3 $OUTDIR/${SAMPLE}$ASR.bam $SEQREG )" properly paired reads in region " >> $OUTDIR/${SAMPLE}$ASR.bam.stats
    fi
   
    BAMREADSRERE=`head -n1 $OUTDIR/${SAMPLE}$ASR.bam.stats | cut -d " " -f 1`
    if [ "$BAMREADSRERE" = "" ]; then let BAMREADSRERE="0"; fi	
    if [[ $BAMREADS -eq $BAMREADSRERE  && ! $BAMREADS -eq 0 ]]; then
        echo "[NOTE] PASS check recalibration and realignment: $BAMREADS == $BAMREADSRERE"
        rm $OUTDIR/${SAMPLE}.real.covar.grp
        rm $OUTDIR/${SAMPLE}.real.recal.covar.grp
        rm $OUTDIR/${SAMPLE}.real.recal.bam
        rm $OUTDIR/${SAMPLE}.real.recal.bam.bai
        rm $OUTDIR/${SAMPLE}.real.recal.bai
        rm $OUTDIR/${SAMPLE}.intervals
        rm $OUTDIR/${SAMPLE}.real.bam.bai
        rm $OUTDIR/${SAMPLE}.real.bam
    else
        echo -e "[ERROR] We are loosing reads during recalibration and realignment-> .bam in $f: bam had $BAMREADS recal+real bam has $BAMREADSRERE"
        exit 1
    fi
    
    # mark checkpoint
    [[ -s $OUTDIR/${SAMPLE}$ASR.bam.stats ]] && NGSANE_CHECKPOINT_CHECK

fi 

################################################################################
[ -e $OUTDIR/${SAMPLE}$ASR.bam.dummy ] && rm $OUTDIR/${SAMPLE}$ASR.bam.dummy
echo ">>>>> recalibration and realignment using GATK - FINISHED"
echo ">>>>> enddate "`date`
