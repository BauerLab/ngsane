#!/bin/bash -e

# Recalibation script; this is for gatk v2.5
# author: Denis C. Bauer
# date: May.2013

# messages to look out for -- relevant for the QC.sh script:
# QCVARIABLES,We are loosing reads,for unmapped read,no such file,file not found,reCalAln.sh: line
# RESULTFILENAME <DIR>/<TASK>/<SAMPLE>.$ASR.bam

echo ">>>>> recalibration and realignment using GATK"
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"


function usage {
echo -e "usage: $(basename $0) -k NGSANE -f bam -r REFERENCE -o OUTDIR -d SNPDB [OPTIONS]

Script running the recalibration and realigment step (GATK)
It expects a bam file (*.$ASD.bam)

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
        -d | --snpdb )          shift; DBROD=$1 ;; # snpdb
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

for MODULE in $MODULE_RECAL; do module load $MODULE; done  # save way to load modules that itself load other modules
export PATH=$PATH_RECAL:$PATH
module list
echo $PATH
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)
PATH_IGVTOOLS=$(dirname $(which igvtools.jar))
PATH_GATK=$(dirname $(which GenomeAnalysisTK.jar))

echo "[NOTE] set java parameters"
JAVAPARAMS="-Xmx"$(python -c "print int($MEMORY_RECAL*0.8)")"g -Djava.io.tmpdir="$TMP" -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1" 
unset _JAVA_OPTIONS
echo "JAVAPARAMS "$JAVAPARAMS

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--JAVA        --\n" $(java -Xmx200m -version 2>&1)
[ -z "$(which java)" ] && echo "[ERROR] no java detected" && exit 1
echo -e "--samtools    --\n "$(samtools 2>&1 | head -n 3 | tail -n-2)
[ -z "$(which samtools)" ] && echo "[ERROR] no samtools detected" && exit 1
echo -e "--R           --\n "$(R --version | head -n 3)
[ -z "$(which R)" ] && echo "[ERROR] no R detected" && exit 1
echo -e "--igvtools    --\n "$(java $JAVAPARAMS -jar $PATH_IGVTOOLS/igvtools.jar version 2>&1)
[ ! -f $PATH_IGVTOOLS/igvtools.jar ] && echo "[ERROR] no igvtools detected" && exit 1
echo -e "--GATK        --\n "$(java $JAVAPARAMS -jar $PATH_GATK/GenomeAnalysisTK.jar --version)
[ ! -f $PATH_GATK/GenomeAnalysisTK.jar ] && echo "[ERROR] no GATK detected" && exit 1

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="parameters"

# get basename of f
f=${f/%.dummy/} #if input came from pipe
n=${f##*/}

BAMREADS=`head -n1 $f.stats | cut -d " " -f 1`

if [ -n "$SEQREG" ]; then REGION="-L $SEQREG"; fi

#is paired ?
p=`grep "paired in " $f.stats | cut -d " " -f 1`
if [ ! "$p" -eq "0" ]; then
    PAIRED="1"
    echo "[NOTE] PAIRED"
else
    echo "[NOTE] SINGLE"
    PAIRED="0"
fi

# delete old bam files unless attempting to recover
if [ -z "$RECOVERFROM" ]; then
    [ -e $OUTDIR/${n/%$ASD.bam/$ASR.bam} ] && rm $OUTDIR/${n/%$ASD.bam/$ASR.bam}
    [ -e $OUTDIR/${n/%$ASD.bam/$ASR.bam}.stats ] && rm $OUTDIR/${n/%$ASD.bam/$ASR.bam}.stats
fi

if [ -z "$DBROD" ] || [ ! -e $DBROD ] ; then
    echo "[ERROR] DBROD parameter not set or data not found"
    exit 1
fi

# bwa/name.$ASD.bam -> /reCalAln/name.$ASD.bam
f2=${f/$TASKBWA/$TASKRCA}
# /reCalAln/name.$ASD.bam -> /reCalAln/name.$ASD.real
f3=${f2/bam/real.bam}

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="recall files from tape"

if [ -n "$DMGET" ]; then
	dmget -a $f
fi
    
echo -e "\n********* $CHECKPOINT\n"    
################################################################################
CHECKPOINT="find intervals to improve"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
    echo "[NOTE] $CHECKPOINT"

    java $JAVAPARAMS -jar $PATH_GATK/GenomeAnalysisTK.jar -l WARN \
        -T RealignerTargetCreator \
        -I $f \
        -R $FASTA \
        -o $f2.intervals \
        -known $DBROD \
        $REGION \
        -nt $CPU_RECAL
        
    # mark checkpoint
    if [ -f $f2.intervals ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi 

################################################################################
CHECKPOINT="realign"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
    echo "[NOTE] $CHECKPOINT"

    java $JAVAPARAMS -jar $PATH_GATK/GenomeAnalysisTK.jar -l WARN \
        -T IndelRealigner \
        -I $f \
        -R $FASTA \
        -targetIntervals $f2.intervals \
        --out ${f2/bam/real.bam} \
        -known $DBROD \
        -compress 0 
    #    -nt $CPU_RECAL

    #samtools sort ${f2/bam/real.fix.bam} $OUTDIR/${n/$ASD.bam/$ASR}
    samtools index $f3

    # mark checkpoint
    if [ -f ${f2/bam/real.bam} ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi 

################################################################################
CHECKPOINT="count covariantes"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
    echo "[NOTE] $CHECKPOINT"
    
    # there seems to be a warning reg. Rscript but it does not affect output
    java $JAVAPARAMS -jar $PATH_GATK/GenomeAnalysisTK.jar -l WARN \
        -T BaseRecalibrator \
        -R $FASTA \
        -knownSites $DBROD \
        -I $f3 \
        -dcov 1000 \
        -cov ReadGroupCovariate \
        -cov QualityScoreCovariate \
        -cov CycleCovariate \
        --plot_pdf_file $OUTDIR/GATKorig.pdf \
        -o ${f3/.bam/.covar.grp} 
    #    -nt $CPU_RECAL

    # mark checkpoint
    if [ -f ${f3/.bam/.covar.grp} ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi 

################################################################################
CHECKPOINT="adjust score"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
    echo "[NOTE] $CHECKPOINT"

    java $JAVAPARAMS -jar $PATH_GATK/GenomeAnalysisTK.jar -l WARN \
        -T PrintReads \
        -R $FASTA \
        -I $f3 \
        -o ${f3/.bam/.recal.bam} \
        -BQSR ${f3/.bam/.covar.grp} \
        -nct $CPU_RECAL
        
    samtools index ${f3/.bam/.recal.bam}

    # mark checkpoint
    if [ -f ${f3/.bam/.recal.bam} ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi 

################################################################################
CHECKPOINT="evaluate performace"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
    echo "[NOTE] $CHECKPOINT"

    java $JAVAPARAMS -jar $PATH_GATK/GenomeAnalysisTK.jar -l WARN \
         -T RecalibrationPerformance \
         -o ${f3/.bam/.recal.performace} \
         -R $FASTA \
         --recal ${f3/.bam/.covar.grp} \
         -nct $CPU_RECAL 
    
    # mark checkpoint
    if [ -f ${f3/.bam/.recal.performace} ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi 

################################################################################
CHECKPOINT="count covariantes after recalibration"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
    echo "[NOTE] $CHECKPOINT"

    # TODO remove the below once it is not experimental anymore
    echo "********* counting covariantes after recalibration"
    java $JAVAPARAMS -jar $PATH_GATK/GenomeAnalysisTK.jar  -l WARN \
        -T BaseRecalibrator \
        -R $FASTA \
        -knownSites $DBROD \
        -I ${f3/.bam/.recal.bam}  \
        -dcov 1000 \
        -cov ReadGroupCovariate \
        -cov QualityScoreCovariate \
        -cov CycleCovariate \
        --plot_pdf_file $OUTDIR/GATKrecal.pdf \
        -o ${f3/.bam/.recal.covar.grp} 
    #    -nt $CPU_RECAL
    
    #echo " plotting both"
    #java $JAVAPARAMS -jar $PATH_GATK/AnalyzeCovariates.jar \
    #    -recalFile ${f3/.bam/.covar.csv} \
    #    -outputDir $OUTDIR/GATKorig/$n \
    #    -ignoreQ 5
    
    #    -Rscript $RSCRIPT \
    #    -resources $GATKHOME/R/ \
    
    
    #java $JAVAPARAMS -jar $PATH_GATK/AnalyzeCovariates.jar \
    #    -recalFile ${f3/.bam/.recal.covar.csv} \
    #    -outputDir $OUTDIR/GATKrcal/$n  \
    #    -ignoreQ 5
    
    # mark checkpoint
    if [ -f ${f3/.bam/.recal.covar.grp} ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi 

################################################################################
CHECKPOINT="sort/index"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
    echo "[NOTE] $CHECKPOINT"

    if [ "$PAIRED" == "1" ]; then
        # fix mates
        samtools sort -@ $CPU_RECAL -n ${f3/bam/recal.bam} ${f3/bam/recal.tmp}
        samtools fixmate ${f3/bam/recal.tmp.bam} ${f3/bam/recal.bam}
        [ -e ${f3/bam/recal.tmp.bam} ] && rm ${f3/bam/recal.tmp.bam}
    fi
    
    samtools sort -@ $CPU_RECAL ${f3/bam/recal.bam} $OUTDIR/${n/%$ASD.bam/$ASR}
    samtools index $OUTDIR/${n/%$ASD.bam/$ASR.bam}

    # mark checkpoint
    if [ -f $OUTDIR/${n/%$ASD.bam/$ASR.bam} ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi 

################################################################################
CHECKPOINT="statistics"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
    echo "[NOTE] $CHECKPOINT"

    samtools flagstat $OUTDIR/${n/%$ASD.bam/$ASR.bam} >> $OUTDIR/${n/%$ASD.bam/$ASR.bam}.stats
    if [ -n "$SEQREG" ]; then
        echo "#custom region " >> $OUTDIR/${n/%$ASD.bam/$ASR.bam}.stats
        echo $(samtools view -@ $CPU_RECAL -c -F 4 $OUTDIR/${n/%$ASD.bam/$ASR.bam} $SEQREG )" total reads in region " >> $OUTDIR/${n/%$ASD.bam/$ASR.bam}.stats
        echo $(samtools view -@ $CPU_RECAL -c -f 3 $OUTDIR/${n/%$ASD.bam/$ASR.bam} $SEQREG )" properly paired reads in region " >> $OUTDIR/${n/%$ASD.bam/$ASR.bam}.stats
    fi

    #f2=/reCalAln/name.$ASD.bam
    #f3=/reCalAln/name.$ASD.real.bam
    
    BAMREADSRERE=`head -n1 $OUTDIR/${n/%$ASD.bam/$ASR.bam}.stats | cut -d " " -f 1`
    if [ "$BAMREADSRERE" = "" ]; then let BAMREADSRERE="0"; fi	
    if [[ $BAMREADS -eq $BAMREADSRERE  && ! $BAMREADS -eq 0 ]]; then
        echo "[NOTE] PASS check recalibration and realignment: $BAMREADS == $BAMREADSRERE"
        rm ${f3/.bam/.covar.grp}
        rm ${f3/.bam/.recal.covar.grp}
        rm ${f3/.bam/.recal.bam}
        rm ${f3/.bam/.recal.bam}.bai
        rm ${f3/.bam/.recal}.bai
        rm $f2.intervals
        rm $f3
        rm ${f2/bam/real.bam}.bai
        rm ${f2/bam/real}.bai
    else
        echo -e "[ERROR] We are loosing reads during recalibration and realignment-> .bam in $f: bam had $BAMREADS recal+real bam has $BAMREADSRERE"
        exit 1
    fi
    
    # mark checkpoint
    if [ -f $OUTDIR/${n/%$ASD.bam/$ASR.bam}.stats ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi 

################################################################################
CHECKPOINT="coverage track"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
    echo "[NOTE] $CHECKPOINT"

    # get the coverage track
    java $JAVAPARAMS -jar $PATH_IGVTOOLS/igvtools.jar count $OUTDIR/${n/%$ASD.bam/$ASR.bam} \
        $OUTDIR/${n/%$ASD.bam/$ASR.bam}.cov.tdf ${FASTA/fasta/genome}
    
    # mark checkpoint
    if [ -f $OUTDIR/${n/%$ASD.bam/$ASR.bam}.cov.tdf ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi

################################################################################
[ -e $OUTDIR/${n/%$ASD.bam/$ASR.bam}.dummy ] && rm $OUTDIR/${n/%$ASD.bam/$ASR.bam}.dummy
echo ">>>>> recalibration and realignment using GATK - FINISHED"
echo ">>>>> enddate "`date`
