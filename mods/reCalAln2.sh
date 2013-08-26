#!/bin/bash

# Recalibation script; this is for gatk v2.5
# author: Denis C. Bauer
# date: May.2013

# messages to look out for -- relevant for the QC.sh script:
# QCVARIABLES,We are loosing reads,for unmapped read,no such file,file not found,reCalAln.sh: line


echo ">>>>> recalibration and realignment using GATK"
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"


function usage {
echo -e "usage: $(basename $0) -k NGSANE -f bam -r REFERENCE -o OUTDIR -d SNPDB [OPTIONS]

Script running the recalibration and realigment step (GATK)
It expects a bam file (*.asd.bam)

required:
  -k | --toolkit <path>     location of the NGSANE repository 
  -f | --fastq <file>       bam file
  -r | --reference <file>   reference genome
  -o | --outdir <path>      output dir
  -d | --snpdb <path>       path to a SNPdb instance (rod)

options:
  -t | --threads <nr>       number of CPUs to use (default: 1)
  -L | --region <ps>        region of specific interest, e.g. targeted reseq
                             format chr:pos-pos
"
exit
}

if [ ! $# -gt 4 ]; then usage ; fi

#DEFAULTS
#THREADS=1
#JAVAPARAMS="-Xmx"$MEMORY"g -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1 -XX:MaxDirectMemorySize=4G"
#module load R
#module load jdk #/1.7.0_03

#INPUTS
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the NGSANE repository
        -t | --threads )        shift; MYTHREADS=$1 ;; # number of CPUs to use
        -f | --fastq )          shift; f=$1 ;; # bam file
        -r | --reference )      shift; FASTA=$1 ;; # reference genome
        -d | --snpdb )          shift; DBROD=$1 ;; # snpdb
        -o | --outdir )         shift; MYOUT=$1 ;; # output dir
        -L | --region )         shift; SEQREG=$1 ;; # (optional) region of specific interest, e.g. targeted reseq
        -h | --help )           usage ;;
        * )                     usage
    esac
    shift
done


#PROGRAMS
. $CONFIG
. ${NGSANE_BASE}/conf/header.sh
. $CONFIG

echo "********** programs"
for MODULE in $MODULE_GATK; do module load $MODULE; done  # save way to load modules that itself load other modules
export PATH=$PATH_GATK:$PATH
module list
echo $PATH
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)
PATH_IGVTOOLS=$(dirname $(which igvtools.jar))
PATH_GATK=$(dirname $(which GenomeAnalysisTK.jar))
echo -e "--JAVA    --\n" $(java -Xmx200m -version 2>&1)
echo -e "--samtools--\n "$(samtools 2>&1 | head -n 3 | tail -n-2)
echo -e "--R       --\n "$(R --version | head -n 3)
echo -e "--igvtools--\n "$(java -jar -Xmx200m $PATH_IGVTOOLS/igvtools.jar version 2>&1)
echo -e "--GATK    --\n "$(java -jar -Xmx200m $PATH_GATK/GenomeAnalysisTK.jar --version)


echo "[NOTE] set java parameters"
JAVAPARAMS="-Xmx"$(python -c "print int($MEMORY_RECAL*0.8)")"g -Djava.io.tmpdir="$TMP"  -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1" 
unset _JAVA_OPTIONS
echo "JAVAPARAMS "$JAVAPARAMS

# get basename of f
n=${f##*/}

BAMREADS=`head -n1 $f.stats | cut -d " " -f 1`

if [ -n "$SEQREG" ]; then REGION="-L $SEQREG"; fi

#is paired ?
p=`grep "paired in " $f.stats | cut -d " " -f 1`
if [ ! "$p" -eq "0" ]; then
    PAIRED="1"
    echo "********* PAIRED"
else
    echo "********* SINGLE"
    PAIRED="0"
fi

# delete old output files
if [ -e $MYOUT/${n/asd/asdrr} ]; then rm $MYOUT/${n/asd/asdrr}; fi
if [ -e $MYOUT/${n/asd/asdrr}.stats ]; then rm $MYOUT/${n/asd/asdrr}.stats; fi


# bwa/name.asd.bam -> /reCalAln/name.asd.bam
f2=${f/$TASKBWA/$TASKRCA}
# /reCalAln/name.asd.bam -> /reCalAln/name.asd.real
f3=${f2/bam/real.bam}

#################
# REALIGMENT
#################


echo "[NOTE] realignment"
echo "********* find intervals to improve"
java $JAVAPARAMS -jar $PATH_GATK/GenomeAnalysisTK.jar -l WARN \
    -T RealignerTargetCreator \
    -I $f \
    -R $FASTA \
    -o $f2.intervals \
    -known $DBROD \
    $REGION \
    -nt $MYTHREADS

echo "********* realine them"
java $JAVAPARAMS -jar $PATH_GATK/GenomeAnalysisTK.jar -l WARN \
    -T IndelRealigner \
    -I $f \
    -R $FASTA \
    -targetIntervals $f2.intervals \
    --out ${f2/bam/real.bam} \
    -known $DBROD \
    -compress 0 
#    -nt $MYTHREADS


echo "********* index after realigment"
#samtools sort ${f2/bam/real.fix.bam} $MYOUT/${n/asd.bam/asdrr}
samtools index $f3

#################
# RECALIBRATION
#################

echo "********* recalibrating"
echo "********* counting covariantes" 
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
    --plot_pdf_file $MYOUT/GATKorig.pdf \
    -o ${f3/.bam/.covar.grp} 
#    -nt $MYTHREADS


echo "********* change score"
java $JAVAPARAMS -jar $PATH_GATK/GenomeAnalysisTK.jar -l WARN \
    -T PrintReads \
    -R $FASTA \
    -I $f3 \
    -o ${f3/.bam/.recal.bam} \
    -BQSR ${f3/.bam/.covar.grp} \
    -nct $MYTHREADS

echo "********* QC Step"
samtools index ${f3/.bam/.recal.bam}

echo "********* evaluate performace" #experimental
java $JAVAPARAMS -jar $PATH_GATK/GenomeAnalysisTK.jar -l WARN \
     -T RecalibrationPerformance \
     -o ${f3/.bam/.recal.performace} \
     -R $FASTA \
     --recal ${f3/.bam/.covar.grp} \
     -nct $MYTHREADS 


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
    --plot_pdf_file $MYOUT/GATKrecal.pdf \
    -o ${f3/.bam/.recal.covar.grp} 
#    -nt $MYTHREADS

#echo " plotting both"
#java $JAVAPARAMS -jar $PATH_GATK/AnalyzeCovariates.jar \
#    -recalFile ${f3/.bam/.covar.csv} \
#    -outputDir $MYOUT/GATKorig/$n \
#    -ignoreQ 5

#    -Rscript $RSCRIPT \
#    -resources $GATKHOME/R/ \


#java $JAVAPARAMS -jar $PATH_GATK/AnalyzeCovariates.jar \
#    -recalFile ${f3/.bam/.recal.covar.csv} \
#    -outputDir $MYOUT/GATKrcal/$n  \
#    -ignoreQ 5


echo "********* sort/index"
samtools sort ${f3/bam/recal.bam} $MYOUT/${n/$ASD.bam/$ASR}

if [ "$PAIRED" == "1" ]; then
    # fix mates
    samtools fixmate $MYOUT/${n/$ASD/$ASR} $MYOUT/${n/$ASD.bam/$ASR}.tmp.bam
    mv $MYOUT/${n/$ASD.bam/$ASR}.tmp.bam $MYOUT/${n/$ASD/$ASR}
fi

samtools index $MYOUT/${n/$ASD/$ASR}


echo "********* statistics"
samtools flagstat $MYOUT/${n/$ASD/$ASR} >> $MYOUT/${n/$ASD/$ASR}.stats
if [ -n $SEQREG ]; then
    echo "#custom region " >> $MYOUT/${n/$ASD/$ASR}.stats
    echo `samtools view $MYOUT/${n/$ASD/$ASR} $SEQREG | wc -l`" total reads in region " >> $MYOUT/${n/$ASD/$ASR}.stats
    echo `samtools view -f 2 $MYOUT/${n/$ASD/$ASR} $SEQREG | wc -l`" properly paired reads in region " >> $MYOUT/${n/$ASD/$ASR}.stats
fi

#f2=/reCalAln/name.asd.bam
#f3=/reCalAln/name.asd.real.bam

BAMREADSRERE=`head -n1 $MYOUT/${n/$ASD/$ASR}.stats | cut -d " " -f 1`
if [ "$BAMREADSRERE" = "" ]; then let BAMREADSRERE="0"; fi	
if [[ $BAMREADS -eq $BAMREADSRERE  && ! $BAMREADS -eq 0 ]]; then
    echo "-----------------> PASS check recalibration and realignment: $BAMREADS == $BAMREADSRERE"
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
    echo -e "***ERROR**** We are loosing reads during recalibration and realignment-> .bam in $f: bam had $BAMREADS recal+real bam has $BAMREADSRERE"
    exit 1
fi

# get the coverage track
echo "********* coverage track"
java $JAVAPARAMS -jar $PATH_IGVTOOLS/igvtools.jar count $MYOUT/${n/$ASD/$ASR} \
    $MYOUT/${n/$ASD/$ASR}.cov.tdf ${FASTA/fasta/genome}

echo ">>>>> recalibration and realignment using GATK - FINISHED"
echo ">>>>> enddate "`date`
