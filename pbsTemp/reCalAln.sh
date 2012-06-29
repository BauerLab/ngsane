#!/bin/bash

# Recalibation script, expects .asd
# author: Denis C. Bauer
# date: Nov.2010

# messages to look out for -- relevant for the QC.sh script:
# QCVARIABLES,We are loosing reads,for unmapped read,no such file,file not found,reCalAln.sh: line


echo ">>>>> recalibration and realignment using GATK"
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> reCalAln.sh $*"


function usage {
echo -e "usage: $(basename $0) -k HISEQINF -f bam -r REFERENCE -o OUTDIR -d SNPDB [OPTIONS]

Script running the recalibration and realigment step (GATK)
It expects a bam file (*.asd.bam)

required:
  -k | --toolkit <path>     location of the HiSeqInf repository 
  -f | --fastq <file>       bam file
  -r | --reference <file>   reference genome
  -o | --outdir <path>      output dir
  -d | --snpdb <path>       path to a SNPdb instance (rod)

options:
  -t | --threads <nr>       number of CPUs to use (default: 1)
  -R | --region <ps>        region of specific interest, e.g. targeted reseq
                             format chr:pos-pos
"
exit
}

if [ ! $# -gt 4 ]; then usage ; fi

#DEFAULTS
THREADS=1
#JAVAPARAMS="-Xmx"$MEMORY"g -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1 -XX:MaxDirectMemorySize=4G"
module load R
module load jdk/1.7.0_03

#INPUTS
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; HISEQINF=$1 ;; # location of the HiSeqInf repository
        -t | --threads )        shift; THREADS=$1 ;; # number of CPUs to use
        -f | --fastq )          shift; f=$1 ;; # bam file
        -r | --reference )      shift; FASTA=$1 ;; # reference genome
	-d | --snpdb )          shift; DBROD=$1 ;; # snpdb
        -o | --outdir )         shift; OUT=$1 ;; # output dir
	-R | --region )         shift; SEQREG=$1 ;; # (optional) region of specific interest, e.g. targeted reseq
        -h | --help )           usage ;;
        * )                     usage
    esac
    shift
done


#PROGRAMS
. $HISEQINF/pbsTemp/header.sh
export PATH=$PATH:$RSCRIPT


n=`basename $f`

BAMREADS=`head -n1 $f.stats | cut -d " " -f 1`

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
if [ -e $OUT/${n/asd/asdrr} ]; then rm $OUT/${n/asd/asdrr}; fi
if [ -e $OUT/${n/asd/asdrr}.stats ]; then rm $OUT/${n/asd/asdrr}.stats; fi

if [ ! -d $OUT/GATKorig ]; then
    mkdir $OUT/GATKorig	
    mkdir $OUT/GATKrcal
fi

if [ ! -d $OUT/GATKorig/$n ]; then
    mkdir $OUT/GATKorig/$n
    mkdir $OUT/GATKrcal/$n
else 
    rm -r $OUT/GATKorig/$n
    rm -r $OUT/GATKrcal/$n
fi


# bwa/name.asd.bam -> /reCalAln/name.asd.bam
f2=${f/$TASKBWA/$TASKRCA}

#################
# REALIGMENT
#################

echo "********* realignment"
echo "********* find intervals to improve"
java -Xmx16g -Djava.io.tmpdir=$TMP -jar $GATKJAR/GenomeAnalysisTK.jar -l WARN \
    -T RealignerTargetCreator \
    -I $f \
    -R $FASTA \
    -o $f2.intervals \
    -known $DBROD \
    -nt $THREADS

echo "********* realine them"
java -Xmx16g -Djava.io.tmpdir=$TMP -jar $GATKJAR/GenomeAnalysisTK.jar -l WARN \
    -T IndelRealigner \
    -I $f \
    -R $FASTA \
    -targetIntervals $f2.intervals \
    --out ${f2/bam/real.bam} \
    -known $DBROD \
    -compress 0 
#    -nt $THREADS

# /reCalAln/name.asd.bam /reCalAln/name.asd.real.bam
f3=${f2/bam/real.bam}


echo "********* index"
#$SAMTOOLS sort ${f2/bam/real.fix.bam} $OUT/${n/asd.bam/asdrr}
$SAMTOOLS index $f3


#################
# RECALIBRATION
#################

echo "********* recalibrating"
echo "********* counting covariantes" 
java -Xmx16g -Djava.io.tmpdir=$TMP -jar $GATKJAR/GenomeAnalysisTK.jar -l WARN \
    -T CountCovariates \
    -R $FASTA \
    -knownSites $DBROD \
    -I $f3 \
    -dcov 1000 \
    -cov ReadGroupCovariate \
    -cov QualityScoreCovariate \
    -cov CycleCovariate \
    -cov DinucCovariate \
    -recalFile ${f3/.bam/.covar.csv} \
    -nt $THREADS


echo "********* change score"
java -Xmx16g -Djava.io.tmpdir=$TMP -jar $GATKJAR/GenomeAnalysisTK.jar -l WARN \
    -T TableRecalibration \
    -R $FASTA \
    -I $f3 \
    --out ${f3/.bam/.recal.bam} \
    -recalFile ${f3/.bam/.covar.csv} 
#    -nt $THREADS



echo "********* QC Step"
echo "********* index"
$SAMTOOLS index ${f3/.bam/.recal.bam}

echo "********* counting covariantes after recalibration"
java -Xmx16g -jar $GATKJAR/GenomeAnalysisTK.jar -l WARN \
    -T CountCovariates \
    -R $FASTA \
    -knownSites $DBROD \
    -I ${f3/.bam/.recal.bam} \
    -dcov 1000 \
    -cov ReadGroupCovariate \
    -cov QualityScoreCovariate \
    -cov CycleCovariate \
    -cov DinucCovariate \
    -recalFile ${f3/.bam/.recal.covar.csv} \
    -nt $THREADS
#    --default_platform illumina

echo "********* plotting both"
java -Xmx4g -jar $GATKJAR/AnalyzeCovariates.jar \
    -recalFile ${f3/.bam/.covar.csv} \
    -outputDir $OUT/GATKorig/$n \
    -ignoreQ 5

#    -Rscript $RSCRIPT \
#    -resources $GATKHOME/R/ \


java -Xmx4g -jar $GATKJAR/AnalyzeCovariates.jar \
    -recalFile ${f3/.bam/.recal.covar.csv} \
    -outputDir $OUT/GATKrcal/$n  \
    -ignoreQ 5



echo "********* sort/index"
$SAMTOOLS sort ${f3/bam/recal.bam} $OUT/${n/$ASD.bam/$ASR}
$SAMTOOLS index $OUT/${n/$ASD/$ASR}

# statistics
echo "********* statistics"
$SAMTOOLS flagstat $OUT/${n/$ASD/$ASR} >> $OUT/${n/$ASD/$ASR}.stats
if [ -n $SEQREG ]; then
    echo "#custom region " >> $OUT/${n/$ASD/$ASR}.stats
    echo `$SAMTOOLS view $OUT/${n/$ASD/$ASR} $SEQREG | wc -l`" total reads in region " >> $OUT/${n/$ASD/$ASR}.stats
    echo `$SAMTOOLS view -f 2 $OUT/${n/$ASD/$ASR} $SEQREG | wc -l`" properly paired reads in region " >> $OUT/${n/$ASD/$ASR}.stats
fi

#f2=/reCalAln/name.asd.bam
#f3=/reCalAln/name.asd.real.bam

BAMREADSRERE=`head -n1 $OUT/${n/$ASD/$ASR}.stats | cut -d " " -f 1`
if [ "$BAMREADSRERE" = "" ]; then let BAMREADSRERE="0"; fi	
if [ $BAMREADS -eq $BAMREADSRERE ]; then
    echo "-----------------> PASS check recalibration and realignment: $BAMREADS == $BAMREADSRERE"
    rm ${f3/.bam/.covar.csv}
    rm ${f3/.bam/.recal.covar.csv}
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
java -Xmx1g -jar $IGVTOOLS count $OUT/${n/$ASD/$ASR} \
    $OUT/${n/$ASD/$ASR}.cov.tdf ${FASTA/fasta/genome}

echo ">>>>> recalibration and realignment using GATK - FINISHED"
echo ">>>>> enddate "`date`