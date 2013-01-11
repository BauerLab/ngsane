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
  -L | --region <ps>        region of specific interest, e.g. targeted reseq
                             format chr:pos-pos
"
exit
}

if [ ! $# -gt 4 ]; then usage ; fi

#DEFAULTS
THREADS=1
#JAVAPARAMS="-Xmx"$MEMORY"g -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1 -XX:MaxDirectMemorySize=4G"
module load R
module load jdk #/1.7.0_03

#INPUTS
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the HiSeqInf repository
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
. $HISEQINF/pbsTemp/header.sh
. $CONFIG
export PATH=$PATH:$RSCRIPT


n=`basename $f`

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

if [ ! -d $MYOUT/GATKorig ]; then
    mkdir $MYOUT/GATKorig	
    mkdir $MYOUT/GATKrcal
fi

if [ ! -d $MYOUT/GATKorig/$n ]; then
    mkdir $MYOUT/GATKorig/$n
    mkdir $MYOUT/GATKrcal/$n
else 
    rm -r $MYOUT/GATKorig/$n
    rm -r $MYOUT/GATKrcal/$n
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
    $REGION \
    -nt $MYTHREADS

echo "********* realine them"
java -Xmx16g -Djava.io.tmpdir=$TMP -jar $GATKJAR/GenomeAnalysisTK.jar -l WARN \
    -T IndelRealigner \
    -I $f \
    -R $FASTA \
    -targetIntervals $f2.intervals \
    --out ${f2/bam/real.bam} \
    -known $DBROD \
    -compress 0 
#    -nt $MYTHREADS

# /reCalAln/name.asd.bam /reCalAln/name.asd.real.bam
f3=${f2/bam/real.bam}


echo "********* index"
#$SAMTOOLS sort ${f2/bam/real.fix.bam} $MYOUT/${n/asd.bam/asdrr}
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
    -nt $MYTHREADS


echo "********* change score"
java -Xmx16g -Djava.io.tmpdir=$TMP -jar $GATKJAR/GenomeAnalysisTK.jar -l WARN \
    -T TableRecalibration \
    -R $FASTA \
    -I $f3 \
    --out ${f3/.bam/.recal.bam} \
    -recalFile ${f3/.bam/.covar.csv} 
#    -nt $MYTHREADS



echo "********* QC Step"
echo "********* index"
$SAMTOOLS index ${f3/.bam/.recal.bam}

echo "********* counting covariantes after recalibration"
java -Xmx16g -Djava.io.tmpdir=$TMP -jar $GATKJAR/GenomeAnalysisTK.jar -l WARN \
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
    -nt $MYTHREADS
#    --default_platform illumina

echo "********* plotting both"
java -Xmx4g -jar $GATKJAR/AnalyzeCovariates.jar \
    -recalFile ${f3/.bam/.covar.csv} \
    -outputDir $MYOUT/GATKorig/$n \
    -ignoreQ 5

#    -Rscript $RSCRIPT \
#    -resources $GATKHOME/R/ \


java -Xmx4g -jar $GATKJAR/AnalyzeCovariates.jar \
    -recalFile ${f3/.bam/.recal.covar.csv} \
    -outputDir $MYOUT/GATKrcal/$n  \
    -ignoreQ 5



echo "********* sort/index"
$SAMTOOLS sort ${f3/bam/recal.bam} $MYOUT/${n/$ASD.bam/$ASR}
$SAMTOOLS index $MYOUT/${n/$ASD/$ASR}

# statistics
echo "********* statistics"
$SAMTOOLS flagstat $MYOUT/${n/$ASD/$ASR} >> $MYOUT/${n/$ASD/$ASR}.stats
if [ -n $SEQREG ]; then
    echo "#custom region " >> $MYOUT/${n/$ASD/$ASR}.stats
    echo `$SAMTOOLS view $MYOUT/${n/$ASD/$ASR} $SEQREG | wc -l`" total reads in region " >> $MYOUT/${n/$ASD/$ASR}.stats
    echo `$SAMTOOLS view -f 2 $MYOUT/${n/$ASD/$ASR} $SEQREG | wc -l`" properly paired reads in region " >> $MYOUT/${n/$ASD/$ASR}.stats
fi

#f2=/reCalAln/name.asd.bam
#f3=/reCalAln/name.asd.real.bam

BAMREADSRERE=`head -n1 $MYOUT/${n/$ASD/$ASR}.stats | cut -d " " -f 1`
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
java -Xmx1g -jar $IGVTOOLS count $MYOUT/${n/$ASD/$ASR} \
    $MYOUT/${n/$ASD/$ASR}.cov.tdf ${FASTA/fasta/genome}

echo ">>>>> recalibration and realignment using GATK - FINISHED"
echo ">>>>> enddate "`date`
