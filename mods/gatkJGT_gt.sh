#!/bin/bash -e

# Script for running joint genotyping using GATK3.5/6
# 
# author: Tim Kahlke
# date: July 2016
# RESULTFILENAME <DIR>/<TASK>/recalibrated_variants.vcf

echo ">>>>> Joint genotyping of gVCF files with subsequent VQSR"
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"

function usage {
echo -e "usage: $(basename $0) 

required:
    -k | --toolkit  location of NGSANE 
    -i | --input <path> directory of input gVCF files
    -o | --output <path>    directory for output files
    -r | --reference <file> path to reference genome

optional:
    -b | --inbreeding     
    "
exit
}

if [ ! $# -gt 3 ]; then usage ; fi

#DEFAULTS
INBRED=""   # Inbreed coefficient set for model training
MAXGAUS="4" # Max gaussian for model

#INPUTS                                                                                                           
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;;     # location of the NGSANE repository                      
        -f | --file )           shift; INPUTFILE=$1;;       # gvcf file
        -o | --outdir )         shift; OUTDIR=$1 ;;     # output dir                                                     
        -r | --reference )      shift; REFERENCE=$1;;   # Reference human genome
        -b | --inbreed )        shift; INBRED=$1 ;; # Inbred coefficient for model training
        --maxGaussians )        shift; MAXGAUS=$1 ;; # MAXGAUS
        --recover-from )        shift; NGSANE_RECOVERFROM=$1 ;; # attempt to recover from log file
        -h | --help )           usage ;;
        * )                     echo "don't understand "$1
    esac
    shift
done

#Default
ALL=0
TS_FILTER=99.9
ANNO=QD,FS,MQRankSum,ReadPosRankSum
EXCLUDE=NC_007605,hs37d5 

#PROGRAMS
. $CONFIG
. ${NGSANE_BASE}/conf/header.sh
. $CONFIG

################################################################################
NGSANE_CHECKPOINT_INIT "programs"

# save way to load modules that itself load other modules
hash module 2>/dev/null && for MODULE in $MODULE_GATKJGT; do module load $MODULE; done && module list

export PATH=$PATH_GATKJGT:$PATH
echo "PATH=$PATH"

PATH_GATK_JAR=$(which GenomeAnalysisTK.jar)
echo "PATH=$PATH"
echo -e "--JAVA        --\n" $(java -Xmx200m -version 2>&1)
[ -z "$(which java)" ] && echo "[ERROR] no java detected" && exit 1
[ ! -f $PATH_GATK_JAR ] && echo "[ERROR] no gatk detected" && exit 1

JAVAPARAMS="-Xmx"$(python -c "print int($MEMORY_GATKJGT_GT*0.75)")"g -Djava.io.tmpdir="$TMP"  -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1"
unset _JAVA_OPTIONS
echo "JAVAPARAMS "$JAVAPARAMS

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--R           --\n "$(R --version | head -n 3)
[ -z "$(which R)" ] && echo "[ERROR] no R detected" && exit 1

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "parameters"


if [ ! -d $OUTDIR ]; then mkdir -p $OUTDIR; fi

# delete old variant files unless attempting to recover
if [ -z "$NGSANE_RECOVERFROM" ]; then
    [ -e $OUTDIR/raw_snps_raw_indels.vcf ] && rm $OUTDIR/raw_snps_raw_indels.vcf
fi

if [ -z "$DBSNPVCF" ] || [ ! -e $DBSNPVCF ]; then
    echo "[ERRPR] DBSNPVCF parameter not specified or file not found"
    exit 1
fi
if [ -z "$HAPMAPVCF" ] || [ ! -e $HAPMAPVCF ]; then
    echo "[ERRPR] HAPMAPVCF parameter not specified or file not found"
    exit 1
fi
if [ -z "$ONEKSNPS" ] || [ ! -e $ONEKSNPS ]; then
    echo "[ERRPR] ONEKGVCF parameter not specified or file not found"
    exit 1
fi
if [ -z "$ONEKOMNI" ] || [ ! -e $ONEKOMNI ]; then
    echo "[ERRPR] ONEKOMNI parameter not specified or file not found"
    exit 1
fi
if [ -z "$REFERENCE" ] || [ ! -e $REFERENCE ]; then
    echo "[ERRPR] REFERENCE parameter not specified or file not found"
    exit 1
fi

if [ -z "$MILLS" ]; then
    echo "[ERRPR] No mills resource file specified"
    exit 1
fi

PARAM_STRING=""
if [ -n $ANNO]; then
    IFS="," read -ra ANNOS <<<$ANNO
    for a in $ANNOS; do
        PARAM_STRING="-an $a "
    done
fi

if [-n $EXCLUDE ]; then
    IFS="," read -ra EXC <<<$EXCLUDE
    for x in $EXC; do
        PARAM_STRING="-XL $x "
    done
fi

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "recall files from tape"

if [ -n "$DMGET" ]; then
	dmget -a $INPUTFILE
    dmget -a $OUTDIR/*
fi
    
NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "joint genotyping"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    echo "Joint genotyping"

    GVCFS=""
    if [ -z $ALL ];then
        GVCFS="--variant $INPUTFILE"
    else
        GVCF_DIR=${INPUTFILE%/*}
        for f in $GVCF_DIR/*.gvcf; do 
            GVCFS=$GVCFS"--variant $f " ;
        done
    fi

    java $JAVAPARAMS -jar $PATH_GATK_JAR -l WARN\
    -T GenotypeGVCFs \
    -nt 20 \
    -R $REFERENCE \
    --dbsnp $DBSNPVCF \
    -o $OUTDIR/raw_snps_raw_indels.vcf \
    $GVCFS

    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/raw_snps_raw_indels.vcf
fi
################################################################################
NGSANE_CHECKPOINT_INIT "training snp model"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    echo "Training SNP model"

    java $JAVAPARAMS -jar $PATH_GATK_JAR -l WARN \
    -T VariantRecalibrator \
    -R $REFERENCE \
    -input $OUTDIR/raw_snps_raw_indels.vcf \
    -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $HAPMAPVCF \
    -resource:omni,known=false,training=true,truth=true,prior=12.0 $ONEKOMNI \
    -resource:1000G,known=false,training=true,truth=false,prior=10.0 $ONEKSNPS \
    -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $DBSNPVCF \
    -nt 20 \
    $PARAM_STRING \
    -mode SNP \
    -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
    -recalFile recalibrate_SNP.recal \
    -tranchesFile recalibrate_SNP.tranches \
    -rscriptFile recalibrate_SNP_plots.R

    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/recalibrate_SNP.recal $OUTDIR/recalibrate_SNP.tranches $OUTDIR/recalibrate_SNP_plots.R

fi    
################################################################################
NGSANE_CHECKPOINT_INIT "recalibrating snps"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then
    echo "recalibrating snps"

    java $JAVAPARAMS -jar $PATH_GATK_JAR -l WARN\
    -T ApplyRecalibration \
    -R $REFERENCE \
    -input $OUTDIR/raw_snps_raw_indels.vcf \
    -nt 20 \ 
    -mode SNP \
    --ts_filter_level $TS_FILTER \
    -recalFile $OUTDIR/recalibrate_SNP.recal \
    -tranchesFile $OUTDIR/recalibrate_SNP.tranches \
    -o $OUTDIR/recalibrated_snps_raw_indels.vcf

    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/recalibrated_snps_raw_indels.vcf

fi
################################################################################
NGSANE_CHECKPOINT_INIT "training indel model"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then
    echo "Training INDEL model"

    java $JAVAPARAMS -jar $PATH_GATK/GenomeAnalysisTK.jar -l WARN \
    -T VariantRecalibrator \
    -R $REFERENCE \
    -input $OUTDIR/recalibrated_snps_raw_indels.vcf \
    -resource:mills,known=true,training=true,truth=true,prior=12.0 $MILLS \
    -nt 20 \
    $PARAM_STRING \
    -mode INDEL \
    -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
    --maxGaussians $MAXGAUS \
    -recalFile $OUTDIR/recalibrate_INDEL.recal \
    -tranchesFile $OUTDIR/recalibrate_INDEL.tranches \
    -rscriptFile $OUTDIR/recalibrate_INDEL_plots.R

    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/recalibrate_INDEL.recal $OUTDIR/recalibrate_INDEL.tranches

fi
################################################################################
NGSANE_CHECKPOINT_INIT "recalibrating indels"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then
    echo "recalibrating indels"

    java $JAVAPARAMS -jar $PATH_GATK_JAR -l WARN \
    -T ApplyRecalibration \
    -R $REFERENCE \
    -input $OUTDIR/recalibrated_snps_raw_indels.vcf \
    -mode INDEL \
    -nt 20 \
    --ts_filter_level $TS_FILTER \
    -recalFile $OUTDIR/recalibrate_INDEL.recal \
    -tranchesFile $OUTDIR/recalibrate_INDEL.tranches \
    -o $OUTDIR/recalibrated_variants.vcf

    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/recalibrated_variants.vcf

fi

################################################################################
[ -e $OUTDIR/recalibrated_variants.vcf.dummy ] && rm $OUTDIR/recalibrated_variants.vcf.dummy
echo ">>>>> GATK jointGenotyping - FINISHED"
echo ">>>>> enddate "`date`
