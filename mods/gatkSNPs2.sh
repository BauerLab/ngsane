#!/bin/bash

# Script running snp calling with GATK version 2.5
# QC:
# author: Denis C. Bauer
# date: May.2013


echo ">>>>> call SNPs using GATK"
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"

function usage {
echo -e "usage: $(basename $0)

required:
  -k | --toolkit <path>     location of the NGSANE repository 
  -i | --input <file>       temp files with all bam files for SNP calls
  -o | --outdir <path>      output dir
  -r | --reference <file>   reference genome

options:
  -t | --threads <nr>       number of CPUs to use (default: 1)
  -R | --region <ps>        region of specific interest, e.g. targeted reseq
                             format chr:pos-pos
  -d | --dbsnp <file> 
"
exit
}


if [ ! $# -gt 3 ]; then usage ; fi

#DEFAULTS
THREADS=1
CALLSNPS="1"
HARDFILTER="1"
VARIANTRECAL="1"

ADDRECAL="" # additional commands for variant recalibrator

#INPUTS
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the NGSANE repository
        -t | --threads )        shift; THREADS=$1 ;; # number of CPUs to use
        -i | --input )          shift; FILES=$1 ;; # temp files with paths to bam file
        -o | --outdir )         shift; MYOUT=$1 ;; # output dir
        -n | --name )           shift; NAME=$1 ;; # name
        -r | --reference )      shift; FASTA=$1 ;; # reference genome
        -g | --refseq )         shift; REFSEQROD=$1 ;; # refseq genome
        -H | --hapmap )         shift; HAPMAPVCF=$1 ;; # hapmap
        -d | --dbsnp )          shift; DBSNPVCF=$1 ;; # dbsnp
        -K | --1kg )            shift; ONEKGVCF=$1 ;; # 1000genomes data
        -L | --region )         shift; SEQREG=$1 ;; # (optional) region of specific interest, e.g. targeted reseq
        --maxGaussians )        shift; ADDRECAL=$ADDRECAL" --maxGaussians "$1 ;; #(additional params for recal)
        --percentBadVariants )  shift; ADDRECAL=$ADDRECAL" --percentBadVariants "$1 ;; #(additional params for recal)
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

for MODULE in $MODULE_GATKSNP; do module load $MODULE; done  # save way to load modules that itself load other modules
export PATH=$PATH_GATKSNP:$PATH
module list
echo $PATH
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)
PATH_GATK=$(dirname $(which GenomeAnalysisTK.jar))
PATH_IGVTOOLS=$(dirname $(which igvtools.jar))

echo "[NOTE] set java parameters"
JAVAPARAMS="-Xmx"$(python -c "print int($MEMORY_VAR*0.8)")"g -Djava.io.tmpdir="$TMP"  -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1" 
unset _JAVA_OPTIONS
echo "JAVAPARAMS "$JAVAPARAMS

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--JAVA        --\n" $(java -Xmx200m -version 2>&1)
[ -z "$(which java)" ] && echo "[ERROR] no java detected" && exit 1
echo -e "--R           --\n "$(R --version | head -n 3)
[ -z "$(which R)" ] && echo "[ERROR] no R detected" && exit 1
echo -e "--igvtools    --\n "$(java $JAVAPARAMS -jar $PATH_IGVTOOLS/igvtools.jar version 2>&1)
[ ! -f $PATH_IGVTOOLS/igvtools.jar ] && echo "[ERROR] no igvtools detected" && exit 1
echo -e "--GATK        --\n "$(java $JAVAPARAMS -jar $PATH_GATK/GenomeAnalysisTK.jar --version)
[ ! -f $PATH_GATK/GenomeAnalysisTK.jar ] && echo "[ERROR] no GATK detected" && exit 1


echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="parameters"

if [ ! -d $MYOUT ]; then mkdir -p $MYOUT; fi

# delete old snp files unless attempting to recover
if [ -z "$RECOVERFROM" ]; then
    [ -e $MYOUT/$NAME.fi.vcf ] && rm $MYOUT/$NAME.fi.vcf
    [ -e $MYOUT/$NAME.fi.vcf.idx ] && rm $MYOUT/$NAME.fi.vcf.idx
    [ -e $MYOUT/gatkSNPcall.tmp ] && rm $MYOUT/gatkSNPcall.tmp
fi
        
if [ -n "$SEQREG" ]; then REGION="-L $SEQREG"; fi
        
echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="recall files from tape"

if [ -n "$DMGET" ]; then 
    dmget -a ${FILES//,/ }; 
fi
    
echo -e "\n********* $CHECKPOINT\n"    
################################################################################
CHECKPOINT="call snps"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
    echo "[NOTE] $CHECKPOINT"
    
    if [ -n "$CALLSNPS" ]; then

        # -nt $THREADS <- it is not parallele (2012)
        # from new versions add --computeSLOD
        # http://seqanswers.com/forums/showthread.php?t=14836
        echo "[NOTE] call SNPs and VariantAnnotation"
        echo "java $JAVAPARAMS -jar $PATH_GATK/GenomeAnalysisTK.jar -l INFO \
             -T UnifiedGenotyper \
             -glm BOTH \
             -R $FASTA \
             --dbsnp $DBSNPVCF \
             -A HomopolymerRun \
             -A MappingQualityRankSumTest \
             -A Coverage \
             -A QualByDepth \
             -A RMSMappingQuality \
             -A SpanningDeletions \
             -A HaplotypeScore \
             -A AlleleBalance \
             -A BaseQualityRankSumTest \
             -A MappingQualityZero \
             --out $MYOUT/$NAME.raw.vcf \
             -stand_call_conf 30.0 \
             $REGION \
             -stand_emit_conf 10.0 \\" > $MYOUT/gatkVarcall.tmp
             
        for f in ${FILES//,/ }; do echo "-I $f \\" >>$MYOUT/gatkVarcall.tmp ; done
        echo " -dcov 1000 " >> $MYOUT/gatkVarcall.tmp
    
        # set up to execute
        echo "rm $MYOUT/gatkVarcall.tmp" >> $MYOUT/gatkVarcall.tmp
        chmod -u=rwx $MYOUT/gatkVarcall.tmp
        $MYOUT/gatkVarcall.tmp
    
    
        echo "[NOTE] get snps only"
        java $JAVAPARAMS -jar $PATH_GATK/GenomeAnalysisTK.jar -l WARN \
    	-T SelectVariants \
    	-R $FASTA \
    	--variant  $MYOUT/$NAME.raw.vcf \
    	--selectTypeToInclude SNP \
    	-o $MYOUT/$NAME.raw.snps.vcf
    
        echo "[NOTE] get indels only"
        java $JAVAPARAMS -jar $PATH_GATK/GenomeAnalysisTK.jar -l WARN \
    	-T SelectVariants \
    	-R $FASTA \
    	--variant  $MYOUT/$NAME.raw.vcf \
    	--selectTypeToInclude INDEL \
    	-o $MYOUT/$NAME.raw.indel.vcf

        echo "[NOTE] SNP call done "`date`

        # mark checkpoint
        if [ -f $MYOUT/$NAME.raw.vcf ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

    fi
fi 

################################################################################
CHECKPOINT="hardfilter"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
    echo "[NOTE] $CHECKPOINT"
        
    if [ -n "$HARDFILTER" ]; then
    
        echo "[NOTE] hard filter SNPs"
        java $JAVAPARAMS -jar $PATH_GATK/GenomeAnalysisTK.jar -l WARN \
    	-T VariantFiltration \
    	-R $FASTA \
    	-o $MYOUT/$NAME.filter.snps.vcf \
    	--variant $MYOUT/$NAME.raw.snps.vcf \
    	$REGION \
    	--mask $MYOUT/$NAME.raw.indel.vcf \
    	--maskName InDel \
    	--clusterWindowSize 10 \
            --filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * (DP+1))) > 0.1)" \
            --filterName "HARD_TO_VALIDATE" \
    	--filterExpression "QUAL < 30.0 || QD < 5.0 || HRun > 5 || SB > -0.10" \
    	--filterName "GATKStandHiCovExomes"
    
    
    #	--filterExpression "AF < 0.2" \
    #	--filterName "AllelFreq" \
    #	--genotypeFilterExpression "DP < 20" \
    #	--genotypeFilterName "Depth" \
    #	--filterExpression "MQ0 > 50" \
    #	--filterName "AccMQ0" \
    #	--filterExpression "SB > -1.0" \
    #	--filterName "StrandBias"
    
    
        echo "[NOTE] hard filter INDELs"
        java $JAVAPARAMS -jar $PATH_GATK/GenomeAnalysisTK.jar -l WARN \
    	-T VariantFiltration \
    	-R $FASTA \
    	-o $MYOUT/$NAME.filter.indel.vcf \
    	--variant $MYOUT/$NAME.raw.indel.vcf \
    	$REGION \
    	--filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * (DP+1))) > 0.1)" \
    	--filterName "HARD_TO_VALIDATE" \
    	--filterExpression "SB >= -1.0" \
    	--filterName "StrandBiasFilter" \
    	--filterExpression "QUAL < 10" \
    	--filterName "QualFilter"
    
        echo "[NOTE] Hard filter eval SNPs"
        java $JAVAPARAMS -jar $PATH_GATK/GenomeAnalysisTK.jar -l WARN \
                -T VariantEval \
                -R $FASTA \
                --dbsnp $DBSNPVCF \
                --eval $MYOUT/$NAME.filter.snps.vcf \
    	       $REGION \
                --evalModule TiTvVariantEvaluator \
                -o $MYOUT/$NAME.filter.snps.eval.txt
    
        echo "[NOTE] Hard filter eval INDELs"
        java $JAVAPARAMS -jar $PATH_GATK/GenomeAnalysisTK.jar -l WARN \
                -T VariantEval \
                -R $FASTA \
                --dbsnp $DBSNPVCF \
                --eval $MYOUT/$NAME.filter.indel.vcf \
    	    $REGION \
                --evalModule TiTvVariantEvaluator \
                -o $MYOUT/$NAME.filter.indel.eval.txt
    
    
        echo "[NOTE] hard filter "`date`

        # mark checkpoint
        if [ -f $MYOUT/$NAME.filter.snps.vcf ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

    fi
fi 

################################################################################
CHECKPOINT="re-calibrate"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
    echo "[NOTE] $CHECKPOINT"
    
    #####################
    # Recalibration
    # http://www.broadinstitute.org/gsa/wiki/index.php/Variant_quality_score_recalibration
    #####################
    
    if [ -n "$VARIANTRECAL" ]; then
    
        echo $MYOUT
        echo $DBSNPVCF
        echo $HAPMAPVCF
        echo $ONEKGVCF
    
    	echo "[NOTE] Recalibrate"
    	if [ ! -e $MYOUT/R ]; then mkdir $MYOUT/R; fi
    
    	echo "[NOTE] train"
    	# maxGaussians 6
    	# no percent bad
    	#-an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an HRun -an InbreedingCoeff\
    	java $JAVAPARAMS -jar $PATH_GATK/GenomeAnalysisTK.jar -l INFO \
    	    -T VariantRecalibrator \
    	    -R $FASTA \
    	    --input $MYOUT/$NAME.raw.vcf \
    	    -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $HAPMAPVCF \
    	    -resource:omni,known=false,training=true,truth=false,prior=12.0 $ONEKGVCF \
    	    -resource:dbsnp,known=true,training=false,truth=false,prior=8.0 $DBSNPVCF \
    	    -an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an MQ \
    	    -mode BOTH \
    	    -recalFile $MYOUT/$NAME.raw.recal \
    	    -tranchesFile $MYOUT/$NAME.raw.tranches \
    	    $ADDRECAL \
    	    -rscriptFile $MYOUT/R/output.plots.R \
    
    #	    -nt $THREADS <- is notparallele
    
    
    	echo "[NOTE] variant cut"
    	java $JAVAPARAMS -jar $PATH_GATK/GenomeAnalysisTK.jar -l WARN \
    	    -T ApplyRecalibration \
    	    -mode Both \
    	    -R $FASTA \
    	    -input $MYOUT/$NAME.raw.vcf \
    	    --ts_filter_level 99.0 \
    	    -tranchesFile $MYOUT/$NAME.raw.tranches \
    	    -recalFile $MYOUT/$NAME.raw.recal \
    	    -o $MYOUT/$NAME.recalfilt.vcf
    
    
    	echo "[NOTE] Recal eval variants"
    	java $JAVAPARAMS -jar $PATH_GATK/GenomeAnalysisTK.jar -l WARN \
                -T VariantEval \
                -R $FASTA \
                --dbsnp $DBSNPVCF \
                --eval $MYOUT/$NAME.recalfilt.snps.vcf \
                --evalModule TiTvVariantEvaluator \
                -o $MYOUT/$NAME.recalfilt.eval.txt
    
        # mark checkpoint
        if [ -f $MYOUT/$NAME.recalfilt.eval.txt ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

    fi
fi 

################################################################################
CHECKPOINT="index for IGV"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
    echo "[NOTE] $CHECKPOINT"
    
    java $JAVAPARAMS -jar $PATH_IGVTOOLS/igvtools.jar index $MYOUT/${n/bam/fi.vcf}

    # mark checkpoint
    echo -e "\n********* $CHECKPOINT\n" && unset RECOVERFROM
fi 
################################################################################
echo ">>>>> call SNPs using GATK - FINISHED"
echo ">>>>> enddate "`date`
