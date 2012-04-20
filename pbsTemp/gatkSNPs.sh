#!/bin/bash

# Script running snp calling with GATK
# QC:
# author: Denis C. Bauer
# date: Nov.2010


echo ">>>>> call SNPs using GATK"
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> gatkSNPs.sh $*"

function usage {
echo -e "usage: $(basename $0)

required:
  -k | --toolkit <path>     location of the HiSeqInf repository 
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
        -k | --toolkit )        shift; HISEQINF=$1 ;; # location of the HiSeqInf repository
        -t | --threads )        shift; THREADS=$1 ;; # number of CPUs to use
        -i | --input )          shift; FILES=$1 ;; # temp files with paths to bam file
        -o | --outdir )         shift; OUT=$1 ;; # output dir
	-n | --name )           shift; NAME=$1 ;; # name
        -r | --reference )      shift; FASTA=$1 ;; # reference genome
	-g | --refseq )         shift; REFSEQROD=$1 ;; # refseq genome
	-H | --hapmap )         shift; HAPMAPVCF=$1 ;; # hapmap
	-d | --dbsnp )          shift; DBSNPVCF=$1 ;; # dbsnp
	-K | --1kg )            shift; ONEKGVCF=$1 ;; # 1000genomes data
	-L | --region )         shift; SEQREG=$1 ;; # (optional) region of specific interest, e.g. targeted reseq
	--maxGaussians )        shift; ADDRECAL=$ADDRECAL" --maxGaussians "$1 ;; #(additional params for recal)
	--percentBadVariants )  shift; ADDRECAL=$ADDRECAL" --percentBadVariants "$1 ;; #(additional params for recal)
        -h | --help )           usage ;;
        * )                     usage
    esac
    shift
done



#PROGRAMS
. $HISEQINF/pbsTemp/header.sh

if [ -n "$CALLSNPS" ]; then

   # delete old snp files
    if [ -e $OUT/$NAME.fi.vcf} ]; then
	rm $OUT/$NAME.fi.vcf}
	rm $OUT/$NAME.fi.vcf.idx
    fi

    if [ -e $OUT/gatkSNPcall.tmp ]; then rm $OUT/gatkSNPcall.tmp; fi
    if [ -n "$SEQREG" ]; then REGION="-L $SEQREG"; fi


    # -nt $THREADS \
    # from new versions add --computeSLOD
    # http://seqanswers.com/forums/showthread.php?t=14836
    echo "********* call SNPs and VariantAnnotation"
    echo "java -Xmx10g -jar $GATKJAR/GenomeAnalysisTK.jar -l INFO \
         -T UnifiedGenotyper \
         -glm BOTH \
         -R $FASTA \
         -B:dbsnp,VCF $DBSNPVCF \
         -A HomopolymerRun \
         -A MappingQualityRankSumTest \
         -A DepthOfCoverage \
         -A QualByDepth \
         -A RMSMappingQuality \
         -A SpanningDeletions \
         -A HaplotypeScore \
         -A AlleleBalance \
         -A BaseQualityRankSumTest \
         -A MappingQualityZero \
         --out $OUT/$NAME.raw.vcf \
         -stand_call_conf 30.0 \
         $REGION \
         -stand_emit_conf 10.0 \\" > $OUT/gatkVarcall.tmp
    for f in $( less $FILES ); do echo "-I $f \\" >>$OUT/gatkVarcall.tmp ; done
    echo " -dcov 1000 " >> $OUT/gatkVarcall.tmp

    # set up to execute
    echo "rm $OUT/gatkVarcall.tmp" >> $OUT/gatkVarcall.tmp
    chmod -u=rwx $OUT/gatkVarcall.tmp
    $OUT/gatkVarcall.tmp


    echo "********* get snps only"
    java -Xmx10g -jar $GATKJAR/GenomeAnalysisTK.jar -l WARN \
	-T SelectVariants \
	-R $FASTA \
	-B:variant,VCF  $OUT/$NAME.raw.vcf \
	-snps \
	-o $OUT/$NAME.raw.snps.vcf

    echo "********* get indels only"
    java -Xmx10g -jar $GATKJAR/GenomeAnalysisTK.jar -l WARN \
	-T SelectVariants \
	-R $FASTA \
	-B:variant,VCF  $OUT/$NAME.raw.vcf \
	-indels \
	-o $OUT/$NAME.raw.indel.vcf


fi

echo ">>>>> SNP call done "`date`


#####################
# Hard filter
#####################

if [ -n "$HARDFILTER" ]; then

    echo "********* HardFilter"


    echo "********* hard filter SNPs"
    java -Xmx10g -jar $GATKJAR/GenomeAnalysisTK.jar -l WARN \
	-T VariantFiltration \
	-R $FASTA \
	-o $OUT/$NAME.filter.snps.vcf \
	-B:variant,VCF $OUT/$NAME.raw.snps.vcf \
	-B:mask,VCF $OUT/$NAME.raw.indel.vcf \
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


    echo "********* hard filter INDELs"
    java -Xmx10g -jar $GATKJAR/GenomeAnalysisTK.jar -l WARN \
	-T VariantFiltration \
	-R $FASTA \
	-o $OUT/$NAME.filter.indel.vcf \
	-B:variant,VCF $OUT/$NAME.raw.indel.vcf \
	--filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * (DP+1))) > 0.1)" \
	--filterName "HARD_TO_VALIDATE" \
	--filterExpression "SB >= -1.0" \
	--filterName "StrandBiasFilter" \
	--filterExpression "QUAL < 10" \
	--filterName "QualFilter"



    echo "********* Hard filter eval SNPs"
    java -Xmx10g -jar $GATKJAR/GenomeAnalysisTK.jar -l WARN \
            -T VariantEval \
            -R $FASTA \
            -B:dbsnp,vcf $DBSNPVCF \
            -B:eval,VCF $OUT/$NAME.filter.snps.vcf \
            --evalModule TiTvVariantEvaluator \
            -o $OUT/$NAME.filter.snps.eval.txt

    echo "********* Hard filter eval INDELs"
    java -Xmx10g -jar $GATKJAR/GenomeAnalysisTK.jar -l WARN \
            -T VariantEval \
            -R $FASTA \
            -B:dbsnp,vcf $DBSNPVCF \
            -B:eval,VCF $OUT/$NAME.filter.indel.vcf \
            --evalModule TiTvVariantEvaluator \
            -o $OUT/$NAME.filter.indel.eval.txt


echo ">>>>> hard filter "`date`

fi

#####################
# Recalibration
# http://www.broadinstitute.org/gsa/wiki/index.php/Variant_quality_score_recalibration
#####################

if [ -n "$VARIANTRECAL" ]; then

    echo $OUT
    echo $DBSNPVCF
    echo $HAPMAPVCF
    echo $ONEKGVCF



	echo "********* Recalibrate"
	if [ ! -e $OUT/R ]; then mkdir $OUT/R; fi

	echo "********* train"
	# maxGaussians 6
	# no percent bad
	#-an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an HRun -an InbreedingCoeff\
	java -Xmx10g -jar $GATKJAR/GenomeAnalysisTK.jar -l INFO \
	    -T VariantRecalibrator \
	    -R $FASTA \
	    -B:input,VCF $OUT/$NAME.raw.snps.vcf \
	    -B:hapmap,VCF,known=false,training=true,truth=true,prior=15.0 $HAPMAPVCF \
	    -B:omni,VCF,known=false,training=true,truth=false,prior=12.0 $ONEKGVCF \
	    -B:dbsnp,VCF,known=true,training=false,truth=false,prior=8.0 $DBSNPVCF \
	    -an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an MQ \
	    -mode SNP \
	    -recalFile $OUT/$NAME.raw.recal \
	    -tranchesFile $OUT/$NAME.raw.tranches \
	    -Rscript $RSCRIPT \
	    $ADDRECAL \
	    -rscriptFile $OUT/R/output.plots.R \
	    -resources $GATKHOME/public/R/ 

#	    -nt $THREADS <- is notparallele


	echo "********* variant cut"
	java -Xmx10g -jar $GATKJAR/GenomeAnalysisTK.jar -l WARN \
	    -T ApplyRecalibration \
	    -R $FASTA \
	    -B:input,VCF $OUT/$NAME.raw.snps.vcf \
	    --ts_filter_level 99.0 \
	    -tranchesFile $OUT/$NAME.raw.tranches \
	    -recalFile $OUT/$NAME.raw.recal \
	    -o $OUT/$NAME.recalfilt.snps.vcf


	echo "********* Recal eval Variants"
	java -Xmx10g -jar $GATKJAR/GenomeAnalysisTK.jar -l WARN \
            -T VariantEval \
            -R $FASTA \
            -B:dbsnp,vcf $DBSNPVCF \
            -B:eval,VCF $OUT/$NAME.recalfilt.snps.vcf \
            --evalModule TiTvVariantEvaluator \
            -o $OUT/$NAME.recalfilt.snps.eval.txt


    

fi



#echo "********* make index for IGV"
#java -Xmx1g -jar $IGVTOOLS index $OUT/${n/bam/fi.vcf}

echo ">>>>> call SNPs using GATK - FINISHED"
echo ">>>>> enddate "`date`