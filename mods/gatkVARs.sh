#!/bin/bash -e

# Script running snp calling with GATK version 2.5
# QC:
# author: Denis C. Bauer
# date: May.2013
# RESULTFILENAME <TASK>/${GATKVAR_AGGREGATE_FOLDER}/<DIR><ADDDUMMY>_snps.flt.vcf

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
HARDFILTER="1"
VARIANTRECAL="1"

ADDRECALALN="" # additional commands for variant recalibrator

#INPUTS
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the NGSANE repository
        -t | --threads )        shift; THREADS=$1 ;; # number of CPUs to use
        -i | --input )          shift; FILES=$1 ;; # temp files with paths to bam file
        -o | --outdir )         shift; OUTDIR=$1 ;; # output dir
        -n | --name )           shift; NAME=$1 ;; # name
        -r | --reference )      shift; FASTA=$1 ;; # reference genome
        -g | --refseq )         shift; REFSEQROD=$1 ;; # refseq genome
        -H | --hapmap )         shift; HAPMAPVCF=$1 ;; # hapmap
        -d | --dbsnp )          shift; DBSNPVCF=$1 ;; # dbsnp
        -K | --1kg )            shift; ONEKGVCF=$1 ;; # 1000genomes data
        -L | --region )         shift; SEQREG=$1 ;; # (optional) region of specific interest, e.g. targeted reseq
        --maxGaussians )        shift; ADDRECALALN=$ADDRECALALN" --maxGaussians "$1 ;; #(additional params for recal)
        --percentBadVariants )  shift; ADDRECALALN=$ADDRECALALN" --percentBadVariants "$1 ;; #(additional params for recal)
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
hash module 2>/dev/null && for MODULE in $MODULE_GATKVAR; do module load $MODULE; done && module list 

export PATH=$PATH_GATKVAR:$PATH
echo $PATH
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)
[ -z "$PATH_GATK" ] && PATH_GATK=$(dirname $(which GenomeAnalysisTK.jar))
[ -z "$PATH_IGVTOOLS" ] && PATH_IGVTOOLS=$(dirname $(which igvtools.jar))

echo "[NOTE] set java parameters"
JAVAPARAMS="-Xmx"$(python -c "print int($MEMORY_GATKVAR*0.75)")"g -Djava.io.tmpdir="$TMP"  -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1" 
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


NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "parameters"

if [ ! -d $OUTDIR ]; then mkdir -p $OUTDIR; fi

# delete old snp files unless attempting to recover
if [ -z "$NGSANE_RECOVERFROM" ]; then
    [ -e $OUTDIR/${NAME}_snps.flt.vcf ] && rm $OUTDIR/${NAME}*
fi
        
if [ -z "$DBSNPVCF" ] || [ ! -e $DBSNPVCF ]; then
    echo "[ERRPR] DBSNPVCF parameter not specified or file not found"
    exit 1
fi
if [ -z "$HAPMAPVCF" ] || [ ! -e $HAPMAPVCF ]; then
    echo "[ERRPR] HAPMAPVCF parameter not specified or file not found"
    exit 1
fi
# ONEKGVCF only optional
#if [ -z "$ONEKGVCF" ] || [ ! -e $ONEKGVCF ]; then
#    echo "[ERRPR] ONEKGVCF parameter not specified or file not found"
#    exit 1
#fi
        
if [ -n "$SEQREG" ]; then REGION="-L $SEQREG"; fi


if [[ $(which GenomeAnalysisTK.jar) =~ "2.8" ]]; then 
        echo "[NOTE] new GATK parallele"
        PARALLELENCT="-nct $CPU_RECALALN"
		PARALLELENT="-nt $CPU_RECALALN"
fi

        
NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "recall files from tape"

if [ -n "$DMGET" ]; then 
    dmget -a ${FILES//,/ };
	dmget -a $OUTDIR/*
fi
    
NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "call snps"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then
    
	for f in ${FILES//,/ }; do LIST=$LIST"-I $f " ; done
	
    # -nt $THREADS <- it is not parallele (2012)
    # from new versions add --computeSLOD
    # http://seqanswers.com/forums/showthread.php?t=14836
    java $JAVAPARAMS -jar $PATH_GATK/GenomeAnalysisTK.jar -l INFO \
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
            --out $OUTDIR/${NAME}.raw.vcf \
            -stand_call_conf 30.0 \
            $REGION \
			$PARALLELENT \
            -stand_emit_conf 10.0 \
			$LIST \
			-dcov 1000

    echo "[NOTE] get snps only"
    java $JAVAPARAMS -jar $PATH_GATK/GenomeAnalysisTK.jar -l WARN \
	-T SelectVariants \
	-R $FASTA \
	--variant  $OUTDIR/${NAME}.raw.vcf \
	--selectTypeToInclude SNP \
	$PARALLELENT \
	-o $OUTDIR/${NAME}_snps.raw.vcf

    echo "[NOTE] get indels only"
    java $JAVAPARAMS -jar $PATH_GATK/GenomeAnalysisTK.jar -l WARN \
	-T SelectVariants \
	-R $FASTA \
	--variant  $OUTDIR/${NAME}.raw.vcf \
	--selectTypeToInclude INDEL \
	$PARALLELENT \
	-o $OUTDIR/${NAME}_indel.raw.vcf

    echo "[NOTE] SNP call done "`date`

    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/${NAME}.raw.vcf

fi 

################################################################################
NGSANE_CHECKPOINT_INIT "hardfilter"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then
        
    java $JAVAPARAMS -jar $PATH_GATK/GenomeAnalysisTK.jar -l WARN \
	-T VariantFiltration \
	-R $FASTA \
	-o $OUTDIR/${NAME}_snps.flt.vcf \
	--variant $OUTDIR/${NAME}_snps.raw.vcf \
	$REGION \
	--mask $OUTDIR/${NAME}_indel.raw.vcf \
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
	-o $OUTDIR/${NAME}_indel.flt.vcf \
	--variant $OUTDIR/${NAME}_indel.raw.vcf \
	$REGION \
	--filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * (DP+1))) > 0.1)" \
	--filterName "HARD_TO_VALIDATE" \
	--filterExpression "SB >= -1.0" \
	--filterName "StrandBiasFilter" \
	--filterExpression "QUAL < 10" \
	--filterName "QualFilter"

    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/${NAME}_snps.flt.vcf $OUTDIR/${NAME}_indel.flt.vcf

    rm $OUTDIR/${NAME}_snps.raw.vcf* $OUTDIR/${NAME}_indel.raw.vcf*


fi 

################################################################################
NGSANE_CHECKPOINT_INIT "index for IGV"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then
    
    java $JAVAPARAMS -jar $PATH_IGVTOOLS/igvtools.jar index $OUTDIR/${NAME}_snps.flt.vcf
    java $JAVAPARAMS -jar $PATH_IGVTOOLS/igvtools.jar index $OUTDIR/${NAME}_indel.flt.vcf

    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/${NAME}_snps.flt.vcf.idx $OUTDIR/${NAME}_indel.flt.vcf.idx
fi 

################################################################################
NGSANE_CHECKPOINT_INIT "evaluate"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    echo "[NOTE] Hard filter eval SNPs"
    java $JAVAPARAMS -jar $PATH_GATK/GenomeAnalysisTK.jar -l WARN \
            -T VariantEval \
            -R $FASTA \
            --dbsnp $DBSNPVCF \
            --eval $OUTDIR/${NAME}_snps.flt.vcf \
	       $REGION \
            --evalModule TiTvVariantEvaluator \
			$PARALLELENT \
            -o $OUTDIR/${NAME}_snps.flt.eval.txt

    echo "[NOTE] Hard filter eval INDELs"
    java $JAVAPARAMS -jar $PATH_GATK/GenomeAnalysisTK.jar -l WARN \
            -T VariantEval \
            -R $FASTA \
            --dbsnp $DBSNPVCF \
            --eval $OUTDIR/${NAME}_indel.flt.vcf \
	       $REGION \
            --evalModule TiTvVariantEvaluator \
			$PARALLELENT \
            -o $OUTDIR/${NAME}_indel.flt.eval.txt

    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/${NAME}_snps.flt.eval.txt $OUTDIR/${NAME}_indel.flt.eval.txt

fi

################################################################################
NGSANE_CHECKPOINT_INIT "re-calibrate"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then
    
    #####################
    # Recalibration
    # http://www.broadinstitute.org/gsa/wiki/index.php/Variant_quality_score_recalibration
    #####################
    
    if [ -n "$VARIANTRECAL" ]; then
    
        echo $OUTDIR
        echo $DBSNPVCF
        echo $HAPMAPVCF
        echo $ONEKGVCF
	    if [ -n "$ONEKGVCF" ]; then 
	        ONEKGPARAMS="-resource:omni,known=false,training=true,truth=false,prior=12.0 $ONEKGVCF"
	    fi

    	if [ ! -e $OUTDIR/R ]; then mkdir $OUTDIR/R; fi
    
    	echo "[NOTE] train"
    	# maxGaussians 6
    	# no percent bad
    	#-an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an HRun -an InbreedingCoeff\
    	java $JAVAPARAMS -jar $PATH_GATK/GenomeAnalysisTK.jar -l INFO \
    	    -T VariantRecalibrator \
    	    -R $FASTA \
    	    --input $OUTDIR/$NAME.raw.vcf \
    	    -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $HAPMAPVCF \
    	    -resource:dbsnp,known=true,training=false,truth=false,prior=8.0 $DBSNPVCF \
	    $ONEKGPARAMS \
    	    -an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an MQ \
    	    -mode BOTH \
    	    -recalFile $OUTDIR/${NAME}.raw.recal \
    	    -tranchesFile $OUTDIR/${NAME}.raw.tranches \
    	    $REGION \
    	    $ADDRECALALN \
			$PARALLELENT \
    	    -rscriptFile $OUTDIR/R/output.plots.R \
    
    #	    -nt $THREADS <- is notparallele
    
    
    	echo "[NOTE] variant cut"
    	java $JAVAPARAMS -jar $PATH_GATK/GenomeAnalysisTK.jar -l WARN \
    	    -T ApplyRecalibration \
    	    -mode Both \
    	    -R $FASTA \
    	    -input $OUTDIR/${NAME}.raw.vcf \
    	    --ts_filter_level 99.0 \
    	    -tranchesFile $OUTDIR/${NAME}.raw.tranches \
    	    -recalFile $OUTDIR/${NAME}.raw.recal \
    	    $REGION \
			$PARALLELENT \
    	    -o $OUTDIR/${NAME}.recalflt.vcf
    
        rm $OUTDIR/${NAME}.raw.recal* $OUTDIR/${NAME}.raw.tranches
    
    	echo "[NOTE] Recal eval variants"
    	java $JAVAPARAMS -jar $PATH_GATK/GenomeAnalysisTK.jar -l WARN \
                -T VariantEval \
                -R $FASTA \
                --dbsnp $DBSNPVCF \
                --eval $OUTDIR/${NAME}.recalflt.vcf \
                --evalModule TiTvVariantEvaluator \
				$PARALLELENT \
                -o $OUTDIR/${NAME}.recalflt.eval.txt
    
        # mark checkpoint
        NGSANE_CHECKPOINT_CHECK $OUTDIR/${NAME}.recalflt.eval.txt
    else
        echo "[NOTE] re-calibrate skipped"
        NGSANE_CHECKPOINT_CHECK
    fi
fi 

################################################################################
[ -e $OUTDIR/${NAME}_snps.flt.vcf.dummy ] && rm $OUTDIR/${NAME}_snps.flt.vcf.dummy
echo ">>>>> call SNPs using GATK - FINISHED"
echo ">>>>> enddate "`date`
