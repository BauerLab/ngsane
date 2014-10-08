#!/bin/bash

# Call variants with pindel
# author: Denis C. Bauer
# date: Dec.2013

# messages to look out for -- relevant for the QC.sh script:
# 

echo ">>>>> Collect Variants after calling with any variant collection tool "
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"


function usage {
echo -e "usage: $(basename $0) -k CONFIG -f BAM -o OUTDIR [OPTIONS]

Variant calling with pindel

required:
  -k | --toolkit <path>     config file
  -f <files>                bam files
  -o <dir>
"
exit
}


if [ ! $# -gt 3 ]; then usage ; fi

#DEFAULTS
TARGET="vcf"

#INPUTS
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the NGSANE repository
        -f             )        shift; FILES=${1//,/ } ;; # bam files
		-i1			   )		shift; STAGEONEINPUT=${1//,/ } ;; # input dir from prev. stage
		-i2			   )		shift; STAGETWOINPUT=${1//,/ } ;; # input from stage two
        -o             )        shift; OUTDIR=$1 ;; # outputdir
		--dummy		   )		shift; DUMMY=$1 ;; # ending of dummy -- only necessary if target is different
		--target	   )		shift; TARGET=$1 ;; # ending of target
#        -L | --region )         shift; SEQREG=$1 ;; # (optional) region of specific interest, e.g. targeted reseq
        --recover-from )        shift; NGSANE_RECOVERFROM=$1 ;; # attempt to recover from log file
        -h | --help    )        usage ;;
        * )                     echo "don't understand "$1
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
hash module 2>/dev/null && for MODULE in $MODULE_VARCOLLECT; do module load $MODULE; done && module list

export PATH=$PATH_VARCOLLECT:$PATH
echo "PATH=$PATH"
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)
PATH_GATK=$(dirname $(which GenomeAnalysisTK.jar))
PATH_IGVTOOLS=$(dirname $(which igvtools.jar))

echo "[NOTE] set java parameters"
JAVAPARAMS="-Xmx"$(python -c "print int($MEMORY_VARCOLLECT*0.75)")"g -Djava.io.tmpdir="$TMP"  -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1" 
unset _JAVA_OPTIONS
echo "JAVAPARAMS "$JAVAPARAMS

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--JAVA        --\n" $(java -Xmx200m -version 2>&1)
[ -z "$(which java)" ] && echo "[ERROR] no java detected" && exit 1
echo -e "--igvtools    --\n "$(java $JAVAPARAMS -jar $PATH_IGVTOOLS/igvtools.jar version 2>&1)
[ ! -f $PATH_IGVTOOLS/igvtools.jar ] && echo "[ERROR] no igvtools detected" && exit 1
echo -e "--GATK        --\n "$(java $JAVAPARAMS -jar $PATH_GATK/GenomeAnalysisTK.jar --version 2>&1)
[ ! -f $PATH_GATK/GenomeAnalysisTK.jar ] && echo "[ERROR] no GATK detected" && exit 1


NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "parameters"

FILES=${FILES//$DUMMY/$TARGET}

if [ ! $CPU_VARCOLLECT -eq 1 ]; then
	echo "[NOTE] parallele"
	PARALLELE="-nt $CPU_VARCOLLECT"
fi

# delete old bam file
if [ -z "$NGSANE_RECOVERFROM" ]; then
    if [ -e $OUTDIR/joined.$TARGET ]; then rm $OUTDIR/joined.$TARGET* ; fi
fi

# ensure dir is there
if [ ! -d $OUTDIR ]; then mkdir -p $OUTDIR; fi

# prep for joining
VARIANTS=""
VARFILES=""
NAMES=""
for f in $FILES; do
	[ ! -e "$f" ] && echo "[NOTE] missing file $f (skipped)" && continue 
    i=${f/$STAGEONEINPUT/$STAGETWOINPUT} #point to var folder
    i=$(echo $i | sed 's/bam$/vcf/g') #${i/bam/"vcf"} # correct ending
	VARFILES=$i
    b=$(basename $i)
    arrIN=(${b//./ })
    name=${arrIN[0]}
    NAMES=$NAMES"$name,"
    VARIANTS=$VARIANTS" --variant:$name $i "
done

echo $VARIANTS

#REGION=""
#if [ -n "$REF" ]; then echo $REF; REGION="-L $REF"; fi

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "recall files from tape"

if [ -n "$DMGET" ]; then
	dmget -a $VARFILES
fi
    
NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "join with GATK"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    java $JAVAPARAMS -jar $PATH_GATK/GenomeAnalysisTK.jar -l INFO \
       -R $FASTA \
       -T CombineVariants \
       $VARIANTS \
       -o $OUTDIR/joined.$TARGET \
       -genotypeMergeOptions PRIORITIZE \
       $ADDCOLLECTREGIONPARAM \
   	   $PARALLELE \
       -priority $NAMES 

    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/joined.$TARGET

fi 

################################################################################
NGSANE_CHECKPOINT_INIT "index for IGV"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    java $JAVAPARAMS -jar $PATH_IGVTOOLS/igvtools.jar index $OUTDIR/joined.$TARGET

    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK 
fi

################################################################################
NGSANE_CHECKPOINT_INIT "evaluate"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    echo "[NOTE] eval variants"
    java $JAVAPARAMS -jar $PATH_GATK/GenomeAnalysisTK.jar -l INFO \
            -T VariantEval \
            -R $FASTA \
            --dbsnp $DBSNPVCF \
            --eval $OUTDIR/joined.$TARGET \
	       $ADDCOLLECTREGIONPARAM \
            --evalModule TiTvVariantEvaluator \
            -o $OUTDIR/joined.$TARGET.eval.txt

    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/joined.$TARGET.eval.txt

fi

################################################################################
echo ">>>>> Join variant after calling with any variant caller - FINISHED"
echo ">>>>> enddate "`date`

