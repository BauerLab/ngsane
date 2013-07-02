#!/bin/bash

# Call variants with samtools
# author: Denis C. Bauer
# date: Feb.2013

# messages to look out for -- relevant for the QC.sh script:
# 

echo ">>>>> Collect Variants after calling with sam "
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> samSNPscollect.sh $*"


function usage {
echo -e "usage: $(basename $0) -k CONFIG -f BAM -o OUTDIR [OPTIONS]

Variant calling with sam

required:
  -k | --toolkit <path>     config file
  -f <files>                 bam files
  -o <dir>
options:

"
exit
}


if [ ! $# -gt 3 ]; then usage ; fi

#DEFAULTS

#INPUTS
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the NGSANE repository
        -f           )          shift; FILES=${1//,/ } ;; # bam files
	-o           )          shift; MYOUT=$1 ;; # outputdir
        -h | --help )           usage ;;
        * )                     echo "don't understand "$1
    esac
    shift
done

#PROGRAMS
. $CONFIG
. ${NGSANE_BASE}/conf/header.sh
. $CONFIG


JAVAPARAMS="-Xmx"$MEMORY_SAMVAR"G -Djava.io.tmpdir="$TMP #-XX:ConcGCThreads=1 -XX:ParallelGCThreads=1 -XX:MaxDirectMemorySize=10G"
echo "JAVAPARAMS "$JAVAPARAMS

echo "********** programs"
for MODULE in $MODULE_SAMVAR; do module load $MODULE; done  # save way to load modules that itself load other modules
export PATH=$PATH_SAMVAR:$PATH
module list
echo "PATH=$PATH"
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)
PATH_GATK=$(dirname $(which GenomeAnalysisTK.jar))
PATH_IGVTOOLS=$(dirname $(which igvtools.jar))
echo -e "--JAVA    --\n" $(java -version 2>&1)
[ -z "$(which java)" ] && echo "[ERROR] no java detected" && exit 1
echo -e "--igvtools--\n "$(java -jar $JAVAPARAMS $PATH_IGVTOOLS/igvtools.jar version 2>&1)
[ ! -f $PATH_IGVTOOLS/igvtools.jar ] && echo "[ERROR] no igvtools detected" && exit 1
echo -e "--GATK  --\n "$(java -jar $JAVAPARAMS $PATH_GATK/GenomeAnalysisTK.jar --version 2>&1)
[ ! -f $PATH_GATK/GenomeAnalysisTK.jar ] && echo "[ERROR] no GATK detected" && exit 1
# get basename of f
n=${f##*/}

# delete old bam file
#if [ -e $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam} ]; then rm $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}; fi
#if [ -e $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}.stats ]; then rm $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}.stats; fi
#if [ -e $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}.dupl ]; then rm $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}.dupl; fi

# ensure dir is there
if [ ! -d $MYOUT ]; then mkdir -p $MYOUT; fi

if [ -n $DMGET ]; then dmget -a $FILES; fi

# prep for joining
VARIANTS=""
NAMES=""
for f in $FILES; do
    i=${f/$TASKBWA/$TASKBWA"-"$TASKSAMVAR} #point to var folder
    i=${i/bam/"clean.vcf"} # correct ending
    b=$(basename $i)
    arrIN=(${b//./ })
    name=${arrIN[0]}
    NAMES=$NAMES"$name,"
    VARIANTS=$VARIANTS" --variant:$name $i "
done

echo $VARIANTS

REGION=""
if [ -n "$REF" ]; then echo $REF; REGION="-L $REF"; fi


echo "********** join with GATK"
java $JAVAPARAMS -jar $PATH_GATK/GenomeAnalysisTK.jar -l INFO \
   -R $FASTA \
   -T CombineVariants \
   $VARIANTS \
   -o $MYOUT/joined.vcf \
   -genotypeMergeOptions PRIORITIZE \
   $REGION \
   -priority $NAMES

echo "********* make index for IGV"
java $JAVAPARAMS -jar $PATH_IGVTOOLS/igvtools.jar index $MYOUT/joined.vcf


echo ">>>>> Join variant after calling with sam"
echo ">>>>> enddate "`date`

