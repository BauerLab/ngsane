#!/bin/bash

# Call variants with samtools
# author: Denis C. Bauer
# date: Feb.2013

# messages to look out for -- relevant for the QC.sh script:
# 

echo ">>>>> Variant calling with sam "
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> samSNPs.sh $*"


function usage {
echo -e "usage: $(basename $0) -k CONFIG -f BAM -o OUTDIR [OPTIONS]

Variant calling with sam

required:
  -k | --toolkit <path>     config file
  -f <file>                 bam file

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
        -f | --file    )        shift; f=$1 ;; # bam file
        -o | --output  )        shift; MYOUT=$1 ;; # output directory
        -h | --help    )        usage ;;
        * )                     echo "don't understand "$1
    esac
    shift
done

#PROGRAMS
. $CONFIG
. ${NGSANE_BASE}/conf/header.sh
. $CONFIG

echo "********** programs"
for MODULE in $MODULE_SAMVAR; do module load $MODULE; done  # save way to load modules that itself load other modules
export PATH=$PATH_SAMVAR:$PATH
module list
echo "PATH=$PATH"
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)
PATH_IGVTOOLS=$(dirname $(which igvtools.jar))
echo -e "--JAVA    --\n" $(java -version 2>&1)
[ -z "$(which java)" ] && echo "[ERROR] no java detected" && exit 1
echo -e "--samtools--\n "$(samtools 2>&1 | head -n 3 | tail -n-2)
[ -z "$(which samtools)" ] && echo "[ERROR] no samtools detected" && exit 1
echo -e "--igvtools--\n "$(java -jar $JAVAPARAMS $PATH_IGVTOOLS/igvtools.jar version 2>&1)
[ ! -f $PATH_IGVTOOLS/igvtools.jar ] && echo "[ERROR] no igvtools detected" && exit 1

echo "[NOTE] set java parameters"
JAVAPARAMS="-Xmx"$(python -c "print int($MEMORY_SAMVAR*0.8)")"g -Djava.io.tmpdir="$TMP"  -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1" 
unset _JAVA_OPTIONS
echo "JAVAPARAMS "$JAVAPARAMS

# get basename of f
n=${f##*/}

# delete old bam file
#if [ -e $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam} ]; then rm $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}; fi
#if [ -e $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}.stats ]; then rm $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}.stats; fi
#if [ -e $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}.dupl ]; then rm $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}.dupl; fi

if [ -n $DMGET ]; then dmget -a $f; fi

echo "********* rm duplicate reads $(date)"
samtools rmdup $f $MYOUT/${n/bam/drm.bam}
samtools index $MYOUT/${n/bam/drm.bam}

echo "********* call variants $(date)"
samtools mpileup -uf $FASTA -q1 -D $MYOUT/${n/bam/drm.bam} |  bcftools view -vcg - >$MYOUT/${n/bam/vcf}

echo "********* convert bcf->vcf; index $(date)"
vcfutils.pl varFilter -D1000 -w0 -e0 $MYOUT/${n/bam/vcf}  >$MYOUT/${n/bam/clean.vcf}

echo "********* index vcf file for viewing in IGV $(date)"
java $JAVAPARAMS -jar $PATH_IGVTOOLS/igvtools.jar index $MYOUT/${n/bam/clean.vcf}


rm $MYOUT/${n/bam/drm.bam} $MYOUT/${n/bam/drm.bam}.bai
rm $MYOUT/${n/bam/vcf}


echo ">>>>> Variant calling with sam"
echo ">>>>> enddate "`date`

