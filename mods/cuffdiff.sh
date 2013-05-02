#!/bin/bash

# cufflinks calling script
# author: Denis C. Bauer
# date: March.2011

# messages to look out for -- relevant for the QC.sh script:
# QCVARIABLES,

echo ">>>>> differential expression with cuffdiff"
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> cuffdiff.sh $*"


function usage {
echo -e "usage: $(basename $0) -k HISEQINF -f  -r REFERENCE -o OUTDIR [OPTIONS]

Script running read mapping for single and paired DNA reads from fastq files
It expects a fastq file, pairdend, reference genome  as input and 
It runs BWA, converts the output to .bam files, adds header information and
writes the coverage information for IGV.

required:
  -k | --toolkit <path>     location of the HiSeqInf repository 
  -b | --basename <b1[,b2]> basename comma separated
  -r | --reference <file>   reference genome
  -o | --outdir <path>      output dir

options:
  -t | --threads <nr>       number of CPUs to use (default: 1)
  -i | --rgid <name>        read group identifier RD ID (default: exp)
  -l | --rglb <name>        read group library RD LB (default: qbi)
  -p | --rgpl <name>        read group platform RD PL (default: illumna)
  -s | --rgsi <name>        read group sample RG SM prefac (default: )
  -R | --region <ps>        region of specific interest, e.g. targeted reseq
                             format chr:pos-pos
  -S | --sam                do not convert to bam file (default confert); not the
                             resulting sam file is not duplicate removed
"
exit
}


if [ ! $# -gt 3 ]; then usage ; fi

#DEFAULTS
THREADS=1

#INPUTS
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; HISEQINF=$1 ;; # location of the HiSeqInf repository
        -t | --threads )        shift; THREADS=$1 ;; # number of CPUs to use
        -b | --basename )       shift; fs=$1 ;; # basename
        -r | --reference )      shift; FASTA=$1 ;; # reference genome
        -o | --outdir )         shift; OUT=$1 ;; # output dir
	-a | --annot )          shift; REFSEQGTF=$1 ;; # refseq annotation
        -h | --help )           usage ;;
        * )                     usage
    esac
    shift
done

#PROGRAMS
. $HISEQINF/conf/header.sh


if [ -d $OUT ]; then rm -r $OUT; fi
mkdir $OUT

n=${fs/,/:}
O=${OUT/$n/}
CUFOUT=${O/$TASKCUFFDIFF/Run\/$TASKCUFF/}
TOPHATOUT=${O/$TASKCUFFDIFF/Run\/$TASKTOPHAT/}

CUFGTFS=""
TOPHATBAM=""
for v in ${fs//,/ }; do
    f=$(basename $v)
    CUFGTFS=$CUFGTFS" "$CUFOUT/$f/transcripts.gtf
    TOPHATBAM=$TOPHATBAM" "$TOPHATOUT/$f/accepted_hits.bam
done

cd $OUT/
# Which transcript reference?
TRASCRIPTS=$REFSEQGTF
if [ -n $TRASCRIPTS ]; then
    echo "********* compare to oneanother"
    $CUFFLINKSHOME/cuffcompare -o comp.txt $CUFGTFS
    
    TRASCRIPTS=$OUT/comp.combined.gtf
    
fi

echo "********* run cuff diff"
$CUFFLINKSHOME/cuffdiff -r $FASTA -p $THREADS -o $OUT $TRASCRIPTS $TOPHATBAM

cd ../../

echo ">>>>> differential expression with cuffdiff - FINISHED"
echo ">>>>> enddate "`date`