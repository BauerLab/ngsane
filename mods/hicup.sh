#!/bin/bash -e

# HiCUP calling script
# author: Fabian Buske
# date: Apr 2013

# messages to look out for -- relevant for the QC.sh script:
# QCVARIABLES

echo ">>>>> HiC readmapping with HiCUP "
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> hicup.sh $*"


function usage {
echo -e "usage: $(basename $0) -k NGSANE -f FASTQ -r REFERENCE -o OUTDIR [OPTIONS]

Script running hicup including reference genome digestion, read mapping for single and 
paired DNA reads with bowtie from fastq files
It expects a fastq file, pairdend, reference genome and digest pattern  as input.

required:
  -k | --toolkit <path>     location of the NGSANE repository 
  -f | --fastq <file>       fastq file
  -r | --reference <file>   reference genome
  -d | --digest <cutsite>   enzyme cutsite pattern, e.g. A^AGCTT,HindIII seperate 2 patterns by ;
  -o | --outdir <path>      output dir

options:
  -t | --threads <nr>       number of CPUs to use (default: 1)
  -m | --memory <nr>        memory available (default: 2)
  --fastqName               name of fastq file ending (fastq.gz)
"
exit
}


if [ ! $# -gt 3 ]; then usage ; fi



#DEFAULTS
MYTHREADS=1
MYMEMORY=2
FASTQNAME=""

#INPUTS
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the NGSANE repository
        -t | --threads )        shift; MYTHREADS=$1 ;; # number of CPUs to use
        -m | --memory )         shift; MYMEMORY=$1 ;; # memory used
        -f | --fastq )          shift; f=$1 ;; # fastq file
        -r | --reference )      shift; FASTA=$1 ;; # reference genome
        -d | --digest )         shift; DIGEST=$1 ;; # digestion patterns
        -o | --outdir )         shift; MYOUT=$1 ;; # output dir
        --fastqName )           shift; FASTQNAME=$1 ;; #(name of fastq or fastq.gz)
        -h | --help )           usage ;;
        * )                     echo "don't understand "$1
    esac
    shift
done

#PROGRAMS
. $CONFIG
. ${NGSANE_BASE}/conf/header.sh
. $CONFIG

JAVAPARAMS="-Xmx"$MYMEMORY"g -Djava.io.tmpdir="$TMP # -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1 -XX:MaxDirectMemorySize=4G"
echo "JAVAPARAMS "$JAVAPARAMS

echo "********** programs"
for MODULE in $MODULE_HICUP; do module load $MODULE; done  # save way to load modules that itself load other modules
export PATH=$PATH_HICUP:$PATH
module list
echo "PATH=$PATH"
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)
echo -e "--JAVA    --\n" $(java -version 2>&1)
[ -z "$(which java)" ] && echo "[ERROR] no java detected" && exit 1
echo -e "--bowtie  --\n "$(bowtie --version | head -n 1 )
[ -z "$(which bowtie)" ] && echo "[ERROR] no bowtie detected" && exit 1
echo -e "--perl    --\n "$(perl -v | grep "version" )
[ -z "$(which perl)" ] && echo "[ERROR] no perl detected" && exit 1
echo -e "--HiCUP   --\n "$(hicup --version )
[ -z "$(which hicup)" ] && echo "[ERROR] no hicup detected" && exit 1
echo -e "--fit-hi-c--\n "$(fit-hi-c.py --version)
[ -z "$(which fit-hi-c.py)" ] && echo "[ERROR] no fit-hi-c detected" && exit 1

# get basename of f
n=${f##*/}


#output for this library
OUTDIR=${n/%$READONE.$FASTQ/}
if [ -d $MYOUT/$OUTDIR ]; then rm -rf $MYOUT/$OUTDIR; fi

# delete old result files 
rm -f $MYOUT/${n/%$READONE.$FASTQ/}*.txt

#is paired ?
if [ "$f" != "${f/$READONE/$READTWO}" ] && [ -e ${f/$READONE/$READTWO} ]; then
    PAIRED="1"
else
    echo "HiCUP requires pair end fastq files"
    exit 1
fi

if [ -n "$DMGET" ]; then
	echo "********** reacall files from tape"
	dmget -a $(dirname $FASTA)/*
	dmget -a $(dirname $(which samtools))/*
	dmget -a $(dirname $(which hicup))/*
	dmget -a ${f/$READONE/"*"}
fi

#is ziped ?
ZCAT="zcat"
if [[ $f != *.gz ]]; then ZCAT="cat"; fi

echo "********** digest reference"
FASTASUFFIX=${FASTA##*.}
FASTABASE=${FASTA##*/}

if [ -z "$DIGEST" ]; then
   echo "[ERROR] No restriction enzyme given!"
   exit 1
fi

ENZYMES=(${DIGEST//;/ })
ENZYME1=(${ENZYMES[0]//,/ })
ENZYME2=(${ENZYMES[1]//,/ })

DIGESTGENOME=""

mkdir -p $MYOUT/$OUTDIR
cd $MYOUT/$OUTDIR
if [ ${#ENZYMES[@]} = 1 ]; then
   echo "Restriction Enzyme 1: ${ENZYME1[1]}:${ENZYME1[0]} "
   DIGESTGENOME=$MYOUT/${FASTABASE/.$FASTASUFFIX/}_${ENZYME1[1]}_None.txt
   hicup_digester -g "${FASTABASE%.*}" -1 ${ENZYME1[0]} $FASTA
   mv Digest_* ${DIGESTGENOME}

elif [ ${#ENZYMES[@]} = 2 ] && [ ! -e $MYOUT/${FASTABASE/.$FASTASUFFIX/}_${ENZYME1[1]}_${ENZYME2[2]}.txt ]; then
   echo "Restriction Enzyme 1: ${ENZYME1[1]}:${ENZYME1[0]} "
   echo "Restriction Enzyme 2: ${ENZYME2[1]}:${ENZYME2[0]} "
   DIGESTGENOME=$MYOUT/${FASTABASE/.$FASTASUFFIX/}_${ENZYME1[1]}_${ENZYME2[2]}.txt
   hicup_digester -g "${FASTABASE%.*}" -1 ${ENZYME1[0]} -2 ${ENZYME2[0]} $FASTA
   mv Digest_* ${DIGESTGENOME}
else
   echo "[ERROR] Invalid number or pattern of enzyme digest patterns."
   exit 1
fi
cd $SOURCE

FULLSAMPLEID=$SAMPLEID"${n/%$READONE.$FASTQ/}"
echo ">>>>> full sample ID "$FULLSAMPLEID

echo "********* create hicup conf script"

HICUP_CONF=$MYOUT/${n/%$READONE.$FASTQ/.conf}

cat /dev/null > $HICUP_CONF
echo "#Number of threads to use" >> $HICUP_CONF
echo "Threads: $MYTHREADS" >> $HICUP_CONF
echo "#Suppress progress updates | 0: off, 1: on" >> $HICUP_CONF
echo "Quiet:0" >> $HICUP_CONF
echo "#Retain all intermediate pipeline files | 0: off, 1: on" >> $HICUP_CONF
echo "Keep:1" >> $HICUP_CONF
echo "#Compress outputfiles | 0: off, 1: on" >> $HICUP_CONF
echo "Zip:1" >> $HICUP_CONF
echo "#Path to the alignment program Bowtie | include the executable Bowtie filename" >> $HICUP_CONF
echo "Bowtie:$(which bowtie)" >> $HICUP_CONF
echo "#Path to the reference genome indices" >> $HICUP_CONF
echo "Index:${BOWTIE_INDEX}"  >> $HICUP_CONF
echo "#Path to the genome digest file" >> $HICUP_CONF
echo "DIGEST:$DIGESTGENOME" >> $HICUP_CONF
echo "#FASTQ file format | phred33-quals, phred64-quals, solexa-quals or solexa1.3-quals" >> $HICUP_CONF
echo "Format:phred33-quals" >> $HICUP_CONF
echo "#Maximum di-tag length | optional parameter" >> $HICUP_CONF
echo "#Longest:" >> $HICUP_CONF
echo "#Minimum di-tag length | optional parameter" >> $HICUP_CONF
echo "#Shortest:" >> $HICUP_CONF
echo "#FASTQ files to be analysed, separating file pairs using the pipe '|' character" >> $HICUP_CONF
echo "$f | ${f/$READONE/$READTWO} " >> $HICUP_CONF

echo "********* execute hicup"
CURDIR=$(pwd)
RUN_COMMAND="$(which perl) $(which hicup) -c $HICUP_CONF"
echo $RUN_COMMAND
cd $MYOUT/$OUTDIR
eval $RUN_COMMAND
cp -f $MYOUT/$OUTDIR/hicup_deduplicater_summary_results_*.txt $MYOUT/${n/%$READONE.$FASTQ/}_hicup_deduplicater_summary_results.txt
cp -f $MYOUT/$OUTDIR/hicup_filter_summary_results_*.txt $MYOUT/${n/%$READONE.$FASTQ/}_hicup_filter_summary_results.txt
cp -f $MYOUT/$OUTDIR/hicup_mapper_summary_*.txt $MYOUT/${n/%$READONE.$FASTQ/}_hicup_mapper_summary.txt
cp -f $MYOUT/$OUTDIR/hicup_truncater_summary_*.txt $MYOUT/${n/%$READONE.$FASTQ/}_hicup_truncater_summary.txt
ln -f -s $OUTDIR/uniques_${n/.$FASTQ/}_trunc_${n/%$READONE.$FASTQ/$READTWO}_trunc.bam $MYOUT/${n/%$READONE.$FASTQ/}_uniques.bam

cd $CURDIR

# copy piecharts
RUNSTATS=$OUT/runStats/$TASKHICUP
mkdir -p $RUNSTATS
cp -f $MYOUT/$OUTDIR/uniques_*_cis-trans.png $RUNSTATS/${n/%$READONE.$FASTQ/}_uniques_cis-trans.png
cp -f $MYOUT/$OUTDIR/*_ditag_classification.png $RUNSTATS/${n/%$READONE.$FASTQ/}_ditag_classification.png

echo "********* fit-hi-c"
python ${NGSANE_BASE}/tools/hicupCountInteractions.py --verbose --genomeFragmentFile=${DIGESTGENOME} --outputDir=$MYOUT/  $MYOUT/${n/%$READONE.$FASTQ/}_uniques.bam
cd $MYOUT
python $(which fit-hi-c.py) --mappabilityThres=2 --fragments=$MYOUT/${n/%$READONE.$FASTQ/}_uniques.bam.fragmentLists --interactions=$MYOUT/${n/%$READONE.$FASTQ/}_uniques.bam.contactCounts --lib=${n/%$READONE.$FASTQ/}
cd $CURDIR

awk '$7<=0.05' $MYOUT/${n/%$READONE.$FASTQ/}.spline_pass1.significances.txt | sort -k7g > $MYOUT/${n/%$READONE.$FASTQ/}.spline_pass1.q05.txt
awk '$7<=0.05' $MYOUT/${n/%$READONE.$FASTQ/}.spline_pass2.significances.txt | sort -k7g > $MYOUT/${n/%$READONE.$FASTQ/}.spline_pass2.q05.txt

$GZIP $MYOUT/${n/%$READONE.$FASTQ/}*.significances.txt

echo ">>>>> readmapping with hicup (bowtie) - FINISHED"
echo ">>>>> enddate "`date`

