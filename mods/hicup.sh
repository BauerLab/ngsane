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
  -i | --rgid <name>        read group identifier RD ID (default: exp)
  -l | --rglb <name>        read group library RD LB (default: qbi)
  -p | --rgpl <name>        read group platform RD PL (default: illumna)
  -s | --rgsi <name>        read group sample RG SM prefac (default: )
  -u | --rgpu <name>        read group platform unit RG PU (default:flowcell )
  -S | --sam                do not convert to bam file (default confert); not the
                             resulting sam file is not duplicate removed
  --noMapping
  --fastqName               name of fastq file ending (fastq.gz)
  --oldIllumina
"
exit
}


if [ ! $# -gt 3 ]; then usage ; fi



#DEFAULTS
MYTHREADS=1
MYMEMORY=2
EXPID="exp"           # read group identifier RD ID
LIBRARY="qbi"         # read group library RD LB
PLATFORM="illumina"   # read group platform RD PL
UNIT="flowcell"       # read group platform unit RG PU
DOBAM=1               # do the bam file
FORCESINGLE=0
NOMAPPING=0
FASTQNAME=""
QUAL="" # standard Sanger


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
	-i | --rgid )           shift; EXPID=$1 ;; # read group identifier RD ID
	-l | --rglb )           shift; LIBRARY=$1 ;; # read group library RD LB
	-p | --rgpl )           shift; PLATFORM=$1 ;; # read group platform RD PL
	-s | --rgsi )           shift; SAMPLEID=$1 ;; # read group sample RG SM (pre)
	-u | --rgpu )           shift; UNIT=$1 ;; # read group platform unit RG PU 
        -S | --sam )            DOBAM=0 ;;
	--fastqName )           shift; FASTQNAME=$1 ;; #(name of fastq or fastq.gz)
	--noMapping )           NOMAPPING=1;;
	--oldIllumina )         QUAL="-S";;   # old illumina encoding 1.3+
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
echo -e "--JAVA    --\n" $(java $JAVAPARAMS -version 2>&1)
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
OUTDIR=${n/'_'$READONE.$FASTQ/}
if [ -d $MYOUT/$OUTDIR ]; then rm -rf $MYOUT/$OUTDIR; fi

# delete old result files 
rm -f $MYOUT/${n/'_'$READONE.$FASTQ/}*.txt

#is paired ?
if [ -e ${f/$READONE/$READTWO} ]; then
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

FULLSAMPLEID=$SAMPLEID"${n/'_'$READONE.$FASTQ/}"
echo ">>>>> full sample ID "$FULLSAMPLEID

echo "********* create hicup conf script"

HICUP_CONF=$MYOUT/${n/'_'$READONE.$FASTQ/.conf}

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
HICUP_CALL="$(which perl) $(which hicup) -c $HICUP_CONF"
echo $HICUP_CALL
cd $MYOUT/$OUTDIR
$($HICUP_CALL)
cp hicup_deduplicater_summary_results_*.txt $MYOUT/${n/'_'$READONE.$FASTQ/}_hicup_deduplicater_summary_results.txt 2>/dev/null
cp hicup_filter_summary_results_*.txt $MYOUT/${n/'_'$READONE.$FASTQ/}_hicup_filter_summary_results.txt 2>/dev/null
cp hicup_mapper_summary_*.txt $MYOUT/${n/'_'$READONE.$FASTQ/}_hicup_mapper_summary.txt 2>/dev/null
cp hicup_truncater_summary_*.txt $MYOUT/${n/'_'$READONE.$FASTQ/}_hicup_truncater_summary.txt 2>/dev/null
ln -f -s $OUTDIR/uniques_${n/.$FASTQ/}_trunc_${n/'_'$READONE.$FASTQ/'_'$READTWO}_trunc.bam $MYOUT/${n/'_'$READONE.$FASTQ/}_uniques.bam

cd $CURDIR

# copy piecharts
RUNSTATS=$OUT/runStats/hicup
mkdir -p $RUNSTATS
cp $MYOUT/$OUTDIR/uniques_*_cis-trans.png $RUNSTATS/${n/'_'$READONE.$FASTQ/}_uniques_cis-trans.png 2>/dev/null
cp $MYOUT/$OUTDIR/*_ditag_classification.png $RUNSTATS/${n/'_'$READONE.$FASTQ/}_ditag_classification.png 2>/dev/null

echo "********* fit-hi-c"
python ${NGSANE_BASE}/tools/hicupCountInteractions.py --verbose --genomeFragmentFile=${DIGESTGENOME} --outputDir=$MYOUT/  $MYOUT/${n/'_'$READONE.$FASTQ/}_uniques.bam
cd $MYOUT
python $(which fit-hi-c.py) --mappabilityThres=2 --fragments=$MYOUT/${n/'_'$READONE.$FASTQ/}_uniques.bam.fragmentLists --interactions=$MYOUT/${n/'_'$READONE.$FASTQ/}_uniques.bam.contactCounts --lib=${n/'_'$READONE.$FASTQ/}
cd $CURDIR

awk 'float($7)<=0.05' $MYOUT/${n/'_'$READONE.$FASTQ/}.spline_pass1.significances.txt | $GZIP > $MYOUT/${n/'_'$READONE.$FASTQ/}.spline_pass1.q05.txt
awk 'float($7)<=0.05' $MYOUT/${n/'_'$READONE.$FASTQ/}.spline_pass2.significances.txt | $GZIP > $MYOUT/${n/'_'$READONE.$FASTQ/}.spline_pass2.q05.txt

$GZIP $MYOUT/${n/'_'$READONE.$FASTQ/}*.significances.txt

echo ">>>>> readmapping with hicup (bowtie) - FINISHED"
echo ">>>>> enddate "`date`

