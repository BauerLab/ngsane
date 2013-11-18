#!/bin/bash -e

# Script to ... 
# It takes a <Run>/*.$FASTQ[.gz] file and returns <Run>_$TASKBLUE/*.$FASTQ[.gz]
#
# author: Denis Bauer
# date: Sept 2013

# messages to look out for -- relevant for the QC.sh script:
# QCVARIABLES,
# RESULTFILENAME fastq/<DIR>_$TASKBLUE/<SAMPLE>$READONE.$FASTQ

echo ">>>>> read correction with Blue"
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"

while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of NGSANE
        -f | --file )           shift; f=$1 ;; # fastq file
        -o | --outdir )         shift; OUTDIR=$1 ;; # output dir                                                     
        --recover-from )        shift; RECOVERFROM=$1 ;; # attempt to recover from log file
        -h | --help )           usage ;;
        * )                     echo "don't understand "$1
    esac
    shift
done

# overwrite defaults
. $CONFIG
. ${NGSANE_BASE}/conf/header.sh
. $CONFIG

################################################################################
CHECKPOINT="programs"

for MODULE in $MODULE_BLUE; do module load $MODULE; done  # save way to load modules that itself load other modules
export PATH=$PATH_BLUE:$PATH;
module list
echo "PATH=$PATH"
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)

BLUE_HOME=$(dirname $(which Blue.exe))

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--MONO        --\n" $(mono --version | head -n 1 )
[ -z "$(which mono)" ] && echo "[ERROR] no mono detected" && exit 1
echo -e "--BLUE        --\n" $(which Blue.exe | gawk '{print $NF}')
[ ! -f $(which Blue.exe) ] && echo "[ERROR] no Blue detected" && exit 1

echo -e "\n********* $CHECKPOINT"
################################################################################
CHECKPOINT="parameters"

#OUTDIR=$(dirname $f)"_blue"
# get basename of f
n=${f##*/}
SAMPLE=${n/$READONE.$FASTQ/}

if [ -z $BLUE_KMER ]; then
    echo "[ERROR] BLUE_KMER not set"; 
    exit 1
fi

if [ -z $GENOMESIZE ]; then
    echo "[ERROR] GENOMESIZE not set"; 
    exit 1
fi

if [ -z ${BLUE_MINREPS} ]; then
    echo "[ERROR] {BLUE_MINREPS} not set"; 
    exit 1
fi

#is paired ?
if [ "$f" != "${f/$READONE/$READTWO}" ] && [ -e ${f/$READONE/$READTWO} ]; then
    PAIRED="1"
else
    PAIRED="0"
fi

TESSELDIR=$OUTDIR/$SAMPLE"_"tessel
mkdir -p $TESSELDIR

FILES=""
THISTMP=$TMP"/"$(whoami)"/"$(echo $OUTDIR | md5sum | cut -d' ' -f1)
mkdir -p $THISTMP

if [[ ${f##*.} != "gz" ]]; then
    FILES="${f/$FASTQ/}"
    if [ $PAIRED = "1" ]; then 
        FILES="$FILES ${f/$READONE.$FASTQ/$READTWO}"
    fi
else
    echo "[NOTE] unzip input"
    zcat $f > $THISTMP/${n/.$FASTQ/.unzipped}
    FILES="$THISTMP/${n/.$FASTQ/.unzipped}"
    if [ $PAIRED = "1" ]; then 
        zcat ${f/$READONE/$READTWO} > $THISTMP/$SAMPLE$READTWO.unzipped
        FILES="$FILES $THISTMP/$SAMPLE$READTWO.unzipped"
    fi
fi

echo "[NOTE] run on $FILES"

echo -e "\n********* $CHECKPOINT"
################################################################################
CHECKPOINT="recall files from tape"

if [ -n "$DMGET" ]; then
    dmget -a ${f/$READONE/"*"}
	dmget -a ${OUTDIR}
fi

echo -e "\n********* $CHECKPOINT"
################################################################################
CHECKPOINT="tessel"    

if [[ -n "$RECOVERFROM" ]] && [[ $(grep "********* $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

	RUN_COMMAND="mono ${BLUE_HOME}/Tessel.exe $TESSELADDPARAM -t ${CPU_BLUE} -tmp $THISTMP -k $BLUE_KMER -g $GENOMESIZE $TESSELDIR/$SAMPLE $FILES"

    echo $RUN_COMMAND && eval $RUN_COMMAND

    # mark checkpoint

    if [ -f $TESSELDIR/$SAMPLE"_"$BLUE_KMER".cbt" ]; then echo -e "\n********* $CHECKPOINT" && unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
fi

################################################################################
CHECKPOINT="blue"    

if [[ -n "$RECOVERFROM" ]] && [[ $(grep "********* $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

	RUN_COMMAND="mono ${BLUE_HOME}/Blue.exe $BLUEADDPARAM -t ${CPU_BLUE} -m ${BLUE_MINREPS} -o $OUTDIR $TESSELDIR/$SAMPLE"*.cbt" $FILES"
    echo $RUN_COMMAND && eval $RUN_COMMAND

    mv $OUTDIR/$SAMPLE$READONE"_corrected_"$BLUE_MINREPS".unzipped" $OUTDIR/$SAMPLE$READONE".fastq"
    if [ $PAIRED = "1" ]; then 
        mv $OUTDIR/$SAMPLE$READTWO"_corrected_"$BLUE_MINREPS".unzipped" $OUTDIR/$SAMPLE$READTWO".fastq"
    fi
	
	# zip
	$GZIP -c $OUTDIR/$SAMPLE$READONE".fastq" > $OUTDIR/$n
	[ -f $OUTDIR/$SAMPLE$READONE".fastq" ] && rm $OUTDIR/$SAMPLE$READONE".fastq"
	if [ $PAIRED = "1" ]; then 
    	$GZIP -c $OUTDIR/$SAMPLE$READTWO".fastq" > $OUTDIR/${n/$READONE.$FASTQ/$READTWO.$FASTQ}
		[ -f $OUTDIR/$SAMPLE$READTWO".fastq" ] && rm $OUTDIR/$SAMPLE$READTWO".fastq"
    fi

    if [ -n "$READONE" ]; then
        mv $OUTDIR/$SAMPLE$READONE"_corrected_"$BLUE_MINREPS"_stats.txt" $OUTDIR/$SAMPLE"_corrected_"$BLUE_MINREPS"_stats.txt"
    fi

    # mark checkpoint
    if [ -f $OUTDIR/$n ]; then echo -e "\n********* $CHECKPOINT" && unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
fi

################################################################################
CHECKPOINT="cleanup"    

for i in $(ls ${f/$READONE/\*}); do
	if [[ ${f##*.} == "gz" ]]; then
		[ -e $THISTMP/$n.unzipped ] && rm $THISTMP/$n.unzipped
	fi
done

if [ -d $TESSELDIR ]; then
    rm -r $TESSELDIR
fi

echo -e "\n********* $CHECKPOINT"
################################################################################
[ -e $OUTDIR/${n}.dummy ] && rm $OUTDIR/${n}.dummy
echo ">>>>> read correction with Blue - FINISHED"
echo ">>>>> enddate "`date`

