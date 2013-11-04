#!/bin/bash -e

# Script to ... 
# It takes a <Run>/*.$FASTQ[.gz] file and returns <Run>_$TASKBLUE/*.$FASTQ[.gz]
#
# author: Denis Bauer
# date: Sept 2013

# messages to look out for -- relevant for the QC.sh script:
# QCVARIABLES,
# RESULTFILENAME fastq/<DIR>_$TASKBLUE/<SAMPLE>$READONE.$FASTQ

echo ">>>>> read screening with FASTQSCREEN"
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
echo -e "--MONO        --\n" $(mono --version 2>&1) | head -n 1
[ -z "$(which mono)" ] && echo "[ERROR] no mono detected" && exit 1
echo -e "--BLUE        --\n" $(which Blue.exe) | gawk '{print $NF}'
[ ! -f $(which Blue.exe) ] && echo "[ERROR] no Blue detected" && exit 1

echo -e "\n********* $CHECKPOINT"
################################################################################
CHECKPOINT="parameters"

#OUTDIR=$(dirname $f)"_blue"
# get basename of f
n=${f##*/}

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
if [ "$f" != "${f/$READONE/$READTWO}" ] && [ -e ${f/$READONE/$READTWO} ] && [ "$FORCESINGLE" = 0 ]; then
    PAIRED="1"
    echo 
else
    PAIRED="0"
fi

mkdir -p $OUTDIR/tessel

FILES=""
THISTMP=$TMP"/"$(whoami)"/"$(echo $OUTDIR | md5sum | cut -d' ' -f1)
echo $THISTMP/blue
mkdir -p $THISTMP
for i in $(ls ${f/$READONE/\*}); do
	if [[ ${f##*.} == "gz" ]]; then
		echo "[NOTE] unzip $i"
		zcat $i > $THISTMP/$n.unzipped
		FILES=$FILES" "$THISTMP/$n.unzipped
	else
		FILES=$FILES" "$i
	fi
done
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

	RUN_COMMAND="mono ${BLUE_HOME}/Tessel.exe $TESSELADDPARAM -t ${CPU_BLUE} -tmp $THISTMP -k $BLUE_KMER -g $GENOMESIZE $OUTDIR/tessel/$n $FILES"

    echo $RUN_COMMAND && eval $RUN_COMMAND

    # mark checkpoint
    [ -f $OUTDIR/tessel/$n"_"$BLUE_KMER.cbt ] && echo -e "\n********* $CHECKPOINT" && unset RECOVERFROM
fi

################################################################################
CHECKPOINT="run blue"    

if [[ -n "$RECOVERFROM" ]] && [[ $(grep "********* $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 


	RUN_COMMAND="mono ${BLUE_HOME}/Blue.exe $BLUEADDPARAM -t ${CPU_BLUE} -m ${BLUE_MINREPS} -o $OUTDIR $OUTDIR/tessel/$n*.cbt $FILES"
    echo $RUN_COMMAND && eval $RUN_COMMAND

	#get back to NGSANE fastq format
	for i in $(ls $OUTDIR/*corrected*.fastq); do 
	   mv $i ${i/_corrected_${BLUE_MINREPS}/}; 
	done
	
	# zip
	$GZIP $OUTDIR/*fastq

    # rename file suffix
    for i in $(ls ${f/$READONE/\*}); do
    	mv ${n/$FASTQ/fastq} $OUTDIR/$n
	done

    # mark checkpoint
    [ -f $OUTDIR/$n ] && echo -e "\n********* $CHECKPOINT" && unset RECOVERFROM
fi

################################################################################
CHECKPOINT="cleanup"    

[ -f $THISTMP ] && rm -r $THISTMP

echo -e "\n********* $CHECKPOINT"
################################################################################
[ -e $OUTDIR/${n}.dummy ] && rm $OUTDIR/${n}.dummy
echo ">>>>> read correction with Blue - FINISHED"
echo ">>>>> enddate "`date`

