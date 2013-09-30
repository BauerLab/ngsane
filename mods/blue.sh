#!/bin/bash -e

# Script to ... 
# It takes a <Run>/*.$FASTQ[.gz] file and returns <Run>_healed/*.$FASTQ[.gz]
#
# author: Denis Bauer
# date: Sept 2013

# messages to look out for -- relevant for the QC.sh script:
# QCVARIABLES,
# RESULTFILENAME fastq/<DIR>_blue/<SAMPLE>$READONE.$FASTQ

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


#defaults
KMER=25
GENOME=200000
CUTOFF=10


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
echo -e "--blue        --\n "  $(which Blue.exe) | gawk '{print $NF}'
[ ! -f $(which Blue.exe) ] && echo "[ERROR] no Blue detected" && exit 1

echo -e "\n********* $CHECKPOINT"
################################################################################
CHECKPOINT="parameters"

#OUTDIR=$(dirname $f)"_blue"
# get basename of f
n=${f##*/}

#is paired ?
#if [ "$f" != "${f/$READONE/$READTWO}" ] && [ -e ${f/$READONE/$READTWO} ]; then
#    echo "[NOTE] PAIRED library"
#    PAIRED="1"
#else
#    echo "[NOTE] SINGLE library"
#    PAIRED="0"
#fi

mkdir -p $OUTDIR/tessel

FILES=""
for i in $(ls ${f/$READONE/\*}); do
	if [[ ${f##*.} == "gz" ]]; then
		echo "[NOTE] unzip $i"
		$GZIP -c $i > $i.unzipped
		FILES=$FILES" "$i.unzipped
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

	RUN_COMMAND="mono ${BLUE_HOME}/Tessel.exe $TESSELADDPARAM -t $CPU_BLUE -tmp $TMP -k $KMER -g $GENOME $OUTDIR/tessel/$n $FILES"

    echo $RUN_COMMAND && eval $RUN_COMMAND

    # mark checkpoint
    [ -f $OUTDIR/tessel/$n"_"$KMER.cbt ] && echo -e "\n********* $CHECKPOINT" && unset RECOVERFROM
fi

################################################################################
CHECKPOINT="run blue"    

if [[ -n "$RECOVERFROM" ]] && [[ $(grep "********* $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 


	RUN_COMMAND="mono ${BLUE_HOME}/Blue.exe $BLUEADDPARAM -t $CPU_BLUE -m $CUTOFF -o $OUTDIR $OUTDIR/tessel/$n*.cbt $FILES"

    echo $RUN_COMMAND && eval $RUN_COMMAND

	#get back to NGSANE fastq format
	for i in $(ls $OUTDIR/*corrected*.fastq); do mv $i ${i/_corrected_$CUTOFF/}; done

    # mark checkpoint
    [ -f $OUTDIR/$n ] && echo -e "\n********* $CHECKPOINT" && unset RECOVERFROM
fi


################################################################################
[ -e $OUTDIR/${n}.dummy ] && rm $OUTDIR/${n}.dummy
echo ">>>>> read correction with Blue - FINISHED"
echo ">>>>> enddate "`date`

