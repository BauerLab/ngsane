#!/bin/bash -e

# Submission script for PBS and SGE


while [ "$1" != "" ]; do
    case $1 in
    -k | --toolkit )        shift; CONFIG=$1 ;; # location of the HiSeqInf repository
    -m | --memory )         shift; SMEMORY=$1 ;; # memory used
    -n | --nodes )          shift; SNODES=$1 ;; # memory used
    -c | --cpu    )         shift; SCPU=$1 ;; # CPU used
    -w | --walltime )       shift; SWALLTIME=$1 ;; # walltime used
    -W | --wait )           shift; JOBIDS=$1 ;; # jobids to wait for
    -j | --jobname   )      shift; SNAME=$1 ;; # name used
    -o | --output )	    shift; SOUTPUT=$1 ;; # pbsoutput
    -a | --additional )	    shift; SADDITIONAL=$1 ;; # additional paramers
    -p | --command )	    shift; SCOMMAND=$1 ;; # Program call
    -t | --tmpdir )	    shift; STMPDIR=$1 ;; # additional paramers	
#	-h | --help )           usage ;;
    * )                     echo "don't understand "$1
    esac
    shift
done

#PROGRAMS
. $CONFIG
. $HISEQINF/pbsTemp/header.sh
. $CONFIG

echo "********** write TMP file"
if [ ! -n "$STMPDIR" ]; then STMPDIR="tmp"; fi
if [ ! -e $STMPDIR ]; then mkdir $STMPDIR; fi
TMPFILE=$STMPDIR/$(cat /dev/urandom | tr -dc A-Za-z0-9 | head -c 9)".tmp"
echo "cd $(pwd)" > $TMPFILE
echo $SCOMMAND >> $TMPFILE
echo "rm $TMPFILE" >>$TMPFILE

if [ "$SUBMISSIONSYSTEM"="PBS" ]; then
#	echo "********** submit with PBS submission system"
	qsub -W after:$JOBIDS -j oe -o $SOUTPUT -w $(pwd) -l $SNODES -l vmem=$SMEMORY \
		-N $SNAME -l walltime=$SWALLTIME $QSUBEXTRA $TMPFILE
elif [ "$SUBMISSIONSYSTEM"="SGE" ]; then
#	echo "********** submit with SGE submission system"
	qsub -hold_jid $JOBIDS -v -S /bin/bash -j y -o $SOUTPUT -cwd -pe smp $SCPU -l h_vmem=$SMEMORY \
		-N $SNAME -l h_rt=$SWALLTIME $QSUBEXTRA $TMPFILE
else
	echo "Submission system, $SUBMISSIONSYSTEM, not implemented; only SGE or PBS work"
	exit
fi
	