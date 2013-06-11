#!/bin/bash -e

# Submission script for PBS and SGE
while [ "$1" != "" ]; do
    case $1 in
    -k | --toolkit )        shift; CONFIG=$1 ;; # location of the NGSANE repository
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
. ${NGSANE_BASE}/conf/header.sh
. $CONFIG

#echo "********** write TMP file"
if [ ! -n "$STMPDIR" ]; then STMPDIR="tmp"; fi
if [ ! -e $STMPDIR ]; then mkdir $STMPDIR; fi
TMPFILE=$STMPDIR/"qsub_"$(date '+%y%m%d%H%m')"_"$(cat /dev/urandom | tr -dc A-Za-z0-9 | head -c 9)".tmp"
echo "cd $(pwd)" > $TMPFILE
echo $SCOMMAND >> $TMPFILE
echo "rm $TMPFILE" >>$TMPFILE

if [ "$SUBMISSIONSYSTEM" == "PBS" ]; then
#	echo "********** submit with PBS submission system"
	JOBIDS=$(echo -e $JOBIDS | sed 's/://')
	command="qsub -W after:$JOBIDS -V -j oe -o $SOUTPUT -w $(pwd) -l $SNODES -l vmem=$SMEMORY \
		-N $SNAME -l walltime=$SWALLTIME $QSUBEXTRA $TMPFILE"
	echo "# $command" >>$TMPFILE
	RECIPT=$($command)
        JOBID=$(echo "$RECIPT" | gawk '{print $(NF-1); split($(NF-1),arr,"."); print arr[1]}' | tail -n 1)
	echo $JOBID

elif [ "$SUBMISSIONSYSTEM" == "SGE" ]; then
#	echo "********** submit with SGE submission system"
	if [ -n "$JOBIDS" ];then JOBIDS=$(echo -e $JOBIDS | sed 's/^://g' | sed 's/:/,/g'); HOLD_JID="-hold_jid $JOBIDS"; fi
	command="qsub $HOLD_JID -V -S /bin/bash -j y -o $SOUTPUT -cwd -pe smp $SCPU -l h_vmem=$SMEMORY \
	    -N $SNAME -l h_rt=$SWALLTIME $QSUBEXTRA $TMPFILE"
	echo "# $command" >>$TMPFILE
	RECIPT=$($command)
	JOBID=$(echo "$RECIPT" | awk '{print $3}')
	echo $JOBID
else
	echo "Submission system, $SUBMISSIONSYSTEM, not implemented; only SGE or PBS work"
	exit
fi
	
