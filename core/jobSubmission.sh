#!/bin/bash -e

# Script to submit jobs to the queue. 
# NOTE: Not meant to be used outside the framework
#
# author: Denis C. Bauer
# date: Nov.2010

function usage {
echo -e "usage: $(basename $0) -k NGSANE [OPTIONS]"
exit
}

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
    -q | --queue )          shift; QUEUE=$1 ;; # nodetype
    -o | --output )	        shift; SOUTPUT=$1 ;; # pbsoutput
    -a | --additional )	    shift; SADDITIONAL=$1 ;; # additional paramers
    -p | --command )	    shift; SCOMMAND=$1 ;; # Program call
    -t | --tmpdir )	        shift; STMPDIR=$1 ;; # additional paramers	
    -h | --help )           usage ;;
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
echo "#!/bin/bash -e" > $TMPFILE
echo "cd $(pwd)" >> $TMPFILE
echo $SCOMMAND >> $TMPFILE
echo "rm $TMPFILE" >> $TMPFILE

# truncate jobname to 64 characters due to job schedule constraints
SNAME="NGs_${SNAME:0:60}"

# split total memory by CPUs if needed
if [ -n "$QUEUEMEMORYPERSLOT" ]; then
    SMEMORY=$(echo $SMEMORY | sed "s|\(^[0-9. ]*\).*|scale=2; \1 / $SCPU |g" | bc -l)$(echo $SMEMORY | sed "s|^[0-9. ]*\(.*\)|\1|g")
fi   

if [ "$SUBMISSIONSYSTEM" == "PBS" ]; then

	# Prepare Prologue Script
	# Emulate SGE's behaviour of appending to the previous $SOUTPUT.sh by cat in the prologue
    if [ -e $SOUTPUT.sh ]; then rm -f $SOUTPUT.sh; fi
	echo -e "#!/bin/sh \n cat $SOUTPUT \n exit 0" > $SOUTPUT.sh
	chmod 500 $SOUTPUT.sh
	# add to the TMPFILE that the prologue scripts needs deleting
	echo "rm -fr $SOUTPUT.sh" >> $TMPFILE

#	echo "********** submit with PBS submission system" 1>&2
	JOBIDS=$QUEUEWAIT${JOBIDS//:/$QUEUEWAITSEP}
	command="qsub $JOBIDS -V -S /bin/bash -j oe -o $SOUTPUT -w $(pwd) -l $SNODES -l vmem=$SMEMORY \
		-N $SNAME -l walltime=$SWALLTIME $TMPFILE $SADDITIONAL -l prologue=$SOUTPUT.sh	"

	echo "# $command" >> $TMPFILE
	RECIPT=$($command)
    JOBID=$(echo "$RECIPT" | gawk '{print $(NF-1); split($(NF-1),arr,"."); print arr[1]}' | tail -n 1)
	echo $JOBID


elif [ "$SUBMISSIONSYSTEM" == "SGE" ]; then
    ## add resource survey for job to logfile
    echo "echo '>>>>> job resources'" >> $TMPFILE
    echo "qstat -j $SNAME 2>/dev/null | egrep 'vmem|parallel environment'" >> $TMPFILE

#   taking care of the SGE module bug
    unset module
#	echo "********** submit with SGE submission system"
	if [ -n "$JOBIDS" ];then JOBIDS=$(echo -e $JOBIDS | sed 's/^://g' | sed 's/:/,/g'); HOLD_JID="-hold_jid $JOBIDS"; fi
	command="qsub $HOLD_JID -V -S /bin/bash -j y -o $SOUTPUT -cwd -pe $QUEUEPARENV $SCPU -l h_vmem=$SMEMORY \
	    -N $SNAME -l h_rt=$SWALLTIME $SADDITIONAL $TMPFILE" 
	echo "# $command" >>$TMPFILE
	RECIPT=$($command)
	JOBID=$(echo "$RECIPT" | awk '{print $3}')
	echo $JOBID

else
	echo "Submission system, $SUBMISSIONSYSTEM, not implemented; only SGE or PBS are currently supported"
	exit
fi
	
