#!/bin/bash -e

#
# General template for submitting a job 
#

#INPUTS
while [ "$1" != "" ]; do
    case $1 in
    -q | --queue )          shift; QUEUE=$1 ;;     # nodetype for qsub
	-k | --toolkit )        shift; CONFIG=$1 ;;    # location of the NGSANE repository
	-i | --input )          shift; ORIGIN=$1 ;;    # subfile in $SOURCE
	-e | --fileending )     shift; ENDING=$1 ;;    # select source files of a specific ending
	-t | --task )           shift; TASK=$1 ;;      # what to do
	-n | --nodes )          shift; NODES=$1;;
	-c | --cpu    )         shift; CPU=$1 ;;       # CPU used
	-m | --memory )         shift; MEMORY=$1;;     # min Memory required
	-w | --walltime )       shift; WALLTIME=$1;;
    -W | --waitfor )        shift; WAITFOR=$1 ;;    # wait for previous TASK to finish (asumes job names follows std. NGSane rules)
	-p | --command )        shift; COMMAND=$1;;
	--postcommand )         shift; POSTCOMMAND=$1;;
	--postnodes )           shift; POSTNODES=$1;;
	--postcpu )             shift; POSTCPU=$1 ;;   # CPU used for postcommand
	--postmemory )          shift; POSTMEMORY=$1;; # Memory used for postcommand
	--postwalltime )        shift; POSTWALLTIME=$1;;
	-r | --reverse )        REV="1";;              # input is fastq
	-d | --nodir )          NODIR="nodir";;
	-a | --armed )          ARMED="armed";;
    -W | --wait )           shift; JOBIDS=$1 ;;    # jobids to wait for
	--force )               FORCE="TRUE";;         # don't double check
	--keep )                KEEP="keep";;
	--new )                 KEEP="new";;
	--recover )             RECOVER="recover";;
	--direct )              DIRECT="direct";;
	--first )               FIRST="first";;
	--postonly )            POSTONLY="postonly" ;;
	--dryrun )              DRYRUN="TRUE" ;;
	-h | --help )           usage ;;
	* )                     echo "prepareJobSubmission.sh: don't understand "$1
    esac
    shift
done

#PROGRAMS
. $CONFIG
. ${NGSANE_BASE}/conf/header.sh
. $CONFIG

echo -e "\e[96m[Task]\e[0m $TASK $NODIR"

if [ ! -d $QOUT/$TASK ]; then mkdir -p $QOUT/$TASK; fi

## Select files in dir to run 
if [[ ! -e $QOUT/$TASK/runnow.tmp || "$DIRECT" || "$KEEP" ]]; then
    echo -e "[NOTE] setup enviroment"
    if [ -e $QOUT/$TASK/runnow.tmp ]; then rm $QOUT/$TASK/runnow.tmp; fi
    
    for dir in ${DIR[@]}; do
        
        #ensure dirs are there
        if [ -z "$NODIR" ]; then
            if [ ! -d $OUT/$dir/$TASK ]; then mkdir -p $OUT/$dir/$TASK; fi
        fi
        
        # add tasks to runnow.tmp
        # search for real files and dummy files, in case both exist only keep real one
        if [ -n "$REV" ]; then
            for f in $( ls $SOURCE/$dir/$ORIGIN/*$ENDING* | grep -P ".$ENDING(.dummy)?\$" | sed 's/.dummy//' | sort -u ); do
                n=${f##*/}
                if [ "$KEEP" = "new" ]; then
                    # check if file has been processed previousely
                	COMMANDARR=(${COMMAND// / })
                	DUMMY="echo "$(grep -P "^# *RESULTFILENAME" ${COMMANDARR[0]} | cut -d " " -f 3- | sed sed "s/<SAMPLE>/$name/" | sed "s/<DIR>/$dir/" | sed "s/<TASK>/$TASK/")
                    D=$(eval $DUMMY)
                	[ -n "$D" ] && [ -f $TASK/$dir/$D ]  && echo -e "\e[34m[SKIP]\e[0m $dir/$D (already processed)" && continue
                fi 
                echo -e "\e[32m[TODO]\e[0m $dir/$n"
                echo $f >> $QOUT/$TASK/runnow.tmp
            done
        
        else
            for f in $( ls $SOURCE/$ORIGIN/$dir/*$ENDING* | grep -P ".$ENDING(.dummy)?\$" | sed 's/.dummy//' | sort -u ); do
                n=${f##*/}
                if [ "$KEEP" = "new" ]; then
                    # check if file has been processed previousely
                	COMMANDARR=(${COMMAND// / })
                	DUMMY="echo "$(grep -P "^# *RESULTFILENAME" ${COMMANDARR[0]} | cut -d " " -f 3- | sed "s/<SAMPLE>/$name/" | sed "s/<DIR>/$dir/" | sed "s/<TASK>/$TASK/")
                    D=$(eval $DUMMY)
                	[ -n "$D" ] && [ -f $dir/$TASK/$D ]  && echo -e "\e[34m[SKIP]\e[0m $dir/$D (already processed)" && continue
                fi 
                echo -e "\e[32m[TODO]\e[0m $dir/$n"
                echo $f >> $QOUT/$TASK/runnow.tmp
            done
        fi
    done
else
    echo -e "[NOTE] previous enviroment setup detected"
fi

if [ -n "$KEEP" ]; then 
    if [ "$KEEP" = "new" ]; then
        echo -e "[NOTE] Data setup finished. Please, start trigger in \e[4marmed\e[24m or \e[4mdirect\e[24m mode."
    else
        echo -e "[NOTE] Data setup finished. Please, inspect/modify runnow.tmp in the qout/TASK folder and then start \e[4marmed\e[24m or \e[4mdirect\e[24m mode."
    fi
    exit ; 
else
    echo -e "[NOTE] proceeding with job scheduling..."
fi

if [[ "$FORCE" != "TRUE" && "$DRYRUN" != "TRUE" ]]; then
    echo -n -e "Double check! Then type \e[4msafetyoff\e[24m and hit enter to launch the job: "
    read safetyoff
    if [ "$safetyoff" != "safetyoff" ];then
        echo -e "Holstering..."
        exit 0
    else
        echo -e "... take cover!"
    fi
fi

MYPBSIDS="" # collect job IDs for postcommand
DIR=""
FILES=""
for i in $(cat $QOUT/$TASK/runnow.tmp); do

    n=$(basename $i) 
    # ending : fastq/xx or xx/bwa
    dir=$(dirname $i | gawk '{n=split($1,arr,"/"); print arr[n]}')
    if [ -n "$REV" ]; then dir=$(dirname $i | gawk '{n=split($1,arr,"/"); print arr[n-1]}'); fi
    name=${n/$ENDING/}
    echo -e "[NOTE] "$dir"/"$name
                
    COMMAND2=${COMMAND//<FILE>/$i} # insert files for which parallele jobs are submitted
    COMMAND2=${COMMAND2//<DIR>/$dir} # insert output dir
    COMMAND2=${COMMAND2//<NAME>/$name} # insert ??

    LOGFILE=$QOUT/$TASK/$dir'_'$name'.out'
    DIR=$DIR" $dir"
    FILES=$FILES" $i"

    # only postcommand 
    if [[ -n "$POSTONLY" || -z "$COMMAND" ]]; then continue; fi

	# create dummy files for the pipe
	COMMANDARR=(${COMMAND// / })
	DUMMY="echo "$(grep -P "^# *RESULTFILENAME" ${COMMANDARR[0]} | cut -d " " -f 3- | sed "s/<SAMPLE>/$name/" | sed "s/<DIR>/$dir/" | sed "s/<TASK>/$TASK/")
	D=$(eval $DUMMY)
	echo "[NOTE] make $D.dummy"
	touch $D.dummy
	ls $D.dummy
	
#	if [ -z "$NODIR" ]; then
#		touch $dir/$TASK/$D.dummy # normal case
#	else
#		if [ -n "$REV" ]; then
#			touch $ORIGIN/$dir"_"$TASK/$D.dummy # no dir
#		else
#			touch $dir/$ORIGIN/$D.dummy # no dir
#		fi
#    fi

    echo -e "\e[97m[JOB]\e[0m  $COMMAND2"

    if [ -n "$DIRECT" ]; then eval $COMMAND2; fi

    if [ -n "$ARMED" ]; then

        echo $ARMED

        if [ -n "$RECOVER" ] && [ -f $LOGFILE ] ; then
            # add log-file for recovery
            COMMAND2="$COMMAND2 --recover-from $LOGFILE"
            
            if [[ $(grep -P "^>{5} .* FINISHED" $LOGFILE | wc -l ) -gt 0 ]] ; then
                echo "[NOTE] Previous $TASK run finished without error - nothing to be done"
                MYPBSIDS=""
                continue
            else
                echo "[NOTE] #########################################################################" >> $LOGFILE
                echo "[NOTE] Recover from logfile: $LOGFILE" >> $LOGFILE
                # mask old errors
                sed -i "s/^\[ERROR\] /[NOTE][PREVIOUS][ERROR] /g" $LOGFILE
                echo "[NOTE] Previous errors masked" >> $LOGFILE
            fi
            
        else
            # remove old submission output logs
            if [ -e $QOUT/$TASK/$dir'_'$name.out ]; then rm $QOUT/$TASK/$dir'_'$name.out; fi
        fi
    
        # record task in log file
        cat $CONFIG ${NGSANE_BASE}/conf/header.sh > $QOUT/$TASK/job.$(date "+%Y%m%d").log
        echo "[NOTE] Jobfile: "$QOUT/$TASK/job.$(date "+%Y%m%d").log >> $LOGFILE

        if [ -n "JOBIDS" ]; then 
            RECIPT=$($BINQSUB -a "$QSUBEXTRA" -W "$JOBIDS" -k $CONFIG -m $MEMORY -n $NODES -c $CPU -w $WALLTIME \
        	   -j $TASK'_'$dir'_'$name -o $LOGFILE --command "$COMMAND2")
        else 
            RECIPT=$($BINQSUB -a "$QSUBEXTRA" -k $CONFIG -m $MEMORY -n $NODES -c $CPU -w $WALLTIME \
        	   -j $TASK'_'$dir'_'$name -o $LOGFILE --command "$COMMAND2")
        fi    	

### Trinity additions - start
    # Add previous task dependencies (need another parameter for array dependencies)
    # .....assumes standard NGsane job name rules
#    if [ -n $WAITFOR ]; then 
#        DEPENDENCIES=${WAITFOR}"_"$dir"_"$name
#        RECIPT=$($BINQSUB -a "$QSUBEXTRA" -k $CONFIG -m $MEMORY -n $NODES -c $CPU -w $WALLTIME -W $DEPENDENCIES \
#            -j $TASK'_'$dir'_'$name -o $LOGFILE --command "$COMMAND2")
#    else
#        RECIPT=$($BINQSUB -a "$QSUBEXTRA" -k $CONFIG -m $MEMORY -n $NODES -c $CPU -w $WALLTIME \
#            -j $TASK'_'$dir'_'$name -o $LOGFILE --command "$COMMAND2")
#    fi

### Trinity additions - end

        echo -e "Jobnumber $RECIPT"
        MYPBSIDS=$MYPBSIDS":"$RECIPT
    
        # if only the first task should be submitted as test
        if [ -n "$FIRST" ]; then exit; fi
    
    fi
done


if [ -n "$POSTCOMMAND" ]; then
   # process for postcommand
    DIR=$(echo -e ${DIR// /\\n} | sort -u | gawk 'BEGIN{x=""};{x=x"_"$0}END{print x}' | sed 's/__//')
    FILES=$(echo -e $FILES | sed 's/ /,/g')
    POSTCOMMAND2=${POSTCOMMAND//<FILE>/$FILES}
    POSTCOMMAND2=${POSTCOMMAND2//<DIR>/$DIR}

    echo "[NOTE] "$DIR" wait for "$MYPBSIDS
    echo $POSTCOMMAND2

    if [[ -n "$DIRECT" || -n "$FIRST" ]]; then eval $POSTCOMMAND2; exit; fi
    if [[ -n "$ARMED" ||  -n "$POSTONLY" ]]; then

    # remove old submission output logs
    if [ -e $QOUT/$TASK/postcommand.out ]; then rm $QOUT/$TASK/postcommand.out; fi

    # record task in log file
    cat $CONFIG ${NGSANE_BASE}/conf/header.sh > $QOUT/$TASK/job.$(date "+%Y%m%d").log

    # unless specified otherwise use HPC parameter from main job 
    if [ -z "$POSTNODES" ];    then POSTNODES=$NODES; fi
    if [ -z "$POSTCPU" ];      then POSTCPU=$CPU; fi
    if [ -z "$POSTMEMORY" ];   then POSTMEMORY=$MEMORY; fi
    if [ -z "$POSTWALLTIME" ]; then POSTWALLTIME=$WALLTIME; fi

    RECIPT=$($BINQSUB -a "$QSUBEXTRA" -W "$MYPBSIDS" -k $CONFIG -m $POSTMEMORY -n $POSTNODES -c $POSTCPU -w $POSTWALLTIME \
            -j $TASK'_postcommand' -o $QOUT/$TASK/postcommand.out --command "$POSTCOMMAND2")

    echo -e "Jobnumber $RECIPT"

    fi
fi

if [ ! -n "$KEEP" ] && [ -e $QOUT/$TASK/runnow.tmp ]; then  rm $QOUT/$TASK/runnow.tmp ; fi
