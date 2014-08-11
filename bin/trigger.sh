#!/bin/bash -e

# Overall trigger pipeline
# author: Denis C. Bauer
# date: Nov.2010

function usage {
echo -e "usage: $(basename $0) CONFIG [TASK]

Main script that interprets the CONFIG file in the project directory and submits
tasks to the HPC queue.

required:
  CONFIG     config.txt file specifying the tasks and location of the resources

options for TASK:
  [empty]        start dry-run: creates dirs, delete old files, print what will be done
  new            detect new data that has not been processed yet.
  fetchdata      get data from remote server (via smbclient)
  pushresult     puts results to remote server (via smbclient)
  armed          submit tasks to the queue
  debug          run task directly (e.g. on node after qrsh) output is written to stdout
  direct         run task directly (e.g. on node after qrsh) output is written to logfiles
  postonly       run only the post analysis steps of a task (if available)
  recover        pick up unfinished business (interrupted jobs)
  recoverdirect  pick up unfinished business (interrupted jobs) on local machine
  html           checks logfiles for errors and creates summary HTML page
  report         alias for html
  trackhubs      generate trackhubs

other options:
  -h         print this help message.
  -v         print version number of ngsane
"
exit
}

function version {

    NGSANE_VERSION=$0
    if [ -e ${NGSANE_VERSION/bin\/trigger.sh/}/.git ]; then 
	    NGSANE_VERSION=`cd ${NGSANE_VERSION/bin\/trigger.sh/}/ && git rev-parse HEAD `" (git hash)"
    else
        NGSANE_VERSION="v0.4.0.3"
    fi
    echo -e "NGSANE version: $NGSANE_VERSION"
    exit
}


if [ ! $# -gt 0 ]; then usage ; fi

while getopts "hv" opt;
do
	case ${opt} in
	h) usage;;
	v) version;;
	\?) print >&2 "$0: error - unrecognized option $1" 
		exit 1;;
	esac
done

CONFIG=$1
ADDITIONALTASK=$2

# convert possibly relative path of CONFIG to absolute path
ABSPATH=`cd \`dirname "$CONFIG"\`; pwd`"/"`basename "$CONFIG"`
CONFIG=$ABSPATH

# check if CONFIG file exists
[ ! -f $CONFIG ] && echo -e "\e[91m[ERROR]\e[0m config file ($CONFIG) not found." && exit 1

# get all the specs defined in the config and defaults from the header (note: sourcing config twice is necessary)
. $CONFIG
# check environment variable
if [ -z ${NGSANE_BASE} ]; then 
    echo -e "[NOTE] NGSANE_BASE environment variable not set. Infering from trigger.sh location";
    TRIGGERDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
    export NGSANE_BASE=${TRIGGERDIR%/bin}
    echo -e "[NOTE] NGSANE_BASE set to $NGSANE_BASE"
fi
. ${NGSANE_BASE}/conf/header.sh
. $CONFIG

################################################################################
#  task fork
################################################################################
if [ -n "$ADDITIONALTASK" ]; then

	if [ "$ADDITIONALTASK" = "fetchdata" ]; then
	    echo -e "\e[35m[NGSANE]\e[0m Trigger mode: \e[4m$ADDITIONALTASK\e[24m"
	    ${NGSANE_BASE}/core/fetchRawDataFromServer.sh -k $CONFIG
	    exit
	    
	elif [ "$ADDITIONALTASK" = "pushresult" ]; then
	    echo -e "\e[35m[NGSANE]\e[0m Trigger mode: \e[4m$ADDITIONALTASK\e[24m"
	    ${NGSANE_BASE}/core/pushResultToServer.sh -k $CONFIG
	    exit
	    
    elif [[ "$ADDITIONALTASK" = "html" || "$ADDITIONALTASK" = "report" ]]; then
        echo -e "\e[35m[NGSANE]\e[0m Trigger mode: \e[4m$ADDITIONALTASK\e[24m"
        ${NGSANE_BASE}/core/makeSummary.sh -k $CONFIG
        exit

    elif [[ "$ADDITIONALTASK" = "trackhubs" ]]; then
        echo -e "\e[35m[NGSANE]\e[0m Trigger mode: \e[4m$ADDITIONALTASK\e[24m"
        ${NGSANE_BASE}/core/makeTrackhubs.sh -k $CONFIG
        exit
        
    elif [ "$ADDITIONALTASK" = "armed" ]; then
	    echo -e "\e[35m[NGSANE]\e[0m Trigger mode: \e[4m$ADDITIONALTASK\e[24m"
	    ARMED="--armed"

        echo -n -e "Double check! Then type \e[4msafetyoff\e[24m and hit enter to launch the job: "
        read safetyoff
        if [ "$safetyoff" != "safetyoff" ];then
            echo -e "Holstering..."
            exit 0
        else
            echo -e "... take cover!"
        fi

    elif [ "$ADDITIONALTASK" = "forcearmed" ]; then
        echo -e "\e[35m[NGSANE]\e[0m Trigger mode: \e[4m$ADDITIONALTASK\e[24m"
        ARMED="--armed"

    elif [ "$ADDITIONALTASK" = "keep" ]; then
        echo -e "\e[35m[NGSANE]\e[0m Trigger mode: \e[4m$ADDITIONALTASK\e[24m"
        ARMED="--keep"

    elif [ "$ADDITIONALTASK" = "new" ]; then
        echo -e "\e[35m[NGSANE]\e[0m Trigger mode: \e[4m$ADDITIONALTASK\e[24m"
        ARMED="--new"

    elif [ "$ADDITIONALTASK" = "debug" ]; then
        echo -e "\e[35m[NGSANE]\e[0m Trigger mode: \e[4m$ADDITIONALTASK\e[24m"
        ARMED="--debug"

    elif [ "$ADDITIONALTASK" = "direct" ]; then
        echo -e "\e[35m[NGSANE]\e[0m Trigger mode: \e[4m$ADDITIONALTASK\e[24m"
        ARMED="--direct"

    elif [ "$ADDITIONALTASK" = "first" ]; then
        echo -e "\e[35m[NGSANE]\e[0m Trigger mode: \e[4m$ADDITIONALTASK\e[24m"
        ARMED="--first --armed"

    elif [ "$ADDITIONALTASK" = "postonly" ]; then
        echo -e "\e[35m[NGSANE]\e[0m Trigger mode: \e[4m$ADDITIONALTASK\e[24m"
        ARMED="--postonly"
        
    elif [ "$ADDITIONALTASK" = "recover" ]; then
        echo -e "\e[35m[NGSANE]\e[0m Trigger mode: \e[4m$ADDITIONALTASK\e[24m"
        ARMED="--recover --armed"

    elif [ "$ADDITIONALTASK" = "recoverdirect" ]; then
        echo -e "\e[35m[NGSANE]\e[0m Trigger mode: \e[4m$ADDITIONALTASK\e[24m"
        ARMED="--recover --direct"

    else
        echo -e "\e[91m[ERROR]\e[0m don't understand $ADDITIONALTASK"
        exit 1
    fi
else
    echo -e "\e[35m[NGSANE]\e[0m Trigger mode: \e[4m[empty]\e[24m (dry run)"
	ARMED="--dryrun"
fi

################################################################################
# test if source data is defined
echo -e "\e[93m[NOTE]\e[0m Folders: ${DIR[@]}"
if [[ -z "${DIR[@]}" ]]; then
    echo -e "\e[91m[ERROR]\e[0m no input directories specified (DIR)."
    exit 1
fi

################################################################################
# create output directories
for dir in ${DIR[@]}; do
    DIRNAME=${dir%%/*} # get (first) folder name
    if [ ! -d $OUT/$DIRNAME ]; then mkdir -p $OUT/$DIRNAME; fi
done

if [ ! -d $QOUT ]; then mkdir -p $QOUT; fi
if [ ! -d $TMP ]; then mkdir -p $TMP; fi

################################################################################
################################################################################
################################################################################
#
#  Pipeline task definitions
#
################################################################################
################################################################################
################################################################################

# source module triggers
for RUN_MODS in ${NGSANE_BASE}/mods/run.d/* ; do source $RUN_MODS; done


