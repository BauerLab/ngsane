#!/bin/bash -e

# transfers the raw data from a different server using smbclient
# author: Fabian Buske
# date: May.2013


echo ">>>>> Transfer data from HPC cluster "
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> $(basename $0) $*"

function usage {
echo -e "usage: $(basename $0) -k CONFIG

Script pushes results from the cluster to a remote server via smbclient
It expects a config file with SOURCE_SERVER SOURCE_FILES and DIR specified.

required:
  -k | --toolkit <path>     location of the NGSANE repository 
"
exit
}

if [ ! $# -gt 1 ]; then usage ; fi

#INPUTS
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the NGSANE repository
        -h | --help )           usage ;;
        * )                     echo "don't understand "$1
    esac
    shift
done

#PROGRAMS
. $CONFIG
. ${NGSANE_BASE}/conf/header.sh
. $CONFIG

################################################################################
CHECKPOINT="programs"

if ! hash smbclient; then
  echo "[ERROR] Could not find smbclient"
  exit 1
fi

echo -e "\n********* $CHECKPOINT"
################################################################################
CHECKPOINT="parameters"

if [ ! -f ~/.smbclient ]; then
  echo "[WARN] ~/.smbclient not configured"
fi

# test if source server is given
if [ -z "$TARGET_SERVER" ]; then
  echo "[ERROR] source server not specified (TARGET_SERVER)."
  exit 1 
fi

# test if source files are given
if [ -z "$TARGET_LOCATION" ]; then
  echo "[ERROR] target location not specified (TARGET_LOCATION)."
  exit 1
fi

echo -e "\n********* $CHECKPOINT"
################################################################################
CHECKPOINT="write results to REMOTE server"

for f in $(ls $SOURCE); do
	fn="${f##*/}" # basename
   	# don't copy the raw files or the tmp folder
   	if [ "$fn" != "fastq" ] && [ "$fn" != "tmp" ]; then 
 		if [ -f  ~/.smbclient ]; then
		   RUN_COMMAND="smbclient ${TARGET_SERVER} -A ~/.smbclient -c \"prompt; recurse; cd ${TARGET_LOCATION}; mput ${fn}\""
		else
		   RUN_COMMAND="smbclient ${TARGET_SERVER} -U `whoami` -c \"prompt; recurse; cd ${TARGET_LOCATION}; mput ${fn}\""
		fi
		echo $RUN_COMMAND
		eval $RUN_COMMAND
	fi
done

################################################################################
echo ">>>>> Transfer data from HPC cluster - FINISHED"
echo ">>>>> enddate "`date`

