#!/bin/bash -e

# transfers the raw data from a different server using smbclient
# author: Fabian Buske
# date: May.2013


echo ">>>>> Transfer data from HPC cluster "
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> pushResultToServer.sh $*"


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

if ! hash smbclient; then
  echo "[ERROR] Could not find smbclient"
  exit 1
fi

if [ ! -f ~/.smbclient ]; then
  echo "[WARN] ~/.smbclient not configured"
fi

# test if source server is given
if [ -z "$SOURCE_SERVER" ]; then
  echo "[ERROR] source server not specified (TARGET_SERVER)."
fi

# test if source files are given
if [ -z "$TARGET_LOCATION" ]; then
  echo "[ERROR] target location not specified (TARGET_LOCATION)."
  exit 1
fi

# test if source data is defined
echo "${DIR[@]}"
if [ -z "${DIR[@]}" ]; then
  echo "[ERROR] no input directories specified (DIR)."
  exit 1
fi

# ensure out directory is there 
for dir in ${DIR[@]}; do
  if [ ! -d $SOURCE/fastq/$dir ]; then 
  	echo "[ERROR] source data not found (DIR)."
  	exit 1
  fi
done

if [ ! -d $QOUT ]; then mkdir -p $QOUT; fi

echo "********* write results to REMOTE server"
for dir in ${DIR[@]}; do
    for f in $(ls $SOURCE); do
   		fn="${f##*/}" # basename
   		# don't copy the raw files or the tmp folder
   		if [ "$fn" != "fastq" ] && [ "$fn" != "tmp" ]; then 
   			echo "smbclient ${SOURCE_SERVER} -A ~/.smbclient -c \"prompt; recurse; cd ${dn}; mput ${dir}\""
#	   		if [ -f  ~/.smbclient ]; then
#			   smbclient ${SOURCE_SERVER} -A ~/.smbclient -c "prompt; recurse; cd ${dn}; mput ${dir}"
#			else
#			   smbclient ${SOURCE_SERVER} -U `whoami` -c "prompt; recurse; cd ${dn}; mput ${dir}"
#			fi
		fi
	done
done
echo ">>>>> enddate "`date`
