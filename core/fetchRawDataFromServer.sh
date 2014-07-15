#!/bin/bash -e

# transfers the raw data from a different server using smbclient
# author: Fabian Buske
# date: May.2013

echo ">>>>> Transfer data to HPC cluster "
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> $(basename $0) $*"

function usage {
echo -e "usage: $(basename $0) -k CONFIG

Script fetches data from a remote server via smbclient
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
if [ -z "${SOURCE_SERVER}" ]; then
    echo "[ERROR] source server not specified (SOURCE_SERVER)."
    exit 1
fi

# test if source files are given
if [ "${#SOURCE_FILES[@]}" -eq "0" ]; then
    echo "[ERROR] no raw data  specified (SOURCE_FILES)."
    exit 1
else
    echo "[NOTE] Number of source file patterns: ${#SOURCE_FILES[@]}"
fi

# test if fastq folder is defined
if [ -z "${DIR[@]}" ]; then
    echo "[ERROR] no input directories specified (DIR)."
    exit 1
else
    echo "[NOTE] place in fastq/${DIR[@]}"
fi

# test if multiple source data is defined
if [ "${#DIR[@]}" -ne "1" ]; then
    echo "[ERROR] multiple input directories specified (DIR)."
    exit 1
fi

if [ -z "${SOURCE_FOLDER}" ]; then
    echo "[ERROR] no path on source server specified (SOURCE_FOLDER)"
    exit 1
else
    echo "[NOTE] will fetch from ${SOURCE_SERVER}/${SOURCE_FOLDER}"
fi

# ensure out directory is there 
for dir in ${DIR[@]}; do
    if [ ! -d $SOURCE/fastq/$dir ]; then mkdir -p $SOURCE/fastq/$dir; fi
done


if [ ! -d $QOUT ]; then mkdir -p $QOUT; fi

if [ -f ~/.smbclient ]; then
    AUTH="-A ~/.smbclient"
else
    AUTH="-U `whoami`"
fi	

CURDIR=$(pwd)

THISTMP=$TMP"/"$(whoami)"/"$(echo ${SOURCE}/fastq/${DIR[0]} | md5sum | cut -d' ' -f1)
[ -d $THISTMP ] && rm -r $THISTMP
mkdir -p $THISTMP

echo -e "\n********* $CHECKPOINT"
################################################################################
CHECKPOINT="fetch files"

cd $THISTMP

for sourcefile in ${SOURCE_FILES[@]}; do

    RUN_COMMAND="smbclient ${SOURCE_SERVER} $AUTH -c 'prompt; recurse; cd ${SOURCE_FOLDER}; mask *$READONE.$FASTQ; mget ${sourcefile}'"
    echo $RUN_COMMAND && eval $RUN_COMMAND
	
	# check second read
	if [ -n "$READTWO" ] && [ "$READONE" != "$READTWO" ]; then  
		echo "[NOTE] get second read"
		RUN_COMMAND="smbclient ${SOURCE_SERVER} $AUTH -c 'prompt; recurse; cd ${SOURCE_FOLDER}; mask *$READTWO.$FASTQ; mget ${sourcefile/%$READONE.$FASTQ/$READTWO.$FASTQ}'"
        	echo $RUN_COMMAND && eval $RUN_COMMAND
	fi
done

# move to get a flat hierarchy
for d in $(find . -mindepth 1  -type f -name "*$READONE.$FASTQ" ); do 
    mv $d ${SOURCE}/fastq/${DIR[0]}
done
if [ -n "$READTWO" ] && [ "$READONE" != "$READTWO" ]; then  
    for d in $(find . -mindepth 1 -type f -name "*$READTWO.$FASTQ" ); do 
        mv $d ${SOURCE}/fastq/${DIR[0]}
    done
fi
# remove dirs
for d in $(find . -mindepth 1 -maxdepth 1 -type d ); do 
    rm -r $d
done

rm -r $THISTMP

echo -e "\n********* $CHECKPOINT"
################################################################################
echo ">>>>> Transfer data to HPC cluster - FINISHED"
echo ">>>>> enddate "`date`
