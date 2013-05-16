#!/bin/bash -e

# transfers the raw data from a different server using smbclient
# author: Fabian Buske
# date: May.2013


echo ">>>>> Transfer data to HPC cluster "
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> getRawDataFromServer.sh $*"


function usage {
echo -e "usage: $(basename $0) -k NGSANE

Script running hicup including reference genome digestion, read mapping for single and 
paired DNA reads with bowtie from fastq files
It expects a fastq file, pairdend, reference genome and digest pattern  as input.

required:
  -k | --toolkit <path>     location of the NGSANE repository 
  -o | --outdir <path>      output dir
"
exit
}

if [ ! $# -gt 3 ]; then usage ; fi

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
  echo "[ERROR] source server not specified (SOURCE_SERVER)."
fi

# test if source files are given
echo "${SOURCE_FILES[@]}"
if [ -z "${SOURCE_FILES[@]}" ]; then
  echo "[ERROR] no raw data  specified (SOURCE_FILES)."
  exit 1
fi

# test if source data is defined
echo "${DIR[@]}"
if [ -z "${DIR[@]}" ]; then
  echo "[ERROR] no input directories specified (DIR)."
  exit 1
fi

# test if multiple source data is defined
if [ "${DIR[@]}" > 1 ]; then
  echo "[ERROR] multiple input directories specified (DIR)."
  exit 1
fi

# ensure out directory is there 
for dir in ${DIR[@]}; do
  if [ ! -d $OUT/$dir ]; then mkdir -p $OUT/$dir; fi
done

if [ ! -d $QOUT ]; then mkdir -p $QOUT; fi


CURDIR=$(pwd)
cd ${DIR[0]}

echo "********* get file from source folder"
for sourcefile in ${SOURCE_FILES[@]}; do
	fn="${sourcefile##*/}" # filename
	dn="${sourcefile%/*}"  # dirname
	
	if [ -f  ~/.smbclient ]; then
	   smbclient ${SOURCE_SERVER} -A ~/.smbclient -c "cd ${dn}; get ${fn}"
	else
	   smbclient ${SOURCE_SERVER} -U `whoami` -c "cd ${dn}; get ${fn}"
	fi
	
	# get second read
	if [ -n "$READTWO" ]; then  
		if [ -f  ~/.smbclient ]; then
		   smbclient ${SOURCE_SERVER} -A ~/.smbclient -c "cd ${dn}; get ${fn/$READONE/$READTWO}"
		else
		   smbclient ${SOURCE_SERVER} -U `whoami` -c "cd ${dn}; get ${fn/$READONE/$READTWO}"
		fi
	fi	
done
echo ">>>>> enddate "`date`
