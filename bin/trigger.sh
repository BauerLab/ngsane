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
  [empty]    start dry-run: creates dirs, delete old files, print what will be done
  new        detect new data that has not been processed yet.
  fetchdata  get data from remote server (via smbclient)
  pushresult puts results to remote server (via smbclient)
  armed      submit tasks to the queue
  debug      run task directly (e.g. on node after qrsh) output is written to stdout
  direct     run task directly (e.g. on node after qrsh) output is written to logfiles
  postonly   run only the post analysis steps of a task (if available)
  recover    pick up unfinished business (interrupted jobs)
  html       checks logfiles for errors and creates summary HTML page
  report     alias for html
  trackhubs  generate trackhubs

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
# define functions 
#
# getJobIds takes 1 parameter
# $1=QSUB output
function waitForJobIds {
    JOBIDS=$(echo -e "$1" | grep "Jobnumber")
    if [ -n "$JOBIDS" ]; then 
		#JOBIDS=$(echo -e $JOBIDS | cut -d " " -f 2 | tr '\n' ':' | sed 's/:$//g' )
        JOBIDS=$(echo -e $JOBIDS | gawk '{ ORS=" "; print; }' | sed 's/Jobnumber //g' | sed 's/ /:/g' )

    fi
    if [ "$JOBIDS" != "" ]; then
        echo "-W $JOBIDS"
    else
        echo ""
    fi
} 

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

    else
        echo -e "\e[91m[ERROR]\e[0m don't understand $ADDITIONALTASK"
        exit 1
    fi
else
    echo -e "\e[35m[NGSANE]\e[0m Trigger mode: \e[4[empty]\e[24m (dry run)"
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


################################################################################
#   FastQC summary of fastq files
#
# IN : $SOURCE/fastq/$dir/*read1.fastq
# OUT: $OUT/$dir/fastQC/*
################################################################################
if [ -n "$RUNFASTQC" ]; then
    if [ -z "$TASK_FASTQC" ] || [ -z "$NODES_FASTQC" ] || [ -z "$CPU_FASTQC" ] || [ -z "$MEMORY_FASTQC" ] || [ -z "$WALLTIME_FASTQC" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi
    
    $QSUB $ARMED -d -k $CONFIG -t $TASK_FASTQC -i $INPUT_FASTQC -e $READONE.$FASTQ -n $NODES_FASTQC \
    	-c $CPU_FASTQC -m $MEMORY_FASTQC"G" -w $WALLTIME_FASTQC \
    	--command "${NGSANE_BASE}/mods/fastQC.sh -k $CONFIG -f <FILE> -o $OUT/<DIR>/$TASK_FASTQC" 
fi

################################################################################
#   FASTQSCREEN 
#
# IN : $SOURCE/fastq/$dir/*read1.fastq
# OUT: $OUT/$dir/fastqscreen/*.
################################################################################

if [ -n "$RUNFASTQSCREEN" ]; then
    $QSUB $ARMED -d -k $CONFIG -t $TASK_FASTQSCREEN -i $INPUT_FASTQSCREEN -e $READONE.$FASTQ -n $NODES_FASTQSCREEN \
        -c $CPU_FASTQSCREEN -m $MEMORY_FASTQSCREEN"G" -w $WALLTIME_FASTQSCREEN \
        --command "$NGSANE_BASE/mods/fastqscreen.sh -k $CONFIG -f <FILE>  -o $OUT/<DIR>/$TASK_FASTQSCREEN"
fi

################################################################################
#   Blue
#
# IN : $SOURCE/fastq/$dir/*read1.fastq
# OUT: $OUT/$dir_healed/*.
################################################################################

if [ -n "$RUNBLUE" ]; then
    if [ -z "$TASK_BLUE" ] || [ -z "$NODES_BLUE" ] || [ -z "$CPU_BLUE" ] || [ -z "$MEMORY_BLUE" ] || [ -z "$WALLTIME_BLUE" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi
    
    $QSUB $ARMED -d -k $CONFIG -t $TASK_BLUE -i $INPUT_BLUE -e $READONE.$FASTQ -n $NODES_BLUE \
	   -c $CPU_BLUE -m $MEMORY_BLUE"G" -w $WALLTIME_BLUE \
	   --command "${NGSANE_BASE}/mods/blue.sh -k $CONFIG -f <FILE> -o $INPUT_BLUE/<DIR>"_"$TASK_BLUE/" 

fi




################################################################################
#   CUTADAPT remove contaminants
#
# IN : $SOURCE/fastq/$dir/*read1.fastq
# OUT: $SOURCE/fastq/$dir_cutadapt/*read1.fastq
################################################################################

if [ -n "$RUNCUTADAPT" ]; then
    if [ -z "$TASK_CUTADAPT" ] || [ -z "$NODES_CUTADAPT" ] || [ -z "$CPU_CUTADAPT" ] || [ -z "$MEMORY_CUTADAPT" ] || [ -z "$WALLTIME_CUTADAPT" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi
    
    $QSUB $ARMED -d -k $CONFIG -t $TASK_CUTADAPT -i $INPUT_CUTADAPT -e $READONE.$FASTQ -n $NODES_CUTADAPT \
        -c $CPU_CUTADAPT -m $MEMORY_CUTADAPT"G" -w $WALLTIME_CUTADAPT \
        --command "${NGSANE_BASE}/mods/cutadapt.sh -k $CONFIG -f <FILE>" 
fi

################################################################################
#   TRIMMOMATIC remove contaminants
#
# IN : $SOURCE/fastq/$dir/*read1.fastq
# OUT: $SOURCE/fastq/$dir_trimmomatic/*read1.fastq
################################################################################

if [ -n "$RUNTRIMMOMATIC" ]; then
    if [ -z "$TASK_TRIMMOMATIC" ] || [ -z "$NODES_TRIMMOMATIC" ] || [ -z "$CPU_TRIMMOMATIC" ] || [ -z "$MEMORY_TRIMMOMATIC" ] || [ -z "$WALLTIME_TRIMMOMATIC" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi
    
    $QSUB $ARMED -d -k $CONFIG -t $TASK_TRIMMOMATIC -i $INPUT_TRIMMOMATIC -e $READONE.$FASTQ -n $NODES_TRIMMOMATIC \
        -c $CPU_TRIMMOMATIC -m $MEMORY_TRIMMOMATIC"G" -w $WALLTIME_TRIMMOMATIC \
        --command "$NGSANE_BASE/mods/trimmomatic.sh -k $CONFIG -f <FILE>"
fi

################################################################################
#   TRIMGALORE remove contaminants
#
# IN : $SOURCE/fastq/$dir/*read1.fastq
# OUT: $SOURCE/fastq/$dir_trimgalore/*read1.fastq
################################################################################ 
if [ -n "$RUNTRIMGALORE" ]; then
    if [ -z "$TASK_TRIMGALORE" ] || [ -z "$NODES_TRIMGALORE" ] || [ -z "$CPU_TRIMGALORE" ] || [ -z "$MEMORY_TRIMGALORE" ] || [ -z "$WALLTIME_TRIMGALORE" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi
    
    $QSUB $ARMED -d -k $CONFIG -t $TASK_TRIMGALORE -i $INPUT_TRIMGALORE -e $READONE.$FASTQ -n $NODES_TRIMGALORE \
        -c $CPU_TRIMGALORE -m $MEMORY_TRIMGALORE"G" -w $WALLTIME_TRIMGALORE \
        --command "${NGSANE_BASE}/mods/trimgalore.sh -k $CONFIG -f <FILE>"
fi

################################################################################ 
# IN : */bwa/*.bam
# OUT: */bwa/*.ann
################################################################################ 
if [ -n "$RUNANNOTATINGBAM" ]; then
    if [ -z "$TASK_BAMANN" ] || [ -z "$NODES_BAMANN" ] || [ -z "$CPU_BAMANN" ] || [ -z "$MEMORY_BAMANN" ] || [ -z "$WALLTIME_BAMANN" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi
    
    $QSUB --nodir -r $ARMED -k $CONFIG -t $TASK_BAMANN -i $INPUT_BAMANN -e .$ASD.bam \
    -n $NODES_BAMANN -c $CPU_BAMANN -m $MEMORY_BAMANN'G' -w $WALLTIME_BAMANN \
        --command "${NGSANE_BASE}/mods/annotateBam.sh -k $CONFIG -f <FILE>"
fi


################################################################################
#   Variance calling
# IN: */bwa/*.bam
# OUT: */bwa_var/*.clean.vcf
################################################################################

if [ -n "$RUNSAMVAR" ]; then
    if [ -z "$TASK_BWA" ] || [ -z "$NODES_SAMVAR" ] || [ -z "$CPU_SAMVAR" ] || [ -z "$MEMORY_SAMVAR" ] || [ -z "$WALLTIME_SAMVAR" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi
    
    $QSUB -r $ARMED -k $CONFIG -t $INPUT_SAMVAR-$TASK_SAMVAR -i $INPUT_SAMVAR -e .$ASD.bam \
       -n $NODES_SAMVAR -c $CPU_SAMVAR -m $MEMORY_SAMVAR'G' -w $WALLTIME_SAMVAR \
		--postnodes $NODES_VARCOLLECT --postcpu $CPU_VARCOLLECT \
		--postwalltime $WALLTIME_VARCOLLECT --postmemory $MEMORY_VARCOLLECT"G" \
        --command "${NGSANE_BASE}/mods/samSNPs.sh -k $CONFIG -f <FILE> -o $OUT/<DIR>/$INPUT_SAMVAR-$TASK_SAMVAR" \
	    --postcommand "${NGSANE_BASE}/mods/variantcollect.sh -k $CONFIG -f <FILE> -i1 $INPUT_SAMVAR \
				-i2 $INPUT_SAMVAR-$TASK_SAMVAR -o $OUT/variant/${INPUT_SAMVAR}-$INPUT_SAMVAR-$TASK_SAMVAR-<DIR> "

	#   --postcommand "${NGSANE_BASE}/mods/samSNPscollect.sh -k $CONFIG -f <FILE> -o $OUT/variant/$TASK_BWA-$TASK_SAMVAR"-"<DIR>"
fi

################################################################################
#   call indels with GATK -- call one vcf file over all folders
################################################################################

if [ -n "$RUNVARCALLS" ]; then
    NAME=$(echo ${DIR[@]}|sed 's/ /_/g')
    $QSUB $ARMED -r -d -k $CONFIG -t $TASK_GATKVAR -i $INPUT_GATKVAR  -e .$ASR.bam -n $NODES_GATKVAR \
        -c $CPU_GATKVAR -m $MEMORY_GATKVAR"G" -w $WALLTIME_GATKVAR \
        --postcommand "${NGSANE_BASE}/mods/gatkVARs.sh -k $CONFIG \
                        -i <FILE> -t $CPU_GATKVAR \
                        -r $FASTA -d $DBSNPVCF -o $OUT/$TASK_GATKVAR/$NAME -n $NAME \
                        -H $HAPMAPVCF" #-K $ONEKGVCF"
fi


################################################################################
#   call indels with GATK -- call one vcf file over all folders but in batches (by chr)
################################################################################

if [ -n "$RUNVARCALLSBATCH" ]; then
	if [ ! -e "$FASTA.fai" ] ; then echo -e "\e[91m[ERROR]\e[0m $FASTA.fai missing"; exit 1; fi
  	BATCHES=$(cut -f 1 $FASTA.fai | grep -v GL | sort -u)
	NAME=$(echo ${DIR[@]} | sort -u |sed 's/ /_/g')

	if [[ ! "$ARMED" -eq "postonly" ]]; then
	  	for i in $BATCHES; do
			echo "[NOTE] Batch $i"
			export ADDDUMMY=$i
	    	JOBID=$( $QSUB $ARMED --postname postcommand$i -r -d -k $CONFIG -t ${TASK_GATKVAR}batch -i $INPUT_GATKVAR  \
				-e .$ASR.bam -n $NODES_GATKVAR -c $CPU_GATKVAR -m $MEMORY_GATKVAR"G" -w $WALLTIME_GATKVAR \
	        	--postcommand "${NGSANE_BASE}/mods/gatkVARs.sh -k $CONFIG \
	                        -i <FILE> -t $CPU_GATKVAR \
	                        -r $FASTA -d $DBSNPVCF -o $OUT/${TASK_GATKVAR}batch/$NAME -n $NAME$ADDDUMMY \
	                        -H $HAPMAPVCF -L $i " 
			) && echo -e "$JOBID"
			if [ -n "$(echo $JOBID | grep Jobnumber)" ]; then JOBIDS=$(waitForJobIds "$JOBID")":"$JOBIDS; fi
	  	done
		JOBIDS=${JOBIDS//-W /}
		JOBIDS=${JOBIDS//::/:}
        [ -n "$JOBIDS" ] && JOBIDS="-W $JOBIDS"
	fi

	echo "[NOTE] filtered SNPs"
   	$QSUB $ARMED --postname joinedSNP --givenDirs $NAME -d -k $CONFIG -t ${TASK_GATKVAR}batch -i ${TASK_GATKVAR}batch -e filter.snps.vcf $JOBIDS \
		-n $NODES_VARCOLLECT -c $CPU_VARCOLLECT -m $MEMORY_VARCOLLECT"G" -w $WALLTIME_VARCOLLECT \
		--postcommand "${NGSANE_BASE}/mods/variantcollect.sh -k $CONFIG -f <FILE> -i1 ${TASK_GATKVAR}batch \
				-i2 ${TASK_GATKVAR}batch -o $OUT/${TASK_GATKVAR}batch/$NAME --dummy filter.snps.vcf --target filter.snps.vcf"

	echo "[NOTE] filtered INDELs"
   	$QSUB $ARMED --postname joinedINDEL --givenDirs $NAME -d -k $CONFIG -t ${TASK_GATKVAR}batch -i ${TASK_GATKVAR}batch -e filter.snps.vcf $JOBIDS \
		-n $NODES_VARCOLLECT -c $CPU_VARCOLLECT -m $MEMORY_VARCOLLECT"G" -w $WALLTIME_VARCOLLECT \
		--postcommand "${NGSANE_BASE}/mods/variantcollect.sh -k $CONFIG -f <FILE> -i1 ${TASK_GATKVAR}batch \
				-i2 ${TASK_GATKVAR}batch -o $OUT/${TASK_GATKVAR}batch/$NAME --dummy filter.snps.vcf --target filter.indel.vcf"

	echo "[NOTE] recal Vars"
   	$QSUB $ARMED --postname joinedRECAL --givenDirs $NAME -d -k $CONFIG -t ${TASK_GATKVAR}batch -i ${TASK_GATKVAR}batch -e filter.snps.vcf $JOBIDS \
		-n $NODES_VARCOLLECT -c $CPU_VARCOLLECT -m $MEMORY_VARCOLLECT"G" -w $WALLTIME_VARCOLLECT \
		--postcommand "${NGSANE_BASE}/mods/variantcollect.sh -k $CONFIG -f <FILE> -i1 ${TASK_GATKVAR}batch \
				-i2 ${TASK_GATKVAR}batch -o $OUT/${TASK_GATKVAR}batch/$NAME --dummy filter.snps.vcf --target recalfilt.vcf"
   
fi	

	
################################################################################
#   Depth of Coverage
# IN: */bwa/*.bam
# OUT: */bwa_var/*.clean.vcf
################################################################################

if [ -n "$RUNDEPTHOFCOVERAGE" ]; then
    $QSUB $ARMED -r -k $CONFIG -t $TASK_GATKDOC -i $INPUT_GATKDOC -e .$ASR.bam \
	-n $NODES_GATKDOC -c $CPU_GATKDOC -m $MEMORY_GATKDOC"G" -w $WALLTIME_GATKDOC \
	--command "${NGSANE_BASE}/mods/gatkDOC.sh -k $CONFIG -f <FILE> -r $FASTA -o $OUT/<DIR>/$TASK_GATKDOC -t $CPU_GATKDOC"
fi

################################################################################
#   Mapping using HiCUP
# IN : $SOURCE/fastq/$dir/*$READONE$FASTQ
# OUT: $OUT/$dir/hicup/*_uniques.bam
################################################################################

if [ -n "$RUNHICUP" ]; then
    if [ -z "$TASK_HICUP" ] || [ -z "$NODES_HICUP" ] || [ -z "$CPU_HICUP" ] || [ -z "$MEMORY_HICUP" ] || [ -z "$WALLTIME_HICUP" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi
    
    if [ ! -f ${FASTA%.*}.1.ebwt ];then
        # submit job for index generation if necessary
        INDEXJOBIDS=$(
            $QSUB $ARMED -k $CONFIG -t $TASK_HICUP -i $INPUT_HICUP -e $READONE.$FASTQ -n $NODES_HICUP -c $CPU_HICUP \
    	   -m $MEMORY_HICUP"G" -w $WALLTIME_HICUP --commontask indexGenome \
            --command "${NGSANE_BASE}/mods/bowtieIndex.sh -k $CONFIG"
        ) && echo -e "$INDEXJOBIDS"
        INDEXJOBIDS=$(waitForJobIds "$INDEXJOBIDS")
    else
        INDEXJOBIDS=""
    fi
    
    if [ ! -f $OUT/common/$TASK_HICUP/Digest_${REFERENCE_NAME}_${HICUP_RENZYME1}_${HICUP_RENZYME2}.txt ];then
        JOBIDS=$( 
        $QSUB $ARMED -k $CONFIG -t $TASK_HICUP -i $INPUT_HICUP -e $READONE.$FASTQ -n $NODES_HICUP -c 1 \
        	-m $MEMORY_HICUP"G" -w $WALLTIME_HICUP $INDEXJOBIDS --commontask digestGenome \
            --command "${NGSANE_BASE}/mods/hicupDigestGenome.sh -k $CONFIG" 
        ) && echo -e "$JOBIDS"
        JOBIDS=$(waitForJobIds "$JOBIDS")

    else
        JOBIDS=""
    fi

    $QSUB $ARMED -k $CONFIG -t $TASK_HICUP -i $INPUT_HICUP -e $READONE.$FASTQ -n $NODES_HICUP -c $CPU_HICUP \
    	-m $MEMORY_HICUP"G" -w $WALLTIME_HICUP $JOBIDS \
        --command "${NGSANE_BASE}/mods/hicup.sh -k $CONFIG -f <FILE> -o $OUT/<DIR>/$TASK_HICUP"

fi

################################################################################
#   Assessing chromatin interactions with fit-hi-c
# IN : $SOURCE/$dir/hicup/*$FRAGMENTLIST
# OUT: $OUT/$dir/hicup/*.spline_pass1.q05.txt.gz
################################################################################

if [ -n "$RUNFITHIC" ]; then
    if [ -z "$TASK_FITHIC" ] || [ -z "$NODES_FITHIC" ] || [ -z "$CPU_FITHIC" ] || [ -z "$MEMORY_FITHIC" ] || [ -z "$WALLTIME_FITHIC" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi

    $QSUB $ARMED -r -k $CONFIG -t $TASK_FITHIC -i $INPUT_FITHIC -e $FRAGMENTLIST -n $NODES_FITHIC -c $CPU_FITHIC \
    	-m $MEMORY_FITHIC"G" -w $WALLTIME_FITHIC \
        --command "${NGSANE_BASE}/mods/fithic.sh -k $CONFIG -f <FILE> -o $OUT/<DIR>/$TASK_FITHIC"

fi

################################################################################
#  Assessing HiC data with hiclib
#
# IN: $SOURCE/$dir/fastq/*read1.fastq
# OUT: $OUT/$dir/hiclib/*.hdf5
################################################################################

if [ -n "$RUNHICLIB" ]; then
    if [ -z "$TASK_HICLIB" ] || [ -z "$NODES_HICLIB" ] || [ -z "$CPU_HICLIB" ] || [ -z "$MEMORY_HICLIB" ] || [ -z "$WALLTIME_HICLIB" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi

    if [ ! -f ${FASTA%.*}.1.bt2 ];then
        # submit job for index generation if necessary
        INDEXJOBIDS=$(
            $QSUB $ARMED -k $CONFIG -t $TASK_HICLIB -i $INPUT_HICLIB -e $READONE.$FASTQ -n $NODES_HICLIB -c $CPU_HICLIB -m $MEMORY_HICLIB"G" \
            -w $WALLTIME_HICLIB --commontask indexGenome \
            --command "${NGSANE_BASE}/mods/bowtie2Index.sh -k $CONFIG"
        ) && echo -e "$INDEXJOBIDS"
        INDEXJOBIDS=$(waitForJobIds "$INDEXJOBIDS")
    else
        INDEXJOBIDS=""
    fi
    
    $QSUB $ARMED -k $CONFIG -t $TASK_HICLIB -i $INPUT_HICLIB -e $READONE.$FASTQ \
    	-n $NODES_HICLIB -c $CPU_HICLIB -m $MEMORY_HICLIB"G" -w $WALLTIME_HICLIB \
    	--postnodes $NODES_HICLIB_POSTCOMMAND --postcpu $CPU_HICLIB_POSTCOMMAND $INDEXJOBIDS \
        --command "${NGSANE_BASE}/mods/hiclibMapping.sh -k $CONFIG --fastq <FILE> --outdir $OUT/<DIR>/$TASK_HICLIB" \
        --postcommand "${NGSANE_BASE}/mods/hiclibCorrelate.sh -f <FILE> -k $CONFIG --outdir $OUT/$TASK_HICLIB/$TASK_HICLIB-<DIR>"
fi

################################################################################
#   Mapping using BWA
#
# IN : $SOURCE/$dir/fastq/*read1.fastq
# OUT: $OUT/$dir/bwa/*.$ASD.bam
################################################################################

if [ -n "$RUNMAPPINGBWA" ]; then
    if [ -z "$TASK_BWA" ] || [ -z "$NODES_BWA" ] || [ -z "$CPU_BWA" ] || [ -z "$MEMORY_BWA" ] || [ -z "$WALLTIME_BWA" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi

    if [ ! -f $FASTA.bwt ];then
        # submit job for index generation if necessary
        JOBIDS=$(
            $QSUB $ARMED -k $CONFIG -t $TASK_BWA -i $INPUT_BWA -e $READONE.$FASTQ -n $NODES_BWA -c $CPU_BWA -m $MEMORY_BWA"G" \
            -w $WALLTIME_BWA --commontask indexGenome \
            --command "${NGSANE_BASE}/mods/bwaIndex.sh -k $CONFIG"
        ) && echo -e "$JOBIDS"
        JOBIDS=$(waitForJobIds "$JOBIDS")
    else
        JOBIDS=""
    fi

    $QSUB $ARMED -k $CONFIG -t $TASK_BWA -i $INPUT_BWA -e $READONE.$FASTQ -n $NODES_BWA -c $CPU_BWA -m $MEMORY_BWA"G" -w $WALLTIME_BWA $JOBIDS \
        --command "${NGSANE_BASE}/mods/bwa.sh -k $CONFIG -f <FILE> -o $OUT/<DIR>/$TASK_BWA --rgsi <DIR>"
fi

################################################################################
#   Mapping using Bowtie v1
#
# IN:$SOURCE/$dir/fastq/*read1.fastq
# OUT: $OUT/$dir/bowtie/*$ASD.bam
################################################################################

if [ -n "$RUNMAPPINGBOWTIE" ]; then
    if [ -z "$TASK_BOWTIE" ] || [ -z "$NODES_BOWTIE" ] || [ -z "$CPU_BOWTIE" ] || [ -z "$MEMORY_BOWTIE" ] || [ -z "$WALLTIME_BOWTIE" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi

    if [ ! -f ${FASTA%.*}.1.ebwt ];then
        # submit job for index generation if necessary
        INDEXJOBIDS=$(
            $QSUB $ARMED -k $CONFIG -t $TASK_BOWTIE -i $INPUT_BOWTIE -e $READONE.$FASTQ -n $NODES_BOWTIE -c $CPU_BOWTIE -m $MEMORY_BOWTIE"G" \
            -w $WALLTIME_BOWTIE --commontask indexGenome \
            --command "${NGSANE_BASE}/mods/bowtieIndex.sh -k $CONFIG"
        ) && echo -e "$INDEXJOBIDS"
        INDEXJOBIDS=$(waitForJobIds "$INDEXJOBIDS")
    else
        INDEXJOBIDS=""
    fi
    
    $QSUB $ARMED -k $CONFIG -t $TASK_BOWTIE -i $INPUT_BOWTIE -e $READONE.$FASTQ -n $NODES_BOWTIE -c $CPU_BOWTIE -m $MEMORY_BOWTIE"G" -w $WALLTIME_BOWTIE $INDEXJOBIDS \
        --command "${NGSANE_BASE}/mods/bowtie.sh -k $CONFIG -f <FILE> -o $OUT/<DIR>/$TASK_BOWTIE --rgsi <DIR>"        
fi


################################################################################
#  Mapping with bowtie v2
#
# IN:  $SOURCE/$dir/fastq/*read1.fastq
# OUT: $OUT/$dir/bowtie2/*$ASD.bam
################################################################################
if [ -n "$RUNMAPPINGBOWTIE2" ]; then
    if [ -z "$TASK_BOWTIE2" ] || [ -z "$NODES_BOWTIE2" ] || [ -z "$CPU_BOWTIE2" ] || [ -z "$MEMORY_BOWTIE2" ] || [ -z "$WALLTIME_BOWTIE2" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi
    
    if [ ! -f ${FASTA%.*}.1.bt2 ];then
        # submit job for index generation if necessary
        INDEXJOBIDS=$(
            $QSUB $ARMED -k $CONFIG -t $TASK_BOWTIE2 -i $INPUT_BOWTIE2 -e $READONE.$FASTQ -n $NODES_BOWTIE2 -c $CPU_BOWTIE2 -m $MEMORY_BOWTIE2"G" \
            -w $WALLTIME_BOWTIE2 --commontask indexGenome \
            --command "${NGSANE_BASE}/mods/bowtie2Index.sh -k $CONFIG"
        ) && echo -e "$INDEXJOBIDS"
        INDEXJOBIDS=$(waitForJobIds "$INDEXJOBIDS")
    else
        INDEXJOBIDS=""
    fi
    
    $QSUB $ARMED -k $CONFIG -t $TASK_BOWTIE2 -i $INPUT_BOWTIE2 -e $READONE.$FASTQ -n $NODES_BOWTIE2 -c $CPU_BOWTIE2 -m $MEMORY_BOWTIE2"G" -w $WALLTIME_BOWTIE2 $INDEXJOBIDS \
	       --command "${NGSANE_BASE}/mods/bowtie2.sh -k $CONFIG -f <FILE> -o $OUT/<DIR>/$TASK_BOWTIE2 --rgsi <DIR>"
      
fi

################################################################################
#   Mapping using Masai
#
# IN:$SOURCE/$dir/fastq/*read1.fastq
# OUT: $OUT/$dir/bowtie/*.bam
################################################################################

# TODO: finish pipeline
#if [ -n "$RUNMASAI" ]; then
#    if [ -z "$TASK_MASAI" ] || [ -z "$NODES_MASAI" ] || [ -z "$CPU_MASAI" ] || [ -z "$MEMORY_MASAI" ] || [ -z "$WALLTIME_MASAI" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi
#
#    $QSUB $ARMED -k $CONFIG -t $TASK_MASAI -i $INPUT_MASAI -e $READONE.$FASTQ -n $NODES_MASAI -c $CPU_MASAI -m $MEMORY_MASAI"G" -w $WALLTIME_MASAI \
#        --command "${NGSANE_BASE}/mods/masai.sh -k $CONFIG -f <FILE> -o $OUT/<DIR>/$TASK_MASAI --rgsi <DIR>"        
#fi


################################################################################
#   Mapping using rrbsmap
#
# IN:$SOURCE/$dir/fastq/*read1.fastq
# OUT: $OUT/$dir/rrbsmap/*.bam
################################################################################

if [ -n "$RUNMAPPINGRRBS" ]; then
    if [ -z "$TASK_RRBSMAP" ] || [ -z "$NODES_RRBSMAP" ] || [ -z "$CPU_RRBSMAP" ] || [ -z "$MEMORY_RRBSMAP" ] || [ -z "$WALLTIME_RRBSMAP" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi

    $QSUB $ARMED -k $CONFIG -t $TASK_RRBSMAP -i $INPUT_RRBSMAP -e $READONE.$FASTQ -n $NODES_RRBSMAP -c $CPU_RRBSMAP -m $MEMORY_RRBSMAP"G" -w $WALLTIME_RRBSMAP \
        --command "${NGSANE_BASE}/mods/rrbsmap.sh -k $CONFIG -f <FILE> -r $FASTA \
            -o $OUT/<DIR>/$TASK_RRBSMAP --rgid $EXPID --rglb $LIBRARY --rgpl $PLATFORM --rgsi <DIR>"
fi


################################################################################
#  Transcript mapping tophat
#
# IN : $SOURCE/$dir/fastq/*read1.fastq
# OUT: $OUT/$dir/tophat/*.bam/
################################################################################       

if [ -n "$RUNTOPHAT" ]; then
    if [ -z "$TASK_TOPHAT" ] || [ -z "$NODES_TOPHAT" ] || [ -z "$CPU_TOPHAT" ] || [ -z "$MEMORY_TOPHAT" ] || [ -z "$WALLTIME_TOPHAT" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi

    if [ ! -f ${FASTA%.*}.1.bt2 ];then
        # submit job for index generation if necessary
        INDEXJOBIDS=$(
            $QSUB $ARMED -k $CONFIG -t $TASK_TOPHAT -i $INPUT_TOPHAT -e $READONE.$FASTQ -n $NODES_TOPHAT -c $CPU_TOPHAT -m $MEMORY_TOPHAT"G" \
            -w $WALLTIME_TOPHAT --commontask indexGenome \
            --command "${NGSANE_BASE}/mods/bowtie2Index.sh -k $CONFIG"
        ) && echo -e "$INDEXJOBIDS"
        INDEXJOBIDS=$(waitForJobIds "$INDEXJOBIDS")
    else
        INDEXJOBIDS=""
    fi

    $QSUB $ARMED -k $CONFIG -t $TASK_TOPHAT -i $INPUT_TOPHAT -e $READONE.$FASTQ -n $NODES_TOPHAT -c $CPU_TOPHAT -m $MEMORY_TOPHAT"G" -w $WALLTIME_TOPHAT $INDEXJOBIDS \
        --command "${NGSANE_BASE}/mods/tophat.sh -k $CONFIG -f <FILE> -o $OUT/<DIR>/$TASK_TOPHAT/<NAME>"

fi

################################################################################
#  Gene expression analysis with cufflinks
#
# IN : $OUT/$dir/tophat/*.bam
# OUT: $OUT/$dir/cufflinks/*_transcript.gtf
################################################################################       

if [ -n "$RUNCUFFLINKS" ]; then
    if [ -z "$TASK_CUFFLINKS" ] || [ -z "$NODES_CUFFLINKS" ] || [ -z "$CPU_CUFFLINKS" ] || [ -z "$MEMORY_CUFFLINKS" ] || [ -z "$WALLTIME_CUFFLINKS" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi

    $QSUB $ARMED -r -k $CONFIG -t $TASK_CUFFLINKS -i $INPUT_CUFFLINKS -e .$ASD.bam -n $NODES_CUFFLINKS -c $CPU_CUFFLINKS -m $MEMORY_CUFFLINKS"G" -w $WALLTIME_CUFFLINKS  \
        --command "${NGSANE_BASE}/mods/cufflinks.sh -k $CONFIG -f <FILE> -o $OUT/<DIR>/$TASK_CUFFLINKS/<NAME>"

fi

################################################################################
#  Feature counting with HTSEQCOUNT
#
# IN : $OUT/$dir/tophat/*.bam
# OUT: $OUT/$dir/htseqcount/*_transcript.gtf
################################################################################       

if [ -n "$RUNHTSEQCOUNT" ]; then
    if [ -z "$TASK_HTSEQCOUNT" ] || [ -z "$NODES_HTSEQCOUNT" ] || [ -z "$CPU_HTSEQCOUNT" ] || [ -z "$MEMORY_HTSEQCOUNT" ] || [ -z "$WALLTIME_HTSEQCOUNT" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi

    $QSUB $ARMED -r -k $CONFIG -t $TASK_HTSEQCOUNT -i $INPUT_HTSEQCOUNT -e .$ASD.bam -n $NODES_HTSEQCOUNT -c $CPU_HTSEQCOUNT -m $MEMORY_HTSEQCOUNT"G" -w $WALLTIME_HTSEQCOUNT \
        --command "${NGSANE_BASE}/mods/htseqcount.sh -k $CONFIG -f <FILE> -o $OUT/<DIR>/$TASK_HTSEQCOUNT/<NAME>" \
        --postcommand "${NGSANE_BASE}/mods/countsTable.sh -f <FILE> -k $CONFIG --outdir $OUT/expression/$TASK_HTSEQCOUNT-<DIR>"
fi


################################################################################
#  Gene expression analysis with tophat + cufflinks + htseqcount
#
# IN : $SOURCE/$dir/fastq/*read1.fastq
# OUT: $OUT/$dir/tophat/*.bam/
################################################################################       

if [ -n "$RUNTOPHATCUFFHTSEQ" ]; then
    if [ -z "$TASK_TOPHAT" ] || [ -z "$NODES_TOPHAT" ] || [ -z "$CPU_TOPHAT" ] || [ -z "$MEMORY_TOPHAT" ] || [ -z "$WALLTIME_TOPHAT" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi

    if [ -z "$TASK_CUFFLINKS" ] || [ -z "$NODES_CUFFLINKS" ] || [ -z "$CPU_CUFFLINKS" ] || [ -z "$MEMORY_CUFFLINKS" ] || [ -z "$WALLTIME_CUFFLINKS" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi

    if [ -z "$TASK_HTSEQCOUNT" ] || [ -z "$NODES_HTSEQCOUNT" ] || [ -z "$CPU_HTSEQCOUNT" ] || [ -z "$MEMORY_HTSEQCOUNT" ] || [ -z "$WALLTIME_HTSEQCOUNT" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi

    if [ -z "$TASK_BAMANN" ] || [ -z "$NODES_BAMANN" ] || [ -z "$CPU_BAMANN" ] || [ -z "$MEMORY_BAMANN" ] || [ -z "$WALLTIME_BAMANN" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi


	# check that if runnow has been altered this is done for all successor jobs of tophat
	if [ -e $QOUT/$TASK_TOPHAT/runnow.tmp ]; then
		echo "[NOTE] checking runnow.tmp"
		LENGTHT=$(wc -l $QOUT/$TASK_TOPHAT/runnow.tmp | cut -d " " -f 1)
		LENGTHC=$(wc -l $QOUT/$TASK_CUFFLINKS/runnow.tmp | cut -d " " -f 1)
		LENGTHH=$(wc -l $QOUT/$TASK_HTSEQCOUNT/runnow.tmp | cut -d " " -f 1)
		LENGTHB=$(wc -l $QOUT/$TASK_BAMANN/runnow.tmp  | cut -d " " -f 1)
		if [[ $LENGTHT != $LENGTHC ]] || [[ $LENGTHT != $LENGTHH ]] || [[ $LENGTHT != $LENGTHB ]] ; then 
			echo "[ERROR] runnow.tmp files of tophat has different length to its successor jobs ($LENGTHT, $LENGTHC, $LENGTHH, $LENGTHB)"; exit -1; 
		fi
	fi
    
    if [ ! -f ${FASTA%.*}.1.bt2 ];then
        # submit job for index generation if necessary
        INDEXJOBIDS=$(
            $QSUB $ARMED -k $CONFIG -t $TASK_TOPHAT -i $INPUT_TOPHAT -e $READONE.$FASTQ -n $NODES_TOPHAT -c $CPU_TOPHAT -m $MEMORY_TOPHAT"G" \
            -w $WALLTIME_TOPHAT --commontask indexGenome \
            --command "${NGSANE_BASE}/mods/bowtie2Index.sh -k $CONFIG"
        ) && echo -e "$INDEXJOBIDS"
        INDEXJOBIDS=$(waitForJobIds "$INDEXJOBIDS")
    else
        INDEXJOBIDS=""
    fi
    
    JOBIDS=$( 
        $QSUB $ARMED -k $CONFIG -t $TASK_TOPHAT -i $INPUT_TOPHAT -e $READONE.$FASTQ -n $NODES_TOPHAT -c $CPU_TOPHAT -m $MEMORY_TOPHAT"G" $INDEXJOBIDS \
            -w $WALLTIME_TOPHAT \
            --command "${NGSANE_BASE}/mods/tophat.sh -k $CONFIG -f <FILE> -o $OUT/<DIR>/$TASK_TOPHAT/<NAME>" 
    ) && echo -e "$JOBIDS"

   JOBIDS=$(waitForJobIds "$JOBIDS")

    # instruct to wait for TOPHATJOBIDS to finish
    $QSUB $ARMED -r -k $CONFIG -t $TASK_CUFFLINKS -i $INPUT_CUFFLINKS -e .$ASD.bam -n $NODES_CUFFLINKS -c $CPU_CUFFLINKS -m $MEMORY_CUFFLINKS"G" -w $WALLTIME_CUFFLINKS $JOBIDS \
        --command "${NGSANE_BASE}/mods/cufflinks.sh -k $CONFIG -f <FILE> -o $OUT/<DIR>/$TASK_CUFFLINKS/<NAME>"

    $QSUB $ARMED -r -k $CONFIG -t $TASK_HTSEQCOUNT -i $INPUT_HTSEQCOUNT -e .$ASD.bam -n $NODES_HTSEQCOUNT -c $CPU_HTSEQCOUNT -m $MEMORY_HTSEQCOUNT"G" -w $WALLTIME_HTSEQCOUNT $JOBIDS \
        --command "${NGSANE_BASE}/mods/htseqcount.sh -k $CONFIG -f <FILE> -o $OUT/<DIR>/$TASK_HTSEQCOUNT/<NAME>"      
    
    $QSUB --nodir -r $ARMED -k $CONFIG -t $TASK_BAMANN -i $INPUT_BAMANN -e .$ASD.bam -n $NODES_BAMANN -c $CPU_BAMANN -m $MEMORY_BAMANN'G' -w $WALLTIME_BAMANN $JOBIDS \
        --command "${NGSANE_BASE}/mods/annotateBam.sh -k $CONFIG -f <FILE>"

fi

################################################################################
#   Create Bigwig from Bam
#
# IN:$SOURCE/$dir/bowtie/*.bam
# OUT: $OUT/$dir/bowtie/*.bw
################################################################################

if [ -n "$RUNBIGWIG" ]; then
    if [ -z "$TASK_BIGWIG" ] || [ -z "$NODES_BIGWIG" ] || [ -z "$CPU_BIGWIG" ] || [ -z "$MEMORY_BIGWIG" ] || [ -z "$WALLTIME_BIGWIG" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi
    
    $QSUB $ARMED --nodir -r -k $CONFIG -t $TASK_BIGWIG -i $INPUT_BIGWIG -e .$ASD.bam -n $NODES_BIGWIG -c $CPU_BIGWIG -m $MEMORY_BIGWIG"G" -w $WALLTIME_BIGWIG \
        --command "${NGSANE_BASE}/mods/bigwig.sh -k $CONFIG -f <FILE> -o $OUT/<DIR>/$INPUT_BIGWIG"
            
fi


################################################################################
#  HiC analysis with homer
#
# IN: $SOURCE/$dir/bwa/*.bam
# OUT: $OUT/$dir/homerhic/
################################################################################
if [ -n "$RUNHOMERHIC" ]; then
    if [ -z "$TASK_HOMERHIC" ] || [ -z "$NODES_HOMERHIC" ] || [ -z "$CPU_HOMERHIC" ] || [ -z "$MEMORY_HOMERHIC" ] || [ -z "$WALLTIME_HOMERHIC" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi
    
    if [ -z "$POOLED_DATA_NAME" ]; then 
        # don't pool data
        $QSUB $ARMED -r -k $CONFIG -t $TASK_HOMERHIC -i $INPUT_HOMERHIC -e $READONE.$ASD.bam -n $NODES_HOMERHIC -c $CPU_HOMERHIC -m $MEMORY_HOMERHIC"G" -w $WALLTIME_HOMERHIC \
	   --command "${NGSANE_BASE}/mods/hicHomer.sh -k $CONFIG -f <FILE> -o $OUT/<DIR>/$TASK_HOMERHIC"
    else
        # pool data
        $QSUB $ARMED --postname $POOLED_DATA_NAME -r -k $CONFIG -t $TASK_HOMERHIC -i $INPUT_HOMERHIC -e $READONE.$ASD.bam -n $NODES_HOMERHIC -c $CPU_HOMERHIC -m $MEMORY_HOMERHIC"G" -w $WALLTIME_HOMERHIC \
	   --postcommand "${NGSANE_BASE}/mods/hicHomer.sh -k $CONFIG -f <FILE> -o $OUT/<DIR>/$TASK_HOMERHIC"
   
    fi
fi


################################################################################
#  ChIP-seq QC with CHANCE
#
# IN: $SOURCE/$dir/bowtie/*.bam
# OUT: $OUT/$dir/chance/
################################################################################
if [ -n "$RUNCHANCE" ]; then
    if [ -z "$TASK_CHANCE" ] || [ -z "$NODES_CHANCE" ] || [ -z "$CPU_CHANCE" ] || [ -z "$MEMORY_CHANCE" ] || [ -z "$WALLTIME_CHANCE" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi

    $QSUB $ARMED -r -k $CONFIG -t $TASK_CHANCE -i $INPUT_CHANCE -e .$ASD.bam -n $NODES_CHANCE -c $CPU_CHANCE -m $MEMORY_CHANCE"G" -w $WALLTIME_CHANCE \
	   --command "${NGSANE_BASE}/mods/chance.sh -k $CONFIG -f <FILE> -o $OUT/<DIR>/$TASK_CHANCE"
fi


################################################################################
#  ChIP-seq analysis with homer
#
# IN: $SOURCE/$dir/bowtie/*.bam
# OUT: $OUT/$dir/homerchipseq/
################################################################################
if [ -n "$RUNHOMERCHIPSEQ" ]; then
    if [ -z "$TASK_HOMERCHIPSEQ" ] || [ -z "$NODES_HOMERCHIPSEQ" ] || [ -z "$CPU_HOMERCHIPSEQ" ] || [ -z "$MEMORY_HOMERCHIPSEQ" ] || [ -z "$WALLTIME_HOMERCHIPSEQ" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi

    if  [ -n "$CHIPINPUT" ];then
        JOBIDS=$( 
        $QSUB $ARMED -r -k $CONFIG -t $TASK_HOMERCHIPSEQ -i $INPUT_HOMERCHIPSEQ -e .$ASD.bam -n $NODES_HOMERCHIPSEQ -c $CPU_HOMERCHIPSEQ \
        	-m $MEMORY_HOMERCHIPSEQ"G" -w $WALLTIME_HOMERCHIPSEQ --commontask ${CONFIG##*/} \
            --command "${NGSANE_BASE}/mods/chipseqHomerInput.sh -k $CONFIG" 
        ) && echo -e "$JOBIDS"
        JOBIDS=$(waitForJobIds "$JOBIDS")

    else
        JOBIDS=""
    fi
    
    $QSUB $ARMED -r -k $CONFIG -t $TASK_HOMERCHIPSEQ -i $INPUT_HOMERCHIPSEQ -e .$ASD.bam -n $NODES_HOMERCHIPSEQ -c $CPU_HOMERCHIPSEQ -m \
        $MEMORY_HOMERCHIPSEQ"G" -w $WALLTIME_HOMERCHIPSEQ $JOBIDS \
        --command "${NGSANE_BASE}/mods/chipseqHomer.sh -k $CONFIG -f <FILE> -o $OUT/<DIR>/$TASK_HOMERCHIPSEQ"
fi

################################################################################
#  ChIP-seq analysis with peakranger
#
# IN: $SOURCE/$dir/bowtie/*.bam
# OUT: $OUT/$dir/peakranger/
################################################################################
if [ -n "$RUNPEAKRANGER" ]; then
    if [ -z "$TASK_PEAKRANGER" ] || [ -z "$NODES_PEAKRANGER" ] || [ -z "$CPU_PEAKRANGER" ] || [ -z "$MEMORY_PEAKRANGER" ] || [ -z "$WALLTIME_PEAKRANGER" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi

    $QSUB $ARMED -r -k $CONFIG -t $TASK_PEAKRANGER -i $INPUT_PEAKRANGER -e .$ASD.bam -n $NODES_PEAKRANGER -c $CPU_PEAKRANGER -m $MEMORY_PEAKRANGER"G" -w $WALLTIME_PEAKRANGER \
	   --command "${NGSANE_BASE}/mods/peakranger.sh -k $CONFIG -f <FILE> -o $OUT/<DIR>/$TASK_PEAKRANGER"
fi

################################################################################
#  ChIP-seq analysis with MACS2
#
# IN: $SOURCE/$dir/bowtie/*.bam
# OUT: $OUT/$dir/macs2/
################################################################################
if [ -n "$RUNMACS2" ]; then
    if [ -z "$TASK_MACS2" ] || [ -z "$NODES_MACS2" ] || [ -z "$CPU_MACS2" ] || [ -z "$MEMORY_MACS2" ] || [ -z "$WALLTIME_MACS2" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi

    $QSUB $ARMED -r -k $CONFIG -t $TASK_MACS2 -i $INPUT_MACS2 -e .$ASD.bam -n $NODES_MACS2 -c $CPU_MACS2 -m $MEMORY_MACS2"G" -w $WALLTIME_MACS2 \
	   --command "${NGSANE_BASE}/mods/macs2.sh -k $CONFIG -f <FILE> -o $OUT/<DIR>/$TASK_MACS2"
fi

################################################################################
#  De-novo motif discovery with memechip
#
# IN: $SOURCE/$dir/peakranger/*.Bedford
# OUT: $OUT/$dir/memechip/
################################################################################
if [ -n "$RUNMEMECHIP" ]; then
    if [ -z "$TASK_MEMECHIP" ] || [ -z "$NODES_MEMECHIP" ] || [ -z "$CPU_MEMECHIP" ] || [ -z "$MEMORY_MEMECHIP" ] || [ -z "$WALLTIME_MEMECHIP" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi

    $QSUB $ARMED -r -k $CONFIG -t $TASK_MEMECHIP -i $INPUT_MEMECHIP -e $BED -n $NODES_MEMECHIP -c $CPU_MEMECHIP -m $MEMORY_MEMECHIP"G" -w $WALLTIME_MEMECHIP \
	   --command "${NGSANE_BASE}/mods/memechip.sh -k $CONFIG -f <FILE> -o $OUT/<DIR>/$TASK_MEMECHIP"
fi


################################################################################
#  Creating normalized (wig) files with wiggler
#
# IN: $SOURCE/$dir/bowtie/*.bam
# OUT: $OUT/$dir/wiggler/
################################################################################
if [ -n "$RUNWIGGLER" ]; then
    if [ -z "$TASK_WIGGLER" ] || [ -z "$NODES_WIGGLER" ] || [ -z "$CPU_WIGGLER" ] || [ -z "$MEMORY_WIGGLER" ] || [ -z "$WALLTIME_WIGGLER" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi

    $QSUB $ARMED -r -k $CONFIG -t $TASK_WIGGLER -i $INPUT_WIGGLER -e .$ASD.bam -n $NODES_WIGGLER -c $CPU_WIGGLER -m $MEMORY_WIGGLER"G" -w $WALLTIME_WIGGLER \
        --command "${NGSANE_BASE}/mods/wiggler.sh -k $CONFIG -f <FILE> -o $OUT/<DIR>/$TASK_WIGGLER" 
fi

################################################################################ 
#   recalibrate quality scores
#   http://www.broadinstitute.org/gsa/wiki/index.php/Base_quality_score_recalibration
#   http://www.broadinstitute.org/gsa/wiki/index.php/Local_realignment_around_indels
#   http://picard.sourceforge.net/command-line-overview.shtml#FixMateInformation
#   full pipe: http://www.broadinstitute.org/gsa/wiki/index.php/Whole_genome,_deep_coverage
# IN:$SOURCE/$dir/fastq/*$READONE.fastq
# OUT: $OUT/$dir/reCal/*.$ASR.bam
################################################################################

if [ -n "$RUNREALRECAL" ]; then

    $QSUB $ARMED -r -k $CONFIG -t $TASK_RECAL -i $INPUT_REALRECAL -e .$ASD.bam \
        -n $NODES_RECAL -c $CPU_RECAL -m $MEMORY_RECAL"G" -w $WALLTIME_RECAL \
        --command "${NGSANE_BASE}/mods/reCalAln.sh -k $CONFIG -f <FILE> -r $FASTA -d $DBSNPVCF -o $OUT/<DIR>/$TASK_RECAL"

fi

################################################################################
#   Pool bam files (e.g. replicates)
#
# IN : $SOURCE/TASK_BOWTIE/PATTERN*$ASD.bam
# OUT: $OUT/TASK_BOWTIE/_pooled*$ASD.bam
################################################################################

if [ -n "$RUNPOOLBAMS" ]; then
    if [ -z "$TASK_POOLBAMS" ] || [ -z "$NODES_POOLBAMS" ] || [ -z "$CPU_POOLBAMS" ] || [ -z "$MEMORY_POOLBAMS" ] || [ -z "$WALLTIME_POOLBAMS" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi
    
    $QSUB $ARMED -r -d -k $CONFIG -t $TASK_POOLBAMS -i $INPUT_POOLBAMS -e .$ASD.bam -n $NODES_POOLBAMS \
    	-c $CPU_POOLBAMS -m $MEMORY_POOLBAMS"G" -w $WALLTIME_POOLBAMS \
    	--postcommand "${NGSANE_BASE}/mods/poolBams.sh -k $CONFIG" 
fi

################################################################################
# downsampling
################################################################################
if [ -n "$RUNDOWNSAMPLING" ]; then

    echo -e "********** Downsampling"
    let ABB=$READNUMBER/100000
    TASK_DOWNSAMPLE="sample"$ABB"K"
    for dir in ${DIR[@]}; do

      if [ ! -d $QOUT/$TASK_DOWNSAMPLE ]; then mkdir -p $QOUT/$TASK_DOWNSAMPLE; fi
      if [ ! -d $OUT/$dir/$TASK_DOWNSAMPLE ]; then mkdir -p $OUT/$dir/$TASK_DOWNSAMPLE; fi


      for f in $( ls $OUT/$dir/$TASK_RECAL/*$ASR.bam ); do
          n=`basename $f`
          NAME=$dir"_"$n
          echo $dir"/"$TASK_RECAL"/"$NAME" -> "$dir/$TASK_DOWNSAMPLE

          if [ -e $QOUT/$TASK_DOWNSAMPLE/$NAME.out ]; then rm $QOUT/$TASK_DOWNSAMPLE/$NAME.out; fi

          if [ -n "$ARMED" ]; then
	  		 qsub $PRIORITY -j y -o $QOUT/$TASK_DOWNSAMPLE/$NAME.out -cwd -b y \
	  		 -N $TASK_DOWNSAMPLE"_"$NAME -l vf=4G \
			${NGSANE_BASE}/mods/downsample.sh -k ${NGSANE_BASE} -i $f -o $OUT/$dir/$TASK_DOWNSAMPLE \
			-s $READNUMBER -r $FASTA

		  fi
      done
  done
fi


################################################################################
# combine flowcell lenes files for recalibration
#
################################################################################

if [ -n "$mergeBWAbams" ]; then
    # make tmp files for first step combining
    ./combined/mergeguide/README
    
    # combine them in lane bam files
    for e in $( ls $OUT/combined/mergeguide/lanes/ ); do
	merge.sh ${NGSANE_BASE} $OUT/combined/mergeguide/lanes/$e $OUT/combined/bwa/ ${e/tmp/$ASD.bam} bam qout/merged/
    done
fi



################################################################################
# combine all into one bamfile
#
################################################################################

if [ -n "$mergeReCalbams" ]; then
    # make tmp files for first step combining
    ls combined/reCalAln/*.bam >combined/mergeguide/combineAll.txt
    
    # combine them in lane bam files
    ${NGSANE_BASE}/mods/merge.sh ${NGSANE_BASE} $OUT/combined/mergeguide/combineAll.txt $OUT/combined/ DISC1_all.bam bam qout/merged/
fi


################################################################################
# downsample
################################################################################

if [ -n "$DOWNSAMPLE" ]
then

    echo -e "********* $TASK_DOWN"

    if [ ! -d $QOUT/$TASK_DOWN ]; then mkdir -p $QOUT/$TASK_DOWN; fi

    for dir in ${DIR[@]}
      do

      for f in $( ls $SOURCE/fastq/$dir/*$READONE.$FASTQ )
	do
	
	n=`basename $f`
	n2=${n/%$READONE.$FASTQ/.$ASD.bam}
	name=${n/%$READONE.$FASTQ/}
	echo -e ">>>>>"$dir$n2

	if [ ! -d $OUT/$dir/$TASK_DOWN ]; then mkdir -p $OUT/$dir/$TASK_DOWN; fi

	# remove old pbs output
	if [ -e $QOUT/$TASK_DOWN/$dir"_"$n2"0.out" ]; then rm $QOUT/$TASK_DOWN/$dir"_"$n2.*; fi

	#check if this is part of the pipe and jobsubmission needs to wait
	if [ -n "$mappingBWA" ]; then HOLD="-hold_jid "$TASK_BWA"_"$dir"_"$name; fi

	#Submit
	if [ -n "$ARMED" ]; then
	    #downsample
	    python ${NGSANE_BASE}/tools/downsample.py -i $OUT/$dir/bwa/$n2 -o $OUT/$dir/$TASK_DOWN/ -t downsample/$dir --region $SEQREG -w 500 -s 500 -q $QOUT/$TASK_DOWN/$dir

	fi

      done
    done

# echo $( ls Hannibal_FC30MEJAAXX/downsample/*500d.bam ) > mergeBamfiles500d.tmp
# qsub -b y -cwd -j y -o qout/$COMBDIR/merged.asd.500d -pe mpich 2 /clusterdata/hiseq_apps/hiSeqInf/mods/merge.sh ${NGSANE_BASE} mergeBamfiles500d.tmp merged merged.$ASD.500d.bam bam

fi

################################################################################
#   call indels with GATK
################################################################################

#TODO run them together with gatkIndelV2.sh

if [ -n "$GATKcallIndelsSeperate" ]
then

    echo -e "********* $TASK_GATKIND"
    if [ ! -d $QOUT/$TASK_GATKIND ]; then mkdir -p $QOUT/$TASK_GATKIND; fi
    if [ -e $QOUT/$TASK_GATKIND/ids.txt ]; then rm $QOUT/$TASK_GATKIND/ids.txt; fi
    if [ -e $QOUT/$TASK_GATKIND/sum.out ]; then rm $QOUT/$TASK_GATKIND/sum.out; fi


    for dir in ${DIR[@]}
      do

      #ensure dirs are there...
      if [ ! -d $OUT/$dir/$TASK_GATKIND ]; then mkdir -p $OUT/$dir/$TASK_GATKIND; fi

      

      #for f in $( ls $SOURCE/$dir/aln2/*$ASR.bam )
      for f in $( ls $SOURCE/fastq/$dir/*$READONE.fastq )
	do
	n=`basename $f`
	n2=${n/%$READONE.$FASTQ/.$ASR.bam}
	name=${n/%$READONE.$FASTQ/}
	echo -e ">>>>>"$dir$n2

	# remove old pbs output
	if [ -e $QOUT/$TASK_GATKIND/$dir"_"$name.out ]; then rm $QOUT/$TASK_GATKIND/$dir"_"$name.*; fi

	#check if this is part of the pipe and jobsubmission needs to wait
	if [ -n "$RUNMAPPINGBWA" ]; then HOLD="-hold_jid "$TASK_BWA"_"$dir"_"$name; fi
	if [ -n "$recalibrateQualScore" ]; then HOLD="-hold_jid "$TASK_RECAL"_"$dir"_"$name; fi

	#Submit
	if [ -n "$ARMED" ]; then
	    qsub $PRIORITY -j y -o $QOUT/$TASK_GATKIND/$dir'_'$name'.out' -cwd -b y \
		-l h_vmem=12G -N $TASK_GATKIND'_'$dir'_'$name $HOLD\
		${NGSANE_BASE}/mods/gatkIndel.sh ${NGSANE_BASE} $OUT/$dir/$TASK_RECAL/$n2 $FASTA $DBSNPVCF \
		$REFSEQROD $OUT/$dir/$TASK_GATKIND
	fi


      done
    done


    #combine into one file
#    qsub -j y -o $QOUT/gatkInd/final.out -cwd -b y -hold_jid `cat $QOUT/gatkInd/ids.txt`\
#	echo -e "merge"
#		merge-vcf `ls $QOUT/dindel/*.dindel.VCF` > genotyping/indels.dindel.vcf
#		merge-vcf A.vcf.gz B.vcf.gz C.vcf.gz | bgzip -c > out.vcf.gz

fi


################################################################################
#   call indels with GATK
################################################################################

#TODO run them together with gatkIndelV2.sh

if [ -n "$GATKcallIndelsCombined" ]
then

    echo -e "********* $TASK_GATKIND"

    if [ ! -d $QOUT/$TASK_GATKIND ]; then mkdir -p $QOUT/$TASK_GATKIND; fi
    if [ -e $QOUT/$TASK_GATKIND/ids.txt ]; then rm $QOUT/$TASK_GATKIND/ids.txt; fi
    if [ ! -d $OUT/genotype ]; then mkdir -p $OUT/genotype; fi

    NAME=$(echo ${DIR[@]}|sed 's/ /_/g')

    for dir in ${DIR[@]};do
      for f in $( ls $SOURCE/fastq/$dir/*$READONE.fastq );do
	n=`basename $f`
	echo $OUT/$dir/$TASK_RECAL/${n/%$READONE.$FASTQ/.$ASR.bam} >> $TASK_GATKIND"bamfiles.tmp"
      done
    done

    # remove old pbs output
    if [ -e $QOUT/$TASK_GATKIND/$NAME.out ]; then rm $QOUT/$TASK_GATKIND/$NAME.out ; fi

    #check if this is part of the pipe and jobsubmission needs to wait
    if [ -n "$recalibrateQualScore" ]; then 
	HOLD="-hold_jid "
	for dir in ${DIR[@]};do
	    HOLD=$HOLD","$TASK_RECAL"_"$dir"*"
        done
	HOLD=${HOLD/,/}
    fi

    #Submit
    if [ -n "$ARMED" ]; then
	qsub -j y -o $QOUT/$TASK_GATKIND/$NAME.out -cwd -b y \
	    -l mem_free=20G -l h_vmem=20G -N $TASK_GATKIND"_"$NAME.out $HOLD\
	    ${NGSANE_BASE}/mods/gatkIndelV2.sh ${NGSANE_BASE} $TASK_GATKIND"bamfiles.tmp" $FASTA $DBSNPVCF \
	    $REFSEQROD $OUT/genotype $SEQREG
    fi

    
fi

########
# run differental expression detection
########
if [ -n "$RUNANNOTATION" ]; then

    # ensure directory is there
    if [ ! -d $QOUT/$TASK_ANNOVAR ]; then mkdir -p $QOUT/$TASK_ANNOVAR; fi
    if [ ! -d $OUT/$TASK_ANNOVAR ]; then mkdir -p $OUT/$TASK_ANNOVAR; fi
    
    #    for dir in $( ls -d $TASK_VAR/* ); do
    for d in ${DIR[@]}; do	
    
        dir=$OUT/$TASK_GATKVAR/$d
        
        n=`basename $dir`
        echo $n
        
        namesnp=$OUT/$TASK_GATKVAR/$n/$n.filter.snps.vcf
        namesnp2=$OUT/$TASK_GATKVAR/$n/$n.recalfilt.snps.vcf
        nameindel=$OUT/$TASK_GATKVAR/$n/$n.filter.indel.vcf
        
        #cleanup old qouts
        if [ -e $QOUT/$TASK_ANNOVAR/$n.out ]; then rm $QOUT/$TASK_ANNOVAR/$n.out; fi
        #ensure subfolder is there
        if [ ! -d $OUT/$TASK_ANNOVAR/$n ]; then mkdir -p $OUT/$TASK_ANNOVAR/$n; fi
        
        
        #check if this is part of the pipe and jobsubmission needs to wait
        if [ -n "$TASK_GATKVAR" ]; then HOLD="-hold_jid "$TASK_GATKVAR"_"$n; fi
        #echo -e "-hold_jid "$TASK_GATKVAR"_"$n
        HOLD="-hold_jid "$TASK_GATKVAR"_"$n
        
            #submit
        if [ -n "$ARMED" ]; then
            qsub $PRIORITY -b y -cwd -j y -o $QOUT/$TASK_ANNOVAR/$n.out \
        	-N $TASK_ANNOVAR'_'$n $HOLD\
        	${NGSANE_BASE}/mods/annovar.sh -k ${NGSANE_BASE} -i1 $namesnp -i2 $namesnp2 -i3 $nameindel -r $FASTA \
        	-o $OUT/$TASK_ANNOVAR/$n 
        fi
    done

fi


########
# run differental expression detection
########
if [ -n "$RUNCUFFDIFF" ]; then

    CPUS=24

    # ensure directory is there
    if [ ! -d $QOUT/$TASK_CUFFDIFF ]; then mkdir -p $QOUT/$TASK_CUFFDIFF; fi
    if [ ! -d $OUT/$TASK_CUFFDIFF ]; then mkdir -p $OUT/$TASK_CUFFDIFF; fi

    for dir in ${DIR[@]}; do
      for f in $( ls $SOURCE/fastq/$dir/*$READONE.$FASTQ );	do
	n=`basename $f`
	name=$name${n/%$READONE.$FASTQ/}","
      done
    done
    name=$(echo $name | sed 's/\(.*\),/\1/')
    n=${name/,/_u_}

    #cleanup old qouts
    if [ -e $QOUT/$TASK_CUFFDIFF/run.out ]; then rm $QOUT/$TASK_CUFFDIFF/run.out; fi

    #check if this is part of the pipe and jobsubmission needs to wait
    if [ -n "$RUNTOPHATCUFF" ]; then HOLD="-hold_jid "$TASK_TOPHAT"_"$dir"_"$name; fi

    #submit
    if [ -n "$ARMED" ]; then
      qsub $PRIORITY -b y -cwd -j y -o $QOUT/$TASK_CUFFDIFF/$n.out \
     	    -N $TASK_CUFFDIFF'_'$n -pe mpich $CPUS $HOLD\
            ${NGSANE_BASE}/mods/cuffdiff.sh -k ${NGSANE_BASE} -b $name -r $FASTA \
	    -t $CPUS -o $OUT/$TASK_CUFFDIFF/$n -a $REFSEQGTF

	fi

fi

################################################################################
#   Transcriptome assembly without a Reference (Trinity)
# IN : $SOURCE/$dir/fastq/*read1.fastq
# OUT: $OUT/$dir/trinity/*.Trinity.fasta.gz
################################################################################

#  echo -e "        _       _     _ _"
#  echo -e "       | |_ ___|_|___|_| |_ _ _"
#  echo -e "       |  _|  _| |   | |  _| | |"
#  echo -e "       |_| |_| |_|_|_|_|_| |_  |"
#  echo -e "   DeNovo Transcriptome    |___|  "
#  echo -e "   Assembly without a reference genome"
#  echo -e ""

if [ -n "$RUNTRINITY" ]; then
    if [ -z "$NODES_INCHWORM" ] || [ -z "$NCPU_INCHWORM" ] || [ -z "$MEMORY_INCHWORM" ] || [ -z "$WALLTIME_INCHWORM" ] || [ -z "$NODETYPE_INCHWORM" ]; then  echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi
    if [ -z "$NODES_CHRYSALIS" ] || [ -z "$NCPU_CHRYSALIS" ] || [ -z "$MEMORY_CHRYSALIS" ] || [ -z "$WALLTIME_CHRYSALIS" ] || [ -z "$NODETYPE_CHRYSALIS" ] ; then  echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi
    if [ -z "$NODES_BUTTERFLY" ] || [ -z "$NCPU_BUTTERFLY" ] || [ -z "$MEMORY_BUTTERFLY" ] || [ -z "$WALLTIME_BUTTERFLY" ] || [ -z "$NODETYPE_BUTTERFLY" ] ; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi
    
    #  ##########   Inchworm   ###########
    JOBIDS=$( \
        $QSUB $ARMED -k $CONFIG -t $TASK_INCHWORM -i $INPUT_INCHWORM -e $READONE.$FASTQ -n $NODES_INCHWORM \
    	         -c $NCPU_INCHWORM -m $MEMORY_INCHWORM"G" -w $WALLTIME_INCHWORM -q $NODETYPE_INCHWORM \
    	         --command "${NGSANE_BASE}/mods/trinity_inchworm.sh -k $CONFIG -f <FILE> -o $OUT/<DIR>/$TASK_INCHWORM/"\
    ) && echo -e "$JOBIDS" 
	JOBIDS=$(waitForJobIds "$JOBIDS")
    
    #  ##########   Chrysalis  ###########
    JOBIDS=$( \
        $QSUB $ARMED -k $CONFIG -t $TASK_CHRYSALIS -i $INPUT_CHRYSALIS -e $READONE.$FASTQ -n $NODES_CHRYSALIS \
         -c $NCPU_CHRYSALIS -m $MEMORY_CHRYSALIS"G" -w $WALLTIME_CHRYSALIS -q $NODETYPE_CHRYSALIS $JOBIDS \
         --command "${NGSANE_BASE}/mods/trinity_chrysalis.sh -k $CONFIG -f <FILE> -o $OUT/<DIR>/$TASK_CHRYSALIS/"\
    ) && echo -e "$JOBIDS" 
	JOBIDS=$(waitForJobIds "$JOBIDS")

    #  ##########   Butterfly  ###########
    $QSUB $ARMED -k $CONFIG -t $TASK_BUTTERFLY -i $INPUT_BUTTERFLY -e $READONE.$FASTQ -n $NODES_BUTTERFLY \
          -c $NCPU_BUTTERFLY -m $MEMORY_BUTTERFLY"G" -w $WALLTIME_BUTTERFLY -q $NODETYPE_BUTTERFLY $JOBIDS \
          --command "${NGSANE_BASE}/mods/trinity_butterfly.sh -k $CONFIG -f <FILE> -o $OUT/<DIR>/$TASK_BUTTERFLY/"
 
fi
################################################################################
# individual calls
#  ##########   Inchworm   ###########
if [ -n "$RUNINCHWORM" ] && [ -z "$RUNTRINITY" ]; then
    if [ -z "$NODES_INCHWORM" ] || [ -z "$NCPU_INCHWORM" ] || [ -z "$MEMORY_INCHWORM" ] || [ -z "$WALLTIME_INCHWORM" ] || [ -z "$NODETYPE_INCHWORM" ]; then  echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi

    $QSUB $ARMED -k $CONFIG -t $TASK_INCHWORM -i $INPUT_INCHWORM -e $READONE.$FASTQ -n $NODES_INCHWORM \
        -c $NCPU_INCHWORM -m $MEMORY_INCHWORM"G" -w $WALLTIME_INCHWORM -q $NODETYPE_INCHWORM \
        --command "${NGSANE_BASE}/mods/trinity_inchworm.sh -k $CONFIG -f <FILE> -o $OUT/<DIR>/$TASK_INCHWORM"
fi
#  ##########   Chrysalis  ###########
if [ -n "$RUNCHRYSALIS" ] && [ -z "$RUNTRINITY" ]; then
    if [ -z "$NODES_CHRYSALIS" ] || [ -z "$NCPU_CHRYSALIS" ] || [ -z "$MEMORY_CHRYSALIS" ] || [ -z "$WALLTIME_CHRYSALIS" ] || [ -z "$NODETYPE_CHRYSALIS" ] ; then  echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi
    
    $QSUB $ARMED -k $CONFIG -t $TASK_CHRYSALIS -i $INPUT_CHRYSALIS -e $READONE.$FASTQ -n $NODES_CHRYSALIS \
        -c $NCPU_CHRYSALIS -m $MEMORY_CHRYSALIS"G" -w $WALLTIME_CHRYSALIS -q $NODETYPE_CHRYSALIS $JOBIDS \
        --command "${NGSANE_BASE}/mods/trinity_chrysalis.sh -k $CONFIG -f <FILE> -o $OUT/<DIR>/$TASK_CHRYSALIS"
fi
#  ##########   Butterfly  ###########
if [ -n "$RUNBUTTERFLY" ] && [ -z "$RUNTRINITY" ]; then
    if [ -z "$NODES_BUTTERFLY" ] || [ -z "$NCPU_BUTTERFLY" ] || [ -z "$MEMORY_BUTTERFLY" ] || [ -z "$WALLTIME_BUTTERFLY" ] || [ -z "$NODETYPE_BUTTERFLY" ] ; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi
    
    $QSUB $ARMED -k $CONFIG -t $TASK_BUTTERFLY -i $INPUT_BUTTERFLY -e $READONE.$FASTQ -n $NODES_BUTTERFLY \
          -c $NCPU_BUTTERFLY -m $MEMORY_BUTTERFLY"G" -w $WALLTIME_BUTTERFLY -q $NODETYPE_BUTTERFLY $JOBIDS \
          --command "${NGSANE_BASE}/mods/trinity_butterfly.sh -k $CONFIG -f <FILE> -o $OUT/<DIR>/$TASK_BUTTERFLY" 
fi

################################################################################ 
#   Pindel
################################################################################

if [ -n "$RUNPINDEL" ]; then

    $QSUB $ARMED -r -k $CONFIG -t $INPUT_PINDEL-$TASK_PINDEL -i $INPUT_PINDEL -e .$ASD.bam \
        -n $NODES_PINDEL -c $CPU_PINDEL -m $MEMORY_PINDEL"G" -w $WALLTIME_PINDEL \
		--postnodes $NODES_VARCOLLECT --postcpu $CPU_VARCOLLECT \
		--postwalltime $WALLTIME_VARCOLLECT --postmemory $MEMORY_VARCOLLECT"G" \
        --command "${NGSANE_BASE}/mods/pindel.sh -k $CONFIG -f <FILE> -o $OUT/<DIR>/$INPUT_PINDEL-$TASK_PINDEL" \
		--postcommand "${NGSANE_BASE}/mods/variantcollect.sh -k $CONFIG -f <FILE> -i1 $INPUT_PINDEL \
				-i2 ${INPUT_PINDEL}-$TASK_PINDEL -o $OUT/variant/${INPUT_PINDEL}-${TASK_PINDEL}-<DIR> "

fi

################################################################################
#   Bigwig generation using fseq
#
# IN:$SOURCE/$dir/$INPUT_FSEQ/*asd.bam
# OUT: $OUT/$dir/fseq/*.bw
################################################################################

if [ -n "$RUNFSEQ" ]; then
    if [ -z "$TASK_FSEQ" ] || [ -z "$NODES_FSEQ" ] || [ -z "$CPU_FSEQ" ] || [ -z "$MEMORY_FSEQ" ] || [ -z "$WALLTIME_FSEQ" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi

    $QSUB $ARMED -r -k $CONFIG -t $TASK_FSEQ -i $INPUT_FSEQ -e .$ASD.bam -n $NODES_FSEQ -c $CPU_FSEQ -m $MEMORY_FSEQ"G" -w $WALLTIME_FSEQ \
        --command "${NGSANE_BASE}/mods/fseq.sh -k $CONFIG -f <FILE> -o $OUT/<DIR>/$TASK_FSEQ"
fi

################################################################################ 
# FEATURE ANNOTATE
################################################################################ 
if [ -n "$RUNANNOTATINGFEATURE" ]; then
    if [ -z "$TASK_FEATANN" ] || [ -z "$NODES_FEATANN" ] || [ -z "$CPU_FEATANN" ] || [ -z "$MEMORY_FEATANN" ] || [ -z "$WALLTIME_FEATANN" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi
    $QSUB --nodir -r $ARMED -k $CONFIG -t ${INPUT_FEATANN}-${TASK_FEATANN} -i $INPUT_FEATANN -e $ENDING \
    -n $NODES_FEATANN -c $CPU_FEATANN -m $MEMORY_FEATANN'G' -w $WALLTIME_FEATANN --postname postcommand-$UPSTREAM+$DOWNSTREAM \
        --postcommand "${NGSANE_BASE}/mods/annotateFeature.sh -k $CONFIG -f <FILE> -o $OUT/${INPUT_FEATANN}-${TASK_FEATANN}-<DIR> "
fi

################################################################################
#  Fusion search tophat
#
# IN : $SOURCE/$dir/fastq/*read1.fastq
# OUT: $OUT/$dir/tophatfusion_out/*
################################################################################       

if [ -n "$RUNFUSION" ]; then
    if [ -z "$TASK_FUSION" ] || [ -z "$NODES_FUSION" ] || [ -z "$CPU_FUSION" ] || [ -z "$MEMORY_FUSION" ] || [ -z "$WALLTIME_FUSION" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi

    $QSUB $ARMED -k $CONFIG -t $TASK_FUSION -i $INPUT_FUSION -e $READONE.$FASTQ -n $NODES_FUSION -c $CPU_FUSION -m $MEMORY_FUSION"G" -w $WALLTIME_FUSION$INDEXJOBIDS \
        --command "${NGSANE_BASE}/mods/fusion.sh -k $CONFIG -f <FILE> -o $OUT/<DIR>/$TASK_FUSION/<NAME>"

fi

