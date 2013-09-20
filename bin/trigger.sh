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
  direct     run task directly (e.g. on node after qrsh)
  postonly   run only the post analysis steps of a task (if available)
  recover    pick up unfinished business (interrupted jobs)
  html       checks logfiles for errors and creates summary HTML page
  report     alias for html

other options:
  -h         print this help message.
  -v         print version number of ngsane
"
exit
}

function version {
    [ -z $(hash git) ] && NGSANE_VERSION=`which trigger.sh` && cd ${NGSANE_VERSION/bin\/trigger.sh/} && NGSANE_VERSION=`git rev-parse HEAD`
    if [ -z "$NGSANE_VERSION" ] || [[ "$NGSANE_VERSION" == *"fatal: Not a git repository"* ]]; then
        NGSANE_VERSION="release v0.2.0.0"
    else
        NGSANE_VERSION="$NGSANE_VERSION (git hash)"
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
        JOBIDS=$(echo $JOBIDS | cut -d " " -f 2 | tr '\n' ':' | sed 's/:$//g' )
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
    if [ ! -d $OUT/$dir ]; then mkdir -p $OUT/$dir; fi
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
# OUT: $OUT/runstats/fastQC/*
################################################################################

if [ -n "$RUNFASTQC" ]; then
    if [ -z "$TASKFASTQC" ] || [ -z "$NODES_FASTQC" ] || [ -z "$CPU_FASTQC" ] || [ -z "$MEMORY_FASTQC" ] || [ -z "$WALLTIME_FASTQC" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi
    
    $QSUB $ARMED -d -k $CONFIG -t $TASKFASTQC -i $INPUT_FASTQC -e $READONE.$FASTQ -n $NODES_FASTQC \
    	-c $CPU_FASTQC -m $MEMORY_FASTQC"G" -w $WALLTIME_FASTQC \
    	--postcommand "${NGSANE_BASE}/mods/fastQC.sh -k $CONFIG" 
fi

################################################################################
#   FASTQSCREEN 
#
# IN : $SOURCE/fastq/$dir/*read1.fastq
# OUT: $OUT/$dir/fastqscreen/*.
################################################################################

if [ -n "$RUNFASTQSCREEN" ]; then
    $QSUB $ARMED -d -k $CONFIG -t $TASKFASTQSCREEN -i $INPUT_FASTQSCREEN -e $READONE.$FASTQ -n $NODES_FASTQSCREEN \
        -c $CPU_FASTQSCREEN -m $MEMORY_FASTQSCREEN"G" -w $WALLTIME_FASTQSCREEN \
        --command "$NGSANE_BASE/mods/fastqscreen.sh -k $CONFIG -f <FILE>  -o $OUT/<DIR>/$TASKFASTQSCREEN"
fi

################################################################################
#   Blue
#
# IN : $SOURCE/fastq/$dir/*read1.fastq
# OUT: $OUT/$dir_healed/*.
################################################################################

if [ -n "$RUNBLUE" ]; then
    if [ -z "$TASKBLUE" ] || [ -z "$NODES_BLUE" ] || [ -z "$CPU_BLUE" ] || [ -z "$MEMORY_BLUE" ] || [ -z "$WALLTIME_BLUE" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi
    
    $QSUB $ARMED -d -k $CONFIG -t $TASKBLUE -i $INPUT_BLUE -e $READONE.$FASTQ -n $NODES_BLUE \
	   -c $CPU_BLUE -m $MEMORY_BLUE"G" -w $WALLTIME_BLUE \
	   --command "${NGSANE_BASE}/mods/blue.sh -k $CONFIG -f <FILE> -o $INPUT_BLUE/<DIR>_$TASKBLUE/" 

fi




################################################################################
#   CUTADAPT remove contaminants
#
# IN : $SOURCE/fastq/$dir/*read1.fastq
# OUT: $SOURCE/fastq/$dir_cutadapt/*read1.fastq
################################################################################

if [ -n "$RUNCUTADAPT" ]; then
    if [ -z "$TASKCUTADAPT" ] || [ -z "$NODES_CUTADAPT" ] || [ -z "$CPU_CUTADAPT" ] || [ -z "$MEMORY_CUTADAPT" ] || [ -z "$WALLTIME_CUTADAPT" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi
    
    $QSUB $ARMED -d -k $CONFIG -t $TASKCUTADAPT -i $INPUT_CUTADAPT -e $READONE.$FASTQ -n $NODES_CUTADAPT \
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
    if [ -z "$TASKTRIMMOMATIC" ] || [ -z "$NODES_TRIMMOMATIC" ] || [ -z "$CPU_TRIMMOMATIC" ] || [ -z "$MEMORY_TRIMMOMATIC" ] || [ -z "$WALLTIME_TRIMMOMATIC" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi
    
    $QSUB $ARMED -d -k $CONFIG -t $TASKTRIMMOMATIC -i $INPUT_TRIMMOMATIC -e $READONE.$FASTQ -n $NODES_TRIMMOMATIC \
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
    if [ -z "$TASKTRIMGALORE" ] || [ -z "$NODES_TRIMGALORE" ] || [ -z "$CPU_TRIMGALORE" ] || [ -z "$MEMORY_TRIMGALORE" ] || [ -z "$WALLTIME_TRIMGALORE" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi
    
    $QSUB $ARMED -d -k $CONFIG -t $TASKTRIMGALORE -i $INPUT_TRIMGALORE -e $READONE.$FASTQ -n $NODES_TRIMGALORE \
        -c $CPU_TRIMGALORE -m $MEMORY_TRIMGALORE"G" -w $WALLTIME_TRIMGALORE \
        --command "${NGSANE_BASE}/mods/trimgalore.sh -k $CONFIG -f <FILE>"
fi

################################################################################ 
# IN : */bwa/*.bam
# OUT: */bwa/*.ann
################################################################################ 
if [ -n "$RUNANNOTATINGBAM" ]; then
    if [ -z "$TASKBAMANN" ] || [ -z "$NODES_BAMANN" ] || [ -z "$CPU_BAMANN" ] || [ -z "$MEMORY_BAMANN" ] || [ -z "$WALLTIME_BAMANN" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi
    
    $QSUB --nodir -r $ARMED -k $CONFIG -t $TASKBAMANN -i $INPUT_BAMANN -e .$ASD.bam \
    -n $NODES_BAMANN -c $CPU_BAMANN -m $MEMORY_BAMANN'G' -w $WALLTIME_BAMANN \
        --command "${NGSANE_BASE}/mods/annotateBam.sh -k $CONFIG -f <FILE>"
fi


################################################################################
#   Variance calling
# IN: */bwa/*.bam
# OUT: */bwa_var/*.clean.vcf
################################################################################

if [ -n "$RUNSAMVAR" ]; then
    if [ -z "$TASKBWA" ] || [ -z "$NODES_SAMVAR" ] || [ -z "$CPU_SAMVAR" ] || [ -z "$MEMORY_SAMVAR" ] || [ -z "$WALLTIME_SAMVAR" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi
    
    $QSUB -r $ARMED -k $CONFIG -t $TASKBWA-$TASKSAMVAR -i $INPUT_SAMVAR -e .$ASD.bam \
       -n $NODES_SAMVAR -c $CPU_SAMVAR -m $MEMORY_SAMVAR'G' -w $WALLTIME_SAMVAR \
       --command "${NGSANE_BASE}/mods/samSNPs.sh -k $CONFIG -f <FILE> -o $OUT/<DIR>/$TASKBWA-$TASKSAMVAR" \
	   --postcommand "${NGSANE_BASE}/mods/samSNPscollect.sh -k $CONFIG -f <FILE> -o $OUT/variant/$TASKBWA-$TASKSAMVAR-<DIR>"
fi

################################################################################
#   call indels with GATK -- call one vcf file over all folders
################################################################################

if [ -n "$RUNVARCALLS" ]; then
    NAME=$(echo ${DIR[@]}|sed 's/ /_/g')
    $QSUB $ARMED -r -d -k $CONFIG -t $TASKVAR -i $INPUT_VAR  -e .$ASR.bam -n $NODES_VAR \
        -c $CPU_VAR -m $MEMORY_VAR"G" -w $WALLTIME_VAR \
        --postcommand "${NGSANE_BASE}/mods/gatkSNPs2.sh -k $CONFIG \
                        -i <FILE> -t $CPU_VAR \
                        -r $FASTA -d $DBROD -o $OUT/$TASKVAR/$NAME -n $NAME \
                        -H $HAPMAPVCF -K $ONEKGVCF"
fi

################################################################################
#   Depth of Coverage
# IN: */bwa/*.bam
# OUT: */bwa_var/*.clean.vcf
################################################################################

if [ -n "$DEPTHOFCOVERAGE2" ]; then

    $QSUB $ARMED -r -k $CONFIG -t $TASKDOC -i $INPUT_GATKDOC -e .$ASR.bam \
	-n $NODES_GATKDOC -c $CPU_GATKDOC -m $MEMORY_GATKDOC"G" -w $WALLTIME_GATKDOC \
	--command "${NGSANE_BASE}/mods/gatkDOC.sh -k ${NGSANE_BASE} -f <FILE> -r $FASTA -o $OUT/<DIR>/$TASKDOC -t $CPU_GATKDOC"

fi

################################################################################
#   Mapping using HiCUP
# IN : $SOURCE/fastq/$dir/*read1.fastq
# OUT: $OUT/$dir/hicup/*.$ASD.bam
################################################################################

if [ -n "$RUNHICUP" ]; then
    if [ -z "$TASKHICUP" ] || [ -z "$NODES_HICUP" ] || [ -z "$CPU_HICUP" ] || [ -z "$MEMORY_HICUP" ] || [ -z "$WALLTIME_HICUP" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi
    
    $QSUB $ARMED -k $CONFIG -t $TASKHICUP -i $INPUT_HICUP -e $READONE.$FASTQ -n $NODES_HICUP -c $CPU_HICUP \
    	-m $MEMORY_HICUP"G" -w $WALLTIME_HICUP \
        --command "${NGSANE_BASE}/mods/hicup.sh -k $CONFIG -f <FILE> -o $OUT/<DIR>/$TASKHICUP"
fi

################################################################################
#  Assessing HiC data with hiclib
#
# IN: $SOURCE/$dir/fastq/*read1.fastq
# OUT: $OUT/$dir/hiclib/*.hdf5
################################################################################

if [ -n "$RUNHICLIB" ]; then
    if [ -z "$TASKHICLIB" ] || [ -z "$NODES_HICLIB" ] || [ -z "$CPU_HICLIB" ] || [ -z "$MEMORY_HICLIB" ] || [ -z "$WALLTIME_HICLIB" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi

    $QSUB $ARMED -k $CONFIG -t $TASKHICLIB -i $INPUT_HICLIB -e $READONE.$FASTQ \
    	-n $NODES_HICLIB -c $CPU_HICLIB -m $MEMORY_HICLIB"G" -w $WALLTIME_HICLIB \
    	--postnodes $NODES_HICLIB_POSTCOMMAND --postcpu $CPU_HICLIB_POSTCOMMAND \
        --command "${NGSANE_BASE}/mods/hiclibMapping.sh -k $CONFIG --fastq <FILE> --outdir $OUT/<DIR>/$TASKHICLIB --fastqName <NAME>" \
        --postcommand "${NGSANE_BASE}/mods/hiclibCorrelate.sh -f <FILE> -k $CONFIG --outdir $OUT/hiclib/$TASKHICLIB-<DIR>"
fi

################################################################################
#   Mapping using BWA
#
# IN : $SOURCE/$dir/fastq/*read1.fastq
# OUT: $OUT/$dir/bwa/*.$ASD.bam
################################################################################

if [ -n "$RUNMAPPINGBWA" ]; then
    if [ -z "$TASKBWA" ] || [ -z "$NODES_BWA" ] || [ -z "$CPU_BWA" ] || [ -z "$MEMORY_BWA" ] || [ -z "$WALLTIME_BWA" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi

    $QSUB $ARMED -k $CONFIG -t $TASKBWA -i $INPUT_BWA -e $READONE.$FASTQ -n $NODES_BWA -c $CPU_BWA -m $MEMORY_BWA"G" -w $WALLTIME_BWA \
        --command "${NGSANE_BASE}/mods/bwa.sh -k $CONFIG -f <FILE> -o $OUT/<DIR>/$TASKBWA --rgsi <DIR>"
fi

################################################################################
#   Mapping using Bowtie v1
#
# IN:$SOURCE/$dir/fastq/*read1.fastq
# OUT: $OUT/$dir/bowtie/*.bam
################################################################################

if [ -n "$RUNMAPPINGBOWTIE" ]; then
    if [ -z "$TASKBOWTIE" ] || [ -z "$NODES_BOWTIE" ] || [ -z "$CPU_BOWTIE" ] || [ -z "$MEMORY_BOWTIE" ] || [ -z "$WALLTIME_BOWTIE" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi

    $QSUB $ARMED -k $CONFIG -t $TASKBOWTIE -i $INPUT_BOWTIE -e $READONE.$FASTQ -n $NODES_BOWTIE -c $CPU_BOWTIE -m $MEMORY_BOWTIE"G" -w $WALLTIME_BOWTIE \
        --command "${NGSANE_BASE}/mods/bowtie.sh -k $CONFIG -f <FILE> -o $OUT/<DIR>/$TASKBOWTIE --rgsi <DIR>"
            
fi


################################################################################
#  Mapping with bowtie v2
#
# IN:  $SOURCE/$dir/fastq/*read1.fastq
# OUT: $OUT/$dir/bowtie/*.bam
################################################################################
if [ -n "$RUNMAPPINGBOWTIE2" ]; then
    if [ -z "$TASKBOWTIE2" ] || [ -z "$NODES_BOWTIE2" ] || [ -z "$CPU_BOWTIE2" ] || [ -z "$MEMORY_BOWTIE2" ] || [ -z "$WALLTIME_BOWTIE2" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi
    
    $QSUB $ARMED -k $CONFIG -t $TASKBOWTIE2 -i $INPUT_BOWTIE2 -e $READONE.$FASTQ -n $NODES_BOWTIE2 -c $CPU_BOWTIE2 -m $MEMORY_BOWTIE2"G" -w $WALLTIME_BOWTIE2 \
	       --command "${NGSANE_BASE}/mods/bowtie2.sh -k $CONFIG -f <FILE> -o $OUT/<DIR>/$TASKBOWTIE2 --rgsi <DIR>"
      
fi


################################################################################
#   Mapping using rrbsmap
#
# IN:$SOURCE/$dir/fastq/*read1.fastq
# OUT: $OUT/$dir/rrbsmap/*.bam
################################################################################

if [ -n "$RUNMAPPINGRRBS" ]; then
    if [ -z "$TASKRRBSMAP" ] || [ -z "$NODES_RRBSMAP" ] || [ -z "$CPU_RRBSMAP" ] || [ -z "$MEMORY_RRBSMAP" ] || [ -z "$WALLTIME_RRBSMAP" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi

    $QSUB $ARMED -k $CONFIG -t $TASKRRBSMAP -i $INPUT_RRBSMAP -e $READONE.$FASTQ -n $NODES_RRBSMAP -c $CPU_RRBSMAP -m $MEMORY_RRBSMAP"G" -w $WALLTIME_RRBSMAP \
        --command "${NGSANE_BASE}/mods/rrbsmap.sh -k $CONFIG -f <FILE> -r $FASTA \
            -o $OUT/<DIR>/$TASKRRBSMAP --rgid $EXPID --rglb $LIBRARY --rgpl $PLATFORM --rgsi <DIR>"
fi


################################################################################
#  Transcript mapping tophat
#
# IN : $SOURCE/$dir/fastq/*read1.fastq
# OUT: $OUT/$dir/tophat/*.bam/
################################################################################       

if [ -n "$RUNTOPHAT" ]; then
    if [ -z "$TASKTOPHAT" ] || [ -z "$NODES_TOPHAT" ] || [ -z "$CPU_TOPHAT" ] || [ -z "$MEMORY_TOPHAT" ] || [ -z "$WALLTIME_TOPHAT" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi

    $QSUB $ARMED -k $CONFIG -t $TASKTOPHAT -i $INPUT_TOPHAT -e $READONE.$FASTQ -n $NODES_TOPHAT -c $CPU_TOPHAT -m $MEMORY_TOPHAT"G" -w $WALLTIME_TOPHAT \
        --command "${NGSANE_BASE}/mods/tophat.sh -k $CONFIG -f <FILE> -o $OUT/<DIR>/$TASKTOPHAT/<NAME>"

fi

################################################################################
#  Gene expression analysis with cufflinks
#
# IN : $OUT/$dir/tophat/*.bam
# OUT: $OUT/$dir/cufflinks/*_transcript.gtf
################################################################################       

if [ -n "$RUNCUFFLINKS" ]; then
    if [ -z "$TASKCUFFLINKS" ] || [ -z "$NODES_CUFFLINKS" ] || [ -z "$CPU_CUFFLINKS" ] || [ -z "$MEMORY_CUFFLINKS" ] || [ -z "$WALLTIME_CUFFLINKS" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi

    $QSUB $ARMED -r -k $CONFIG -t $TASKCUFFLINKS -i $INPUT_CUFFLINKS -e .$ASD.bam -n $NODES_CUFFLINKS -c $CPU_CUFFLINKS -m $MEMORY_CUFFLINKS"G" -w $WALLTIME_CUFFLINKS  \
        --command "${NGSANE_BASE}/mods/cufflinks.sh -k $CONFIG -f <FILE> -o $OUT/<DIR>/$TASKCUFFLINKS/<NAME>"

fi

################################################################################
#  Feature counting with HTSEQCOUNT
#
# IN : $OUT/$dir/tophat/*.bam
# OUT: $OUT/$dir/htseqcount/*_transcript.gtf
################################################################################       

if [ -n "$RUNHTSEQCOUNT" ]; then
    if [ -z "$TASKHTSEQCOUNT" ] || [ -z "$NODES_HTSEQCOUNT" ] || [ -z "$CPU_HTSEQCOUNT" ] || [ -z "$MEMORY_HTSEQCOUNT" ] || [ -z "$WALLTIME_HTSEQCOUNT" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi

    $QSUB $ARMED -r -k $CONFIG -t $TASKHTSEQCOUNT -i $INPUT_HTSEQCOUNT -e .$ASD.bam -n $NODES_HTSEQCOUNT -c $CPU_HTSEQCOUNT -m $MEMORY_HTSEQCOUNT"G" -w $WALLTIME_HTSEQCOUNT \
        --command "${NGSANE_BASE}/mods/htseqcount.sh -k $CONFIG -f <FILE> -o $OUT/<DIR>/$TASKHTSEQCOUNT/<NAME>"

fi


################################################################################
#  Gene expression analysis with tophat + cufflinks + htseqcount
#
# IN : $SOURCE/$dir/fastq/*read1.fastq
# OUT: $OUT/$dir/tophat/*.bam/
################################################################################       

if [ -n "$RUNTOPHATCUFFHTSEQ" ]; then
    if [ -z "$TASKTOPHAT" ] || [ -z "$NODES_TOPHAT" ] || [ -z "$CPU_TOPHAT" ] || [ -z "$MEMORY_TOPHAT" ] || [ -z "$WALLTIME_TOPHAT" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi

    if [ -z "$TASKCUFFLINKS" ] || [ -z "$NODES_CUFFLINKS" ] || [ -z "$CPU_CUFFLINKS" ] || [ -z "$MEMORY_CUFFLINKS" ] || [ -z "$WALLTIME_CUFFLINKS" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi

    if [ -z "$TASKHTSEQCOUNT" ] || [ -z "$NODES_HTSEQCOUNT" ] || [ -z "$CPU_HTSEQCOUNT" ] || [ -z "$MEMORY_HTSEQCOUNT" ] || [ -z "$WALLTIME_HTSEQCOUNT" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi

    if [ -z "$TASKBAMANN" ] || [ -z "$NODES_BAMANN" ] || [ -z "$CPU_BAMANN" ] || [ -z "$MEMORY_BAMANN" ] || [ -z "$WALLTIME_BAMANN" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi
    
    JOBIDS=$( 
        $QSUB $ARMED -k $CONFIG -t $TASKTOPHAT -i $INPUT_TOPHAT -e $READONE.$FASTQ -n $NODES_TOPHAT -c $CPU_TOPHAT -m $MEMORY_TOPHAT"G" \
            -w $WALLTIME_TOPHAT \
            --command "${NGSANE_BASE}/mods/tophat.sh -k $CONFIG -f <FILE> -o $OUT/<DIR>/$TASKTOPHAT/<NAME>" 
    ) && echo -e "$JOBIDS"

   JOBIDS=$(waitForJobIds "$JOBIDS")

    # instruct to wait for TOPHATJOBIDS to finish
    $QSUB $ARMED -r -k $CONFIG -t $TASKCUFFLINKS -i $INPUT_CUFFLINKS -e .$ASD.bam -n $NODES_CUFFLINKS -c $CPU_CUFFLINKS -m $MEMORY_CUFFLINKS"G" -w $WALLTIME_CUFFLINKS $JOBIDS \
        --command "${NGSANE_BASE}/mods/cufflinks.sh -k $CONFIG -f <FILE> -o $OUT/<DIR>/$TASKCUFFLINKS/<NAME>"

    $QSUB $ARMED -r -k $CONFIG -t $TASKHTSEQCOUNT -i $INPUT_HTSEQCOUNT -e .$ASD.bam -n $NODES_HTSEQCOUNT -c $CPU_HTSEQCOUNT -m $MEMORY_HTSEQCOUNT"G" -w $WALLTIME_HTSEQCOUNT $JOBIDS \
        --command "${NGSANE_BASE}/mods/htseqcount.sh -k $CONFIG -f <FILE> -o $OUT/<DIR>/$TASKHTSEQCOUNT/<NAME>"      
    
    $QSUB --nodir -r $ARMED -k $CONFIG -t $TASKBAMANN -i $INPUT_BAMANN -e .$ASD.bam -n $NODES_BAMANN -c $CPU_BAMANN -m $MEMORY_BAMANN'G' -w $WALLTIME_BAMANN $JOBIDS \
        --command "${NGSANE_BASE}/mods/annotateBam.sh -k $CONFIG -f <FILE>"

fi

################################################################################
#   Create Bigwig from Bam
#
# IN:$SOURCE/$dir/bowtie/*.bam
# OUT: $OUT/$dir/bowtie/*.bw
################################################################################

if [ -n "$RUNBIGWIG" ]; then
    if [ -z "$TASKBIGWIG" ] || [ -z "$NODES_BIGWIG" ] || [ -z "$CPU_BIGWIG" ] || [ -z "$MEMORY_BIGWIG" ] || [ -z "$WALLTIME_BIGWIG" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi
    
    $QSUB $ARMED --nodir -r -k $CONFIG -t $TASKBIGWIG -i $INPUT_BIGWIG -e .$ASD.bam -n $NODES_BIGWIG -c $CPU_BIGWIG -m $MEMORY_BIGWIG"G" -w $WALLTIME_BIGWIG \
        --command "${NGSANE_BASE}/mods/bigwig.sh -k $CONFIG -f <FILE> -o $OUT/<DIR>/$INPUT_BIGWIG"
            
fi


################################################################################
#  HiC analysis with homer
#
# IN: $SOURCE/$dir/bwa/*.bam
# OUT: $OUT/$dir/homerhic/
################################################################################
if [ -n "$RUNHOMERHIC" ]; then
    if [ -z "$TASKHOMERHIC" ] || [ -z "$NODES_HOMERHIC" ] || [ -z "$CPU_HOMERHIC" ] || [ -z "$MEMORY_HOMERHIC" ] || [ -z "$WALLTIME_HOMERHIC" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi
    
    $QSUB $ARMED -r -k $CONFIG -t $TASKHOMERHIC -i $INPUT_HOMERHIC -e .$ASD.bam -n $NODES_HOMERHIC -c $CPU_HOMERHIC -m $MEMORY_HOMERHIC"G" -w $WALLTIME_HOMERHIC \
	   --command "${NGSANE_BASE}/mods/hicHomer.sh -k $CONFIG -f <FILE> -o $OUT/<DIR>/$TASKHOMERHIC"
fi


################################################################################
#  ChIP-seq QC with CHANCE
#
# IN: $SOURCE/$dir/bowtie/*.bam
# OUT: $OUT/$dir/chance/
################################################################################
if [ -n "$RUNCHANCE" ]; then
    if [ -z "$TASKCHANCE" ] || [ -z "$NODES_CHANCE" ] || [ -z "$CPU_CHANCE" ] || [ -z "$MEMORY_CHANCE" ] || [ -z "$WALLTIME_CHANCE" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi

    $QSUB $ARMED -r -k $CONFIG -t $TASKCHANCE -i $INPUT_CHANCE -e .$ASD.bam -n $NODES_CHANCE -c $CPU_CHANCE -m $MEMORY_CHANCE"G" -w $WALLTIME_CHANCE \
	   --command "${NGSANE_BASE}/mods/chance.sh -k $CONFIG -f <FILE> -o $OUT/<DIR>/$TASKCHANCE"
fi


################################################################################
#  ChIP-seq analysis with homer
#
# IN: $SOURCE/$dir/bowtie/*.bam
# OUT: $OUT/$dir/homerchipseq/
################################################################################
if [ -n "$RUNHOMERCHIPSEQ" ]; then
    if [ -z "$TASKHOMERCHIPSEQ" ] || [ -z "$NODES_HOMERCHIPSEQ" ] || [ -z "$CPU_HOMERCHIPSEQ" ] || [ -z "$MEMORY_HOMERCHIPSEQ" ] || [ -z "$WALLTIME_HOMERCHIPSEQ" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi
    
    $QSUB $ARMED -r -k $CONFIG -t $TASKHOMERCHIPSEQ -i $INPUT_HOMERCHIPSEQ -e .$ASD.bam -n $NODES_HOMERCHIPSEQ -c $CPU_HOMERCHIPSEQ -m $MEMORY_HOMERCHIPSEQ"G" -w $WALLTIME_HOMERCHIPSEQ \
	   --command "${NGSANE_BASE}/mods/chipseqHomer.sh -k $CONFIG -f <FILE> -o $OUT/<DIR>/$TASKHOMERCHIPSEQ"
fi

################################################################################
#  ChIP-seq analysis with peakranger
#
# IN: $SOURCE/$dir/bowtie/*.bam
# OUT: $OUT/$dir/peakranger/
################################################################################
if [ -n "$RUNPEAKRANGER" ]; then
    if [ -z "$TASKPEAKRANGER" ] || [ -z "$NODES_PEAKRANGER" ] || [ -z "$CPU_PEAKRANGER" ] || [ -z "$MEMORY_PEAKRANGER" ] || [ -z "$WALLTIME_PEAKRANGER" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi

    $QSUB $ARMED -r -k $CONFIG -t $TASKPEAKRANGER -i $INPUT_PEAKRANGER -e .$ASD.bam -n $NODES_PEAKRANGER -c $CPU_PEAKRANGER -m $MEMORY_PEAKRANGER"G" -w $WALLTIME_PEAKRANGER \
	   --command "${NGSANE_BASE}/mods/peakranger.sh -k $CONFIG -f <FILE> -o $OUT/<DIR>/$TASKPEAKRANGER"
fi

################################################################################
#  ChIP-seq analysis with MACS2
#
# IN: $SOURCE/$dir/bowtie/*.bam
# OUT: $OUT/$dir/macs2/
################################################################################
if [ -n "$RUNMACS2" ]; then
    if [ -z "$TASKMACS2" ] || [ -z "$NODES_MACS2" ] || [ -z "$CPU_MACS2" ] || [ -z "$MEMORY_MACS2" ] || [ -z "$WALLTIME_MACS2" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi

    $QSUB $ARMED -r -k $CONFIG -t $TASKMACS2 -i $INPUT_MACS2 -e .$ASD.bam -n $NODES_MACS2 -c $CPU_MACS2 -m $MEMORY_MACS2"G" -w $WALLTIME_MACS2 \
	   --command "${NGSANE_BASE}/mods/macs2.sh -k $CONFIG -f <FILE> -o $OUT/<DIR>/$TASKMACS2"
fi

################################################################################
#  De-novo motif discovery with memechip
#
# IN: $SOURCE/$dir/peakranger/*.Bedford
# OUT: $OUT/$dir/memechip/
################################################################################
if [ -n "$RUNMEMECHIP" ]; then
    if [ -z "$TASKMEMECHIP" ] || [ -z "$NODES_MEMECHIP" ] || [ -z "$CPU_MEMECHIP" ] || [ -z "$MEMORY_MEMECHIP" ] || [ -z "$WALLTIME_MEMECHIP" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi

    $QSUB $ARMED -r -k $CONFIG -t $TASKMEMECHIP -i $INPUT_MEMECHIP -e $BED -n $NODES_MEMECHIP -c $CPU_MEMECHIP -m $MEMORY_MEMECHIP"G" -w $WALLTIME_MEMECHIP \
	   --command "${NGSANE_BASE}/mods/memechip.sh -k $CONFIG -f <FILE> -o $OUT/<DIR>/$TASKMEMECHIP"
fi


################################################################################
#  Creating normalized (wig) files with wiggler
#
# IN: $SOURCE/<DIR>/bowtie/*.bam
# OUT: $OUT/$dir/wiggler/
################################################################################
if [ -n "$RUNWIGGLER" ]; then
    if [ -z "$TASKWIGGLER" ] || [ -z "$NODES_WIGGLER" ] || [ -z "$CPU_WIGGLER" ] || [ -z "$MEMORY_WIGGLER" ] || [ -z "$WALLTIME_WIGGLER" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi

    $QSUB $ARMED -r -k $CONFIG -t $TASKWIGGLER -i $INPUT_WIGGLER -e .$ASD.bam -n $NODES_WIGGLER -c $CPU_WIGGLER -m $MEMORY_WIGGLER"G" -w $WALLTIME_WIGGLER \
        --postcommand "${NGSANE_BASE}/mods/wiggler.sh -k $CONFIG -f <FILE> -o $OUT/<DIR>/$TASKWIGGLER" 
fi

################################################################################ 
#   recalibrate quality scores
#   http://www.broadinstitute.org/gsa/wiki/index.php/Base_quality_score_recalibration
#   http://www.broadinstitute.org/gsa/wiki/index.php/Local_realignment_around_indels
#   http://picard.sourceforge.net/command-line-overview.shtml#FixMateInformation
#   full pipe: http://www.broadinstitute.org/gsa/wiki/index.php/Whole_genome,_deep_coverage
# IN:$SOURCE/$dir/fastq/*$READONE.fastq
# OUT: $OUT/$dir/reCal/*.$$ASR.bam
################################################################################

if [ -n "$RUNREALRECAL" ]; then

    $QSUB $ARMED -r -k $CONFIG -t $TASKRCA -i $INPUT_REALRECAL -e .$ASD.bam \
        -n $NODES_RECAL -c $CPU_RECAL -m $MEMORY_RECAL"G" -w $WALLTIME_RECAL \
        --command "${NGSANE_BASE}/mods/reCalAln2.sh -k $CONFIG -f <FILE> -r $FASTA -d $DBROD -o $OUT/<DIR>/$TASKRCA"

fi

################################################################################
#   Pool bam files (e.g. replicates)
#
# IN : $SOURCE/TASKBOWTIE/PATTERN*$ASD.bam
# OUT: $OUT/TASKBOWTIE/_pooled*$ASD.bam
################################################################################

if [ -n "$RUNPOOLBAMS" ]; then
    if [ -z "$TASKPOOLBAMS" ] || [ -z "$NODES_POOLBAMS" ] || [ -z "$CPU_POOLBAMS" ] || [ -z "$MEMORY_POOLBAMS" ] || [ -z "$WALLTIME_POOLBAMS" ]; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi
    
    $QSUB $ARMED -r -d -k $CONFIG -t $TASKPOOLBAMS -i $INPUT_POOLBAMS -e .$ASD.bam -n $NODES_POOLBAMS \
    	-c $CPU_POOLBAMS -m $MEMORY_POOLBAMS"G" -w $WALLTIME_POOLBAMS \
    	--postcommand "${NGSANE_BASE}/mods/poolBams.sh -k $CONFIG" 
fi

################################################################################
# downsampling
################################################################################
if [ -n "$RUNDOWNSAMPLING" ]; then

    echo -e "********** Downsampling"
    let ABB=$READNUMBER/100000
    TASKDOWNSAMPLE="sample"$ABB"K"
    for dir in ${DIR[@]}; do

      if [ ! -d $QOUT/$TASKDOWNSAMPLE ]; then mkdir -p $QOUT/$TASKDOWNSAMPLE; fi
      if [ ! -d $OUT/$dir/$TASKDOWNSAMPLE ]; then mkdir -p $OUT/$dir/$TASKDOWNSAMPLE; fi


      for f in $( ls $OUT/$dir/$TASKRCA/*$ASR.bam ); do
          n=`basename $f`
          NAME=$dir"_"$n
          echo $dir"/"$TASKRCA"/"$NAME" -> "$dir/$TASKDOWNSAMPLE

          if [ -e $QOUT/$TASKDOWNSAMPLE/$NAME.out ]; then rm $QOUT/$TASKDOWNSAMPLE/$NAME.out; fi

          if [ -n "$ARMED" ]; then
	  		 qsub $PRIORITY -j y -o $QOUT/$TASKDOWNSAMPLE/$NAME.out -cwd -b y \
	  		 -N $TASKDOWNSAMPLE"_"$NAME -l vf=4G \
			${NGSANE_BASE}/mods/downsample.sh -k ${NGSANE_BASE} -i $f -o $OUT/$dir/$TASKDOWNSAMPLE \
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
#   recalibrate quality scores OLD
#   http://www.broadinstitute.org/gsa/wiki/index.php/Base_quality_score_recalibration
#   http://www.broadinstitute.org/gsa/wiki/index.php/Local_realignment_around_indels
#   http://picard.sourceforge.net/command-line-overview.shtml#FixMateInformation
#   full pipe: http://www.broadinstitute.org/gsa/wiki/index.php/Whole_genome,_deep_coverage
# IN:$SOURCE/$dir/fastq/*$READONE.fastq
# OUT: $OUT/$dir/reCal/*.$$ASR.bam
################################################################################
if [ -n "$recalibrateQualScore" ]; then

    echo -e "********* $TASKRCA"

    # generates the .dict file
    #java -Xmx4g -jar /home/Software/picard-tools-1.22/CreateSequenceDictionary.jar R= $FASTA O= ${$FASTA/.fasta/.dict}
    #sort dbSNP rod file according to fasta
    #/home/Software/Sting/perl/sortByRef.pl --k 2 $GATKSUP/tmp.dbsnp.txt $FASTA.fai > $DBROD

    if [ ! -d $QOUT/$TASKRCA ]; then mkdir -p $QOUT/$TASKRCA; fi
    #ensure dirs are there...
    if [ ! -d $OUT/combined/$TASKRCA ]; then mkdir -p $OUT/combined/$TASKRCA; fi

    for f in $( ls $OUT/combined/bwa/*$ASD.bam )
#    for f in $( ls combined/bwa/Nina_FC3056JAAXX_l6.$ASD.bam )
	do
	n=`basename $f`
	name=${n/.bam/}
	echo -e ">>>>>"$n


	# wait on pipeline steps
	HOLD="-hold_jid *mergeBams*"
	# remove old pbs output
	if [ -e $QOUT/$TASKRCA/$name.out ]; then rm $QOUT/$TASKRCA/$name.*; fi

	#Sumit (ca. 80 min on AV 1297205 reads) -l h_rt=20:00:00
	if [ -n "$ARMED" ]; then
	    qsub $PRIORITY -j y -o $QOUT/$TASKRCA/$name'.out' -cwd -b y -l vf=20G \
		-l mem_free=20G -l h_vmem=20G -N $TASKRCA'_'$name $HOLD\
		${NGSANE_BASE}/mods/reCalAln.sh ${NGSANE_BASE} $f $FASTA $DBROD \
		$OUT/combined/$TASKRCA $SEQREG
	fi
	
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
# DepthOfCoverage
# expects to be run fom <dir>/<TASKRCA>/<name>.<ASR>.bam
# e.g. Run/reCalAln/name.ashrr.bam
# change that by setting TASKRCA=TASKBWA and ASR=ASD
################################################################################

if [ -n "$DEPTHOFCOVERAGE" ]
then
    CPUS=24

    echo -e "********* $TASKDOC"

    if [ ! -d $QOUT/$TASKDOC ]; then mkdir -p $QOUT/$TASKDOC; fi

    for dir in ${DIR[@]}
      do

      for f in $( ls $SOURCE/fastq/$dir/*$READONE.$FASTQ )
	do
	
	n=`basename $f`
	n2=${n/%$READONE.$FASTQ/.$ASR.bam}
	name=${n/%$READONE.$FASTQ/}
	echo -e ">>>>>"$dir/$TASKRCA/$n2

	if [ ! -d $OUT/$dir/$TASKDOC ]; then mkdir -p $OUT/$dir/$TASKDOC; fi

	# remove old pbs output
	if [ -e $QOUT/$TASKDOC/$dir'_'$name'.out' ]; then rm $QOUT/$TASKDOC/$dir'_'$name'.out'; fi

	#check if this is part of the pipe and jobsubmission needs to wait
#	if [ -n "$$RUNREALRECAL" ]; then HOLD="-hold_jid "$TASKRCA"_"$dir"_"$name; fi

	#Submit
	if [ -n "$ARMED" ]; then
	    qsub $PRIORITY -j y -o $QOUT/$TASKDOC/$dir'_'$name'.out' -cwd -b y -pe mpich $CPUS \
		-l mem_free=11G -l h_vmem=11G -l vf=500K -N $TASKDOC'_'$dir'_'$name $HOLD\
		${NGSANE_BASE}/mods/gatkDOC.sh -k ${NGSANE_BASE} -f $OUT/$dir/$TASKRCA/$n2 -r $FASTA \
		-o $OUT/$dir/$TASKDOC -t $CPUS
	fi

      done
    done

fi




################################################################################
# downsample
################################################################################

if [ -n "$DOWNSAMPLE" ]
then

    echo -e "********* $TASKDOWN"

    if [ ! -d $QOUT/$TASKDOWN ]; then mkdir -p $QOUT/$TASKDOWN; fi

    for dir in ${DIR[@]}
      do

      for f in $( ls $SOURCE/fastq/$dir/*$READONE.$FASTQ )
	do
	
	n=`basename $f`
	n2=${n/%$READONE.$FASTQ/.$ASD.bam}
	name=${n/%$READONE.$FASTQ/}
	echo -e ">>>>>"$dir$n2

	if [ ! -d $OUT/$dir/$TASKDOWN ]; then mkdir -p $OUT/$dir/$TASKDOWN; fi

	# remove old pbs output
	if [ -e $QOUT/$TASKDOWN/$dir"_"$n2"0.out" ]; then rm $QOUT/$TASKDOWN/$dir"_"$n2.*; fi

	#check if this is part of the pipe and jobsubmission needs to wait
	if [ -n "$mappingBWA" ]; then HOLD="-hold_jid "$TASKBWA"_"$dir"_"$name; fi

	#Submit
	if [ -n "$ARMED" ]; then
	    #downsample
	    python ${NGSANE_BASE}/tools/downsample.py -i $OUT/$dir/bwa/$n2 -o $OUT/$dir/$TASKDOWN/ -t downsample/$dir --region $SEQREG -w 500 -s 500 -q $QOUT/$TASKDOWN/$dir

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

    echo -e "********* $TASKIND"
    if [ ! -d $QOUT/$TASKIND ]; then mkdir -p $QOUT/$TASKIND; fi
    if [ -e $QOUT/$TASKIND/ids.txt ]; then rm $QOUT/$TASKIND/ids.txt; fi
    if [ -e $QOUT/$TASKIND/sum.out ]; then rm $QOUT/$TASKIND/sum.out; fi


    for dir in ${DIR[@]}
      do

      #ensure dirs are there...
      if [ ! -d $OUT/$dir/$TASKIND ]; then mkdir -p $OUT/$dir/$TASKIND; fi

      

      #for f in $( ls $SOURCE/$dir/aln2/*$ASR.bam )
      for f in $( ls $SOURCE/fastq/$dir/*$READONE.fastq )
	do
	n=`basename $f`
	n2=${n/%$READONE.$FASTQ/.$ASR.bam}
	name=${n/%$READONE.$FASTQ/}
	echo -e ">>>>>"$dir$n2

	# remove old pbs output
	if [ -e $QOUT/$TASKIND/$dir"_"$name.out ]; then rm $QOUT/$TASKIND/$dir"_"$name.*; fi

	#check if this is part of the pipe and jobsubmission needs to wait
	if [ -n "$RUNMAPPINGBWA" ]; then HOLD="-hold_jid "$TASKBWA"_"$dir"_"$name; fi
	if [ -n "$recalibrateQualScore" ]; then HOLD="-hold_jid "$TASKRCA"_"$dir"_"$name; fi

	#Submit
	if [ -n "$ARMED" ]; then
	    qsub $PRIORITY -j y -o $QOUT/$TASKIND/$dir'_'$name'.out' -cwd -b y \
		-l h_vmem=12G -N $TASKIND'_'$dir'_'$name $HOLD\
		${NGSANE_BASE}/mods/gatkIndel.sh ${NGSANE_BASE} $OUT/$dir/$TASKRCA/$n2 $FASTA $DBROD \
		$REFSEQROD $OUT/$dir/$TASKIND
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

    echo -e "********* $TASKIND"

    if [ ! -d $QOUT/$TASKIND ]; then mkdir -p $QOUT/$TASKIND; fi
    if [ -e $QOUT/$TASKIND/ids.txt ]; then rm $QOUT/$TASKIND/ids.txt; fi
    if [ ! -d $OUT/genotype ]; then mkdir -p $OUT/genotype; fi

    NAME=$(echo ${DIR[@]}|sed 's/ /_/g')

    for dir in ${DIR[@]};do
      for f in $( ls $SOURCE/fastq/$dir/*$READONE.fastq );do
	n=`basename $f`
	echo $OUT/$dir/$TASKRCA/${n/%$READONE.$FASTQ/.$ASR.bam} >> $TASKIND"bamfiles.tmp"
      done
    done

    # remove old pbs output
    if [ -e $QOUT/$TASKIND/$NAME.out ]; then rm $QOUT/$TASKIND/$NAME.out ; fi

    #check if this is part of the pipe and jobsubmission needs to wait
    if [ -n "$recalibrateQualScore" ]; then 
	HOLD="-hold_jid "
	for dir in ${DIR[@]};do
	    HOLD=$HOLD","$TASKRCA"_"$dir"*"
        done
	HOLD=${HOLD/,/}
    fi

    #Submit
    if [ -n "$ARMED" ]; then
	qsub -j y -o $QOUT/$TASKIND/$NAME.out -cwd -b y \
	    -l mem_free=20G -l h_vmem=20G -N $TASKIND"_"$NAME.out $HOLD\
	    ${NGSANE_BASE}/mods/gatkIndelV2.sh ${NGSANE_BASE} $TASKIND"bamfiles.tmp" $FASTA $DBROD \
	    $REFSEQROD $OUT/genotype $SEQREG
    fi

    
fi

########
# run differental expression detection
########
if [ -n "$RUNANNOTATION" ]; then

    # ensure directory is there
    if [ ! -d $QOUT/$TASKANNOVAR ]; then mkdir -p $QOUT/$TASKANNOVAR; fi
    if [ ! -d $OUT/$TASKANNOVAR ]; then mkdir -p $OUT/$TASKANNOVAR; fi

    
#    for dir in $( ls -d $TASKVAR/* ); do
     for d in ${DIR[@]}; do	

	dir=$OUT/$TASKVAR/$d

	n=`basename $dir`
	echo $n

	namesnp=$OUT/$TASKVAR/$n/$n.filter.snps.vcf
	namesnp2=$OUT/$TASKVAR/$n/$n.recalfilt.snps.vcf
	nameindel=$OUT/$TASKVAR/$n/$n.filter.indel.vcf

        #cleanup old qouts
        if [ -e $QOUT/$TASKANNOVAR/$n.out ]; then rm $QOUT/$TASKANNOVAR/$n.out; fi
	#ensure subfolder is there
	if [ ! -d $OUT/$TASKANNOVAR/$n ]; then mkdir -p $OUT/$TASKANNOVAR/$n; fi


        #check if this is part of the pipe and jobsubmission needs to wait
        if [ -n "$TASKVAR" ]; then HOLD="-hold_jid "$TASKVAR"_"$n; fi
	#echo -e "-hold_jid "$TASKVAR"_"$n
	HOLD="-hold_jid "$TASKVAR"_"$n

        #submit
	if [ -n "$ARMED" ]; then
	    qsub $PRIORITY -b y -cwd -j y -o $QOUT/$TASKANNOVAR/$n.out \
		-N $TASKANNOVAR'_'$n $HOLD\
		${NGSANE_BASE}/mods/annovar.sh -k ${NGSANE_BASE} -i1 $namesnp -i2 $namesnp2 -i3 $nameindel -r $FASTA \
		-o $OUT/$TASKANNOVAR/$n 

	fi
   done

fi


################################################################################
#   call snps with GATK
################################################################################

if [ -n "$GATKcallSNPS" ]
then

    echo -e "********* $TASKSNP"

    if [ ! -d $QOUT/$TASKSNP ]; then mkdir -p $QOUT/$TASKSNP; fi
    if [ -e $QOUT/$TASKSNP/ids.txt ]; then rm $QOUT/$TASKSNP/ids.txt; fi
    if [ -e $QOUT/$TASKSNP/sum.out ]; then rm $QOUT/$TASKSNP/sum.out; fi


    for dir in ${DIR[@]}
      do

      #ensure dirs are there...
      if [ ! -d $OUT/$dir/$TASKSNP ]; then mkdir -p $OUT/$dir/$TASKSNP; fi

      

      #for f in $( ls $SOURCE/$dir/aln2/*$ASR.bam )
      for f in $( ls $SOURCE/fastq/$dir/*$READONE.$FASTQ )
	do
	n=`basename $f`
	n2=${n/%$READONE.$FASTQ/.$ASR.bam}
	name=${n/%$READONE.$FASTQ/}
	echo -e ">>>>>"$dir$n2

	# remove old pbs output
	if [ -e $QOUT/$TASKSNP/$dir"_"$name.out ]; then rm $QOUT/$TASKSNP/$dir"_"$name.*; fi

	#check if this is part of the pipe and jobsubmission needs to wait
	if [ -n "$GATKcallIndelsSeperate" ] || [ -n "$GATKcallIndelsCombined" ]
	    then HOLD="-hold_jid "$TASKIND"_"$dir"_"$name; fi

	#Submit
	if [ -n "$ARMED" ]; then
	    qsub $PRIORITY -j y -o $QOUT/$TASKSNP/$dir'_'$name'.out' -cwd -b y \
		-l h_vmem=12G -N $TASKSNP'_'$dir'_'$name $HOLD\
		${NGSANE_BASE}/mods/gatkSNPs.sh $CONFIG $OUT/$dir/$TASKRCA/$n2 $FASTA $DBVCF \
		$REFSEQROD $OUT/$dir/$TASKSNP $OUT/$dir/$TASKIND | cut -d "" -f 2 >>$QOUT/$TASKSNP/ids.txt
	fi


      done
    done


    #combine into one file
#    qsub -j y -o $QOUT/gatkInd/final.out -cwd -b y -hold_jid `cat $QOUT/gatkInd/ids.txt`\
#	echo -e "merge"
#		merge-vcf `ls $QOUT/dindel/*.dindel.VCF` > genotyping/indels.dindel.vcf
#		merge-vcf A.vcf.gz B.vcf.gz C.vcf.gz | bgzip -c > out.vcf.gz

fi




########
# run differental expression detection
########
if [ -n "$RUNCUFFDIFF" ]; then

    CPUS=24

    # ensure directory is there
    if [ ! -d $QOUT/$TASKCUFFDIFF ]; then mkdir -p $QOUT/$TASKCUFFDIFF; fi
    if [ ! -d $OUT/$TASKCUFFDIFF ]; then mkdir -p $OUT/$TASKCUFFDIFF; fi

    for dir in ${DIR[@]}; do
      for f in $( ls $SOURCE/fastq/$dir/*$READONE.$FASTQ );	do
	n=`basename $f`
	name=$name${n/%$READONE.$FASTQ/}","
      done
    done
    name=$(echo $name | sed 's/\(.*\),/\1/')
    n=${name/,/_u_}

    #cleanup old qouts
    if [ -e $QOUT/$TASKCUFFDIFF/run.out ]; then rm $QOUT/$TASKCUFFDIFF/run.out; fi

    #check if this is part of the pipe and jobsubmission needs to wait
    if [ -n "$RUNTOPHATCUFF" ]; then HOLD="-hold_jid "$TASKTOPHAT"_"$dir"_"$name; fi

    #submit
    if [ -n "$ARMED" ]; then
      qsub $PRIORITY -b y -cwd -j y -o $QOUT/$TASKCUFFDIFF/$n.out \
     	    -N $TASKCUFFDIFF'_'$n -pe mpich $CPUS $HOLD\
            ${NGSANE_BASE}/mods/cuffdiff.sh -k ${NGSANE_BASE} -b $name -r $FASTA \
	    -t $CPUS -o $OUT/$TASKCUFFDIFF/$n -a $REFSEQGTF

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
        $QSUB $ARMED -k $CONFIG -t $TASKINCHWORM -i $INPUT_INCHWORM -e $READONE.$FASTQ -n $NODES_INCHWORM \
    	         -c $NCPU_INCHWORM -m $MEMORY_INCHWORM"G" -w $WALLTIME_INCHWORM -q $NODETYPE_INCHWORM \
    	         --command "${NGSANE_BASE}/mods/trinity_inchworm.sh -k $CONFIG -f <FILE> -o $OUT/<DIR>/$TASKINCHWORM/"\
    ) && echo -e "$JOBIDS" 
	JOBIDS=$(waitForJobIds "$JOBIDS")
    
    #  ##########   Chrysalis  ###########
    JOBIDS=$( \
        $QSUB $ARMED -k $CONFIG -t $TASKCHRYSALIS -i $INPUT_CHRYSALIS -e $READONE.$FASTQ -n $NODES_CHRYSALIS \
         -c $NCPU_CHRYSALIS -m $MEMORY_CHRYSALIS"G" -w $WALLTIME_CHRYSALIS -q $NODETYPE_CHRYSALIS $JOBIDS \
         --command "${NGSANE_BASE}/mods/trinity_chrysalis.sh -k $CONFIG -f <FILE> -o $OUT/<DIR>/$TASKCHRYSALIS/"\
    ) && echo -e "$JOBIDS" 
	JOBIDS=$(waitForJobIds "$JOBIDS")

    #  ##########   Butterfly  ###########
    $QSUB $ARMED -k $CONFIG -t $TASKBUTTERFLY -i $INPUT_BUTTERFLY -e $READONE.$FASTQ -n $NODES_BUTTERFLY \
          -c $NCPU_BUTTERFLY -m $MEMORY_BUTTERFLY"G" -w $WALLTIME_BUTTERFLY -q $NODETYPE_BUTTERFLY $JOBIDS \
          --command "${NGSANE_BASE}/mods/trinity_butterfly.sh -k $CONFIG -f <FILE> -o $OUT/<DIR>/$TASKBUTTERFLY/"
 
fi
################################################################################
# individual calls
#  ##########   Inchworm   ###########
if [ -n "$RUNINCHWORM" ] && [ -z "$RUNTRINITY" ]; then
    if [ -z "$NODES_INCHWORM" ] || [ -z "$NCPU_INCHWORM" ] || [ -z "$MEMORY_INCHWORM" ] || [ -z "$WALLTIME_INCHWORM" ] || [ -z "$NODETYPE_INCHWORM" ]; then  echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi

    $QSUB $ARMED -k $CONFIG -t $TASKINCHWORM -i $INPUT_INCHWORM -e $READONE.$FASTQ -n $NODES_INCHWORM \
        -c $NCPU_INCHWORM -m $MEMORY_INCHWORM"G" -w $WALLTIME_INCHWORM -q $NODETYPE_INCHWORM \
        --command "${NGSANE_BASE}/mods/trinity_inchworm.sh -k $CONFIG -f <FILE> -o $OUT/<DIR>/$TASKINCHWORM"
fi
#  ##########   Chrysalis  ###########
if [ -n "$RUNCHRYSALIS" ] && [ -z "$RUNTRINITY" ]; then
    if [ -z "$NODES_CHRYSALIS" ] || [ -z "$NCPU_CHRYSALIS" ] || [ -z "$MEMORY_CHRYSALIS" ] || [ -z "$WALLTIME_CHRYSALIS" ] || [ -z "$NODETYPE_CHRYSALIS" ] ; then  echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi
    
    $QSUB $ARMED -k $CONFIG -t $TASKCHRYSALIS -i $INPUT_CHRYSALIS -e $READONE.$FASTQ -n $NODES_CHRYSALIS \
        -c $NCPU_CHRYSALIS -m $MEMORY_CHRYSALIS"G" -w $WALLTIME_CHRYSALIS -q $NODETYPE_CHRYSALIS $JOBIDS \
        --command "${NGSANE_BASE}/mods/trinity_chrysalis.sh -k $CONFIG -f <FILE> -o $OUT/<DIR>/$TASKCHRYSALIS"
fi
#  ##########   Butterfly  ###########
if [ -n "$RUNBUTTERFLY" ] && [ -z "$RUNTRINITY" ]; then
    if [ -z "$NODES_BUTTERFLY" ] || [ -z "$NCPU_BUTTERFLY" ] || [ -z "$MEMORY_BUTTERFLY" ] || [ -z "$WALLTIME_BUTTERFLY" ] || [ -z "$NODETYPE_BUTTERFLY" ] ; then echo -e "\e[91m[ERROR]\e[0m Server misconfigured"; exit 1; fi
    
    $QSUB $ARMED -k $CONFIG -t $TASKBUTTERFLY -i $INPUT_BUTTERFLY -e $READONE.$FASTQ -n $NODES_BUTTERFLY \
          -c $NCPU_BUTTERFLY -m $MEMORY_BUTTERFLY"G" -w $WALLTIME_BUTTERFLY -q $NODETYPE_BUTTERFLY $JOBIDS \
          --command "${NGSANE_BASE}/mods/trinity_butterfly.sh -k $CONFIG -f <FILE> -o $OUT/<DIR>/$TASKBUTTERFLY" 
fi
