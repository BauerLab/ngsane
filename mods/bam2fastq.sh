#!/bin/bash -e

# Script to convert bam files back to fastq
# It takes a <Run>/<TASK>/*.bam file and returns fastq/<Run>/*.$FASTQ[.gz]
#
# author: Denis Bauer
# date: Nov 2013

# messages to look out for -- relevant for the QC.sh script:
# QCVARIABLES,
# RESULTFILENAME fastq/<DIR>/<SAMPLE>$READONE.$FASTQ

echo ">>>>> bam2fastq"
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"

while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of NGSANE
        -f | --file )           shift; f=$1 ;; # fastq file
        -o | --outdir )         shift; OUTDIR=$1 ;; # output dir                                                     
        --recover-from )        shift; RECOVERFROM=$1 ;; # attempt to recover from log file
        -h | --help )           usage ;;
        * )                     echo "don't understand "$1
    esac
    shift
done

# overwrite defaults
. $CONFIG
. ${NGSANE_BASE}/conf/header.sh
. $CONFIG

################################################################################
CHECKPOINT="programs"

for MODULE in $MODULE_BAM2FASTQ; do module load $MODULE; done  # save way to load modules that itself load other modules
export PATH=$PATH_BAM2FASTQ:$PATH;
module list
echo "PATH=$PATH"
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)
PATH_PICARD=$(dirname $(which MarkDuplicates.jar))
#PATH_RESOLVEPAIR=$(dirname $(which resolvepair))

echo "[NOTE] set java parameters"
JAVAPARAMS="-Xmx"$(python -c "print int($MEMORY_BAM2FASTQ*0.8)")"g -Djava.io.tmpdir="$TMP" -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1" 
unset _JAVA_OPTIONS
echo "JAVAPARAMS "$JAVAPARAMS


echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--JAVA        --\n" $(java -Xmx200m -version 2>&1)
[ -z "$(which java)" ] && echo "[ERROR] no java detected" && exit 1
echo -e "--PICARD      --\n "$(java $JAVAPARAMS -jar $PATH_PICARD/MarkDuplicates.jar --version 2>&1)
[ ! -f $PATH_PICARD/MarkDuplicates.jar ] && echo "[ERROR] no picard detected" && exit 1

echo -e "\n********* $CHECKPOINT"
################################################################################
CHECKPOINT="parameters"

n=${f##*/}

#DIR=$(basename $(basename $(dirname $f)))
#mkdir -p $OUTDIR/fastq/$DIR
f2=$OUTDIR/${n/.$ASD.bam/$READONE.$FASTQ}
GZEND="."${FASTQ/*./}

echo -e "\n********* $CHECKPOINT"
################################################################################
CHECKPOINT="recall files from tape"

if [ -n "$DMGET" ]; then
    dmget -a $f
fi

echo -e "\n********* $CHECKPOINT"
################################################################################
CHECKPOINT="findout pairedness"    

if [[ -n "$RECOVERFROM" ]] && [[ $(grep "********* $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

	if [[ ! -e $f.stats ||  -z $(cat $f.stats) ]]; then	echo "[NOTE] generate stats file"; samtools flagstat $f > $f.stats; fi
	
	#is paired ?
	p=$(grep "paired in " $f.stats | cut -d " " -f 1)
	if [ ! "$p" -eq "0" ]; then
	    PAIRED="1"
	    echo "[NOTE] PAIRED"
	else
	    echo "[NOTE] SINGLE"
	fi

    # mark checkpoint
    if [ -f $f.stats ]; then echo -e "\n********* $CHECKPOINT" && unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
fi


################################################################################
CHECKPOINT="convert bam to fastq"    

if [[ -n "$RECOVERFROM" ]] && [[ $(grep "********* $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

	if [ -n "$PAIRED" ]; then 
#		SECONDREAD="SECOND_END_FASTQ=${f2/$READONE/$READTWO}";
		SECONDREAD="${f2/$READONE/$READTWO}"; 
	fi

	#TODO : resolvpair just gives unique IDs to reads, however, they were multimappers and all but one instance of them should be removed. 
	#python ${PATH_RESOLVEPAIR}/resolvepair
	
	# note $f needs to be there eventhough it is ignored because of -o
	#RUN_COMMAND="samtools view -h $f | head -n 1000 | samtools view -S -b - | samtools sort -o -n - $f | java $JAVAPARAMS -jar ${PATH_PICARD}/SamToFastq.jar $BAM2FASTQADDPARAM VALIDATION_STRINGENCY=SILENT INPUT=/dev/stdin INCLUDE_NON_PF_READS=true FASTQ=$f2 $SECONDREAD"
	#RUN_COMMAND="samtools sort -o -n $f $f | java $JAVAPARAMS -jar ${PATH_PICARD}/SamToFastq.jar $BAM2FASTQADDPARAM VALIDATION_STRINGENCY=SILENT INPUT=/dev/stdin INCLUDE_NON_PF_READS=true FASTQ=$f2 $SECONDREAD"
	#RUN_COMMAND="
	THISTMP=$TMP/$n$RANDOM
	mkdir -p $THISTMP
	htscmd bamshuf -uOn 128 $f $THISTMP | htscmd bam2fq -aO - | gawk -v F0=${f2/$GZEND/} -v F1=${SECONDREAD/$GZEND/} 'NR%4==1{if ($0 ~ "/1$"){x=F0}else{x=F1}}{print > x}'
	#| awk 'NR%4==1{x="F"i++%2}{print > x}' 
    #echo $RUN_COMMAND && eval $RUN_COMMAND
	$GZIP -c ${f2/$GZEND/} > $f2
	[ -e ${SECONDREAD/$GZEND/} ] && $GZIP -c ${SECONDREAD/$GZEND/} > ${SECONDREAD}

    # mark checkpoint
    if [ -f $f2 ]; then echo -e "\n********* $CHECKPOINT" && unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
	rm -f ${f2/$GZEND/} ${SECONDREAD/$GZEND/}
fi


################################################################################
CHECKPOINT="check"    

if [[ -n "$RECOVERFROM" ]] && [[ $(grep "********* $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

	# no singletons
	#BAMREADS=$(grep "with itself and mate mapped" $f.stats | cut -f 1 -d " " )
	# all reads
	BAMREADS=$( cat $f.stats | head -n 1 | cut -f 1 -d " " ) 
	FASTQREADS=$( echo $(zcat ${f2/$READONE.$FASTQ/}*$FASTQ | wc -l | cut -f 1 -d " ")" / 4" | bc)

	if [ $BAMREADS -eq $FASTQREADS ]; then
	    echo "[NOTE] PASS check mapping: $BAMREADS == $FASTQREADS"
	else
	    echo -e "[ERROR] We are loosing reads from .bam -> .fastq in $f: \nFastq had $FASTQREADS Bam has $BAMREADS"
	    exit 1 
	fi
	echo -e "\n********* $CHECKPOINT\n"
fi




################################################################################
[ -e $f2.dummy ] && rm $f2.dummy
echo ">>>>> bam2fastq - FINISHED"
echo ">>>>> enddate "`date`

