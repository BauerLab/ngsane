#!/bin/bash -e

# Annotation of bam files
# author: Denis C. Bauer
# date: Jan.2013

# messages to look out for -- relevant for the QC.sh script:
# QCVARIABLES,Resource temporarily unavailable
# RESULTFILENAME <SAMPLE>.merg.anno.bed

echo ">>>>> Annotate BAM file "
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"


function usage {
echo -e "usage: $(basename $0) -k CONFIG -f BAM -o OUTDIR [OPTIONS]

Annotating BAM file with annotations in a folder specified in the config
file

required:
  -k | --toolkit <path>     config file
  -f <file>                 bam file

options:
"
exit
}


if [ ! $# -gt 3 ]; then usage ; fi

#INPUTS
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the NGSANE repository
        -f           )          shift; f=$1 ;; # bam file
        --recover-from )        shift; RECOVERFROM=$1 ;; # attempt to recover from log file
        -h | --help )           usage ;;
        * )                     echo "don't understand "$1
    esac
    shift
done

#PROGRAMS
. $CONFIG
. ${NGSANE_BASE}/conf/header.sh
. $CONFIG


# get basename of f
n=${f##*/}

CHECKPOINT="programs"
for MODULE in $MODULE_BAMANN; do module load $MODULE; done  # save way to load modules that itself load other modules
export PATH=$PATH_BAMANN:$PATH
module list
echo "PATH=$PATH"

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--bedtools    --\n "$(bedtools -version)
[ -z "$(which bedtools)" ] && echo "[WARN] bedtools not detected" && exit 1

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="parameters"

# get basename of f
n=${f##*/}

# check library variables are set
if [[ -z "$BAMANNLIB" ]]; then
    echo "[ERROR] library info not set (BAMANNLIB)"
    exit 1;
else
    echo "[NOTE] using libraries in $BAMANNLIB"
fi

# delete old bam files unless attempting to recover
if [ -z "$RECOVERFROM" ]; then
    [ -e $f.merg.anno.bed ] && rm $f.merg.anno.bed
fi

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="recall files from tape"
	
if [ -n "$DMGET" ]; then
    dmget -a $BAMANNLIB/*
	dmget -a ${f}
fi
    
echo -e "\n********* $CHECKPOINT\n"

################################################################################
CHECKPOINT="bedmerge"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

	bedtools bamtobed -i $f | bedtools merge -i - -n  > $f.merg.bed

	# mark checkpoint
    if [ -f $f.merg.bed ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi
 

################################################################################
CHECKPOINT="run annotation"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

	NAMES=""
	NUMBER=$( ls $BAMANNLIB/*.gtf | wc -l )
	for i in $(ls $BAMANNLIB/*.gtf ); do
		NAME=$(basename $i)
		NAMES=$NAMES" "${NAME/.gtf/}
	done

	echo "[NOTE] annotate with $NAMES"

   	# strandedness does not work for RRNA
   	COMMAND="bedtools annotate -counts -i $f.merg.bed -files $( ls $BAMANNLIB/*.gtf ) -names $NAMES | python ${NGSANE_HOME}/tools/addUnannotated.py >$f.merg.anno.bed"
	echo $COMMAND && eval $COMMAND
   
#   	ALLREADS=$(head -n 1 $f.stats | cut -d" " -f1)


	# mark checkpoint
    if [ -f $f.merg.anno.bed ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
	rm $f.merg.bed

	
fi

################################################################################
CHECKPOINT="summarize"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

   	COMMAND="python ${NGSANE_HOME}/tools/sumRows.py -i $f.merg.anno.bed -l 3 -s 4 -e $(expr $NUMBER + 5) -n $(echo $NAMES | sed 's/ /,/g'),unannotated >$f.anno.stats"
	echo $COMMAND && eval $COMMAND
	
#   echo "total Pgenes rRNA tRNA lincRNA miRNA snoRNA snRNA miscRNA PolyA other HiSeq ucsc_rRNA segDups unannotated unmapped" >$f.merg.anno.stats
#   gawk -v all=$ALLREADS 'BEGIN{a=0; b=0; c=0; d=0; e=0; f=0; g=0; h=0; i=0; j=0; k=0; l=0; m=0; n=0; o=0; x=0}{
#        a=a+$4; 
#        if( $5 !=0){b=b+$4; next}; 
#        if( $6 !=0){c=c+$4; next};
#        if( $7 !=0){d=d+$4; next}; 
#        if( $8 !=0){e=e+$4; next}; 
#        if( $9 !=0){f=f+$4; next};
#        if( $10!=0){g=g+$4; next};
#        if( $11!=0){h=h+$4; next};
#        if( $12!=0){i=i+$4; next};
#        if( $13!=0){j=j+$4; next};
#        if( $14!=0){k=k+$4; next};
#        if( $15!=0){l=l+$4; next};
#        if( $16!=0){m=m+$4; next};
#        if( $17!=0){n=n+$4; next};
#        if( $5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17==0){x=x+$4}}
#        END{print all" "b" "c" "d" "e" "f" "g" "h" "i" "j" "k" "l" "m" "n" "x" "all-a}' $f.merg.anno.bed >> $f.merg.anno.stats

#    head -n 1 $f.merg.anno.bed >> $f.merg.anno.stats
    cat $f.merg.anno.bed | sort -k4gr | head -n 20 >> $f.anno.stats  
    
	# mark checkpoint
    if [ -f $f.anno.stats ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi

################################################################################

################################################################################
echo ">>>>> Annotate BAM file - FINISHED"
echo ">>>>> enddate "`date`

