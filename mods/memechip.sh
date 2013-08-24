#!/bin/bash -e

echo ">>>>> Motif discovery with memechip"
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*""

function usage {
echo -e "usage: $(basename $0) -k NGSANE -f FASTQ -r REFERENCE -o OUTDIR [OPTIONS]"
exit
}

# Script for de-novo motif discovery using meme-chip
# It takes bed regions that are enriched for the ChIPed molecule.
# It produces enriched DNA binding motifs and run the most enriched motif on the input bed file
# author: Fabian Buske
# date: August 2013

# QCVARIABLES,Resource temporarily unavailable
if [ ! $# -gt 3 ]; then usage ; fi

#INPUTS
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the NGSANE repository
        -f | --bed )            shift; f=$1 ;; # bed file containing enriched regions (peaks)
        -o | --outdir )         shift; MYOUT=$1 ;; # output dir 
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

###################################################################################################
CHECKPOINT="programs"

for MODULE in $MODULE_MEMECHIP; do module load $MODULE; done  # save way to load modules that itself load other modules

export PATH=$PATH_MEMECHIP:$PATH
module list
echo "PATH=$PATH"
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)
echo -e "--bedtools --\n "$(bedtools --version)
[ -z "$(which bedtools)" ] && echo "[ERROR] no bedtools detected" && exit 1
echo -e "--meme-chip         --\n "$(cat `which meme`.bin | strings | grep -A 2 "MEME - Motif discovery tool" | tail -n 1)
[ -z "$(which meme-chip)" ] && echo "[ERROR] meme-chip not detected" && exit 1

echo -n "********* $CHECKPOINT"
###################################################################################################
CHECKPOINT="parameters"

# get basename of f
n=${f##*/}

if [ -z "$FASTA" ] || [ ! -f $FASTA ]; then
    echo "[ERROR] genome not provided" && exit 1
fi
if [ -z "$CHROMSIZES" ] || [ ! -f $CHROMSIZES ]; then
    echo "[ERROR] chromosome sizes not provided" && exit 1
fi

echo -n "********* $CHECKPOINT"
###################################################################################################
CHECKPOINT="recall files from tape"

if [ -n "$DMGET" ]; then
	dmget -a ${f}
	dmls -l ${f}
fi

echo -n "********* $CHECKPOINT"
###################################################################################################
CHECKPOINT="get sequence data"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep "********* $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo -n "::::::::: passed $CHECKPOINT"
else

    if [ -n "$SLOPBEDADDPARAM" ]; then
        echo "[NOTE] extend bed regions: $EXTENDREGION"
    
        RUN_COMMAND="bedtools slop -i $f -g $CHROMSIZES $SLOPBEDADDPARAM  > $MYOUT/$n"
        echo $RUN_COMMAND && eval $RUN_COMMAND
        f=$MYOUT/$n
    fi
    
    bedtools getfasta -name -fi $FASTA -bed $f -fo $MYOUT/${n/$BED/.fasta}
    echo "Peak regions: `wc -l $f | awk '{print $1}'`" > $MYOUT/${n/$BED/.summary.txt}

    # mark checkpoint
    [ -f $MYOUT/${n/$BED/.fasta} ] && echo -n "********* $CHECKPOINT"
fi

###################################################################################################
CHECKPOINT="create background model"    

if [[ -n "$RECOVERFROM" ]] && [[ $(grep "********* $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo -n "::::::::: passed $CHECKPOINT"
else 

    # create background from bed file unless provided
    if [ -z $MEMEBACKGROUND ]; then
        fasta-get-markov -nostatus $FASTAGETMARKOVADDPARAM < $MYOUT/${n/$BED/.fasta} > $MYOUT/${n/$BED/.bg}
        MEMEBACKGROUND=$MYOUT/${n/$BED/.bg}
    fi
    # mark checkpoint
    [ -f $MYOUT/${n/$BED/.bg} ] && echo -n "********* $CHECKPOINT"
fi

###################################################################################################
CHECKPOINT="run meme-chip"    

if [[ -n "$RECOVERFROM" ]] && [[ $(grep "********* $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo -n "::::::::: passed $CHECKPOINT"
else 
    
    RUN_COMMAND="meme-chip $MEMECHIPADDPARAM -oc $MYOUT/${n/$BED/} -bfile $MEMEBACKGROUND -desc ${n/$BED/} -db $MEMECHIPDATABASES -meme-p $CPU_MEMECHIP $MYOUT/${n/$BED/.fasta}"
    echo $RUN_COMMAND && eval $RUN_COMMAND
    
    # mark checkpoint
    [ -f $MYOUT/${n/$BED/}/combined.meme} ] && echo -n "********* $CHECKPOINT"
fi

###################################################################################################
CHECKPOINT="classify bound regions"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep "********* $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo -n "::::::::: passed $CHECKPOINT"
else 
   
    RUN_COMMAND="fimo $FIMOADDPARAM --bgfile $MEMEBACKGROUND --oc $MYOUT/${n/$BED/_fimo} $MYOUT/${n/$BED/}/combined.meme $MYOUT/${n/$BED/.fasta}"
    echo $RUN_COMMAND && eval $RUN_COMMAND

    sort -k4,4 -k1,1 -k2,2g $f > $MYOUT/${n/$BED/_sorted.bed}
    for PATTERN in $(tail -n+2 $MYOUT/${n/$BED/_fimo}/fimo.txt | awk '{print $1}' | sort -u); do
        echo "[NOTE] Motif: $PATTERN"
    
        grep -P "^${PATTERN}\t" $MYOUT/${n/$BED/_fimo}/fimo.txt | cut -f2-4,6 | tail -n+2 | sort -k1,1 > $MYOUT/${n/$BED/_fimo}/$PATTERN.txt
    
        join -1 1 -2 4 $MYOUT/${n/$BED/_fimo}/$PATTERN.txt $MYOUT/${n/$BED/_sorted.bed} | awk '{OFS="\t"; print $5,$6+$2,$6+$3,$1,$4,$9}' > $MYOUT/${n/$BED/_motif}_${PATTERN}.direct.bed
        
        comm -13 <(awk '{print $1}' $MYOUT/${n/$BED/_fimo}/$PATTERN.txt | sort -u ) <(awk '{print $4}' $f | sort -u ) > $MYOUT/${n/$BED/_fimo}/${PATTERN}_tmp.txt
    
        join -1 4 -2 1 $MYOUT/${n/$BED/_sorted.bed} $MYOUT/${n/$BED/_fimo}/${PATTERN}_tmp.txt | awk '{OFS="\t"; print 2,$3,$4,$1,$5,$6}' > $MYOUT/${n/$BED/_motif}_${PATTERN}.indirect.bed
        
        echo "Motif $PATTERN bound directely (strong site): $(cat $MYOUT/${n/$BED/_motif}_${PATTERN}.direct.bed | awk '{print $4}' | sort -u | wc -l | awk '{print $1}')" >> $MYOUT/${n/$BED/.summary.txt}
        echo "Motif $PATTERN bound indirectely (weak or no site): $(cat $MYOUT/${n/$BED/_motif}_${PATTERN}.indirect.bed | awk '{print $4}' | sort -u | wc -l | awk '{print $1}')" >> $MYOUT/${n/$BED/.summary.txt}
    done
    
    # mark checkpoint
    [ -f $MYOUT/${n/$BED/_motif}_${PATTERN}.direct.bed ] && [ -f $MYOUT/${n/$BED/_motif}_${PATTERN}.indirect.bed ] && echo -n "********* $CHECKPOINT"
fi

###################################################################################################
CHECKPOINT="cleanup"    

rm -rf $MYOUT/${n/$BED/.fasta} $MYOUT/${n/$BED/_fimo} $MYOUT/${n/$BED/_sorted.bed} $MYOUT/${n/$BED/.bg} $MYOUT/$n

echo -n "********* $CHECKPOINT"
###################################################################################################
echo ">>>>> Motif discovery with memechip - FINISHED"
echo ">>>>> enddate "`date`

