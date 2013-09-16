#!/bin/bash -e

echo ">>>>> Motif discovery with memechip"
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"

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
# RESULTFILENAME <SAMPLE>.summary.txt

if [ ! $# -gt 3 ]; then usage ; fi

#INPUTS
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the NGSANE repository
        -f | --bed )            shift; f=$1 ;; # bed file containing enriched regions (peaks)
        -o | --outdir )         shift; OUTDIR=$1 ;; # output dir 
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

################################################################################
CHECKPOINT="programs"

for MODULE in $MODULE_MEMECHIP; do module load $MODULE; done  # save way to load modules that itself load other modules

export PATH=$PATH_MEMECHIP:$PATH
module list
echo "PATH=$PATH"
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--bedtools    --\n "$(bedtools --version)
[ -z "$(which bedtools)" ] && echo "[ERROR] no bedtools detected" && exit 1
echo -e "--meme-chip   --\n "$(cat `which meme`.bin | strings | grep -A 2 "MEME - Motif discovery tool" | tail -n 1)
[ -z "$(which meme-chip)" ] && echo "[ERROR] meme-chip not detected" && exit 1

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="parameters"

# get basename of f
n=${f##*/}

if [ -z "$FASTA" ] || [ ! -f $FASTA ]; then
    echo "[ERROR] genome not provided" && exit 1
fi

GENOME_CHROMSIZES=${FASTA%%.*}.chrom.sizes
[ ! -f $GENOME_CHROMSIZES ] && echo "[ERROR] GENOME_CHROMSIZES not found. Excepted at $GENOME_CHROMSIZES" && exit 1

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="recall files from tape"

if [ -n "$DMGET" ]; then
	dmget -a ${f}
	dmls -l ${f}
fi

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="get sequence data"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else

    if [ -n "$SLOPBEDADDPARAM" ]; then
        echo "[NOTE] extend bed regions: $EXTENDREGION"
    
        RUN_COMMAND="bedtools slop -i $f -g $GENOME_CHROMSIZES $SLOPBEDADDPARAM  > $OUTDIR/$n"
        echo $RUN_COMMAND && eval $RUN_COMMAND
        f=$OUTDIR/$n
    fi
    
    bedtools getfasta -name -fi $FASTA -bed $f -fo $OUTDIR/${n/$BED/.fasta}
    echo "Peak regions: `wc -l $f | awk '{print $1}'`" > $OUTDIR/${n/$BED/.summary.txt}

    # mark checkpoint
    if [ -f $OUTDIR/${n/$BED/.fasta} ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi

################################################################################
CHECKPOINT="create background model"    

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

    # create background from bed file unless provided
    if [ -z $MEMEBACKGROUND ]; then
        fasta-get-markov -nostatus $FASTAGETMARKOVADDPARAM < $OUTDIR/${n/$BED/.fasta} > $OUTDIR/${n/$BED/.bg}
        MEMEBACKGROUND=$OUTDIR/${n/$BED/.bg}
    fi
    # mark checkpoint
    if [ -f $OUTDIR/${n/$BED/.bg} ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi

################################################################################
CHECKPOINT="meme-chip"    

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
    
    RUN_COMMAND="meme-chip $MEMECHIPADDPARAM -oc $OUTDIR/${n/$BED/} -bfile $MEMEBACKGROUND -desc ${n/$BED/} -db $MEMECHIPDATABASES -meme-p $CPU_MEMECHIP $OUTDIR/${n/$BED/.fasta}"
    echo $RUN_COMMAND && eval $RUN_COMMAND
    
    # mark checkpoint
    if [ -f $OUTDIR/${n/$BED/}/combined.meme ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi

################################################################################
CHECKPOINT="classify bound regions"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
   
    RUN_COMMAND="fimo $FIMOADDPARAM --bgfile $MEMEBACKGROUND --oc $OUTDIR/${n/$BED/_fimo} $OUTDIR/${n/$BED/}/combined.meme $OUTDIR/${n/$BED/.fasta}"
    echo $RUN_COMMAND && eval $RUN_COMMAND

    sort -k4,4 -k1,1 -k2,2g $f > $OUTDIR/${n/$BED/_sorted.bed}
    for PATTERN in $(tail -n+2 $OUTDIR/${n/$BED/_fimo}/fimo.txt | awk '{print $1}' | sort -u); do
        echo "[NOTE] Motif: $PATTERN"
    
        grep -P "^${PATTERN}\t" $OUTDIR/${n/$BED/_fimo}/fimo.txt | cut -f2-4,6 | tail -n+2 | sort -k1,1 > $OUTDIR/${n/$BED/_fimo}/$PATTERN.txt
    
        join -1 1 -2 4 $OUTDIR/${n/$BED/_fimo}/$PATTERN.txt $OUTDIR/${n/$BED/_sorted.bed} | awk '{OFS="\t"; print $5,$6+$2,$6+$3,$1,$4,$9}' > $OUTDIR/${n/$BED/_motif}_${PATTERN}.direct.bed
        
        comm -13 <(awk '{print $1}' $OUTDIR/${n/$BED/_fimo}/$PATTERN.txt | sort -u ) <(awk '{print $4}' $f | sort -u ) > $OUTDIR/${n/$BED/_fimo}/${PATTERN}_tmp.txt
    
        join -1 4 -2 1 $OUTDIR/${n/$BED/_sorted.bed} $OUTDIR/${n/$BED/_fimo}/${PATTERN}_tmp.txt | awk '{OFS="\t"; print 2,$3,$4,$1,$5,$6}' > $OUTDIR/${n/$BED/_motif}_${PATTERN}.indirect.bed
        
        echo "Motif $PATTERN bound directely (strong site): $(cat $OUTDIR/${n/$BED/_motif}_${PATTERN}.direct.bed | awk '{print $4}' | sort -u | wc -l | awk '{print $1}')" >> $OUTDIR/${n/$BED/.summary.txt}
        echo "Motif $PATTERN bound indirectely (weak or no site): $(cat $OUTDIR/${n/$BED/_motif}_${PATTERN}.indirect.bed | awk '{print $4}' | sort -u | wc -l | awk '{print $1}')" >> $OUTDIR/${n/$BED/.summary.txt}
    done
    
    # mark checkpoint
    if [ -f $OUTDIR/${n/$BED/_motif}_${PATTERN}.direct.bed ] && [ -f $OUTDIR/${n/$BED/_motif}_${PATTERN}.indirect.bed ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi

################################################################################
CHECKPOINT="cleanup"    

[ -e $OUTDIR/${n/$BED/.fasta} ] && rm $OUTDIR/${n/$BED/.fasta}
[ -d $OUTDIR/${n/$BED/_fimo} ] && rm -r $OUTDIR/${n/$BED/_fimo}
[ -e $OUTDIR/${n/$BED/_sorted.bed} ] && rm $OUTDIR/${n/$BED/_sorted.bed}
[ -e $OUTDIR/${n/$BED/.bg} ] && rm $OUTDIR/${n/$BED/.bg} 

echo -e "\n********* $CHECKPOINT\n"
################################################################################
[ -e $OUTDIR/${n/$BED/.summary.txt}.dummy ] && rm $OUTDIR/${n/$BED/.summary.txt}.dummy
echo ">>>>> Motif discovery with memechip - FINISHED"
echo ">>>>> enddate "`date`

