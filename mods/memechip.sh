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
# RESULTFILENAME <DIR>/<TASK>/<SAMPLE>.summary.txt

if [ ! $# -gt 3 ]; then usage ; fi

#INPUTS
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the NGSANE repository
        -f | --bed )            shift; f=$1 ;; # bed file containing enriched regions (peaks)
        -o | --outdir )         shift; OUTDIR=$1 ;; # output dir 
        --recover-from )        shift; NGSANE_RECOVERFROM=$1 ;; # attempt to recover from log file
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
NGSANE_CHECKPOINT_INIT "programs"

# save way to load modules that itself loads other modules
hash module 2>/dev/null && for MODULE in $MODULE_MEMECHIP; do module load $MODULE; done && module list 

export PATH=$PATH_MEMECHIP:$PATH
echo "PATH=$PATH"
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--bedtools    --\n "$(bedtools --version)
[ -z "$(which bedtools)" ] && echo "[ERROR] no bedtools detected" && exit 1
echo -e "--perl        --\n "$(perl -v | grep "This is perl" )
[ -z "$(which perl)" ] && echo "[ERROR] no perl detected" && exit 1
echo -e "--meme-chip   --\n "$(cat `which meme` | strings | grep -A 2 "MEME - Motif discovery tool" | tail -n 1)
[ -z "$(which meme-chip)" ] && echo "[ERROR] meme-chip not detected" && exit 1

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "parameters"

# get basename of f
n=${f##*/}
SAMPLE=${n/$BED/}

if [ -z "$FASTA" ] || [ ! -f $FASTA ]; then
    echo "[ERROR] no reference provided (FASTA)"
    exit 1
else
    echo "[NOTE] Reference: $FASTA"
fi

GENOME_CHROMSIZES=${FASTA%.*}.chrom.sizes
if [ ! -f $GENOME_CHROMSIZES ]; then
    echo "[ERROR] GENOME_CHROMSIZES not found. Excepted at $GENOME_CHROMSIZES"
    exit 1
else
    echo "[NOTE] Chromosome size: $GENOME_CHROMSIZES"
fi

if [ -n "$MEMEBACKGROUNDFILE" ]; then
    MEMEBACKGROUND=$MEMEBACKGROUNDFILE
else
    MEMEBACKGROUND=$OUTDIR/$SAMPLE.bg
fi

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "recall files from tape"

if [ -n "$DMGET" ]; then
	dmget -a ${f}
	dmget -a $OUTDIR/*
fi

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "get sequence data"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then


    if [ -n "$SLOPBEDADDPARAM" ]; then
        echo "[NOTE] extend bed regions: $EXTENDREGION"
    
        RUN_COMMAND="bedtools slop -i $f -g $GENOME_CHROMSIZES $SLOPBEDADDPARAM  > $OUTDIR/$n"
        echo $RUN_COMMAND && eval $RUN_COMMAND
        f=$OUTDIR/$n
    fi
    
    bedtools getfasta -name -fi $FASTA -bed $f -fo $OUTDIR/$SAMPLE.fasta

    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE.fasta

fi

################################################################################
NGSANE_CHECKPOINT_INIT "create background model"    

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    # create background from bed file unless provided
    if [ -z "$MEMEBACKGROUNDFILE" ]; then
        echo "[NOTE] create background from peak regions"
        fasta-get-markov -nostatus $FASTAGETMARKOVADDPARAM < $OUTDIR/$SAMPLE.fasta > $MEMEBACKGROUND
    fi
    
    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK "$MEMEBACKGROUND"

fi

################################################################################
NGSANE_CHECKPOINT_INIT "meme-chip"    

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then
    
    RUN_COMMAND="meme-chip $MEMECHIPADDPARAM -oc $OUTDIR/$SAMPLE -bfile $MEMEBACKGROUND -desc $SAMPLE -db $MEMECHIPDATABASES -meme-p $CPU_MEMECHIP $OUTDIR/$SAMPLE.fasta"
    echo $RUN_COMMAND && eval $RUN_COMMAND
    
    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE/combined.meme

fi

################################################################################
NGSANE_CHECKPOINT_INIT "classify bound regions"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    echo "Peak regions: $(wc -l $f |  cut -d' ' -f 1)" > $OUTDIR/$SAMPLE.summary.txt
    
    # extract top motif from combined motif
    sed -n '/MEME/,/MOTIF 1/p'  $OUTDIR/$SAMPLE/combined.meme | head -n-1 >  $OUTDIR/$SAMPLE/topmotif.meme
    # get motif from meme that best matches top hit
    RUN_COMMAND="tomtom -verbosity 1 -text -thresh 0.01 $OUTDIR/$SAMPLE/topmotif.meme $OUTDIR/$SAMPLE/meme_out/meme.xml | tail -n+2 1> $OUTDIR/$SAMPLE/topmotif_memehit.txt 2> /dev/null"
    echo $RUN_COMMAND && eval $RUN_COMMAND
    
    if [ -s $OUTDIR/$SAMPLE/topmotif_memehit.txt ]; then
        # get target id motif from hit
        MEMEMOTIFNUM=$(cut -f 2 $OUTDIR/$SAMPLE/topmotif_memehit.txt | head -n 1)
        MEMECONSENSUS=$(cut -f 9 $OUTDIR/$SAMPLE/topmotif_memehit.txt | head -n 1)
        MOTIFNUM=$(awk -v CONS=$MEMECONSENSUS  '{if ($9 == CONS){print $0; exit 1}}' $OUTDIR/$SAMPLE/motif_alignment.txt | cut -f 2)
        if [ -z "$MOTIFNUM" ]; then 
            # try reverse complement
            MEMECONSENSUS=$(echo "$MEMECONSENSUS" | tr [GCTAgcta] [CGATcgat] | rev)
            MOTIFNUM=$(awk -v CONS="$MEMECONSENSUS"  '{if ($9 == CONS){print $0; exit 1}}' $OUTDIR/$SAMPLE/motif_alignment.txt | cut -f 2)
        fi
    else
        # otherwise take longest motif that clusters with the top motif (avoid running fimo on dreme results
        MOTIFNUM=$(awk '{if ($1==0 && $6<0.01){OFS="\t";print $0,length($9)}}' $OUTDIR/$SAMPLE/motif_alignment.txt | sort -k11,11gr | cut -f 2 | head -n 1)
        MEMECONSENSUS=$(sed -n '2,2p' $OUTDIR/$SAMPLE/motif_alignment.txt | cut -f 8)
    fi
    echo "MEME motif: $MEMEMOTIFNUM" >> $OUTDIR/$SAMPLE.summary.txt
    cp $OUTDIR/$SAMPLE/meme_out/logo$MEMEMOTIFNUM.png $OUTDIR/${SAMPLE}_topmememotif.png
    echo "Combined motif: $MOTIFNUM" >> $OUTDIR/$SAMPLE.summary.txt
    echo "Query consensus: $MEMECONSENSUS" >> $OUTDIR/$SAMPLE.summary.txt
       
    RUN_COMMAND="fimo $FIMOADDPARAM --motif $MOTIFNUM --bgfile $MEMEBACKGROUND --oc $OUTDIR/$SAMPLE"_"fimo $OUTDIR/$SAMPLE/combined.meme $OUTDIR/$SAMPLE.fasta"
    echo $RUN_COMMAND && eval $RUN_COMMAND
    
    if [[ "$(wc -l $OUTDIR/$SAMPLE"_"fimo/fimo.txt | cut -d' ' -f 1)" -le 1 ]]; then
        echo "[NOTE] no motif occurences enriched using fimo with given cutoff"
    
    else

        for PATTERN in $(tail -n+2 $OUTDIR/$SAMPLE"_"fimo/fimo.txt | awk '{print $1}' | sort -u); do

            sort -k4,4 -k1,1 -k2,2g $f > $OUTDIR/$SAMPLE"_"sorted.bed
            grep -P "^${PATTERN}\t" $OUTDIR/$SAMPLE"_"fimo/fimo.txt | cut -f2-4,6 | tail -n+2 | sort -k1,1 > $OUTDIR/$SAMPLE"_"fimo/$PATTERN.txt
        
            join -1 1 -2 4 $OUTDIR/$SAMPLE"_"fimo/$PATTERN.txt $OUTDIR/$SAMPLE"_"sorted.bed | awk '{OFS="\t"; print $5,$6+$2,$6+$3,$1,$4,$9}' > $OUTDIR/${n/$BED/_motif}_${PATTERN}.direct.bed
            
            comm -13 <(awk '{print $1}' $OUTDIR/$SAMPLE"_"fimo/$PATTERN.txt | sort -u ) <(awk '{print $4}' $f | sort -u ) > $OUTDIR/$SAMPLE"_"fimo/${PATTERN}_tmp.txt
        
            join -1 4 -2 1 $OUTDIR/$SAMPLE"_"sorted.bed $OUTDIR/$SAMPLE"_"fimo/${PATTERN}_tmp.txt | awk '{OFS="\t"; print 2,$3,$4,$1,$5,$6}' > $OUTDIR/${n/$BED/_motif}_${PATTERN}.indirect.bed
            
            echo "Motif $PATTERN bound directly (strong site): $(cat $OUTDIR/${n/$BED/_motif}_${PATTERN}.direct.bed | awk '{print $4}' | sort -u | wc -l | awk '{print $1}')" >> $OUTDIR/$SAMPLE.summary.txt
            echo "Motif $PATTERN bound indirectly (weak or no site): $(cat $OUTDIR/${n/$BED/_motif}_${PATTERN}.indirect.bed | awk '{print $4}' | sort -u | wc -l | awk '{print $1}')" >> $OUTDIR/$SAMPLE.summary.txt
            EVALUE=$(grep "letter-probability matrix" $OUTDIR/$SAMPLE/combined.meme | awk -v MOTIF=`expr $MOTIFNUM + 1` '{;if (NR==MOTIF){IFS="=";print $10}}')
            echo "E-value: $EVALUE">> $OUTDIR/$SAMPLE.summary.txt
            TOMTOM=$(awk -v MEME=$MEMEMOTIFNUM '{if ($1 == MEME){print $0; exit}}' $OUTDIR/$SAMPLE/meme_tomtom_out/tomtom.txt)
            echo "Most similar known motif: "$(echo "$TOMTOM" | awk '{print $2}') >> $OUTDIR/$SAMPLE.summary.txt
            echo "Q-value: "$(echo "$TOMTOM" | awk '{print $6}') >> $OUTDIR/$SAMPLE.summary.txt
            echo "Target consensus: "$(echo "$TOMTOM" | awk '{print $9}') >> $OUTDIR/$SAMPLE.summary.txt
        done
        
        [ -e $OUTDIR/$SAMPLE"_"sorted.bed ] && rm $OUTDIR/$SAMPLE"_"sorted.bed
    fi
    
    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE.summary.txt
    
fi

################################################################################
NGSANE_CHECKPOINT_INIT "cleanup"    

[ -e $OUTDIR/$SAMPLE.fasta ] && rm $OUTDIR/$SAMPLE.fasta
[ -e $OUTDIR/$SAMPLE.bg ] && rm $OUTDIR/$SAMPLE.bg 
[ -d $OUTDIR/$SAMPLE"_"fimo ] && rm -r $OUTDIR/$SAMPLE"_"fimo
[ -f $OUTDIR/$SAMPLE/topmotif.meme ] && rm $OUTDIR/$SAMPLE/topmotif.meme $OUTDIR/$SAMPLE/topmotif_memehit.txt

NGSANE_CHECKPOINT_CHECK
################################################################################
[ -e $OUTDIR/$SAMPLE.summary.txt.dummy ] && rm $OUTDIR/$SAMPLE.summary.txt.dummy
echo ">>>>> Motif discovery with memechip - FINISHED"
echo ">>>>> enddate "`date`

