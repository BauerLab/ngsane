#!/bin/bash -e

# Script running HIC HOMER pipeline
# It expects bam files, paired end, as input.
# author: Fabian Buske
# date: August 2013

echo ">>>>> HiC analysis with homer"
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"

function usage {
echo -e "usage: $(basename $0) -k NGSANE -f bam -o OUTDIR [OPTIONS]"
exit
}
# QCVARIABLES,Resource temporarily unavailable
# RESULTFILENAME <DIR>/<TASK>/<SAMPLE>_significantInteractions.log

if [ ! $# -gt 3 ]; then usage ; fi

#INPUTS                                                                                                           
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the NGSANE repository
        -f | --bam )            shift; f=$1 ;; # bam file
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
hash module 2>/dev/null && for MODULE in $MODULE_HOMERHIC; do module load $MODULE; done && module list 

export PATH=$PATH_HOMERHIC:$PATH
echo "PATH=$PATH"
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--samtools    --\n "$(samtools 2>&1 | head -n 3 | tail -n-2)
[ -z "$(which samtools)" ] && echo "[ERROR] no samtools detected" && exit 1
echo -e "--R           --\n "$(R --version | head -n 3)
[ -z "$(which R)" ] && echo "[ERROR] no R detected" && exit 1
echo -e "--homer       --\n "$(which makeTagDirectory)
[ -z "$(which makeTagDirectory)" ] && echo "[ERROR] homer not detected" && exit 1
echo -e "--circos      --\n "$(circos --version)
[ -z "$(which circos)" ] && echo "[WARN] circos not detected"
echo -e "--wigToBigWig --\n "$(wigToBigWig 2>&1 | tee | head -n 1)
[ -z "$(which wigToBigWig)" ] && echo "[WARN] wigToBigWig not detected"

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "parameters"

# get basename of f
if [ -z "$POOLED_DATA_NAME" ]; then 
    n=${f##*/}
    SAMPLE=${n/%$READONE$ASD.bam/}
    
    #is paired ?                                                                                                      
    if [ "$f" != "${f/%$READONE$ASD.bam/$READTWO$ASD.bam}" ] && [ -e ${f/%$READONE$ASD.bam/$READTWO$ASD.bam} ]; then
        PAIRED="1"
    else
        PAIRED="0"
        echo "[ERROR] paired library required for HIC analysis" && exit 1
    fi

else
    SAMPLE=$POOLED_DATA_NAME
    array=(${f//,/ })
    #is paired ?                                                                                                      
    if [ "${array[i]}" != "${array[i]/%$READONE$ASD.bam/$READTWO$ASD.bam}" ] && [ -e ${array[i]/%$READONE$ASD.bam/$READTWO$ASD.bam} ]; then
        PAIRED="1"
    else
        PAIRED="0"
        echo "[ERROR] paired library required for HIC analysis" && exit 1
    fi
    
fi

if [ -z "$FASTA" ]; then
    echo "[ERROR] no reference provided (FASTA)"
    exit 1
fi

GENOME_CHROMSIZES=${FASTA%.*}.chrom.sizes
if [ ! -f $GENOME_CHROMSIZES ]; then
    echo "[WARN] GENOME_CHROMSIZES not found. Excepted at $GENOME_CHROMSIZES. Will not create bigBed file"
else
    echo "[NOTE] Chromosome size: $GENOME_CHROMSIZES"
fi


THISTMP=$TMP"/"$(whoami)"/"$(echo $OUTDIR/$SAMPLE | md5sum | cut -d' ' -f1)
[ -d $THISTMP ] && rm -r $THISTMP
mkdir -p $THISTMP

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "recall files from tape"

if [ -n "$DMGET" ]; then
	dmget -a ${f/$READONE/"*"}
	dmget -a $OUTDIR/*
fi

NGSANE_CHECKPOINT_CHECK
################################################################################

#homer likes to write in the current directory, so change to target
CURDIR=$(pwd)
cd $OUTDIR

################################################################################
NGSANE_CHECKPOINT_INIT "create unfiltered tagdirectory"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then


    if [ -z "$POOLED_DATA_NAME" ]; then 
        cp $f ${f/%$READONE$ASD.bam/$READTWO$ASD.bam} $THISTMP
        RUN_COMMAND="makeTagDirectory $THISTMP/$SAMPLE"_tagdir_unfiltered" $THISTMP/$n,$THISTMP/${n/%$READONE$ASD.bam/$READTWO$ASD.bam} -format sam -illuminaPE -tbp 1"
    else
        # pool data
        RUN_COMMAND="makeTagDirectory $THISTMP/$SAMPLE"_tagdir_unfiltered
        array=(${f//,/ })
        for i in "${!array[@]}"; do
            cp ${array[i]} ${array[i]/%$READONE$ASD.bam/$READTWO$ASD.bam} $THISTMP
            n=${array[i]##*/} 
            RUN_COMMAND="$RUN_COMMAND $THISTMP/$n,$THISTMP/${n/%$READONE$ASD.bam/$READTWO$ASD.bam}"
        done
        RUN_COMMAND="$RUN_COMMAND -format sam -illuminaPE -tbp 1"
    fi
    echo $RUN_COMMAND && eval $RUN_COMMAND
    
    mv $THISTMP/$SAMPLE"_tagdir_unfiltered" $OUTDIR
    
    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE"_tagdir_unfiltered"/tagInfo.txt
fi

################################################################################
NGSANE_CHECKPOINT_INIT "filter tagdirectory"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then


    count=`ls -1 $OUTDIR/$SAMPLE"_tagdir_unfiltered"/*tags.tsv.gz 2>/dev/null | wc -l`
    if [ $count != 0 ]; then
        gunzip $OUTDIR/$SAMPLE"_tagdir_unfiltered"/*tags.tsv.gz
    fi
    
    cp -r $OUTDIR/$SAMPLE"_tagdir_unfiltered" $THISTMP
    $GZIP $OUTDIR/$SAMPLE"_tagdir_unfiltered"/*tags.tsv
        
    RUN_COMMAND="makeTagDirectory $THISTMP/$SAMPLE"_tagdir_unfiltered" -update $HOMER_HIC_TAGDIR_ADDPARAM"
    echo $RUN_COMMAND && eval $RUN_COMMAND
       
    mv $THISTMP/$SAMPLE"_tagdir_unfiltered" $OUTDIR/$SAMPLE"_tagdir_filtered"

    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE"_tagdir_filtered"/tagInfo.txt

    [ -d $OUTDIR/$SAMPLE"_tagdir_unfiltered" ] && rm -r $OUTDIR/$SAMPLE"_tagdir_unfiltered"

fi

################################################################################
NGSANE_CHECKPOINT_INIT "output matrices"    

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    # normalize matrix for each chromosome
    LASTCHR=""
    for CHR in $(grep -v "_" $GENOME_CHROMSIZES | cut -f 1); do
        
        RUN_COMMAND="analyzeHiC $OUTDIR/"$SAMPLE"_tagdir_filtered $HOMER_HIC_NORMALIZE_ADDPARAM -chr $CHR | $GZIP -9 > $OUTDIR/"$SAMPLE"_"$CHR"_"matrix.txt.gz
        echo $RUN_COMMAND && eval $RUN_COMMAND
        LASTCHR=$CHR
    done
    
    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE"_"$LASTCHR"_"matrix.txt.gz
fi

################################################################################
NGSANE_CHECKPOINT_INIT "Significant high-res cis interactions"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    RUN_COMMAND="findHiCInteractionsByChr.pl $OUTDIR/${SAMPLE}_tagdir_filtered $HOMER_HIC_CISINTERACTIONS_ADDPARAM -cpu $CPU_HOMERHIC > $OUTDIR/${SAMPLE}_significantCisInteractions.txt"
    echo $RUN_COMMAND && eval $RUN_COMMAND
    
    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/${SAMPLE}_significantCisInteractions.txt
fi

################################################################################
NGSANE_CHECKPOINT_INIT "Significant low-res interactions"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    cat /dev/null > "$SAMPLE"_significantInteractions.txt
        
    RUN_COMMAND="analyzeHiC $OUTDIR/${SAMPLE}_tagdir_filtered $HOMER_HIC_LOWRESINTERACTIONS_ADDPARAM -interactions $OUTDIR/${SAMPLE}_significantInteractions.txt -nomatrix -cpu $CPU_HOMERHIC "
    echo $RUN_COMMAND && eval $RUN_COMMAND

    echo "$OUTDIR/${SAMPLE}_significantInteractions.txt all" > $OUTDIR/${SAMPLE}_significantInteractions.log

    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/${SAMPLE}_significantInteractions.log
fi

################################################################################
NGSANE_CHECKPOINT_INIT "Annotate interactions"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    RUN_COMMAND="annotateInteractions.pl $OUTDIR/${SAMPLE}_significantCisInteractions.txt $HOMER_HIC_ANNOTATE_ADDPARAM $OUTDIR/${SAMPLE}_annotations_cisInteractions"
    echo $RUN_COMMAND && eval $RUN_COMMAND

    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/${SAMPLE}_annotations_cisInteractions/interactionAnnotation.txt

fi

################################################################################
NGSANE_CHECKPOINT_INIT "PCA clustering"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    RUN_COMMAND="runHiCpca.pl $SAMPLE $OUTDIR/${SAMPLE}_tagdir_filtered $HOMER_HIC_PCA_ADDPARAM -cpu $CPU_HOMERHIC"
    echo $RUN_COMMAND && eval $RUN_COMMAND
    
    if hash wigToBigWig && [ -f $GENOME_CHROMSIZES ]; then   
          $OUTDIR/$SAMPLE.PC1.bedGraph ${GENOME_CHROMSIZES} $OUTDIR/$SAMPLE.PC1.bw
          rm $OUTDIR/$SAMPLE.PC1.bedGraph
    fi
    
    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE.PC1.txt

fi

################################################################################
NGSANE_CHECKPOINT_INIT "SIMA"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    if [[ "$HOMER_HIC_SIMA_ADDPARAM" =~ -p ]]; then

        RUN_COMMAND="findHiCCompartments.pl <( tail -n+2 $OUTDIR/$SAMPLE.PC1.txt) > $OUTDIR/$SAMPLE.activeDomains.txt"
        echo $RUN_COMMAND && eval $RUN_COMMAND
        
        RUN_COMMAND="findHiCCompartments.pl <( tail -n+2 $OUTDIR/$SAMPLE.PC1.txt) -opp > $OUTDIR/$SAMPLE.inactiveDomains.txt"
        echo $RUN_COMMAND && eval $RUN_COMMAND
        
        tail -n+2 $OUTDIR/$SAMPLE.activeDomains.txt > $OUTDIR/$SAMPLE.activeAndInactiveDomains.txt
        tail -n+2 $OUTDIR/$SAMPLE.inactiveDomains.txt >> $OUTDIR/$SAMPLE.activeAndInactiveDomains.txt
        
        awk '{if ($5<0){OFS="\t"; print $2,$3,$4,$1,"Active-"$6,$5}else{OFS="\t"; print $2,$3,$4,$1,"Inactive-"$6,$5}}' $OUTDIR/$SAMPLE.activeAndInactiveDomains.txt > $OUTDIR/$SAMPLE.activeAndInactiveDomains.bed
        
#        rm $OUTDIR/$SAMPLE.activeDomains.txt
#        rm $OUTDIR/$SAMPLE.inactiveDomains.txt
        
        RUN_COMMAND="SIMA.pl $OUTDIR/${SAMPLE}_tagdir_filtered -d $OUTDIR/$SAMPLE.activeAndInactiveDomains.txt $HOMER_HIC_SIMA_ADDPARAM -cpu $CPU_HOMERHIC > $SAMPLE.sima.txt"
        echo $RUN_COMMAND && eval $RUN_COMMAND
        
        # mark checkpoint
        NGSANE_CHECKPOINT_CHECK $SAMPLE.sima.txt

    else
        echo "[NOTE] no peaks provided for SIMA analysis. Skipping it"
        # mark checkpoint
        NGSANE_CHECKPOINT_CHECK
        
    fi
fi

################################################################################
NGSANE_CHECKPOINT_INIT "Circos plots (optional)"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then
    if hash circos 2>&- ; then

        RUN_COMMAND="analyzeHiC $OUTDIR/${SAMPLE}_tagdir_filtered -res 1000000 -pvalue 1e-7 -cpu $CPU_HOMERHIC -circos $SAMPLE -minDist 2000000000 -nomatrix"
        echo $RUN_COMMAND && eval $RUN_COMMAND
    fi

    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK
fi

################################################################################
# back to where we came from
cd $CURDIR

################################################################################
[ -e $OUTDIR/${SAMPLE}_significantInteractions.log.dummy ] && rm $OUTDIR/${SAMPLE}_significantInteractions.log.dummy
echo ">>>>> HiC analysis with homer - FINISHED"
echo ">>>>> enddate "`date`

