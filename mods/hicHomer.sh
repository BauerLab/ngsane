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

for MODULE in $MODULE_HOMERHIC; do module load $MODULE; done  # save way to load modules that itself load other modules

export PATH=$PATH_HOMERHIC:$PATH
module list
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

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="parameters"

# get basename of f
if [ -z "$POOLED_DATA_NAME" ]; then 
    n=${f##*/}
    SAMPLE=${n/%$READONE.$ASD.bam/}
    
    #is paired ?                                                                                                      
    if [ "$f" != "${f/%$READONE.$ASD.bam/$READTWO.$ASD.bam}" ] && [ -e ${f/%$READONE.$ASD.bam/$READTWO.$ASD.bam} ]; then
        PAIRED="1"
    else
        PAIRED="0"
        echo "[ERROR] paired library required for HIC analysis" && exit 1
    fi

else
    SAMPLE=$POOLED_DATA_NAME
    array=(${f//,/ })
    #is paired ?                                                                                                      
    if [ "${array[i]}" != "${array[i]/%$READONE.$ASD.bam/$READTWO.$ASD.bam}" ] && [ -e ${array[i]/%$READONE.$ASD.bam/$READTWO.$ASD.bam} ]; then
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

THISTMP=$TMP"/"$(whoami)"/"$(echo $OUTDIR/$SAMPLE | md5sum | cut -d' ' -f1)
[ -d $THISTMP ] && rm -r $THISTMP
mkdir -p $THISTMP

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="recall files from tape"

if [ -n "$DMGET" ]; then
	dmget -a ${f/$READONE/"*"}
	dmget -a $OUTDIR/*
fi

echo -e "\n********* $CHECKPOINT\n"
################################################################################

#homer likes to write in the current directory, so change to target
CURDIR=$(pwd)
cd $OUTDIR

################################################################################
CHECKPOINT="create unfiltered tagdirectory"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else

    if [ -z "$POOLED_DATA_NAME" ]; then 
        cp $f ${f/%$READONE.$ASD.bam/$READTWO.$ASD.bam} $THISTMP
        RUN_COMMAND="makeTagDirectory $THISTMP/$SAMPLE"_tagdir_unfiltered" $THISTMP/$n,$THISTMP/${n/%$READONE.$ASD.bam/$READTWO.$ASD.bam} -format sam -illuminaPE -tbp 1"
    else
        # pool data
        RUN_COMMAND="makeTagDirectory $THISTMP/$SAMPLE"_tagdir_unfiltered
        array=(${f//,/ })
        for i in "${!array[@]}"; do
            cp ${array[i]} ${array[i]/%$READONE.$ASD.bam/$READTWO.$ASD.bam} $THISTMP
            n=${array[i]##*/} 
            RUN_COMMAND="$RUN_COMMAND $THISTMP/$n,$THISTMP/${n/%$READONE.$ASD.bam/$READTWO.$ASD.bam}"
        done
        RUN_COMMAND="$RUN_COMMAND -format sam -illuminaPE -tbp 1"
    fi
    echo $RUN_COMMAND && eval $RUN_COMMAND
    
    mv $THISTMP/$SAMPLE"_tagdir_unfiltered" $OUTDIR
    
    # mark checkpoint
    if [ -e $OUTDIR/$SAMPLE"_tagdir_unfiltered"/tagInfo.txt ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
fi

################################################################################
CHECKPOINT="filter tagdirectory"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else

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
    if [ -e $OUTDIR/$SAMPLE"_tagdir_filtered"/tagInfo.txt ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

    [ -d $OUTDIR/$SAMPLE"_tagdir_unfiltered" ] && rm -r $OUTDIR/$SAMPLE"_tagdir_unfiltered"

fi

################################################################################
CHECKPOINT="normalize matrices"    

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

    if [ "$HOMER_HIC_INTERACTIONS" == "all" ]; then
        RUN_COMMAND="analyzeHiC $OUTDIR/"$SAMPLE"_tagdir_filtered $HOMER_HIC_NORMALIZE_ADDPARAM  > $OUTDIR/"$SAMPLE"_matrix.txt"
        echo $RUN_COMMAND && eval $RUN_COMMAND
    
    [ ! -f $FASTA.fai ] && samtools faidx $FASTA
    
    elif [ "$HOMER_HIC_INTERACTIONS" == "cis" ]; then
    
        for CHR in $(awk '{print $1'} $FASTA.fai); do
    	    RUN_COMMAND="analyzeHiC $OUTDIR/"$SAMPLE"_tagdir_filtered $HOMER_HIC_NORMALIZE_ADDPARAM -chr $CHR > $OUTDIR/"$SAMPLE"_${CHR}_matrix.txt"
    	    echo $RUN_COMMAND && eval $RUN_COMMAND
        done
    elif [ "$HOMER_HIC_INTERACTIONS" == "trans" ]; then
    
        for CHR1 in $(awk '{print $1'} $FASTA.fai); do
            for CHR2 in $(awk '{print $1'} $FASTA.fai); do
                if [ "$CHR1" != "$CHR2" ]; then
                    RUN_COMMAND="analyzeHiC $OUTDIR/"$SAMPLE"_tagdir_filtered $HOMER_HIC_NORMALIZE_ADDPARAM -chr $CHR1 -chr2 $CHR2 > $OUTDIR/"$SAMPLE"_${CHR1}-${CHR2}_matrix.txt"
                    echo $RUN_COMMAND && eval $RUN_COMMAND
                fi
            done
        done
    fi

    # mark checkpoint
    echo -e "\n********* $CHECKPOINT\n"
fi

################################################################################
CHECKPOINT="Significant high-res cis interactions"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

    RUN_COMMAND="findHiCInteractionsByChr.pl $OUTDIR/${SAMPLE}_tagdir_filtered $HOMER_HIC_CISINTERACTIONS_ADDPARAM -cpu $CPU_HOMERHIC > $OUTDIR/${SAMPLE}_significantCisInteractions.txt"
    echo $RUN_COMMAND && eval $RUN_COMMAND
    
    # mark checkpoint
    if [ -f $OUTDIR/${SAMPLE}_significantCisInteractions.txt ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
fi

################################################################################
CHECKPOINT="Significant low-res interactions"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
    cat /dev/null > "$SAMPLE"_significantInteractions.txt
        
    RUN_COMMAND="analyzeHiC $OUTDIR/${SAMPLE}_tagdir_filtered $HOMER_HIC_LOWRESINTERACTIONS_ADDPARAM -interactions $OUTDIR/${SAMPLE}_significantInteractions.txt -nomatrix -cpu $CPU_HOMERHIC "
    echo $RUN_COMMAND && eval $RUN_COMMAND

    echo "$OUTDIR/${SAMPLE}_significantInteractions.txt all" > $OUTDIR/${SAMPLE}_significantInteractions.log

    # mark checkpoint
    if [[ $(wc -l $OUTDIR/${SAMPLE}_significantInteractions.log | cut -d' ' -f 1) -ge 1 ]];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
fi

################################################################################
CHECKPOINT="Annotate interactions"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

    RUN_COMMAND="annotateInteractions.pl $OUTDIR/${SAMPLE}_significantCisInteractions.txt $HOMER_HIC_ANNOTATE_ADDPARAM $OUTDIR/${SAMPLE}_annotations_cisInteractions"
    echo $RUN_COMMAND && eval $RUN_COMMAND

    # mark checkpoint
    if [ -f $OUTDIR/${SAMPLE}_annotations_cisInteractions/interactionAnnotation.txt ]; then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi

################################################################################
CHECKPOINT="PCA clustering"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

    RUN_COMMAND="runHiCpca.pl $SAMPLE $OUTDIR/${SAMPLE}_tagdir_filtered $HOMER_HIC_PCA_ADDPARAM -cpu $CPU_HOMERHIC"
    echo $RUN_COMMAND && eval $RUN_COMMAND
    
    # mark checkpoint
    if [ -f $OUTDIR/$SAMPLE.PC1.txt ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi


################################################################################
CHECKPOINT="SIMA"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

    if [[ "$HOMER_HIC_SIMA_ADDPARAM" ~= "-p" ]]; then

        RUN_COMMAND="findHiCCompartments.pl $OUTDIR/$SAMPLE.PC1.txt > $OUTDIR/$SAMPLE.activeDomains.txt"
        echo $RUN_COMMAND && eval $RUN_COMMAND
        
        RUN_COMMAND="findHiCCompartments.pl -opp $OUTDIR/$SAMPLE.PC1.txt > $OUTDIR/$SAMPLE.inactiveDomains.txt"
        echo $RUN_COMMAND && eval $RUN_COMMAND
        
        RUN_COMMAND="SIMA.pl $OUTDIR/${SAMPLE}_tagdir_filtered -d <(cat $OUTDIR/$SAMPLE.activeDomains.txt $OUTDIR/$SAMPLE.inactiveDomains.txt)$HOMER_HIC_SIMA_ADDPARAM -cpu $CPU_HOMERHIC > $SAMPLE.sima.txt"
        echo $RUN_COMMAND && eval $RUN_COMMAND
        
        # mark checkpoint
        if [ -f $SAMPLE.sima.txt ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

    else
        echo "[NOTE] no peaks provided for SIMA analysis. Skipping it"
        # mark checkpoint
        echo -e "\n********* $CHECKPOINT\n"
        
    fi
fi

################################################################################
CHECKPOINT="Circos plots (optional)"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
    if hash circos 2>&- ; then

        RUN_COMMAND="analyzeHiC $OUTDIR/${SAMPLE}_tagdir_filtered -res 1000000 -pvalue 1e-7 -cpu $CPU_HOMERHIC -circos $SAMPLE -minDist 2000000000 -nomatrix"
        echo $RUN_COMMAND && eval $RUN_COMMAND
    fi

    # mark checkpoint
    echo -e "\n********* $CHECKPOINT\n"
fi

################################################################################
# back to where we came from
cd $CURDIR

################################################################################
[ -e $OUTDIR/${SAMPLE}_significantInteractions.log.dummy ] && rm $OUTDIR/${SAMPLE}_significantInteractions.log.dummy
echo ">>>>> HiC analysis with homer - FINISHED"
echo ">>>>> enddate "`date`

