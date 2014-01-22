#!/bin/bash -e

# Script for ChIP-seq peak calling using MACS v2.
# It takes read alignments in .bam format.
# It produces output files: peak regions in bed format
# author: Fabian Buske
# date: August 2013

# QCVARIABLES,Resource temporarily unavailable
# RESULTFILENAME <DIR>/<TASK>/<SAMPLE>_refinepeak.bed


echo ">>>>> ChIPseq analysis with MACS2"
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"

function usage {
echo -e "usage: $(basename $0) -k NGSANE -f FASTQ -o OUTDIR [OPTIONS]"
exit
}


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

for MODULE in $MODULE_MACS2; do module load $MODULE; done  # save way to load modules that itself load other modules

export PATH=$PATH_MACS2:$PATH
module list
echo "PATH=$PATH"
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--macs2       --\n "$(macs2 --version 2>&1)
[ -z "$(which macs2)" ] && echo "[ERROR] macs2 not detected" && exit 1
echo -e "--R           --\n "$(R --version | head -n 3)
[ -z "$(which R)" ] && echo "[ERROR] no R detected" && exit 1
echo -e "--convert     --\n "$(convert -version | head -n 1)
[ -z "$(which convert)" ] && echo "[WARN] imagemagick convert not detected"
echo -e "--bedToBigBed --\n "$(bedToBigBed 2>&1 | tee | head -n 1 )
[ -z "$(which bedToBigBed)" ] && echo "[WARN] bedToBigBed not detected, cannot compress bedgraphs"

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="parameters"

# get basename of f
f=${f/%.dummy/} #if input came from pip
n=${f##*/}
SAMPLE=${n/%.$ASD.bam/}
c=${CHIPINPUT##*/}
CONTROL=${c/%.$ASD.bam/}

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
# set default method to ppois unless specified
if [ -z "$MACS2_BDGCMP_METHOD" ]; then
    echo "[NOTE] no method provided for MACS2_BDGCMP_METHOD, defaulting to ppois"
    MACS2_BDGCMP_METHOD="ppois"
fi

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="recall files from tape"

if [ -n "$DMGET" ]; then
	dmget -a $(dirname $FASTA)/*
	dmget -a ${f}
	dmget -a $OUTDIR/*
	[ -n "$CHIPINPUT" ] && dmget -a $CHIPINPUT
fi

if [ -z "$CHIPINPUT" ] || [ ! -f $CHIPINPUT ]; then
    echo "[WARN] input control not provided or invalid (CHIPINPUT)"
    unset CHIPINPUT
else
	echo "[NOTE] CHIPINPUT $CHIPINPUT"
    CHIPINPUT="--control $SOURCE/$CHIPINPUT"
fi

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="macs 2 - call peaks "

cd $OUTDIR
if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else

    if [ -n "$MACS2_MAKEBIGBEDS" ]; then
        MAKEBEDGRAPH="--bdg "
    fi
    
    RUN_COMMAND="macs2 callpeak $MACS2_CALLPEAK_ADDPARAM $MACS2_MAKEBIGBED $MAKEBEDGRAPH --treatment $f $CHIPINPUT --gsize $MACS2_GENOMESIZE --name $SAMPLE > $SAMPLE.summary.txt 2>&1"    
    echo $RUN_COMMAND && eval $RUN_COMMAND

    if [ -f $SAMPLE"_"model.r ];then 
        Rscript $SAMPLE"_"model.r
        if hash convert; then 
            convert -format png $SAMPLE"_"model.pdf $SAMPLE"_"model.png
        fi
    fi
    
    # mark checkpoint
    if [ -f $SAMPLE"_"peaks.xls ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
fi
################################################################################
CHECKPOINT="macs 2 - check fragment length "

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else

    if [ "$(grep '#2 predicted fragment length is -' $SAMPLE.summary.txt)" > 0 ]; then
        echo "[NOTE] Sample has negative fragment length, adjusting shiftsize with alternatives if available"

        [ -f $SAMPLE"_"peaks.xls ] && rm $SAMPLE"_"peaks.xls
    
        if [ "$(grep '#2 alternative fragment length(s) may be' $SAMPLE.summary.txt | egrep -o '(-[0-9]{2,4},[0-9]{2,4})' | cut -d ',' -f 2)" > 0 ]; then
            SHIFTSIZE="--shiftsize $(grep '#2 alternative fragment length(s) may be' $SAMPLE.summary.txt | egrep -o '(-[0-9]{2,4},[0-9]{2,4})' | cut -d ',' -f 2)"           
        fi
        
        RUN_COMMAND="macs2 callpeak $MACS2_CALLPEAK_ADDPARAM $MACS2_MAKEBIGBED $MAKEBEDGRAPH --nomodel $SHIFTSIZE --treatment $f $CHIPINPUT --gsize $MACS2_GENOMESIZE --name $SAMPLE > $SAMPLE.summary.txt 2>&1"    
        echo $RUN_COMMAND && eval $RUN_COMMAND

    fi
    
    # mark checkpoint
    if [ -f $SAMPLE"_"peaks.xls ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
fi
################################################################################
CHECKPOINT="macs 2 - convert bedgraph to bigbed"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else
    
    if [ -z "$MACS2_MAKEBIGBEDS" ]; then
        echo -e "[NOTE] skipping bigbed generation"
        echo -e "\n********* $CHECKPOINT\n"
        
    else
        if [ -n "$CHIPINPUT" ]; then
    
            RUN_COMMAND="macs2 bdgcmp $MACS2_BDGCMP_ADDPARAM --method $MACS2_BDGCMP_METHOD --tfile $SAMPLE"_"treat_pileup.bdg --cfile $SAMPLE"_"control_lambda.bdg --output $SAMPLE >> $SAMPLE.summary.txt 2>&1"
            echo $RUN_COMMAND && eval $RUN_COMMAND
    
        	if hash bedToBigBed ; then 
                bedToBigBed -type=bed4 $SAMPLE"_"treat_pileup.bdg $GENOME_CHROMSIZES $SAMPLE"_"treat_pileup.bb
                bedToBigBed -type=bed4 ${n/.$ASD.bam/"_"$MACS2_BDGCMP_METHOD.bdg} $GENOME_CHROMSIZES ${n/.$ASD.bam/"_"$MACS2_BDGCMP_METHOD.bb}
                # don't create another one for the control in case if already exists
                if [ ! -f $SAMPLE"_"control_lambda.bb ]; then 
                    BBTMP=$RANDOM
                    bedToBigBed -type=bed4 $SAMPLE"_"control_lambda.bdg $GENOME_CHROMSIZES $TMP/$SAMPLE"_"control_lambda.bb$BBTMP
                    mv $TMP/$SAMPLE"_"control_lambda.bb$BBTMP $SAMPLE"_"control_lambda.bb
                fi
            fi
    
        else
        	if hash bedToBigBed ; then 
                bedToBigBed -type=bed4 $SAMPLE"_"treat_pileup.bdg $GENOME_CHROMSIZES $SAMPLE"_"treat_pileup.bb
            fi
        fi
        
        # mark checkpoint
        if [ -f $SAMPLE"_"treat_pileup.bb ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
    fi 
    
    # cleanup
    [ -e $SAMPLE"_"treat_pileup.bdg ] && rm $SAMPLE"_"treat_pileup.bdg 
    [ -e $SAMPLE"_"control_lambda.bdg ] && rm $SAMPLE"_"control_lambda.bdg 
    [ -e $SAMPLE"_"$MACS2_BDGCMP_METHOD.bdg ] && rm $SAMPLE"_"$MACS2_BDGCMP_METHOD.bdg

fi
################################################################################
CHECKPOINT="macs 2 - refine peaks "

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else

    if [[ "$MACS2_CALLPEAK_ADDPARAM" == *--broad* ]]; then
        echo "[NOTE] convert broadpeaks to 6 bed file for refining"
        cat $SAMPLE"_"peaks.broadPeak | cut -f1-6  > $SAMPLE"_"peaks.bed
    else
        echo "[NOTE] convert narrowpeaks to 6 bed file for refining"
        cat $SAMPLE"_"peaks.narrowPeak | cut -f1-6 > $SAMPLE"_"peaks.bed
    fi

    RUN_COMMAND="macs2 refinepeak $MACS2_REFINEPEAK_ADDPARAM -b $SAMPLE"_"peaks.bed -i $f --o-prefix $SAMPLE  >> $SAMPLE.summary.txt 2>&1"
    echo $RUN_COMMAND && eval $RUN_COMMAND
    [ -e $SAMPLE"_"peaks.bed ] && rm $SAMPLE"_"peaks.bed    

    echo "Final number of refined peaks: $(wc -l ${SAMPLE}_refinepeak.bed )" >> $SAMPLE.summary.txt
 
    # mark checkpoint
    if [ -f ${SAMPLE}_refinepeak.bed ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
fi
################################################################################
# back to where we came from
cd $SOURCE
################################################################################
[ -e $OUTDIR/$SAMPLE"_"refinepeak.bed.dummy ] && rm $OUTDIR/$SAMPLE"_"refinepeak.bed.dummy
echo ">>>>> ChIPseq analysis with MACS2 - FINISHED"
echo ">>>>> enddate "`date`

