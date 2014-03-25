#!/bin/bash -e

# Script for ChIP-seq peak calling using MACS v2.
# It takes read alignments in .bam format.
# It produces output files: peak regions in bed format
# author: Fabian Buske
# date: August 2013

# QCVARIABLES,Resource temporarily unavailable
# RESULTFILENAME <DIR>/<TASK>/<SAMPLE>.refinepeak.bed


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
CHECKPOINT="macs 2 - predictd "

cd $OUTDIR
if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else
    
    if [ -n "$MACS2_FRAGMENTSIZE" ]; then
        echo "[NOTE] Skip modeling"
    
    else
    
        RUN_COMMAND="macs2 predictd $MACS2_PREDICTD_ADDPARAM --ifile $f --gsize $MACS2_GENOMESIZE --rfile $SAMPLE > $SAMPLE.summary.txt 2>&1"    
        echo $RUN_COMMAND && eval $RUN_COMMAND
    
        Rscript $SAMPLE"_"model.R
        if hash convert; then 
            convert -format png $SAMPLE"_"model.pdf $SAMPLE"_"model.png
        fi
    
    fi
    
    # mark checkpoint
    if [ -f $SAMPLE"_"model.R ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
fi
################################################################################
CHECKPOINT="macs 2 - call peaks "

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else

    if [ -z "$MACS2_FRAGMENTSIZE" ]; then
        MACS2_FRAGMENTSIZE=$(grep 'alternative fragment length(s) may be' $SAMPLE.summary.txt | sed 's/.* be //' | egrep -o '((-[0-9]{2,4},)?[0-9]{2,4})')
    fi
    
    RUN_COMMAND="macs2 callpeak $MACS2_CALLPEAK_ADDPARAM --nomodel --extsize $MACS2_FRAGMENTSIZE --treatment $f $CHIPINPUT --gsize $MACS2_GENOMESIZE --name $SAMPLE >> $SAMPLE.summary.txt 2>&1"    
    echo $RUN_COMMAND && eval $RUN_COMMAND
    
    # mark checkpoint
    if [ -f $SAMPLE"_"peaks.xls ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
fi
################################################################################
CHECKPOINT="macs 2 - refine peaks "

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else

    if [[ "$MACS2_CALLPEAK_ADDPARAM" == *--broad* ]]; then
        echo "[NOTE] convert broadpeaks to 6 bed file for refining"
        cat $SAMPLE"_"peaks.broadPeak | awk '{OFS="\t"; print $1,$2,$3,$4,$9,$6}'  > $SAMPLE"_"peaks.bed
    else
        echo "[NOTE] convert narrowpeaks to 6 bed file for refining"
        cat $SAMPLE"_"peaks.narrowPeak | awk '{OFS="\t"; print $1,$2,$3,$4,$9,$6}' > $SAMPLE"_"peaks.bed
    fi
    
    RUN_COMMAND="macs2 refinepeak $MACS2_REFINEPEAK_ADDPARAM -b $SAMPLE"_"peaks.bed -i $f --o-prefix $SAMPLE >> $SAMPLE.summary.txt 2>&1"
    echo $RUN_COMMAND && eval $RUN_COMMAND
    mv ${SAMPLE}_refinepeak.bed ${SAMPLE}.refinepeak.bed

	if hash bedToBigBed ; then 
        echo "[NOTE] create bigbed from peaks" 
        awk '{OFS="\t"; print $1,$2,$3,$5}' $SAMPLE"_"peaks.bed > $SAMPLE"_"peaks.tmp
        bedToBigBed -type=bed4 $SAMPLE"_"peaks.tmp $GENOME_CHROMSIZES $SAMPLE.bb
        rm $SAMPLE"_"peaks.tmp
    fi
    [ -e $SAMPLE"_"peaks.bed ] && rm $SAMPLE"_"peaks.bed    

    echo "Final number of refined peaks: $(wc -l ${SAMPLE}.refinepeak.bed )" >> $SAMPLE.summary.txt

    # mark checkpoint
    if [ -f ${SAMPLE}.refinepeak.bed ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
fi
################################################################################
# back to where we came from
cd $SOURCE
################################################################################
[ -e $OUTDIR/$SAMPLE"_"refinepeak.bed.dummy ] && rm $OUTDIR/$SAMPLE"_"refinepeak.bed.dummy
echo ">>>>> ChIPseq analysis with MACS2 - FINISHED"
echo ">>>>> enddate "`date`

