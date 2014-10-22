#!/bin/bash -e

# Script for ChIP-seq peak calling using MACS v2.
# It takes read alignments in .bam format.
# It produces output files: peak regions in bed format
# author: Fabian Buske
# date: August 2013

# QCVARIABLES,Resource temporarily unavailable
# RESULTFILENAME <DIR>/<TASK>/<SAMPLE>_summits.bed


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
hash module 2>/dev/null && for MODULE in $MODULE_MACS2; do module load $MODULE; done && module list 

export PATH=$PATH_MACS2:$PATH
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

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "parameters"

# get basename of f
n=${f##*/}
SAMPLE=${n/%$ASD.bam/}
c=${CHIPINPUT##*/}
CONTROL=${c/%$ASD.bam/}

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

if [[ "$MACS2_CALLPEAK_ADDPARAM" == *--broad* ]]; then
    PEAKSUFFIX="_peaks.broadPeak"
else
    PEAKSUFFIX="_peaks.narrowPeak"
fi

if [ -z "$NGSANE_RECOVERFROM" ]; then
    [ -e $OUTDIR/$SAMPLE$PEAKSUFFIX ] && rm $OUTDIR/$SAMPLE*
fi

cd $OUTDIR  

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "recall files from tape"

if [ -n "$DMGET" ]; then
	dmget -a $(dirname $FASTA)/*
	dmget -a ${f}
	dmget -a $OUTDIR/*
	[ -n "$CHIPINPUT" ] && dmget -a $CHIPINPUT
fi

if [[ -z "$CHIPINPUT" || ! -f $CHIPINPUT ]]; then
    echo "[WARN] input control not provided or invalid ($CHIPINPUT)"
    unset CHIPINPUT
else
	echo "[NOTE] CHIPINPUT $CHIPINPUT"
    CHIPINPUT="--control $CHIPINPUT"
fi

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "macs 2 - predictd "

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then
 
    if [ -n "$MACS2_FRAGMENTSIZE" ]; then
        echo "[NOTE] Skip modeling"
        cat /dev/null > $SAMPLE.summary.txt
        # mark checkpoint
        NGSANE_CHECKPOINT_CHECK
        
    else
    
        RUN_COMMAND="macs2 predictd $MACS2_PREDICTD_ADDPARAM --ifile $f --gsize $MACS2_GENOMESIZE --rfile $SAMPLE.R > $SAMPLE.summary.txt 2>&1"    
        echo $RUN_COMMAND && eval $RUN_COMMAND
        
        Rscript $SAMPLE.R
        if hash convert; then 
            convert -format png $SAMPLE.R"_"model.pdf $SAMPLE.R"_"model.png
        fi
    
        # mark checkpoint
        NGSANE_CHECKPOINT_CHECK $SAMPLE.R

    fi    
fi
################################################################################
NGSANE_CHECKPOINT_INIT "macs 2 - call peaks "

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    if [ -z "$MACS2_FRAGMENTSIZE" ]; then
        # estimated fragment sizes up to the tagsize are usually artifacts
        TAGSIZE=$(fgrep '# tag size =' $SAMPLE.summary.txt | cut -d'=' -f 2 | tr -d ' \t\r\f')
        
        MACS2_FRAGMENTSIZE=$(grep 'alternative fragment length(s) may be' $SAMPLE.summary.txt | sed 's/.* be //' | cut -d' ' -f 1 | tr ',' '\n' | awk -v t=$TAGSIZE '($1>t){print $1; exit;}')
        
        if [[ -z "$MACS2_FRAGMENTSIZE" || $MACS2_FRAGMENTSIZE -le 0 ]]; then
            echo "[ERROR] could not determine fragments size, please provide the parameter MACS2_FRAGMENTSIZE."
            exit 1
        fi
    fi
    echo "Determined fragment size: $MACS2_FRAGMENTSIZE" >> $SAMPLE.summary.txt
    echo "ChIP input: $CONTROL" >> $SAMPLE.summary.txt

    RUN_COMMAND="macs2 callpeak $MACS2_CALLPEAK_ADDPARAM --nomodel --extsize $MACS2_FRAGMENTSIZE --treatment $f $CHIPINPUT --gsize $MACS2_GENOMESIZE --name $SAMPLE >> $SAMPLE.summary.txt 2>&1"    
    echo $RUN_COMMAND && eval $RUN_COMMAND

    # remove excel stuff
    [ -f $SAMPLE"_peaks.xls" ] && rm $SAMPLE"_peaks.xls"

    # slim down peak file by removing samplenames from column 4
    sed -i "s/${SAMPLE}_peak/peak/g" $SAMPLE$PEAKSUFFIX
    
    echo "Nucleotides covered: $(awk '{sum+=$3-$2}END{sum>0?sum=sum:sum=0; print sum}' $SAMPLE$PEAKSUFFIX)" >> $SAMPLE.summary.txt

    # create bigbed making sure score is between 0 and 1000
	if hash bedToBigBed; then 
	   if [ -s $SAMPLE$PEAKSUFFIX ]; then
            echo "[NOTE] create bigbed from peaks" 
            sort -k1,1 -k2,2n $SAMPLE$PEAKSUFFIX | awk '{OFS="\t";$5>1000?$5=1000:$5=$5; print $0}' > $SAMPLE$PEAKSUFFIX.tmp
            if [[ "$MACS2_CALLPEAK_ADDPARAM" == *--broad* ]]; then
                bedToBigBed -type=bed6+3 $SAMPLE$PEAKSUFFIX.tmp $GENOME_CHROMSIZES $SAMPLE.bb
            else
                bedToBigBed -type=bed6+4 $SAMPLE$PEAKSUFFIX.tmp $GENOME_CHROMSIZES $SAMPLE.bb
            fi
            rm $SAMPLE$PEAKSUFFIX.tmp
        fi
    fi

    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $SAMPLE$PEAKSUFFIX
fi
################################################################################
NGSANE_CHECKPOINT_INIT "macs 2 - refine peaks "

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    RUN_COMMAND="macs2 refinepeak $MACS2_REFINEPEAK_ADDPARAM -b $SAMPLE$PEAKSUFFIX -i $f --o-prefix $SAMPLE >> $SAMPLE.summary.txt 2>&1"
    echo $RUN_COMMAND && eval $RUN_COMMAND
    
    # overwrite old summits
    mv $SAMPLE"_refinepeak.bed" $SAMPLE"_summits.bed"
        
    echo "Final number of refined summits: $(wc -l ${SAMPLE}_summits.bed | awk '{print $1}')" >> $SAMPLE.summary.txt
        
    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK 
fi
################################################################################
[ -e $OUTDIR/${SAMPLE}_summits.bed.dummy ] && rm $OUTDIR/${SAMPLE}_summits.bed.dummy
echo ">>>>> ChIPseq analysis with MACS2 - FINISHED"
echo ">>>>> enddate "`date`

