#!/bin/bash -e

# Script for ChIP-seq peak calling using MACS v2.
# It takes read alignments in .bam format.
# It produces output files: peak regions in bed format
# author: Fabian Buske
# date: August 2013

echo ">>>>> ChIPseq analysis with MACS2"
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"

function usage {
echo -e "usage: $(basename $0) -k NGSANE -f FASTQ -r REFERENCE -o OUTDIR [OPTIONS]"
exit
}

# QCVARIABLES,Resource temporarily unavailable
if [ ! $# -gt 3 ]; then usage ; fi

#INPUTS                                                                                                           
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the NGSANE repository
        -f | --bam )            shift; f=$1 ;; # bam file
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

################################################################################
CHECKPOINT="programs"

for MODULE in $MODULE_MACS2; do module load $MODULE; done  # save way to load modules that itself load other modules

export PATH=$PATH_MACS2:$PATH
module list
echo "PATH=$PATH"
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)

echo -e "--macs2      --\n "$(macs2 --version 2>&1)
[ -z "$(which macs2)" ] && echo "[ERROR] macs2 not detected" && exit 1
echo -e "--R          --\n "$(R --version | head -n 3)
[ -z "$(which R)" ] && echo "[ERROR] no R detected" && exit 1
echo -e "--convert     --\n "$(convert -version | head -n 1)
[ -z "$(which convert)" ] && echo "[WARN] imagemagick convert not detected"
echo -e "--bedToBigBed --\n "$(bedToBigBed 2>&1 | tee | head -n 1 )
[ -z "$(which bedToBigBed)" ] && echo "[WARN] bedToBigBed not detected, cannot compress bedgraphs"

echo -e "\n********* $CHECKPOINT"
################################################################################
CHECKPOINT="parameters"

# get basename of f
n=${f##*/}
c=${CHIPINPUT##*/}

if [ -z "$CHIPINPUT" ] || [ ! -f $CHIPINPUT ]; then
    echo "[WARN] input control not provided or invalid (CHIPINPUT)"
    unset CHIPINPUT
else
    CHIPINPUT="--control $SOURCE/$CHIPINPUT"
fi

echo -e "\n********* $CHECKPOINT"
################################################################################
CHECKPOINT="recall files from tape"

if [ -n "$DMGET" ]; then
	dmget -a ${f}
	[ -n $CHIPINPUT ] && dmget -a $CHIPINPUT
fi

cd $MYOUT

echo -e "\n********* $CHECKPOINT"

################################################################################
CHECKPOINT="macs 2 - call peaks "

if [[ -n "$RECOVERFROM" ]] && [[ $(grep "********* $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else

    macs2 callpeak $MACS2_CALLPEAK_ADDPARAM --treatment $f $CHIPINPUT --gsize $MACS2_GENOMESIZE --name ${n/.$ASD.bam/} > ${n/.$ASD.bam/}.summary.txt 2>&1

    Rscript ${n/.$ASD.bam/}_model.r
    cat ${n/.$ASD.bam/}_model.r | sed 's/pdf(/png(/' | sed 's/=6/=200/g' > ${n/.$ASD.bam/}_model.png.r
    Rscript ${n/.$ASD.bam/}_model.png.r && rm ${n/.$ASD.bam/}_model.png.r
    if [ -n "$(which convert)" ]; then convert -size 200x200 ${n/.$ASD.bam/}_model.pdf ${n/.$ASD.bam/}_model.png; fi
   
    if [ -n "$CHIPINPUT" ]; then
        macs2 bdgcmp $MACS2_BDGCMP_ADDPARAM --treatment ${n/.$ASD.bam/_treat_pileup.bdg} ${n/.$ASD.bam/_control_lambda.bdg} --output ${n/.$ASD.bam/}
	if [ -n "$(which bedToBigBed)" ]; then 
            bedToBigBed -type=bed4 ${n/.$ASD.bam/_treat_pileup.bdg} $CHROM_SIZES ${n/.$ASD.bam/_treat_pileup.bb}
            bedToBigBed -type=bed4 ${n/.$ASD.bam/_control_lambda.bdg} $CHROM_SIZES ${n/.$ASD.bam/_control_lambda.bdg}
            rm ${n/.$ASD.bam/_treat_pileup.bdg} ${n/.$ASD.bam/_control_lambda.bdg}
        fi

    else
        if [ -n "$(which bedToBigBed)" ]; then
            bedToBigBed -type=bed4 ${n/.$ASD.bam/_treat_pileup.bdg} $CHROM_SIZES ${n/.$ASD.bam/_treat_pileup.bb}
            rm ${n/.$ASD.bam/_treat_pileup.bdg}
        fi
    fi

    # mark checkpoint
    [ -f ${n/.$ASD.bam/}_model.r ] && echo -e "\n********* $CHECKPOINT" && unset RECOVERFROM
fi

################################################################################
CHECKPOINT="macs 2 - refine peaks "

if [[ -n "$RECOVERFROM" ]] && [[ $(grep "********* $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else

    if [[ "$MACS2_CALLPEAK_ADDPARAM" == *--broad* ]]; then
        echo "[NOTE] convert broadpeaks to 6 bed file for refining"
        cat ${n/.$ASD.bam/}_peaks.broadPeak | cut -f1-6  > ${n/.$ASD.bam/}_peaks.bed
    else
        echo "[NOTE] convert narrowpeaks to 6 bed filei for refining"
        cat ${n/.$ASD.bam/}_peaks.narrowPeak | cut -f1-6 > ${n/.$ASD.bam/}_peaks.bed
    fi

    macs2 refinepeak $MACS2_REFINEPEAK_ADDPARAM -b ${n/.$ASD.bam/}_peaks.bed -i $f --o-prefix ${n/.$ASD.bam/}_refined  >> ${n/.$ASD.bam/}.summary.txt 2>&1
    
    # mark checkpoint
    [ -f ${n/.$ASD.bam/}_refined_peaks.bed ] && echo -e "\n********* $CHECKPOINT" && unset RECOVERFROM
fi

################################################################################
CHECKPOINT="zip"

# $GZIP $MYOUT/${n/.$ASD.bam/}-${c/.$ASD.bam/}_details
echo "Peaks: `cat ${n/.$ASD.bam/}_peaks.bed | wc -l | awk '{print $1}'`" >> ${n/.$ASD.bam/}.summary.txt
echo "Summits: `cat ${n/.$ASD.bam/}_summits.bed | wc -l | awk '{print $1}'`" >> ${n/.$ASD.bam/}.summary.txt

cd $SOURCE

echo -e "\n********* $CHECKPOINT"
################################################################################
echo ">>>>> ChIPseq analysis with MACS2 - FINISHED"
echo ">>>>> enddate "`date`

