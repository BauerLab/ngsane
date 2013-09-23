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
n=${f##*/}
c=${CHIPINPUT##*/}

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

    RUN_COMMAND="macs2 callpeak $MACS2_CALLPEAK_ADDPARAM $MACS2_MAKEBIGBED --bdg --treatment $f $CHIPINPUT --gsize $MACS2_GENOMESIZE --name ${n/.$ASD.bam/} > ${n/.$ASD.bam/}.summary.txt 2>&1"
    echo $RUN_COMMAND && eval $RUN_COMMAND

    if [ -f ${n/.$ASD.bam/}_model.r ];then 
    
        Rscript ${n/.$ASD.bam/}_model.r
        if [ hash convert ]; then 
            convert -format png ${n/.$ASD.bam/_model.pdf} ${n/.$ASD.bam/_model.png}
        fi
    fi
    
    # mark checkpoint
    if [ -f ${n/.$ASD.bam/}_peaks.xls ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi

################################################################################
CHECKPOINT="macs 2 - convert bedgraph to bigbed"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else
    
    if [ -n "$CHIPINPUT" ]; then

        RUN_COMMAND="macs2 bdgcmp $MACS2_BDGCMP_ADDPARAM --method $MACS2_BDGCMP_METHOD --tfile ${n/.$ASD.bam/_treat_pileup.bdg} --cfile ${n/.$ASD.bam/_control_lambda.bdg} --output ${n/.$ASD.bam/} >> ${n/.$ASD.bam/}.summary.txt 2>&1"
        echo $RUN_COMMAND && eval $RUN_COMMAND

    	if hash bedToBigBed ; then 
            bedToBigBed -type=bed4 ${n/.$ASD.bam/_treat_pileup.bdg} $GENOME_CHROMSIZES ${n/.$ASD.bam/_treat_pileup.bb}
            bedToBigBed -type=bed4 ${n/.$ASD.bam/_control_lambda.bdg} $GENOME_CHROMSIZES ${n/.$ASD.bam/_control_lambda.bb}
            bedToBigBed -type=bed4 ${n/.$ASD.bam/_$MACS2_BDGCMP_METHOD.bdg} $GENOME_CHROMSIZES ${n/.$ASD.bam/_$MACS2_BDGCMP_METHOD.bb}
        fi

    else
    	if hash bedToBigBed ; then 
            bedToBigBed -type=bed4 ${n/.$ASD.bam/_treat_pileup.bdg} $GENOME_CHROMSIZES ${n/.$ASD.bam/_treat_pileup.bb}
        fi
    fi

    [ -e ${n/.$ASD.bam/_treat_pileup.bdg} ] && rm ${n/.$ASD.bam/_treat_pileup.bdg} 
    [ -e ${n/.$ASD.bam/_control_lambda.bdg} ] && rm ${n/.$ASD.bam/_control_lambda.bdg} 
    [ -e ${n/.$ASD.bam/_$MACS2_BDGCMP_METHOD.bdg} ] && rm ${n/.$ASD.bam/_$MACS2_BDGCMP_METHOD.bdg}

    # mark checkpoint
    if [ -f ${n/.$ASD.bam/_treat_pileup.bb} ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi

################################################################################
CHECKPOINT="macs 2 - refine peaks "

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else

    if [[ "$MACS2_CALLPEAK_ADDPARAM" == *--broad* ]]; then
        echo "[NOTE] convert broadpeaks to 6 bed file for refining"
        cat ${n/.$ASD.bam/}_peaks.broadPeak | cut -f1-6  > ${n/.$ASD.bam/}_peaks.bed
    else
        echo "[NOTE] convert narrowpeaks to 6 bed file for refining"
        cat ${n/.$ASD.bam/}_peaks.narrowPeak | cut -f1-6 > ${n/.$ASD.bam/}_peaks.bed
    fi

    RUN_COMMAND="macs2 refinepeak $MACS2_REFINEPEAK_ADDPARAM -b ${n/.$ASD.bam/}_peaks.bed -i $f --o-prefix ${n/.$ASD.bam/}  >> ${n/.$ASD.bam/}.summary.txt 2>&1"
    echo $RUN_COMMAND && eval $RUN_COMMAND
    [ -e ${n/.$ASD.bam/}_peaks.bed ] && rm ${n/.$ASD.bam/}_peaks.bed    

 
    # mark checkpoint
    if [ -f ${n/.$ASD.bam/}_refinepeak.bed ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
fi

################################################################################
cd $SOURCE

################################################################################
[ -e $OUTDIR/${n/.$ASD.bam/_refinepeak.bed}.dummy ] && rm $OUTDIR/${n/.$ASD.bam/_refinepeak.bed}.dummy
echo ">>>>> ChIPseq analysis with MACS2 - FINISHED"
echo ">>>>> enddate "`date`

