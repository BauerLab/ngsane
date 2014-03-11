#!/bin/bash -e

# Script to generate bigwigs from bam files.
# author: Fabian Buske
# date: August 2013

# QCVARIABLES,Resource temporarily unavailable
# RESULTFILENAME <DIR>/$INPUT_BIGWIG/<SAMPLE>.bw

echo ">>>>> bigwig generation"
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"

function usage {
echo -e "usage: $(basename $0) -k NGSANE -f bam -o OUTDIR [OPTIONS]"
exit
}

if [ ! $# -gt 3 ]; then usage ; fi

FORCESINGLE=0

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

for MODULE in $MODULE_BIGWIG; do module load $MODULE; done  # save way to load modules that itself load other modules
export PATH=$PATH_BIGWIG:$PATH
module list
echo "PATH=$PATH"
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--bedtools --\n "$(bedtools --version)
[ -z "$(which bedtools)" ] && echo "[ERROR] no bedtools detected" && exit 1
echo -e "--wigToBigWig --\n "$(wigToBigWig 2>&1 | tee | head -n 1)
[ -z "$(which wigToBigWig)" ] && echo "[ERROR] wigToBigWig not detected" && exit 1


echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="parameters"

[ ! -f $f ] && echo "[ERROR] input file not found: $f" && exit 1

# get basename of f
n=${f##*/}
SAMPLE=${n/%.$ASD.bam/}

# delete old bw files unless attempting to recover
if [ -z "$RECOVERFROM" ]; then
    [ -e $OUTDIR/$SAMPLE.bw ] && rm $OUTDIR/$SAMPLE.bw
fi

if [ -z "$FASTA" ] || [ ! -f $FASTA ]; then
    echo "[ERROR] no reference provided (FASTA)"
    exit 1
else
    echo "[NOTE] Reference: $FASTA"
fi

GENOME_CHROMSIZES=${FASTA%.*}.chrom.sizes
if [ ! -f $GENOME_CHROMSIZES ]; then
    echo "[ERROR] GENOME_CHROMSIZES not found. Expected at $GENOME_CHROMSIZES"
    exit 1
else
    echo "[NOTE] Chromosome size: $GENOME_CHROMSIZES"

fi

if [ -z "$FRAGMENTLENGTH" ]; then
    FRAGMENTLENGTH=0
fi

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="recall files from tape"

if [ -n "$DMGET" ]; then
	dmget -a $(dirname $FASTA)/*
	dmget -a ${f}*
	dmget -a ${OUTDIR}
fi
    
echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="generate bigwigs"    

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
   
    # run flagstat if no stats available for bam file
    [ ! -e $f.stats ] && samtools flagstat > $f.stats
    # check "paired in sequencing" entry to detect library
    if [[ $(cat $f.stats | head -n 4 | tail -n 1 | cut -d' ' -f 1) -gt 0 ]]; then
        PAIRED=1
        echo "[NOTE] paired library detected"
    else 
        PAIRED=0
        echo "[NOTE] single-end library detected"
    fi

	if [ -z "$NORMALIZETO" ]; then NORMALIZETO=1000000; fi
    NUMBEROFREADS=$(samtools view -c -F 1028 $f )
	if [ -z "$SCALEFACTOR" ]; then 
		SCALEFACTOR=`echo "scale=3; $NORMALIZETO/$NUMBEROFREADS" | bc`; 
	else
		NORMALIZETO="NA"
	fi
    
    echo "library size: $NUMBEROFREADS" > $OUTDIR/$SAMPLE.bw.stats
    echo "normalize to: $NORMALIZETO" >> $OUTDIR/$SAMPLE.bw.stats
    echo "scale factor: $SCALEFACTOR" >> $OUTDIR/$SAMPLE.bw.stats
        
    if [ "$PAIRED" = "1" ]; then 
    	echo "[NOTE] Paired end libaray detected, ignore fragment length parameter"

        samtools sort -@ $CPU_BIGWIG -n $f $OUTDIR/$SAMPLE.tmp
        echo "[NOTE] generate bigwig for properly paired reads on the same chromosomes"
            
        samtools view -@ $CPU_BIGWIG -b -F 1028 -f 0x2 $OUTDIR/$SAMPLE.tmp.bam | bamToBed -bedpe | awk '($1 == $4){OFS="\t"; print $1,$2,$6,$7,$8,$9}' | sort -k1,1 -k2,3g | genomeCoverageBed -scale $SCALEFACTOR -g ${GENOME_CHROMSIZES} -i stdin -bg | wigToBigWig stdin  ${GENOME_CHROMSIZES} $OUTDIR/$SAMPLE.bw

            if [ "$BIGWIGSTRANDS" = "strand-specific" ]; then 
                echo "[NOTE] generate strand-specific bigwigs too"
                samtools view -@ $CPU_BIGWIG -b -F 1028 -f 0x2 $OUTDIR/$SAMPLE.tmp.bam | bamToBed -bedpe | awk '($1 == $4){OFS="\t"; print $1,$2,$6,$7,$8,$9}' | sort -k1,1 -k2,3g | genomeCoverageBed -strand "+" -scale $SCALEFACTOR -g ${GENOME_CHROMSIZES} -i stdin -bg | wigToBigWig stdin  ${GENOME_CHROMSIZES} $OUTDIR/$SAMPLE.+.bw
                
                samtools view -@ $CPU_BIGWIG -b -F 1028 -f 0x2 $OUTDIR/$SAMPLE.tmp.bam | bamToBed -bedpe | awk '($1 == $4){OFS="\t"; print $1,$2,$6,$7,$8,$9}' | sort -k1,1 -k2,3g | genomeCoverageBed -strand "-" -scale -$SCALEFACTOR -g ${GENOME_CHROMSIZES} -i stdin -bg | wigToBigWig stdin  ${GENOME_CHROMSIZES} $OUTDIR/$SAMPLE.-.bw
        	fi
            
        [ -e $OUTDIR/$SAMPLE.tmp.bam ] && rm $OUTDIR/$SAMPLE.tmp.bam

    else # single end libraries
        echo "[NOTE] fragment length: $FRAGMENTLENGTH"

        if [[ $FRAGMENTLENGTH -lt 0 ]]; then
    	   # leave fragments as is
            echo "[NOTE] generate bigwig from original reads."

            samtools view -@ $CPU_BIGWIG -b -F 1028 $f | bamToBed | genomeCoverageBed -scale $SCALEFACTOR -g ${GENOME_CHROMSIZES} -i stdin -bg | wigToBigWig stdin  ${GENOME_CHROMSIZES} $OUTDIR/$SAMPLE.bw

            if [ "$BIGWIGSTRANDS" = "strand-specific" ]; then 
                echo "[NOTE] generate strand-specific bigwigs too"
                samtools view -@ $CPU_BIGWIG -b -F 1028 $f | bamToBed | genomeCoverageBed -strand "+" -scale $SCALEFACTOR -g ${GENOME_CHROMSIZES} -i stdin -bg | wigToBigWig stdin  ${GENOME_CHROMSIZES} $OUTDIR/$SAMPLE.+.bw
                
                samtools view -@ $CPU_BIGWIG -b -F 1028 $f | bamToBed | genomeCoverageBed -strand "-" -scale $SCALEFACTOR -g ${GENOME_CHROMSIZES} -i stdin -bg | wigToBigWig stdin  ${GENOME_CHROMSIZES} $OUTDIR/$SAMPLE.-.bw   
            fi
            	
        else
            # adjust fragment length to given size
            echo "[NOTE] generate bigwig with specified fragment length"
        
            samtools view -@ $CPU_BIGWIG -b -F 1028 $f | bamToBed | awk 'BEGIN{OFS="\t"}{if($6=="+"){print $1,$2,$2+1,$4,$5,$6}else if($6=="-"){print $1,$3-1,$3,$4,$5,$6}else{print $0}}' | slopBed -s -r $FRAGMENTLENGTH -l 0 -i stdin -g ${GENOME_CHROMSIZES}  | genomeCoverageBed -scale $SCALEFACTOR -g ${GENOME_CHROMSIZES} -i stdin -bg | wigToBigWig stdin ${GENOME_CHROMSIZES} $OUTDIR/$SAMPLE.bw

            if [ "$BIGWIGSTRANDS" = "strand-specific" ]; then 
                echo "[NOTE] generate strand-specific bigwigs too"
                samtools view -@ $CPU_BIGWIG -b -F 1028 $f | bamToBed | awk 'BEGIN{OFS="\t";}{if($6=="+"){print $1,$2,$2+1,$4,$5,$6}}' | slopBed -s -r $FRAGMENTLENGTH -l 0 -i stdin -g ${GENOME_CHROMSIZES}  | genomeCoverageBed -strand "+" -scale $SCALEFACTOR -g ${GENOME_CHROMSIZES} -i stdin -bg | wigToBigWig stdin  ${GENOME_CHROMSIZES} $OUTDIR/$SAMPLE.+.bw
                
                samtools view -@ $CPU_BIGWIG -b -F 1028 $f | bamToBed | awk 'BEGIN{OFS="\t";}{if($6=="-"){print $1,$3-1,$3,$4,$5,$6}}' | slopBed -s -r $FRAGMENTLENGTH -l 0 -i stdin -g ${GENOME_CHROMSIZES}  | genomeCoverageBed -strand "-" -scale -$SCALEFACTOR -g ${GENOME_CHROMSIZES} -i stdin -bg | wigToBigWig stdin  ${GENOME_CHROMSIZES} $OUTDIR/$SAMPLE.-.bw
        	fi
    	fi
    fi 
      
    # mark checkpoint
    if [ -f $OUTDIR/$SAMPLE.bw ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
 
fi

################################################################################
[ -e $OUTDIR/$SAMPLE.bw.dummy ] && rm $OUTDIR/$SAMPLE.bw.dummy
echo ">>>>> bigwig generation - FINISHED"
echo ">>>>> enddate "`date`
