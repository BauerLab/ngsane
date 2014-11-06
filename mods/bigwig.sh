#!/bin/bash -e

# Script to generate bigwigs from bam files.
# author: Fabian Buske
# date: August 2013

# QCVARIABLES,Resource temporarily unavailable
# RESULTFILENAME <DIR>/$TASK_BIGWIG/<SAMPLE>.bw

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
hash module 2>/dev/null && for MODULE in $MODULE_BIGWIG; do module load $MODULE; done && module list 

export PATH=$PATH_BIGWIG:$PATH
echo "PATH=$PATH"
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--bedtools --\n "$(bedtools --version)
[ -z "$(which bedtools)" ] && echo "[ERROR] no bedtools detected" && exit 1
echo -e "--wigToBigWig --\n "$(wigToBigWig 2>&1 | tee | head -n 1)
[ -z "$(which wigToBigWig)" ] && echo "[ERROR] wigToBigWig not detected" && exit 1


NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "parameters"

[ ! -f $f ] && echo "[ERROR] input file not found: $f" && exit 1

# get basename of f
n=${f##*/}
SAMPLE=${n/%$ASD.bam/}

# delete old bw files unless attempting to recover
if [ -z "$NGSANE_RECOVERFROM" ]; then
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

DUPLICATEFILTERFLAG=1028  # filter unmapped and duplicate reads
if [ -n "$KEEPDUPLICATES" ]; then
    DUPLICATEFILTERFLAG=4 # filter unmapped reads only
fi

THISTMP=$TMP"/"$(whoami)"/"$(echo $OUTDIR/$SAMPLE | md5sum | cut -d' ' -f1)
[ -d $THISTMP ] && rm -r $THISTMP
mkdir -p $THISTMP

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "recall files from tape"

if [ -n "$DMGET" ]; then
	dmget -a $(dirname $FASTA)/*
	dmget -a ${f}*
	dmget -a ${OUTDIR}
fi
    
NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "generate bigwigs"    

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then
   
    # run flagstat if no stats available for bam file
    [ ! -e $f.stats ] && samtools flagstat > $f.stats
    # check "paired in sequencing" entry to detect library
    if [[ $(cat $f.stats | fgrep -w 'paired in sequencing' | cut -d' ' -f1) -gt 0 ]]; then
        PAIRED=1
        echo "[NOTE] paired library detected"
    else 
        PAIRED=0
        echo "[NOTE] single-end library detected"
    fi

    NUMBEROFREADS=$(samtools view -c -F $DUPLICATEFILTERFLAG $f )
	if [[ -n "$NORMALIZETO" ]] && [[ -z "$SCALEFACTOR" ]]; then 
		SCALEFACTOR=`echo "scale=3; $NORMALIZETO/$NUMBEROFREADS" | bc`; 
	else
    	SCALEFACTOR="1"
		NORMALIZETO="NA"
	fi
    
    echo "library size: $NUMBEROFREADS" > $OUTDIR/$SAMPLE.bw.stats
    echo "normalize to: $NORMALIZETO" >> $OUTDIR/$SAMPLE.bw.stats
    echo "scale factor: $SCALEFACTOR" >> $OUTDIR/$SAMPLE.bw.stats
    
        
    if [ "$PAIRED" = "1" ]; then 
    	echo "[NOTE] Paired end libaray detected, ignore fragment length parameter"

        samtools sort -@ $CPU_BIGWIG -n $f $THISTMP/$SAMPLE.tmp
        echo "[NOTE] generate bigwig for properly paired reads on the same chromosomes"
        
        if [ -z "$FRAGMENTMIDPOINT" ]; then    
 
            echo "[NOTE] signal spans readpair fragment"
            samtools view -@ $CPU_BIGWIG -b -F $DUPLICATEFILTERFLAG -f 0x2 $THISTMP/$SAMPLE.tmp.bam | bamToBed -bedpe | awk '($1 == $4){OFS="\t"; print $1,$2,$6,$7,$8,$9}' | sort -k1,1 -k2,3g | genomeCoverageBed -scale $SCALEFACTOR -g ${GENOME_CHROMSIZES} -i stdin -bg | wigToBigWig stdin  ${GENOME_CHROMSIZES} $OUTDIR/$SAMPLE.bw
            
            if [ "$BIGWIGSTRANDS" = "strand-specific" ]; then 
                echo "[NOTE] generate strand-specific bigwigs too"
                samtools view -@ $CPU_BIGWIG -b -F $DUPLICATEFILTERFLAG -f 0x2 $THISTMP/$SAMPLE.tmp.bam | bamToBed -bedpe | awk '($1 == $4){OFS="\t"; print $1,$2,$6,$7,$8,$9}' | sort -k1,1 -k2,3g | genomeCoverageBed -strand "+" -scale $SCALEFACTOR -g ${GENOME_CHROMSIZES} -i stdin -bg | wigToBigWig stdin  ${GENOME_CHROMSIZES} $OUTDIR/$SAMPLE.watson.bw
                
                samtools view -@ $CPU_BIGWIG -b -F $DUPLICATEFILTERFLAG -f 0x2 $THISTMP/$SAMPLE.tmp.bam | bamToBed -bedpe | awk '($1 == $4){OFS="\t"; print $1,$2,$6,$7,$8,$9}' | sort -k1,1 -k2,3g | genomeCoverageBed -strand "-" -scale -$SCALEFACTOR -g ${GENOME_CHROMSIZES} -i stdin -bg | wigToBigWig stdin  ${GENOME_CHROMSIZES} $OUTDIR/$SAMPLE.crick.bw
        	fi
        	
        else
            echo "[NOTE] considered midpoint of fragment only"
            samtools view -@ $CPU_BIGWIG -b -F $DUPLICATEFILTERFLAG -f 0x2 $THISTMP/$SAMPLE.tmp.bam | bamToBed -bedpe | awk '($1 == $4){OFS="\t"; print $1,int(($6+$2)/2),int(($6+$2)/2)+1,$7,$8,$9}' | sort -k1,1 -k2,3g | genomeCoverageBed -scale $SCALEFACTOR -g ${GENOME_CHROMSIZES} -i stdin -bg | wigToBigWig stdin  ${GENOME_CHROMSIZES} $OUTDIR/$SAMPLE.bw
            
            if [ "$BIGWIGSTRANDS" = "strand-specific" ]; then 
                echo "[NOTE] generate strand-specific bigwigs too"
                samtools view -@ $CPU_BIGWIG -b -F $DUPLICATEFILTERFLAG -f 0x2 $THISTMP/$SAMPLE.tmp.bam | bamToBed -bedpe | awk '($1 == $4){OFS="\t"; print $1,int(($6+$2)/2),int(($6+$2)/2)+1,$7,$8,$9}' | sort -k1,1 -k2,3g | genomeCoverageBed -strand "+" -scale $SCALEFACTOR -g ${GENOME_CHROMSIZES} -i stdin -bg | wigToBigWig stdin  ${GENOME_CHROMSIZES} $OUTDIR/$SAMPLE.watson.bw
                
                samtools view -@ $CPU_BIGWIG -b -F $DUPLICATEFILTERFLAG -f 0x2 $THISTMP/$SAMPLE.tmp.bam | bamToBed -bedpe | awk '($1 == $4){OFS="\t"; print $1,int(($6+$2)/2),int(($6+$2)/2)+1,$7,$8,$9}' | sort -k1,1 -k2,3g | genomeCoverageBed -strand "-" -scale -$SCALEFACTOR -g ${GENOME_CHROMSIZES} -i stdin -bg | wigToBigWig stdin  ${GENOME_CHROMSIZES} $OUTDIR/$SAMPLE.crick.bw
        	fi
        fi
        
        if [ "$BIGWIGSTRANDS" = "strand-specific" ]; then 
            echo "[NOTE] generate strand-specific bigwigs too"
            samtools view -@ $CPU_BIGWIG -b -F $DUPLICATEFILTERFLAG -f 0x2 $THISTMP/$SAMPLE.tmp.bam | bamToBed -bedpe | awk '($1 == $4){OFS="\t"; print $1,$2,$6,$7,$8,$9}' | sort -k1,1 -k2,3g | genomeCoverageBed -strand "+" -scale $SCALEFACTOR -g ${GENOME_CHROMSIZES} -i stdin -bg | wigToBigWig stdin  ${GENOME_CHROMSIZES} $OUTDIR/$SAMPLE.watson.bw
            
            samtools view -@ $CPU_BIGWIG -b -F $DUPLICATEFILTERFLAG -f 0x2 $THISTMP/$SAMPLE.tmp.bam | bamToBed -bedpe | awk '($1 == $4){OFS="\t"; print $1,$2,$6,$7,$8,$9}' | sort -k1,1 -k2,3g | genomeCoverageBed -strand "-" -scale -$SCALEFACTOR -g ${GENOME_CHROMSIZES} -i stdin -bg | wigToBigWig stdin  ${GENOME_CHROMSIZES} $OUTDIR/$SAMPLE.crick.bw
    	fi
            
        [ -e $THISTMP/$SAMPLE.tmp.bam ] && rm $THISTMP/$SAMPLE.tmp.bam

    else # single end libraries
        echo "[NOTE] fragment length: $FRAGMENTLENGTH"

        if [[ $FRAGMENTLENGTH -lt 0 ]]; then
    	   # leave fragments as is
            echo "[NOTE] generate bigwig from original reads."

            samtools view -@ $CPU_BIGWIG -b -F $DUPLICATEFILTERFLAG $f | bamToBed | genomeCoverageBed -scale $SCALEFACTOR -g ${GENOME_CHROMSIZES} -i stdin -bg | wigToBigWig stdin  ${GENOME_CHROMSIZES} $OUTDIR/$SAMPLE.bw

            if [ "$BIGWIGSTRANDS" = "strand-specific" ]; then 
                echo "[NOTE] generate strand-specific bigwigs too"
                samtools view -@ $CPU_BIGWIG -b -F $DUPLICATEFILTERFLAG $f | bamToBed | genomeCoverageBed -strand "+" -scale $SCALEFACTOR -g ${GENOME_CHROMSIZES} -i stdin -bg | wigToBigWig stdin  ${GENOME_CHROMSIZES} $OUTDIR/$SAMPLE.watson.bw
                
                samtools view -@ $CPU_BIGWIG -b -F $DUPLICATEFILTERFLAG $f | bamToBed | genomeCoverageBed -strand "-" -scale $SCALEFACTOR -g ${GENOME_CHROMSIZES} -i stdin -bg | wigToBigWig stdin  ${GENOME_CHROMSIZES} $OUTDIR/$SAMPLE.crick.bw   
            fi
            	
        else
            # adjust fragment length to given size
            echo "[NOTE] generate bigwig with specified fragment length"
        
            samtools view -@ $CPU_BIGWIG -b -F $DUPLICATEFILTERFLAG $f | bamToBed | awk 'BEGIN{OFS="\t"}{if($6=="+"){print $1,$2,$2+1,$4,$5,$6}else if($6=="-"){print $1,$3-1,$3,$4,$5,$6}else{print $0}}' | slopBed -s -r $FRAGMENTLENGTH -l 0 -i stdin -g ${GENOME_CHROMSIZES}  | genomeCoverageBed -scale $SCALEFACTOR -g ${GENOME_CHROMSIZES} -i stdin -bg | wigToBigWig stdin ${GENOME_CHROMSIZES} $OUTDIR/$SAMPLE.bw

            if [ "$BIGWIGSTRANDS" = "strand-specific" ]; then 
                echo "[NOTE] generate strand-specific bigwigs too"
                samtools view -@ $CPU_BIGWIG -b -F $DUPLICATEFILTERFLAG $f | bamToBed | awk 'BEGIN{OFS="\t";}{if($6=="+"){print $1,$2,$2+1,$4,$5,$6}}' | slopBed -s -r $FRAGMENTLENGTH -l 0 -i stdin -g ${GENOME_CHROMSIZES}  | genomeCoverageBed -strand "+" -scale $SCALEFACTOR -g ${GENOME_CHROMSIZES} -i stdin -bg | wigToBigWig stdin  ${GENOME_CHROMSIZES} $OUTDIR/$SAMPLE.watson.bw
                
                samtools view -@ $CPU_BIGWIG -b -F $DUPLICATEFILTERFLAG $f | bamToBed | awk 'BEGIN{OFS="\t";}{if($6=="-"){print $1,$3-1,$3,$4,$5,$6}}' | slopBed -s -r $FRAGMENTLENGTH -l 0 -i stdin -g ${GENOME_CHROMSIZES}  | genomeCoverageBed -strand "-" -scale -$SCALEFACTOR -g ${GENOME_CHROMSIZES} -i stdin -bg | wigToBigWig stdin  ${GENOME_CHROMSIZES} $OUTDIR/$SAMPLE.crick.bw
        	fi
    	fi
    fi 
      
    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE.bw 
 
fi

################################################################################
[ -e $OUTDIR/$SAMPLE.bw.dummy ] && rm $OUTDIR/$SAMPLE.bw.dummy
echo ">>>>> bigwig generation - FINISHED"
echo ">>>>> enddate "`date`
