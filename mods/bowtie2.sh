#!/bin/bash -e

# Script to run bowtie program.
# It takes comma-seprated list of files containing short sequence reads in fasta or fastq format and bowtie index files as input.
# It produces output files: read alignments in .bam format and other files.
# author: Denis Bauer
# date: June 2012
# modified: August 2013 Fabian Buske

# QCVARIABLES,Resource temporarily unavailable
# RESULTFILENAME <SAMPLE>.$ASD.bam

echo ">>>>> readmapping with Bowtie2 "
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"

function usage {
echo -e "usage: $(basename $0) -k NGSANE -f FASTQ -r REFERENCE -o OUTDIR [OPTIONS]"
exit
}

if [ ! $# -gt 3 ]; then usage ; fi

FORCESINGLE=0

#INPUTS                                                                                                           
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the NGSANE repository                       
        -f | --fastq )          shift; f=$1 ;; # fastq file                                                       
        -r | --reference )      shift; FASTA=$1 ;; # reference genome                                             
        -o | --outdir )         shift; MYOUT=$1 ;; # output dir                                                     
        -i | --rgid )           shift; EXPID=$1 ;; # read group identifier RD ID                                  
        -l | --rglb )           shift; LIBRARY=$1 ;; # read group library RD LB                                   
        -p | --rgpl )           shift; PLATFORM=$1 ;; # read group platform RD PL                                 
        -s | --rgsi )           shift; SAMPLEID=$1 ;; # read group sample RG SM (pre)                             
        -u | --rgpu )           shift; UNIT=$1 ;; # read group platform unit RG PU
        --forceSingle )         FORCESINGLE=1;;
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

echo "********** programs"
for MODULE in $MODULE_BOWTIE2; do module load $MODULE; done  # save way to load modules that itself load other modules
export PATH=$PATH_BOWTIE2:$PATH
module list
echo "PATH=$PATH"
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)
PATH_IGVTOOLS=$(dirname $(which igvtools.jar))
PATH_PICARD=$(dirname $(which MarkDuplicates.jar))

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--JAVA        --\n" $(java -version 2>&1)
[ -z "$(which java)" ] && echo "[ERROR] no java detected" && exit 1
echo -e "--bowtie2     --\n "$(bowtie2 --version)
[ -z "$(which bowtie2)" ] && echo "[ERROR] no bowtie2 detected" && exit 1
echo -e "--samtools    --\n "$(samtools 2>&1 | head -n 3 | tail -n-2)
[ -z "$(which samtools)" ] && echo "[ERROR] no samtools detected" && exit 1
echo -e "--R           --\n "$(R --version | head -n 3)
[ -z "$(which R)" ] && echo "[ERROR] no R detected" && exit 1
echo -e "--igvtools    --\n "$(java -jar $JAVAPARAMS $PATH_IGVTOOLS/igvtools.jar version 2>&1)
[ ! -f $PATH_IGVTOOLS/igvtools.jar ] && echo "[ERROR] no igvtools detected" && exit 1
echo -e "--PICARD      --\n "$(java -jar $JAVAPARAMS $PATH_PICARD/MarkDuplicates.jar --version 2>&1)
[ ! -f $PATH_PICARD/MarkDuplicates.jar ] && echo "[ERROR] no picard detected" && exit 1
echo -e "--samstat     --\n "$(samstat -h | head -n 2 | tail -n1)
[ -z "$(which samstat)" ] && echo "[ERROR] no samstat detected" && exit 1
echo -e "--convert     --\n "$(convert -version | head -n 1)
[ -z "$(which convert)" ] && echo "[WARN] imagemagick convert not detected" && exit 1
echo -e "--wigToBigWig --\n "$(wigToBigWig 2>&1 | tee | head -n 1)
[ -z "$(which wigToBigWig)" ] && echo "[WARN] wigToBigWig not detected - no bigwigs will be generated"

echo "[NOTE] set java parameters"
JAVAPARAMS="-Xmx"$(python -c "print int($MEMORY_BOWTIE2*0.8)")"g -Djava.io.tmpdir="$TMP"  -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1" 
unset _JAVA_OPTIONS
echo "JAVAPARAMS "$JAVAPARAMS

echo -e "\n********* $CHECKPOINT"
################################################################################
CHECKPOINT="parameters"

# get basename of f
n=${f##*/}

# check library variables are set
if [[ -z "$EXPID" || -z "$LIBRARY" || -z "$PLATFORM" ]]; then
    echo "[ERROR] library info not set (EXPID, LIBRARY, and PLATFORM): free text needed"
    exit 1;
else
    echo "[NOTE] EXPID $EXPID; LIBRARY $LIBRARY; PLATFORM $PLATFORM"
fi


# delete old bam files unless attempting to recover
if [ -z "$RECOVERFROM" ]; then
    [ -e $MYOUT/${n/%$READONE.$FASTQ/.$ASD.bam} ] && rm $MYOUT/${n/%$READONE.$FASTQ/.$ASD.bam}
    [ -e $MYOUT/${n/%$READONE.$FASTQ/.$ASD.bam}.stats ] && rm $MYOUT/${n/%$READONE.$FASTQ/.$ASD.bam}.stats
    [ -e $MYOUT/${n/%$READONE.$FASTQ/.$ASD.bam}.dupl ] && rm $MYOUT/${n/%$READONE.$FASTQ/.$ASD.bam}.dupl
fi

#is ziped ?
ZCAT="zcat"
if [[ ${f##*.} != "gz" ]]; then ZCAT="cat"; fi

#is paired ?                                                                                                      
if [ "$f" != "${f/$READONE/$READTWO}" ] && [ -e ${f/$READONE/$READTWO} ] && [ "$FORCESINGLE" = 0 ]; then
    PAIRED="1"
    READS="-1 $f -2 ${f/$READONE/$READTWO}"
    READ1=`$ZCAT $f | wc -l | gawk '{print int($1/4)}' `
    READ2=`$ZCAT ${f/$READONE/$READTWO} | wc -l | gawk '{print int($1/4)}' `
    let FASTQREADS=$READ1+$READ2
else
    PAIRED="0"
    READS="-U $f"
    let FASTQREADS=`$ZCAT $f | wc -l | gawk '{print int($1/4)}' `
fi

# get encoding
FASTQ_ENCODING=$($ZCAT $f |  awk 'NR % 4 ==0' | python $NGSANE_BASE/tools/GuessFastqEncoding.py |  tail -n 1)
if [[ "$FASTQ_ENCODING" == *Phred33* ]]; then
    FASTQ_PHRED="--phred33"    
elif [[ "$FASTQ_ENCODING" == *Illumina* ]]; then
    FASTQ_PHRED="--phred64"
elif [[ "$FASTQ_ENCODING" == *Solexa* ]]; then
    FASTQ_PHRED="--solexa-quals"
else
    echo "[NOTE] cannot detect/don't understand fastq format: $FASTQ_ENCODING - using default"
fi
echo "[NOTE] $FASTQ_ENCODING fastq format detected"

#readgroup
FULLSAMPLEID=$SAMPLEID"${n/%$READONE.$FASTQ/}"
RG="--sam-rg \"ID:$EXPID\" --sam-rg \"SM:$FULLSAMPLEID\" --sam-rg \"LB:$LIBRARY\" --sam-rg \"PL:$PLATFORM\""


echo -e "\n********* $CHECKPOINT"
################################################################################
CHECKPOINT="recall files from tape"
	
if [ -n "$DMGET" ]; then
    dmget -a $FASTA*
	dmget -a ${f/$READONE/"*"}
fi
    
echo -e "\n********* $CHECKPOINT"
################################################################################
CHECKPOINT="generating the index files"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep "********* $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
    
    FASTASUFFIX=${FASTA##*.}
    if [ ! -e ${FASTA/.${FASTASUFFIX}/}.1.bt2 ]; then echo ">>>>> make .bt2"; bowtie2-build $FASTA ${FASTA/.${FASTASUFFIX}/}; fi
    if [ ! -e $FASTA.fai ]; then echo ">>>>> make .fai"; samtools faidx $FASTA; fi

    # mark checkpoint
    [ -f ${FASTA/.${FASTASUFFIX}/}.1.bt2 ] && echo -e "\n********* $CHECKPOINT" && unset RECOVERFROM
fi 

################################################################################
CHECKPOINT="run bowtie2"
if [[ -n "$RECOVERFROM" ]] && [[ $(grep "********* $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
    
    RUN_COMMAND="bowtie2 $RG $BOWTIE2ADDPARAM $FASTQ_PHRED -t -x ${FASTA/.${FASTASUFFIX}/} -p $CPU_BOWTIE2 $READS --un $MYOUT/${n/%$READONE.$FASTQ/.$ALN.un.sam} | samtools view -bS -t $FASTA.fai - > $MYOUT/${n/%$READONE.$FASTQ/.$ALN.bam}"
    echo $RUN_COMMAND && eval $RUN_COMMAND
    
    # mark checkpoint
    [ -f $MYOUT/${n/%$READONE.$FASTQ/.$ALN.bam} ] && echo -e "\n********* $CHECKPOINT" && unset RECOVERFROM
fi

################################################################################
CHECKPOINT="bam conversion and sorting"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep "********* $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
       
    samtools sort $MYOUT/${n/%$READONE.$FASTQ/.$ALN.bam} $MYOUT/${n/%$READONE.$FASTQ/.map}
    samtools view -bt $FASTA.fai $MYOUT/${n/%$READONE.$FASTQ/.$ALN.un.sam} | samtools sort - $MYOUT/${n/%$READONE.$FASTQ/.unm}
    # merge mappend and unmapped
    samtools merge -f $MYOUT/${n/%$READONE.$FASTQ/.ash}.bam $MYOUT/${n/%$READONE.$FASTQ/.map}.bam $MYOUT/${n/%$READONE.$FASTQ/.unm}.bam 
    # remove sam files 
    [ -e $MYOUT/${n/%$READONE.$FASTQ/.$ALN.bam} ] && rm $MYOUT/${n/%$READONE.$FASTQ/.$ALN.bam}
    [ -e $MYOUT/${n/%$READONE.$FASTQ/.$ALN.un.sam} ] && rm $MYOUT/${n/%$READONE.$FASTQ/.$ALN.un.sam}
    
    if [ "$PAIRED" = "1" ]; then
        # fix mates
        samtools sort -n $MYOUT/${n/%$READONE.$FASTQ/.ash}.bam $MYOUT/${n/%$READONE.$FASTQ/.ash}.bam.tmp
        samtools fixmate $MYOUT/${n/%$READONE.$FASTQ/.ash}.bam.tmp.bam - | samtools sort - $MYOUT/${n/%$READONE.$FASTQ/.ash}
        [ -e $MYOUT/${n/%$READONE.$FASTQ/.ash}.bam.tmp.bam ] && rm $MYOUT/${n/%$READONE.$FASTQ/.ash}.bam.tmp.bam
    fi
    
    # mark checkpoint
    [ -f $MYOUT/${n/%$READONE.$FASTQ/.ash}.bam ] && echo -e "\n********* $CHECKPOINT" && unset RECOVERFROM

fi

################################################################################
CHECKPOINT="mark duplicates"
# create bam files for discarded reads and remove fastq files
if [[ -n "$RECOVERFROM" ]] && [[ $(grep "********* $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
   
    if [ ! -e $MYOUT/metrices ]; then mkdir -p $MYOUT/metrices ; fi
    THISTMP=$TMP/$n$RANDOM #mk tmp dir because picard writes none-unique files                                        
    mkdir -p $THISTMP
    java $JAVAPARAMS -jar $PATH_PICARD/MarkDuplicates.jar \
        INPUT=$MYOUT/${n/%$READONE.$FASTQ/.ash.bam} \
        OUTPUT=$MYOUT/${n/%$READONE.$FASTQ/.$ASD.bam} \
        METRICS_FILE=$MYOUT/metrices/${n/%$READONE.$FASTQ/.$ASD.bam}.dupl AS=true \
        VALIDATION_STRINGENCY=LENIENT \
        TMP_DIR=$THISTMP
    [ -d $THISTMP ] && rm -r $THISTMP
    samtools index $MYOUT/${n/%$READONE.$FASTQ/.$ASD.bam}
    
    # mark checkpoint
    [ -f $MYOUT/${n/%$READONE.$FASTQ/.$ASD.bam} ] && echo -e "\n********* $CHECKPOINT" && unset RECOVERFROM
fi

################################################################################
CHECKPOINT="statistics"                                                                                                

if [[ -n "$RECOVERFROM" ]] && [[ $(grep "********* $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
    
    STATSOUT=$MYOUT/${n/%$READONE.$FASTQ/.$ASD.bam}.stats
    samtools flagstat $MYOUT/${n/%$READONE.$FASTQ/.$ASD.bam} > $STATSOUT
    if [ -n $SEQREG ]; then
        echo "#custom region" >> $STATSOUT
        echo `samtools view $MYOUT/${n/%$READONE.$FASTQ/.ash.bam} $SEQREG | wc -l`" total reads in region " >> $STATSOUT
        echo `samtools view -f 2 $MYOUT/${n/%$READONE.$FASTQ/.ash.bam} $SEQREG | wc -l`" properly paired reads in region " >> $STATSOUT
    fi

    # mark checkpoint
    [ -e $STATSOUT ] && echo -e "\n********* $CHECKPOINT" && unset RECOVERFROM
fi

################################################################################
CHECKPOINT="calculate inner distance"                                                                                                

if [[ -n "$RECOVERFROM" ]] && [[ $(grep "********* $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
    

    THISTMP=$TMP/$n$RANDOM #mk tmp dir because picard writes none-unique files
    mkdir $THISTMP
    java $JAVAPARAMS -jar $PATH_PICARD/CollectMultipleMetrics.jar \
        INPUT=$MYOUT/${n/%$READONE.$FASTQ/.$ASD.bam} \
        REFERENCE_SEQUENCE=$FASTA \
        OUTPUT=$MYOUT/metrices/${n/%$READONE.$FASTQ/.$ASD.bam} \
        VALIDATION_STRINGENCY=LENIENT \
        PROGRAM=CollectAlignmentSummaryMetrics \
        PROGRAM=CollectInsertSizeMetrics \
        PROGRAM=QualityScoreDistribution \
        TMP_DIR=$THISTMP
    for im in $( ls $MYOUT/metrices/*.pdf ); do
        convert $im ${im/pdf/jpg}
    done
    rm -r $THISTMP

    # mark checkpoint
    [ -f $MYOUT/metrices/${n/%$READONE.$FASTQ/.$ASD.bam}.alignment_summary_metrics ] && echo -e "\n********* $CHECKPOINT" && unset RECOVERFROM
fi

################################################################################
CHECKPOINT="coverage track"    

if [[ -n "$RECOVERFROM" ]] && [[ $(grep "********* $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

    
    java $JAVAPARAMS -jar $PATH_IGVTOOLS/igvtools.jar count $MYOUT/${n/%$READONE.$FASTQ/.$ASD.bam} $MYOUT/${n/%$READONE.$FASTQ/.$ASD.bam.cov.tdf} ${FASTA/.$FASTASUFFIX/}.genome
    # mark checkpoint
    [ -f $MYOUT/${n/%$READONE.$FASTQ/.$ASD.bam.cov.tdf} ] && echo -e "\n********* $CHECKPOINT" && unset RECOVERFROM
fi

################################################################################
CHECKPOINT="samstat"    

if [[ -n "$RECOVERFROM" ]] && [[ $(grep "********* $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
    
    samstat $MYOUT/${n/%$READONE.$FASTQ/.$ASD.bam}

    # mark checkpoint
    [ -f $MYOUT/${n/%$READONE.$FASTQ/.$ASD.bam}.stats ] && echo -e "\n********* $CHECKPOINT" && unset RECOVERFROM
    
fi

###############################################################################
CHECKPOINT="verify"    

BAMREADS=`head -n1 $MYOUT/${n/%$READONE.$FASTQ/.$ASD.bam}.stats | cut -d " " -f 1`
if [ "$BAMREADS" = "" ]; then let BAMREADS="0"; fi
if [ $BAMREADS -eq $FASTQREADS ]; then
    echo "-----------------> PASS check mapping: $BAMREADS == $FASTQREADS"
    [ -e $MYOUT/${n/%$READONE.$FASTQ/.ash.bam} ] && rm $MYOUT/${n/%$READONE.$FASTQ/.ash.bam}
    [ -e $MYOUT/${n/%$READONE.$FASTQ/.unm}.bam ] && rm $MYOUT/${n/%$READONE.$FASTQ/.unm}.bam
    [ -e $MYOUT/${n/%$READONE.$FASTQ/.map}.bam ] && rm $MYOUT/${n/%$READONE.$FASTQ/.map}.bam
else
    echo -e "[ERROR] We are loosing reads from .fastq -> .bam in $f: \nFastq had $FASTQREADS Bam has $BAMREADS"
    exit 1
fi

echo "********* $CHECKPOINT"
################################################################################
CHECKPOINT="generate bigwigs"    

FRAGMENTLENGTH=0
GENOME_CHROMSIZES=$FASTA.chrom.size
. $CONFIG # overwrite defaults

echo "$FASTQ $MYOUT/${n/%$READONE.$FASTQ/.$ASD.bam}"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep "********* $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
        
    if [ -z "$(which wigToBigWig)" ]; then
        echo "[NOTE] Skip bigwig generation due to missing software: wigToBigWig"
        
    else
        NORMALIZETO=1000000
        NUMBEROFREADS=$(samtools view -c -F 1028 $MYOUT/${n/%$READONE.$FASTQ/.$ASD.bam})
        SCALEFACTOR=`echo "scale=3; $NORMALIZETO/$NUMBEROFREADS" | bc`
      
        echo "[NOTE] library size (mapped reads): $NORMALIZETO" 
        echo "[NOTE] scale factor: $SCALEFACTOR"
        echo "[NOTE] fragment length: $FRAGMENTLENGTH"
    
        if [ "$PAIRED" = "1" ] && [[ $FRAGMENTLENGTH -le 0 ]]; then
            echo "[NOTE] generate bigwig for properly paired reads on the same chromosomes"
            samtools sort -n $MYOUT/${n/%$READONE.$FASTQ/.$ASD.bam} $MYOUT/${n/%$READONE.$FASTQ/.$ASD.tmp}
            samtools view -b -F 1028 -f 0x2 $MYOUT/${n/%$READONE.$FASTQ/.$ASD.tmp.bam} | bamToBed -bedpe | awk '($1 == $4){OFS="\t"; print $1,$2,$6,$7,$8,$9}' | genomeCoverageBed -scale $SCALEFACTOR -g ${GENOME_CHROMSIZES} -i stdin -bg | wigToBigWig stdin  ${GENOME_CHROMSIZES} $MYOUT/${n/%$READONE.$FASTQ/.bw}
            [ -e $MYOUT/${n/%$READONE.$FASTQ/.$ASD.tmp.bam} ] && rm $MYOUT/${n/%$READONE.$FASTQ/.$ASD.tmp.bam}
    	
        else
            
            if [[ $FRAGMENTLENGTH -ge 0 ]]; then
            	echo "[NOTE] Skip bigwig generation due to invalid value for parameter: FRAGMENTLENGTH"

            else

                if [ "$BIGWIGSTRANDS" = "strand-specific" ]; then 
                    echo "[NOTE] generate strand-specific bigwigs considering single reads"
                    samtools view -b -F 1028 $MYOUT/${n/%$READONE.$FASTQ/.$ASD.bam} | bamToBed | slopBed -s -r $FRAGMENTLENGTH -l 0 -i stdin -g ${GENOME_CHROMSIZES}  | genomeCoverageBed -strand "+" -scale $SCALEFACTOR -g ${GENOME_CHROMSIZES} -i stdin -bg | wigToBigWig stdin  ${GENOME_CHROMSIZES} $MYOUT/${n/%$READONE.$FASTQ/.+.bw}
                    
                    samtools view -b -F 1028 $MYOUT/${n/%$READONE.$FASTQ/.$ASD.bam} | bamToBed | slopBed -s -r $FRAGMENTLENGTH -l 0 -i stdin -g ${GENOME_CHROMSIZES}  | genomeCoverageBed -strand "-" -scale $SCALEFACTOR -g ${GENOME_CHROMSIZES} -i stdin -bg | wigToBigWig stdin  ${GENOME_CHROMSIZES} $MYOUT/${n/%$READONE.$FASTQ/.-.bw}
        
                else
                    echo "[NOTE] generate bigwig considering single reads"
                    samtools view -b -F 1028 $MYOUT/${n/%$READONE.$FASTQ/.$ASD.bam} | bamToBed | slopBed -s -r $FRAGMENTLENGTH -l 0 -i stdin -g ${GENOME_CHROMSIZES}  | genomeCoverageBed -scale $SCALEFACTOR -g ${GENOME_CHROMSIZES} -i stdin -bg | wigToBigWig stdin  ${GENOME_CHROMSIZES} $MYOUT/${n/%$READONE.$FASTQ/.bw}
                    
            	fi
        	fi
        fi
    
        # mark checkpoint
        [ -f $MYOUT/${n/%$READONE.$FASTQ/.$ASD.bam}.stats ] && echo -e "\n********* $CHECKPOINT" && unset RECOVERFROM

    fi   
fi

################################################################################
[ -e $MYOUT/${n/%$READONE.$FASTQ/.$ASD.bam}.dummy ] && rm $MYOUT/${n/%$READONE.$FASTQ/.$ASD.bam}.dummy
echo ">>>>> readmapping with Bowtie2 - FINISHED"
echo ">>>>> enddate "`date`
