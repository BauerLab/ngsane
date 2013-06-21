#!/bin/bash

echo ">>>>> HiC analysis with homer"
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> hicHomer.sh $*"

function usage {
echo -e "usage: $(basename $0) -k NGSANE -f FASTQ -r REFERENCE -e ENZYMES -o OUTDIR [OPTIONS]

Script running HICHOMER pipeline tapping into bowtie2
It expects a fastq file, paired end, reference genome and digest pattern  as input.

required:
  -k | --toolkit <path>     location of the NGSANE repository 
  -f | --fastq <file>       fastq file
  -e | --enzymes <name>     restriction enzyme (one per library) seperated by comma 
                            see http://biopython.org/DIST/docs/api/Bio.Restriction-module.html
  -o | --outdir <path>      output dir

options:
  -t | --threads <nr>       number of CPUs to use (default: 8)
  --fastqName               name of fastq file ending (fastq.gz)
  --oldIllumina
"
exit
}

# Script to run a hic analysis based on the HICHOMER framework.
# It takes comma-seprated list of files containing short sequence reads in fasta or fastq format and bowtie index files as input.
# It produces output files: read alignments in .bam format and other files.
# author: Fabian Buske
# date: April 2013

# QCVARIABLES,Resource temporarily unavailable

if [ ! $# -gt 3 ]; then usage ; fi

#DEFAULTS
THREADS=8
EXPID="exp"           # read group identifier RD ID
LIBRARY="tkcc"        # read group library RD LB
PLATFORM="illumina"   # read group platform RD PL
UNIT="flowcell"       # read group platform unit RG PU
FASTQNAME=""
ENZYME=""
QUAL="" # standard Sanger

#INPUTS                                                                                                           
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the NGSANE repository                       
        -t | --threads )        shift; THREADS=$1 ;; # number of CPUs to use                                      
        -f | --fastq )          shift; f=$1 ;; # fastq file                                                       
        -e | --enzymes )        shift; ENZYME=$1 ;; # digestion patterns
        -o | --outdir )         shift; MYOUT=$1 ;; # output dir                                                     
        -i | --rgid )           shift; EXPID=$1 ;; # read group identifier RD ID                                  
        -l | --rglb )           shift; LIBRARY=$1 ;; # read group library RD LB                                   
        -p | --rgpl )           shift; PLATFORM=$1 ;; # read group platform RD PL                                 
        -s | --rgsi )           shift; SAMPLEID=$1 ;; # read group sample RG SM (pre)                             
        -u | --rgpu )           shift; UNIT=$1 ;; # read group platform unit RG PU
        --fastqName )           shift; FASTQNAME=$1 ;; #(name of fastq or fastq.gz)
        -h | --help )           usage ;;
        * )                     echo "don't understand "$1
    esac
    shift
done

if [ -z "$ENZYME" ]; then
	echo "[ERROR] restriction enzyme not specified"
	exit 1
fi

if [ -z "BOWTIE2INDEX" ]; then
	echo "[ERROR] bowtie index not specified"
	exit 1
fi

#PROGRAMS
. $CONFIG
. ${NGSANE_BASE}/conf/header.sh
. $CONFIG

echo "********** programs"
for MODULE in $MODULE_HICHOMER; do module load $MODULE; done  # save way to load modules that itself load other modules

export PATH=$PATH_HICHOMER:$PATH
module list
echo "PATH=$PATH"
echo -e "--bowtie2 --\n "$(bowtie2 --version)
[ -z "$(which bowtie2)" ] && echo "[ERROR] no bowtie2 detected" && exit 1
echo -e "--samtools--\n "$(samtools 2>&1 | head -n 3 | tail -n-2)
[ -z "$(which samtools)" ] && echo "[ERROR] no samtools detected" && exit 1
echo -e "--Python      --\n" $(python --version)
[ -z "$(which python)" ] && echo "[ERROR] no python detected" && exit 1
echo -e "--Python libs --\n "$(yolk -l)

# get basename of f
n=${f##*/}

# delete old bam file                                                                                             
if [ -e $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam} ]; then rm $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}; fi
if [ -e $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}.stats ]; then rm $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}.stats; fi
if [ -e $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}.dupl ]; then rm $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}.dupl; fi

#is paired ?                                                                                                      
if [ -e ${f/$READONE/$READTWO} ] && [ "$FORCESINGLE" = 0 ]; then
    PAIRED="1"
else
    PAIRED="0"
fi

#is ziped ?                                                                                                       
ZCAT="zcat"
if [[ ${f##*.} != "gz" ]]; then ZCAT="cat"; fi

echo "********* generating the index files"
FASTASUFFIX=${FASTA##*.}
if [ ! -e ${FASTA/.${FASTASUFFIX}/}.1.bt2 ]; then echo ">>>>> make .bt2"; bowtie2-build $FASTA ${FASTA/.${FASTASUFFIX}/}; fi
if [ ! -e $FASTA.fai ]; then echo ">>>>> make .fai"; samtools faidx $FASTA; fi

if [ -n "$DMGET" ]; then
	echo "********** reacall files from tape"
	dmget -a $(dirname $FASTA)/*
	dmls -l $FASTA*
	dmget -a ${f/$READONE/"*"}
	dmls -l ${f/$READONE/"*"}
fi

echo "********* bowtie" 
if [ $PAIRED == "0" ]; then 
    echo "[ERROR] paired library required for HIC analysis"
    exit 1
else 
    READ1=`$ZCAT $f | wc -l | gawk '{print int($1/4)}' `
    READ2=`$ZCAT ${f/$READONE/$READTWO} | wc -l | gawk '{print int($1/4)}' `
    let FASTQREADS=$READ1+$READ2
fi

#readgroup
FULLSAMPLEID=$SAMPLEID"${n/'_'$READONE.$FASTQ/}"
RG="--sam-rg \"ID:$EXPID\" --sam-rg \"SM:$FULLSAMPLEID\" --sam-rg \"LB:$LIBRARY\" --sam-rg \"PL:$PLATFORM\""

for READLIB in $READONE $READTWO; do
    
    RUN_COMMAND="bowtie2 $RG -t -x ${FASTA/.${FASTASUFFIX}/} -p $THREADS  ${f/$READONE/$READLIB} -S $MYOUT/${n/'_'$READONE.$FASTQ/'_'$READLIB.$ALN.sam} --un $MYOUT/${n/'_'$READONE.$FASTQ/'_'$READLIB.$ALN.un.sam}"
    echo $RUN_COMMAND
    eval $RUN_COMMAND
    
    # continue for normal bam file conversion                                                                         
    samtools view -bt $FASTA.fai $MYOUT/${n/'_'$READONE.$FASTQ/'_'$READLIB.$ALN.sam} | samtools sort - $MYOUT/${n/'_'$READONE.$FASTQ/'_'$READLIB.map}
    samtools view -bt $FASTA.fai $MYOUT/${n/'_'$READONE.$FASTQ/'_'$READLIB.$ALN.un.sam} | samtools sort - $MYOUT/${n/'_'$READONE.$FASTQ/'_'$READLIB.unm}
    
    # merge mappend and unmapped
    samtools merge -f $MYOUT/${n/'_'$READONE.$FASTQ/'_'$READLIB.ash}.bam $MYOUT/${n/'_'$READONE.$FASTQ/'_'$READLIB.map}.bam $MYOUT/${n/'_'$READONE.$FASTQ/'_'$READLIB.unm}.bam 
    
    if [ ! -e $MYOUT/metrices ]; then mkdir -p $MYOUT/metrices ; fi
    THISTMP=$TMP/$n$RANDOM #mk tmp dir because picard writes none-unique files                                        
    mkdir -p $THISTMP
    java $JAVAPARAMS -jar $PATH_PICARD/MarkDuplicates.jar \
        INPUT=$MYOUT/${n/'_'$READONE.$FASTQ/'_'$READLIB.ash.bam} \
        OUTPUT=$MYOUT/${n/'_'$READONE.$FASTQ/'_'$READLIB.$ASD.bam} \
        METRICS_FILE=$MYOUT/metrices/${n/'_'$READONE.$FASTQ/'_'$READLIB.$ASD.bam}.dupl AS=true \
        VALIDATION_STRINGENCY=LENIENT \
        TMP_DIR=$THISTMP
    rm -rf $THISTMP
    samtools index $MYOUT/${n/'_'$READONE.$FASTQ/'_'$READLIB.$ASD.bam}
    
    
    # statistics                                                                                                      
    STATSOUT=$MYOUT/${n/'_'$READONE.$FASTQ/'_'$READLIB.$ASD.bam}.stats
    samtools flagstat $MYOUT/${n/'_'$READONE.$FASTQ/'_'$READLIB.$ASD.bam} > $STATSOUT
    if [ -n $SEQREG ]; then
        echo "#custom region" >> $STATSOUT
        echo `samtools view $MYOUT/${n/'_'$READONE.$FASTQ/'_'$READLIB.ash.bam} $SEQREG | wc -l`" total reads in region " >> $STAT\
    SOUT
        echo `samtools view -f 2 $MYOUT/${n/'_'$READONE.$FASTQ/'_'$READLIB.ash.bam} $SEQREG | wc -l`" properly paired reads in re\
    gion " >> $STATSOUT
    fi
    
    #verify
    BAMREADS=`head -n1 $MYOUT/${n/'_'$READONE.$FASTQ/'_'$READLIB.$ASD.bam}.stats | cut -d " " -f 1`
    if [ "$BAMREADS" = "" ]; then let BAMREADS="0"; fi
    if [ $BAMREADS -eq $FASTQREADS ]; then
        echo "-----------------> PASS check mapping: $BAMREADS == $FASTQREADS"
        rm $MYOUT/${n/'_'$READONE.$FASTQ/'_'$READLIB.$ALN.sam}
        rm $MYOUT/${n/'_'$READONE.$FASTQ/'_'$READLIB.$ALN.un.sam}
        rm $MYOUT/${n/'_'$READONE.$FASTQ/'_'$READLIB.ash.bam}
        rm $MYOUT/${n/'_'$READONE.$FASTQ/'_'$READLIB.unm}.bam
        rm $MYOUT/${n/'_'$READONE.$FASTQ/'_'$READLIB.map}.bam
    else
        echo -e "[ERROR] We are loosing reads from .fastq -> .bam in $f: \nFastq had $FASTQREADS Bam has $BAMREADS"
        exit 1
    fi
    
    #samstat
    samstat $MYOUT/${n/'_'$READONE.$FASTQ/'_'$READLIB.$ASD.bam}
done

echo "********* makeTagDirectory" 
RUN_COMMAND="makeTagDirectory $MYOUT/${n/'_'$READONE.$FASTQ/_tagdir_unfiltered} $MYOUT/${n/'_'$READONE.$FASTQ/'_'$READONE.$ASD.bam} $MYOUT/${n/'_'$READONE.$FASTQ/'_'$READTWO.$ASD.bam} $HOMER_HIC_TAGDIR_OPTIONS"
echo $RUN_COMMAND
eval $RUN_COMMAND


cp -r $MYOUT/${n/'_'$READONE.$FASTQ/_tagdir_unfiltered} $MYOUT/${n/'_'$READONE.$FASTQ/_tagdir_filtered}

RUN_COMMAND="makeTagDirectory $MYOUT/${n/'_'$READONE.$FASTQ/_tagdir_filtered} -update $HOMER_HIC_TAGDIR_OPTIONS"

echo $RUN_COMMAND
eval $RUN_COMMAND

echo "********* analyzeHiC" 




echo ">>>>> HiC analysis with homer - FINISHED"
echo ">>>>> enddate "`date`

