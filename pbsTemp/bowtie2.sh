#!/bin/bash

echo ">>>>> alignment with with bowtie "
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> bowtie.sh $*"

function usage {
echo -e "usage: $(basename $0) -k HISEQINF -f FASTQ -r REFERENCE -o OUTDIR [OPTIONS]"
exit
}

# Script to run bowtie program.
# It takes comma-seprated list of files containing short sequence reads in fasta or fastq format and bowtie index files as input.
# It produces output files: read alignments in .bam format and other files.
# author: Denis Bauer
# date: June. 2012

# QCVARIABLES,Resource temporarily unavailable

if [ ! $# -gt 3 ]; then usage ; fi

THREADS=1
MEMORY=2
EXPID="exp"           # read group identifier RD ID                                                               
LIBRARY="qbi"         # read group library RD LB                                                                  
PLATFORM="illumina"   # read group platform RD PL                                                                 
UNIT="flowcell"       # read group platform unit RG PU                                                            
FASTQNAME=""
FORCESINGLE=0

#INPUTS                                                                                                           
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the HiSeqInf repository                       
        -t | --threads )        shift; THREADS=$1 ;; # number of CPUs to use                                      
        -m | --memory )         shift; MEMORY=$1 ;; # memory used 
        -f | --fastq )          shift; f=$1 ;; # fastq file                                                       
        -r | --reference )      shift; FASTA=$1 ;; # reference genome                                             
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

#PROGRAMS
. $CONFIG
. $HISEQINF/pbsTemp/header.sh
. $CONFIG
#if [ -n "$FASTQNAME" ]; then FASTQ=$FASTQNAME ; fi

module load R
module load jdk #/1.7.0_03
#module load $BOWTIETWO
JAVAPARAMS="-Xmx"$MEMORY"g -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1 -XX:MaxDirectMemorySize=4G"

n=`basename $f`

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
if [[ $f != *.fastq.gz ]]; then ZCAT="cat"; fi

# generating the index files                                                                                      
if [ ! -e ${FASTA/.fasta/}.1.bt2 ]; then echo ">>>>> make .bt2"; $BOWTIETWO/bowtie2-build $FASTA ${FASTA/.fasta/}; fi
if [ ! -e $FASTA.fai ]; then echo ">>>>> make .fai"; $SAMTOOLS faidx $FASTA; fi

dmget -a $(dirname $FASTA)/*

#run bowtie command -v $MISMATCH -m 1
echo "********* bowtie" 
if [ $PAIRED == "0" ]; then 
    READS="-U $f"
    dmget -a $f
    let FASTQREADS=`$ZCAT $f | wc -l | gawk '{print int($1/4)}' `
else 
    READS="-1 $f -2 ${f/$READONE/$READTWO}"
    dmget -a $f
    dmget -a ${f/$READONE/$READTWO}
    READ1=`$ZCAT $f | wc -l | gawk '{print int($1/4)}' `
    READ2=`$ZCAT ${f/$READONE/$READTWO} | wc -l | gawk '{print int($1/4)}' `
    let FASTQREADS=$READ1+$READ2
fi

#readgroup
FULLSAMPLEID=$SAMPLEID"${n/'_'$READONE.$FASTQ/}"
RG='--sam-rg "ID:$EXPID" --sam-rg "SM:$FULLSAMPLEID" --sam-rg "LB:$LIBRARY" --sam-rg "PL:$PLATFORM"'

$BOWTIETWO/bowtie2 $RG -t -x ${FASTA/.fasta/} -p $THREADS $READS -S $MYOUT/${n/'_'$READONE.$FASTQ/.$ALN.sam} --un $MYOUT/${n/'_'$READONE.$FASTQ/.$ALN.un.sam}
#$BOWTIETWO/bowtie2 -t -x ~/Documents/datahome/ErrorCorrection/ExCap/bowtie_index_konsta/human_g1k_v37 -p $THREADS $READS -S $MYOUT/${n/'_'$READONE.$FASTQ/.$ALN.sam} --un $MYOUT/${n/'_'$READONE.$FASTQ/.$ALN.un.sam}


# continue for normal bam file conversion                                                                         
echo "********* sorting and bam-conversion"
$SAMTOOLS view -bt $FASTA.fai $MYOUT/${n/'_'$READONE.$FASTQ/.$ALN.sam} | $SAMTOOLS sort - $MYOUT/${n/'_'$READONE.$FASTQ/.map}
$SAMTOOLS view -bt $FASTA.fai $MYOUT/${n/'_'$READONE.$FASTQ/.$ALN.un.sam} | $SAMTOOLS sort - $MYOUT/${n/'_'$READONE.$FASTQ/.unm}
# merge mappend and unmapped
$SAMTOOLS merge -f $MYOUT/${n/'_'$READONE.$FASTQ/.ash}.bam $MYOUT/${n/'_'$READONE.$FASTQ/.map}.bam $MYOUT/${n/'_'$READONE.$FASTQ/.unm}.bam 


echo "********* mark duplicates"
if [ ! -e $MYOUT/metrices ]; then mkdir $MYOUT/metrices ; fi
THISTMP=$TMP/$n$RANDOM #mk tmp dir because picard writes none-unique files                                        
mkdir $THISTMP
java $JAVAPARAMS -jar $PICARD/MarkDuplicates.jar \
    INPUT=$MYOUT/${n/'_'$READONE.$FASTQ/.ash.bam} \
    OUTPUT=$MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam} \
    METRICS_FILE=$MYOUT/metrices/${n/'_'$READONE.$FASTQ/.$ASD.bam}.dupl AS=true \
    VALIDATION_STRINGENCY=LENIENT \
    TMP_DIR=$THISTMP
rm -r $THISTMP
$SAMTOOLS index $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}



# statistics                                                                                                      
echo "********* statistics"
STATSOUT=$MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}.stats
$SAMTOOLS flagstat $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam} > $STATSOUT
if [ -n $SEQREG ]; then
    echo "#custom region" >> $STATSOUT
    echo `$SAMTOOLS view $MYOUT/${n/'_'$READONE.$FASTQ/.ash.bam} $SEQREG | wc -l`" total reads in region " >> $STAT\
SOUT
    echo `$SAMTOOLS view -f 2 $MYOUT/${n/'_'$READONE.$FASTQ/.ash.bam} $SEQREG | wc -l`" properly paired reads in re\
gion " >> $STATSOUT
fi

echo "********* verify"
BAMREADS=`head -n1 $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}.stats | cut -d " " -f 1`
if [ "$BAMREADS" = "" ]; then let BAMREADS="0"; fi
if [ $BAMREADS -eq $FASTQREADS ]; then
    echo "-----------------> PASS check mapping: $BAMREADS == $FASTQREADS"
    rm $MYOUT/${n/'_'$READONE.$FASTQ/.$ALN.sam}
    rm $MYOUT/${n/'_'$READONE.$FASTQ/.$ALN.un.sam}
    rm $MYOUT/${n/'_'$READONE.$FASTQ/.ash.bam}
    rm $MYOUT/${n/'_'$READONE.$FASTQ/.unm}.bam
    rm $MYOUT/${n/'_'$READONE.$FASTQ/.map}.bam
else
    echo -e "***ERROR**** We are loosing reads from .fastq -> .bam in $f: \nFastq had $FASTQREADS Bam has $BAMREA\
DS"
    exit 1

fi

echo "********* coverage track"
GENOME=$(echo $FASTA| sed 's/.fasta/.genome/' | sed 's/.fa/.genome/' )
java $JAVAPARAMS -jar $IGVTOOLS count $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam} \
    $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam.cov.tdf} $GENOME


echo "********* samstat"
$SAMSTAT $MYOUT/${n/'_'$READONE.$FASTQ/.$ASD.bam}

echo ">>>>> readmapping with BWA - FINISHED"
echo ">>>>> enddate "`date`

