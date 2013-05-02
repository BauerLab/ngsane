#!/bin/bash


# Script running dindel on a bam file containing the reads from several individuals
# It requires the candidate indels from a prior dindel run on the bamfiles for each
# individual seperately.
# author: Denis C. Bauer
# date: Jan.2011


#INPUTS
HISEQINF=$1   # location of the HiSeqInf repository
FILE=$2       # bam files
CANDIDATES=$3 # candidate regions from the seperate run
FASTA=$4      # reference genome
OUT=$5        # output dir

#PROGRAMS
. $HISEQINF/conf/header.sh

#PARAMETERS
WINDOWS=1000


echo ">>>>> indel calling with DINDEL simultaneously"
echo ">>>>> startdate "`date`
echo ">>>>> dindelMerge.sh $HISEQINF $FILE $CANDIDATES $FASTA $OUT"


# delete privious output files
if [ -e $OUT/genotype.dindel.vcf ]; then 
    rm $OUT/genotype.dindel.vcf
fi

# make dir
if [ ! -d $OUT/dindelWindow ]; then 
    mkdir $OUT/dindelWindow
fi

samp=12


# get library size
echo "********* get library size"
$DINDELHOME/binaries/dindel-1.01-linux-64bit --analysis getCIGARindels --bamFile $FILE \
    --outputFile $FILE.dindel --ref $FASTA --maxRead 100000 #--region 229994688-230071581 --tid chr1

# #combine the variants from the previous run
# cat $FILE.dindel.variants.txt $CANDIDATES > $FILE.dindel.ext.variants.txt

# make windows
echo "********* make windows"
$DINDELHOME/dindel-1.01-python/makeWindows.py \
    --inputVarFile $FILE.dindel.variants.txt \
    --windowFilePrefix $OUT/dindelWindow/genotype.realWind \
    --numWindowsPerFile $WINDOWS

#get indels
#
#@Kees: this is where the problem occurs
echo "********* get indels"
for i in $( ls $OUT/dindelWindow/genotype.realWind* )
  do
  echo "prozess $i"
  $DINDELHOME/binaries/dindel-1.01-linux-64bit --analysis indels --doPooled \
      --bamFile $FILE --ref $FASTA --quiet --maxRead 100000 \
      --varFile $FILE.dindel.libraries.txt \
      --libFile $OUT/ \
      --outputFile ${i/realWind/stage2_realWind}
done
ls $OUT/dindelWindow/genotype*.glf.txt > $OUT/dindelWindow/genotype.stage2.outputfiles.txt

exit

# merge into a VCF file
echo "********* merge"
$DINDELHOME/dindel-1.01-python/mergeOutputPooled.py --numSamples $samp --numBAMFiles $samp \
    --inputFiles $OUT/dindelWindow/genotype.stage2.outputfiles.txt\
    --outputFile $OUT/genotype.dindel.vcf --ref $FASTA

# get individual genotypes
echo "********* get individual genotypelikelyhoods"
$DINDELHOME/dindel-1.01-python/makeGenotypeLikelihoodFilePooled.py \
    --inputGLFFiles ${i/realWind/stage2_realWind} \
    --callFile $OUT/genotype.dindel.vcf \
    --bamFiles $FILES \
    --outputFile $OUT/genotype.dindel.glf


#clean up
#rm $OUT/genotype.dindel.variants.txt
#rm $OUT/genotype.dindel.libraries.txt
#rm $OUT/dindelWindow/genotype.stage2"_"realWind.*
#rm $OUT/dindelWindow/genotype.realWind*
#rm $OUT/dindelWindow/genotype".staget2.outputfiles.txt"
#rm bamFilePaths.txt


echo ">>>>> indel calling with DINDEL - FINISHED"
echo ">>>>> enddate "`date`