#!/bin/bash

# Script running dindel
# It expects an aligned bam file, the reference genome as input and requires the
# directory *ind* to be present.


# cause of errors:
# mergeOutputDiploid.py needs absolute path to util/ dir


#INPUTS
NGSANE_BASE=$1   # location of the NGSANE repository
f=$1          #fastq file
FASTA=$2      #reference genome
OUT=$3        #output dir
CANDIDATES=$4

#PROGRAMS
. ${NGSANE_BASE}/conf/header.sh

#PARAMETERS
WINDOWS=1000


n=`basename $f`
echo ">>>>> indel calling with DINDEL"
echo ">>>>> startdate "`date`
echo ">>>>> dindel.sh $f $FASTA $OUT $MERGE"


# delete privious output files
if [ -e $OUT/$n.dindel.VCF ]; then 
    rm $OUT/$n.dindel.VCF
fi

# make dir
if [ ! -d $OUT/dindelWindow ]; then 
    mkdir $OUT/dindelWindow
fi

# index bam file
if [ ! -e $f.bai ]; then
    $SAMTOOLS index $f
fi


# get candidates
echo "********* get candidates"
$DINDELHOME/binaries/dindel-1.01-linux-64bit --analysis getCIGARindels --bamFile $f \
    --outputFile $OUT/$n.dindel --ref $FASTA --quiet --maxRead 100000 


# make windows
echo "********* make windows"
$DINDELHOME/dindel-1.01-python/makeWindows.py \
    --inputVarFile $OUT/$n.dindel.variants.txt \
    --windowFilePrefix $OUT/dindelWindow/$n.realWind \
    --numWindowsPerFile $WINDOWS

#realign windows
echo "********* realign windows"
for i in $( ls $OUT/dindelWindow/$n.realWind* )
  do
  echo "prozess $i"
  $DINDELHOME/binaries/dindel-1.01-linux-64bit --analysis indels --doDiploid \
      --bamFile $f --ref $FASTA --quiet --maxRead 100000 \
      --varFile $i \
      --libFile $OUT/$n.dindel.libraries.txt \
      --outputFile ${i/realWind/stage2_realWind}
done
ls $OUT/dindelWindow/$n*.glf.txt > $OUT/dindelWindow/$n.staget2.outputfiles.txt

#merge
echo "********* merge"
$DINDELHOME/dindel-1.01-python/mergeOutputDiploid.py \
    --inputFiles $OUT/dindelWindow/$n.staget2.outputfiles.txt\
    --outputFile $OUT/$n.dindel.VCF --ref $FASTA

#clean up
rm $OUT/$n.dindel.variants.txt
rm $OUT/$n.dindel.libraries.txt
rm $OUT/dindelWindow/$n.stage2"_"realWind.*
rm $OUT/dindelWindow/$n.realWind*
rm $OUT/dindelWindow/$n".staget2.outputfiles.txt"


echo ">>>>> indel calling with DINDEL - FINISHED"
echo ">>>>> enddate "`date`