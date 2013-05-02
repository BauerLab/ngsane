#!/bin/bash

# Script running dindel
# It expects an aligned bam file, the reference genome as input and requires the
# directory *ind* to be present.


# cause of errors:
# mergeOutputDiploid.py needs absolute path to util/ dir


#INPUTS
HISEQINF=$1   # location of the HiSeqInf repository
c=$2
FASTA=$3      #reference genome
OUT=$4        #output dir



#PROGRAMS
. $HISEQINF/conf/header.sh

#PARAMETERS
WINDOWS=100


echo ">>>>> indel calling with DINDEL"
echo ">>>>> startdate "`date`
echo ">>>>> dindel.sh $HISEQINF $c $FASTA $OUT"

n=`basename $c`
d=`dirname $c`

# make dir
if [ ! -d $d/dindelWindow ]; then 
    mkdir $d/dindelWindow
fi


BAM=${c/$TASKDIN/$TASKBWA}
BAM=${BAM/".dindel.VCF"/}
BAM=${BAM/ashrr/asd}


# get candidates
echo "********* get candidates"
$DINDELHOME/binaries/dindel-1.01-linux-64bit --analysis getCIGARindels --bamFile $BAM \
    --outputFile $c.dindel --ref $FASTA --quiet --maxRead 1000000 

# concatinate them with previous ones
#cat $OUT/dindel.comVvariants.txt $c.dindel.variants.txt | grep -w chr1 | sort -u > $c.dindel.combVariants.txt 
#ASSESSPOS=$c.dindel.combVariants.txt

ASSESSPOS=$OUT/dindel.comVariants.txt


#realign them
echo "********* realign them $ASSESSPOS"
$DINDELHOME/binaries/dindel-1.01-linux-64bit --analysis realignCandidates \
    --varFile $ASSESSPOS --outputFile $c.dindel.combVariants.real.txt --ref $FASTA


# make windows
echo "********* make windows"
$DINDELHOME/dindel-1.01-python/makeWindows.py \
    --inputVarFile $c.dindel.combVariants.real.txt.variants.txt \
    --windowFilePrefix $d/dindelWindow/$n.realWind \
    --numWindowsPerFile $WINDOWS
    #--inputVarFile $c.dindel.variants.txt \


#call indels
echo "********* call indels"
for i in $( ls $d/dindelWindow/$n.realWind* )
  do
  echo "prozess $i"
  $DINDELHOME/binaries/dindel-1.01-linux-64bit --analysis indels --doDiploid \
      --bamFile $BAM --ref $FASTA --quiet --maxRead 1000000 \
      --varFile $i \
      --libFile $c.dindel.libraries.txt \
      --outputFile ${i/realWind/stage2_realWind}
done
ls $d/dindelWindow/$n*.glf.txt > $d/dindelWindow/$n.staget2.outputfiles.txt

#merge
echo "********* merge "${c/"dindel.VCF"/"dindelP2.VCF"}
$DINDELHOME/dindel-1.01-python/mergeOutputDiploid.py \
    --inputFiles $d/dindelWindow/$n.staget2.outputfiles.txt \
    --outputFile ${c/"dindel.VCF"/"dindelP2.VCF"} --ref $FASTA


#clean up
rm $d/dindelWindow/$n.stage2"_"realWind.*
rm $d/dindelWindow/$n.realWind*
rm $d/dindelWindow/$n".staget2.outputfiles.txt"
rm $c.pos
rm $c.dindel.variants.txt
rm $c.dindel.libraries.txt
rm $c.dindel.combVariants.real.txt.variants.txt
#rm $c.dindel.combVariants.txt

echo ">>>>> indel calling with DINDEL - FINISHED"
echo ">>>>> enddate "`date`