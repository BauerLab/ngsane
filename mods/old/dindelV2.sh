#!/bin/bash

# Script running dindel
# It expects an aligned bam file, the reference genome as input and requires the
# directory *ind* to be present.


# cause of errors:
# mergeOutputDiploid.py needs absolute path to util/ dir


#INPUTS
HISEQINF=$1   # location of the HiSeqInf repository
CANDIDATES=$2
FASTA=$3      #reference genome
OUT=$4        #output dir

QOUT=qout/

PRIORITY="-l hp=TRUE"

#PROGRAMS
. $HISEQINF/conf/header.sh

#PARAMETERS
WINDOWS=100


echo ">>>>> indel calling with DINDEL"
echo ">>>>> startdate "`date`
echo ">>>>> dindel.sh $HISEQINF $CANDIDATES $FASTA $OUT"

CANDIDATES=${CANDIDATES//,/ }

echo $CADIDATES


# delete privious output files
if [ -e $OUT/dindelMerge.VCF ]; then 
    rm $OUT/dindelMerge.VCF
    touch 
fi

if [ -e merge.tmp ]; then
    rm merge.tmp
fi

#get candidates from previous run
LIST=""
for c in $CANDIDATES; do
    echo $c

    touch ${c/"dindel.VCF"/"dindelP2.VCF"} # make for the merge 

    python $DINDELHOME/dindel-1.01-python/convertVCFToDindel.py --inputFile $c \
	--outputFile $c.pos --refFile $FASTA
    LIST=$LIST" "$c.pos
done

#get in one file
cat $LIST | grep -w chr1 | sort -u > $OUT/dindel.comVariants.txt


if [ -e merge.tmp ]; then rm merge.tmp; fi

#align them to the right
for c in $CANDIDATES; do
    name=`basename $c`
    qsub $PRIORITY -j y -o $QOUT/$TASKDINC/$name.out -cwd -b y -l h_vmem=12G \
    	-N $TASKDINC"_"$name $HOLD\
    	$HISEQINF/mods/dindelV2p2.sh $HISEQINF $c $FASTA $OUT

    gawk -v OFS='\t' '{$5=toupper($5);$4=toupper($4); print $_}' ${c/"dindel.VCF"/"dindelP2.VCF"} >${c/"dindel.VCF"/"dindelP3.VCF"}

    echo ${c/"dindel.VCF"/"dindelP3.VCF"} >>merge.tmp

done


$HISEQINF/mods/merge.sh $HISEQINF merge.tmp $OUT dindelMerge.vcf vcf $QOUT/$TASKDINC $TASKDINC"*"

echo ">>>>> indel calling with DINDEL - FINISHED"
echo ">>>>> enddate "`date`