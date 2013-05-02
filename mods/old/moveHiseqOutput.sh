
# moves the essential files from a hiSeq run
#
# author: Denis C. Bauer
# date: Dec.2010


function usage {
    echo "moveHiseqOutput.sh <Run dir> <root of destination> <lane1(,lane2)>"
    exit
}


if [ ! $# -gt 2 ]; then usage ; fi


WD=$1
DEST=$2
LANES=$3

LATESTBASECALLDIR=`ls $WD/Data/Intensities/ | grep "BaseCalls_" | tail -n1`
LATESTGERALDDIR=`ls $WD/Data/Intensities/$LATESTBASECALLDIR | grep "GERALD_" | tail -n1`

SOURCE="$WD/Data/Intensities/$LATESTBASECALLDIR/$LATESTGERALDDIR/"

RUN=`basename $WD`

if [ ! -d $DEST/$RUN ]; then
    mkdir $DEST/$RUN
    mkdir $DEST/$RUN/fastq
    mkdir $DEST/$RUN/runStats
else
    echo "Folder already there"
   # exit 1
fi


echo "copy $SOURCE/*_sequences.txt to $DEST/$RUN/fastq/ $LANES"

for f in ${LANES//,/ }; do
  echo $f
  qsub -l hp=TRUE -b y -j y -o $WD/MyLogs/move"_s_"$f".out" -N "s_"$f \
      cp $SOURCE/"s_"$f"_1_sequence.txt" $DEST/$RUN/fastq/"s_"$f"_read1.fastq"
  if [ -e $SOURCE/"s_"$f"_2_sequence.txt" ]; then
      echo "second"
      qsub -l hp=TRUE -b y -j y -o $WD/MyLogs/move"_s_"$f".out" -N "s_"$f \
	  cp $SOURCE/"s_"$f"_2_sequence.txt" $DEST/$RUN/fastq/"s_"$f"_read2.fastq"
  fi
done

qsub -l hp=TRUE -b y -j y -o $WD/MyLogs/move"_s_"$f".out" -hold_jid "s_*" \
    chmod 444 $DEST/$RUN/fastq/"s_"$f*.fastq

echo "copy other stuff over"
for i in InterOp RunInfo.xml runParameters.xml /Data/Intensities/$LATESTBASECALLDIR/BustardSummary.xml /Data/Intensities/$LATESTBASECALLDIR/$LATESTGERALDDIR/Summary.xml
do
  echo $i
  cp -r $WD/$i $DEST/$RUN/runStats
done

/opt/rocks/bin/chgrp -R hiseq $DEST/$RUN


echo 'AuthType Basic
AuthName "Restricted Access" 
AuthUserFile /home/secure/.htpasswd
Require user cbguser public' > $DEST/$RUN/.htaccess
