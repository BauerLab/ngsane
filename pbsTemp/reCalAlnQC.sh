
# author: Denis C. Bauer
# date: Dec.2010

QOUT=$1

echo ""
echo "###################################################"
echo "# recalibrating and realigning "
echo "###################################################"


files=`ls $QOUT/*.out | wc -l`

# FINISHED ?
finished=`grep "FINISHED" $QOUT/* | cut -d ":" -f 1 | sort -u | wc -l`
if [ "$finished" = "$files" ]; then
    echo "QC_PASS .. finished are $finished/$files"
else
    echo "**_FAIL .. finished are $finished/$files"
fi

#stats
for i in "PAIRED" "SINGLE"
do
  var=`grep "$i" $QOUT/* | cut -d ":" -f 1 | sort -u | wc -l`
  echo "           $var $i"
done


# Errors
for i in "We are loosing reads" "error" "no such file" "file not found" "reCalAln.sh: line"
do
  var=`grep -i "$i" $QOUT/* | grep -v "Writing out data tables for read group:" | cut -d ":" -f 1 | sort -u | wc -l`
  if [ "$var" = "0" ]; then
    echo "QC_PASS .. $var have $i/$files"
  else
    echo "**_FAIL .. $var have $i/$files"
  fi
done


# Steps
echo ">>>>>>>>>> Progress"
for i in "recalibrating" "counting covariantes" "change score" "index" "QC Step" "counting covariantes after recalibration" "plotting both" "realignment" "find intervals to improve" "realine them" "sort/index" "fix mates" "statistics" "coverage track"
do
  var=`grep "********* $i" $QOUT/* | cut -d ":" -f 1 | sort -u | wc -l`
  verd="**_FAIL"
  if [ "$var" = "$files" ]; then
      verd="QC_PASS"
  fi
  echo "$verd .. $i $var/$files"
done