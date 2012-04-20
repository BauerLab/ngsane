

QOUT=$1

files=`ls $QOUT/*.err | wc -l`

# FINISHED ?
finished=`grep "FINISHED" $QOUT* | cut -d ":" -f 1 | sort -u | wc -l`
if [ "$finished" = "$files" ]; then
    echo "QC_PASS .. finished are $finished"
else
    echo "QC_FAIL .. finished are $finished"
fi


# Errors
for i in "error_above_read_count_threshold" "matepos inconsistency" "Exception: Cannot open variant file"
do
  var=`grep "$i" $QOUT* | cut -d ":" -f 1 | sort -u | wc -l`
  if [ "$var" = "0" ]; then
    echo "QC_PASS .. $var have $i/$files"
  else
    echo "QC_FAIL .. $var have $i/$files"
  fi
done


# Steps
echo "** Progress **"
for i in "get candidates" "make windows" "realign windows" "merge"
do
  var=`grep "********* $i" $QOUT* | cut -d ":" -f 1 | sort -u | wc -l`
  verd="QC_FAIL"
  if [ "$var" = "$files" ]; then
      verd="QC_PASS"
  fi
  echo "$verd .. $i $var/$files"
done