
# author: Denis C. Bauer
# date: Dec.2010

QOUT=$1

echo ""
echo "###################################################"
echo "# GATK indel calling "$QOUT
echo "###################################################"


files=`ls $QOUT/*.out | wc -l`

# FINISHED ?
finished=`grep "FINISHED" $QOUT/* | cut -d ":" -f 1 | sort -u | wc -l`
if [ "$finished" = "$files" ]; then
    echo "QC_PASS .. finished are $finished/$files"
else
    echo "**_FAIL .. finished are $finished/$files"
fi

# Errors
for i in "error"
do
  var=`grep "$i" $QOUT/* | cut -d ":" -f 1 | sort -u | wc -l`
  if [ "$var" = "0" ]; then
    echo "QC_PASS .. $var have $i/$files"
  else
    echo "**_FAIL .. $var have $i/$files"
  fi
done


# Steps
echo ">>>>>>>>>> Progress"
for i in "call indels" "make index for IGV"
do
  var=`grep "********* $i" $QOUT/* | cut -d ":" -f 1 | sort -u | wc -l`
  verd="**_FAIL"
  if [ "$var" = "$files" ]; then
      verd="QC_PASS"
  fi
  echo "$verd .. $i $var/$files"
done