# author: Denis C. Bauer
# date: Feb.2011
# modified by Fabian Buske, Jul 2013
SCRIPT=$1
QOUT=$2

if [ -n "$DMGET" ]; then
    dmget $QOUT/*.out
fi

echo ""
echo "###################################################"
echo "# NGSANE ${SCRIPT/.sh/} "
echo "###################################################"

files=$(ls $QOUT/*.out | wc -l)

# FINISHED ?
finished=$(grep "FINISHED" $QOUT/*.out | cut -d ":" -f 1 | sort -u | wc -l)
if [ "$finished" = "$files" ]; then
    echo "QC_PASS .. finished are $finished/$files"
else
    echo "**_FAIL .. finished are $finished/$files"
fi

IFS=','

echo ">>>>>>>>>> Errors"

# Errors
ERROR=$(grep QCVARIABLES $1)
ERROR=${ERROR/"# QCVARIABLES,"/}
for i in $ERROR
do
  var=$(grep -i "$i" $QOUT/*.out | cut -d ":" -f 1 | sort -u | wc -l)
  if [ "$var" = "0" ]; then
    echo "QC_PASS .. $var have $i/$files"
  else
    echo "**_FAIL .. $var have $i/$files"
  fi
done

echo ">>>>>>>>>> CheckPoints "

PROGRESS=`grep -P '^CHECKPOINT="' $SCRIPT | awk -F'"' '{print $2}' | tr '\n' ','`

for i in $PROGRESS
do
  var=$(grep -i "\*\*\*\*\*\* $i" $QOUT/*.out | cut -d ":" -f 1 | sort -u | wc -l)
  if [ ! "$var" = $files ]; then
    echo "**_FAIL .. $var have $i/$files"
  else
    echo "QC_PASS .. $var have $i/$files"
  fi
done
