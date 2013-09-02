#!/bin/sh

# author: Denis C. Bauer
# date: Feb.2011
# modified by Fabian Buske, Aug 2013

function usage {
echo -e "usage: $(basename $0) [OPTIONS] -m MOD.SCRIPT - l LOG.QOUT"
exit
}

if [ ! $# -gt 2 ]; then usage ; fi

while [ "$1" != "" ]; do
    case $1 in
        -m | --mod )            shift; SCRIPT=$1 ;; # location of mod script                       
        -l | --log )            shift; QOUT=$1 ;;   # location of log output
        -o | --output-html )    HTMLOUTPUT='TRUE';; # flag                                                     
        -h | --help )           usage ;;
        * )                     echo "don't understand "$1
    esac
    shift
done

################################################################################

if [ -n "$DMGET" ]; then
    dmget $QOUT/*.out
fi


if [ "$HTMLOUTPUT" == 'TRUE' ]; then
    echo "<pre>"
fi

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
ERROR=$(grep QCVARIABLES $SCRIPT)
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

if [ "$HTMLOUTPUT" = 'TRUE' ]; then
    echo "</pre>"
fi

################################################################################

