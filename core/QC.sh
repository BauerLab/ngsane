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
        -t | --task )           shift; TASK=$1 ;;   # location of log output
        -o | --output-html )    HTMLOUTPUT='TRUE';; # flag                                                     
        -h | --help )           usage ;;
        * )                     echo "don't understand "$1
    esac
    shift
done

################################################################################

if [ -n "$DMGET" ]; then
    dmget $QOUT/$TASK/*.out
fi

LOGFOLDER=$(dirname $QOUT)

if [ "$HTMLOUTPUT" == 'TRUE' ]; then
    echo "<div id='${TASK}_checklist'><pre>"
fi

echo "###################################################"
echo "# NGSANE ${SCRIPT/.sh/} "
echo "###################################################"

files=$(ls $QOUT/$TASK/*.out | wc -l)

CHECKPOINTS_PASSED=0
CHECKPOINTS_FAILED=0

# FINISHED ?
finished=$(tail -n 2  $QOUT/$TASK/*.out | grep "FINISHED" | cut -d ":" -f 1 | sort -u | wc -l)
if [ "$finished" = "$files" ]; then
    echo "QC_PASS .. finished are $finished/$files"
    CHECKPOINTS_PASSED=`expr $CHECKPOINTS_PASSED + 1`
else
    echo "**_FAIL .. finished are $finished/$files"
    CHECKPOINTS_FAILED=`expr $CHECKPOINTS_FAILED + 1`
fi

IFS=','

echo ">>>>>>>>>> Errors"
ERROR=$(grep QCVARIABLES $SCRIPT)
ERROR=${ERROR/"# QCVARIABLES,"/}
for i in $ERROR; do
  var=$(grep -i "$i" $QOUT/$TASK/*.out | cut -d ":" -f 1 | sort -u | wc -l)
  if [ "$var" = "0" ]; then
    echo "QC_PASS .. $var have $i/$files"
    CHECKPOINTS_PASSED=`expr $CHECKPOINTS_PASSED + 1`
  else
    echo "**_FAIL .. $var have $i/$files"
    CHECKPOINTS_FAILED=`expr $CHECKPOINTS_FAILED + 1`
  fi
done

echo ">>>>>>>>>> CheckPoints "
PROGRESS=`grep -P '^CHECKPOINT="' $SCRIPT | awk -F'"' '{print $2}' | tr '\n' ','`
for i in $PROGRESS; do
  var=$(grep -i "\*\*\*\*\*\* $i" $QOUT/$TASK/*.out | cut -d ":" -f 1 | sort -u | wc -l)
  if [ ! "$var" = $files ]; then
    echo "**_FAIL .. $var have $i/$files"
    CHECKPOINTS_FAILED=`expr $CHECKPOINTS_FAILED + 1`
  else
    echo "QC_PASS .. $var have $i/$files"
    CHECKPOINTS_PASSED=`expr $CHECKPOINTS_PASSED + 1`
  fi
done

if [ "$HTMLOUTPUT" = 'TRUE' ]; then
    echo "</pre></div>"
    echo "<div id='${TASK}_notes'><pre>"
fi

echo ">>>>>>>>>> Notes"
SUMNOTES=0
for i in $QOUT/$TASK/*.out;do
    echo ${i/$LOGFOLDER\//}
    NOTELIST=$(grep -P "^\[NOTE\]" $i)
    echo $NOTELIST
    SUMNOTES=`expr $SUMNOTES + $(echo $NOTELIST | awk 'BEGIN{count=0} NF != 0 {++count} END {print count}' )`
done

echo ">>>>>>>>>> Errors"
SUMERRORS=0
for i in $QOUT/$TASK/*.out;do
    echo ${i/$LOGFOLDER\//}
    ERRORLIST=$(grep -P "^\[ERROR\]" $i)
    echo $ERRORLIST
    SUMERRORS=`expr $SUMERRORS + $(echo $ERRORLIST | awk 'BEGIN{count=0} NF != 0 {++count} END {print count}' )`
done

if [ "$HTMLOUTPUT" = 'TRUE' ]; then
    echo "</pre></div>"
    echo "<script type='text/javascript'> if (typeof jQuery === 'undefined') {console.log('jquery not loaded');} else {\$(\"#${TASK}_counter_notes\").text('$SUMNOTES');\$('#${TASK}_counter_errors').text('$SUMERRORS');\$('#${TASK}_counter_checkpoints_passed').text('$CHECKPOINTS_PASSED');\$('#${TASK}_counter_checkpoints_failed').text('$CHECKPOINTS_FAILED'); if ($SUMERRORS==0){\$('#${TASK}_counter_errors').toggleClass('errors neutral');}; if($CHECKPOINTS_FAILED==0){\$('#${TASK}_counter_checkpoints_failed').toggleClass('failed nofailed');};}</script>"
    
    echo "<div id='${TASK}_logfiles'>"
    for i in $QOUT/$TASK/*.out ;do
        FN=${i/$LOGFOLDER\//}
        echo "<a href='$FN'>${i/$QOUT\/$TASK\//}</a><br/>"
    done
    echo "</div>"
fi


################################################################################

