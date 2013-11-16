#!/bin/bash

# author: Denis C. Bauer
# date: Feb.2011
# modified by Fabian Buske, Aug 2013

function usage {
echo -e "usage: $(basename $0) [OPTIONS] -m MOD.SCRIPT -l QOUT -t TASK"
exit
}

if [ ! $# -gt 2 ]; then usage ; fi

HTMLOUTPUT=""
RESULTSUFFIX=""

while [ "$1" != "" ]; do
    case $1 in
        -m | --modscript )      shift; SCRIPT=$1 ;; # location of mod script                       
        -l | --log )            shift; QOUT=$1 ;;   # location of log output
        -t | --task )           shift; TASK=$1 ;;   # task at hand
        -r | --results-dir )    shift; OUTDIR=$1 ;; # location of the output 
        -s | --results-task )   shift; OUTTASK=$1 ;; # task the results are put in 
        -f | --filesuffix )     shift; RESULTSUFFIX=$1 ;; # suffix of the result file
        -o | --html-file )      shift; HTMLOUTPUT=$1;; # where the output will be place in the end
        -h | --help )           usage ;;
        * )                     echo "don't understand "$1
    esac
    shift
done

# default output location corresponds to the task at hand
if [ -z "$OUTTASK" ]; then
    OUTTASK=$TASK
fi

################################################################################

if [ -n "$DMGET" ]; then
    dmget $QOUT/$TASK/*.out
fi

LOGFOLDER=$(dirname $QOUT)

if [ -n "$HTMLOUTPUT" ]; then
    echo "<div id='${TASK}_checklist'><pre>"
fi

echo "###################################################"
echo "# NGSANE ${SCRIPT/.sh/} "
echo "###################################################"

files=$(ls $QOUT/$TASK/*.out | wc -l)

CHECKPOINTS_PASSED=0
CHECKPOINTS_FAILED=0

# FINISHED ?
finished=$(grep -P "^>{5} .* FINISHED" $QOUT/$TASK/*.out | cut -d ":" -f 1 | sort -u | wc -l)
if [ "$finished" = "$files" ]; then
    echo "QC_PASS .. finished are $finished/$files"
    CHECKPOINTS_PASSED=`expr $CHECKPOINTS_PASSED + 1`
else
    echo "**_FAIL .. finished are $finished/$files"
    CHECKPOINTS_FAILED=`expr $CHECKPOINTS_FAILED + 1`
fi


echo ">>>>>>>>>> Errors"
DEFAULTSEP=$IFS
IFS=','
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
PROGRESS=$(grep -P '^CHECKPOINT="' $SCRIPT | awk -F'"' '{print $2}' | tr '\n' ',')
for i in $PROGRESS; do
  var=$(grep -P "^\*{9} $i" $QOUT/$TASK/*.out | cut -d ":" -f 1 | sort -u | wc -l)
  if [ ! "$var" = $files ]; then
    echo "**_FAIL .. $var have $i/$files"
    CHECKPOINTS_FAILED=`expr $CHECKPOINTS_FAILED + 1`
  else
    echo "QC_PASS .. $var have $i/$files"
    CHECKPOINTS_PASSED=`expr $CHECKPOINTS_PASSED + 1`
  fi
done


if [ -n "$HTMLOUTPUT" ]; then
    echo "</pre></div>"
    echo "<div id='${TASK}_notes'><pre>"
else
    echo ">>>>>>>>>> Notes"
fi

SUMNOTES=0
for i in $QOUT/$TASK/*.out;do
    echo -e "\n${i/$LOGFOLDER\//}"
    NOTELIST=$(grep -P "^\[NOTE\]" $i)
    echo $NOTELIST
    SUMNOTES=`expr $SUMNOTES + $(echo $NOTELIST | awk 'BEGIN{count=0} NF != 0 {++count} END {print count}' )`
done

if [ -n "$HTMLOUTPUT" ]; then
    echo "</pre></div>"
    echo "<div id='${TASK}_errors'><pre>"
else
    echo ">>>>>>>>>> Errors"
fi

SUMERRORS=0
for i in $QOUT/$TASK/*.out;do
    echo -e "\n${i/$LOGFOLDER\//}"
    ERRORLIST=$(grep -P "^\[ERROR\]" $i)
    if [ -n "$ERRORLIST" ]; then 
        echo $ERRORLIST
    else
        echo "-- all good, no errors"
    fi
    SUMERRORS=$(expr $SUMERRORS + $(echo $ERRORLIST | awk 'BEGIN{count=0} NF != 0 {++count} END {print count}' ))
done

IFS=$DEFAULTSEP

if [ -n "$HTMLOUTPUT" ]; then

    echo "</pre></div>"
    echo "<div id='${TASK}_logfiles'><div class='box'>"
    for i in $( ls $QOUT/$TASK/*.out ); do
        FN=$(python -c "import os.path; print os.path.relpath(os.path.realpath('$i'),os.path.realpath('$(dirname $HTMLOUTPUT)'))")
        echo "<a href='$FN'>${i/$QOUT\/$TASK\//}</a><br/>"
    done
    echo "</div></div>"
    if [ -n "$RESULTSUFFIX" ]; then
        echo "<div id='${TASK}_files'><div class='box'>"
        for i in $(find $OUTDIR/*/$OUTTASK/ -maxdepth 2 -type f -name *$RESULTSUFFIX ); do
            FN=$(python -c "import os.path; print os.path.relpath(os.path.realpath('$i'),os.path.realpath('$(dirname $HTMLOUTPUT)'))")
            echo "<a href='$FN'>${i/$OUTDIR\/*\/$OUTTASK\//}</a><br/>"
        done
        echo "</div></div>"
    fi
    echo "<script type='text/javascript'> 
        if (typeof jQuery === 'undefined') {
            console.log('jquery not loaded');
        } else {
            \$('#${TASK}_counter_notes').text('$SUMNOTES');
            \$('#${TASK}_counter_errors').text('$SUMERRORS');
            \$('#${TASK}_counter_checkpoints_passed').text('$CHECKPOINTS_PASSED');
            \$('#${TASK}_counter_checkpoints_failed').text('$CHECKPOINTS_FAILED'); 
            if ($SUMERRORS==0){
                \$('#${TASK}_counter_errors').toggleClass('errors neutral');
            }; 
            if($CHECKPOINTS_FAILED==0){
                \$('#${TASK}_counter_checkpoints_failed').toggleClass('failed nofailed');
            };

            \$('#${TASK}_panelback h2').click(function() {
                \$(this).parent().children().addClass('inactive');
                \$(this).removeClass('inactive');
                console.log();
                var panel=\$(this).attr('id').replace('_h_','_');
                \$('#${TASK}_panel div.wrapper div.display').html(\$('#'+ panel).html());
            });
        }
    </script>"    
fi



################################################################################

