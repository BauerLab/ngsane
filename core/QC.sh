#!/bin/bash

# author: Denis C. Bauer
# date: Feb.2011
# modified by Fabian Buske, Aug 2013

function usage {
echo -e "usage: $(basename $0) [OPTIONS] -m MOD.SCRIPT -l QOUT -t TASK -k CONFIG"
exit
}

if [ ! $# -gt 2 ]; then usage ; fi

HTMLOUTPUT=""
RESULTSUFFIX=""

while [ "$1" != "" ]; do
    case $1 in
        -m | --modscript )      shift; SCRIPT=$1 ;; # location of mod script (comma-sep if postcommand)
        -l | --log )            shift; QOUT=$1 ;;   # location of log output
        -t | --task )           shift; TASK=$1 ;;   # task at hand
        -r | --results-dir )    shift; OUTDIR=$1 ;; # location of the output 
        -s | --results-task )   shift; OUTTASK=$1 ;; # task the results are put in 
        -f | --filesuffix )     shift; RESULTSUFFIX=$1 ;; # suffix of the result file
        -o | --html-file )      shift; HTMLOUTPUT=$1;; # where the output will be place in the end
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the NGSANE repository
        -h | --help )           usage ;;
        * )                     echo "don't understand "$1
    esac
    shift
done

#PROGRAMS
. $CONFIG
. ${NGSANE_BASE}/conf/header.sh
. $CONFIG

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
    echo "<div class='tabContent_hide' id='DC_${TASK}_checklist'><div><pre>"
fi


# pack normal and post command in a long string for the subsequent loops
# file1*file2*file3:task.sh postcommand:posttask.sh
for dir in ${DIR[@]}; do
    for LOGFILE in $(find $QOUT/$TASK/ -maxdepth 1 -type f -name "${dir%%/*}*.out" 2>/dev/null ); do
        LOGFILES="${LOGFILES} ${LOGFILE}"
    done
done
# remove potential duplicates
LOGFILES=$(echo "${LOGFILES}" | tr ' ' '\n' | sort -u | tr '\n' ' ')
echo $LOGFILES

if [ ! -z "$LOGFILES" ]; then # if there is really more than just the postcommand
	SCRIPTFILES=$(echo $LOGFILES | tr ' ' '*')":"${NGSANE_BASE}/mods/${SCRIPT/%,*/}
fi
# normal postcommand
if [ -e $QOUT/$TASK/postcommand.out ]; then
	SCRIPTFILES=$SCRIPTFILES" "$QOUT/$TASK/postcommand.out":"${NGSANE_BASE}/mods/${SCRIPT/*,/}
fi

#echo $SCRIPTFILES 1>&2


#########################################################################################
# FIRST TAB
#########################################################################################

for T in $SCRIPTFILES; do

	# unpack script name and files
	SCRIPT=${T/*:/}
	files=$(echo ${T/%:*/} | tr "*" " "files=${FILES/*/ })
	nrfiles=$(echo $files | wc -w | sed -e 's/^[ \t]*//')

	echo "###################################################"
	echo "# NGSANE ${SCRIPT/.sh/} "
	echo "###################################################"

	CHECKPOINTS_PASSED=0
	CHECKPOINTS_FAILED=0
	
	echo ">>>>>>>>>> Finished"
	finished=$(egrep "^>{5} .* FINISHED" $(echo $files) | cut -d ":" -f 1 | sort -u | wc -l | sed -e 's/^[ \t]*//')
	if [ "$finished" = "$nrfiles" ]; then
	    echo "QC_PASS .. $finished/$nrfiles are finished"
	    CHECKPOINTS_PASSED=`expr $CHECKPOINTS_PASSED + 1`
	else
	    echo "**_FAIL .. $finished/$nrfiles are finished"
	    CHECKPOINTS_FAILED=`expr $CHECKPOINTS_FAILED + 1`
	fi

	echo ">>>>>>>>>> Errors"
	ERROR=$(egrep QCVARIABLES $SCRIPT | tr ' ' '_' | tr ',' ' ')
	ERROR=${ERROR/"# QCVARIABLES"/}
	for i in $ERROR; do
	  i=${i//_/ }
	  var=$(grep -i "$i" $(echo $files) | cut -d ":" -f 1 | sort -u | wc -l | sed -e 's/^[ \t]*//')
	  if [ "$var" = "0" ]; then
	    echo "QC_PASS .. $var/$nrfiles have $i"
	    CHECKPOINTS_PASSED=`expr $CHECKPOINTS_PASSED + 1`
	  else
	    echo "**_FAIL .. $var/$nrfiles have $i"
	    CHECKPOINTS_FAILED=`expr $CHECKPOINTS_FAILED + 1`
	  fi
	done

	echo ">>>>>>>>>> CheckPoints "
	PROGRESS=$(grep '^NGSANE_CHECKPOINT_INIT "' $SCRIPT | awk -F'"' '{print $2}' | tr ' ' '_')
	for i in $PROGRESS; do
	  i=${i//_/ }
	  var=$(egrep "^\*{9} $i" $files | cut -d ":" -f 1 | sort -u | wc -l | sed -e 's/^[ \t]*//')
	  if [ ! "$var" = "$nrfiles" ]; then
	    echo "**_FAIL .. $var/$nrfiles have $i"
	    CHECKPOINTS_FAILED=`expr $CHECKPOINTS_FAILED + 1`
	  else
	    echo "QC_PASS .. $var/$nrfiles have $i"
	    CHECKPOINTS_PASSED=`expr $CHECKPOINTS_PASSED + 1`
	  fi
	done

	echo ""

done

#########################################################################################
# Second TAB
#########################################################################################


if [ -n "$HTMLOUTPUT" ]; then
    echo "</pre></div></div>"
    echo "<div class='tabContent_hide' id='DC_${TASK}_notes'><div class='box scroll'><pre>"
else
    echo ">>>>>>>>>> Notes"
fi

SUMNOTES=0
#for i in $( ls $QOUT/$TASK/*.out ) ; do
for i in $LOGFILES ;do
    echo -e "\n${i/$LOGFOLDER\//}"
    NOTELIST=$(egrep "^\[NOTE\]" $i)
    echo -e "$NOTELIST"
    SUMNOTES=`expr $SUMNOTES + $(echo -e "$NOTELIST" | awk 'BEGIN{count=0} NF != 0 {++count} END {print count}' )`
done


#########################################################################################
# Third TAB
#########################################################################################

if [ -n "$HTMLOUTPUT" ]; then
    echo "</pre></div></div>"
    echo "<div class='tabContent_hide' id='DC_${TASK}_errors'><div class='box scroll'><pre>"
else
    echo ">>>>>>>>>> Errors"
fi

SUMERRORS=0
for i in $LOGFILES ;do
    echo -e "\n${i/$LOGFOLDER\//}"
    ERRORLIST=$(egrep "^\[ERROR\]" $i)
    if [ -n "$ERRORLIST" ]; then 
        echo -e "$ERRORLIST"
    else
        echo "-- all good, no errors"
    fi
    SUMERRORS=$(expr $SUMERRORS + $(echo -e "$ERRORLIST" | awk 'BEGIN{count=0} NF != 0 {++count} END {print count}' ))
done


#########################################################################################
# Third TAB
#########################################################################################
if [ -n "$HTMLOUTPUT" ]; then

    echo "</pre></div></div>"
    echo "<div class='tabContent_hide' id='DC_${TASK}_logfiles'><div><div class='box scroll'>"
    for i in $LOGFILES ;do
        FN=$(python -c "import os.path; print os.path.relpath(os.path.realpath('$i'),os.path.realpath('$(dirname $HTMLOUTPUT)'))")
        echo "<a href='$FN'>${i/$QOUT\/$TASK\//}</a><br/>"
    done
    echo "</div></div></div>"
    if [ -n "$RESULTSUFFIX" ]; then
        echo "<div class='tabContent_hide' id='DC_${TASK}_PRIMARY_RESULT_FILES'><div><div class='box scroll'>"
        for dir in ${DIR[@]}; do
            if [ ! -d $OUTDIR/${dir%%/*}/$OUTTASK/ ]; then
                echo "[NOTE] No result detected: $dir" 1>&2
            else
                for i in $(find $OUTDIR/${dir%%/*}/$OUTTASK/ -maxdepth 2 -type f -name *$RESULTSUFFIX | sort -n ); do
                    FN=$(python -c "import os.path; print os.path.relpath(os.path.realpath('$i'),os.path.realpath('$(dirname $HTMLOUTPUT)'))")
                    echo "<a href='$FN'>${i/$OUTDIR\/*\/$OUTTASK\//}</a><br/>"
                done
            fi
        done
        echo "</div></div></div>"
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
