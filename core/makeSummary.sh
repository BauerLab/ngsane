#!/bin/bash

# author: Fabian Buske
# date: Aug 2014

echo ">>>>> Generate HTML report "
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> $(basename $0) $*"

function usage {
echo -e "usage: $(basename $0) -k CONFIG

Generates the html report

required:
  -k | --toolkit <path>     location of the NGSANE repository 
"
exit
}

if [ ! $# -gt 1 ]; then usage ; fi

#INPUTS
while [ "$1" != "" ]; do
    case $1 in
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

################################################################################
hash module 2>/dev/null && for MODULE in $MODULE_SUMMARY; do module load $MODULE; done  # save way to load modules that itself load other modules
hash module 2>/dev/null && for MODULE in $MODULE_R; do module load $MODULE; done  # save way to load modules that itself load other modules


echo -e "--R           --\n "$(R --version | head -n 3)
[ -z "$(which R)" ] && echo "[ERROR] no R detected" && exit 1
echo -e "--Python      --\n "$(python --version 2>&1 | tee | head -n 1)
[ -z "$(which python)" ] && echo "[ERROR] no Python detected" && exit 1

################################################################################

SUMMARYTMP=$HTMLOUT".tmp"
SUMMARYFILE=$HTMLOUT".html"
SUMMARYCITES=$HTMLOUT".cites"
NGSANE_COMPILE_REPORT="TRUE"

mkdir -p $(dirname $SUMMARYTMP) && cat /dev/null > $SUMMARYTMP && cat /dev/null > $SUMMARYCITES # clean temporary content

PROJECT_RELPATH=$(python -c "import os.path; print os.path.relpath('$(pwd -P)',os.path.realpath('$(dirname $SUMMARYTMP)'))")
[ -z "$PROJECT_RELPATH" ] && PROJECT_RELPATH="."

# source module requested for report generation (in order provided by the config)
for RUN_MODS in $(cat $CONFIG | egrep '^RUN.*=' | cut -d'=' -f 1); do 
    eval "if [ -n $RUN_MODS ]; then source ${NGSANE_BASE}/mods/run.d/${RUN_MODS/RUN}; fi"; 
done
    
################################################################################
cat >> $SUMMARYFILE.tmp <<EOF
<!DOCTYPE html><html>
<head>
<meta charset="windows-1252">
<title>NGSANE project card</title>
<link rel="shortcut icon" href="data:image/png;base64,$(cat $NGSANE_BASE/core/includes/images/favicon.ico | base64)">
<link rel="icon" type="image/png" href="data:image/png;base64,$(cat $NGSANE_BASE/core/includes/images/favicon-196x196.png | base64)" sizes="196x196">
<link rel="icon" type="image/png" href="data:image/png;base64,$(cat $NGSANE_BASE/core/includes/images/favicon-160x160.png | base64)" sizes="160x160">
<link rel="icon" type="image/png" href="data:image/png;base64,$(cat $NGSANE_BASE/core/includes/images/favicon-96x96.png | base64)" sizes="96x96">
<link rel="icon" type="image/png" href="data:image/png;base64,$(cat $NGSANE_BASE/core/includes/images/favicon-16x16.png | base64)" sizes="16x16">
<link rel="icon" type="image/png" href="data:image/png;base64,$(cat $NGSANE_BASE/core/includes/images/favicon-32x32.png | base64)" sizes="32x32">
<meta name="msapplication-TileColor" content="#da532c">
<meta name="msapplication-TileImage" content="data:text/xml;base64,$(cat $NGSANE_BASE/core/includes/images/mstile-144x144.png | base64)">
<meta name="msapplication-config" content="data:image/png;base64,$(cat $NGSANE_BASE/core/includes/images/browserconfig.xml | base64)">
<script type="text/javascript" src="data:text/javascript;base64,$(cat $NGSANE_BASE/core/includes/js/jquery.js | base64)"></script>
<script type="text/javascript" src="data:text/javascript;base64,$(cat $NGSANE_BASE/core/includes/js/jquery.dataTables.min.js | base64)"></script>
<link rel="stylesheet" type="text/css" href="data:text/css;base64,$(cat $NGSANE_BASE/core/includes/css/jquery.dataTables.min.css | base64)">
<link rel="stylesheet" type="text/css" href="data:text/css;base64,$(cat $NGSANE_BASE/core/includes/css/ngsane.css | base64)">
</head>
<body>
<script type="text/javascript">
//<![CDATA[
    var datatable_array=[];
    var quicklink_array=[];
    $(cat $NGSANE_BASE/core/includes/js/ngsane.js)
//]]>
</script>
<div id="center">
<div class='panel' id='quicklinks'><div style="display:table; width: 100%"><h2>Quicklinks</h2><div id='quicklink' style="float:left"></div>
<div id="right"><label>Aggregate: <input id="showAggregation" type="checkbox" /></label> <label>Global Filter: <input type="search" id="search" /></label></div>
</div></div><!-- panel -->
EOF



################################################################################
# citation list
################################################################################    

cat >> $SUMMARYTMP <<EOF
<div class='panel'><div class='headbagb'><a name='References_panel'><h2 class='sub'>References</h2></a></div>
<h3>Please cite the following publications/resources</h3>
<table class="data"><thead><tr><th><div style="width:10px"></div></th><th><div></div></th></tr></thead><tbody>
EOF

COUNT=1
while read -r line; do
    echo "<tr class='citation'><td class='top'>[$COUNT]</td><td class='left' >${line//\"/}</td></tr>" >>$SUMMARYTMP
    COUNT=$(( $COUNT + 1 ))
done <<< "$(cat $SUMMARYCITES | cut -d']' -f 2 | sort -u )"

cat >> $SUMMARYTMP <<EOF
</tbody></table>
</div></div>

<hr><div class="footerline"><img style="float:left;padding-right:10px;" src="data:image/png;base64,$(cat $NGSANE_BASE/core/includes/images/favicon-32x32.png | base64)" /> Report generated with $($NGSANE_BASE/bin/trigger.sh -v)<br/>Last modified: $(date)</div>
</div><!-- center --></body>
<script type="text/javascript">
    var tables = [];
    \$(document).ready(function() { 
        for (var i = 0; i < datatable_array.length; i++) { 
            eval("var "+datatable_array[i].html+" = "+datatable_array[i].json);
            eval("tables.push(\$(\"#"+datatable_array[i].html+"\").dataTable("+datatable_array[i].json+"));");
        }
        \$("#search").keyup(function(){
            for (var i=0;i<tables.length;i++){
                tables[i].fnDraw();
            }
        });
        \$("#showAggregation").click(function(){
            for (var i=0;i<tables.length;i++){
                tables[i].fnDraw();
            }
        });

        var quicklinks_str=""
        for (var i = 0; i < quicklink_array.length; i++) {
            quicklinks_str = quicklinks_str+"<a href='#"+ quicklink_array[i]+"_panel' class='quicklinks'>"+ capitaliseFirstLetter(quicklink_array[i])+"</a> | "
        }
        // add references
        quicklinks_str = quicklinks_str+"<a href='#References_panel' class='quicklinks'>References</a>"
        \$("#quicklink").html(quicklinks_str);
    }); 
    $(cat $NGSANE_BASE/core/includes/js/genericJavaScript.js)
</script>
EOF
rm $SUMMARYCITES

################################################################################
cat $SUMMARYFILE.tmp $SUMMARYTMP > $SUMMARYFILE

rm $SUMMARYTMP
rm $SUMMARYFILE.tmp

################################################################################
# make tar containing all files smaller than 80k
if [ -n "$SUMMARYTAR" ];then

    [ -f Summary_files.tmp ] && rm Summary_files.tmp
    find . -size -$SUMMARYTAR"k" > Summary_files.tmp
    echo "$SUMMARYFILE" >> Summary_files.tmp
    tar -czf ${SUMMARYFILE%.*}.tar.gz -T Summary_files.tmp --no-recursion
    rm Summary_files.tmp
fi

################################################################################
echo ">>>>> Generate HTML report - FINISHED"
echo ">>>>> enddate "`date`
