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

# add bamannotate section
# $1=vali 
# $2=TASK (e.g.TASK_BWA)
# $3=output file ($SUMMARYTMP)
function bamAnnotate {
	echo "<h3 class='overall'>Reads overlapping annotated regions</h3>" >>$3
	python ${NGSANE_BASE}/core/Summary.py "$2" ${1} .anno.stats annostats >> $3
	BAMANNOUT=runStats/bamann/$(echo ${DIR[@]} | sed 's/ /_/g' | cut -c 1-60 )_${2}.ggplot
	BAMANNIMAGE=${BAMANNOUT/ggplot/pdf}
	if [ ! -f $BAMANNOUT ]; then mkdir -p $( dirname $BAMANNOUT); fi
	
	find ${1} -type f -name *anno.stats | xargs -d"\n" cat | head -n 1 | gawk '{print "type "$0" sample"}' > $BAMANNOUT
    for i in $(find ${1} -type f -name *anno.stats); do
        name=$(basename $i)
        arrIN=(${name//$ASD/ })
        grep --no-messages sum $i | gawk -v x=${arrIN[0]} '{print $0" "x}';
	done >> $BAMANNOUT
	sed -i -r 's/\s+/ /g' $BAMANNOUT
	Rscript ${NGSANE_BASE}/tools/bamann.R $BAMANNOUT $BAMANNIMAGE "Genome Features ${2}"
	convert $BAMANNIMAGE ${BAMANNIMAGE/pdf/jpg}
	echo "<h3>Annotation of mapped reads</h3>" >> $3
	echo "<div><a href=$PROJECT_RELPATH/$BAMANNIMAGE><img src=\""$PROJECT_RELPATH/${BAMANNIMAGE/.pdf/}"-0.jpg\" width='250px' style='float:left;'><img src=\""$PROJECT_RELPATH/${BAMANNIMAGE/.pdf/}"-1.jpg\" width='250px' style='float:left;'></a></div>">>$3

#	    python ${NGSANE_BASE}/tools/makeBamHistogram.py "${PROJECT}" $ROUTH >>$3

}

# source module requested for report generation (in order provided by the config)
for RUN_MODS in $(cat $CONFIG | egrep '^RUN.*=' | cut -d'=' -f 1); do 
    eval "if [ -n $RUN_MODS ]; then source ${NGSANE_BASE}/mods/run.d/${RUN_MODS/RUN}; fi"; 
done

################################################################################
# citation list
################################################################################    
PIPELINE="References"
PIPELINK="References"

LINKS="$LINKS $PIPELINK"
echo "<div class='panel'><div class='headbagb'><a name='"$PIPELINK"_panel'><h2 class='sub'>$PIPELINE</h2></a></div>" >>$SUMMARYTMP

echo '<h3>Please cite the following publications/resources</h3>' >>$SUMMARYTMP
echo '<table class="data"><thead><tr><th><div style="width:10px"></div></th><th><div></div></th></tr></thead><tbody>' >>$SUMMARYTMP

COUNT=1
while read -r line; do
    echo "<tr class='citation'><td class='top'>[$COUNT]</td><td class='left' >${line//\"/}</td></tr>" >>$SUMMARYTMP
    COUNT=$(( $COUNT + 1 ))
done <<< "$(cat $SUMMARYCITES | cut -d']' -f 2 | sort -u )"

echo "</tbody></table>" >>$SUMMARYTMP

echo "</div></div>">>$SUMMARYTMP
rm $SUMMARYCITES
    
################################################################################
echo '''
<!DOCTYPE html><html>
<head>
<meta charset="windows-1252">
<title>NGSANE project card</title>
<link rel="shortcut icon" href="includes/images/favicon.ico">
<link rel="icon" type="image/png" href="includes/images/favicon-196x196.png" sizes="196x196">
<link rel="icon" type="image/png" href="includes/images/favicon-160x160.png" sizes="160x160">
<link rel="icon" type="image/png" href="includes/images/favicon-96x96.png" sizes="96x96">
<link rel="icon" type="image/png" href="includes/images/favicon-16x16.png" sizes="16x16">
<link rel="icon" type="image/png" href="includes/images/favicon-32x32.png" sizes="32x32">
<meta name="msapplication-TileColor" content="#da532c">
<meta name="msapplication-TileImage" content="includes/images/mstile-144x144.png">
<meta name="msapplication-config" content="includes/images/browserconfig.xml">
<script type="text/javascript" src="http://ajax.googleapis.com/ajax/libs/jquery/1.11.1/jquery.min.js"></script>
<script type="text/javascript">
if (!window.jQuery) {
    document.write("<script src=includes/js/jquery.js><\/script>");
}
</script>
<script type="text/javascript" src="http://cdn.datatables.net/1.10.2/js/jquery.dataTables.min.js"></script>
<script type="text/javascript">
if (!$.fn.dataTableExt) { 
    document.write("<script src=includes/js/jquery.dataTables.min.js><\/script>");
}
</script>
<link rel="stylesheet" type="text/css" href="http://cdn.datatables.net/1.10.2/css/jquery.dataTables.css">
<style>''' >> $SUMMARYFILE.tmp
cat $NGSANE_BASE/core/includes/css/ngsane.css >> $SUMMARYFILE.tmp
echo '''
</style>
</head>
<body>
<script type="text/javascript">
//<![CDATA[
    var datatable_array=[];
''' >> $SUMMARYFILE.tmp
cat $NGSANE_BASE/core/includes/js/ngsane.js >> $SUMMARYFILE.tmp
echo '''
//]]>
</script>
<div id="center">
<div class='panel' id='quicklinks'><h2>Quicklinks</h2><div>
''' >> $SUMMARYFILE.tmp
declare -a LINKSET=( )
for i in $LINKS; do
    LINKSET=("${LINKSET[@]}" "<a href='#"$i"_panel' class='quicklinks'>$i</a>")
done
echo $(IFS='|' ; echo "${LINKSET[*]}") >> $SUMMARYFILE.tmp

echo '''
<div id ="right"><label>Aggregate: <input id="showAggregation" type="checkbox" /></label> <label>Global Filter: <input type="search" id="search" /></label></div>
</div><!-- Links -->
</div><!-- panel -->''' >>$SUMMARYFILE.tmp

echo "<hr><div class="footerline"><img style='float:left;padding-right:10px;' src='data:image/png;base64,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' /> Report generated with "`$NGSANE_BASE/bin/trigger.sh -v`"<br/>Last modified: "`date`"</div>" >> $SUMMARYTMP
echo "</div><!-- center --></body>" >> $SUMMARYTMP

echo '''
<script type='text/javascript'>
    var tables = [];
    $(document).ready(function() { 
        for (var i = 0; i < datatable_array.length; i++) { 
            eval("var "+datatable_array[i].html+" = "+datatable_array[i].json);
            eval("tables.push($(\"#"+datatable_array[i].html+"\").dataTable("+datatable_array[i].json+"));");
        }
        $("#search").keyup(function(){
            for (var i=0;i<tables.length;i++){
                tables[i].fnDraw();
            }
        });
        $("#showAggregation").click(function(){
            for (var i=0;i<tables.length;i++){
                tables[i].fnDraw();
            }
        });
    }); 
''' >> $SUMMARYTMP
cat $NGSANE_BASE/core/includes/js/genericJavaScript.js >> $SUMMARYTMP
echo '''
</script>
''' >> $SUMMARYTMP

#echo "<script type='text/javascript' src='includes/js/genericJavaScript.js'></script></html>" >> $SUMMARYTMP

# copy includes folder ommitting hidden files and folders
rsync -av --exclude=".*" $NGSANE_BASE/core/includes $(dirname $SUMMARYTMP)
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
