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
<head><meta charset="windows-1252"><title>NGSANE project card</title>
<script type="text/javascript" src="includes/js/jquery.js"></script>
<script type="text/javascript" charset="utf8" src="includes/js/jquery.dataTables.min.js"></script>
<script type="text/javascript" charset="utf8" src="includes/js/dataTables.responsive.js"></script>
<link rel="stylesheet" type="text/css" href="includes/css/jquery.dataTables.min.css">
<link rel="stylesheet" type="text/css" href="includes/css/dataTables.responsive.css">
<link rel="stylesheet" type="text/css" href="includes/css/ngsane.css">
<script type='text/javascript' src='includes/js/ngsane.js'></script>
<link rel="shortcut icon" href="includes/images/favicon.ico">
<link rel="apple-touch-icon" sizes="57x57" href="includes/images/apple-touch-icon-57x57.png">
<link rel="apple-touch-icon" sizes="114x114" href="includes/images/apple-touch-icon-114x114.png">
<link rel="apple-touch-icon" sizes="72x72" href="includes/images/apple-touch-icon-72x72.png">
<link rel="apple-touch-icon" sizes="144x144" href="includes/images/apple-touch-icon-144x144.png">
<link rel="apple-touch-icon" sizes="60x60" href="includes/images/apple-touch-icon-60x60.png">
<link rel="apple-touch-icon" sizes="120x120" href="includes/images/apple-touch-icon-120x120.png">
<link rel="apple-touch-icon" sizes="76x76" href="includes/images/apple-touch-icon-76x76.png">
<link rel="apple-touch-icon" sizes="152x152" href="includes/images/apple-touch-icon-152x152.png">
<link rel="icon" type="image/png" href="includes/images/favicon-196x196.png" sizes="196x196">
<link rel="icon" type="image/png" href="includes/images/favicon-160x160.png" sizes="160x160">
<link rel="icon" type="image/png" href="includes/images/favicon-96x96.png" sizes="96x96">
<link rel="icon" type="image/png" href="includes/images/favicon-16x16.png" sizes="16x16">
<link rel="icon" type="image/png" href="includes/images/favicon-32x32.png" sizes="32x32">
<meta name="msapplication-TileColor" content="#da532c">
<meta name="msapplication-TileImage" content="includes/images/mstile-144x144.png">
<meta name="msapplication-config" content="includes/images/browserconfig.xml">
</head><body>
<script type="text/javascript">
//<![CDATA[
    var datatable_array=[];
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
<div id ="right"><label>Global Filter: <input type="search" id="search" /></label></div>
</div><!-- Links -->
</div><!-- panel -->''' >>$SUMMARYFILE.tmp

echo "<hr><span>Report generated with "`$NGSANE_BASE/bin/trigger.sh -v`"</span><span style='float:right;'>Last modified: "`date`"</span>" >> $SUMMARYTMP
echo "</div><!-- center --></body>" >> $SUMMARYTMP


echo '''
<script type='text/javascript'>
    $(document).ready(function() { 
        var tables = [];
        for (var i = 0; i < datatable_array.length; i++) { 
            eval("var "+datatable_array[i]["html"]+" = "+datatable_array[i]["json"]);
            eval("$(\"#"+datatable_array[i]["html"]+"\").dataTable("+datatable_array[i]["json"]+")");
            eval("tables.push("+datatable_array[i]["html"]+")");
        }
        $("#search").keyup(function(){
            for (var i=0;i<tables.length;i++){
                tables[i].fnDraw();
            }
        });
    }); 
</script>
''' >> $SUMMARYTMP

#for i in $LINKS; do
#    if [ "$i" != "References" ]; then 
#	    echo "if (typeof ${i}_Table_json !== 'undefined'){var ${i}_Table = \$("\'"#${i}_id_Table"\'").dataTable(${i}_Table_json);}" >> $SUMMARYTMP
#     fi
#done
#echo "\$('#search').keyup(function(){" >> $SUMMARYTMP
#for i in $LINKS; do
#	if [ "$i" != "References" ]; then
#		echo "    if (typeof ${i}_Table !== 'undefined'){${i}_Table.fnDraw();}" >> $SUMMARYTMP
#	fi
#done
#echo "});} ); </script>" >> $SUMMARYTMP


echo "<script type='text/javascript' src='includes/js/genericJavaScript.js'></script></html>" >> $SUMMARYTMP

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
