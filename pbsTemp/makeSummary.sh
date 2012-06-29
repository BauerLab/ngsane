#!/bin/bash

HISEQINF=$1   # location of the HiSeqInf repository
CONFIG=$2


#PROGRAMS
. $HISEQINF/pbsTemp/header.sh
. $CONFIG


SUMMARYTMP="Summary.tmp"
SUMMARYFILE="Summary.html"


echo "Last modified "`date` >$SUMMARYTMP


if [ -n "$fastQC" ]; then
    echo "<h2>Read biases (FASTQC)</h2>">>$SUMMARYTMP
    $HISEQINF/bin/makeFastQCplot.sh $(pwd)/runStats/fastQC/ $(pwd)/runStats/ fastQCSummary.pdf $CONFIG > /dev/null #2>&1
    echo "<table><tr><td valign=top>" >>$SUMMARYTMP
    for f in $( ls runStats/fastQC/*.zip ); do
	n=`basename $f`
	n=${n/"_fastqc.zip"/}
	ICO="<img height=15px src=\"runStats/fastQC/"$n"_fastqc/Icons/"
        P=$(grep "PASS" -c runStats/fastQC/$n"_fastqc/summary.txt")
	W=$(grep "WARN" -c runStats/fastQC/$n"_fastqc/summary.txt")
	F=$(grep "FAIL" -c runStats/fastQC/$n"_fastqc/summary.txt")
	CHART=$ICO"tick.png\" title=\"$P\"\>"
	if [ "$W" -ne "0" ]; then CHART=$CHART""$ICO"warning.png\"\>"$W; fi
	if [ "$F" -ne "0" ]; then CHART=$CHART""$ICO"error.png\"\>"$F; fi
	echo "<a href=\"runStats/fastQC/"$n"_fastqc/fastqc_report.html\">$n.fastq</a>$CHART<br>" >>$SUMMARYTMP
    done
    echo "</td><td>">>$SUMMARYTMP
    echo "<img src=\"runStats/fastQCSummary.jpg\" alt=\"Quality scores for all reads\"/>" >>$SUMMARYTMP
    echo "</td></tr></table>">>$SUMMARYTMP
fi



if [ -n "$RUNMAPPINGBWA" ]; then
    LINKS=$LINKS" mapping"
    echo "<a name=\"mapping\"><h2>BWA Mapping</h2>">>$SUMMARYTMP
    echo "<pre>" >>$SUMMARYTMP
    $HISEQINF/pbsTemp/QC.sh $HISEQINF/pbsTemp/bwa.sh $QOUT/$TASKBWA >>$SUMMARYTMP
    for dir in ${DIR[@]}; do
	vali=$vali" $OUT/$dir/$TASKBWA"
    done
    echo "</pre><h3>Result</h3><pre>">>$SUMMARYTMP
    python $HISEQINF/bin/Summary.py "$vali" $ASD.bam.stats samstats >>$SUMMARYTMP
    echo "</pre>" >>$SUMMARYTMP
fi



if [ -n "$RUNREALRECAL" ]; then 
    LINKS=$LINKS" recal"
    echo "<a name=\"recal\"><h2>RECAL Mapping</h2>">>$SUMMARYTMP
    echo "<pre>" >>$SUMMARYTMP
    $HISEQINF/pbsTemp/QC.sh $HISEQINF/pbsTemp/reCalAln.sh $QOUT/$TASKRCA >>$SUMMARYTMP
    vali=""
    for dir in ${DIR[@]}; do
	vali=$vali" $OUT/$dir/$TASKRCA"
    done
    echo "</pre><h3>Result</h3><pre>">>$SUMMARYTMP
    python $HISEQINF/bin/Summary.py "$vali" $ASR".bam.stats" samstatsrecal >>$SUMMARYTMP
    echo "</pre>" >>$SUMMARYTMP
fi


if [ -n "$RUNTOPHATCUFF" ]; then
    vali=""
    LINKS=$LINKS" tophat"
    echo "<a name=\"tophat\"><h2>tophat Mapping</h2><br>Note, the duplication rate is not calculated by tophat and hence zero.">>$SUMMARYTMP
    echo "<pre>" >>$SUMMARYTMP
    $HISEQINF/pbsTemp/QC.sh $HISEQINF/pbsTemp/tophatcuff.sh $QOUT/$TASKTOPHAT >>$SUMMARYTMP
    for dir in ${DIR[@]}; do
	vali=$vali" $OUT/$dir/$TASKTOPHAT"
#	vali=$vali" "$(ls -d $OUT/$dir/$TASKTOPHAT/*/)
    done
    echo "</pre><h3>Result</h3><pre>">>$SUMMARYTMP
    python $HISEQINF/bin/Summary.py "$vali" bam.stats tophat >>$SUMMARYTMP
    echo "</pre>" >>$SUMMARYTMP
fi



if [ -n "$DEPTHOFCOVERAGE" ]; then
    LINKS=$LINKS" coverage"
    echo "<a name=\"coverage\"><h2>Coverage</h2>">>$SUMMARYTMP
    vali=""
    for dir in ${DIR[@]}; do
	vali=$vali" $OUT/$dir/$TASKDOC"
    done
    echo "<h3>Average coverage</h3><pre>">>$SUMMARYTMP
    python $HISEQINF/bin/Summary.py "$vali" $ASR".bam.doc.sample_summary" coverage >>$SUMMARYTMP
    echo "</pre><h3>Base pair coverage over all intervals</h3><pre>" >>$SUMMARYTMP
    python $HISEQINF/bin/Summary.py "$vali" $ASR".bam.doc.sample_cumulative_coverage_counts" coverage --p >>$SUMMARYTMP
    echo "</pre><h3>Intervals covered</h3><pre>" >>$SUMMARYTMP
    python $HISEQINF/bin/Summary.py "$vali" $ASR".bam.doc.sample_interval_statistics" coverage --p >>$SUMMARYTMP
    echo "</pre><h3>On Target</h3><pre>" >>$SUMMARYTMP
    python $HISEQINF/bin/Summary.py "$vali" $ASR".bam.stats" target >>$SUMMARYTMP
    echo "</pre>" >>$SUMMARYTMP
    
fi

if [ -n "$RUNVARCALLS" ]; then 
    LINKS=$LINKS" varcalls"
    echo "<a name=\"varcalls\"><h2>Variant calling</h2><pre>">>$SUMMARYTMP
    $HISEQINF/pbsTemp/QC.sh $HISEQINF/pbsTemp/gatkSNPs.sh $QOUT/$TASKVAR >> $SUMMARYTMP

    vali=""
    for dir in ${DIR[@]}; do
	vali=$vali" $OUT/$TASKVAR/$dir"
    done

    echo "</pre><h3>SNPs</h3><pre>">>$SUMMARYTMP
    python $HISEQINF/bin/Summary.py "$vali" "filter.snps.eval.txt" variant --n --l>>$SUMMARYTMP
    python $HISEQINF/bin/Summary.py "$vali" "recalfilt.snps.eval.txt" variant --n --l>>$SUMMARYTMP
    echo "</pre><h3>INDELs</h3><pre>" >>$SUMMARYTMP
    python $HISEQINF/bin/Summary.py "$vali" "filter.indel.eval.txt" variant --n --l>>$SUMMARYTMP
    echo "</pre>" >>$SUMMARYTMP
fi


if [ -n "$RUNANNOTATION" ]; then
    LINKS=$LINKS" annotation"
    echo "<a name=\"annotation\"><h2>Variant annotation</h2></a><pre>">>$SUMMARYTMP
    $HISEQINF/pbsTemp/QC.sh $HISEQINF/pbsTemp/annovar.sh $QOUT/$TASKANNOVAR >> $SUMMARYTMP

    vali=""
    for dir in ${DIR[@]}; do
	vali=$vali" "$( ls $OUT/$TASKANNOVAR/$dir/*.csv)
    done

    echo "</pre><h3>Annotation Files</h3>">>$SUMMARYTMP
    echo "Please right click the link and choose \"Save as...\" to download.<br><br>">> $SUMMARYTMP
    for v in $vali; do
	name=`basename $v`
	address=${v/\/illumina/http:\/\/cluster-vm.qbi.uq.edu.au}
	echo "<a href=\"$address\">$name</a><br>" >> $SUMMARYTMP
    done
    echo "<br>More information about the columns can be found on the <a target=new href=\"http://redmine.qbi.uq.edu.au/knowledgebase/articles/12\">Project Server</a> (uqlogin). and the description of the <a href=\"http://www.broadinstitute.org/gsa/wiki/index.php/Understanding_the_Unified_Genotyper%27s_VCF_files\">vcf file</a>">> $SUMMARYTMP


    echo "</pre>" >>$SUMMARYTMP

fi

#TODO add IGV

echo "<h2>Illumina Stuff</h2>">>$SUMMARYTMP
for f in BustardSummary.xml  RunInfo.xml  runParameters.xml  Summary.xml; do
echo "<a href=\"runStats/"$f"\">"$f"</a><br>" >>$SUMMARYTMP
done


echo "<h3>Quicklink</h3>" >$SUMMARYFILE.tmp
for i in $LINKS; do
    echo "<a href=#$i>$i</a> | ">>$SUMMARYFILE.tmp
done
echo "<br><br>" >>$SUMMARYFILE.tmp

cat $SUMMARYFILE.tmp  $SUMMARYTMP> $SUMMARYFILE

#echo "</body>" >>$SUMMARYFILE

rm $SUMMARYTMP
rm $SUMMARYFILE.tmp