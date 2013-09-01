#!/bin/bash

echo "run"

NGSANE_BASE=$1   # location of the NGSANE repository
CONFIG=$2


#PROGRAMS
. $CONFIG
. ${NGSANE_BASE}/conf/header.sh
. $CONFIG

################################################################################
for MODULE in $MODULE_SUMMARY; do module load $MODULE; done  # save way to load modules that itself load other modules
for MODULE in $MODULE_R; do module load $MODULE; done  # save way to load modules that itself load other modules

################################################################################

SUMMARYTMP=$HTMLOUT".tmp"
SUMMARYFILE=$HTMLOUT".html"


echo "Last modified "`date` >$SUMMARYTMP

################################################################################
if [ -n "$RUNFASTQC" ]; then
    echo "fastqc"
    echo "<h2>Read biases (FASTQC)</h2>">>$SUMMARYTMP
#    echo ${NGSANE_BASE}/tools/makeFastQCplot.sh
#    if [ -n "$RUNFASTQC" ]; then
#	   ${NGSANE_BASE}/tools/makeFastQCplot.sh $(pwd)/runStats/$TASKFASTQC/ $(pwd)/runStats/ fastQCSummary.pdf $CONFIG > /dev/null #2>&1
#    fi
#    echo "done"
    echo "<table>" >>$SUMMARYTMP
    echo "<thead><th>Libary</th><th>Chart</th><th>Encoding</th><th>Library size</th><th>Read</th><th>Read length</th><th>%GC</th><th> Read Qualities</th><thead>" >>$SUMMARYTMP

    if [[ -e runStats/ && -e runStats/$TASKFASTQC/ ]]; then
        for f in $( ls runStats/$TASKFASTQC/*.zip ); do
            # get basename of f
            n=${f##*/}
            n=${n/"_fastqc.zip"/}
            ICO="<img height=15px src=\"runStats/$TASKFASTQC/"$n"_fastqc/Icons/"
            P=$(grep "PASS" -c runStats/$TASKFASTQC/$n"_fastqc/summary.txt")
            W=$(grep "WARN" -c runStats/$TASKFASTQC/$n"_fastqc/summary.txt")
            F=$(grep "FAIL" -c runStats/$TASKFASTQC/$n"_fastqc/summary.txt")
            CHART=$ICO"tick.png\" title=\"$P\"\>$P"
            if [ "$W" -ne "0" ]; then CHART=$CHART""$ICO"warning.png\"\>"$W; fi
            if [ "$F" -ne "0" ]; then CHART=$CHART""$ICO"error.png\"\>"$F; fi
            ENCODING=$(grep "Encoding" runStats/$TASKFASTQC/$n"_fastqc/fastqc_data.txt" | head -n 1 | cut -f 2)
            LIBRARYSIZE=$(grep "Total Sequences" runStats/$TASKFASTQC/$n"_fastqc/fastqc_data.txt" | head -n 1 | cut -f 2)
            READLENGTH=$(grep "Sequence length" runStats/$TASKFASTQC/$n"_fastqc/fastqc_data.txt" | head -n 1 | cut -f 2)
            GCCONTENT=$(grep "\%GC" runStats/$TASKFASTQC/$n"_fastqc/fastqc_data.txt" | head -n 1 | cut -f 2)
            if [[ "$f" == *$READTWO* ]] && [ "$f" != "${f/$READTWO/$READONE}" ]; then
                READ=2
            else
                READ=1
            fi
            echo "<tr><td><a href=\"runStats/$TASKFASTQC/"$n"_fastqc/fastqc_report.html\">$n.fastq</a></td><td>$CHART</td><td>$ENCODING</td><td>$LIBRARYSIZE</td><td>$READ</td><td>$READLENGTH</td><td>$GCCONTENT</td><td>" >>$SUMMARYTMP

            if [[ "$f" == *$READONE* ]]; then
                echo "<a href=\"runStats/$TASKFASTQC/${n}_fastqc/fastqc_report.html\"><img src=\"runStats/$TASKFASTQC/${n}_fastqc/Images/per_base_quality.png\" width=200 alt=\"Quality scores for all first reads\"/></a>" >>$SUMMARYTMP
                if [ -e ${f/$READONE/$READTWO} ] && [ "$f" != "${f/$READONE/$READTWO}" ]; then
                     echo "<a href=\"runStats/$TASKFASTQC/${n/$READONE/$READTWO}_fastqc/fastqc_report.html\"><img src=\"runStats/$TASKFASTQC/${n/$READONE/$READTWO}_fastqc/Images/per_base_quality.png\" width=200 alt=\"Quality scores for all second reads\"/></a>" >>$SUMMARYTMP
                fi
            
            elif [[ "$f" == *$READTWO* ]] && [ "$f" != "${f/$READTWO/$READONE}" ]; then
                echo "<a href=\"runStats/$TASKFASTQC/${n/$READTWO/$READONE}_fastqc/fastqc_report.html\"><img src=\"runStats/$TASKFASTQC/${n/$READTWO/$READONE}_fastqc/Images/per_base_quality.png\" width=200 alt=\"Quality scores for all first reads\"/></a>" >>$SUMMARYTMP
                echo "<a href=\"runStats/$TASKFASTQC/${n}_fastqc/fastqc_report.html\"><img src=\"runStats/$TASKFASTQC/${n}_fastqc/Images/per_base_quality.png\" width=200 alt=\"Quality scores for all second reads\"/></a>" >>$SUMMARYTMP
            fi
            echo "</td></tr><br>" >>$SUMMARYTMP
        done
    fi
    echo "</table>">>$SUMMARYTMP
fi


################################################################################
if [[ -n "$RUNMAPPINGBWA" || -n "$RUNMAPPINGBWA2" ]]; then
    echo "bwa"
    LINKS=$LINKS" mapping"
    echo "<a name=\"mapping\"><h2>BWA Mapping</h2>">>$SUMMARYTMP
    echo "<pre>" >>$SUMMARYTMP
    echo "QC"
    ${NGSANE_BASE}/mods/QC.sh ${NGSANE_BASE}/mods/bwa.sh $QOUT/$TASKBWA >>$SUMMARYTMP
    echo "gather dirs"
    for dir in ${DIR[@]}; do
	   vali=$vali" $OUT/$dir/$TASKBWA/"
    done
    echo "</pre><h3>Result</h3><pre>">>$SUMMARYTMP
    python ${NGSANE_BASE}/tools/Summary.py "$vali" $ASD.bam.stats samstats >>$SUMMARYTMP
    echo "</pre>" >>$SUMMARYTMP
    if [ -n "$RUNANNOTATINGBAM" ]; then
	echo "anno"
	echo "</pre><h3>Anno</h3><pre>">>$SUMMARYTMP
	python ${NGSANE_BASE}/tools/Summary.py "$vali" merg.anno.stats annostats >>$SUMMARYTMP
	echo "</pre>" >>$SUMMARYTMP
	ROUTH=runStats/$(echo ${DIR[@]}|sed 's/ /_/g')
	if [ ! -e $ROUTH ]; then mkdir $ROUTH; fi
	python ${NGSANE_BASE}/tools/makeBamHistogram.py "$vali" $ROUTH >>$SUMMARYTMP
    fi
fi


################################################################################
if [[ -n "$RUNREALRECAL" || -n "$RUNREALRECAL2" || -n "$RUNREALRECAL3" ]]; then 
    LINKS=$LINKS" recal"
    echo "<a name=\"recal\"><h2>RECAL Mapping</h2>">>$SUMMARYTMP
    echo "<pre>" >>$SUMMARYTMP
    ${NGSANE_BASE}/mods/QC.sh ${NGSANE_BASE}/mods/reCalAln.sh $QOUT/$TASKRCA >>$SUMMARYTMP
    vali=""
    for dir in ${DIR[@]}; do
	   vali=$vali" $OUT/$dir/$TASKRCA/"
    done
    echo "</pre><h3>Result</h3><pre>">>$SUMMARYTMP
    python ${NGSANE_BASE}/tools/Summary.py "$vali" $ASR".bam.stats" samstatsrecal >>$SUMMARYTMP
    echo "</pre>" >>$SUMMARYTMP
fi

################################################################################
if [[ -n "$RUNMAPPINGBOWTIE" ]]; then
    LINKS=$LINKS" mapping"
    echo "<a name=\"mapping\"><h2>BOWTIE v1 Mapping</h2>">>$SUMMARYTMP
    echo "<pre>" >>$SUMMARYTMP
    echo "QC"
    ${NGSANE_BASE}/mods/QC.sh ${NGSANE_BASE}/mods/bowtie.sh $QOUT/$TASKBOWTIE/ >>$SUMMARYTMP
    echo "gather dirs"
    vali=""
    for dir in ${DIR[@]}; do
        vali=$vali" $OUT/$dir/$TASKBOWTIE/"
    done
    echo "</pre><h3>Result</h3><pre>">>$SUMMARYTMP
    python ${NGSANE_BASE}/tools/Summary.py "$vali" $ASD.bam.stats samstats >>$SUMMARYTMP
    echo "</pre>" >>$SUMMARYTMP
fi

################################################################################
if [[ -n "$RUNFASTQSCREEN" ]]; then
    LINKS=$LINKS" screening"
    echo "<a name=\"screening\"><h2>Fastq screen</h2>">>$SUMMARYTMP
    echo "<pre>" >>$SUMMARYTMP
    echo "QC"
    ${NGSANE_BASE}/mods/QC.sh ${NGSANE_BASE}/mods/fastqscreen.sh $QOUT/$TASKFASTQSCREEN/ >>$SUMMARYTMP
    echo "gather dirs"
    vali=""
    for dir in ${DIR[@]}; do
        vali=$vali" $OUT/$dir/$TASKFASTQSCREEN/"
    done
    echo "</pre><h3>Result</h3><pre>">>$SUMMARYTMP
    python ${NGSANE_BASE}/tools/Summary.py "$vali" _screen.txt fastqscreen --noSummary --noOverallSummary >>$SUMMARYTMP

    row0=""
    row1=""
    for dir in $vali; do
        for f in $(ls $dir/*_screen.png); do
            n=${f##*/}
            n=${n/"_screen.png"/}

            row0+="<td>$n</td>"
            row1+="<td><a href=\"$dir/"$n"_screen.png\"><img src=\"$dir/"$n"_screen.png\" width=\"300px\"/></a></td>"
        done
    done
    echo "<table><tr>$row0</tr><tr>$row1</tr></table>" >> $SUMMARYTMP

    echo "</pre>" >>$SUMMARYTMP
fi

################################################################################
if [[ -n "$RUNMAPPINGBOWTIE2" ]]; then
    LINKS=$LINKS" mapping"
    echo "<a name=\"mapping\"><h2>BOWTIE v2 Mapping</h2>">>$SUMMARYTMP
    echo "<pre>" >>$SUMMARYTMP
    echo "QC"
    ${NGSANE_BASE}/mods/QC.sh ${NGSANE_BASE}/mods/bowtie2.sh $QOUT/$TASKBOWTIE2/ >>$SUMMARYTMP
    echo "gather dirs"
    vali=""
    for dir in ${DIR[@]}; do
        vali=$vali" $OUT/$dir/$TASKBOWTIE2/"
    done
    echo "</pre><h3>Result</h3><pre>">>$SUMMARYTMP
    python ${NGSANE_BASE}/tools/Summary.py "$vali" $ASD.bam.stats samstats >>$SUMMARYTMP
    echo "</pre>" >>$SUMMARYTMP
fi

################################################################################
if [[ -n "$RUNTOPHATCUFF" || -n "$RUNTOPHATCUFF2" ]]; then
    vali=""
    LINKS=$LINKS" tophat"
    echo "<a name=\"tophat\"><h2>tophat Mapping</h2><br>Note, the duplication rate is not calculated by tophat and hence zero.">>$SUMMARYTMP
    echo "<pre>" >>$SUMMARYTMP
    ${NGSANE_BASE}/mods/QC.sh ${NGSANE_BASE}/mods/tophatcuff.sh $QOUT/$TASKTOPHAT/ >>$SUMMARYTMP
    echo "</pre><h3>Result</h3><pre>">>$SUMMARYTMP
    CURDIR=$(pwd)
    for dir in ${DIR[@]}; do
    	vali=$vali" $OUT/$dir/$TASKTOPHAT/"
    #	vali=$vali" "$(ls -d $OUT/$dir/$TASKTOPHAT/*/)
    	cd $OUT/$dir/$TASKTOPHAT
    	for d in $(find . -maxdepth 1 -mindepth 1 -type d -exec basename '{}' \; | grep "RNASeQC"); do
            echo "<a href=\"$dir/$TASKTOPHAT/$d/index.html\">RNAseq-QC for $dir/$d</a><br/>" >> $CURDIR/$SUMMARYTMP
	done
    done
    cd $CURDIR
    python ${NGSANE_BASE}/tools/Summary.py "$vali" bam.stats tophat >>$SUMMARYTMP
    echo "</pre>" >>$SUMMARYTMP
fi


################################################################################
if [[ -n "$DEPTHOFCOVERAGE"  || -n "$DEPTHOFCOVERAGE2" ]]; then
    LINKS=$LINKS" coverage"
    echo "<a name=\"coverage\"><h2>Coverage</h2>">>$SUMMARYTMP
    vali=""
    for dir in ${DIR[@]}; do
	   vali=$vali" $OUT/$dir/$TASKDOC/"
    done
    echo "<h3>Average coverage</h3><pre>">>$SUMMARYTMP
    python ${NGSANE_BASE}/tools/Summary.py "$vali" $ASR".bam.doc.sample_summary" coverage >>$SUMMARYTMP
    echo "</pre><h3>Base pair coverage over all intervals</h3><pre>" >>$SUMMARYTMP
    python ${NGSANE_BASE}/tools/Summary.py "$vali" $ASR".bam.doc.sample_cumulative_coverage_counts" coverage --p >>$SUMMARYTMP
    echo "</pre><h3>Intervals covered</h3><pre>" >>$SUMMARYTMP
    python ${NGSANE_BASE}/tools/Summary.py "$vali" $ASR".bam.doc.sample_interval_statistics" coverage --p >>$SUMMARYTMP
    echo "</pre><h3>On Target</h3><pre>" >>$SUMMARYTMP
    python ${NGSANE_BASE}/tools/Summary.py "$vali" $ASR".bam.stats" target >>$SUMMARYTMP
    echo "</pre>" >>$SUMMARYTMP
    
fi

################################################################################
if [ -n "$RUNVARCALLS" ]; then 
    LINKS=$LINKS" varcalls"
    echo "<a name=\"varcalls\"><h2>Variant calling</h2><pre>">>$SUMMARYTMP
    ${NGSANE_BASE}/mods/QC.sh ${NGSANE_BASE}/mods/gatkSNPs.sh $QOUT/$TASKVAR >> $SUMMARYTMP

    vali=""
    for dir in ${DIR[@]}; do
	   vali=$vali" $OUT/$TASKVAR/$dir/"
    done
    echo "</pre><h3>SNPs</h3><pre>">>$SUMMARYTMP
    python ${NGSANE_BASE}/tools/Summary.py "$vali" "filter.snps.eval.txt" variant --n --l>>$SUMMARYTMP
    python ${NGSANE_BASE}/tools/Summary.py "$vali" "recalfilt.snps.eval.txt" variant --n --l>>$SUMMARYTMP
    echo "</pre><h3>INDELs</h3><pre>" >>$SUMMARYTMP
    python ${NGSANE_BASE}/tools/Summary.py "$vali" "filter.indel.eval.txt" variant --n --l>>$SUMMARYTMP
    echo "</pre>" >>$SUMMARYTMP
fi

################################################################################
if [ -n "$RUNANNOTATION" ]; then
    LINKS=$LINKS" annotation"
    echo "<a name=\"annotation\"><h2>Variant annotation</h2></a><pre>">>$SUMMARYTMP
    ${NGSANE_BASE}/mods/QC.sh ${NGSANE_BASE}/mods/annovar.sh $QOUT/$TASKANNOVAR >> $SUMMARYTMP
    
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

################################################################################
if [ -n "$RUNTRIMGALORE" ];then
    LINKS=$LINKS" trimgalore"
    echo "<a name=\"trimgalore\"><h2>Trim-Galore results</h2></a><pre>">>$SUMMARYTMP
    ${NGSANE_BASE}/mods/QC.sh ${NGSANE_BASE}/mods/trimgalore.sh $QOUT/$TASKTRIMGALORE >> $SUMMARYTMP

    echo "</pre><h3>trimgalore</h3><pre>">>$SUMMARYTMP
    vali=""
    for dir in ${DIR[@]}; do
        vali=$vali" $OUT/fastq/${dir/_$TASKTRIMGALORE/}_$TASKTRIMGALORE/"
    done
    python ${NGSANE_BASE}/tools/Summary.py "$vali" "_trimming_report.txt" trimgalore --noSummary >> $SUMMARYTMP
    echo "</pre>" >>$SUMMARYTMP
fi

################################################################################
if [ -n "$RUNTRIMMOMATIC" ];then
    LINKS=$LINKS" trimmomatic"
    echo "<a name=\"trimmomatic\"><h2>Trimmomatic results</h2></a><pre>">>$SUMMARYTMP
    ${NGSANE_BASE}/mods/QC.sh ${NGSANE_BASE}/mods/trimmomatic.sh $QOUT/$TASKTRIMMOMATIC >> $SUMMARYTMP

    echo "</pre><h3>trimmomatic</h3><pre>">>$SUMMARYTMP
    vali=""
    for dir in ${DIR[@]}; do
        vali=$vali" $OUT/fastq/${dir/_$TASKTRIMMOMATIC/}_$TASKTRIMMOMATIC/"
    done
    python ${NGSANE_BASE}/tools/Summary.py "$vali" ".log" trimmomatic --noSummary >> $SUMMARYTMP
    echo "</pre>" >>$SUMMARYTMP
fi

################################################################################
if [ -n "$RUNCUTADAPT" ];then
    LINKS=$LINKS" cutadapt"
    echo "<a name=\"cutadapt\"><h2>Cutadapt results</h2></a><pre>">>$SUMMARYTMP
    ${NGSANE_BASE}/mods/QC.sh ${NGSANE_BASE}/mods/cutadapt.sh $QOUT/$TASKCUTADAPT >> $SUMMARYTMP

    echo "</pre><h3>cutadapt</h3><pre>">>$SUMMARYTMP
    vali=""
    for dir in ${DIR[@]}; do
        vali=$vali" $OUT/fastq/${dir/_$TASKCUTADAPT/}_$TASKCUTADAPT/"
    done
    python ${NGSANE_BASE}/tools/Summary.py "$vali" ".stats" cutadapt --noSummary >> $SUMMARYTMP
    echo "</pre>" >>$SUMMARYTMP
fi

################################################################################
if [ -n "$RUNHICLIB" ];then
    LINKS=$LINKS" hiclib"
    echo "<a name=\"hiclib\"><h2>HiC results</h2></a><pre>">>$SUMMARYTMP
    ${NGSANE_BASE}/mods/QC.sh ${NGSANE_BASE}/mods/hiclibMapping.sh $QOUT/$TASKHICLIB >> $SUMMARYTMP

    echo "</pre><h3>hiclib</h3><pre>">>$SUMMARYTMP
    vali=""
    for dir in ${DIR[@]}; do
        vali=$vali" $OUT/$dir/$TASKHICLIB/"
    done
    python ${NGSANE_BASE}/tools/Summary.py "$vali" ".log" hiclibMapping >> $SUMMARYTMP
    echo "</pre>" >>$SUMMARYTMP
    for dir in $vali; do
        for pdf in $(ls -f $dir/*.pdf 2>/dev/null ); do
            echo "<a href='$pdf'>${pdf##*/}</a> " >> $SUMMARYTMP
        done
    done
fi

################################################################################
if [ -n "$RUNHICUP" ];then
    LINKS=$LINKS" hicup"
    echo "<a name=\"hicup\"><h2>HiC results</h2></a><pre>">>$SUMMARYTMP
    ${NGSANE_BASE}/mods/QC.sh ${NGSANE_BASE}/mods/hicup.sh $QOUT/$TASKHICUP >> $SUMMARYTMP

    vali=""
    for dir in ${DIR[@]}; do
	   vali=$vali" $OUT/$dir/$TASKHICUP/"
    done
    
    echo "</pre><h3>hicup</h3>">>$SUMMARYTMP
    echo "<h4>truncater</h4><pre>">>$SUMMARYTMP
    python ${NGSANE_BASE}/tools/Summary.py "$vali" "hicup_truncater_summary.txt" hicup --noSummary --noOverallSummary >> $SUMMARYTMP
    echo "</pre><h4>mapper</h4><pre>">>$SUMMARYTMP
    python ${NGSANE_BASE}/tools/Summary.py "$vali" "hicup_mapper_summary.txt" hicup --noSummary --noOverallSummary >> $SUMMARYTMP
    echo "</pre><h4>filter</h4><pre>">>$SUMMARYTMP
    python ${NGSANE_BASE}/tools/Summary.py "$vali" "hicup_filter_summary_results.txt" hicup --noSummary --noOverallSummary >> $SUMMARYTMP
    echo "</pre><h4>deduplicator</h4><pre>">>$SUMMARYTMP
    python ${NGSANE_BASE}/tools/Summary.py "$vali" "hicup_deduplicater_summary_results.txt" hicup --noSummary --noOverallSummary >> $SUMMARYTMP
    echo "</pre>" >>$SUMMARYTMP
    
    row0=""
    row1=""
    row2=""
    for f in $(ls runStats/$TASKHICUP/*_ditag_classification.png); do
    	n=${f##*/}
	    n=${n/"_ditag_classification.png"/}
	    
	    row0+="<td>$n</td>"
	    row1+="<td><a href=\"runStats/$TASKHICUP/"$n"_ditag_classification.png\"><img src=\"runStats/$TASKHICUP/"$n"_ditag_classification.png\" width=\"200px\"/></a></td>"
   	    row2+="<td><a href=\"runStats/$TASKHICUP/"$n"_uniques_cis-trans.png\"><img src=\"runStats/$TASKHICUP/"$n"_uniques_cis-trans.png\" width=\"200px\"/></a></td>"

    done
    echo "<table><tr>$row0</tr><tr>$row1</tr><tr>$row2</tr></table>" >> $SUMMARYTMP
fi

################################################################################
if [ -n "$RUNHOMERCHIPSEQ" ];then
    LINKS=$LINKS" Homer-ChIP-seq"
    echo "<a name=\"Homer-ChIP-seq\"><h2>Homer ChIP-seq results</h2></a><pre>">>$SUMMARYTMP
    ${NGSANE_BASE}/mods/QC.sh ${NGSANE_BASE}/mods/chipseqHomer.sh $QOUT/$TASKHOMERCHIPSEQ >> $SUMMARYTMP

    echo "</pre><h3>Homer ChIP-seq</h3><pre>">>$SUMMARYTMP
    vali=""
    for dir in ${DIR[@]}; do
        vali=$vali" $OUT/$dir/$TASKHOMERCHIPSEQ/"
    done
    python ${NGSANE_BASE}/tools/Summary.py "$vali" ".summary.txt" homerchipseq >> $SUMMARYTMP
    echo "</pre>" >>$SUMMARYTMP
fi

################################################################################
if [ -n "$RUNPEAKRANGER" ];then
    LINKS=$LINKS" peakranger"
    echo "<a name=\"peakranger\"><h2>Peakranger results</h2></a><pre>">>$SUMMARYTMP
    ${NGSANE_BASE}/mods/QC.sh ${NGSANE_BASE}/mods/peakranger.sh $QOUT/$TASKPEAKRANGER >> $SUMMARYTMP

    echo "</pre><h3>Peakranger</h3><pre>">>$SUMMARYTMP
    vali=""
    for dir in ${DIR[@]}; do
        vali=$vali" $OUT/$dir/$TASKPEAKRANGER/"
    done
    python ${NGSANE_BASE}/tools/Summary.py "$vali" ".summary.txt" peakranger >> $SUMMARYTMP
    echo "</pre>" >>$SUMMARYTMP
fi

################################################################################
if [ -n "$RUNMACS2" ];then
    LINKS=$LINKS" MACS2"
    echo "<a name=\"MACS2\"><h2>MACS2 results</h2></a><pre>">>$SUMMARYTMP
    ${NGSANE_BASE}/mods/QC.sh ${NGSANE_BASE}/mods/macs2.sh $QOUT/$TASKMACS2 >> $SUMMARYTMP

    echo "</pre><h3>MACS2</h3><pre>">>$SUMMARYTMP
    vali=""
    for dir in ${DIR[@]}; do
        vali=$vali" $OUT/$dir/$TASKMACS2/"
    done
    python ${NGSANE_BASE}/tools/Summary.py "$vali" ".summary.txt" macs2 >> $SUMMARYTMP
    echo "</pre>" >>$SUMMARYTMP

    row0=""
    row1=""
    row2=""
    for dir in $vali; do
    for f in $(ls $dir/*model-0.png); do
        n=${f##*/}
            n=${n/"_model-0.png"/}
            row0+="<td>$n</td>"
            row1+="<td><a href=\"${f/-0.png/.pdf/}\"><img src=\"$f\" width=\"200px\"/></a></td>"
            row2+="<td><a href=\"${f/-1.png/.pdf/}\"><img src=\"${f/model-0.png/model-1.png}\" width=\"200px\"/></a></td>"
    done
    done
    echo "<table><tr>$row0</tr><tr>$row1</tr><tr>$row2</tr></table>" >> $SUMMARYTMP

fi

################################################################################
if [ -n "$RUNMEMECHIP" ];then
    LINKS=$LINKS" meme-chip"
    echo "<a name=\"meme-chip\"><h2>MEME-chip results</h2></a><pre>">>$SUMMARYTMP
    ${NGSANE_BASE}/mods/QC.sh ${NGSANE_BASE}/mods/memechip.sh $QOUT/$TASKMEMECHIP >> $SUMMARYTMP

    echo "</pre><h3>MEME-chip</h3><pre>">>$SUMMARYTMP
    vali=""
    CURDIR=$(pwd)
    for dir in ${DIR[@]}; do
        vali=$vali" $OUT/$dir/$TASKMEMECHIP/"
        cd $OUT/$dir/$TASKMEMECHIP
        for d in $(find . -maxdepth 1 -mindepth 1 -type d -exec basename '{}' \; ); do
                echo "<a href=\"$dir/$TASKMEMECHIP/$d/index.html\">$dir/$d</a> " >> $CURDIR/$SUMMARYTMP
        done
    done
    cd $CURDIR
    python ${NGSANE_BASE}/tools/Summary.py "$vali" ".summary.txt" memechip >> $SUMMARYTMP
    echo "</pre>" >>$SUMMARYTMP
fi

################################################################################
# Old code ...
################################################################################
if [ -n "$RUNANNOTATINGBAM3" ]; then
    echo "ANNOTATION"

    for typ in bam; do # dups.bam
	echo $typ
	DIR=runStats/$typ"_Annotate"
	if [ ! -e $DIR ]; then mkdir $DIR; fi
	for i in $(ls */*/*$typ.merg.anno.stats); do echo $i $(head -n 2 $i | tail -1 ) ; done > $DIR/distribution$typ.txt
	gawk 'BEGIN{print "sample type feature number"}{
                split($1,arr,"[/.]"); print arr[3]" "arr[1]" genes "$3"\n" arr[3]" "arr[1]" rRNA "$4"\n" arr[3]" "arr[1]" tRNA "$5"\n" arr[3]" "arr[1]" lincRNA "$6"\n" arr[3]" "arr[1]" miRNA "$7"\n" arr[3]" "arr[1]" snoRNA "$8"\n" arr[3]" "arr[1]" snRNA "$9"\n" arr[3]" "arr[1]" miscRNA "$10"\n" arr[3]" "arr[1]" PolyA "$11"\n" arr[3]" "arr[1]" other "$12"\n" arr[3]" "arr[1]" HiSeq "$13"\n" arr[3]" "arr[1]" UCSC_rRNA "$14"\n" arr[3]" "arr[1]" SegDups "$15"\n" arr[3]" "arr[1]" unannotated "$16"\n" arr[3]" "arr[1]" unmapped "$17}' $DIR/distribution$typ.txt > $DIR/distribution$typ.ggplot

	RSCRIPT=$DIR/"distribution$typ.ggplot".R
	P=$(pwd)
	DESCRIPT=$(basename $P)
	IMAGE=$DIR/"distribution$type.pdf"
	echo 'library("ggplot2")' > $RSCRIPT
	echo 'library("reshape")' >>$RSCRIPT
	echo 'pdf(file = "'$IMAGE'")' >> $RSCRIPT
	echo 'distribution <- read.table("'$DIR'/distribution'$typ'.ggplot", header=T, quote="\"")' >> $RSCRIPT
	echo 'distribution$feature <- factor(distribution$feature, levels = c("genes","lincRNA" ,"miRNA","snoRNA","snRNA", "miscRNA", "rRNA", "UCSC_rRNA", "tRNA", "PolyA", "other", "HiSeq", "SegDups", "unannotated", "unmapped"))' >> $RSCRIPT
	echo 'ggplot(distribution, aes(x = sample, y=number)) + geom_bar(aes(fill = feature), position = "fill") + scale_y_continuous("fraction") + opts(axis.text.x=theme_text(angle=-90, hjust=0),title = expression("'$DESCRIPT'")) + facet_grid(. ~ type , space = "free", scales = "free_x") ' >> $RSCRIPT
	echo 'ggplot(distribution, aes(x = sample, y=number)) + geom_bar(aes(fill = feature)) + opts(axis.text.x=theme_text(angle=-90, hjust=0),title = expression("'$DESCRIPT'")) + ylim(0,9e+07) +facet_grid(. ~ type , space = "free", scales = "free_x")' >> $RSCRIPT
	echo "dev.off()"  >> $RSCRIPT

	Rscript --vanilla $RSCRIPT
	convert $IMAGE ${IMAGE/pdf/jpg}
	echo "</pre><h3>Annotation of mapped reads</h3><pre>">>$SUMMARYTMP
	echo "<table><tr><td><a href=$IMAGE><img src=\""${IMAGE/.pdf/}"-0.jpg\"></a></td>">>$SUMMARYTMP
	echo "<td><a href=$IMAGE><img src=\""${IMAGE/.pdf/}"-1.jpg""\"></a></td></tr></table>">>$SUMMARYTMP

done

fi




#TODO add IGV

#echo "<h2>Illumina Stuff</h2>">>$SUMMARYTMP
#for f in BustardSummary.xml  RunInfo.xml  runParameters.xml  Summary.xml; do
#echo "<a href=\"runStats/"$f"\">"$f"</a><br>" >>$SUMMARYTMP
#done


################################################################################
echo "<h3>Quicklink</h3>" >$SUMMARYFILE.tmp
for i in $LINKS; do
    echo "<a href=#$i>$i</a> | ">>$SUMMARYFILE.tmp
done
echo "<br><br>" >>$SUMMARYFILE.tmp

cat $SUMMARYFILE.tmp  $SUMMARYTMP > $SUMMARYFILE

#echo "</body>" >>$SUMMARYFILE

rm $SUMMARYTMP
rm $SUMMARYFILE.tmp

################################################################################
# convert html to pdf
#if [ "$(hash prince)" == "" ]; then
#    prince $SUMMARYFILE -o ${HTMLOUT}.pdf
#fi

