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
    PIPELINE="FASTQC"
    PIPELINK="fastqc"
    
    LINKS=$LINKS" $PIPELINK"
    echo "<div class='panel'><div class='headbagb'><a name='$PIPELINK'><h2 class='sub'>Read biases (FASTQC) </h2></a></div>" >>$SUMMARYTMP

    echo "<table class="data">" >>$SUMMARYTMP
    echo "<thead><tr><th>Libary</th><th>Chart</th><th>Encoding</th><th>Library size</th><th>Read</th><th>Read length</th><th>%GC</th><th> Read Qualities</th></tr><thead><tbody>" >>$SUMMARYTMP

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
    echo "</tbody></table>">>$SUMMARYTMP

    echo "</div></div>">>$SUMMARYTMP
fi

################################################################################
if [[ -n "$RUNFASTQSCREEN" ]]; then
    PIPELINE="FASTQ screen"
    PIPELINK="fastqscreen"
    
    LINKS=$LINKS" $PIPELINK"
    echo "<div class='panel'><div class='headbagb'><a name='$PIPELINK'><h2 class='sub'>$PIPELINE</h2></a></div><div class='wrapper'><div class='results'>" >>$SUMMARYTMP

    echo "QC"
    ${NGSANE_BASE}/mods/QC.sh -o -m ${NGSANE_BASE}/mods/fastqscreen.sh -l $QOUT/$TASKFASTQSCREEN/ >>$SUMMARYTMP
    echo "gather dirs"
    vali=""
    for dir in ${DIR[@]}; do
        vali=$vali" $OUT/$dir/$TASKFASTQSCREEN/"
    done
    echo "<h3>Result</h3>">>$SUMMARYTMP
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

    echo "</div></div></div>" >>$SUMMARYTMP
fi


################################################################################
if [[ -n "$RUNMAPPINGBWA" || -n "$RUNMAPPINGBWA2" ]]; then
    PIPELINE="BWA mapping"
    PIPELINK="bwa"
   
    LINKS=$LINKS" $PIPELINK"
    echo "<div class='panel'><div class='headbagb'><a name='$PIPELINK'><h2 class='sub'>$PIPELINE</h2></a></div><div class='wrapper'><div class='results'>" >>$SUMMARYTMP

    echo "QC"
    ${NGSANE_BASE}/mods/QC.sh -o -m ${NGSANE_BASE}/mods/bwa.sh -l $QOUT/$TASKBWA >>$SUMMARYTMP
    echo "gather dirs"
    for dir in ${DIR[@]}; do
	   vali=$vali" $OUT/$dir/$TASKBWA/"
    done
    echo "<h3>Result</h3>">>$SUMMARYTMP
    python ${NGSANE_BASE}/tools/Summary.py "$vali" $ASD.bam.stats samstats >>$SUMMARYTMP

    if [ -n "$RUNANNOTATINGBAM" ]; then
    	echo "<h3>Annotation</h3>" >>$SUMMARYTMP
    	python ${NGSANE_BASE}/tools/Summary.py "$vali" merg.anno.stats annostats >>$SUMMARYTMP
    	ROUTH=runStats/$(echo ${DIR[@]}|sed 's/ /_/g')
    	if [ ! -e $ROUTH ]; then mkdir $ROUTH; fi
	   python ${NGSANE_BASE}/tools/makeBamHistogram.py "$vali" $ROUTH >>$SUMMARYTMP
    fi
    
    echo "</div></div></div>" >>$SUMMARYTMP
fi


################################################################################
if [[ -n "$RUNREALRECAL" || -n "$RUNREALRECAL2" || -n "$RUNREALRECAL3" ]]; then 
    PIPELINE="RECAL mapping"
    PIPELINK="recal"
   
    LINKS=$LINKS" $PIPELINK"
    echo "<div class='panel'><div class='headbagb'><a name='$PIPELINK'><h2 class='sub'>$PIPELINE</h2></a></div><div class='wrapper'><div class='results'>" >>$SUMMARYTMP

    ${NGSANE_BASE}/mods/QC.sh -o -m ${NGSANE_BASE}/mods/reCalAln.sh -l $QOUT/$TASKRCA >>$SUMMARYTMP
    vali=""
    for dir in ${DIR[@]}; do
	   vali=$vali" $OUT/$dir/$TASKRCA/"
    done
    echo "<h3>Result</h3>">>$SUMMARYTMP
    python ${NGSANE_BASE}/tools/Summary.py "$vali" $ASR".bam.stats" samstatsrecal >>$SUMMARYTMP

    echo "</div></div></div>" >>$SUMMARYTMP
fi

################################################################################
if [[ -n "$RUNMAPPINGBOWTIE" ]]; then
    PIPELINE="BOWTIE v1 mapping"
    PIPELINK="bowtie"
   
    LINKS=$LINKS" $PIPELINK"
    echo "<div class='panel'><div class='headbagb'><a name='$PIPELINK'><h2 class='sub'>$PIPELINE</h2></a></div><div class='wrapper'><div class='results'>" >>$SUMMARYTMP

    echo "QC"
    ${NGSANE_BASE}/mods/QC.sh -o -m ${NGSANE_BASE}/mods/bowtie.sh -l $QOUT/$TASKBOWTIE/ >>$SUMMARYTMP
    echo "gather dirs"
    vali=""
    for dir in ${DIR[@]}; do
        vali=$vali" $OUT/$dir/$TASKBOWTIE/"
    done
    echo "<h3>Result</h3>">>$SUMMARYTMP
    python ${NGSANE_BASE}/tools/Summary.py "$vali" $ASD.bam.stats samstats >>$SUMMARYTMP

    echo "</div></div></div>" >>$SUMMARYTMP
fi

################################################################################
if [[ -n "$RUNMAPPINGBOWTIE2" ]]; then
    PIPELINE="BOWTIE v2 mapping"
    PIPELINK="bowtie"
   
    LINKS=$LINKS" $PIPELINK"
    echo "<div class='panel'><div class='headbagb'><a name='$PIPELINK'><h2 class='sub'>$PIPELINE</h2></a></div><div class='wrapper'><div class='results'>" >>$SUMMARYTMP

    echo "QC"
    ${NGSANE_BASE}/mods/QC.sh -o -m ${NGSANE_BASE}/mods/bowtie2.sh -l $QOUT/$TASKBOWTIE2/ >>$SUMMARYTMP
    echo "gather dirs"
    vali=""
    for dir in ${DIR[@]}; do
        vali=$vali" $OUT/$dir/$TASKBOWTIE2/"
    done
    echo "<h3>Result</h3>">>$SUMMARYTMP
    python ${NGSANE_BASE}/tools/Summary.py "$vali" $ASD.bam.stats samstats >>$SUMMARYTMP

    echo "</div></div></div>" >>$SUMMARYTMP
fi

################################################################################
if [[ -n "$RUNTOPHATCUFF" || -n "$RUNTOPHATCUFF2" ]]; then
    PIPELINE="TOPHAT + Cufflinks"
    PIPELINK="tophat"
   
    LINKS=$LINKS" $PIPELINK"
    echo "<div class='panel'><div class='headbagb'><a name='$PIPELINK'><h2 class='sub'>$PIPELINE</h2></a></div><div class='wrapper'><div class='results'>" >>$SUMMARYTMP

    echo "<br>Note, the duplication rate is not calculated by tophat and hence zero." >>$SUMMARYTMP
    ${NGSANE_BASE}/mods/QC.sh -o -m ${NGSANE_BASE}/mods/tophatcuff.sh -l $QOUT/$TASKTOPHAT/ >>$SUMMARYTMP
    echo "<h3>Result</h3>">>$SUMMARYTMP
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

    echo "</div></div></div>" >>$SUMMARYTMP
fi


################################################################################
if [[ -n "$DEPTHOFCOVERAGE"  || -n "$DEPTHOFCOVERAGE2" ]]; then
    PIPELINE="COVERAGE"
    PIPELINK="coverage"
   
    LINKS=$LINKS" $PIPELINK"
     echo "<div class='panel'><div class='headbagb'><a name='$PIPELINK'><h2 class='sub'>$PIPELINE</h2></a></div><div class='wrapper'><div class='results'>" >>$SUMMARYTMP
 
    ${NGSANE_BASE}/mods/QC.sh -o -m ${NGSANE_BASE}/mods/gatkSNPs.sh -l $QOUT/$TASKVAR >> $SUMMARYTMP
    vali=""
    for dir in ${DIR[@]}; do
	   vali=$vali" $OUT/$dir/$TASKDOC/"
    done
    echo "<h3>Average coverage</h3>">>$SUMMARYTMP
    python ${NGSANE_BASE}/tools/Summary.py "$vali" $ASR".bam.doc.sample_summary" coverage >>$SUMMARYTMP
    echo "<h3>Base pair coverage over all intervals</h3>" >>$SUMMARYTMP
    python ${NGSANE_BASE}/tools/Summary.py "$vali" $ASR".bam.doc.sample_cumulative_coverage_counts" coverage --p >>$SUMMARYTMP
    echo "<h3>Intervals covered</h3>" >>$SUMMARYTMP
    python ${NGSANE_BASE}/tools/Summary.py "$vali" $ASR".bam.doc.sample_interval_statistics" coverage --p >>$SUMMARYTMP
    echo "<h3>On Target</h3>" >>$SUMMARYTMP
    python ${NGSANE_BASE}/tools/Summary.py "$vali" $ASR".bam.stats" target >>$SUMMARYTMP

    echo "</div></div></div>" >>$SUMMARYTMP
    
fi

################################################################################
if [ -n "$RUNVARCALLS" ]; then 
    PIPELINE="Variant calling"
    PIPELINK="varcall"
   
    LINKS=$LINKS" $PIPELINK"
    echo "<div class='panel'><div class='headbagb'><a name='$PIPELINK'><h2 class='sub'>$PIPELINE</h2></a></div><div class='wrapper'><div class='results'>" >>$SUMMARYTMP

    ${NGSANE_BASE}/mods/QC.sh -o -l ${NGSANE_BASE}/mods/gatkSNPs.sh -m $QOUT/$TASKVAR >> $SUMMARYTMP

    vali=""
    for dir in ${DIR[@]}; do
	   vali=$vali" $OUT/$TASKVAR/$dir/"
    done
    echo "<h3>SNPs</h3>">>$SUMMARYTMP
    python ${NGSANE_BASE}/tools/Summary.py "$vali" "filter.snps.eval.txt" variant --n --l>>$SUMMARYTMP
    python ${NGSANE_BASE}/tools/Summary.py "$vali" "recalfilt.snps.eval.txt" variant --n --l>>$SUMMARYTMP
    echo "<h3>INDELs</h3>" >>$SUMMARYTMP
    python ${NGSANE_BASE}/tools/Summary.py "$vali" "filter.indel.eval.txt" variant --n --l>>$SUMMARYTMP

    echo "</div></div></div>" >>$SUMMARYTMP
fi

################################################################################
if [ -n "$RUNANNOTATION" ]; then
    PIPELINE="Variant annotation"
    PIPELINK="varanno"
   
    LINKS=$LINKS" $PIPELINK"
    echo "<div class='panel'><div class='headbagb'><a name='$PIPELINK'><h2 class='sub'>$PIPELINE</h2></a></div><div class='wrapper'><div class='results'>" >>$SUMMARYTMP

    ${NGSANE_BASE}/mods/QC.sh -o -m ${NGSANE_BASE}/mods/annovar.sh -l $QOUT/$TASKANNOVAR >> $SUMMARYTMP
    
    vali=""
    for dir in ${DIR[@]}; do
        vali=$vali" "$( ls $OUT/$TASKANNOVAR/$dir/*.csv)
    done
    echo "<h3>Annotation Files</h3>">>$SUMMARYTMP
    echo "Please right click the link and choose \"Save as...\" to download.<br><br>">> $SUMMARYTMP
    for v in $vali; do
        name=`basename $v`
        address=${v/\/illumina/http:\/\/cluster-vm.qbi.uq.edu.au}
        echo "<a href=\"$address\">$name</a><br>" >> $SUMMARYTMP
    done
    echo "<br>More information about the columns can be found on the <a target=new href=\"http://redmine.qbi.uq.edu.au/knowledgebase/articles/12\">Project Server</a> (uqlogin). and the description of the <a href=\"http://www.broadinstitute.org/gsa/wiki/index.php/Understanding_the_Unified_Genotyper%27s_VCF_files\">vcf file</a>">> $SUMMARYTMP
    
    echo "</div></div></div>" >>$SUMMARYTMP
fi

################################################################################
if [ -n "$RUNTRIMGALORE" ];then
    PIPELINE="Trimgalore trimming"
    PIPELINK="trimgalore"
   
    LINKS=$LINKS" $PIPELINK"
    echo "<div class='panel'><div class='headbagb'><a name='$PIPELINK'><h2 class='sub'>$PIPELINE</h2></a></div><div class='wrapper'><div class='results'>" >>$SUMMARYTMP

    ${NGSANE_BASE}/mods/QC.sh -o -m ${NGSANE_BASE}/mods/trimgalore.sh -l $QOUT/$TASKTRIMGALORE >> $SUMMARYTMP

    echo "<h3>trimgalore</h3>">>$SUMMARYTMP
    vali=""
    for dir in ${DIR[@]}; do
        vali=$vali" $OUT/fastq/${dir/_$TASKTRIMGALORE/}_$TASKTRIMGALORE/"
    done
    python ${NGSANE_BASE}/tools/Summary.py "$vali" "_trimming_report.txt" trimgalore --noSummary >> $SUMMARYTMP

    echo "</div></div></div>" >>$SUMMARYTMP
fi

################################################################################
if [ -n "$RUNTRIMMOMATIC" ];then
    PIPELINE="Trimmomatic trimming"
    PIPELINK="trimmomatic"
   
    LINKS=$LINKS" $PIPELINK"
    echo "<div class='panel'><div class='headbagb'><a name='$PIPELINK'><h2 class='sub'>$PIPELINE</h2></a></div><div class='wrapper'><div class='results'>" >>$SUMMARYTMP

    ${NGSANE_BASE}/mods/QC.sh -o -m ${NGSANE_BASE}/mods/trimmomatic.sh -l $QOUT/$TASKTRIMMOMATIC >> $SUMMARYTMP

    echo "<h3>trimmomatic</h3>">>$SUMMARYTMP
    vali=""
    for dir in ${DIR[@]}; do
        vali=$vali" $OUT/fastq/${dir/_$TASKTRIMMOMATIC/}_$TASKTRIMMOMATIC/"
    done
    python ${NGSANE_BASE}/tools/Summary.py "$vali" ".log" trimmomatic --noSummary >> $SUMMARYTMP

    echo "</div></div></div>" >>$SUMMARYTMP
fi

################################################################################
if [ -n "$RUNCUTADAPT" ];then
    PIPELINE="Cutadapt trimming"
    PIPELINK="cutadapt"
   
    LINKS=$LINKS" $PIPELINK"
    echo "<div class='panel'><div class='headbagb'><a name='$PIPELINK'><h2 class='sub'>$PIPELINE</h2></a></div><div class='wrapper'><div class='results'>" >>$SUMMARYTMP

    ${NGSANE_BASE}/mods/QC.sh -o -m ${NGSANE_BASE}/mods/cutadapt.sh -l $QOUT/$TASKCUTADAPT >> $SUMMARYTMP

    echo "<h3>cutadapt</h3>">>$SUMMARYTMP
    vali=""
    for dir in ${DIR[@]}; do
        vali=$vali" $OUT/fastq/${dir/_$TASKCUTADAPT/}_$TASKCUTADAPT/"
    done
    python ${NGSANE_BASE}/tools/Summary.py "$vali" ".stats" cutadapt --noSummary >> $SUMMARYTMP

    echo "</div></div></div>" >>$SUMMARYTMP
fi

################################################################################
if [ -n "$RUNHICLIB" ];then
    PIPELINE="HiClib"
    PIPELINK="hiclib"
   
    LINKS=$LINKS" $PIPELINK"
    echo "<div class='panel'><div class='headbagb'><a name='$PIPELINK'><h2 class='sub'>$PIPELINE</h2></a></div><div class='wrapper'><div class='results'>" >>$SUMMARYTMP

    ${NGSANE_BASE}/mods/QC.sh -o -l ${NGSANE_BASE}/mods/hiclibMapping.sh -m $QOUT/$TASKHICLIB >> $SUMMARYTMP

    echo "<h3>hiclib</h3>">>$SUMMARYTMP
    vali=""
    for dir in ${DIR[@]}; do
        vali=$vali" $OUT/$dir/$TASKHICLIB/"
    done
    python ${NGSANE_BASE}/tools/Summary.py "$vali" ".log" hiclibMapping >> $SUMMARYTMP
    for dir in $vali; do
        for pdf in $(ls -f $dir/*.pdf 2>/dev/null ); do
            echo "<a href='$pdf'>${pdf##*/}</a> " >> $SUMMARYTMP
        done
    done

    echo "</div></div></div>" >>$SUMMARYTMP
fi

################################################################################
if [ -n "$RUNHICUP" ];then
    PIPELINE="HiCUP + fit-hi-C"
    PIPELINK="hicup"
   
    LINKS=$LINKS" $PIPELINK"
    echo "<div class='panel'><div class='headbagb'><a name='$PIPELINK'><h2 class='sub'>$PIPELINE</h2></a></div><div class='wrapper'><div class='results'>" >>$SUMMARYTMP

    ${NGSANE_BASE}/mods/QC.sh -o -l ${NGSANE_BASE}/mods/hicup.sh -m $QOUT/$TASKHICUP >> $SUMMARYTMP

    vali=""
    for dir in ${DIR[@]}; do
	   vali=$vali" $OUT/$dir/$TASKHICUP/"
    done
    
    echo "<h3>hicup</h3>">>$SUMMARYTMP
    echo "<h4>truncater</h4>">>$SUMMARYTMP
    python ${NGSANE_BASE}/tools/Summary.py "$vali" "hicup_truncater_summary.txt" hicup --noSummary --noOverallSummary >> $SUMMARYTMP
    echo "<h4>mapper</h4>">>$SUMMARYTMP
    python ${NGSANE_BASE}/tools/Summary.py "$vali" "hicup_mapper_summary.txt" hicup --noSummary --noOverallSummary >> $SUMMARYTMP
    echo "<h4>filter</h4>">>$SUMMARYTMP
    python ${NGSANE_BASE}/tools/Summary.py "$vali" "hicup_filter_summary_results.txt" hicup --noSummary --noOverallSummary >> $SUMMARYTMP
    echo "<h4>deduplicator</h4>">>$SUMMARYTMP
    python ${NGSANE_BASE}/tools/Summary.py "$vali" "hicup_deduplicater_summary_results.txt" hicup --noSummary --noOverallSummary >> $SUMMARYTMP
    
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
    
    echo "</div></div></div>" >>$SUMMARYTMP
fi

################################################################################
if [ -n "$RUNHOMERCHIPSEQ" ];then
    PIPELINE="Homer ChIP-Seq"
    PIPELINK="homerchipseq"
   
    LINKS=$LINKS" $PIPELINK"
    echo "<div class='panel'><div class='headbagb'><a name='$PIPELINK'><h2 class='sub'>$PIPELINE</h2></a></div><div class='wrapper'><div class='results'>" >>$SUMMARYTMP

    ${NGSANE_BASE}/mods/QC.sh -o -m ${NGSANE_BASE}/mods/chipseqHomer.sh -l $QOUT/$TASKHOMERCHIPSEQ >> $SUMMARYTMP

    echo "<h3>Homer ChIP-seq</h3>">>$SUMMARYTMP
    vali=""
    for dir in ${DIR[@]}; do
        vali=$vali" $OUT/$dir/$TASKHOMERCHIPSEQ/"
    done
    python ${NGSANE_BASE}/tools/Summary.py "$vali" ".summary.txt" homerchipseq >> $SUMMARYTMP

    echo "</div></div></div>" >>$SUMMARYTMP
fi

################################################################################
if [ -n "$RUNPEAKRANGER" ];then
    PIPELINE="Peakranger"
    PIPELINK="peakranger"
   
    LINKS=$LINKS" $PIPELINK"
    echo "<div class='panel'><div class='headbagb'><a name='$PIPELINK'><h2 class='sub'>$PIPELINE</h2></a></div><div class='wrapper'><div class='results'>" >>$SUMMARYTMP

    ${NGSANE_BASE}/mods/QC.sh -o -m ${NGSANE_BASE}/mods/peakranger.sh -l $QOUT/$TASKPEAKRANGER >> $SUMMARYTMP

    echo "<h3>Peakranger</h3>">>$SUMMARYTMP
    vali=""
    for dir in ${DIR[@]}; do
        vali=$vali" $OUT/$dir/$TASKPEAKRANGER/"
    done
    python ${NGSANE_BASE}/tools/Summary.py "$vali" ".summary.txt" peakranger >> $SUMMARYTMP
    echo "</div>" >>$SUMMARYTMP
fi

################################################################################
if [ -n "$RUNMACS2" ];then
    PIPELINE="MACS2"
    PIPELINK="mac2"
   
    LINKS=$LINKS" $PIPELINK"
    echo "<div class='panel'><div class='headbagb'><a name='$PIPELINK'><h2 class='sub'>$PIPELINE</h2></a></div><div class='wrapper'><div class='results'>" >>$SUMMARYTMP

    ${NGSANE_BASE}/mods/QC.sh -o -m ${NGSANE_BASE}/mods/macs2.sh -l $QOUT/$TASKMACS2 >> $SUMMARYTMP

    echo "<h3>MACS2</h3>">>$SUMMARYTMP
    vali=""
    for dir in ${DIR[@]}; do
        vali=$vali" $OUT/$dir/$TASKMACS2/"
    done
    python ${NGSANE_BASE}/tools/Summary.py "$vali" ".summary.txt" macs2 >> $SUMMARYTMP

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

    echo "</div></div></div>" >>$SUMMARYTMP
fi

################################################################################
if [ -n "$RUNMEMECHIP" ];then
    PIPELINE="MEME-chip Motif discovery"
    PIPELINK="meme"
   
    LINKS=$LINKS" $PIPELINK"
    echo "<div class='panel'><div class='headbagb'><a name='$PIPELINK'><h2 class='sub'>$PIPELINE</h2></a></div><div class='wrapper'><div class='results'>" >>$SUMMARYTMP

    ${NGSANE_BASE}/mods/QC.sh -o -m ${NGSANE_BASE}/mods/memechip.sh -l $QOUT/$TASKMEMECHIP >> $SUMMARYTMP

    echo "<h3>MEME-chip</h3>">>$SUMMARYTMP
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

    echo "</div></div></div>" >>$SUMMARYTMP
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
	echo "<h3>Annotation of mapped reads</h3>">>$SUMMARYTMP
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
echo '''
<html>
<head>
<style type="text/css">
body {
	font: 12px/19px "Lucida Sans", "Lucida Grande", "Lucida Sans Unicode",
		Verdana, sans-serif;
}

h2 {
	display: block;
	font-size: 1.17em;
	-webkit-margin-before: 1em;
	-webkit-margin-after: 1em;
	-webkit-margin-start: 0px;
	-webkit-margin-end: 0px;
	font-weight: bold;
}

.panel h2 {
	background-color: #999;
	border-radius: 10px;
	color: white;
	font-size: 100%;
	font-weight: normal;
	letter-spacing: 0.2em;
	text-transform: uppercase;
	padding: 2px 10px;
	margin: 0px;
}

.panel h3 {
	background-color: #69c;
	color: white;
	font-size: 100%;
	font-weight: normal;
	letter-spacing: 0.2em;
	text-transform: uppercase;
	padding: 2px 10px;
	margin: 0px;
}

div.panel h2,div.panel h2.sub.inactive {
	background: #888;
	background-image: -webkit-gradient(linear, left top, left bottom, from(#69c),
		to(#669));
	background-image: -webkit-linear-gradient(top, #69c, #669);
	background-image: -moz-linear-gradient(top, #69c, #669);
	background-image: -ms-linear-gradient(top, #69c, #669);
	background-image: -o-linear-gradient(top, #69c, #669);
	background-image: linear-gradient(top, #69c, #669);
}

div.panel h2.sub.inactive {
	color: #dedede;
	cursor: pointer;
	border-top-color: #ccc;
	border-left-color: #aaa;
}

div.headbagb {
	width: 100%;
	border-radius: 10px 10px 0 0;
	background: #555;
	background-image: -webkit-gradient(linear, left top, left bottom, from(#69c),
		to(#669));
	background-image: -webkit-linear-gradient(top, #69c, #669);
	background-image: -moz-linear-gradient(top, #69c, #669);
	background-image: -ms-linear-gradient(top, #69c, #669);
	background-image: -o-linear-gradient(top, #69c, #669);
	background-image: linear-gradient(top, #69c, #669);
}

div.panel {
	background: #fafaff;
	border: 1px solid #999;
	border-radius: 12px;
	padding: 2px;
	margin-bottom: 25px;
}

#controls div,#gallery div {
	margin: 2px 0;
	width: 240px;
}

#quicklinks a {
	display: inline-block;
	padding: 2px;
	margin: 2px;
	outline: 0;
	color: #333;
	transition-duration: 0.25s;
	transition-property: transform;
	transform: scale(1) rotate(0);
}

#quicklinks a:hover {
	background: #aaa;
	text-decoration: none;
	color: #000;
	border-radius: 10px;
	transform: scale(1.05) rotate(-1deg);
}

div.panel h2.sub {
	margin-left: 2px;
	margin-top: 2px;
	background: #fafaff;
	border-radius: 10px 10px 0 0;
	border: 1px solid #666;
	border-top-color: #EEE;
	border-left-color: #DDD;
	color: #444;
	border-bottom: none;
	display: inline-block;
	vertical-align: top;
	width: 400px;
}

div.wrapper {
	margin: auto; 
	width: 100%;
	padding-top: 5px;
}

div.results {
	overflow-x: scroll;
	overflow-y: hidden;
	margin-bottom: 10px;
}

p {
	color: blue;
}

table {
	border-collapse: collapse;
	border-spacing: 0;
}

table.data {
        table-layout: fixed;
	font-size: 12px;
	text-align: left;
	border-collapse: collapse;
	border: 1px solid #69c;
	margin: 20px 5px 20px 5px;
}

table,caption,tbody,tfoot,thead,tr,th,td {
	font-size: 100%;
	vertical-align: baseline;
	margin: 0;
	padding: 0;
	outline: 0;
	border: 0;
	background: transparent;
}

table.data th {
	font-weight: bold;
	font-size: 12px;
	color: #039;
	border-bottom: 1px dashed #69c;
	padding: 12px 5px;
	text-align:right;
}

table.data th div{
	width: 100px;
}

table.data td {
	color: #669;
	padding: 5px 5px;
	text-align:right;
	white-space: nowrap;
}

table.data td.left, table.data  th.left{
	text-align:left;
	width: auto;
}

table.data tfoot {
	border-top: 1px dashed #69c;
}

div.library {
	background: #fff;
	margin: 5px;
	padding: 5px;
	border-radius: 10px;
	border: 1px solid #666;
	max-width: 1200px;
}

hr {
	border: 0;
	border-top: 1px solid #ccc;
}
</style>
</head>
<body>

<div id="center">
''' > $SUMMARYFILE.tmp
echo "<div class='panel' id='quicklinks'><h2>Quicklink</h2><div>" >> $SUMMARYFILE.tmp
for i in $LINKS; do
    echo "<a href=#$i>$i</a> | ">>$SUMMARYFILE.tmp
done
echo "</div></div>" >>$SUMMARYFILE.tmp

cat $SUMMARYFILE.tmp  $SUMMARYTMP > $SUMMARYFILE

echo "</div><!-- center --></body></html>" >>$SUMMARYFILE

rm $SUMMARYTMP
rm $SUMMARYFILE.tmp

################################################################################
# convert html to pdf
#if [ "$(hash prince)" == "" ]; then
#    prince $SUMMARYFILE -o ${HTMLOUT}.pdf
#fi

