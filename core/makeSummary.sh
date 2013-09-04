#!/bin/bash

# author: Fabian Buske
# date: Aug 2013

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
for MODULE in $MODULE_SUMMARY; do module load $MODULE; done  # save way to load modules that itself load other modules
for MODULE in $MODULE_R; do module load $MODULE; done  # save way to load modules that itself load other modules

echo -e "--R           --\n "$(R --version | head -n 3)
[ -z "$(which R)" ] && echo "[ERROR] no R detected" && exit 1
echo -e "--Python      --\n "$(python --version 2>&1 | tee | head -n 1)
[ -z "$(which python)" ] && echo "[ERROR] no Python detected" && exit 1

################################################################################

SUMMARYTMP=$HTMLOUT".tmp"
SUMMARYFILE=$HTMLOUT".html"

echo "Last modified "`date` >$SUMMARYTMP

################################################################################
# define functions for generating summary scaffold
#
# summaryHeader takes 4 parameters
# $1=PIPELINE name
# $2=TASK (e.g. $TASKBWA)
# $3=pipeline mod script (e.g. bwa.sh)
# $4=output file ($SUMMARYTMP)
function summaryHeader {
    LINKS=$LINKS" $2"
    echo "<div class='panel'><div class='headbagb'><a name='$2'><h2 class='sub'>$1</h2></a></div><div class='wrapper'><div class='results'>" >> $4
    echo "QC"
    ${NGSANE_BASE}/core/QC.sh -o -m ${NGSANE_BASE}/mods/$3 -l $QOUT/$2 >> $4
} 

# summaryFooter takes 1 parameter
# $1=output file ($SUMMARYTMP)
function summaryFooter {
    echo "</div></div></div>" >> $1
}
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
    summaryHeader "FASTQ screen" "TASKFASTQSCREEN" "fastqscreen.sh" "$SUMMARYTMP"

    echo "gather dirs"
    vali=""
    for dir in ${DIR[@]}; do
        vali=$vali" $OUT/$dir/$TASKFASTQSCREEN/"
    done
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

    summaryFooter "$SUMMARYTMP"
fi


################################################################################
if [[ -n "$RUNMAPPINGBWA" || -n "$RUNMAPPINGBWA2" ]]; then
    summaryHeader "BWA mapping" "$TASKBWA" "bwa.sh" "$SUMMARYTMP"

    echo "gather dirs"
    for dir in ${DIR[@]}; do
	   vali=$vali" $OUT/$dir/$TASKBWA/"
    done
    python ${NGSANE_BASE}/tools/Summary.py "$vali" $ASD.bam.stats samstats >>$SUMMARYTMP

    if [ -n "$RUNANNOTATINGBAM" ]; then
    	echo "<h3>Annotation</h3>" >>$SUMMARYTMP
    	python ${NGSANE_BASE}/tools/Summary.py "$vali" merg.anno.stats annostats >>$SUMMARYTMP
    	ROUTH=runStats/$(echo ${DIR[@]}|sed 's/ /_/g')
    	if [ ! -e $ROUTH ]; then mkdir $ROUTH; fi
	   python ${NGSANE_BASE}/tools/makeBamHistogram.py "$vali" $ROUTH >>$SUMMARYTMP
    fi
    
    summaryFooter "$SUMMARYTMP"
fi


################################################################################
if [[ -n "$RUNREALRECAL" || -n "$RUNREALRECAL2" || -n "$RUNREALRECAL3" ]]; then 
    summaryHeader "RECAL mapping" "$TASKRCA" "reCalAln.sh" "$SUMMARYTMP"

    vali=""
    for dir in ${DIR[@]}; do
	   vali=$vali" $OUT/$dir/$TASKRCA/"
    done

    python ${NGSANE_BASE}/tools/Summary.py "$vali" $ASR".bam.stats" samstatsrecal >>$SUMMARYTMP

    summaryFooter "$SUMMARYTMP"
fi

################################################################################
if [[ -n "$RUNMAPPINGBOWTIE" ]]; then
    summaryHeader "BOWTIE v1 mapping" "$TASKBOWTIE" "bowtie.sh" "$SUMMARYTMP"

    vali=""
    for dir in ${DIR[@]}; do
        vali=$vali" $OUT/$dir/$TASKBOWTIE/"
    done

    python ${NGSANE_BASE}/tools/Summary.py "$vali" $ASD.bam.stats samstats >>$SUMMARYTMP

    summaryFooter "$SUMMARYTMP"
fi

################################################################################
if [[ -n "$RUNMAPPINGBOWTIE2" ]]; then
    summaryHeader "BOWTIE v2 mapping" "$TASKBOWTIE2" "bowtie2.sh" "$SUMMARYTMP"

    echo "gather dirs"
    vali=""
    for dir in ${DIR[@]}; do
        vali=$vali" $OUT/$dir/$TASKBOWTIE2/"
    done

    python ${NGSANE_BASE}/tools/Summary.py "$vali" $ASD.bam.stats samstats >>$SUMMARYTMP

    summaryFooter "$SUMMARYTMP"
fi

################################################################################
if [[ -n "$RUNTOPHATCUFF" || -n "$RUNTOPHATCUFF2" ]]; then
    summaryHeader "TOPHAT + Cufflinks" "$TASKTOPHAT" "tophatcuff.sh" "$SUMMARYTMP"


    echo "<br>Note, the duplication rate is not calculated by tophat and hence zero." >>$SUMMARYTMP
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

    summaryFooter "$SUMMARYTMP"
fi


################################################################################
if [[ -n "$DEPTHOFCOVERAGE"  || -n "$DEPTHOFCOVERAGE2" ]]; then
    summaryHeader "COVERAGE" "$TASKVAR" "gatkSNPs.sh" "$SUMMARYTMP"

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

    summaryFooter "$SUMMARYTMP"
    
fi

################################################################################
if [ -n "$RUNVARCALLS" ]; then 
    summaryHeader "Variant calling" "$TASKVAR" "gatkSNPs.sh" "$SUMMARYTMP"

    vali=""
    for dir in ${DIR[@]}; do
	   vali=$vali" $OUT/$TASKVAR/$dir/"
    done
    echo "<h3>SNPs</h3>">>$SUMMARYTMP
    python ${NGSANE_BASE}/tools/Summary.py "$vali" "filter.snps.eval.txt" variant --n --l>>$SUMMARYTMP
    python ${NGSANE_BASE}/tools/Summary.py "$vali" "recalfilt.snps.eval.txt" variant --n --l>>$SUMMARYTMP
    echo "<h3>INDELs</h3>" >>$SUMMARYTMP
    python ${NGSANE_BASE}/tools/Summary.py "$vali" "filter.indel.eval.txt" variant --n --l>>$SUMMARYTMP

    summaryFooter "$SUMMARYTMP"
fi

################################################################################
if [ -n "$RUNANNOTATION" ]; then
    summaryHeader "Variant annotation" "$TASKANNOVAR" "annovar.sh" "$SUMMARYTMP"

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
    
    summaryFooter "$SUMMARYTMP"
fi

################################################################################
if [ -n "$RUNTRIMGALORE" ];then
    summaryHeader "Trimgalore trimming" "$TASKTRIMGALORE" "trimgalore.sh" "$SUMMARYTMP"


    vali=""
    for dir in ${DIR[@]}; do
        vali=$vali" $OUT/fastq/${dir/_$TASKTRIMGALORE/}_$TASKTRIMGALORE/"
    done
    python ${NGSANE_BASE}/tools/Summary.py "$vali" "_trimming_report.txt" trimgalore --noSummary >> $SUMMARYTMP

    summaryFooter "$SUMMARYTMP"
fi

################################################################################
if [ -n "$RUNTRIMMOMATIC" ];then
    summaryHeader "Trimmomatic trimming" "$TASKTRIMMOMATIC" "trimmomatic.sh" "$SUMMARYTMP"

    echo "<h3>trimmomatic</h3>">>$SUMMARYTMP
    vali=""
    for dir in ${DIR[@]}; do
        vali=$vali" $OUT/fastq/${dir/_$TASKTRIMMOMATIC/}_$TASKTRIMMOMATIC/"
    done
    python ${NGSANE_BASE}/tools/Summary.py "$vali" ".log" trimmomatic --noSummary >> $SUMMARYTMP

    summaryFooter "$SUMMARYTMP"
fi

################################################################################
if [ -n "$RUNCUTADAPT" ];then
    summaryHeader "Cutadapt trimming" "$TASKCUTADAPT" "cutadapt.sh" "$SUMMARYTMP"

    echo "<h3>cutadapt</h3>">>$SUMMARYTMP
    vali=""
    for dir in ${DIR[@]}; do
        vali=$vali" $OUT/fastq/${dir/_$TASKCUTADAPT/}_$TASKCUTADAPT/"
    done
    python ${NGSANE_BASE}/tools/Summary.py "$vali" ".stats" cutadapt --noSummary >> $SUMMARYTMP

    summaryFooter "$SUMMARYTMP"
fi

################################################################################
if [ -n "$RUNHICLIB" ];then
    summaryHeader "HiClib" "$TASKHICLIB" "hiclibMapping.sh" "$SUMMARYTMP"

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

    summaryFooter "$SUMMARYTMP"
fi

################################################################################
if [ -n "$RUNHICUP" ];then
    summaryHeader "HiCUP + fit-hi-C" "$TASKHICUP" "hicup.sh" "$SUMMARYTMP"

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
    
    summaryFooter "$SUMMARYTMP"
fi

################################################################################
if [ -n "$RUNHOMERCHIPSEQ" ];then
    summaryHeader "Homer ChIP-Seq" "$TASKHOMERCHIPSEQ" "chipseqHomer.sh" "$SUMMARYTMP"

    echo "<h3>Homer ChIP-seq</h3>">>$SUMMARYTMP
    vali=""
    for dir in ${DIR[@]}; do
        vali=$vali" $OUT/$dir/$TASKHOMERCHIPSEQ/"
    done
    python ${NGSANE_BASE}/tools/Summary.py "$vali" ".summary.txt" homerchipseq >> $SUMMARYTMP

    summaryFooter "$SUMMARYTMP"
fi

################################################################################
if [ -n "$RUNPEAKRANGER" ];then
    summaryHeader "Peakranger" "$TASKPEAKRANGER" "peakranger.sh" "$SUMMARYTMP"

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
    summaryHeader "MACS2" "$TASKMACS2" "macs2.sh" "$SUMMARYTMP"

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

    summaryFooter "$SUMMARYTMP"
fi

################################################################################
if [ -n "$RUNMEMECHIP" ];then
    summaryHeader "MEME-chip Motif discovery" "$TASKMEMECHIP" "memechip.sh" "$SUMMARYTMP"

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

    summaryFooter "$SUMMARYTMP"
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
echo '''<html><head>''' > $SUMMARYFILE.tmp
cat ${NGSANE_BASE}/core/Summary.css >> $SUMMARYFILE.tmp

echo '<script type="javascript">' >> $SUMMARYFILE.tmp
cat ${NGSANE_BASE}/core/jquery-1.9.1.min.js >> $SUMMARYFILE.tmp
cat ${NGSANE_BASE}/core/Summary.js >> $SUMMARYFILE.tmp
echo '''</script></head><body>

<div id="center">
''' >> $SUMMARYFILE.tmp
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

################################################################################
echo ">>>>> Generate HTML report - FINISHED"
echo ">>>>> enddate "`date`
