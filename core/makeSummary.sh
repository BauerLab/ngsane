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

mkdir -p $(dirname $SUMMARYTMP) && cat /dev/null > $SUMMARYTMP # clean temporary content

PROJECT_RELPATH=$(python -c "import os.path; print os.path.relpath('$(pwd)',os.path.abspath('$(dirname $SUMMARYTMP)'))")

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
    echo "<div class='panel' id='$2_panel'>
        <div class='headbagb' id='$2_panelback'><a name='$2'></a>
            <h2 id='$2_h_results' class='sub'>$1</h2>
            <h2 id='$2_h_checklist' class='sub inactive' rel='checklist'>Checklist<span class='counter'><span class='passed' id='$2_counter_checkpoints_passed'></span><span class='failed' id='$2_counter_checkpoints_failed'></span></span></h2>
            <h2 id='$2_h_notes' class='sub inactive' rel='notes'>Notes<span class='counter'><span class='neutral' id='$2_counter_notes'></span></span></span></h2>
            <h2 id='$2_h_errors' class='sub inactive' rel='notes'>Errors<span class='counter'><span class='errors' id='$2_counter_errors'></span></span></h2>
            <h2 id='$2_h_logfiles' class='sub inactive' rel='errors'>Logfiles</h2>
        </div>
        <div class='wrapper'><div class='hidden'>" >> $4
    echo "QC - $2"
    ${NGSANE_BASE}/core/QC.sh -o $4 -m ${NGSANE_BASE}/mods/$3 -l $QOUT -t $2 >> $4
    echo "<div id='$2_results'>" >> $4
} 

# summaryFooter takes 2 parameters
# $1=TASK (e.g. $TASKBWA)
# $2=output file ($SUMMARYTMP)
function summaryFooter {
    echo "</div></div><div class='display'></div></div></div>" >> $2
    echo "<script type='text/javascript'> 
        if (typeof jQuery === 'undefined') {
            console.log('jquery not loaded');
        } else {
            \$('#$1_panel div.wrapper div.display').html(\$('#$1_results').html());
        }
        </script>" >> $2
    
}

# gatherDirs takes 1 parameter
# $1=TASK (e.g. $TASKBWA)
function gatherDirs {
    vali=""
    for dir in ${DIR[@]}; do
        [ -d $OUT/$dir/$1/ ] && vali=$vali" $OUT/$dir/$1/"
    done
	echo $vali

}

# gatherDirsAggregate takes 1 parameter
# $1=TASK (e.g. $TASKVAR)
function gatherDirsAggregate {
    vali=""
	for dir in ${DIR[@]}; do
	   [ -d $OUT/$1/$dir/ ] && vali=$vali" $OUT/$1/$dir/"
    done
	echo $vali

}


################################################################################
if [ -n "$RUNFASTQC" ]; then
    PIPELINE="FASTQC"
    PIPELINK="fastqc"
    
    LINKS=$LINKS" $PIPELINK"
    echo "<div class='panel'><div class='headbagb'><a name='$PIPELINK'><h2 class='sub'>$PIPELINE</h2></a></div>" >>$SUMMARYTMP

    echo "<table class='data'>" >>$SUMMARYTMP
    echo "<thead><tr><th class='left'>Libary</th><th><div style='width:100px'>Chart</div></th><th><div style='width:140px'>Encoding</div></th><th><div style='width:120px'>Library size</div></th><th><div style='width:50px'>Read</div></th><th><div style='width:80px'>Read length</div></th><th><div style='width:50px'>%GC</div></th><th><div style='width:120px'>Read Qualities</th></tr><thead><tbody>" >>$SUMMARYTMP

    if [[ -e runStats/ && -e runStats/$TASKFASTQC/ ]]; then
        for f in $( ls runStats/$TASKFASTQC/*.zip ); do
            # get basename of f
            n=${f##*/}
            n=${n/"_fastqc.zip"/}
            ICO="<img height='15px' class='noborder' src='$PROJECT_RELPATH/runStats/$TASKFASTQC/"$n"_fastqc/Icons/"
            P=$(grep "PASS" -c runStats/$TASKFASTQC/$n"_fastqc/summary.txt")
            W=$(grep "WARN" -c runStats/$TASKFASTQC/$n"_fastqc/summary.txt")
            F=$(grep "FAIL" -c runStats/$TASKFASTQC/$n"_fastqc/summary.txt")
            CHART=$ICO"tick.png' title='$P'\>$P"
            if [ "$W" -ne "0" ]; then CHART=$CHART""$ICO"warning.png'\>"$W; fi
            if [ "$F" -ne "0" ]; then CHART=$CHART""$ICO"error.png'\>"$F; fi
            ENCODING=$(grep "Encoding" runStats/$TASKFASTQC/$n"_fastqc/fastqc_data.txt" | head -n 1 | cut -f 2)
            LIBRARYSIZE=$(grep "Total Sequences" runStats/$TASKFASTQC/$n"_fastqc/fastqc_data.txt" | head -n 1 | cut -f 2)
            READLENGTH=$(grep "Sequence length" runStats/$TASKFASTQC/$n"_fastqc/fastqc_data.txt" | head -n 1 | cut -f 2)
            GCCONTENT=$(grep "\%GC" runStats/$TASKFASTQC/$n"_fastqc/fastqc_data.txt" | head -n 1 | cut -f 2)
            if [[ "$f" == *$READTWO* ]] && [ "$f" != "${f/$READTWO/$READONE}" ]; then
                READ=2
            else
                READ=1
            fi
            echo "<tr style='vertical-align: middle;'><td class='left'><a href='$PROJECT_RELPATH/runStats/$TASKFASTQC/"$n"_fastqc/fastqc_report.html'>$n.fastq</a></td><td>$CHART</td><td>$ENCODING</td><td>$LIBRARYSIZE</td><td>$READ</td><td>$READLENGTH</td><td>$GCCONTENT</td><td>" >>$SUMMARYTMP

            if [[ "$f" == *$READONE* ]]; then
                echo "<a href='$PROJECT_RELPATH/runStats/$TASKFASTQC/${n}_fastqc/fastqc_report.html'><img src='$PROJECT_RELPATH/runStats/$TASKFASTQC/${n}_fastqc/Images/per_base_quality.png' height=75 alt='Quality scores for all first reads'/></a>" >>$SUMMARYTMP
                
                if [ -e ${f/$READONE/$READTWO} ] && [ "$f" != "${f/$READONE/$READTWO}" ]; then
                     echo "<a href='$PROJECT_RELPATH/runStats/$TASKFASTQC/${n/$READONE/$READTWO}_fastqc/fastqc_report.html'><img src='$PROJECT_RELPATH/runStats/$TASKFASTQC/${n/$READONE/$READTWO}_fastqc/Images/per_base_quality.png' height=75 alt='Quality scores for all second reads'/></a>" >>$SUMMARYTMP
                fi
            
            elif [[ "$f" == *$READTWO* ]] && [ "$f" != "${f/$READTWO/$READONE}" ]; then
                echo "<a href='$PROJECT_RELPATH/runStats/$TASKFASTQC/${n/$READTWO/$READONE}_fastqc/fastqc_report.html'><img src='$PROJECT_RELPATH/runStats/$TASKFASTQC/${n/$READTWO/$READONE}_fastqc/Images/per_base_quality.png' height=75 alt='Quality scores for all first reads'/></a>" >>$SUMMARYTMP
                echo "<a href='$PROJECT_RELPATH/runStats/$TASKFASTQC/${n}_fastqc/fastqc_report.html'><img src='$PROJECT_RELPATH/runStats/$TASKFASTQC/${n}_fastqc/Images/per_base_quality.png' height=75 alt='Quality scores for all second reads'/></a>" >>$SUMMARYTMP
			else
				echo "[ERROR] no fastq files $f"
            fi
            echo "</td></tr>" >>$SUMMARYTMP
        done
    fi
    echo "</tbody></table>">>$SUMMARYTMP

    echo "</div></div>">>$SUMMARYTMP
fi

################################################################################
if [[ -n "$RUNFASTQSCREEN" ]]; then
    summaryHeader "FASTQ screen" "$TASKFASTQSCREEN" "fastqscreen.sh" "$SUMMARYTMP"

    vali=$(gatherDirs $TASKFASTQSCREEN)
    python ${NGSANE_BASE}/core/Summary.py "$vali" _screen.txt fastqscreen --noSummary --noOverallSummary >>$SUMMARYTMP

    row0=""
    row1=""
    for dir in $vali; do
        for f in $(ls $dir/*_screen.png); do
            n=${f##*/}
            n=${n/"_screen.png"/}

            row0+="<td>$n</td>"
            row1+="<td><a href=\"$PROJECT_RELPATH/$dir/"$n"_screen.png\"><img src=\"$PROJECT_RELPATH/$dir/"$n"_screen.png\" width=\"300px\"/></a></td>"
        done
    done
    echo "<table><tr>$row0</tr><tr>$row1</tr></table>" >> $SUMMARYTMP

    summaryFooter "$TASKFASTQSCREEN" "$SUMMARYTMP"
fi


################################################################################
if [[ -n "$RUNMAPPINGBWA" || -n "$RUNMAPPINGBWA2" ]]; then
    summaryHeader "BWA mapping" "$TASKBWA" "bwa.sh" "$SUMMARYTMP"

    vali=$(gatherDirs $TASKBWA)
    python ${NGSANE_BASE}/core/Summary.py "$vali" $ASD.bam.stats samstats >>$SUMMARYTMP

    if [ -n "$RUNANNOTATINGBAM" ]; then
    	echo "<h3>Annotation</h3>" >>$SUMMARYTMP
    	python ${NGSANE_BASE}/core/Summary.py "$vali" merg.anno.stats annostats >>$SUMMARYTMP
    	ROUTH=runStats/$(echo ${DIR[@]}|sed 's/ /_/g')
    	if [ ! -e $ROUTH ]; then mkdir $ROUTH; fi
	   python ${NGSANE_BASE}/tools/makeBamHistogram.py "$vali" $ROUTH >>$SUMMARYTMP
    fi
    
    summaryFooter "$TASKBWA" "$SUMMARYTMP"
fi


################################################################################
if [[ -n "$RUNREALRECAL" || -n "$RUNREALRECAL2" || -n "$RUNREALRECAL3" ]]; then 
    summaryHeader "RECAL mapping" "$TASKRCA" "reCalAln.sh" "$SUMMARYTMP"

    python ${NGSANE_BASE}/core/Summary.py "$(gatherDirs $TASKRCA)" $ASR".bam.stats" samstatsrecal >>$SUMMARYTMP

    summaryFooter "$TASKRCA" "$SUMMARYTMP"
fi

################################################################################
if [[ -n "$RUNMAPPINGBOWTIE" ]]; then
    summaryHeader "BOWTIE v1 mapping" "$TASKBOWTIE" "bowtie.sh" "$SUMMARYTMP"

    python ${NGSANE_BASE}/core/Summary.py "$(gatherDirs $TASKBOWTIE)" $ASD.bam.stats samstats >>$SUMMARYTMP

    summaryFooter "$TASKBOWTIE" "$SUMMARYTMP"
fi

################################################################################
if [[ -n "$RUNMAPPINGBOWTIE2" ]]; then
    summaryHeader "BOWTIE v2 mapping" "$TASKBOWTIE2" "bowtie2.sh" "$SUMMARYTMP"

    python ${NGSANE_BASE}/core/Summary.py "$(gatherDirs $TASKBOWTIE2)" $ASD.bam.stats samstats >>$SUMMARYTMP

    summaryFooter "$TASKBOWTIE2" "$SUMMARYTMP"
fi

################################################################################
if [[ -n "$RUNTOPHATCUFF" || -n "$RUNTOPHATCUFF2" ]]; then
    summaryHeader "TOPHAT + Cufflinks" "$TASKTOPHAT" "tophatcuff.sh" "$SUMMARYTMP"

	vali=""
    echo "<br>Note, the duplication rate is not calculated by tophat and hence zero.<br>" >>$SUMMARYTMP
    CURDIR=$(pwd)
    for dir in ${DIR[@]}; do
    	vali=$vali" $OUT/$dir/$TASKTOPHAT/"
    	cd $OUT/$dir/$TASKTOPHAT
    	for d in $(find . -maxdepth 1 -mindepth 1 -type d -exec basename '{}' \; | grep "RNASeQC"); do
            echo "<a href=\"$PROJECT_RELPATH/$dir/$TASKTOPHAT/$d/index.html\">RNAseq-QC for $dir/$d</a><br/>" >> $CURDIR/$SUMMARYTMP
		done
    done
    cd $CURDIR
    python ${NGSANE_BASE}/core/Summary.py "$vali" bam.stats tophat >>$SUMMARYTMP

    summaryFooter "$TASKTOPHAT" "$SUMMARYTMP"
fi


################################################################################
if [[ -n "$DEPTHOFCOVERAGE"  || -n "$DEPTHOFCOVERAGE2" ]]; then
    summaryHeader "COVERAGE" "$TASKVAR" "gatkSNPs.sh" "$SUMMARYTMP"

    vali=$(gatherDirs $TASKDOC)
    echo "<h3>Average coverage</h3>">>$SUMMARYTMP
    python ${NGSANE_BASE}/core/Summary.py "$vali" $ASR".bam.doc.sample_summary" coverage >>$SUMMARYTMP
    echo "<h3>Base pair coverage over all intervals</h3>" >>$SUMMARYTMP
    python ${NGSANE_BASE}/core/Summary.py "$vali" $ASR".bam.doc.sample_cumulative_coverage_counts" coverage --p >>$SUMMARYTMP
    echo "<h3>Intervals covered</h3>" >>$SUMMARYTMP
    python ${NGSANE_BASE}/core/Summary.py "$vali" $ASR".bam.doc.sample_interval_statistics" coverage --p >>$SUMMARYTMP
    echo "<h3>On Target</h3>" >>$SUMMARYTMP
    python ${NGSANE_BASE}/core/Summary.py "$vali" $ASR".bam.stats" target >>$SUMMARYTMP

    summaryFooter "$TASKVAR" "$SUMMARYTMP" 
fi

################################################################################
if [ -n "$RUNVARCALLS" ]; then 
    summaryHeader "Variant calling" "$TASKVAR" "gatkSNPs.sh" "$SUMMARYTMP"

	vali=$(gatherDirs $TASKVAR)
    echo "<h3>SNPs</h3>">>$SUMMARYTMP
    python ${NGSANE_BASE}/core/Summary.py "$vali" "filter.snps.eval.txt" variant --n --l>>$SUMMARYTMP
    python ${NGSANE_BASE}/core/Summary.py "$vali" "recalfilt.snps.eval.txt" variant --n --l>>$SUMMARYTMP
    echo "<h3>INDELs</h3>" >>$SUMMARYTMP
    python ${NGSANE_BASE}/core/Summary.py "$vali" "filter.indel.eval.txt" variant --n --l>>$SUMMARYTMP

    summaryFooter "$TASKVAR" "$SUMMARYTMP"
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
        address=${v/\/illumina/http:\/\/hpsc.csiro.au}
        echo "<a href=\"$address\">$name</a><br>" >> $SUMMARYTMP
    done
    echo "<br>More information about the columns can be found on the <a target=new href=\"http://redmine.qbi.uq.edu.au/knowledgebase/articles/12\">Project Server</a> (uqlogin). and the description of the <a href=\"http://www.broadinstitute.org/gsa/wiki/index.php/Understanding_the_Unified_Genotyper%27s_VCF_files\">vcf file</a>">> $SUMMARYTMP
    
    summaryFooter "$TASKANNOVAR" "$SUMMARYTMP"
fi

################################################################################
if [ -n "$RUNTRIMGALORE" ];then
    summaryHeader "Trimgalore trimming" "$TASKTRIMGALORE" "trimgalore.sh" "$SUMMARYTMP"

    vali=""
    for dir in ${DIR[@]}; do
        vali=$vali" $OUT/fastq/${dir/_$TASKTRIMGALORE/}_$TASKTRIMGALORE/"
    done
    python ${NGSANE_BASE}/core/Summary.py "$vali" "_trimming_report.txt" trimgalore --noSummary >> $SUMMARYTMP

    summaryFooter "$TASKTRIMGALORE" "$SUMMARYTMP"
fi

################################################################################
if [ -n "$RUNTRIMMOMATIC" ];then
    summaryHeader "Trimmomatic trimming" "$TASKTRIMMOMATIC" "trimmomatic.sh" "$SUMMARYTMP"

    vali=""
    for dir in ${DIR[@]}; do
        vali=$vali" $OUT/fastq/${dir/_$TASKTRIMMOMATIC/}_$TASKTRIMMOMATIC/"
    done
    python ${NGSANE_BASE}/core/Summary.py "$vali" ".log" trimmomatic --noSummary >> $SUMMARYTMP

    summaryFooter "$TASKTRIMMOMATIC" "$SUMMARYTMP"
fi

################################################################################
if [ -n "$RUNCUTADAPT" ];then
    summaryHeader "Cutadapt trimming" "$TASKCUTADAPT" "cutadapt.sh" "$SUMMARYTMP"

    vali=""
    for dir in ${DIR[@]}; do
        vali=$vali" $OUT/fastq/${dir/_$TASKCUTADAPT/}_$TASKCUTADAPT/"
    done
    python ${NGSANE_BASE}/core/Summary.py "$vali" ".stats" cutadapt --noSummary >> $SUMMARYTMP

    summaryFooter "$TASKCUTADAPT" "$SUMMARYTMP"
fi

################################################################################
if [ -n "$RUNHICLIB" ];then
    summaryHeader "HiClib" "$TASKHICLIB" "hiclibMapping.sh" "$SUMMARYTMP"

    vali=$(gatherDirs $TASKHICLIB)
    python ${NGSANE_BASE}/core/Summary.py "$vali" ".log" hiclibMapping >> $SUMMARYTMP
    for dir in $vali; do
        for pdf in $(ls -f $dir/*.pdf 2>/dev/null ); do
            echo "<a href='$pdf'>${pdf##*/}</a> " >> $SUMMARYTMP
        done
    done

    summaryFooter "$TASKHICLIB" "$SUMMARYTMP"
fi

################################################################################
if [ -n "$RUNHICUP" ];then
    summaryHeader "HiCUP + fit-hi-C" "$TASKHICUP" "hicup.sh" "$SUMMARYTMP"

    vali=$(gatherDirs $TASKHICUP)
    echo "<h4>truncater</h4>">>$SUMMARYTMP
    python ${NGSANE_BASE}/core/Summary.py "$vali" "hicup_truncater_summary.txt" hicup --noSummary --noOverallSummary >> $SUMMARYTMP
    echo "<h4>mapper</h4>">>$SUMMARYTMP
    python ${NGSANE_BASE}/core/Summary.py "$vali" "hicup_mapper_summary.txt" hicup --noSummary --noOverallSummary >> $SUMMARYTMP
    echo "<h4>filter</h4>">>$SUMMARYTMP
    python ${NGSANE_BASE}/core/Summary.py "$vali" "hicup_filter_summary_results.txt" hicup --noSummary --noOverallSummary >> $SUMMARYTMP
    echo "<h4>deduplicator</h4>">>$SUMMARYTMP
    python ${NGSANE_BASE}/core/Summary.py "$vali" "hicup_deduplicater_summary_results.txt" hicup --noSummary --noOverallSummary >> $SUMMARYTMP
    
    row0=""
    row1=""
    row2=""
    for f in $(ls runStats/$TASKHICUP/*_ditag_classification.png); do
    	n=${f##*/}
	    n=${n/"_ditag_classification.png"/}
	    
	    row0+="<td>$n</td>"
	    row1+="<td><a href=\"$PROJECT_RELPATH/runStats/$TASKHICUP/"$n"_ditag_classification.png\"><img src=\"$PROJECT_RELPATH/runStats/$TASKHICUP/"$n"_ditag_classification.png\" width=\"200px\"/></a></td>"
   	    row2+="<td><a href=\"$PROJECT_RELPATH/runStats/$TASKHICUP/"$n"_uniques_cis-trans.png\"><img src=\"$PROJECT_RELPATH/runStats/$TASKHICUP/"$n"_uniques_cis-trans.png\" width=\"200px\"/></a></td>"

    done
    echo "<table><tr>$row0</tr><tr>$row1</tr><tr>$row2</tr></table>" >> $SUMMARYTMP
    
    summaryFooter "$TASKHICUP" "$SUMMARYTMP"
fi

################################################################################
if [ -n "$RUNHOMERCHIPSEQ" ];then
    summaryHeader "Homer ChIP-Seq" "$TASKHOMERCHIPSEQ" "chipseqHomer.sh" "$SUMMARYTMP"

    python ${NGSANE_BASE}/core/Summary.py "$(gatherDirs $TASKHOMERCHIPSEQ)" ".summary.txt" homerchipseq >> $SUMMARYTMP

    summaryFooter "$TASKHOMERCHIPSEQ" "$SUMMARYTMP"
fi

################################################################################
if [ -n "$RUNPEAKRANGER" ];then
    summaryHeader "Peakranger" "$TASKPEAKRANGER" "peakranger.sh" "$SUMMARYTMP"

    python ${NGSANE_BASE}/core/Summary.py "$(gatherDirs $TASKPEAKRANGER)" ".summary.txt" peakranger >> $SUMMARYTMP

    summaryFooter "$TASKPEAKRANGER" "$SUMMARYTMP"
fi

################################################################################
if [ -n "$RUNMACS2" ];then
    summaryHeader "MACS2" "$TASKMACS2" "macs2.sh" "$SUMMARYTMP"

    vali=$(gatherDirs $TASKMACS2)
    python ${NGSANE_BASE}/core/Summary.py "$vali" ".summary.txt" macs2 >> $SUMMARYTMP

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

    summaryFooter "$TASKMACS2" "$SUMMARYTMP"
fi

################################################################################
if [ -n "$RUNMEMECHIP" ];then
    summaryHeader "MEME-chip Motif discovery" "$TASKMEMECHIP" "memechip.sh" "$SUMMARYTMP"

    vali=""
    CURDIR=$(pwd)
    for dir in ${DIR[@]}; do
        vali=$vali" $OUT/$dir/$TASKMEMECHIP/"
        cd $OUT/$dir/$TASKMEMECHIP
        for d in $(find . -maxdepth 1 -mindepth 1 -type d -exec basename '{}' \; ); do
                echo "<a href=\"$PROJECT_RELPATH/$dir/$TASKMEMECHIP/$d/index.html\">$dir/$d</a> " >> $CURDIR/$SUMMARYTMP
        done
    done
    cd $CURDIR
    python ${NGSANE_BASE}/core/Summary.py "$vali" ".summary.txt" memechip >> $SUMMARYTMP

    summaryFooter "$TASKMEMECHIP" "$SUMMARYTMP"
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


################################################################################
echo '<!DOCTYPE html><html><title>NGSANE project card</title><head>' > $SUMMARYFILE.tmp
cat ${NGSANE_BASE}/core/Summary.css >> $SUMMARYFILE.tmp

echo "<script type='text/javascript'>" >> $SUMMARYFILE.tmp
cat ${NGSANE_BASE}/core/jquery-1.9.1.min.js >> $SUMMARYFILE.tmp
echo '''</script></head><body>

<div id="center">
''' >> $SUMMARYFILE.tmp
echo "<div class='panel' id='quicklinks'><h2>Quicklinks</h2><div>" >> $SUMMARYFILE.tmp
declare -a LINKSET=( )
for i in $LINKS; do
    LINKSET=("${LINKSET[@]}" "<a href='#$i'>$i</a>")
done
echo $(IFS='|' ; echo "${LINKSET[*]}") >> $SUMMARYFILE.tmp

echo "</div><!-- Links --></div><!-- panel -->" >>$SUMMARYFILE.tmp

echo "<hr><span>Report generated with "`trigger.sh -v`"</span><span style='float:right;'>Last modified: "`date`"</span>" >> $SUMMARYTMP
echo "</div><!-- center --></body></html>" >> $SUMMARYTMP
################################################################################
cat $SUMMARYFILE.tmp $SUMMARYTMP > $SUMMARYFILE

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
