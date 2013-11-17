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

PROJECT_RELPATH=$(python -c "import os.path; print os.path.relpath('$(pwd -P)',os.path.realpath('$(dirname $SUMMARYTMP)'))")
[ -z "$PROJECT_RELPATH" ] && PROJECT_RELPATH="."


################################################################################
# define functions for generating summary scaffold
#
# summaryHeader takes 4 parameters
# $1=PIPELINE name
# $2=TASK (e.g. $TASKBWA)
# $3=pipeline mod script (e.g. bwa.sh)
# $4=html report file ($SUMMARYTMP)
# $5=result file suffix (e.g. asb.bwa)
# $6=output location if different to the default ($TASK)
function summaryHeader {
    LINKS=$LINKS" $2"   
    echo "<div class='panel' id='$2_panel'>
        <div class='headbagb' id='$2_panelback'><a name='$2'></a>
            <h2 id='$2_h_results' class='sub'>$1</h2>
            <h2 id='$2_h_checklist' class='sub inactive' rel='checklist'>Checklist<span class='counter'><span class='passed' id='$2_counter_checkpoints_passed'></span><span class='failed' id='$2_counter_checkpoints_failed'></span></span></h2>
            <h2 id='$2_h_notes' class='sub inactive' rel='notes'>Notes<span class='counter'><span class='neutral' id='$2_counter_notes'></span></span></span></h2>
            <h2 id='$2_h_errors' class='sub inactive' rel='notes'>Errors<span class='counter'><span class='errors' id='$2_counter_errors'></span></span></h2>
            <h2 id='$2_h_logfiles' class='sub inactive' rel='errors'>Log files</h2>" >> $4
    if [ -n "$5" ]; then 
        SUFFIX="--filesuffix $5"; 
        echo "<h2 id='$2_h_files' class='sub inactive' rel='files'>Result files</h2>" >> $4
    fi
    echo "</div><div class='wrapper'><div class='hidden'>" >> $4
    echo "QC - $2"
    #check of resultfiles are redirected into to a different folder
    if [ -n "$6" ]; then 
       RESULTLOCATION="--results-task $6"
    else 
        RESULTLOCATION=""
    fi
    ${NGSANE_BASE}/core/QC.sh --results-dir $OUT --html-file $4 --modscript ${NGSANE_BASE}/mods/$3 --log $QOUT --task $2 $RESULTLOCATION $SUFFIX >> $4    
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

# add bamannotate section
# $1=vali 
# $2=TASK (e.g.TASKBWA)
# $3=output file ($SUMMARYTMP)
function bamAnnotate {
	echo "<h3 class='overall'>Reads overlapping annotated regions</h3>" >>$3
	python ${NGSANE_BASE}/core/Summary.py ${1} .anno.stats annostats >> $3
	BAMANNOUT=runStats/bamann/$(echo ${DIR[@]}|sed 's/ /_/g')_${2}.ggplot
	BAMANNIMAGE=${BAMANNOUT/ggplot/pdf}
	if [ ! -f $BAMANNOUT ]; then mkdir -p $( dirname $BAMANNOUT); fi
	
	find ${1} -type f -name *anno.stats |  xargs -d"\n" cat | head -n 1 | gawk '{print "type "$0" sample"}' > $BAMANNOUT
    for i in $(find ${1} -type f -name *anno.stats); do
        name=$(basename $i)
        arrIN=(${name//.$ASD/ })
        grep --no-messages sum $i | gawk -v x=${arrIN[0]} '{print $0" "x}';
	done >> $BAMANNOUT
	sed -i -r 's/\s+/ /g' $BAMANNOUT
	Rscript ${NGSANE_BASE}/tools/bamann.R $BAMANNOUT $BAMANNIMAGE "Genome Features ${2}"
	convert $BAMANNIMAGE ${BAMANNIMAGE/pdf/jpg}
	echo "<h3>Annotation of mapped reads</h3>" >> $3
	echo "<div><a href=$PROJECT_RELPATH/$BAMANNIMAGE><img src=\""$PROJECT_RELPATH/${BAMANNIMAGE/.pdf/}"-0.jpg\" width='250px' style='float:left;'><img src=\""$PROJECT_RELPATH/${BAMANNIMAGE/.pdf/}"-1.jpg\" width='250px' style='float:left;'></a></div>">>$3

#	    python ${NGSANE_BASE}/tools/makeBamHistogram.py "${PROJECT}" $ROUTH >>$3

}

################################################################################
if [ -n "$RUNFASTQC" ]; then
    PIPELINE="FASTQC"
    PIPELINK="fastqc"
    
    LINKS=$LINKS" $PIPELINK"
    echo "<div class='panel'><div class='headbagb'><a name='$PIPELINK'><h2 class='sub'>$PIPELINE</h2></a></div>" >>$SUMMARYTMP

    for dir in ${DIR[@]}; do
    echo "<h3>${OUT##*/}/fastq/$dir/</h3>" >>$SUMMARYTMP
        echo "<table class='data'>" >>$SUMMARYTMP
        echo "<thead><tr><th><div style='width:100px'>Chart</div></th><th><div style='width:140px'>Encoding</div></th><th><div style='width:120px'>Library size</div></th><th><div style='width:50px'>Read</div></th><th><div style='width:80px'>Read length</div></th><th><div style='width:50px'>%GC</div></th><th><div style='width:120px'>Read qualities</th><th class='left'>Library</th></tr></thead><tbody>" >>$SUMMARYTMP
        
        if [[ -e $dir/$TASKFASTQC/ ]]; then
            for fasta in $(ls $OUT/fastq/$dir/*.$FASTQ); do
                fastan=${fasta##*/}
                f="$dir/$TASKFASTQC/${fastan/.$FASTQ/}_fastqc.zip"
                # get basename of f
                n=${f##*/}
                n=${n/"_fastqc.zip"/}
                ICO=" <img height='15px' class='noborder' style='vertical-align:middle' src='$PROJECT_RELPATH/$dir/$TASKFASTQC/"$n"_fastqc/Icons/"
                P=$(grep "PASS" -c $dir/$TASKFASTQC/$n"_fastqc/summary.txt")
                W=$(grep "WARN" -c $dir/$TASKFASTQC/$n"_fastqc/summary.txt")
                F=$(grep "FAIL" -c $dir/$TASKFASTQC/$n"_fastqc/summary.txt")
                CHART=$ICO"tick.png' title='$P'\>$P"
                if [ "$W" -ne "0" ]; then CHART=$CHART""$ICO"warning.png'\>"$W; fi
                if [ "$F" -ne "0" ]; then CHART=$CHART""$ICO"error.png'\>"$F; fi
                ENCODING=$(grep "Encoding" $dir/$TASKFASTQC/$n"_fastqc/fastqc_data.txt" | head -n 1 | cut -f 2)
                LIBRARYSIZE=$(grep "Total Sequences" $dir/$TASKFASTQC/$n"_fastqc/fastqc_data.txt" | head -n 1 | cut -f 2)
                READLENGTH=$(grep "Sequence length" $dir/$TASKFASTQC/$n"_fastqc/fastqc_data.txt" | head -n 1 | cut -f 2)
                GCCONTENT=$(grep "\%GC" $dir/$TASKFASTQC/$n"_fastqc/fastqc_data.txt" | head -n 1 | cut -f 2)
                if [[ "$f" == *$READTWO* ]] && [ "$f" != "${f/$READTWO/$READONE}" ]; then
                    READ=2
                else
                    READ=1
                fi
                echo "<tr style='vertical-align: middle;'><td>$CHART</td><td>$ENCODING</td><td>$LIBRARYSIZE</td><td>$READ</td><td>$READLENGTH</td><td>$GCCONTENT</td><td>" >>$SUMMARYTMP
        
                if [[ "$f" == *$READONE* ]]; then
                    echo "<a href='$PROJECT_RELPATH/$dir/$TASKFASTQC/${n}_fastqc/fastqc_report.html'><img src='$PROJECT_RELPATH/$dir/$TASKFASTQC/${n}_fastqc/Images/per_base_quality.png' height=75 alt='Quality scores for all first reads'/></a>" >>$SUMMARYTMP
                    
                    if [ -e ${f/$READONE/$READTWO} ] && [ "$f" != "${f/$READONE/$READTWO}" ]; then
                        echo "<a href='$PROJECT_RELPATH/$dir/$TASKFASTQC/${n/$READONE/$READTWO}_fastqc/fastqc_report.html'><img src='$PROJECT_RELPATH/$dir/$TASKFASTQC/${n/$READONE/$READTWO}_fastqc/Images/per_base_quality.png' height=75 alt='Quality scores for all second reads'/></a>" >>$SUMMARYTMP
                    fi
                
                elif [[ "$f" == *$READTWO* ]] && [ "$f" != "${f/$READTWO/$READONE}" ]; then
                    echo "<a href='$PROJECT_RELPATH/$dir/$TASKFASTQC/${n/$READTWO/$READONE}_fastqc/fastqc_report.html'><img src='$PROJECT_RELPATH/$dir/$TASKFASTQC/${n/$READTWO/$READONE}_fastqc/Images/per_base_quality.png' height=75 alt='Quality scores for all first reads'/></a>" >>$SUMMARYTMP
                    echo "<a href='$PROJECT_RELPATH/$dir/$TASKFASTQC/${n}_fastqc/fastqc_report.html'><img src='$PROJECT_RELPATH/$dir/$TASKFASTQC/${n}_fastqc/Images/per_base_quality.png' height=75 alt='Quality scores for all second reads'/></a>" >>$SUMMARYTMP
        		else
        			echo "[ERROR] no fastq files $f"
                fi
                echo "</td><td class='left'><a href='$PROJECT_RELPATH/$dir/$TASKFASTQC/"$n"_fastqc/fastqc_report.html'>${fastan/.$FASTQ/}</a></td></tr>" >>$SUMMARYTMP

            done
        fi
        echo "</tbody></table>">>$SUMMARYTMP
    done
    
    echo "</div></div>">>$SUMMARYTMP
fi

################################################################################
if [[ -n "$RUNFASTQSCREEN" ]]; then
    summaryHeader "FASTQ screen" "$TASKFASTQSCREEN" "fastqscreen.sh" "$SUMMARYTMP"

    vali=$(gatherDirs $TASKFASTQSCREEN)
    python ${NGSANE_BASE}/core/Summary.py "$vali" _screen.txt fastqscreen --noSummary --noOverallSummary >>$SUMMARYTMP

    imgs=""
    for dir in ${DIR[@]}; do
        for f in $(ls $dir/$TASKFASTQSCREEN/*_screen.png 2> /dev/null); do
            n=${f##*/}
            imgs+="<div class='inset_image'>${n/"_screen.png"/}<br/><a href=\"$PROJECT_RELPATH/$dir/$TASKFASTQSCREEN/${n}\"><img src=\"$PROJECT_RELPATH/$dir/$TASKFASTQSCREEN/$n\" width=\"200px\"/></a></div>"
        done
    done
    echo "<div>$imgs</div>" >> $SUMMARYTMP
    
    summaryFooter "$TASKFASTQSCREEN" "$SUMMARYTMP"
fi


################################################################################
if [[ -n "$RUNBLUE" ]]; then
    summaryHeader "BLUE error correction" "$TASKBLUE" "blue.sh" "$SUMMARYTMP"

    vali=""
    for dir in ${DIR[@]}; do
        vali=$vali" $OUT/fastq/${dir/_$TASKBLUE/}_$TASKBLUE/"
    done

    python ${NGSANE_BASE}/core/Summary.py "$vali" stats.txt blue >>$SUMMARYTMP

	mkdir -p runstats/blue
	BLUEOUT=runstats/blue/$(echo ${DIR[@]}|sed 's/ /_/g').ggplot
	IMAGE=runstats/blue/$(echo ${DIR[@]}|sed 's/ /_/g').pdf
	echo -e "copy\tcount\tvalue\tperc\tsample" > $BLUEOUT
	for i in $(ls $vali/tessel/*histo*); do
		name=$(basename $i)
		arrIN=(${name//$READONE/ })
		head -n -10 $i | tail -n +4 | gawk -v x=${arrIN[0]} '{print $0"\t"x}'; 
	done >> $BLUEOUT
	Rscript ${NGSANE_BASE}/tools/blue.R $BLUEOUT $IMAGE
	convert $IMAGE ${IMAGE/pdf/jpg}
	echo "<h3>Annotation of mapped reads</h3>" >> $SUMMARYTMP
	echo "<a href=$PROJECT_RELPATH/$IMAGE><img src=\""$PROJECT_RELPATH/${IMAGE/.pdf/}".jpg\"></a>">>$SUMMARYTMP

    summaryFooter "$TASKBLUE" "$SUMMARYTMP"
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
if [[ -n "$RUNMAPPINGBWA" || -n "$RUNMAPPINGBWA2" ]]; then
    summaryHeader "BWA mapping" "$TASKBWA" "bwa.sh" "$SUMMARYTMP"

    vali=$(gatherDirs $TASKBWA)
    python ${NGSANE_BASE}/core/Summary.py "$vali" .$ASD.bam.stats samstats >>$SUMMARYTMP

    if [ -n "$RUNANNOTATINGBAM" ]; then
        bamAnnotate "$vali" $TASKBWA  $SUMMARYTMP
    fi

    summaryFooter "$TASKBWA" "$SUMMARYTMP"
fi


################################################################################
if [[ -n "$RUNREALRECAL" || -n "$RUNREALRECAL2" || -n "$RUNREALRECAL3" ]]; then 
    summaryHeader "Recalibrate + Realign" "$TASKRCA" "reCalAln2.sh" "$SUMMARYTMP"

    python ${NGSANE_BASE}/core/Summary.py "$(gatherDirs $TASKRCA)" .$ASR".bam.stats" samstatsrecal >>$SUMMARYTMP

    summaryFooter "$TASKRCA" "$SUMMARYTMP"
fi

################################################################################
if [[ -n "$RUNMAPPINGBOWTIE" ]]; then
    summaryHeader "Bowtie v1 mapping" "$TASKBOWTIE" "bowtie.sh" "$SUMMARYTMP" ".$ASD.bam"

    vali=$(gatherDirs $TASKBOWTIE)
    python ${NGSANE_BASE}/core/Summary.py "$(gatherDirs $TASKBOWTIE)" .$ASD.bam.stats samstats >>$SUMMARYTMP

    if [ -n "$RUNANNOTATINGBAM" ]; then
        bamAnnotate "$vali" $TASKBOWTIE  $SUMMARYTMP
    fi
    
    summaryFooter "$TASKBOWTIE" "$SUMMARYTMP"
fi

################################################################################
if [[ -n "$RUNMAPPINGBOWTIE2" ]]; then
    summaryHeader "Bowtie v2 mapping" "$TASKBOWTIE2" "bowtie2.sh" "$SUMMARYTMP" ".$ASD.bam"

    python ${NGSANE_BASE}/core/Summary.py "$(gatherDirs $TASKBOWTIE2)" .$ASD.bam.stats samstats >>$SUMMARYTMP

    if [ -n "$RUNANNOTATINGBAM" ]; then
        bamAnnotate "$vali" $TASKBOWTIE2 $SUMMARYTMP
    fi
    
    summaryFooter "$TASKBOWTIE2" "$SUMMARYTMP"
fi

################################################################################
if [[ -n "$RUNTOPHAT" || -n "$RUNTOPHATCUFFHTSEQ" ]]; then
    summaryHeader "Tophat" "$TASKTOPHAT" "tophat.sh" "$SUMMARYTMP" ".$ASD.bam"

	vali=""
    echo "<br>Note, the duplication rate is not calculated by tophat and hence zero.<br>" >>$SUMMARYTMP
    CURDIR=$(pwd -P)
    for dir in ${DIR[@]}; do
    	vali=$vali" $OUT/$dir/$TASKTOPHAT/"
    	cd $OUT/$dir/$TASKTOPHAT
    	for d in $(find . -maxdepth 1 -mindepth 1 -type d -exec basename '{}' \; | grep --no-messages "RNASeQC"); do
            echo "<a href=\"$PROJECT_RELPATH/$dir/$TASKTOPHAT/$d/index.html\">RNAseq-QC for $dir/$d</a><br/>" >> $CURDIR/$SUMMARYTMP
		done
    done
    cd $CURDIR
    python ${NGSANE_BASE}/core/Summary.py "$vali" .$ASD.bam.stats tophat >>$SUMMARYTMP
    
    if [ -n "$RUNANNOTATINGBAM" ] || [ -n "$RUNTOPHATCUFFHTSEQ" ]; then
        bamAnnotate "$vali" $TASKTOPHAT  $SUMMARYTMP
    fi
    
    summaryFooter "$TASKTOPHAT" "$SUMMARYTMP"

fi

################################################################################
if [[ -n "$RUNCUFFLINKS" || -n "$RUNTOPHATCUFFHTSEQ" ]]; then
    summaryHeader "Cufflinks" "$TASKCUFFLINKS" "cufflinks.sh" "$SUMMARYTMP" "_transcripts.gtf"

    python ${NGSANE_BASE}/core/Summary.py "$(gatherDirs $TASKCUFFLINKS)" .summary.txt cufflinks >>$SUMMARYTMP

    summaryFooter "$TASKCUFFLINKS" "$SUMMARYTMP"
fi


################################################################################
if [[ -n "$RUNHTSEQCOUNT" || -n "$RUNTOPHATCUFFHTSEQ" ]]; then
    summaryHeader "Htseq-count" "$TASKHTSEQCOUNT" "htseqcount.sh" "$SUMMARYTMP" ".RPKM.csv"

    python ${NGSANE_BASE}/core/Summary.py "$(gatherDirs $TASKHTSEQCOUNT)" summary.txt htseqcount >>$SUMMARYTMP

    summaryFooter "$TASKHTSEQCOUNT" "$SUMMARYTMP"
fi

################################################################################
if [[ -n "$DEPTHOFCOVERAGE"  || -n "$DEPTHOFCOVERAGE2" ]]; then
    summaryHeader "Coverage" "$TASKVAR" "gatkSNPs.sh" "$SUMMARYTMP"

    vali=$(gatherDirs $TASKDOC)
    echo "<h3>Average coverage</h3>">>$SUMMARYTMP
    python ${NGSANE_BASE}/core/Summary.py "$vali" .$ASR".bam.doc.sample_summary" coverage >>$SUMMARYTMP
    echo "<h3>Base pair coverage over all intervals</h3>" >>$SUMMARYTMP
    python ${NGSANE_BASE}/core/Summary.py "$vali" .$ASR".bam.doc.sample_cumulative_coverage_counts" coverage --p >>$SUMMARYTMP
    echo "<h3>Intervals covered</h3>" >>$SUMMARYTMP
    python ${NGSANE_BASE}/core/Summary.py "$vali" .$ASR".bam.doc.sample_interval_statistics" coverage --p >>$SUMMARYTMP
    echo "<h3>On Target</h3>" >>$SUMMARYTMP
    python ${NGSANE_BASE}/core/Summary.py "$vali" .$ASR".bam.stats" target >>$SUMMARYTMP

    summaryFooter "$TASKVAR" "$SUMMARYTMP" 
fi

################################################################################
if [ -n "$RUNVARCALLS" ]; then 
    summaryHeader "Variant calling" "$TASKVAR" "gatkSNPs.sh" "$SUMMARYTMP"

	vali=$OUT/$TASKVAR/$(echo ${DIR[@]}|sed 's/ /_/g')/
    echo "<h3>SNPs</h3>">>$SUMMARYTMP
    python ${NGSANE_BASE}/core/Summary.py "$vali" "filter.snps.eval.txt" variant --n --l "../$PROJECT_RELPATH" >>$SUMMARYTMP
    python ${NGSANE_BASE}/core/Summary.py "$vali" "recalfilt.snps.eval.txt" variant --n --l "../$PROJECT_RELPATH" >>$SUMMARYTMP
    echo "<h3>INDELs</h3>" >>$SUMMARYTMP
    python ${NGSANE_BASE}/core/Summary.py "$vali" "filter.indel.eval.txt" variant --n --l "../$PROJECT_RELPATH" >>$SUMMARYTMP

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
    
    imgs=""
    for dir in ${DIR[@]}; do
        for f in $(ls $dir/$TASKHICUP/*_ditag_classification.png 2> /dev/null); do
            n=${f##*/}
            n=${n/"_ditag_classification.png"/}
            imgs+="<div class='inset_image'>$n<br/><a href=\"$PROJECT_RELPATH/runStats/$TASKHICUP/"$n"_ditag_classification.png\"><img src=\"$PROJECT_RELPATH/runStats/$TASKHICUP/"$n"_ditag_classification.png\" width=\"200px\"/></a><br/><a href=\"$PROJECT_RELPATH/runStats/$TASKHICUP/"$n"_uniques_cis-trans.png\"><img src=\"$PROJECT_RELPATH/runStats/$TASKHICUP/"$n"_uniques_cis-trans.png\" width=\"200px\"/></a></div>"
        done
    done
    echo "<div>$imgs</div>" >> $SUMMARYTMP
    
    summaryFooter "$TASKHICUP" "$SUMMARYTMP"
fi


################################################################################
if [ -n "$RUNCHANCE" ];then
    summaryHeader "Chance" "$TASKCHANCE" "chance.sh" "$SUMMARYTMP"

    vali=$(gatherDirs $TASKCHANCE)
    python ${NGSANE_BASE}/core/Summary.py "$vali" ".IPstrength" chance >> $SUMMARYTMP

    imgs=""
    for dir in ${DIR[@]}; do
        for f in $(ls $dir/$TASKCHANCE/*.png 2> /dev/null); do
            n=${f##*/}
            imgs+="<div class='inset_image'>${n/".png"/}<br/><a href=\"$PROJECT_RELPATH/$dir/$TASKCHANCE/${n/.png/.pdf}\"><img src=\"$PROJECT_RELPATH/$dir/$TASKCHANCE/$n\" width=\"200px\"/></a></div>"
        done
    done
    echo "<div>$imgs</div>" >> $SUMMARYTMP

    summaryFooter "$TASKCHANCE" "$SUMMARYTMP"
fi


################################################################################
if [ -n "$RUNBIGWIG" ];then
    summaryHeader "BigWig" "$TASKBIGWIG" "bigwig.sh" "$SUMMARYTMP" ".bw" $INPUT_BIGWIG

    summaryFooter "$TASKBIGWIG" "$SUMMARYTMP"
fi

################################################################################
if [ -n "$RUNHOMERCHIPSEQ" ];then
    summaryHeader "Homer ChIP-Seq" "$TASKHOMERCHIPSEQ" "chipseqHomer.sh" "$SUMMARYTMP" ".bed"

    python ${NGSANE_BASE}/core/Summary.py "$(gatherDirs $TASKHOMERCHIPSEQ)" ".summary.txt" homerchipseq >> $SUMMARYTMP

    summaryFooter "$TASKHOMERCHIPSEQ" "$SUMMARYTMP"
fi

################################################################################
if [ -n "$RUNPEAKRANGER" ];then
    summaryHeader "Peakranger" "$TASKPEAKRANGER" "peakranger.sh" "$SUMMARYTMP" "_region.bed"

    python ${NGSANE_BASE}/core/Summary.py "$(gatherDirs $TASKPEAKRANGER)" ".summary.txt" peakranger >> $SUMMARYTMP

    summaryFooter "$TASKPEAKRANGER" "$SUMMARYTMP"
fi

################################################################################
if [ -n "$RUNMACS2" ];then
    summaryHeader "MACS2" "$TASKMACS2" "macs2.sh" "$SUMMARYTMP" "_refinepeak.bed"

    vali=$(gatherDirs $TASKMACS2)
    python ${NGSANE_BASE}/core/Summary.py "$vali" ".summary.txt" macs2 >> $SUMMARYTMP

    imgs=""
    for dir in ${DIR[@]}; do
        for f in $(ls $dir/$TASKMACS2/*_model-0.png 2> /dev/null); do
            n=${f##*/}
            imgs+="<div class='inset_image'>${n/_model-0.png/}<br/><a href=\"$PROJECT_RELPATH/$dir/$TASKMACS2/${n/-0.png/.pdf}\"><img src=\"$PROJECT_RELPATH/$dir/$TASKMACS2/$n\" width=\"200px\"/><img src=\"$PROJECT_RELPATH/$dir/$TASKMACS2/${n/model-0.png/model-1.png}\" width=\"200px\"/></a></div>"
        done
    done
    echo "<div>$imgs</div>" >> $SUMMARYTMP

    summaryFooter "$TASKMACS2" "$SUMMARYTMP"
fi

################################################################################
if [ -n "$RUNMEMECHIP" ]; then
    summaryHeader "MEME-chip Motif discovery" "$TASKMEMECHIP" "memechip.sh" "$SUMMARYTMP"

    for dir in ${DIR[@]}; do
        if [ ! -d $dir/$TASKMEMECHIP ]; then
            continue
        fi
        
        echo "<h3>${OUT##*/}/$dir/$TASKMEMECHIP/</h3>" >>$SUMMARYTMP
        echo "<table class='data'>" >>$SUMMARYTMP
        echo "<thead><tr><th><div style='width:200px'>Logo</div></th><th><div style='width:140px'>Consensus motif</div></th><th><div style='width:80px'>q-value</div></th><th><div style='width:100px'>Similar to</div></th><th><div style='width:120px'>Peaks</div></th><th><div style='width:120px'>With strong sites</div></th><th><div style='width:40px'>%</div></th><th><div style='width:120px'>With weak sites</div><th><div style='width:40px'>%</div></th><th class='left'>Library</th></thead><tbody>" >>$SUMMARYTMP

        
        for summary in $(ls $OUT/$dir/$TASKMEMECHIP/*.summary.txt); do
            SAMPLE=${summary##*/}
            SAMPLE=${SAMPLE/.summary.txt/}

            MEMEMOTIF=$(grep "Query consensus:" $summary | cut -d':' -f 2)
            MEMEQVALUE=$(grep "Q-value:" $summary | cut -d':' -f 2)
            TOMTOMKNOWNMOTIF=$(grep "Most similar known motif:" $summary | cut -d':' -f 2)
            PEAKS=$(grep "Peak regions:" $summary | cut -d':' -f 2)
            FIMODIRECT=$(grep "bound directly" $summary | cut -d':' -f 2)
            FIMODIRECTP=$(echo "scale=2;100 * $FIMODIRECT / $PEAKS" | bc)
            FIMOINDIRECT=$(grep "bound indirectly" $summary | cut -d':' -f 2)
            FIMOINDIRECTP=$(echo "scale=2;100 * $FIMOINDIRECT / $PEAKS" | bc)
            echo "<tr style='vertical-align: middle;'>"  >>$SUMMARYTMP
            echo "<td><a href='$PROJECT_RELPATH/${dir/$OUT/}/$TASKMEMECHIP/$SAMPLE/index.html'><img src='$PROJECT_RELPATH/${dir/$OUT/}/$TASKMEMECHIP/$SAMPLE/meme_out/logo1.png' height=75 alt='Meme Motif LOGO'/></a></td>" >>$SUMMARYTMP            
            echo "<td>$MEMEMOTIF</td><td>$MEMEQVALUE</td><td>$TOMTOMKNOWNMOTIF</td><td>$PEAKS</td><td>$FIMODIRECT</td><td>$FIMODIRECTP</td><td>$FIMOINDIRECT</td><td>$FIMOINDIRECTP</td><td class='left'><a href='$PROJECT_RELPATH/${dir/$OUT/}/$TASKMEMECHIP/$SAMPLE/index.html'>$SAMPLE</a></td></tr>" >>$SUMMARYTMP
       
        done
        echo "</tbody></table>">>$SUMMARYTMP
    done

    summaryFooter "$TASKMEMECHIP" "$SUMMARYTMP"
fi


################################################################################
if [ -n "$RUNTRINITY" ] || [ -n "$RUNINCHWORM" ];then
    summaryHeader "Trinity - Inchworm" "$TASKINCHWORM" "trinity_inchworm.sh" "$SUMMARYTMP"

    python ${NGSANE_BASE}/core/Summary.py "$(gatherDirs $TASKINCHWORM)" .summary.txt "trinity_inchworm" --noSummary  >>$SUMMARYTMP

    summaryFooter "$TASKINCHWORM" "$SUMMARYTMP"
fi  
if [ -n "$RUNTRINITY" ] || [ -n "$RUNCHRYSALIS" ];then
    summaryHeader "Trinity - Chrysalis" "$TASKCHRYSALIS" "trinity_chrysalis.sh" "$SUMMARYTMP"

    python ${NGSANE_BASE}/core/Summary.py "$(gatherDirs $TASKCHRYSALIS)" .summary.txt "trinity_chrysalis" --noSummary  >>$SUMMARYTMP

    summaryFooter "$TASKCHRYSALIS" "$SUMMARYTMP"
fi    
if [ -n "$RUNTRINITY" ] || [ -n "$RUNBUTTERFLY" ];then
    summaryHeader "Trinity - Butterfly" "$TASKBUTTERFLY" "trinity_butterfly.sh" "$SUMMARYTMP"

    python ${NGSANE_BASE}/core/Summary.py "$(gatherDirs $TASKBUTTERFLY)" .summary.txt "trinity_butterfly" --noSummary  >>$SUMMARYTMP

    summaryFooter "$TASKBUTTERFLY" "$SUMMARYTMP"
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
	P=$(pwd -P)
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
