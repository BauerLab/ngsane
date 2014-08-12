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


################################################################################
if [[ -n "$RUNFASTQSCREEN" ]]; then
    NGSANE_REPORT_HEADER "FASTQ screen" "$TASK_FASTQSCREEN" "fastqscreen.sh" "$SUMMARYTMP"

    vali=$(NGSANE_REPORT_GATHERDIRS $TASK_FASTQSCREEN)
    python ${NGSANE_BASE}/core/Summary.py "$TASK_FASTQSCREEN" "$vali" _screen.txt fastqscreen --noSummary --noOverallSummary >>$SUMMARYTMP

    imgs=""
    for dir in ${DIR[@]}; do
        for f in $(ls ${dir%%/*}/$TASK_FASTQSCREEN/*_screen.png 2> /dev/null); do
            n=${f##*/}
            imgs+="<div class='inset_image'>${n/"_screen.png"/}<br/><a href=\"$PROJECT_RELPATH/${dir%%/*}/$TASK_FASTQSCREEN/${n}\"><img src=\"$PROJECT_RELPATH/${dir%%/*}/$TASK_FASTQSCREEN/$n\" width=\"200px\"/></a></div>"
        done
    done
    echo "<div>$imgs</div>" >> $SUMMARYTMP
    
    NGSANE_REPORT_FOOTER "$SUMMARYTMP"
fi


################################################################################
if [[ -n "$RUNBLUE" ]]; then
    NGSANE_REPORT_HEADER "BLUE error correction" "$TASK_BLUE" "blue.sh" "$SUMMARYTMP"

    vali=""
    for dir in ${DIR[@]}; do
        vali=$vali" $OUT/fastq/${dir}"_"$TASK_BLUE/"
    done

    python ${NGSANE_BASE}/core/Summary.py "$TASK_BLUE" "$vali" stats.txt blue >>$SUMMARYTMP

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

    NGSANE_REPORT_FOOTER "$SUMMARYTMP"
fi


################################################################################
if [ -n "$RUNTRIMGALORE" ];then
    NGSANE_REPORT_HEADER "Trimgalore trimming" "$TASK_TRIMGALORE" "trimgalore.sh" "$SUMMARYTMP"

    vali=""
    for dir in ${DIR[@]}; do
        vali=$vali" $OUT/fastq/${dir}"_"$TASK_TRIMGALORE/"
    done
    python ${NGSANE_BASE}/core/Summary.py "$TASK_TRIMGALORE" "$vali" "_trimming_report.txt" trimgalore --noSummary >> $SUMMARYTMP

    NGSANE_REPORT_FOOTER "$SUMMARYTMP"
fi

################################################################################
if [ -n "$RUNTRIMMOMATIC" ];then
    NGSANE_REPORT_HEADER "Trimmomatic trimming" "$TASK_TRIMMOMATIC" "trimmomatic.sh" "$SUMMARYTMP"

    vali=""
    for dir in ${DIR[@]}; do
        vali=$vali" $OUT/fastq/${dir}"_"$TASK_TRIMMOMATIC/"
    done
    python ${NGSANE_BASE}/core/Summary.py "$TASK_TRIMMOMATIC" "$vali" ".log" trimmomatic --noSummary >> $SUMMARYTMP

    NGSANE_REPORT_FOOTER "$SUMMARYTMP"
fi

################################################################################
if [ -n "$RUNCUTADAPT" ];then
    NGSANE_REPORT_HEADER "Cutadapt trimming" "$TASK_CUTADAPT" "cutadapt.sh" "$SUMMARYTMP"

    vali=""
    for dir in ${DIR[@]}; do
        vali=$vali" $OUT/fastq/${dir}"_"$TASK_CUTADAPT/"
    done
    python ${NGSANE_BASE}/core/Summary.py "$TASK_CUTADAPT" "$vali" ".stats" cutadapt --noSummary >> $SUMMARYTMP

    NGSANE_REPORT_FOOTER "$SUMMARYTMP"
fi


################################################################################
if [[ -n "$RUNBWA" ]]; then
    NGSANE_REPORT_HEADER "BWA mapping" "$TASK_BWA" "bwa.sh" "$SUMMARYTMP"

    vali=$(NGSANE_REPORT_GATHERDIRS $TASK_BWA)
    python ${NGSANE_BASE}/core/Summary.py "$TASK_BWA" "$vali" $ASD.bam.stats samstats >>$SUMMARYTMP

    if [ -n "$RUNANNOTATINGBAM" ]; then
        bamAnnotate "$vali" $TASK_BWA  $SUMMARYTMP
    fi

    NGSANE_REPORT_FOOTER "$SUMMARYTMP"
fi


################################################################################
if [[ -n "$RUNREALRECAL" ]]; then 
    NGSANE_REPORT_HEADER "Recalibrate + Realign" "$TASK_RECAL" "reCalAln.sh" "$SUMMARYTMP"

    python ${NGSANE_BASE}/core/Summary.py "$TASK_RECAL" "$(NGSANE_REPORT_GATHERDIRS $TASK_RECAL)" $ASR.bam.stats samstatsrecal >>$SUMMARYTMP

    NGSANE_REPORT_FOOTER "$SUMMARYTMP"
fi

################################################################################
if [[ -n "$RUNBOWTIE" ]]; then
    NGSANE_REPORT_HEADER "Bowtie v1 mapping" "$TASK_BOWTIE" "bowtie.sh" "$SUMMARYTMP" "$ASD.bam"

    vali=$(NGSANE_REPORT_GATHERDIRS $TASK_BOWTIE)
    python ${NGSANE_BASE}/core/Summary.py "$TASK_BOWTIE" "$(NGSANE_REPORT_GATHERDIRS $TASK_BOWTIE)" $ASD.bam.stats samstats >>$SUMMARYTMP

    if [ -n "$RUNANNOTATINGBAM" ]; then
        bamAnnotate "$vali" $TASK_BOWTIE  $SUMMARYTMP
    fi
    
    NGSANE_REPORT_FOOTER "$SUMMARYTMP"
fi

################################################################################
if [[ -n "$RUNBOWTIE2" ]]; then
    NGSANE_REPORT_HEADER "Bowtie v2 mapping" "$TASK_BOWTIE2" "bowtie2.sh" "$SUMMARYTMP" "$ASD.bam"

    python ${NGSANE_BASE}/core/Summary.py "$TASK_BOWTIE2" "$(NGSANE_REPORT_GATHERDIRS $TASK_BOWTIE2)" $ASD.bam.stats samstats >>$SUMMARYTMP

    if [ -n "$RUNANNOTATINGBAM" ]; then
        bamAnnotate "$vali" $TASK_BOWTIE2 $SUMMARYTMP
    fi
    
    NGSANE_REPORT_FOOTER "$SUMMARYTMP"
fi

################################################################################
if [[ -n "$RUNTOPHAT" || -n "$RUNTOPHATCUFFHTSEQ" ]]; then
    NGSANE_REPORT_HEADER "Tophat" "$TASK_TOPHAT" "tophat.sh" "$SUMMARYTMP" "$ASD.bam"

	vali=$(NGSANE_REPORT_GATHERDIRS $TASK_TOPHAT)
    echo "<br>[NOTE] the duplication rate is not calculated by tophat and hence zero.<br>" >>$SUMMARYTMP
    python ${NGSANE_BASE}/core/Summary.py "$TASK_TOPHAT" "$vali" $ASD.bam.stats tophat >>$SUMMARYTMP
    
    if [ -n "$RUNANNOTATINGBAM" ]; then
        bamAnnotate "$vali" $TASK_TOPHAT  $SUMMARYTMP
    fi
    
    NGSANE_REPORT_FOOTER "$SUMMARYTMP"

fi

################################################################################
if [[ -n "$RUNRNASEQC" ]]; then
    NGSANE_REPORT_HEADER "RNA-SeQC" "$TASK_RNASEQC" "rnaseqc.sh" "$SUMMARYTMP" "$ASD.bam"

	vali=""
    CURDIR=$(pwd -P)
    for dir in ${DIR[@]}; do
    	vali=$vali" $OUT/${dir%%/*}/$TASK_RNASEQC/"
    	cd $OUT/${dir%%/*}/$TASK_RNASEQC
    	for d in $(find . -maxdepth 1 -mindepth 1 -type d -exec basename '{}' \; ); do
            echo "<a href=\"$PROJECT_RELPATH/${dir%%/*}/$TASK_RNASEQC/$d/index.html\">RNAseq-QC for ${dir%%/*}/$d</a><br/>" >> $CURDIR/$SUMMARYTMP
		done
    done
    cd $CURDIR
    python ${NGSANE_BASE}/core/Summary.py "$TASK_RNASEQC" "$vali" $ASD.bam.stats tophat >>$SUMMARYTMP
    
    NGSANE_REPORT_FOOTER "$SUMMARYTMP"
fi

################################################################################
if [[ -n "$RUNCUFFLINKS" || -n "$RUNTOPHATCUFFHTSEQ" ]]; then
    NGSANE_REPORT_HEADER "Cufflinks" "$TASK_CUFFLINKS" "cufflinks.sh,cuffpost.sh" "$SUMMARYTMP" "_transcripts.gtf"

    python ${NGSANE_BASE}/core/Summary.py "$TASK_CUFFLINKS" "$(NGSANE_REPORT_GATHERDIRS $TASK_CUFFLINKS)" .summary.txt cufflinks >>$SUMMARYTMP

    NGSANE_REPORT_FOOTER "$SUMMARYTMP"
fi

################################################################################
if [[ -n "$RUNHTSEQCOUNT" || -n "$RUNTOPHATCUFFHTSEQ" ]]; then
    NGSANE_REPORT_HEADER "Htseq-count" "$TASK_HTSEQCOUNT" "htseqcount.sh,countsTable.sh" "$SUMMARYTMP"

    python ${NGSANE_BASE}/core/Summary.py "$TASK_HTSEQCOUNT" "$(NGSANE_REPORT_GATHERDIRS $TASK_HTSEQCOUNT)" .summary.txt htseqcount >>$SUMMARYTMP

    NGSANE_REPORT_FOOTER "$SUMMARYTMP"
fi

################################################################################
if [[ -n "$DEPTHOFCOVERAGE"  || -n "$DEPTHOFCOVERAGE2" ]]; then
    NGSANE_REPORT_HEADER "Coverage" "$TASK_VAR" "gatkSNPs.sh" "$SUMMARYTMP"

    vali=$(NGSANE_REPORT_GATHERDIRS $TASK_GATKDOC)
    echo "<h3>Average coverage</h3>">>$SUMMARYTMP
    python ${NGSANE_BASE}/core/Summary.py "$TASK_VAR" "$vali" $ASR".bam.doc.sample_summary" coverage >>$SUMMARYTMP
    echo "<h3>Base pair coverage over all intervals</h3>" >>$SUMMARYTMP
    python ${NGSANE_BASE}/core/Summary.py "$TASK_VAR" "$vali" $ASR".bam.doc.sample_cumulative_coverage_counts" coverage --p >>$SUMMARYTMP
    echo "<h3>Intervals covered</h3>" >>$SUMMARYTMP
    python ${NGSANE_BASE}/core/Summary.py "$TASK_VAR" "$vali" $ASR".bam.doc.sample_interval_statistics" coverage --p >>$SUMMARYTMP
    echo "<h3>On Target</h3>" >>$SUMMARYTMP
    python ${NGSANE_BASE}/core/Summary.py "$TASK_VAR" "$vali" $ASR".bam.stats" target >>$SUMMARYTMP

    NGSANE_REPORT_FOOTER "$SUMMARYTMP" 
fi

################################################################################
if [ -n "$RUNVARCALLS" ]; then 
    NGSANE_REPORT_HEADER "Variant calling" "$TASK_GATKVAR" "gatkVARs.sh" "$SUMMARYTMP"

	vali=$OUT/$TASK_GATKVAR/$(echo ${DIR[@]}|sed 's/ /_/g')/
    echo "<h3>SNPs</h3>">>$SUMMARYTMP
    python ${NGSANE_BASE}/core/Summary.py "$TASK_GATKVAR" "$vali" "filter.snps.eval.txt" variant --n --l "../$PROJECT_RELPATH" >>$SUMMARYTMP
    python ${NGSANE_BASE}/core/Summary.py "$TASK_GATKVAR" "$vali" "recalfilt.snps.eval.txt" variant --n --l "../$PROJECT_RELPATH" >>$SUMMARYTMP
    echo "<h3>INDELs</h3>" >>$SUMMARYTMP
    python ${NGSANE_BASE}/core/Summary.py "$TASK_GATKVAR" "$vali" "filter.indel.eval.txt" variant --n --l "../$PROJECT_RELPATH" >>$SUMMARYTMP

    NGSANE_REPORT_FOOTER "$SUMMARYTMP"
fi

################################################################################
if [ -n "$RUNANNOTATION" ]; then
    NGSANE_REPORT_HEADER "Variant annotation" "$TASK_ANNOVAR" "annovar.sh" "$SUMMARYTMP"

    vali=""
    for dir in ${DIR[@]}; do
        vali=$vali" "$( ls $OUT/$TASK_ANNOVAR/${dir%%/*}/*.csv)
    done
    echo "<h3>Annotation Files</h3>">>$SUMMARYTMP
    echo "Please right click the link and choose \"Save as...\" to download.<br><br>">> $SUMMARYTMP
    for v in $vali; do
        name=`basename $v`
        address=${v/\/illumina/http:\/\/hpsc.csiro.au}
        echo "<a href=\"$address\">$name</a><br>" >> $SUMMARYTMP
    done
    echo "<br>More information about the columns can be found on the <a target=new href=\"http://redmine.qbi.uq.edu.au/knowledgebase/articles/12\">Project Server</a> (uqlogin). and the description of the <a href=\"http://www.broadinstitute.org/gsa/wiki/index.php/Understanding_the_Unified_Genotyper%27s_VCF_files\">vcf file</a>">> $SUMMARYTMP
    
    NGSANE_REPORT_FOOTER "$SUMMARYTMP"
fi

################################################################################
if [ -n "$RUNHICLIB" ];then
    NGSANE_REPORT_HEADER "HiClib" "$TASK_HICLIB" "hiclibMapping.sh" "$SUMMARYTMP"

    vali=$(NGSANE_REPORT_GATHERDIRS $TASK_HICLIB)
    python ${NGSANE_BASE}/core/Summary.py "$TASK_HICLIB" "$vali" ".log" hiclibMapping >> $SUMMARYTMP
    for dir in $vali; do
        for pdf in $(ls -f ${dir%%/*}/*.pdf 2>/dev/null ); do
            echo "<a href='$pdf'>${pdf##*/}</a> " >> $SUMMARYTMP
        done
    done

    NGSANE_REPORT_FOOTER "$SUMMARYTMP"
fi

################################################################################
if [ -n "$RUNHICUP" ];then
    NGSANE_REPORT_HEADER "HiCUP" "$TASK_HICUP" "hicup.sh" "$SUMMARYTMP"

    vali=$(NGSANE_REPORT_GATHERDIRS $TASK_HICUP)
    echo "<h4>truncater</h4>">>$SUMMARYTMP
    python ${NGSANE_BASE}/core/Summary.py "$TASK_HICUP" "$vali" "_truncater_summary.txt" hicup --noSummary --noOverallSummary >> $SUMMARYTMP
    echo "<h4>mapper</h4>">>$SUMMARYTMP
    python ${NGSANE_BASE}/core/Summary.py "$TASK_HICUP" "$vali" "_mapper_summary.txt" hicup --noSummary --noOverallSummary >> $SUMMARYTMP
    echo "<h4>filter</h4>">>$SUMMARYTMP
    python ${NGSANE_BASE}/core/Summary.py "$TASK_HICUP" "$vali" "_filter_summary.txt" hicup --noSummary --noOverallSummary >> $SUMMARYTMP
    echo "<h4>deduplicator</h4>">>$SUMMARYTMP
    python ${NGSANE_BASE}/core/Summary.py "$TASK_HICUP" "$vali" "_deduplicator_summary.txt" hicup --noSummary --noOverallSummary >> $SUMMARYTMP
    
    imgs=""
#    for f in $(ls runStats/$TASK_HICUP/*_ditag_classification.png 2> /dev/null); do
    for dir in ${DIR[@]}; do
        echo ${dir%%/*}
        for f in $(ls $OUT/${dir%%/*}/$TASK_HICUP/*_ditag_classification.png 2> /dev/null); do
            echo $f
            n=${f/"_ditag_classification.png"/}
            
            imgs+="<div class='inset_image'>${n##*/}<br/><a href=\"${n}_ditag_classification.png\"><img src=\"${n}_ditag_classification.png\" width=\"200px\"/></a><br/><a href=\"${n}_uniques_cis-trans.png\"><img src=\"${n}_uniques_cis-trans.png\" width=\"200px\"/></a><br/><a href=\"${n}_ditag_size_distribution.png\"><img src=\"${n}_ditag_size_distribution.png\" width=\"200px\"/></a></div>"
        done
    done
    echo "<div>$imgs</div>" >> $SUMMARYTMP
    
    NGSANE_REPORT_FOOTER "$SUMMARYTMP"
fi

################################################################################
if [ -n "$RUNFITHIC" ];then
    NGSANE_REPORT_HEADER "Fit-hi-c" "$TASK_FITHIC" "fithic.sh" "$SUMMARYTMP"

    python ${NGSANE_BASE}/core/Summary.py "$2" "$(NGSANE_REPORT_GATHERDIRS $TASK_FITHIC)" ".log" fithic >> $SUMMARYTMP
    
    NGSANE_REPORT_FOOTER "$SUMMARYTMP"
fi

################################################################################
if [ -n "$RUNCHANCE" ];then
    NGSANE_REPORT_HEADER "Chance" "$TASK_CHANCE" "chance.sh" "$SUMMARYTMP"

    vali=$(NGSANE_REPORT_GATHERDIRS $TASK_CHANCE)
    python ${NGSANE_BASE}/core/Summary.py "$TASK_CHANCE" "$vali" ".IPstrength" chance >> $SUMMARYTMP

    imgs=""
    for dir in ${DIR[@]}; do
        for f in $(ls ${dir%%/*}/$TASK_CHANCE/*.png 2> /dev/null); do
            n=${f##*/}
            imgs+="<div class='inset_image'>${n/".png"/}<br/><a href=\"$PROJECT_RELPATH/${dir%%/*}/$TASK_CHANCE/${n/.png/.pdf}\"><img src=\"$PROJECT_RELPATH/${dir%%/*}/$TASK_CHANCE/$n\" width=\"200px\"/></a></div>"
        done
    done
    echo "<div>$imgs</div>" >> $SUMMARYTMP

    NGSANE_REPORT_FOOTER "$SUMMARYTMP"
fi

################################################################################
if [ -n "$RUNBIGWIG" ];then
    NGSANE_REPORT_HEADER "BigWig" "$TASK_BIGWIG" "bigwig.sh" "$SUMMARYTMP" ".bw"

    python ${NGSANE_BASE}/core/Summary.py "$TASK_BIGWIG" "$(NGSANE_REPORT_GATHERDIRS $TASK_BIGWIG)" ".bw.stats" bigwig  --noSummary --noOverallSummary >> $SUMMARYTMP
    
    NGSANE_REPORT_FOOTER "$SUMMARYTMP"
fi

################################################################################
if [ -n "$RUNFSEQ" ];then
    NGSANE_REPORT_HEADER "Fseq" "$TASK_FSEQ" "fseq.sh" "$SUMMARYTMP" ".narrowPeak"

    python ${NGSANE_BASE}/core/Summary.py "$(NGSANE_REPORT_GATHERDIRS $TASK_FSEQ)" ".narrowPeak" fseq >> $SUMMARYTMP

    NGSANE_REPORT_FOOTER "$SUMMARYTMP"
fi
################################################################################
if [ -n "$RUNHOMERCHIPSEQ" ];then
    NGSANE_REPORT_HEADER "Homer ChIP-Seq" "$TASK_HOMERCHIPSEQ" "chipseqHomer.sh" "$SUMMARYTMP" ".bed"

    python ${NGSANE_BASE}/core/Summary.py "$TASK_HOMERCHIPSEQ" "$(NGSANE_REPORT_GATHERDIRS $TASK_HOMERCHIPSEQ)" ".summary.txt" homerchipseq >> $SUMMARYTMP

    NGSANE_REPORT_FOOTER "$SUMMARYTMP"
fi

################################################################################
if [ -n "$RUNPEAKRANGER" ];then
    NGSANE_REPORT_HEADER "Peakranger" "$TASK_PEAKRANGER" "peakranger.sh" "$SUMMARYTMP" "_region.bed"

    python ${NGSANE_BASE}/core/Summary.py "$TASK_PEAKRANGER" "$(NGSANE_REPORT_GATHERDIRS $TASK_PEAKRANGER)" ".summary.txt" peakranger >> $SUMMARYTMP

    NGSANE_REPORT_FOOTER "$SUMMARYTMP"
fi

################################################################################
if [ -n "$RUNMACS2" ];then
    NGSANE_REPORT_HEADER "MACS2" "$TASK_MACS2" "macs2.sh" "$SUMMARYTMP" "_refinepeak.bed"


    vali=$(NGSANE_REPORT_GATHERDIRS $TASK_MACS2)
    python ${NGSANE_BASE}/core/Summary.py "$TASK_MACS2" "$vali" ".summary.txt" macs2 >> $SUMMARYTMP

    imgs=""
    for dir in ${DIR[@]}; do
        for f in $(ls ${dir%%/*}/$TASK_MACS2/*_model-0.png 2> /dev/null); do
            n=${f##*/}
            imgs+="<div class='inset_image'>${n/_model-0.png/}<br/><a href=\"$PROJECT_RELPATH/${dir%%/*}/$TASK_MACS2/${n/-0.png/.pdf}\"><img src=\"$PROJECT_RELPATH/${dir%%/*}/$TASK_MACS2/$n\" width=\"200px\"/><img src=\"$PROJECT_RELPATH/${dir%%/*}/$TASK_MACS2/${n/model-0.png/model-1.png}\" width=\"200px\"/></a></div>"
        done
    done
    echo "<div>$imgs</div>" >> $SUMMARYTMP

    NGSANE_REPORT_FOOTER "$SUMMARYTMP"
fi

################################################################################
if [ -n "$RUNMEMECHIP" ]; then
    NGSANE_REPORT_HEADER "MEME-chip Motif discovery" "$TASK_MEMECHIP" "memechip.sh" "$SUMMARYTMP"

echo "<table id='memechip_id_Table' class='data'>" >>$SUMMARYTMP
echo "</table>" >>$SUMMARYTMP
echo "<script type=\"text/javascript\">
//<![CDATA[
memechip_Table_json ={
	\"bAutoWidth\": false,
    \"bFilter\": true,
	\"aoColumns\": [
        {\"sWidth\": \"15%\", \"sTitle\": \"Library\"},
		{\"sWidth\": \"10%\", \"sTitle\": \"Experiment\"},
		{\"sWidth\": \"15%\", \"sTitle\": \"Logo\"},
		{\"sWidth\": \"10%\", \"sTitle\": \"Consensus motif\"},
		{\"sWidth\": \"7.5%\", \"sTitle\": \"E-value\"},
		{\"sWidth\": \"7.5%\", \"sTitle\": \"Similar to\"},
		{\"sWidth\": \"10%\", \"sTitle\": \"TOMTOM q-value\"},
		{\"sWidth\": \"5%\", \"sTitle\": \"Peaks\"},
		{\"sWidth\": \"5%\", \"sTitle\": \"With strong sites\"},
		{\"sWidth\": \"5%\", \"sTitle\": \"%\"},
        {\"sWidth\": \"5%\", \"sTitle\": \"With weak sites\"},
		{\"sWidth\": \"5%\", \"sTitle\": \"%\"},
		
			],	
	\"aaData\": [" >>$SUMMARYTMP 

    for dir in ${DIR[@]}; do
        if [ ! -d ${dir%%/*}/$TASK_MEMECHIP ]; then
            continue
        fi
               
        for summary in $(ls $OUT/${dir%%/*}/$TASK_MEMECHIP/*.summary.txt); do
            SAMPLE=${summary##*/}
            SAMPLE=${SAMPLE/.summary.txt/}

	    EXPERIMENT=${dir%%/*}

            CONSENSUSMOTIF=$(grep "Query consensus:" $summary | cut -d':' -f 2 | tr -d ' ')
            MEMEEVALUE=$(grep "E-value:" $summary | cut -d':' -f 2 | tr -d ' ')
            
            TOMTOMQVALUE=$(grep "Q-value:" $summary | cut -d':' -f 2 | tr -d ' ')
            TOMTOMKNOWNMOTIF=$(grep "Most similar known motif:" $summary | cut -d':' -f 2)
            PEAKS=$(grep "Peak regions:" $summary | cut -d':' -f 2)
            FIMODIRECT=$(grep "bound directly" $summary | cut -d':' -f 2 | tr -d ' ')
            [ -n "$FIMODIRECT" ] && FIMODIRECTP=$(echo "scale=2;100 * $FIMODIRECT / $PEAKS" | bc) || FIMODIRECTP=""
            FIMOINDIRECT=$(grep "bound indirectly" $summary | cut -d':' -f 2 | tr -d ' ')
            [ -n "$FIMOINDIRECT" ] && FIMOINDIRECTP=$(echo "scale=2;100 * $FIMOINDIRECT / $PEAKS" | bc) || FIMOINDIRECTP=""
            MEMEMOTIF=$(grep "MEME motif:" $summary | cut -d':' -f 2 | tr -d ' ')
            #echo "<tr style='vertical-align: middle;'>"  >>$SUMMARYTMP
            
            if [ -n "$MEMEMOTIF" ]; then
                #echo "<td><a href='$PROJECT_RELPATH/${dir/$OUT/}/$TASK_MEMECHIP/$SAMPLE/index.html'><img src='$PROJECT_RELPATH/${dir/$OUT/}/$TASK_MEMECHIP/$SAMPLE/meme_out/logo$MEMEMOTIF.png' height=75 alt='Meme Motif LOGO'/></a></td>" >>$SUMMARYTMP

		MEMECHIP=$MEMECHIP"[[\"<a href='$PROJECT_RELPATH/${dir/$OUT/}/$TASK_MEMECHIP/$SAMPLE/index.html'>'$SAMPLE'</a>\"], '$EXPERIMENT', [\"<a href='$PROJECT_RELPATH/${dir/$OUT/}/$TASK_MEMECHIP/$SAMPLE/index.html'><img src='$PROJECT_RELPATH/${dir/$OUT/}/$TASK_MEMECHIP/$SAMPLE/meme_out/logo$MEMEMOTIF.png' height=75 alt='Meme Motif LOGO'/></a>\"], '$CONSENSUSMOTIF', '$MEMEEVALUE', '$TOMTOMKNOWNMOTIF', '$TOMTOMQVALUE', $PEAKS, $FIMODIRECT, $FIMODIRECTP, '$FIMOINDIRECT', $FIMOINDIRECTP] "
            else
                #echo "<td><a href='$PROJECT_RELPATH/${dir/$OUT/}/$TASK_MEMECHIP/$SAMPLE/index.html'><div style='height:75px; width:150px; background-color:#ffffff;border: 1px dotted #999999; '></div></a></td>" >>$SUMMARYTMP

		MEMECHIP=$MEMECHIP"[[\"<a href='$PROJECT_RELPATH/${dir/$OUT/}/$TASK_MEMECHIP/$SAMPLE/index.html'>'$SAMPLE'</a>\"], '$TASK', [\"<a href='$PROJECT_RELPATH/${dir/$OUT/}/$TASK_MEMECHIP/$SAMPLE/index.html'></a>\"], '$CONSENSUSMOTIF', '$MEMEEVALUE', '$TOMTOMKNOWNMOTIF', '$TOMTOMQVALUE', $PEAKS, $FIMODIRECT, $FIMODIRECTP, $FIMOINDIRECT, $FIMOINDIRECTP] "
            fi
            #echo "<td>$CONSENSUSMOTIF</td><td>$MEMEEVALUE</td><td>$TOMTOMKNOWNMOTIF</td><td>$TOMTOMQVALUE</td><td>$PEAKS</td><td>$FIMODIRECT</td><td>$FIMODIRECTP</td><td>$FIMOINDIRECT</td><td>$FIMOINDIRECTP</td><td class='left'><a href='$PROJECT_RELPATH/${dir/$OUT/}/$TASK_MEMECHIP/$SAMPLE/index.html'>$SAMPLE</a></td></tr>" >>$SUMMARYTMP
       
        done
       
    done

MEMECHIP=${MEMECHIP//] [/], 
[}
    echo $MEMECHIP >>$SUMMARYTMP
    echo "]} </script>" >>$SUMMARYTMP
    NGSANE_REPORT_FOOTER "$SUMMARYTMP"
fi


################################################################################
if [ -n "$RUNTRINITY" ] || [ -n "$RUNINCHWORM" ];then
    NGSANE_REPORT_HEADER "Trinity - Inchworm" "$TASK_INCHWORM" "trinity_inchworm.sh" "$SUMMARYTMP"

    python ${NGSANE_BASE}/core/Summary.py "$TASK_INCHWORM" "$(NGSANE_REPORT_GATHERDIRS $TASK_INCHWORM)" .summary.txt "trinity_inchworm" --noSummary  >>$SUMMARYTMP

    NGSANE_REPORT_FOOTER "$SUMMARYTMP"
fi  
if [ -n "$RUNTRINITY" ] || [ -n "$RUNCHRYSALIS" ];then
    NGSANE_REPORT_HEADER "Trinity - Chrysalis" "$TASK_CHRYSALIS" "trinity_chrysalis.sh" "$SUMMARYTMP"

    python ${NGSANE_BASE}/core/Summary.py "$TASK_INCHWORM" "$(NGSANE_REPORT_GATHERDIRS $TASK_CHRYSALIS)" .summary.txt "trinity_chrysalis" --noSummary  >>$SUMMARYTMP

    NGSANE_REPORT_FOOTER "$SUMMARYTMP"
fi    
if [ -n "$RUNTRINITY" ] || [ -n "$RUNBUTTERFLY" ];then
    NGSANE_REPORT_HEADER "Trinity - Butterfly" "$TASK_BUTTERFLY" "trinity_butterfly.sh" "$SUMMARYTMP"

    python ${NGSANE_BASE}/core/Summary.py "$TASK_INCHWORM" "$(NGSANE_REPORT_GATHERDIRS $TASK_BUTTERFLY)" .summary.txt "trinity_butterfly" --noSummary  >>$SUMMARYTMP

    NGSANE_REPORT_FOOTER "$SUMMARYTMP"
fi 


################################################################################
# pindel
################################################################################
if [ -n "$RUNPINDEL" ]; then 
    NGSANE_REPORT_HEADER "Structural Variants" "$INPUT_PINDEL-$TASK_PINDEL" "pindel.sh,variantcollect.sh" "$SUMMARYTMP" 
#$OUT/variant/${INPUT_PINDEL}-${TASK_PINDEL}-$(echo ${DIR[@]}|sed 's/ /_/g')/ "joined.eval.txt"

	vali=$OUT/variant/${INPUT_PINDEL}-${TASK_PINDEL}-$(echo ${DIR[@]}|sed 's/ /_/g')/
    echo "<h3>Variants</h3>">>$SUMMARYTMP
    python ${NGSANE_BASE}/core/Summary.py "$INPUT_PINDEL-$TASK_PINDEL" "$vali" "eval.txt" variant --n --l "../$PROJECT_RELPATH" >>$SUMMARYTMP

    NGSANE_REPORT_FOOTER "$SUMMARYTMP"
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



################################################################################
# citation list
################################################################################    
PIPELINE="References"
PIPELINK="References"

LINKS=$LINKS" $PIPELINK"
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
<link rel="stylesheet" type="text/css" href="includes/css/jquery.dataTables.min.css">
<link rel="stylesheet" type="text/css" href="includes/css/ngsane.css">
</head><body>
<div id="center">
<div class='panel' id='quicklinks'><h2>Quicklinks</h2><div>
''' >> $SUMMARYFILE.tmp
declare -a LINKSET=( )
for i in $LINKS; do
    LINKSET=("${LINKSET[@]}" "<a href='#"$i"_panel'>$i</a>")
done
echo $(IFS='|' ; echo "${LINKSET[*]}") >> $SUMMARYFILE.tmp

echo '''
<div id ="Right" style="float:right;width:285px;padding:2px 2px 2px 2px;"> Search: <input type="text" id="search"></div>
</div><!-- Links -->
</div><!-- panel -->''' >>$SUMMARYFILE.tmp

echo "<hr><span>Report generated with "`$NGSANE_BASE/bin/trigger.sh -v`"</span><span style='float:right;'>Last modified: "`date`"</span>" >> $SUMMARYTMP
echo "</div><!-- center --></body>" >> $SUMMARYTMP


echo "<script type='text/javascript'>" >> $SUMMARYTMP
echo '$(document).ready(function() {' >> $SUMMARYTMP

for i in $LINKS; do
    if [ "$i" != "References" ]; then #[ "$i" != "fastQC" ] && [ "$i" != "memechip" ] && 
	    echo "var $i""_Table = $""("\'"#$i""_id_Table"\'").dataTable(       
		$i""_Table_json
	     );" >> $SUMMARYTMP
     fi
done
echo "$""("\'"#search"\'").keyup(function(){" >> $SUMMARYTMP
for i in $LINKS; do
	if [ "$i" != "References" ]; then #[ "$i" != "fastQC" ] && [ "$i" != "memechip" ] && 
		echo "    $i""_Table.fnDraw();" >> $SUMMARYTMP
	fi
done
echo "}); 
} ); </script>" >> $SUMMARYTMP

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
