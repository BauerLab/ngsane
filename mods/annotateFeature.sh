#!/bin/bash -e

# Coverage of bam/bed files around genomic regions
# author: Denis C. Bauer
# date: Mar.2014

# messages to look out for -- relevant for the QC.sh script:
# QCVARIABLES,Resource temporarily unavailable

echo ">>>>> Coverage at genomic regions"
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"


function usage {
echo -e "usage: $(basename $0) -k CONFIG -f BAM

Coverage at genomic regions 

required:
  -k | --toolkit <path>     config file
  -f <file>                 bam file
  -o | --outdir <path>      location for results

"
exit
}


if [ ! $# -gt 3 ]; then usage ; fi

#INPUTS
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the NGSANE repository
        -f | --file )           shift; FILES=$1 ;; # files
        -o | --outdir )         shift; OUTDIR=$1 ;; # output dir                                                       
        --recover-from )        shift; NGSANE_RECOVERFROM=$1 ;; # attempt to recover from log file
        -h | --help )           usage ;;
        * )                     echo "don't understand "$1
    esac
    shift
done

#PROGRAMS
export CONFIG=$CONFIG
export OUTDIR=$OUTDIR
. $CONFIG
. ${NGSANE_BASE}/conf/header.sh
. $CONFIG

if [ -n "$BIN" ]; then export PYBIN="--bin $BIN"; fi
if [ -n "$STRANDETNESS" ]; then export STRAND="-s"; fi

################################################################################
NGSANE_CHECKPOINT_INIT "programs"

# save way to load modules that itself loads other modules
hash module 2>/dev/null && for MODULE in $MODULE_FEATANN; do module load $MODULE; done && module list 

export PATH=$PATH_FEATANN:$PATH
echo "PATH=$PATH"

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--bedtools    --\n "$(bedtools -version)
[ -z "$(which bedtools)" ] && echo "[WARN] bedtools not detected" && exit 1

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "parameters"

FILES=${FILES//,/ /}

# get basename of f
RESULTFILES=""
for i in $FILES; do
    n=$(basename $i)
    RESULTFILES=$RESULTFILES" "$OUTDIR/${n/%$ENDING/}"-"${UPSTREAM}"+"${DOWNSTREAM}"_"$METRIC${STRAND/-/_}".txt"
done

# check library variables are set
if [[ -z "$FEATUREFILE" ]]; then
    echo "[ERROR] library info not set (FEATUREFILE)"
    exit 1;
else
    echo "[NOTE] using libraries in $FEATUREFILE"
fi

# delete old bam files unless attempting to recover
if [ -z "$NGSANE_RECOVERFROM" ]; then
    [ -e $(echo $RESULTFILES | cut -d " " -f 1) ] && rm -f ${RESULTFILES//.txt/*}
    [ -e $OUTDIR/joined"-"${UPSTREAM}"+"${DOWNSTREAM}"_"$METRIC${STRAND/-/_}.png ] && rm -f $OUTDIR/joined"-"${UPSTREAM}"+"${DOWNSTREAM}"_"$METRIC${STRAND/-/_}*
fi

if [[ -z "$FEATURE_START" && -z "$FEATURE_END" ]]; then
    echo "[ERROR] Either the feature start and/or the end column need to be specified"
    exit 1
fi

if [ -n "$FEATURE_START" && $(head -n 1 $FEATUREFILE | awk '{FS = "\t"};{print NF}') -lt $FEATURE_START ]; then
    echo "[ERROR] Feature bed file contains too few columns"
    exit 1
fi
if [[ -n "$FEATURE_END" && $(head -n 1 $FEATUREFILE | awk '{FS = "\t"};{print NF}') -lt $FEATURE_END ]]; then
    echo "[ERROR] Feature bed file contains too few columns"
    exit 1
fi
    
GENOMESIZE=${FASTA%.*}.chrom.sizes
FEATURENAME=${FEATUREFILE##*/}
export REGIONS=$OUTDIR/${FEATURENAME/bed/}feat-$UPSTREAM"+"$DOWNSTREAM"_"$METRIC${STRAND/-/_}.bed

echo "[NOTE] Padding featues with $UPSTREAM and $DOWNSTREAM bps up- and downstream, respectively"

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "recall files from tape"
	
if [ -n "$DMGET" ]; then
    dmget -a $FEATUREFILE
    dmget -a ${FILES}
fi
    
NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "pad and clean features"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    head $FEATUREFILE
    if [ -z "$FEATURE_END" ];then    
        awk -v c=$FEATURE_START 'BEGIN{OFS=FS="\t"}{if(NF >= 6){print $1,$c,$c+1,$4,$5,$6}else{print $1,$c,$c+1,$4,0,"."}}' $FEATUREFILE | sort -u -k1,1 -k2,2g | bedtools slop $STRAND -r $DOWNSTREAM -l $UPSTREAM -g $GENOMESIZE -i - | bedtools sort > $REGIONS
    elif [ -z "$FEATURE_START" ];then    
        awk -v c=$FEATURE_END 'BEGIN{OFS=FS="\t"}{if(NF >= 6){print $1,$c-1,$c,$4,$5,$6}else{print $1,$c-1,$c,$4,0,"."}}' $FEATUREFILE | sort -u -k1,1 -k2,2g | bedtools slop $STRAND -r $DOWNSTREAM -l $UPSTREAM -g $GENOMESIZE -i - | bedtools sort > $REGIONS
    else
        awk -v s=$FEATURE_START -v e=$FEATURE_END 'BEGIN{OFS=FS="\t"}{if(NF >= 6){print $1,$s,$e,$4,$5,$6}else{print $1,$s,$e,$4,0,"."}}' $FEATUREFILE | sort -u -k1,1 -k2,2g | bedtools slop $STRAND -r $DOWNSTREAM -l $UPSTREAM -g $GENOMESIZE -i - | bedtools sort > $REGIONS
    fi 

#    if [[ $(head -n 1 $FEATUREFILE | awk '{FS = "\t"};{print NF}') -ge 7 ]]; then
#    	gawk '{OFS="\t"; print $1,$7,$7+1,$4,$5,$6}' $FEATUREFILE | sort -u -k1,1 -k2,2g | bedtools slop $STRAND -r $DOWNSTREAM -l $UPSTREAM -g $GENOMESIZE -i - | bedtools sort > $REGIONS
#    else
#        cat $FEATUREFILE | sort -u -k1,1 -k2,2g | bedtools slop $STRAND -r $DOWNSTREAM -l $UPSTREAM -g $GENOMESIZE -i - | bedtools sort > $REGIONS
#    fi
    
#    head $REGIONS

    if [ -n "$REMOVEMULTIFEATURE" ]; then
        bedtools merge $STRAND -n -i $REGIONS | awk '{if ($4==1){$5==""?$5=".":$5=$5;OFS="\t"; print $1,$2,$3,".",".",$5}}' | bedtools sort > $REGIONS.tmp
        echo "[NOTE] drop "$(echo $(wc -l $REGIONS | cut -d " " -f 1 )-$(wc -l $REGIONS.tmp | cut -d " " -f 1) | bc)" multi-feature regions"
        mv $REGIONS.tmp $REGIONS
    fi

    head -n 2 $REGIONS
 
    # mark checkpoint
    [ -f $REGIONS ] && NGSANE_CHECKPOINT_CHECK 
fi
################################################################################

# get number of regions
export LENGTH=$(wc -l $REGIONS | cut -d " " -f 1)
if [ $LENGTH == 0 ]; then
    echo "[ERROR] No features left in $REGIONS $STRAND"
    exit 1
else
    echo "[NOTE] $LENGTH features in $REGIONS $STRAND"
fi

################################################################################
NGSANE_CHECKPOINT_INIT "run annotation"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    get_coverage() {
        . $CONFIG
        . ${NGSANE_BASE}/conf/header.sh
        . $CONFIG
        name=$(basename $1)
        name=${name/%$ENDING/}"-"$UPSTREAM"+"$DOWNSTREAM"_"$METRIC${STRAND/-/_}
        arrIN=(${1//\// })
        fileloc=$(echo ${arrIN[@]:(-3)} | sed 's/ /\//g')
        mark=$(grep $fileloc $CONFIG | cut -d " " -f 3)
        if [ -z "$mark" ]; then echo "[WARN] skipped undefined mark (provide FEATANN_LAB definition in the config file if you want to consider this data)"; exit; fi
        
        echo "[NOTE] process $fileloc with mark $mark"
        if [[ "$ENDING" =~ ".bam" ]]; then
            if [ -n "$NORMALIZE" ]; then
                if [ ! -e $1.stats ]; then
                    echo "[NOTE] flagstat"
                    samtools flagstat $1 >$1.stats
                fi
                TOTALREADS="--normalize "$(head -n 1 $1.stats | cut -d " " -f 1)
            fi
            A="-abam"
        elif [[ "$ENDING" =~ ".bed" ]]; then
            if [ -n "$NORMALIZE" ]; then TOTALREADS="--normalize "$( wc -l $1 | cut -d " " -f 1); fi
            A="-a"
        else
            echo "[ERROR] input file format ($ENDING) not recognized"
            exit 1
        fi

#        if [ -n "$STRAND" ]; then IND=7 ;VAL=8; STR=6; else IND=7 ;VAL=8; STR=10; fi
#       TODO seems to be the same
        if [ -z "$STRAND" ]; then IND=7 ;VAL=8; STR=6; else IND=7 ;VAL=8; STR=6; fi

        echo "[NOTE] coverage $STRAND"
#        bedtools coverage -d $STRAND $A $1 -b $REGIONS | head -n 2

        EXPREG=$(bedtools coverage $STRAND $A $1 -b $REGIONS | gawk -v v=$VAL '{if ($v!=0) {print $0}}' | wc -l)  # non-zero covered features
        if [ $EXPREG != 0 ]; then
            echo "[NOTE] nonzero regions $EXPREG"
            if [ -n "$BIN" ]; then
                echo "bin with $BIN"
                bedtools coverage $STRAND -d $A $1 -b $REGIONS | gawk -v i=$IND -v v=$VAL -v s=$STR '{OFS="\t";print $i,$v,$s}' | gawk -v bin=$BIN 'BEGIN{sum=0;len=1}{if ($1%bin==0){if(len==bin){print $1"\t"sum/len"\t"$3}; sum=0;len=1}else{if($1<len){sum=0;len=1};sum=sum+$2;len=len+1}}' > $OUTDIR/$name.bed
            else
                bedtools coverage $STRAND -d $A $1 -b $REGIONS | gawk -v i=$IND -v v=$VAL -v s=$STR '{OFS="\t";print $i,$v,$s}' > $OUTDIR/$name.bed
            fi
    
#            head $OUTDIR/$name.bed

            echo "[NOTE] process file"
            RUNCOMMAND="python ${NGSANE_BASE}/tools/coverageAtFeature.py -f $OUTDIR/$name.bed $PYBIN -u $UPSTREAM -d $DOWNSTREAM -l $LENGTH -n $mark -o $OUTDIR/$name $IGNOREUNCOVERED $REMOVEOUTLIER $TOTALREADS --metric $METRIC"
            echo $RUNCOMMAND && eval $RUNCOMMAND

        else
            echo "[NOTE] no regions overlap features"
            touch $OUTDIR/$name.txt $OUTDIR/$name.bed;
        fi
        head $OUTDIR/$name.bed

    }
    export -f get_coverage

    echo "[NOTE] Files $FILES"
    parallel --gnu -env get_coverage ::: $FILES
    #for i in $FILES; do get_coverage $i; done

	# mark checkpoint
    [ -f $(echo $RESULTFILES | cut -d " " -f 1) ] && NGSANE_CHECKPOINT_CHECK 
    
    rm ${RESULTFILES//.txt/.bed}
fi
################################################################################
NGSANE_CHECKPOINT_INIT "summarize"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    JOINED=$OUTDIR/"joined-"$UPSTREAM"+"$DOWNSTREAM"_"$METRIC${STRAND/-/_}
    echo "[NOTE] plot $JOINED.pdf"
    head -n 1 $(echo $RESULTFILES | cut -d " " -f 1) > $JOINED.txt
    cat $RESULTFILES | grep -v "x" >> $JOINED.txt
    python ${NGSANE_BASE}/tools/coverageAtFeature.py -o $JOINED $PYBIN -u $UPSTREAM -d $DOWNSTREAM -l $LENGTH -g "TSS" -i $JOINED --metric $METRIC
    Rscript $JOINED.R
    convert $JOINED.pdf $JOINED.png
    
    ls $JOINED.png

	# mark checkpoint
    [ -f $JOINED.pdf ] && NGSANE_CHECKPOINT_CHECK 
    
    rm $RESULTFILES

fi
################################################################################
#[ -e $f.merg.anno.bed.dummy ] && rm $f.merg.anno.bed.dummy
echo ">>>>> Coverage at genomic regions - FINISHED"
echo ">>>>> enddate "`date`

