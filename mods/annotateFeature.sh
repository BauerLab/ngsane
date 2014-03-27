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

options:
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
        --recover-from )        shift; RECOVERFROM=$1 ;; # attempt to recover from log file
        -h | --help )           usage ;;
        * )                     echo "don't understand "$1
    esac
    shift
done

#PROGRAMS
export CONFIG=$CONFIG
export OUTDIR=$OUTDIR
#export NGSANE_BASE=${NGSANE_BASE}
. $CONFIG
. ${NGSANE_BASE}/conf/header.sh
. $CONFIG

if [ -n "$BIN" ]; then export PYBIN="--bin $BIN"; fi
if [ -n "$STRANDETNESS" ]; then export STRAND="-s"; fi

################################################################################
CHECKPOINT="programs"

for MODULE in $MODULE_FEATANN; do module load $MODULE; done  # save way to load modules that itself load other modules
export PATH=$PATH_FEATANN:$PATH
module list
echo "PATH=$PATH"

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--bedtools    --\n "$(bedtools -version)
[ -z "$(which bedtools)" ] && echo "[WARN] bedtools not detected" && exit 1

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="parameters"

FILES=${FILES//,/ /}

# get basename of f
RESULTFILES=""
for i in $FILES; do
    n=$(basename $i)
    ending=${n/*./}
    RESULTFILES=$RESULTFILES" "$OUTDIR/${n/.$ending/}"-"${UPSTREAM}"+"${DOWNSTREAM}"_"$METRIC${STRAND/-/_}".txt"
done


# check library variables are set
if [[ -z "$FEATUREFILE" ]]; then
    echo "[ERROR] library info not set (FEATUREFILE)"
    exit 1;
else
    echo "[NOTE] using libraries in $FEATUREFILE"
fi

# delete old bam files unless attempting to recover
if [ -z "$RECOVERFROM" ]; then
    [ -e $(echo $RESULTFILES | cut -d " " -f 1) ] && rm ${RESULTFILES//.txt/*}
    [ -e $OUTDIR/joined.png ] && rm $OUTDIR/joined*
fi

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="recall files from tape"
	
if [ -n "$DMGET" ]; then
    dmget -a $FEATUREFILE
	dmget -a ${FILES}
fi
    
echo -e "\n********* $CHECKPOINT\n"

################################################################################
CHECKPOINT="bedmerge"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

    ending=${FASTA/*./}
    GENOMESIZE=${FASTA/.$ending/}.chrom.sizes
    name=$(basename $FEATUREFILE)
    export REGIONS=$OUTDIR/${name/bed/}feat-$UPSTREAM"+"$DOWNSTREAM"_"$METRIC${STRAND/-/_}.bed
	gawk '{print $1"\t"$7"\t"$7+1"\t"$4"\t"$5"\t"$6}' $FEATUREFILE | sort -u -k1,1 -k2,2g | bedtools slop $STRAND -r $DOWNSTREAM -l $UPSTREAM -g $GENOMESIZE -i - | bedtools sort > $REGIONS

    if [ -n "$REMOVEMULTIFEATURE" ]; then
        bedtools merge $STRAND -n -i $REGIONS | gawk '{if ($4==1){print $1"\t"$2"\t"$3"\t"$5}}' > $REGIONS.tmp
        echo "drop "$(echo $(wc -l $REGIONS | cut -d " " -f 1 )-$(wc -l $REGIONS.tmp | cut -d " " -f 1) | bc)" multi-feature regions"
        mv $REGIONS.tmp $REGIONS
    fi

    export LENGTH=$(wc -l $REGIONS | cut -d " " -f 1)
    echo "[NOTE] $LENGTH features in $REGIONS $STRAND"

	# mark checkpoint
    if [ -f $REGIONS ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi

fi
 

################################################################################
CHECKPOINT="run annotation"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

    get_coverage() {
        . $CONFIG
        . ${NGSANE_BASE}/conf/header.sh
        . $CONFIG
        name=$(basename $1)
        ending=${name/*./}
        name=${name/.$ending/}"-"$UPSTREAM"+"$DOWNSTREAM"_"$METRIC${STRAND/-/_}
        arrIN=(${1//\// })
        fileloc=$(echo ${arrIN[@]:(-3)} | sed 's/ /\//g')
        mark=$(grep $fileloc $CONFIG | cut -d " " -f 3)
        if [ -z "$mark" ]; then echo "error no mark definition (please provide FEATANN_LAB in the config file)"; exit; fi
        echo "[NOTE] process $fileloc with mark $mark"
        if [[ "$ending" == "bam" ]]; then
            if [ -n "$NORMALIZE" ]; then
                if [ ! -e $1.stats ]; then
                    echo "[NOTE] flagstat"
                    samtools flagstat $1 >$1.stats
                fi
                TOTALREADS="--normalize "$(head -n 1 $1.stats | cut -d " " -f 1)
            fi
            if [ -n "$STRAND" ]; then IND=5 ;VAL=6; STR=4; else IND=4 ;VAL=5; STR=10; fi
            A="-abam"
        elif [[ "$ending" == "bed" ]]; then
            if [ -n "$NORMALIZE" ]; then TOTALREADS="--normalize "$( wc -l $1 | cut -d " " -f 1); fi
            if [ -n "$STRAND" ]; then IND=5 ;VAL=6; STR=4; else IND=5 ;VAL=5; STR=10; fi
            A="-a"
        else
            echo "input file format not recognized"
            exit
        fi
        echo "[NOTE] coverage $STRAND"
#	    bedtools coverage $STRAND -d $A $1 -b $REGIONS | head
        EXPREG=$(bedtools coverage $STRAND $A $1 -b $REGIONS | gawk '{if ($5!=0) {print $0}}' | wc -l)  # non-zero covered features
        echo "nonzero regions $EXPREG"
        if [ $EXPREG == 0 ]; then touch $OUTDIR/$name.txt;  continue; fi
        if [ -n "$BIN" ]; then
            echo "bin with $BIN"
            bedtools coverage -d $A $1 -b $REGIONS | gawk -v i=$IND -v v=$VAL -v s=$STR '{print $i"\t"$v"\t"$s}' | gawk -v bin=$BIN 'BEGIN{sum=0;len=1}{if ($1%bin==0){if(len==bin){print $1"\t"sum/len}; sum=0;len=1}else{if($1<len){sum=0;len=1};sum=sum+$2;len=len+1}}' > $OUTDIR/$name.bed
        else
            bedtools coverage -d $A $1 -b $REGIONS | gawk -v i=$IND -v v=$VAL -v s=$STR '{print $i"\t"$v"\t"$s}' > $OUTDIR/$name.bed
        fi

        echo "[NOTE] process file"
        python ${NGSANE_BASE}/tools/coverageAtFeature.py -f $OUTDIR/$name.bed $PYBIN -u $UPSTREAM -d $DOWNSTREAM -l $LENGTH -n $mark -o $OUTDIR/$name $IGNOREUNCOVERED $REMOVEOUTLIER $TOTALREADS --metric $METRIC
    }
    export -f get_coverage

    echo "[NOTE] Files $FILES"
    parallel --gnu -env get_coverage ::: $FILES
#    for i in $FILES; do get_coverage $i; done

	# mark checkpoint
    if [ -f $(echo $RESULTFILES | cut -d " " -f 1) ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
    
    rm ${RESULTFILES//.txt/.bed}
fi


################################################################################
CHECKPOINT="summarize"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

    JOINED=$OUTDIR/"joined-"$UPSTREAM"+"$DOWNSTREAM"_"$METRIC${STRAND/-/_}
    echo "[NOTE] plot $JOINED.pdf"
    head -n 1 $(echo $RESULTFILES | cut -d " " -f 1) > $JOINED.txt
    cat $RESULTFILES | grep -v "x" >> $JOINED.txt
    python ${NGSANE_BASE}/tools/coverageAtFeature.py -o $JOINED $PYBIN -u $UPSTREAM -d $DOWNSTREAM -l $LENGTH -g "TSS" -i $JOINED --metric $METRIC
    Rscript $JOINED.R
    convert $JOINED.pdf $JOINED.png
    
	# mark checkpoint
    if [ -f $JOINED.pdf ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
    
    rm $RESULTFILES

fi

################################################################################
#[ -e $f.merg.anno.bed.dummy ] && rm $f.merg.anno.bed.dummy
echo ">>>>> Coverage at genomic regions - FINISHED"
echo ">>>>> enddate "`date`

