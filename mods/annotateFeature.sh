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

FILES=${FILES//,/ }

if [[ ! "$ENDING" =~ ".bed" ]] && [[ ! "$ENDING" =~ ".bam" ]]; then
    echo "[ERROR] input file format ($ENDING) not recognized"
    exit 1
fi

if [ -z "$FASTA" ]; then
    echo "[ERROR] no reference provided (FASTA)"
    exit 1
fi
GENOMESIZE=${FASTA%.*}.chrom.sizes

# check library variables are set
if [[ -z "$FEATUREFILE" ]]; then
    echo "[ERROR] library info not set (FEATUREFILE)"
    exit 1;
else
    echo "[NOTE] using libraries in $FEATUREFILE"
fi

if [[ -z "$EXPERIMENTNAME" ]]; then
    echo "[ERROR] please specify EXPERIMENTNAME"
    exit 1
fi

FEATURENAME=${FEATUREFILE##*/}
#is ziped ?
CAT="cat"
if [[ ${FEATUREFILE##*.} == "gz" ]]; then 
    CAT="zcat"; 
    export FEATURENAME=${FEATURENAME/.bed.gz/}
elif [[ ${FEATUREFILE##*.} == "bz2" ]]; 
    then CAT="bzcat"; 
    export FEATURENAME=${FEATURENAME/.bed.bz2/}
else
    export FEATURENAME=${FEATURENAME/.bed/}
fi
[ -d $OUTDIR/$FEATURENAME ] && rm -r $OUTDIR/$FEATURENAME 
mkdir -p $OUTDIR/$FEATURENAME

# get basename of f
RESULTFILES=""
for i in $FILES; do
    n=$(basename $i)
    RESULTFILE="$OUTDIR/$EXPERIMENTNAME${n/%$ENDING/}${STRAND/-/_}.txt"
    [ -f $RESULTFILE ] && rm $RESULTFILE
    RESULTFILES="$RESULTFILES $RESULTFILE" 
done

# delete old bam files unless attempting to recover
if [ -z "$NGSANE_RECOVERFROM" ]; then
    [ -e $(echo $RESULTFILES | cut -d " " -f 1) ] && rm -f ${RESULTFILES//.txt/*}
    [ -e $OUTDIR/$EXPERIMENTNAME"_joined"${STRAND/-/_}".png" ] && rm -f $OUTDIR/$EXPERIMENTNAME"_joined"${STRAND/-/_}*
fi

if [[ -z "$FEATURE_START" && -z "$FEATURE_END" ]]; then
    echo "[ERROR] Either the feature start and/or the end column need to be specified"
    exit 1
fi

if [[ -n "$FEATURE_START" && $(head -n 1 $FEATUREFILE | awk '{FS = "\t"};{print NF}') -lt $FEATURE_START ]]; then
    echo "[ERROR] Feature bed file contains too few columns"
    exit 1
fi
if [[ -n "$FEATURE_END" && $(head -n 1 $FEATUREFILE | awk '{FS = "\t"};{print NF}') -lt $FEATURE_END ]]; then
    echo "[ERROR] Feature bed file contains too few columns"
    exit 1
fi
    
echo "[NOTE] Padding featues with $UPSTREAM and $DOWNSTREAM bps up- and downstream, respectively"
if [ -n "$BIN" ]; then
    echo "[NOTE] bin with $BIN"
fi

CENTERBINS=100

THISTMP=$TMP"/"$(whoami)"/"$(echo $OUTDIR/$FEATURENAME | md5sum | cut -d' ' -f1)
[ -d $THISTMP ] && rm -r $THISTMP
mkdir -p $THISTMP

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "recall files from tape"
	
if [ -n "$DMGET" ]; then
    dmget -a $FEATUREFILE
    dmget -a ${FILES}
fi
    
NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "split features by category"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    if [ -n "$FEATANN_SEGREGATEBY" ]; then
        $CAT $FEATUREFILE | awk -v C=$FEATANN_SEGREGATEBY -v O=$OUTDIR/$FEATURENAME/ '{print > O"/"$C".bed.tmp"}' 

    else
        if [ "$CAT" != "cat" ]; then
            $CAT $FEATUREFILE > $OUTDIR/$FEATURENAME/$FEATURENAME.bed.tmp
        else
            cp $FEATUREFILE $OUTDIR/$FEATURENAME/$FEATURENAME.bed.tmp
        fi
    fi
    NGSANE_CHECKPOINT_CHECK
fi

################################################################################
NGSANE_CHECKPOINT_INIT "pad and clean features"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    for FEATURES in $(ls $OUTDIR/$FEATURENAME/*.tmp 2> /dev/null); do
        F=${FEATURES##*/}
        FEATURE_REGIONS=$OUTDIR/$FEATURENAME/${F/.bed.tmp/}${STRAND/-/_}.bed
#        head  -n 2 $FEATURES
        
        if [ -z "$FEATURE_END" ];then    
            awk -v c=$FEATURE_START 'BEGIN{OFS=FS="\t"}{if(NF >= 6){print $1,$c,$c+1,$4,$5,$6}else{print $1,$c,$c+1,$4,0,"."}}' $FEATURES | sort -u -k1,1 -k2,2g | bedtools slop $STRAND -r $DOWNSTREAM -l $UPSTREAM -g $GENOMESIZE -i - | bedtools sort > $FEATURE_REGIONS
        elif [ -z "$FEATURE_START" ];then    
            awk -v c=$FEATURE_END 'BEGIN{OFS=FS="\t"}{if(NF >= 6){print $1,$c-1,$c,$4,$5,$6}else{print $1,$c-1,$c,$4,0,"."}}' $FEATURES | sort -u -k1,1 -k2,2g | bedtools slop $STRAND -r $DOWNSTREAM -l $UPSTREAM -g $GENOMESIZE -i - | bedtools sort > $FEATURE_REGIONS
        else
            # split into up center and downstream
            awk -v s=$FEATURE_START -v e=$FEATURE_END 'BEGIN{OFS=FS="\t"}{if(NF >= 6){print $1,$s,$e,$4,$5,$6}else{print $1,$s,$e,$4,0,"."}}' $FEATURES | sort -u -k1,1 -k2,2g | bedtools flank $STRAND -r 0 -l $UPSTREAM -g $GENOMESIZE -i - | bedtools sort > ${FEATURE_REGIONS}_upstream
            awk -v s=$FEATURE_START -v e=$FEATURE_END 'BEGIN{OFS=FS="\t"}{if(NF >= 6){print $1,$s,$e,$4,$5,$6}else{print $1,$s,$e,$4,0,"."}}' $FEATURES | sort -u -k1,1 -k2,2g | bedtools flank $STRAND -r $DOWNSTREAM -l 0 -g $GENOMESIZE -i - | bedtools sort > ${FEATURE_REGIONS}_downstream
            cp $FEATURES ${FEATURE_REGIONS}

        fi 
#        head -n 2 $FEATURE_REGIONS       
    
        if [ -n "$REMOVEMULTIFEATURE" ]; then
            bedtools merge $STRAND -n -i $FEATURE_REGIONS | awk '{if ($4==1){$5==""?$5=".":$5=$5;OFS="\t"; print $1,$2,$3,".",".",$5}}' | bedtools sort > $FEATURE_REGIONS.tmp2
            echo "[NOTE] drop "$(echo $(wc -l $FEATURE_REGIONS | cut -d " " -f 1 )-$(wc -l $FEATURE_REGIONS.tmp2 | cut -d " " -f 1) | bc)" multi-feature regions from "$(echo $(wc -l $FEATURE_REGIONS | cut -d " " -f 1 ))
            mv $FEATURE_REGIONS.tmp2 $FEATURE_REGIONS
        fi
#        head -n 2 $FEATURE_REGIONS
    done 
    
    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $FEATURE_REGIONS 
    
    # cleanup 
    rm $OUTDIR/$FEATURENAME/*.tmp 2> /dev/null
fi
################################################################################
NGSANE_CHECKPOINT_INIT "run annotation"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    get_coverage_onesided() {
        . $CONFIG
        . ${NGSANE_BASE}/conf/header.sh
        . $CONFIG
        name=$(basename $1)
        name=$EXPERIMENTNAME${name/%$ENDING/}${STRAND/-/_}
        
        [ -f $OUTDIR/$name.bed ] && rm $OUTDIR/$name.bed
        [ -f $OUTDIR/$name.txt ] && rm $OUTDIR/$name.txt
        
        arrIN=(${1//\// })
        fileloc=$(echo ${arrIN[@]:(-3)} | sed 's/ /\//g')
        fileloc=$(basename $fileloc)
        echo $fileloc
        fgrep -w "$fileloc" $CONFIG
        mark=$(fgrep "$fileloc" $CONFIG | cut -d " " -f 3)
        if [ -z "$mark" ]; then echo "[WARN] skipped undefined mark (provide FEATANN_LAB definition in the config file if you want to consider this data)"; exit; fi
        
        for FEATURES in $(ls $OUTDIR/$FEATURENAME/*.bed 2> /dev/null); do
            F=${FEATURES##*/}
            F=${F/.bed/}
            FEATURE_REGIONS=$OUTDIR/$FEATURENAME/${F}${STRAND/-/_}.bed
#            head  -n 2 $FEATURES    
            echo $FEATURE_REGIONS
            # get number of regions
            FEATURE_LENGTH=$(wc -l $FEATURE_REGIONS | cut -d " " -f 1)
            echo "[NOTE] process $fileloc with mark $mark"
            if [[ "$ENDING" =~ ".bam" ]]; then
                if [ "$NORMALIZE" == "genome" ]; then
                    if [ ! -e $1.stats ]; then
                        echo "[NOTE] flagstat"
                        samtools flagstat $1 >$1.stats
                    fi
                    # normalize by mapped reads
                    TOTALREADS="--normalize "$(head -n 3 $1.stats | tail -n 1 | cut -d " " -f 1)
                    echo "[NOTE] $TOTALREADS"
                elif [ "$NORMALIZE" == "features" ]; then 
                    TOTALREADS="--normalize "$( bedtools coverage -counts -hist -abam $1 -b $FEATURE_REGIONS | awk '{ sum+=$NF} END {print sum}'); 
                    echo "[NOTE] $TOTALREADS"
                fi
                A="-abam"
            elif [[ "$ENDING" =~ ".bed" ]]; then
                if [ "$NORMALIZE" == "genome" ]; then 
                    TOTALREADS="--normalize "$( wc -l $1 | cut -d " " -f 1); 
                elif [ "$NORMALIZE" == "features" ]; then 
                    TOTALREADS="--normalize "$( bedtools coverage -counts -hist -a $1 -b $FEATURE_REGIONS | awk '{ sum+=$NF} END {print sum}'); 
                fi
                A="-a"
            fi
    
    #        if [ -n "$STRAND" ]; then IND=7 ;VAL=8; STR=6; else IND=7 ;VAL=8; STR=10; fi
    #       TODO seems to be the same
            if [ -z "$STRAND" ]; then IND=7 ;VAL=8; STR=6; else IND=7 ;VAL=8; STR=6; fi
    
            echo "[NOTE] coverage $STRAND"
    #        bedtools coverage -d $STRAND $A $1 -b $FEATURE_REGIONS | head -n 2

    
            EXPREG=$(bedtools coverage $STRAND $A $1 -b $FEATURE_REGIONS | gawk -v v=$VAL '{if ($v!=0) {print $0}}' | wc -l)  # non-zero covered features
            if [ $EXPREG != 0 ]; then
                echo "[NOTE] nonzero FEATURE_REGIONS $EXPREG"
                if [ -n "$BIN" ]; then
                    bedtools coverage $STRAND -d $A $1 -b $FEATURE_REGIONS | gawk -v i=$IND -v v=$VAL -v s=$STR '{OFS="\t";print $i,$v,$s}' | gawk -v bin=$BIN  'BEGIN{sum=0;len=1}{if ($1%bin==0){if(len==bin){OFS="\t";print $1,sum/len,$3}; sum=0;len=1}else{if($1<len){sum=0;len=1};sum=sum+$2;len=len+1}}' > $OUTDIR/$name.bed

                else
                    bedtools coverage $STRAND -d $A $1 -b $FEATURE_REGIONS | gawk -v i=$IND -v v=$VAL -v s=$STR '{OFS="\t";print $i,$v,$s}' >> $OUTDIR/$name.bed
                fi
        
    #            head $OUTDIR/$name.bed
    
                echo "[NOTE] process file"
                RUNCOMMAND="python ${NGSANE_BASE}/tools/coverageAtFeature.py -f $OUTDIR/$name.bed -C $F $PYBIN -u $UPSTREAM -d $DOWNSTREAM -l $FEATURE_LENGTH -n $mark -o $OUTDIR/$name $IGNOREUNCOVERED $REMOVEOUTLIER $TOTALREADS --metric $METRIC"
                echo $RUNCOMMAND && eval $RUNCOMMAND
    
            else
                echo "[NOTE] no regions overlap features"
                touch $OUTDIR/$name.txt $OUTDIR/$name.bed;
            fi
            head $OUTDIR/$name.bed
        done
    }
    
    get_coverage_twosided() {
        . $CONFIG
        . ${NGSANE_BASE}/conf/header.sh
        . $CONFIG
        name=$(basename $1)
        name=$EXPERIMENTNAME${name/%$ENDING/}${STRAND/-/_}

        [ -f $OUTDIR/$name.bed ] && rm $OUTDIR/$name.bed
        [ -f $OUTDIR/$name.txt ] && rm $OUTDIR/$name.txt
        
        arrIN=(${1//\// })
        fileloc=$(echo ${arrIN[@]:(-3)} | sed 's/ /\//g')
        fileloc=$(basename $fileloc)
        echo $fileloc
        fgrep -w "$fileloc" $CONFIG
        mark=$(fgrep "$fileloc" $CONFIG | cut -d " " -f 3)
        if [ -z "$mark" ]; then echo "[WARN] skipped undefined mark (provide FEATANN_LAB definition in the config file if you want to consider this data)"; exit; fi
        
        for FEATURES in $(ls $OUTDIR/$FEATURENAME/*.bed 2> /dev/null); do
            F=${FEATURES##*/}
            F=${F/.bed/}
            FEATURE_REGIONS=$OUTDIR/$FEATURENAME/${F}${STRAND/-/_}.bed
#            head  -n 2 $FEATURES    
            echo $FEATURE_REGIONS
            # get number of regions
            FEATURE_LENGTH=$(wc -l $FEATURE_REGIONS | cut -d " " -f 1)
            echo "[NOTE] process $fileloc with mark $mark"
            if [[ "$ENDING" =~ ".bam" ]]; then
                if [ "$NORMALIZE" == "genome" ]; then
                    if [ ! -e $1.stats ]; then
                        echo "[NOTE] flagstat"
                        samtools flagstat $1 >$1.stats
                    fi
                    # normalize by mapped reads
                    TOTALREADS="--normalize "$(head -n 3 $1.stats | tail -n 1 | cut -d " " -f 1)
                    echo "[NOTE] $TOTALREADS"
                elif [ "$NORMALIZE" == "features" ]; then 
                    TOTALREADS="--normalize "$( bedtools coverage -counts -hist -abam $1 -b $FEATURE_REGIONS | awk '{ sum+=$NF} END {print sum}'); 
                    echo "[NOTE] $TOTALREADS"
                fi
                A="-abam"
            elif [[ "$ENDING" =~ ".bed" ]]; then
                if [ "$NORMALIZE" == "genome" ]; then 
                    TOTALREADS="--normalize "$( wc -l $1 | cut -d " " -f 1); 
                elif [ "$NORMALIZE" == "features" ]; then 
                    TOTALREADS="--normalize "$( bedtools coverage -counts -hist -a $1 -b $FEATURE_REGIONS | awk '{ sum+=$NF} END {print sum}'); 
                fi
                A="-a"
            fi
    
    #        if [ -n "$STRAND" ]; then IND=7 ;VAL=8; STR=6; else IND=7 ;VAL=8; STR=10; fi
    #       TODO seems to be the same
            if [ -z "$STRAND" ]; then IND=7 ;VAL=8; STR=6; else IND=7 ;VAL=8; STR=6; fi
    
            echo "[NOTE] coverage $STRAND"
    #        bedtools coverage -d $STRAND $A $1 -b $FEATURE_REGIONS | head -n 2

    
            EXPREGUP=$(bedtools coverage $STRAND $A $1 -b ${FEATURE_REGIONS}_upstream | gawk -v v=$VAL '{if ($v!=0) {print $0}}' | wc -l)  
            EXPREGCENTER=$(bedtools coverage $STRAND $A $1 -b ${FEATURE_REGIONS} | gawk -v v=$VAL '{if ($v!=0) {print $0}}' | wc -l)  
            EXPREGDOWN=$(bedtools coverage $STRAND $A $1 -b ${FEATURE_REGIONS}_downstream | gawk -v v=$VAL '{if ($v!=0) {print $0}}' | wc -l)  
            # non-zero covered features
            if [[ $EXPREGUP != 0 || $EXPREGCENTER != 0 || $EXPREGDOWN != 0 ]]; then
                echo "[NOTE] nonzero FEATURE_REGIONS $EXPREG"
                if [ -n "$BIN" ]; then

                    bedtools coverage $STRAND -d $A $1 -b ${FEATURE_REGIONS}_upstream | \
                        gawk -v i=$IND -v v=$VAL -v s=$STR '{OFS="\t";print $i,$v,$s}' | \
                        gawk -v bin=$BIN  'BEGIN{OFS="\t";sum=0;len=1}{if ($1%bin==0){if(len==bin){print $1,sum/len,$3}; sum=0;len=1}else{if($1<len){sum=0;len=1};sum=sum+$2;len=len+1}}' > $OUTDIR/$name.bed

                    bedtools coverage $STRAND -d $A $1 -b $FEATURE_REGIONS | \
                    gawk -v i=$IND -v v=$VAL -v s=$STR '{OFS="\t";print $3-$2,$i,$v,$s}' | \
                    gawk -v upstream=$UPSTREAM -v bins=$CENTERBINS 'BEGIN{OFS="\t";sum=0;len=1;bin=0; lastbin=0}
                        {bin=$1/bins; 
                        if(lastbin < int($2/bin)){
                            for (i=lastbin;i<int($2/bin);i++){print upstream+i,(sum+$3)/len,$4};
                            lastbin=int($2/bin);
                            len=1;
                            sum=0;
                        } else {
                            len = len+1;
                            sum=sum+$3;
                        }
                        if (lastbin >= bins){
                            lastbin=0;
                        }                             
                    }' >> $OUTDIR/$name.bed >> $OUTDIR/$name.bed

                    bedtools coverage $STRAND -d $A $1 -b ${FEATURE_REGIONS}_downstream | \
                        gawk -v i=$IND -v v=$VAL -v s=$STR '{OFS="\t";print $i,$v,$s}' | 
                        gawk -v bin=$BIN -v upstream=$UPSTREAM -v center=$CENTERBINS 'BEGIN{OFS="\t";sum=0;len=1}{if ($1%bin==0){if(len==bin){print upstream+center+$1,sum/len,$3}; sum=0;len=1}else{if($1<len){sum=0;len=1};sum=sum+$2;len=len+1}}' >> $OUTDIR/$name.bed

                else
                
                    bedtools coverage $STRAND -d $A $1 -b ${FEATURE_REGIONS}_upstream | \
                        gawk -v i=$IND -v v=$VAL -v s=$STR '{OFS="\t";print $i,$v,$s}' > $OUTDIR/$name.bed


                    bedtools coverage $STRAND -d $A $1 -b $FEATURE_REGIONS | \
                       gawk -v i=$IND -v v=$VAL -v s=$STR '{OFS="\t";print $3-$2,$i,$v,$s}' | \
                       gawk -v upstream=$UPSTREAM -v bins=$CENTERBINS 'BEGIN{OFS="\t";sum=0;len=1;bin=0;lastbin=0}
                           {bin=$1/bins; 
                           if(lastbin < int($2/bin)){
                               for (i=lastbin;i<int($2/bin);i++){print upstream+i,(sum+$3)/len,$4};
                               lastbin=int($2/bin);
                               len=1;
                               sum=0;
                           } else {
                               len = len+1;
                               sum=sum+$3;
                           }
                           if (lastbin >= bins){
                               lastbin=0;
                           }                             
                       }' >> $OUTDIR/$name.bed

                    bedtools coverage $STRAND -d $A $1 -b ${FEATURE_REGIONS}_downstream | \
                    gawk -v i=$IND -v v=$VAL -v s=$STR -v upstream=$UPSTREAM -v center=$CENTERBINS '{OFS="\t";print upstream+center+$i,$v,$s}' >> $OUTDIR/$name.bed

                fi
        
    #            head $OUTDIR/$name.bed
    
                echo "[NOTE] process file"
                RUNCOMMAND="python ${NGSANE_BASE}/tools/coverageAtFeature.py -f $OUTDIR/$name.bed -C $F $PYBIN -u $UPSTREAM -d $DOWNSTREAM -c $CENTERBINS -l $FEATURE_LENGTH -n $mark -o $OUTDIR/$name $IGNOREUNCOVERED $REMOVEOUTLIER $TOTALREADS --metric $METRIC"
                echo $RUNCOMMAND && eval $RUNCOMMAND
    
            else
                echo "[NOTE] no regions overlap features"
                touch $OUTDIR/$name.txt $OUTDIR/$name.bed;
            fi
            head $OUTDIR/$name.bed
        done
    }
    export -f get_coverage_onesided
    export -f get_coverage_twosided
    
    if [[ -z "$FEATURE_END" || -z "$FEATURE_START" ]];then    
        echo "[NOTE] one-sided coverage"
    #    parallel --gnu -env get_coverage_onesided ::: $FILES
        for i in $FILES; do get_coverage_onesided $i; done
    else
        echo "[NOTE] two-sided coverage"
    #    parallel --gnu -env get_coverage_twosided ::: $FILES
         for i in $FILES; do get_coverage_twosided $i; done

    fi
    echo "[NOTE] Files $FILES"

	# mark checkpoint
    NGSANE_CHECKPOINT_CHECK $RESULTFILES
    
#    rm ${RESULTFILES//.txt/.bed}
fi
################################################################################
NGSANE_CHECKPOINT_INIT "summarize"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    JOINED="$OUTDIR/${EXPERIMENTNAME}_joined${STRAND/-/_}"
    echo "[NOTE] plot $JOINED.pdf"
    head -n 1 $(echo $RESULTFILES | cut -d " " -f 1) > $JOINED.txt
    cat $RESULTFILES | grep -v "x" >> $JOINED.txt
    
    if [[ -z "$FEATURELABEL" ]]; then
        if [[ -z "$FEATURE_END" ]];then    
            FEATURELABEL="Feature start"
        elif [[ -z "$FEATURE_START" ]];then  
            FEATURELABEL="Feature end"
        else 
            FEATURELABEL="Feature"
        fi
    fi
        
    if [[ -n "$FEATURE_END"  && -n "$FEATURE_END" ]]; then
        python ${NGSANE_BASE}/tools/coverageAtFeature.py -o $JOINED $PYBIN -u $UPSTREAM -c $CENTERBINS -d $DOWNSTREAM -g "$FEATURELABEL" -i $JOINED --metric $METRIC
    else
        python ${NGSANE_BASE}/tools/coverageAtFeature.py -o $JOINED $PYBIN -u $UPSTREAM -d $DOWNSTREAM -g "$FEATURELABEL" -i $JOINED --metric $METRIC
    fi
        
    Rscript $JOINED.R
    convert $JOINED.pdf $JOINED.png
    
    ls $JOINED.png

	# mark checkpoint
    NGSANE_CHECKPOINT_CHECK $JOINED.pdf 
    
    rm $RESULTFILES

fi
################################################################################
#[ -e $f.merg.anno.bed.dummy ] && rm $f.merg.anno.bed.dummy
echo ">>>>> Coverage at genomic regions - FINISHED"
echo ">>>>> enddate "`date`

