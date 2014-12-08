#!/bin/bash -e

# Script to run Sailfish
# It takes comma-seprated list of files containing short sequence reads in fasta or fastq format and bowtie index files as input.
# It produces quant.sf & quant_biasquant.sf count files.
# author: Boris Guennewig
# date: 2014

# QCVARIABLES,Resource temporarily unavailable
# RESULTFILENAME <DIR>/<TASK>/reads.count_info


echo ">>>>> feature counting with Sailfish"
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"

#INPUTS
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of NGSANE
        -f | --file )           shift; f=$1 ;; # fastq file
        -o | --outdir )         shift; OUTDIR=$1 ;; # output dir
        --recover-from )        shift; RECOVERFROM=$1 ;; # attempt to recover from log file
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
NGSANE_CHECKPOINT_INIT "programs"

for MODULE in $MODULE_SAILFISH; do module load $MODULE; done  # save way to load modules that itself load other modules
export PATH=$PATH_SAILFISH:$PATH
module list
echo "PATH=$PATH"

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--sailfish --\n "$(sailfish --version  | grep version  | tr -d ' ')
[ -z "$(which sailfish)" ] && echo "[ERROR] no sailfish detected" && exit 1

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "parameters"

# get basename of input file f
n=${f##*/}
SAMPLE=${n/%$READONE.$FASTQ/}

if [ -z "$FASTA" ]; then
    echo "[ERROR] no reference provided (FASTA)"
    exit 1
fi

#is paired ?
if [ "$f" != "${f/%$READONE.$FASTQ/$READTWO.$FASTQ}" ] && [ -e ${f/%$READONE.$FASTQ/$READTWO.$FASTQ} ]; then
    echo "[NOTE] PAIRED library"
    PAIRED="1"
else
    echo "[NOTE] SINGLE library"
    PAIRED="0"
fi

## is ziped ?
ZCAT="cat" # always cat
if [[ ${f##*.} == "gz" ]]; then # unless its zipped
    ZCAT="zcat" 
elif [[ ${f##*.} == "bz2" ]]; then 
    ZCAT="bzcat"; 
fi

#remove old files
if [ -z "$RECOVERFROM" ]; then
    if [ -d $OUTDIR ]; then rm -r $OUTDIR; fi
fi

# unique temp folder that should be used to store temporary files
THISTMP=$TMP"/"$(whoami)"/"$(echo $OUTDIR | md5sum | cut -d' ' -f1)
mkdir -p $THISTMP

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "recall files from tape"

if [ -n "$DMGET" ]; then
    dmget -a $(dirname $FASTA)/*
    dmget -a ${f/$READONE/"*"}
    dmget -a $OUTDIR/*
fi

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "Sailfish run"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    if [ "$PAIRED" == 1 ]; then 

        sailfish quant \
        -i $FASTA \
        -o $OUTDIR/$SAMPLE \
        -p $CPU_SAILFISH \
        -l $SAILFISHLIBRARY \
        -1 <($ZCAT $f) \
        -2 <($ZCAT ${f/%$READONE.$FASTQ/$READTWO.$FASTQ})
    
    else
        sailfish quant \
        -i $FASTA \
        -o $OUTDIR/$SAMPLE \
        -p $CPU_SAILFISH \
        -l $SAILFISHLIBRARY \
        -1 <($ZCAT $f)
    fi

    # mark checkpoint
     NGSANE_CHECKPOINT_CHECK $OUTDIR/$SAMPLE/quant_bias_corrected.sf $OUTDIR/$SAMPLE/quant.sf 
fi
################################################################################
NGSANE_CHECKPOINT_INIT "postprocess sailfish reads"
    
if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then
  
    if [ ! -d $OUTDIR/Combined/ ]; then
        mkdir -p ${OUTDIR}/Combined; else echo " Couldn't generate Combined."
    fi
    
    if [ ! -f $OUTDIR/"$SAMPLE"_quant.sf_cut ]; then
        sed '1,5d' $OUTDIR/$SAMPLE/quant.sf > $OUTDIR/"$SAMPLE"_quant.sf_cut; else echo " Couldn't delete first 5 lines of quant_sf."
    fi

    if [ ! -f $OUTDIR/"$SAMPLE"_quant.sf_sort ]; then
        sort $OUTDIR/"$SAMPLE"_quant.sf_cut > $OUTDIR/"$SAMPLE"_quant.sf_sort; else echo " Couldn't sort quant_sf."
    fi

    if [ ! -f $OUTDIR/Combined/"$SAMPLE"_biassf_value ]; then
        (echo "$SAMPLE"; awk '{print $7}' $OUTDIR/"$SAMPLE"_quant.sf_sort) > $OUTDIR/Combined/"$SAMPLE"_quant.sf_value; else echo " Couldn't generate final quant.sf."
    fi


    if [ ! -f $OUTDIR/"$SAMPLE"_quant_bias_corrected.sf_cut ]; then
        sed '1,5d' $OUTDIR/$SAMPLE/quant_bias_corrected.sf > $OUTDIR/"$SAMPLE"_quant_bias_corrected.sf_cut; else echo " Couldn't delete first 5 lines of bias_sf."
    fi

    if [ ! -f $OUTDIR/"$SAMPLE"_biassf_value_sort ]; then
        sort $OUTDIR/"$SAMPLE"_quant_bias_corrected.sf_cut > $OUTDIR/"$SAMPLE"_biassf_value_sort; else echo " Couldn't sort bias_sf."
    fi

    if [ ! -f $OUTDIR/Combined/"$SAMPLE"_biassf_value ]; then
        (echo "$SAMPLE"; awk '{print $7}' $OUTDIR/"$SAMPLE"_biassf_value_sort) > $OUTDIR/Combined/"$SAMPLE"_biassf_value; else echo " Couldn't generate final bias_sf."
    fi

    if [ ! -f $OUTDIR/Combined/sailfish_transcripts ]; then
        (echo transcripts; awk '{print $1}' $OUTDIR/"$SAMPLE"_quant.sf_cut) > $OUTDIR/Combined/sailfish_transcripts; else echo " Couldn't extract Transcripts."
    fi

# mark checkpoint
   NGSANE_CHECKPOINT_CHECK $OUTDIR/Combined/$SAMPLE_biassf_value $OUTDIR/Combined/"$SAMPLE"_quant.sf_value
   
fi
# ################################################################################
NGSANE_CHECKPOINT_INIT "summarize"

cat $OUTDIR/Combined/$SAMPLE_biassf_value | awk '{print "masked",$0}' >> $OUTDIR/$SAMPLE.summary.txt
   
NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "cleanup"    

rm $OUTDIR/"$SAMPLE"_quant.sf_cut 
rm $OUTDIR/"$SAMPLE"_quant.sf_sort 
rm $OUTDIR/"$SAMPLE"_quant_bias_corrected.sf_cut 
rm $OUTDIR/"$SAMPLE"_biassf_value_sort 
rm $OUTDIR/$SAMPLE/reads.sfc

NGSANE_CHECKPOINT_CHECK
################################################################################
[ -e $OUTDIR/Combined/$SAMPLE_biassf_value.dummy ] && rm $OUTDIR/Combined/$SAMPLE_biassf_value.dummy 
echo ">>>>> Counting with Sailfish - FINISHED"
echo ">>>>> enddate "`date`
