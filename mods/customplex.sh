#!/bin/bash

# Custom de multiplexing
# author: Denis Bauer 
# date: Nov. 2011

# messages to look out for -- relevant for the QC.sh script:
# QCVARIABLES,

echo ">>>>> customplex with fastXtoolkit "
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"

function usage {
echo -e "usage: $(basename $0) -k NGSANE -f FASTA -r REFERENCE -o OUTDIR [OPTIONS]

Script

required:
  -k | --toolkit <path>     config file 
  -f | --fastq <file>       fastq file
  -b | --barcode <file>     barcode
  -o | --outdir <path>      output dir
"
exit
}

if [ ! $# -gt 3 ]; then usage ; fi

#DEFAULTS
THREADS=1

#INPUTS
while [ "$1" != "" ]; do
    case $1 in
        -k | toolkit )          shift; CONFIG=$1 ;; # ENSURE NO VARIABLE NAMES FROM CONFIG
        -t | --threads )        shift; THREADS=$1 ;; # number of CPUs to use
        -f | --fastq )          shift; f=$1 ;; # fastq file
        -b | --barcode )        shift; BARCODR=$1 ;; # reference genome
        -o | --outdir )         shift; OUTDIR=$1 ;; # output dir
        -p | --prefix )         shift; PREFIX=$1 ;; # prefix for the line
        --recover-from )        shift; RECOVERFROM=$1 ;; # attempt to recover from log file                                                  
        -h | --help )           usage ;;
        * )                     usage
    esac
    shift
done


#PROGRAMS (note, both configs are necessary to overwrite the default, here:e.g.  TASKTOPHAT)
. $CONFIG
. ${NGSANE_BASE}/conf/header.sh
. $CONFIG

################################################################################
CHECKPOINT="programs"

for MODULE in $MODULE_DEMULTIPLEX; do module load $MODULE; done  # save way to load modules that itself load other modules
export PATH=$PATH_DEMULTIPLEX:$PATH
module list
echo "PATH=$PATH"
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)

echo -e "--NGSANE                     --\n" $(trigger.sh -v 2>&1)
echo -e "--perl                       --\n "$(perl -v  | head -n 2 | tail -n 1)
[ -z "$(which perl)" ] && echo "[ERROR] no perl detected" && exit 1
echo -e "--fastx_barcode_splitter.pl  --\n "$(fastx_barcode_splitter.pl | head -n 1)
[ -z "$(which fastx_barcode_splitter.pl)" ] && echo "[ERROR] no fastx_barcode_splitter.pl detected" && exit 1
echo -e "--fastx_trimmer              --\n "$(fastx_trimmer -h | head -n 2 | tail -n 1)
[ -z "$(which fastx_trimmer)" ] && echo "[ERROR] no fastx_trimmer detected" && exit 1


echo -e "\n********* $CHECKPOINT"
################################################################################
CHECKPOINT="parameters"

#echo $OUTDIR

if [[ ${f##*.} == "gz" ]]; then
    echo "[NOTE] unzip first"
    f=${f/.gz/}
    if [ ! -e $f ]; then gunzip -c $f.gz >$f; fi
    if [ ! -e ${f/$READONE/$READTWO} ]; then gunzip -c ${f/$READONE/$READTWO}.gz >${f/$READONE/$READTWO} ; fi
fi

# get basename of f
n=${f##*/}

################################################################################
CHECKPOINT="read1-read2"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep "********* $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

    # put read1 and read2 side by side
    perl ${NGSANE_BASE}/bin/shuffleSequences_fastq_sidebyside.pl $f ${f/$READONE/$READTWO} \
        $OUTDIR/${n/$READONE/"sidebyside$READONE"}

    # mark checkpoint
    [ -e $OUTDIR/${n/$READONE/"sidebyside$READONE"} ] && echo -e "\n********* $CHECKPOINT" && unset RECOVERFROM
fi 

################################################################################
CHECKPOINT="read1-read2 demult"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep "********* $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

    cat $OUTDIR/${n/$READONE/"sidebyside$READONE"} | fastx_barcode_splitter.pl  \
        --bcfile $CUSTOMBARCODE --bol --mismatches 1 --prefix $OUTDIR/$PREFIX$READONE"_" \
        --suffix "_seq.fastq" > $OUTDIR/$PREFIX$READONE"_read_counts"
    
    [ -e $OUTDIR/$PREFIX$READONE"_read_counts" ] && echo -e "\n********* $CHECKPOINT" && unset RECOVERFROM
fi 

################################################################################
CHECKPOINT="read1-read2 unmatched"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep "********* $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

    # find unmatched
    perl ${NGSANE_BASE}/bin/splitintoforandrevreads.pl $OUTDIR/$PREFIX$READONE"_unmatched"

    [ -e $OUTDIR/$PREFIX$READONE"_unmatched" ] && echo -e "\n********* $CHECKPOINT" && unset RECOVERFROM
fi 

################################################################################
CHECKPOINT="read2-read1"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep "********* $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

    # put for the unmached read one first
    perl ${NGSANE_BASE}/bin/shuffleSequences_fastq_sidebyside.pl \
        $OUTDIR/$PREFIX$READONE"_unmatched_2_seq.fastq" \
        $OUTDIR/$PREFIX$READONE"_unmatched_1_seq.fastq" \
        $OUTDIR/${n/$READONE/"sidebyside$READTWO"}

    [ -e $OUTDIR/${n/$READONE/"sidebyside$READTWO"} ] && echo -e "\n********* $CHECKPOINT" && unset RECOVERFROM
fi 

################################################################################
CHECKPOINT="read2-read1 demult"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep "********* $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 
    
    # demultiplex with read2/read1
    cat $OUTDIR/${n/$READONE/"sidebyside$READTWO"} | fastx_barcode_splitter.pl \
        --bcfile $CUSTOMBARCODE --bol --mismatches 1 --prefix $OUTDIR/$PREFIX$READTWO"_" \
        --suffix "_seq.fastq" > $OUTDIR/$PREFIX$READTWO"_read_counts"

    [ -e $OUTDIR/$PREFIX$READTWO"_read_counts" ] && echo -e "\n********* $CHECKPOINT" && unset RECOVERFROM
fi 

################################################################################
CHECKPOINT="cat and trim"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep "********* $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

    for p in $( ls $OUTDIR/$PREFIX"$READONE"*.fastq ); do
    
        # skip the unmatched ones
        if ( $p == *unmatched* ); then continue; fi
    
        RONE=${p/_seq.fastq/}
        RTWO=${RONE/$READONE/$READTWO}
    
        echo $RONE
        echo $RTWO
    
        # unconcatinate again
        perl ${NGSANE_BASE}/bin/splitintoforandrevreads.pl $RONE
        perl ${NGSANE_BASE}/bin/splitintoforandrevreads_2readfirst.pl $RTWO
    
        # get read1 and read2
        cat $RONE"_1_seq.fastq" $RTWO"_1_seq.fastq" > ${RONE/"$READONE"/}"_read1untr_seq.fastq"
        cat $RONE"_2_seq.fastq" $RTWO"_2_seq.fastq" > ${RTWO/"$READTWO"/}"_read2untr_seq.fastq"
    
        # trim 
        #echo "trimming with sanger quality score (-Q 33) "${RONE/"$READONE"/}"_read1untr_seq.fastq"
        fastx_trimmer -Q 33 -f 7 -z -i ${RONE/"$READONE"/}"_read1untr_seq.fastq" -o ${RONE/"$READONE"/}"_"$READONE.fastq.gz
        #echo "trimming with sanger quality score (-Q 33) "${RTWO/"$READTWO"/}"_read2untr_seq.fastq"
        fastx_trimmer -Q 33 -f 7 -z -i ${RTWO/"$READTWO"/}"_read2untr_seq.fastq" -o ${RTWO/"$READTWO"/}"_"$READTWO.fastq.gz
    
        rm $RONE"_seq.fastq" $RTWO"_seq.fastq"
        rm $RONE"_1_seq.fastq" $RTWO"_1_seq.fastq"
        rm $RONE"_2_seq.fastq" $RTWO"_2_seq.fastq"
        rm ${RONE/$READONE/}"_read1untr_seq.fastq"
        rm ${RTWO/$READTWO/}"_read2untr_seq.fastq"
     
    done
    
    [ -e $OUTDIR/$PREFIX$READTWO"_read_counts" ] && echo -e "\n********* $CHECKPOINT" && unset RECOVERFROM
fi 

################################################################################
CHECKPOINT="cleanup"

[ -e $OUTDIR/$PREFIX$READONE"_unmatched_seq.fastq" ] && rm $OUTDIR/$PREFIX$READONE"_unmatched_seq.fastq"
[ -e $OUTDIR/$PREFIX$READTWO"_unmatched_seq.fastq" ] && rm $OUTDIR/$PREFIX$READTWO"_unmatched_seq.fastq"
[ -e $OUTDIR/${n/$READONE/"sidebyside$READONE"} ] && rm $OUTDIR/${n/$READONE/"sidebyside$READONE"}
[ -e $OUTDIR/${n/$READONE/"sidebyside$READTWO"} ] && rm $OUTDIR/${n/$READONE/"sidebyside$READTWO"}

echo -e "\n********* $CHECKPOINT"
################################################################################
echo ">>>>> customplex with fastXtoolkit - FINISHED"
echo ">>>>> enddate "`date`

exit

#######3
# this is just with one read
######

# need to deal with read one and read two
# split according to barcodes
less $f | fastx_barcode_splitter.pl --bcfile $CUSTOMBARCODE --bol --mismatches 2 \
    --prefix $OUTDIR/$PREFIX --suffix ".fastq"

# trim barcode sequence and gzip
for f in $( ls $OUTDIR/*fastq ); do
    #less $f | fastx_trimmer -f 7 -z -o ${f/.fastq/.r1} -Q 33
    fastx_trimmer -f 7 -z -i $f -o ${f/.fastq/.fastq.gz} -Q 33
done


