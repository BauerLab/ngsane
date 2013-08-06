#!/bin/bash

echo ">>>>> HiC analysis with homer"
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> hicHomer.sh $*"

function usage {
echo -e "usage: $(basename $0) -k NGSANE -f bam -o OUTDIR [OPTIONS]

Script running HIC HOMER pipeline tapping into bowtie2
It expects bam files, paired end, as input.

required:
  -k | --toolkit <path>     location of the NGSANE repository 
  -f | --bam <file>         bam file
  -o | --outdir <path>      output dir

options:
  -t | --threads <nr>       number of CPUs to use (default: 8)
"
exit
}
# QCVARIABLES,Resource temporarily unavailable

if [ ! $# -gt 3 ]; then usage ; fi

#DEFAULTS
THREADS=8

#INPUTS                                                                                                           
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the NGSANE repository
        -t | --threads )        shift; THREADS=$1 ;; # number of CPUs to use
        -f | --bam )            shift; f=$1 ;; # bam file
        -o | --outdir )         shift; MYOUT=$1 ;; # output dir 
        -h | --help )           usage ;;
        * )                     echo "don't understand "$1
    esac
    shift
done

#PROGRAMS
. $CONFIG
. ${NGSANE_BASE}/conf/header.sh
. $CONFIG

echo "********** programs"
for MODULE in $MODULE_HOMERHIC; do module load $MODULE; done  # save way to load modules that itself load other modules

export PATH=$PATH_HOMERHIC:$PATH
module list
echo "PATH=$PATH"
#this is to get the full path (modules should work but for path we need the full path and this is the\
# best common denominator)
echo -e "--R       --\n "$(R --version | head -n 3)
[ -z "$(which R)" ] && echo "[ERROR] no R detected" && exit 1
echo -e "--homer   --\n "$(which makeTagDirectory)
[ -z "$(which makeTagDirectory)" ] && echo "[ERROR] homer not detected" && exit 1
echo -e "--circos  --\n "$(circos --version)
[ -z "$(which circos)" ] && echo "[WARN] circos not detected"

if [ "$HOMER_HIC_INTERACTIONS" != "all" ] && [ "$HOMER_HIC_INTERACTIONS" != "cis" ] && [ "$HOMER_HIC_INTERACTIONS" != "trans" ]; then
    echo "[ERROR] HiC interactions not specified (all, cis or trans) : $HOMER_HIC_INTERACTIONS"
fi

# get basename of f
n=${f##*/}

#is paired ?                                                                                                      
if [ -e ${f/$READONE/$READTWO} ]; then
    PAIRED="1"
else
    PAIRED="0"
fi

FASTASUFFIX=${FASTA##*.}

if [ -n "$DMGET" ]; then
	echo "********** reacall files from tape"
	dmget -a ${f/$READONE/"*"}
	dmls -l ${f/$READONE/"*"}
fi

if [ $PAIRED == "0" ]; then 
    echo "[ERROR] paired library required for HIC analysis"
    exit 1
fi

#homer likes to write in the current directory, so change to target
CURDIR=$(pwd)
cd $MYOUT

echo "********* makeTagDirectory" 
RUN_COMMAND="makeTagDirectory $MYOUT/${n/'_'$READONE.$ASD.bam/_tagdir_unfiltered} $f,${f/'_'$READONE/'_'$READTWO} $HOMER_HIC_TAGDIR_OPTIONS"
echo $RUN_COMMAND
eval $RUN_COMMAND

cp -r $MYOUT/${n/'_'$READONE.$ASD.bam/_tagdir_unfiltered} $MYOUT/${n/'_'$READONE.$ASD.bam/_tagdir_filtered}

RUN_COMMAND="makeTagDirectory $MYOUT/${n/'_'$READONE.$ASD.bam/_tagdir_filtered} -update $HOMER_HIC_TAGDIR_OPTIONS"
echo $RUN_COMMAND
eval $RUN_COMMAND

echo "********* create background model" 

RUN_COMMAND="analyzeHiC $MYOUT/${n/'_'$READONE.$ASD.bam/_tagdir_filtered} $HOMER_HIC_BACKGROUND_OPTIONS -createModel $MYOUT/${n/'_'$READONE.$ASD.bam/_background.txt} active.model.txt -cpu $THREADS"
echo $RUN_COMMAND
eval $RUN_COMMAND

echo "********* normalize matrices"

if [ "$HOMER_HIC_INTERACTIONS" == "all" ]; then
    RUN_COMMAND="analyzeHiC $MYOUT/${n/'_'$READONE.$ASD.bam/_tagdir_filtered} $HOMER_HIC_NORMALIZE_OPTIONS -model $MYOUT/${n/'_'$READONE.$ASD.bam/_background.txt}  > $MYOUT/${n/'_'$READONE.$ASD.bam/_matrix.txt}"
    echo $RUN_COMMAND
    eval $RUN_COMMAND

elif [ "$HOMER_HIC_INTERACTIONS" == "cis" ]; then
    [ ! -f $FASTA.fai ] && samtools faidx $FASTA

    for CHR in $(awk '{print $1'} $FASTA.fai); do
	    RUN_COMMAND="analyzeHiC $MYOUT/${n/'_'$READONE.$ASD.bam/_tagdir_filtered} $HOMER_HIC_NORMALIZE_OPTIONS -chr $CHR -model $MYOUT/${n/'_'$READONE.$ASD.bam/_background.txt}  > $MYOUT/${n/'_'$READONE.$ASD.bam/_${CHR}_matrix.txt}"
	    echo $RUN_COMMAND
	    eval $RUN_COMMAND
    done
elif [ "$HOMER_HIC_INTERACTIONS" == "trans" ]; then
   [ ! -f $FASTA.fai ] && samtools faidx $FASTA

    for CHR1 in $(awk '{print $1'} $FASTA.fai); do
        for CHR2 in $(awk '{print $1'} $FASTA.fai); do
            if [ "$CHR1" != "$CHR2" ]; then
                RUN_COMMAND="analyzeHiC $MYOUT/${n/'_'$READONE.$ASD.bam/_tagdir_filtered} $HOMER_HIC_NORMALIZE_OPTIONS -chr $CHR1 -chr2 $CHR2 -model $MYOUT/${n/'_'$READONE.$ASD.bam/_background.txt}  > $MYOUT/${n/'_'$READONE.$ASD.bam/_${CHR1}-${CHR2}_matrix.txt}"
                echo $RUN_COMMAND
                eval $RUN_COMMAND
            fi
        done
    done
fi



echo $RUN_COMMAND
eval $RUN_COMMAND

echo "********* PCA clustering"

RUN_COMMAND="runHiCpca.pl ${n/'_'$READONE.$ASD.bam/} $MYOUT/${n/'_'$READONE.$ASD.bam/_tagdir_filtered} $HOMER_HIC_PCA_OPTIONS -cpu $THREADS "
echo $RUN_COMMAND
eval $RUN_COMMAND

echo "********* Significant interactions"

if [ "$HOMER_HIC_INTERACTIONS" == "all" ]; then
    RUN_COMMAND="analyzeHiC $MYOUT/${n/'_'$READONE.$ASD.bam/_tagdir_filtered} $HOMER_HIC_INTERACTION_OPTIONS -interactions $MYOUT/${n/'_'$READONE.$ASD.bam/_significantInteractions.txt} -nomatrix -cpu $THREADS "
    echo $RUN_COMMAND
    eval $RUN_COMMAND

elif [ "$HOMER_HIC_INTERACTIONS" == "cis" ]; then
    for CHR in $(awk '{print $1'} $FASTA.fai); do
        RUN_COMMAND="analyzeHiC $MYOUT/${n/'_'$READONE.$ASD.bam/_tagdir_filtered} $HOMER_HIC_INTERACTION_OPTIONS -chr $CHR -interactions $MYOUT/${n/'_'$READONE.$ASD.bam/_significantInteractions.txt} -nomatrix -cpu $THREADS "
        echo $RUN_COMMAND
        eval $RUN_COMMAND
    done

elif [ "$HOMER_HIC_INTERACTIONS" == "trans" ]; then
    for CHR1 in $(awk '{print $1'} $FASTA.fai); do
        for CHR2 in $(awk '{print $1'} $FASTA.fai); do
            if [ "$CHR1" != "$CHR2" ]; then
                RUN_COMMAND="analyzeHiC $MYOUT/${n/'_'$READONE.$ASD.bam/_tagdir_filtered} $HOMER_HIC_INTERACTION_OPTIONS -chr $CHR1 -chr2 $CHR2 -interactions $MYOUT/${n/'_'$READONE.$ASD.bam/_significantInteractions.txt} -nomatrix -cpu $THREADS "
                echo $RUN_COMMAND
                eval $RUN_COMMAND
            fi
        done
    done
fi

echo "********* Annotate interactions"

RUN_COMMAND="annotateInteractions.pl $MYOUT/${n/'_'$READONE.$ASD.bam/_significantInteractions.txt} $HOMER_HIC_ANNOTATE_OPTIONS $MYOUT/${n/'_'$READONE.$ASD.bam/_annotations}"
echo $RUN_COMMAND
eval $RUN_COMMAND

if hash ${CIRCOS} 2>&- ; then
    echo "********* Circos plots"
    RUN_COMMAND="analyzeHiC $MYOUT/${n/'_'$READONE.$ASD.bam/_tagdir_filtered} -res 1000000 -pvalue 1e-7 -cpu $THREADS -circos ${n/'_'$READONE.$ASD.bam/} -minDist 2000000000 -nomatrix"
    echo $RUN_COMMAND
    eval $RUN_COMMAND
fi

# and back to where we used to be
cd $CURDIR

echo ">>>>> HiC analysis with homer - FINISHED"
echo ">>>>> enddate "`date`

