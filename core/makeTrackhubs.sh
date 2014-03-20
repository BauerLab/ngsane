#!/bin/bash

# author: Fabian Buske
# date: Mar 2014

echo ">>>>> Generate UCSC trackhubs "
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

if [ -z "UCSC_GENOMEASSEMBLY" ]; then
    echo "[ERROR] genome assembly not specified"
    exit 1
fi

if [ -z "$TRACKHUB_DIR" ]; then
    echo "[ERROR] output folder not specified"
    exit 1
else
    mkdir -p $TRACKHUB_DIR/$UCSC_GENOMEASSEMBLY
fi

if [ -z "$TRACKHUB_NAME" ]; then
    echo "[ERROR] trackhub name not specified"
    exit 1
fi

# get trackhub colors
COLORS=$(grep '^TRACKHUB_COLOR ' $CONFIG | cut -d ' ' -f 2,3 )

################################################################################

echo "
genome $UCSC_GENOMEASSEMBLY
trackDb $UCSC_GENOMEASSEMBLY/trackDb.txt
" > $TRACKHUB_DIR/genomes.txt

echo "
hub $(echo $TRACKHUB_NAME | sed -e 's/ /_/g')
shortLabel $TRACKHUB_NAME
longLabel $TRACKHUB_NAME
genomesFile genomes.txt
email $TRACKHUB_EMAIL
" > $TRACKHUB_DIR/hub.txt

PROJECT_RELPATH=$(python -c "import os.path; print os.path.relpath('$(pwd -P)',os.path.realpath('$(dirname $UCSC_GENOMEASSEMBLY/trackDb.txt)'))")
[ -z "$PROJECT_RELPATH" ] && PROJECT_RELPATH="."


################################################################################
# define functions for generating summary scaffold
#
# gatherDirs takes 1 parameter
# $1=TASK (e.g. $TASK_BWA)
function gatherDirs {
    vali=""
    for dir in ${DIR[@]}; do
        [ -d $OUT/${dir%%/*}/$1/ ] && vali=$vali" $OUT/${dir%%/*}/$1/"
    done
	echo $vali
}

# makeCompositeTrack takes 2 parameters
# $1=trackDb.txt folder
# $2=Libary folder
function makeCompositeTrack {
    cat >> $1/trackDb.txt <<DELIM
    
track $2
shortLabel $2
longLabel $2
compositeTrack on
showSubtrackColorOnUi on
viewLimits 0:3
configurable on
visibility full
priority 7
maxHeightPixels 64:64:11
type bed 3
subGroup1 view Views PK=Peaks SIG=Signals READS=Reads VC=VCF
subGroup2 mod Module\ #TODO aggregate onces
    mod=bowtie
sortOrder software=+ view=+
DELIM
} 

# makeTrack takes 2 parameters
# $1=trackDb.txt folder
# $2=Libary folder
# $3=subGroup1
# $4=type
# $5=visibility
function makeSubTrack {
    cat >> $1/trackDb.txt <<DELIM
    
    track $2_$3
    parent $2
    shortLabel $2 $3
    view $3
    type $4
    visibility $5
    viewUi on
DELIM
} 

# makeTrack takes 2 parameters
# $1=trackDb.txt folder
# $2=Libary folder
# $3=subGroup1  (e.g. LNCaP)
# $4=subGroup2  (e.g. bowtie)
# $5=type
# $6=sample pattern
# $7=filesuffix (e.g: .asd.bam)
# $8=additional track infos
function makeTracks {
    echo "ls $SOURCE/$2/$4/$6*$7"
    for f in $( ls $SOURCE/$2/$4/$6*$7); do
        RELPATH=$(python -c "import os.path; print os.path.relpath(os.path.realpath('$(dirname $f)'),'$1/$2_$3')")
        TRACKNAME=${f##*/}          # remove folders
        TRACKNAME=${TRACKNAME/%$7/} # remove file suffix
        TRACKNAME=$(echo $TRACKNAME | sed -n "s/$TRACKHUB_NAMEPATTERN/\1/p") # extract pattern

        COLOR=$(awk -v FILE=$f 'BEGIN{IGNORECASE=1} FILE~$1 {print $2}' <(echo "$COLORS") | tail -n 1)
        # default color
        if [ -z "$COLOR" ]; then COLOR="100,100,100"; fi
        
        mkdir -p $1/$2_$3
        [ -f $1/$2_$3/$4/${f##*/} ] && rm $1/$2_$3/$4/${f##*/}
        ln -s $RELPATH/${f##*/} $1/$2_$3/$4/${f##*/}
        
        cat >> $1/trackDb.txt <<DELIM
            track $2_$3_$4_${f##*/}
            shortLabel $2 $TRACKNAME 
            longLabel $TRACKNAME ($2 $3 $4)
            parent $2_$3 on
            type $5
            color $COLOR
            bigDataUrl $2_$3/$4/${f##*/}
            subGroups view=$3 mod=$4
            $9
            
DELIM
    done
} 

################################################################################
################################################################################
################################################################################

################################################################################
# make all subtracks first
SAMPLESETS=""
for dir in ${DIR[@]}; do
    # separate folder from sample pattern
    DIRNAME=${dir%%/*} 
    if [ -z "$SAMPLESETS" ]; then 
        SAMPLESETS="${DIRNAME}"    
    else
        SAMPLESETS="${SAMPLESETS[@]} ${DIRNAME}"
    fi
done

cat /dev/null > $TRACKHUB_DIR/$UCSC_GENOMEASSEMBLY/trackDb.txt

for SET in "$( echo $SAMPLESETS | tr ' ' '\n' | sort -u)"; do
    makeCompositeTrack $TRACKHUB_DIR/$UCSC_GENOMEASSEMBLY $SET
done

for dir in ${DIR[@]}; do

    # separate folder from sample pattern
    DIRNAME=${dir%%/*} # TODO ONLY ONCE PER 
    SAMPLEPATTERN=${dir/$DIRNAME/}
    
    ############################################################################
    # make signal composite (all bigwig and wig tracks)
    SUBGROUP1="SIG"
    makeSubTrack $TRACKHUB_DIR/$UCSC_GENOMEASSEMBLY $DIR $SUBGROUP1 "bigWig 0 100" "full"
    
    if [[ -n "$RUNBIGWIG" ]]; then        
        makeTracks $TRACKHUB_DIR/$UCSC_GENOMEASSEMBLY $DIR $SUBGROUP1 $INPUT_BIGWIG "bigWig 0 100" "$SAMPLEPATTERN" ".bw" ""
    fi

    if [[ -n "$RUNWIGGLER" ]]; then        
        makeTracks $TRACKHUB_DIR/$UCSC_GENOMEASSEMBLY $DIR $SUBGROUP1 $TASK_WIGGLER "bigWig 0 100" "$SAMPLEPATTERN" ".bw" ""
    fi
    

    ############################################################################
    # make peak composite (all bed, bg, bb)
    makeSubTrack $TRACKHUB_DIR/$UCSC_GENOMEASSEMBLY $DIR "PK" "bigBed 4" "dense"   
    SUBGROUP1="PK"
        
    if [[ -n "$RUNMACS2" ]]; then        
        makeTracks $TRACKHUB_DIR/$UCSC_GENOMEASSEMBLY $DIR $SUBGROUP1 $TASK_MACS2 "bigBed 4" "$SAMPLEPATTERN" ".bb" ""
    fi
    
    if [[ -n "$RUNPEAKRANGER" ]]; then        
        makeTracks $TRACKHUB_DIR/$UCSC_GENOMEASSEMBLY $DIR $SUBGROUP1 $TASK_PEAKRANGER "bed 4" "$SAMPLEPATTERN" "_region.bed" ""
    fi

    if [[ -n "$RUNHOMERCHIPSEQ" ]]; then        
        makeTracks $TRACKHUB_DIR/$UCSC_GENOMEASSEMBLY $DIR $SUBGROUP1 $TASK_HOMERCHIPSEQ "bed 4" "$SAMPLEPATTERN" ".bed" "" 
    fi
    
    ############################################################################
    # make read (all bams)
    makeSubTrack $TRACKHUB_DIR/$UCSC_GENOMEASSEMBLY $DIR "READ" "bam" "squish"   
    SUBGROUP1="READ"
        
    if [[ -n "$RUNMAPPINGBOWTIE" ]]; then        
        makeTracks $TRACKHUB_DIR/$UCSC_GENOMEASSEMBLY $DIR $SUBGROUP1 $TASK_BOWTIE "bam" "$SAMPLEPATTERN" ".$ASD.bam" "bamColorMode=strand"
    fi

    if [[ -n "$RUNMAPPINGBOWTIE2" ]]; then        
        makeTracks $TRACKHUB_DIR/$UCSC_GENOMEASSEMBLY $DIR $SUBGROUP1 $TASK_BOWTIE2 "bam" "$SAMPLEPATTERN" ".$ASD.bam" "bamColorMode=strand"
    fi

    if [[ -n "$RUNTOPHAT" ]]; then        
        makeTracks $TRACKHUB_DIR/$UCSC_GENOMEASSEMBLY $DIR $SUBGROUP1 $TASK_TOPHAT "bam" "$SAMPLEPATTERN" ".$ASD.bam" "bamColorMode=strand"
    fi

    if [[ -n "$RUNMAPPINGBWA" ]]; then        
        makeTracks $TRACKHUB_DIR/$UCSC_GENOMEASSEMBLY $DIR $SUBGROUP1 $TASK_BWA "bam" "$SAMPLEPATTERN" ".$ASD.bam" "bamColorMode=strand"
    fi

    if [[ -n "$RUNHICLIB" ]]; then        
        makeTracks $TRACKHUB_DIR/$UCSC_GENOMEASSEMBLY $DIR $SUBGROUP1 $TASK_HICLIB "bam" "$SAMPLEPATTERN" ".$ASD.bam" "bamColorMode=strand"
    fi

    if [[ -n "$RUNHICUP" ]]; then        
        makeTracks $TRACKHUB_DIR/$UCSC_GENOMEASSEMBLY $DIR $SUBGROUP1 $TASK_HICUP "bam" "$SAMPLEPATTERN" ".bam" "bamColorMode=strand"
    fi

    if [[ -n "$RUNREALRECAL" ]]; then        
        makeTracks $TRACKHUB_DIR/$UCSC_GENOMEASSEMBLY $DIR $SUBGROUP1 $TASK_RECAL "bam" "$SAMPLEPATTERN" ".$ASR.bam" "bamColorMode=strand"
    fi

    ############################################################################
    # make read (all vcg)
    makeSubTrack $TRACKHUB_DIR/$UCSC_GENOMEASSEMBLY $DIR "VC" "vcfTabix" "squish"   
    SUBGROUP1="VC"
    
    if [[ -n "$RUNPINDEL" ]]; then        
        makeTracks $TRACKHUB_DIR/$UCSC_GENOMEASSEMBLY $DIR $SUBGROUP1 $TASK_PINDEL "vcfTabix" "$SAMPLEPATTERN" "vcf.gz" ""
    fi

        
done
################################################################################


################################################################################
echo ">>>>> Generate UCSC trackhubs - FINISHED"
echo ">>>>> enddate "`date`
