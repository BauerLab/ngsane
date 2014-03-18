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

################################################################################

echo '''
genome $UCSC_GENOMEASSEMBLY
trackDb $UCSC_GENOMEASSEMBLY/trackDb.txt
''' > $TRACKHUB_DIR/genomes.txt

echo ''' 
hub $(echo $$TRACKHUB_NAME | sed -e 's/[ \t]*/d')
shortLabel $TRACKHUB_NAME
longLabel $TRACKHUB_NAME
genomesFile genomes.txt
email $TRACKHUB_EMAIL
''' > $TRACKHUB_DIR/hub.txt

mkdir -p $(dirname $SUMMARYTMP) && cat /dev/null > $SUMMARYTMP && cat /dev/null > $SUMMARYCITES # clean temporary content

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
    echo '''
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
''' >> $1/trackDb.txt
} 

# makeTrack takes 2 parameters
# $1=trackDb.txt folder
# $2=Libary folder
# $3=subGroup1
# $4=type
# $5=visibility
function makeSubTrack {
    echo '''
    track $2"_"$3
    parent $2
    shortLabel $2 $3
    view $3
    type $4
    visibility $5
    viewUi on
''' >> $1/trackDb.txt
} 

# makeTrack takes 2 parameters
# $1=trackDb.txt folder
# $2=Libary folder
# $3=subGroup1  (e.g. LNCaP)
# $4=subGroup2  (e.g. bowtie)
# $5=type
# $6=color
# $7=sample pattern
# $8=filesuffix (e.g: .asd.bam)
# $9=additional track infos
function makeTracks {
    for f in $( ls $SOURCE/$2/$4/$$7*$8); do
        RELPATH=$(python -c "import os.path; print os.path.relpath('$1',os.path.realpath('$(dirname $f)'))")
        TRACKNAME=${f##*/}          # remove folders
        TRACKNAME=${TRACKNAME/%$8/} # remove file suffix
        TRACKNAME=$(echo $TRACKNAME | sed -n "s/$TRACKHUB_NAMEPATTERN/\1/p") # extract pattern
        
        
        mkdir -p $1/$2_$3
        [ -f $1/$2_$3/${f##*/} ] && rm $1/$2_$3/${f##*/}
        ln -s $RELPATH/${f##*/} $1/$2_$3/${f##*/} $1/$2_$3/${f##*/}
        
        echo '''
            track $2"_"$3"_"$4
            shortLabel $2 $3 $4
            longLabel $TRACKNAME
            parent $2"_"$3 on
            type $5
            color $6
            bigDataUrl $1/$2_$3/${f##*/}
            subGroups view=$3 mod=$4
            $9
    ''' >> $1/trackDb.txt
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
    DATASETS="${SAMPLESETS[@]} ${DIRNAME}"
done

for SET in "$( echo $SAMPLESETS | sort -u); do
    makeCompositeTrack $UCSC_GENOMEASSEMBLY "$SET" 
done

for dir in ${DIR[@]}; do

    # separate folder from sample pattern
    DIRNAME=${dir%%/*} # TODO ONLY ONCE PER 
    SAMPLEPATTERN=${dir/$DIRNAME/}
    
    ############################################################################
    # make signal composite (all bigwig and wig tracks)
    SUBGROUP1="SIG"
    makeSubTrack $UCSC_GENOMEASSEMBLY $DIR $SUBGROUP1 "bigWig 0 100" "full"
    
    if [[ -n "$RUNBIGWIG" ]]; then        
        makeTracks $UCSC_GENOMEASSEMBLY $DIR $SUBGROUP1 $TASK_BIGWIG "bigWig 0 100" "255,214,79" $SAMPLEPATTERN ".bw" ""
    fi

    if [[ -n "$RUNWIGGLER" ]]; then        
        makeTracks $UCSC_GENOMEASSEMBLY $DIR $SUBGROUP1 $TASK_WIGGLER "bigWig 0 100" "255,214,79" $SAMPLEPATTERN ".bw" ""
    fi
    

    ############################################################################
    # make peak composite (all bed, bg, bb)
    makeSubTrack $UCSC_GENOMEASSEMBLY $DIR "PK" "bigBed 4" "dense"   
    
    if [[ -n "$RUNMACS2" ]]; then        
        makeTracks $UCSC_GENOMEASSEMBLY $DIR $SUBGROUP1 $TASK_MACS2 "bigBed 4" "255,214,79" $SAMPLEPATTERN ".bb" ""
    fi
    
    if [[ -n "$RUNPEAKRANGER" ]]; then        
        makeTracks $UCSC_GENOMEASSEMBLY $DIR $SUBGROUP1 $TASK_PEAKRANGER "bed 4" "255,214,79" $SAMPLEPATTERN "_region.bed" ""
    fi

    if [[ -n "$RUNHOMERCHIPSEQ" ]]; then        
        makeTracks $UCSC_GENOMEASSEMBLY $DIR $SUBGROUP1 $TASK_HOMERCHIPSEQ "bed 4" "255,214,79" $SAMPLEPATTERN ".bed" "" 
    fi
    
    ############################################################################
    # make read (all bams)
    makeSubTrack $UCSC_GENOMEASSEMBLY $DIR "READ" "bam" "squish"   
    
    if [[ -n "$RUNMAPPINGBOWTIE" ]]; then        
        makeTracks $UCSC_GENOMEASSEMBLY $DIR $SUBGROUP1 $TASK_BOWTIE "bam" "0,0,0" $SAMPLEPATTERN ".$ASD.bam" "bamColorMode=strand"
    fi

    if [[ -n "$RUNMAPPINGBOWTIE2" ]]; then        
        makeTracks $UCSC_GENOMEASSEMBLY $DIR $SUBGROUP1 $TASK_BOWTIE2 "bam" "0,0,0" $SAMPLEPATTERN ".$ASD.bam" "bamColorMode=strand"
    fi

    if [[ -n "$RUNTOPHAT" ]]; then        
        makeTracks $UCSC_GENOMEASSEMBLY $DIR $SUBGROUP1 $TASK_TOPHAT "bam" "0,0,0" $SAMPLEPATTERN ".$ASD.bam" "bamColorMode=strand"
    fi

    if [[ -n "$RUNMAPPINGBWA" ]]; then        
        makeTracks $UCSC_GENOMEASSEMBLY $DIR $SUBGROUP1 $TASK_BWA "bam" "0,0,0" $SAMPLEPATTERN ".$ASD.bam" "bamColorMode=strand"
    fi

    if [[ -n "$RUNHICLIB" ]]; then        
        makeTracks $UCSC_GENOMEASSEMBLY $DIR $SUBGROUP1 $TASK_HICLIB "bam" "0,0,0" $SAMPLEPATTERN ".$ASD.bam" "bamColorMode=strand"
    fi

    if [[ -n "$RUNHICUP" ]]; then        
        makeTracks $UCSC_GENOMEASSEMBLY $DIR $SUBGROUP1 $TASK_HICUP "bam" "0,0,0" $SAMPLEPATTERN ".bam" "bamColorMode=strand"
    fi

    if [[ -n "$RUNREALRECAL" ]]; then        
        makeTracks $UCSC_GENOMEASSEMBLY $DIR $SUBGROUP1 $TASK_RECAL "bam" "0,0,0" $SAMPLEPATTERN ".$ASR.bam" "bamColorMode=strand"
    fi

    ############################################################################
    # make read (all vcg)
    makeSubTrack $UCSC_GENOMEASSEMBLY $DIR "VC" "vcfTabix" "squish"   
    
    if [[ -n "$RUNPINDEL" ]]; then        
        makeTracks $UCSC_GENOMEASSEMBLY $DIR $SUBGROUP1 $TASK_PINDEL "vcfTabix" "0,0,0" $SAMPLEPATTERN "vcf.gz" ""
    fi

        
done
################################################################################


################################################################################
echo ">>>>> Generate UCSC trackhubs - FINISHED"
echo ">>>>> enddate "`date`
