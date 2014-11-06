#!/bin/bash -e

# Annotation of bam files
# author: Denis C. Bauer
# date: Jan.2013

# messages to look out for -- relevant for the QC.sh script:
# QCVARIABLES,Resource temporarily unavailable

echo ">>>>> Annotate BAM file "
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"


function usage {
echo -e "usage: $(basename $0) -k CONFIG -f BAM -o OUTDIR [OPTIONS]

Annotating BAM file with annotations in a folder specified in the config
file

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
        -f | --bam )            shift; FILES=$1 ;; # bam file
        -o | --outdir )         shift; OUTDIR=$1 ;; # output dir             
        --recover-from )        shift; NGSANE_RECOVERFROM=$1 ;; # attempt to recover from log file
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

# save way to load modules that itself loads other modules
hash module 2>/dev/null && for MODULE in $MODULE_BAMANN; do module load $MODULE; done && module list 

export PATH=$PATH_BAMANN:$PATH
echo "PATH=$PATH"

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
echo -e "--R           --\n "$(R --version | head -n 3)
[ -z "$(which R)" ] && echo "[ERROR] no R detected" && exit 1

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "parameters"

echo "[NOTE] Files: $FILES"
OLDFS=$IFS
IFS=","
DATASETS=""
for f in $FILES; do
    # get basename of f
    
    n=${f/%$ASD.bam/.anno.stats}
    FILE=${n/$INPUT_BAMANN/$TASK_BAMANN}
    # get directory
    d=$(dirname $f)
    d=${d##*/}    # add to dataset
    if [ -n "$FILE" ]; then 
        DATASETS="${DATASETS[@]} ${FILE[@]}"
    fi
done
IFS=" "
read -a DATASETS <<< "$DATASETS"

echo $DATASETS

BAMANNOUT=$OUTDIR/$(echo ${DIR[@]} | sed 's/ /_/g' | cut -c 1-60 )_${INPUT_BAMANN}.ggplot
if [ ! -f $BAMANNOUT ]; then mkdir -p $( dirname $BAMANNOUT); fi
BAMANNIMAGE=${BAMANNOUT/ggplot/pdf}


NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "recall files from tape"
	
if [ -n "$DMGET" ]; then
	dmget -a $TASK_BAMANN/*
fi
    
NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "aggregate data"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    # code for plotting
    cat > $OUTDIR/bamann.R <<EOF
library(ggplot2)
library(reshape)

args <- commandArgs(trailingOnly = TRUE)
pdf(file=args[2], width=15,height=7)

file=args[1]
x<-read.table(file, row.names = NULL ,header=T, quote="\"")
head(x)
distribution <- melt(x[,c(2:length(x))], id.vars="sample")
colnames(distribution)=c("sample","feature","value")

ggplot(distribution, aes(x = sample, y=value)) + 
  geom_bar(stat="identity", aes(fill = feature), position = "fill") + 
  scale_y_continuous("fraction") + 
  labs(title=args[3]) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+
  coord_flip()
#  facet_grid(. ~ type , space = "free", scales = "free_x")

ggplot(distribution, aes(x = sample, y=value)) + 
  geom_bar(stat="identity", aes(fill = feature)) +
  labs(title=args[3]) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  coord_flip()  
#  facet_grid(. ~ type , space = "free", scales = "free_x")


sink(type = "message")
sessionInfo()
EOF

    # add data
	head -n 1 ${DATASETS[0]} | gawk '{OFS=" "; print "type",$0,"sample"}' > $BAMANNOUT
    for i in ${DATASETS[@]}; do
        SAMPLE=${i/%.anno.stats/}
        SAMPLE=${SAMPLE/$OUT\//}
        SAMPLE=${SAMPLE/$TASK_BAMANN\//}
        grep --no-messages sum $i | gawk -F\\t -v x=$SAMPLE '{OFS=" "; print $0,x}';
	done >> $BAMANNOUT
	sed -i -r 's/\s+/ /g' $BAMANNOUT
	# remove mergedReads column
	cut -d' ' -f1,3- $BAMANNOUT > $BAMANNOUT.tmp && mv $BAMANNOUT.tmp $BAMANNOUT
	
	# plot
	RUNCOMMAND="Rscript $OUTDIR/bamann.R $BAMANNOUT $BAMANNIMAGE 'Genome Features $INPUT_BAMANN'"
	echo $RUNCOMMAND && eval $RUNCOMMAND

    convert $BAMANNIMAGE ${BAMANNIMAGE/pdf/jpg}

	# mark checkpoint
    NGSANE_CHECKPOINT_CHECK $BAMANNIMAGE
fi
################################################################################
echo ">>>>> Annotate BAM file - FINISHED"
echo ">>>>> enddate "`date`

