#!/bin/bash -e

# Script running annovar
# QC:
# author: Denis C. Bauer
# date: Sept.2011

echo ">>>>> Annotate variants"
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"

function usage {
echo -e "usage: $(basename $0)

required:
  -k | --toolkit <path>     location of the NGSANE repository 
  -i | --input <file>       vcf file
  -o | --outdir <path>      output dir
  -r | --reference <file>   reference genome

options:
  -t | --threads <nr>       number of CPUs to use (default: 1)
"
exit
}


if [ ! $# -gt 3 ]; then usage ; fi

#DEFAULTS
THREADS=1
ADDRECALALN="" # additional commands for variant recalibrator

#INPUTS
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the NGSANE repository
        -t | --threads )        shift; THREADS=$1 ;; # number of CPUs to use
        -f | --file )           shift; INPUTFILE=$1 ;;  # input file
        -o | --outdir )         shift; OUTDIR=$1 ;;     # output dir 
        -n | --name )           shift; NAME=$1 ;; # name
        --maxGaussians )        shift; ADDRECALALN=$ADDRECALALN" --maxGaussians "$1 ;; #(additional params for recal)
        --percentBadVariants )  shift; ADDRECALALN=$ADDRECALALN" --percentBadVariants "$1 ;; #(additional params for recal)
        --recover-from )        shift; NGSANE_RECOVERFROM=$1 ;; # attempt to recover from log file                                                  
        -h | --help )           usage ;;
        * )                     usage
    esac
    shift
done



#PROGRAMS
. $CONFIG
. ${NGSANE_BASE}/conf/header.sh
. $CONFIG

################################################################################
NGSANE_CHECKPOINT_INIT "programs"

for MODULE in $MODULE_ANNOVAR; do module load $MODULE; done  # save way to load modules that itself load other modules
export PATH=$PATH_ANNOVAR:$PATH
module list
echo "PATH=$PATH"

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)


NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "parameters"

f=$INPUTFILE

NGSANE_CHECKPOINT_CHECK
################################################################################
NGSANE_CHECKPOINT_INIT "recall files from tape"

if [ -n "$DMGET" ]; then
	dmget -a $INPUTFILE
fi
    
NGSANE_CHECKPOINT_CHECK



# get basename of f
n=${f##*/}
mkdir -p $OUTDIR

################################################################################
NGSANE_CHECKPOINT_INIT "convert to annovar format"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

	#head -n 1000 $f >test

    command="convert2annovar.pl $f --format vcf4old --filter PASS --includeinfo >$OUTDIR/${n/vcf/txt}"
    echo $command && eval $command
    
    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/${n/vcf/txt}
fi

################################################################################
NGSANE_CHECKPOINT_INIT "autofilter"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    #summarize_annovar.pl --buildver $LIBRARY $OUTDIR/${n/vcf/txt} $DATABASE --remove -outfile $OUTDIR/${n/vcf/sum}
	command="table_annovar.pl $OUTDIR/${n/vcf/txt} $DATABASE --buildver $LIBRARY --protocol refGene,knownGene,ensGene,wgEncodeGencodeManualV4,gerp++elem,phastConsElements46way,genomicSuperDups,esp6500si_all,1000g2012apr_all,1000g2012apr_eur,1000g2012apr_amr,1000g2012apr_asn,1000g2012apr_afr,cg46,cosmic64,snp129,snp132,snp138,avsift,ljb2_all --operation g,g,g,g,r,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f --outfile $OUTDIR/${n/vcf/sum} --remove --otherinfo"
	#command="table_annovar.pl $OUTDIR/${n/vcf/txt} $DATABASE --buildver $LIBRARY --protocol refGene,genomicSuperDups,cosmic64 --operation g,r,f --outfile $OUTDIR/${n/vcf/sum} --remove --otherinfo"
	echo $command && eval $command

    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/${n/vcf/sum}.$LIBRARY"_multianno.txt" 
fi

################################################################################
NGSANE_CHECKPOINT_INIT "summarize over all samples"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

	# add the genotype information to the sample and make it then an csv file (cannot do it in the previous step because the genotype info is tab separated...)
    REP=$(grep "CHROM" $f) # | cut -f 10- | sed 's/\t/,/g')
    sed "s/Otherinfo/$REP/" $OUTDIR/${n/vcf/sum}.$LIBRARY"_multianno.txt" |	gawk -F \\t -v OFS="\",\"" '{$1=$1; print "\""$0"\""}' > $OUTDIR/${n/vcf/sum}.$LIBRARY"_multianno.csv"

# mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/${n/vcf/sum}.$LIBRARY"_multianno.csv" 
fi

################################################################################
NGSANE_CHECKPOINT_INIT "make statistics"

if [[ $(NGSANE_CHECKPOINT_TASK) == "start" ]]; then

    for i in nonsynonymous synonymous splicing; do
        echo $i $(grep $i $OUTDIR/${n/vcf/sum}.$LIBRARY"_multianno.txt" | wc -l) >> $OUTDIR/${n/vcf/sum}.$LIBRARY"_multianno.stats"
    done
    
    # mark checkpoint
    NGSANE_CHECKPOINT_CHECK $OUTDIR/${n/vcf/sum}.$LIBRARY"_multianno.stats" 
fi

#clean
rm $OUTDIR/${n/vcf/txt}

echo ">>>>> Annotate variants - FINISHED"
echo ">>>>> enddate "`date`


